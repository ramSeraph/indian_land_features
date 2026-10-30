#!/usr/bin/env python3
"""Download CartoDEM inputs and prepare EGM2008 ocean-masked sheets.

All paths are relative to --workdir. Existing verified inputs can be reused.
GLO-30 supplies ocean and river classification across the DEM footprint; elevation
observations remain CartoDEM. Paths are independent of the execution environment.
"""
import os
os.environ.setdefault("GDAL_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("OMP_NUM_THREADS", "1")

import argparse
import csv
import hashlib
import json
import math
import shutil
import time
import urllib.error
import urllib.request
import xml.etree.ElementTree as ET
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

import numpy as np
import rasterio
from rasterio.transform import Affine

from cartodem_to_egm2008 import GEOID_SHA256, convert, sha256, windows

RELEASE = "https://github.com/ramSeraph/indian_land_features/releases/download/cartodem-30m-v3r1-orig/"
GEOID_URL = "https://cdn.proj.org/us_nga_egm08_25.tif"
WBM_BASE = "https://copernicus-dem-30m.s3.amazonaws.com/"
SOURCE_COUNT = 543
SOURCE_BYTES = 14466784647
SOURCE_LIST_SHA256 = '233b812f228fa8c96ddb1ceca6f5d3e1a0eaeecec3526413bebac4bf1369de88'
NODATA = -255.0


def atomic_npz(path,**arrays):
    path=Path(path)
    path.parent.mkdir(parents=True,exist_ok=True)
    temp=path.with_name(path.name+'.part')
    with temp.open('wb') as stream:
        np.savez_compressed(stream,**arrays)
    temp.replace(path)


def atomic_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    data = json.dumps(value, indent=2, sort_keys=True, allow_nan=False) + "\n"
    temp = path.with_name(path.name + ".part")
    temp.write_text(data)
    temp.replace(path)


def digest_json(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, allow_nan=False).encode()).hexdigest()


def emit(event, **fields):
    print(json.dumps(dict(event=event, **fields), allow_nan=False), flush=True)


def read_json(path):
    return json.loads(Path(path).read_text())


def setup(root):
    root = Path(root).absolute()
    for name in ("inputs/cartodem", "inputs/glo30-wbm", "inputs/geoid", "sheets/coastal",
                 "mosaics", "tiles", "output", "reports/sheets", "logs", "tmp"):
        (root / name).mkdir(parents=True, exist_ok=True)
    return root


def download(url, path, *, expected_size=None, expected_sha=None):
    """Download atomically, checking exact bytes; never replace a bad existing input."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    headers = {}
    reused = path.exists()
    if not reused:
        temp = path.with_name(path.name + ".download")
        for attempt in range(4):
            try:
                request = urllib.request.Request(url, headers={"User-Agent": "CartoDEM-EGM2008-build/1"})
                with urllib.request.urlopen(request, timeout=90) as response, temp.open("wb") as out:
                    shutil.copyfileobj(response, out, length=1024*1024)
                    out.flush()
                    headers = {k: response.headers.get(k) for k in ("ETag", "Last-Modified")}
                    length = response.headers.get("Content-Length")
                    if length is not None and temp.stat().st_size != int(length):
                        raise IOError("Incomplete HTTP response")
                if expected_size is not None and temp.stat().st_size != expected_size:
                    raise ValueError(f"Downloaded size differs for {path.name}")
                if expected_sha is not None and sha256(temp) != expected_sha:
                    raise ValueError(f"Downloaded checksum differs for {path.name}")
                temp.replace(path)
                break
            except urllib.error.HTTPError as exc:
                if exc.code in (403, 404) or attempt == 3:
                    temp.unlink(missing_ok=True)
                    raise
                time.sleep(2**attempt)
            except Exception:
                if attempt == 3:
                    temp.unlink(missing_ok=True)
                    raise
                time.sleep(2**attempt)
    size = path.stat().st_size
    digest = sha256(path)
    if expected_size is not None and size != expected_size:
        raise ValueError(f"Existing input size differs: {path.name}")
    if expected_sha is not None and digest != expected_sha:
        raise ValueError(f"Existing input checksum differs: {path.name}")
    return dict(url=url, bytes=size, sha256=digest, **{k: v for k, v in headers.items() if v})


def mask_name(lon, lat):
    lat, lon = math.floor(lat), math.floor(lon)
    ns, ew = ("N" if lat >= 0 else "S"), ("E" if lon >= 0 else "W")
    stem = f"Copernicus_DSM_COG_10_{ns}{abs(lat):02d}_00_{ew}{abs(lon):03d}_00"
    filename = stem + "_WBM.tif"
    return filename, WBM_BASE + stem + "_DEM/AUXFILES/" + filename


def fetch_inputs(root, workers=4):
    previous_path = root / "reports/inputs.json"
    pinned_path = Path(__file__).resolve().parent / 'provenance/inputs.lock.json'
    previous = (read_json(previous_path) if previous_path.exists() else
                read_json(pinned_path) if pinned_path.exists() else {})
    old_sources = {x["name"]: x for x in previous.get("sheets", [])}
    old_masks = {x["name"]: x for x in previous.get("masks", [])}
    index_path = root / "inputs/cartodem/listing_files.csv"
    index = download(RELEASE + "listing_files.csv", index_path,
                     expected_sha=SOURCE_LIST_SHA256)
    with index_path.open(newline="") as stream:
        rows = [r for r in csv.DictReader(stream) if r["name"].endswith(".tif")]
    if len(rows) != SOURCE_COUNT or sum(int(r["size"]) for r in rows) != SOURCE_BYTES:
        raise ValueError("Release inventory changed; inspect before rebuilding")
    if len({r['name'] for r in rows}) != len(rows):
        raise ValueError("Duplicate source names in release inventory")
    geoid_path = root / "inputs/geoid/us_nga_egm08_25.tif"
    geoid = download(GEOID_URL, geoid_path, expected_sha=GEOID_SHA256)
    sheets = []
    masks_needed = {}
    for i, row in enumerate(rows):
        name = row["name"]
        if Path(name).name != name or row["url"] != RELEASE + name:
            raise ValueError("Unexpected release asset path or URL")
        source = root / "inputs/cartodem" / name
        entry = download(row["url"], source, expected_size=int(row["size"]),
                         expected_sha=old_sources.get(name, {}).get("sha256"))
        with rasterio.open(source) as src:
            if (src.crs != rasterio.crs.CRS.from_epsg(4326) or src.count != 1
                    or src.dtypes != ("int16",)):
                raise ValueError(f"Unexpected original raster: {name}")
            center = src.xy(src.height // 2, src.width // 2)
            mask, url = mask_name(*center)
            entry.update(name=name, width=src.width, height=src.height,
                         transform=list(src.transform)[:6], bounds=list(src.bounds), mask=mask)
        masks_needed[mask] = url
        sheets.append(entry)
        if i % 20 == 0 or i + 1 == len(rows):
            emit("inventory", complete=i+1, total=len(rows))

    mask_entries = []
    missing = []
    failures = []
    def get_mask(name, url):
        try:
            entry = download(url, root / "inputs/glo30-wbm" / name,
                             expected_sha=old_masks.get(name, {}).get("sha256"))
            # Preserve HTTP provenance across reuse.
            entry = {**old_masks.get(name, {}), **entry, "name": name}
            return entry
        except urllib.error.HTTPError as exc:
            if exc.code == 404:
                return dict(name=name, url=url, unavailable=True)
            raise
    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = {pool.submit(get_mask, name, url):name for name, url in masks_needed.items()}
        for i, future in enumerate(as_completed(futures)):
            try:
                entry = future.result()
                (missing if entry.get("unavailable") else mask_entries).append(entry)
            except Exception as exc:
                failure=dict(name=futures[future],error=str(exc))
                failures.append(failure)
                emit('download_error',**failure)
            if i % 20 == 0 or i + 1 == len(futures):
                emit("masks", complete=i+1, total=len(futures), unavailable=len(missing),failed=len(failures))
    if failures:
        atomic_json(root/'reports/download-failures.json',failures)
        raise ValueError(f'{len(failures)} mask downloads failed; see reports/download-failures.json')
    # A missing product is not proof of ocean. It is harmless only if CartoDEM
    # has no usable elevation in that entire sheet, so nothing will be changed.
    missing_names = {e['name'] for e in missing}
    for sheet in sheets:
        if sheet['mask'] not in missing_names:
            continue
        with rasterio.open(root / 'inputs/cartodem' / sheet['name']) as src:
            for window in windows(src.width, src.height):
                a = src.read(1, window=window)
                if np.any((src.read_masks(1, window=window) != 0) & ~np.isin(a, [-255,-1024])):
                    raise ValueError(f"Missing GLO-30 mask at valid CartoDEM sheet {sheet['name']}")
        sheet['mask'] = None
        sheet['mask_note'] = 'Unavailable product; entire source sheet is already nodata'
    result = dict(schema=1, release=RELEASE, listing=index, geoid=geoid,
                  sheets=sorted(sheets, key=lambda x:x['name']),
                  masks=sorted(mask_entries, key=lambda x:x['name']),
                  unavailable_masks=sorted(missing, key=lambda x:x['name']))
    atomic_json(previous_path, result)
    return result


def prepare(root, names=None):
    inputs = read_json(root / "reports/inputs.json")
    masks = {e['name']:e for e in inputs['masks']}
    recipe = dict(converter_sha256=sha256(Path(__file__).with_name('cartodem_to_egm2008.py')),
                  output_nodata=NODATA, ocean_class=1, river_flattening=False,
                  source_nodata=[-255,-1024], gdal=rasterio.__gdal_version__,
                  proj=rasterio.__proj_version__)
    recipe_id = digest_json(recipe)
    selected = [s for s in inputs['sheets'] if names is None or s['name'] in names]
    if names and len(selected) != len(set(names)):
        raise ValueError('A requested sheet is absent from the input inventory')
    reports = []
    for i, sheet in enumerate(selected):
        started = time.monotonic()
        name = sheet['name']
        source = root / 'inputs/cartodem' / name
        mask = root / 'inputs/glo30-wbm' / sheet['mask'] if sheet['mask'] else None
        output = root / 'sheets/coastal' / name
        report_path = root / 'reports/sheets' / (Path(name).stem + '.json')
        mask_hash = masks[sheet['mask']]['sha256'] if mask else None
        expected = dict(source_sha256=sheet['sha256'], water_mask_sha256=mask_hash,
                        geoid_sha256=inputs['geoid']['sha256'], recipe_sha256=recipe_id)
        if sha256(source) != expected['source_sha256']:
            raise ValueError(f'Source changed: {name}')
        if mask and sha256(mask) != mask_hash:
            raise ValueError(f'Mask changed: {mask.name}')
        if output.exists():
            if not report_path.exists():
                raise ValueError(f'Output without completion report: {name}')
            report = read_json(report_path)
            if any(report.get(k) != v for k,v in expected.items()):
                raise ValueError(f'Existing sheet was built with different inputs/settings: {name}')
            if sha256(output) != report['output_sha256']:
                raise ValueError(f'Prepared sheet checksum changed: {name}')
            reused = True
        else:
            report = convert(source, root/'inputs/geoid/us_nga_egm08_25.tif', output,
                             water_mask_path=mask, output_nodata=NODATA)
            report.update(expected)
            report.update(source=f'inputs/cartodem/{name}', output=f'sheets/coastal/{name}',
                          water_mask=f'inputs/glo30-wbm/{mask.name}' if mask else None,
                          recipe=recipe, output_sha256=sha256(output),
                          output_bytes=output.stat().st_size,
                          elapsed_seconds=round(time.monotonic()-started,3))
            atomic_json(report_path, report)
            reused = False
        reports.append(report)
        emit('prepare', sheet=name, complete=i+1, total=len(selected), reused=reused,
             seconds=round(time.monotonic()-started,2),
             valid=report['valid_pixels'], ocean_removed=report['newly_masked_ocean_pixels'])
        atomic_json(root/'reports/prepare-status.json', dict(complete=i+1,total=len(selected),
                    last_sheet=name,state='complete' if i+1 == len(selected) else 'running'))
    summary = dict(recipe=recipe, sheets=len(reports),
                   valid_pixels=sum(r['valid_pixels'] for r in reports),
                   source_nodata_pixels=sum(r['source_nodata_pixels'] for r in reports),
                   newly_masked_ocean_pixels=sum(r['newly_masked_ocean_pixels'] for r in reports),
                   source_nodata_on_glo30_land_pixels=sum(r['source_nodata_on_glo30_land_pixels'] for r in reports),
                   output_bytes=sum(r['output_bytes'] for r in reports),
                   source_fill_counts={k:sum(r['source_fill_counts'][k] for r in reports)
                                       for k in ('-255','-1024')})
    atomic_json(root/('reports/prepare-summary.json' if names is None else 'reports/pilot-prepare-summary.json'), summary)
    return reports


def create_mosaic(paths, height_path, validity_path):
    """Write portable VRTs on the exact common pixel grid, with 0/255 validity."""
    paths = [Path(p) for p in paths]
    records = []
    for p in paths:
        with rasterio.open(p) as src:
            if (src.count != 1 or src.crs != rasterio.crs.CRS.from_epsg(4326)
                    or src.dtypes != ('float32',) or src.nodata != NODATA):
                raise ValueError(f'Unexpected prepared sheet {p.name}')
            records.append((p, src.transform, src.width, src.height))
    if not records:
        raise ValueError('No sheets to mosaic')
    dx,dy=records[0][1].a,records[0][1].e
    west=min(r[1].c for r in records)
    north=max(r[1].f for r in records)
    placed=[]
    footprints=set()
    for p,t,w,h in records:
        if abs(t.a-dx)>1e-12 or abs(t.e-dy)>1e-12 or t.b or t.d:
            raise ValueError('Prepared sheets do not share a north-up resolution')
        x,y=(t.c-west)/dx,(t.f-north)/dy
        if abs(x-round(x))>1e-5 or abs(y-round(y))>1e-5:
            raise ValueError('Prepared sheets are not aligned on the same pixel grid')
        x,y=round(x),round(y)
        if (x,y,w,h) in footprints:
            raise ValueError('Duplicate mosaic footprint')
        footprints.add((x,y,w,h))
        placed.append((p,x,y,w,h))
    width=max(x+w for p,x,y,w,h in placed)
    height=max(y+h for p,x,y,w,h in placed)
    transform=Affine(dx,0,west,0,dy,north)
    # Materialize native 0/255 validity as DATA, with no nodata or overviews.
    # An implicit nodata mask exposed through SourceBand=mask,1 can inherit
    # numeric overview/resampling behaviour and erase small holes. The explicit
    # mask is small, compresses well and has unambiguous nearest-pixel semantics.
    mask_dir=Path(validity_path).parent/'validity'
    mask_dir.mkdir(parents=True,exist_ok=True)
    mask_paths={}
    for p,x,y,w,h in placed:
        mask_path=mask_dir/(p.stem+'-validity.tif')
        temp=mask_path.with_name(mask_path.name+'.part')
        with rasterio.open(p) as src:
            profile=dict(driver='GTiff',width=w,height=h,count=1,dtype='uint8',
                         crs=src.crs,transform=src.transform,tiled=True,
                         blockxsize=512,blockysize=512,compress='DEFLATE')
            with rasterio.open(temp,'w',**profile) as dst:
                for window in windows(w,h):
                    dst.write(src.read_masks(1,window=window),1,window=window)
        temp.replace(mask_path)
        mask_paths[p]=mask_path
    for output,validity in ((Path(height_path),False),(Path(validity_path),True)):
        output.parent.mkdir(parents=True,exist_ok=True)
        vrt=ET.Element('VRTDataset',rasterXSize=str(width),rasterYSize=str(height))
        ET.SubElement(vrt,'SRS',dataAxisToSRSAxisMapping='2,1').text='EPSG:4326'
        ET.SubElement(vrt,'GeoTransform').text=','.join(map(str,transform.to_gdal()))
        band=ET.SubElement(vrt,'VRTRasterBand',dataType='Byte' if validity else 'Float32',band='1')
        if not validity:
            ET.SubElement(band,'NoDataValue').text=str(NODATA)
        for p,x,y,w,h in placed:
            source=ET.SubElement(band,'SimpleSource' if validity else 'ComplexSource')
            source_path=mask_paths[p] if validity else p
            ET.SubElement(source,'SourceFilename',relativeToVRT='1').text=os.path.relpath(source_path,output.parent)
            ET.SubElement(source,'SourceBand').text='1'
            ET.SubElement(source,'SrcRect',xOff='0',yOff='0',xSize=str(w),ySize=str(h))
            ET.SubElement(source,'DstRect',xOff=str(x),yOff=str(y),xSize=str(w),ySize=str(h))
            if not validity:
                ET.SubElement(source,'NODATA').text=str(NODATA)
        ET.indent(vrt)
        output.write_bytes(ET.tostring(vrt,encoding='utf-8',xml_declaration=True))
    return dict(width=width,height=height,transform=list(transform)[:6],
                bounds=[west,north+height*dy,west+width*dx,north],sheets=len(records))


def main():
    import fcntl
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('stage',choices=['download','prepare'])
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--download-workers',type=int,default=4)
    parser.add_argument('--sheets',nargs='+',help='Source filenames; omit for the full inventory')
    args=parser.parse_args()
    root=setup(args.workdir)
    with (root/'logs/build.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        if args.stage=='download':fetch_inputs(root,args.download_workers)
        else:prepare(root,args.sheets)


if __name__=='__main__':main()
