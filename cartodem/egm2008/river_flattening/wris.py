"""Pinned WRIS river polygons, excluding sand and islands on the native DEM grid.

ArcGIS's NWIC/River_GangaBasin/MapServer/11 renderer defines ripcode labels.
The national WaterBodies/MapServer/9 renderer alone does not distinguish them.
An on-disk R-tree avoids repeatedly loading the national geometry into memory.
Invalid linework can be normalized only after proving identical membership on
the native raster and less than one pixel of area change. No simplification,
reprojection or height substitution is done.
"""
from collections import Counter
import hashlib
import inspect
import json
from pathlib import Path
import sqlite3

import numpy as np
import pyarrow
import pyarrow.parquet as pq
import rasterio
from rasterio.features import rasterize
from rasterio.transform import array_bounds
import shapely
from shapely.geometry import box

from build import atomic_json, digest_json, download, read_json, sha256


CLASS_NAMES = {'01': 'River', '02': 'Sand', '03': 'Island', ' ': 'River', None: 'River'}


def prepare_source(root):
    lock = read_json(Path(__file__).with_name('wris.lock.json'))
    folder = root / 'inputs/wris'
    source = folder / 'WRIS_River_Polygons.parquet'
    record = download(lock['url'], source, expected_size=lock['bytes'], expected_sha=lock['sha256'])
    record['legend_url'] = lock['legend_url']
    # Reusing verified bytes should not erase the original HTTP provenance.
    report = folder / 'source.json'
    previous = read_json(report) if report.exists() else {}
    atomic_json(report, {**previous, **record})
    return source


def build_index(source, target):
    """Index full WGS84 WKB once; preserve attributes and interior rings."""
    source, target = Path(source), Path(target)
    recipe = dict(source_sha256=sha256(source), code_sha256=digest_json(inspect.getsource(build_index)),
                  pyarrow=pyarrow.__version__,
                  shapely=shapely.__version__)
    report_path = target.with_suffix('.json')
    if target.exists():
        report = read_json(report_path)
        if report['recipe'] != recipe or report['sha256'] != sha256(target):
            raise ValueError('WRIS index recipe or bytes changed; use a new index path')
        return report
    target.parent.mkdir(parents=True, exist_ok=True)
    temporary = target.with_name(target.name + '.part')
    temporary.unlink(missing_ok=True)
    parquet = pq.ParquetFile(source)
    geo = json.loads(parquet.schema_arrow.metadata[b'geo'])
    column = geo['primary_column']
    metadata = geo['columns'][column]
    # GeoParquet's absent CRS defaults to OGC:CRS84. An explicitly different
    # CRS must be implemented and tested, never silently interpreted as degrees.
    if metadata.get('encoding') != 'WKB' or 'crs' in metadata:
        raise ValueError('Expected the pinned WKB GeoParquet with implicit CRS84')
    counts = Counter(); invalid = []; total = 0
    with sqlite3.connect(temporary) as db:
        db.execute('CREATE TABLE features (id INTEGER PRIMARY KEY, objectid INTEGER UNIQUE, code TEXT, geometry BLOB)')
        db.execute('CREATE VIRTUAL TABLE bounds USING rtree(id, west, east, south, north)')
        for batch in parquet.iter_batches(batch_size=128, columns=['objectid', 'ripcode', column]):
            rows = batch.to_pydict()
            for oid, code, wkb in zip(rows['objectid'], rows['ripcode'], rows[column]):
                if code not in CLASS_NAMES:
                    raise ValueError(f'Unrecognised WRIS ripcode {code!r} at OBJECTID {oid}')
                geometry = shapely.from_wkb(wkb)
                if geometry is None or geometry.is_empty or geometry.geom_type not in ('Polygon', 'MultiPolygon'):
                    raise ValueError(f'Missing/non-polygon WRIS geometry at OBJECTID {oid}')
                west, south, east, north = geometry.bounds
                if not np.isfinite(geometry.bounds).all() or not (-180 <= west <= east <= 180 and -90 <= south <= north <= 90):
                    raise ValueError(f'Invalid WGS84 bounds at OBJECTID {oid}')
                if not geometry.is_valid:
                    invalid.append(int(oid))
                total += 1; counts[code] += 1
                db.execute('INSERT INTO features VALUES (?,?,?,?)', (total, int(oid), code, wkb))
                db.execute('INSERT INTO bounds VALUES (?,?,?,?,?)', (total, west, east, south, north))
        db.commit()
    temporary.replace(target)
    report = dict(recipe=recipe, sha256=sha256(target), features=total,
                  classes=[dict(code=code, label=CLASS_NAMES[code], features=count) for code, count in counts.items()],
                  invalid_geometry_objectids=invalid)
    atomic_json(report_path, report)
    return report


def normalized_polygon(geometry, shape, transform):
    """Accept only a linework repair that leaves native pixel membership intact.

    Self-touching rings occur in the pinned release. GEOS linework make_valid
    keeps their edges; polygonal parts are retained and collapsed lines dropped.
    Both the area check and the pixel-centre comparison must pass. This is not
    permission to repair arbitrary source topology or move water boundaries.
    """
    repaired = shapely.make_valid(geometry, method='linework')
    polygons = []
    def collect(part):
        if part.geom_type == 'Polygon':
            polygons.append(part)
        elif hasattr(part, 'geoms'):
            for child in part.geoms:
                collect(child)
    collect(repaired)
    if not polygons:
        raise ValueError('Invalid WRIS polygon has no polygonal repair')
    repaired = shapely.union_all(polygons)
    pixel_area = abs(transform.a*transform.e-transform.b*transform.d)
    if not repaired.is_valid or abs(repaired.area-geometry.area) > pixel_area:
        raise ValueError('Invalid WRIS repair changes more than one native pixel of area')
    original_mask = rasterize([(geometry, 1)], out_shape=shape, transform=transform,
                              fill=0, all_touched=False, dtype='uint8')
    repaired_mask = rasterize([(repaired, 1)], out_shape=shape, transform=transform,
                              fill=0, all_touched=False, dtype='uint8')
    if not np.array_equal(original_mask, repaired_mask):
        raise ValueError('Invalid WRIS repair changes native pixel membership')
    return repaired, dict(reason=shapely.is_valid_reason(geometry),
                          original_wkb_sha256=hashlib.sha256(geometry.wkb).hexdigest(),
                          repaired_wkb_sha256=hashlib.sha256(repaired.wkb).hexdigest(),
                          native_membership_changes=0,
                          area_change_pixels=float(abs(repaired.area-geometry.area)/pixel_area))


def river_mask(features, shape, transform, *, normalize_invalid=False):
    """Pixel-centre membership. Explicit sand/island exclusions override water."""
    water, excluded = [], []
    counts = Counter(); invalid = []; normalized = []
    region = box(*array_bounds(*shape, transform))
    for oid, code, geometry in features:
        if code not in CLASS_NAMES:
            raise ValueError(f'Unrecognised WRIS ripcode {code!r}')
        if geometry.is_empty or not geometry.intersects(region):
            continue
        if not geometry.is_valid:
            if not normalize_invalid:
                invalid.append(int(oid))
                continue
            geometry, detail = normalized_polygon(geometry, shape, transform)
            normalized.append(dict(objectid=int(oid), **detail))
        counts[CLASS_NAMES[code]] += 1
        (water if CLASS_NAMES[code] == 'River' else excluded).append((geometry, 1))
    if invalid:
        raise ValueError(f'Invalid WRIS polygons intersect this sheet: {invalid}; inspect source geometry')
    def burn(shapes):
        return (rasterize(shapes, out_shape=shape, transform=transform, fill=0,
                          all_touched=False, dtype='uint8').astype(bool)
                if shapes else np.zeros(shape, dtype=bool))
    selected = burn(water); protected = burn(excluded)
    overlap = int((selected & protected).sum())
    selected &= ~protected
    return selected, dict(features=dict(counts), excluded_overlap_pixels=overlap,
                          normalized_invalid_polygons=normalized,
                          river_pixels=int(selected.sum()))


def sheet_mask(index, glo_path, output):
    with rasterio.open(glo_path) as src:
        original = src.read(1)
        shape, transform, profile = original.shape, src.transform, src.profile.copy()
        west, south, east, north = src.bounds
    with sqlite3.connect(f'file:{index}?mode=ro', uri=True) as db:
        rows = db.execute('SELECT f.objectid,f.code,f.geometry FROM features f JOIN bounds b ON f.id=b.id '
                          'WHERE b.west<=? AND b.east>=? AND b.south<=? AND b.north>=? ORDER BY f.id',
                          (east, west, north, south))
        river, report = river_mask(((oid, code, shapely.from_wkb(wkb)) for oid, code, wkb in rows),
                                  shape, transform, normalize_invalid=True)
    # GLO-30 supplies only the ocean boundary. Its river/lake classifications do
    # not add to or subtract from WRIS River. Everything outside the selected
    # river footprint retains the existing CartoDEM height, including lakes.
    classes = np.zeros_like(original)
    classes[original == 1] = 1
    classes[original > 3] = 255
    eligible = (original != 1) & (original <= 3)
    classes[river & eligible] = 3
    report.update(ocean_overlap_pixels=int((river & (original == 1)).sum()),
                  diagnostic_glo_lake_overlap_pixels=int((river & (original == 2)).sum()),
                  missing_classification_overlap_pixels=int((river & (original > 3)).sum()),
                  selected_river_pixels=int((classes == 3).sum()))
    if report['missing_classification_overlap_pixels']:
        raise ValueError('WRIS river intersects missing ocean classification')
    profile.update(driver='GTiff', count=1, dtype='uint8', nodata=255, compress='deflate',
                   tiled=True, blockxsize=256, blockysize=256)
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = output.with_name(output.name + '.part')
    with rasterio.open(temporary, 'w', **profile) as dst:
        dst.write(classes, 1)
        dst.update_tags(RIVER_MASK='WRIS: 01, blank and null; Sand 02 and Island 03 excluded',
                        OCEAN_MASK='GLO-30 class 1', NONRIVER_ELEVATIONS='Unchanged')
    temporary.replace(output)
    return report


def prepare_masks(root, work, selected):
    source = prepare_source(root)
    index_id=digest_json(dict(source=sha256(source),code=inspect.getsource(build_index),
                             pyarrow=pyarrow.__version__,shapely=shapely.__version__))[:16]
    index = root / 'inputs/wris' / f'polygons-{index_id}.sqlite'
    source_report = build_index(source, index)
    lock = read_json(Path(__file__).with_name('wris.lock.json'))
    expected = {r['code']: r['features'] for r in lock['classes']}
    if {r['code']: r['features'] for r in source_report['classes']} != expected:
        raise ValueError('WRIS class inventory differs from the pinned release')
    reports = []; paths = {}
    for sheet in selected:
        if sheet['mask'] is None:
            raise ValueError('WRIS mask requires verified ocean classification for each sheet')
        source_mask = root / 'inputs/glo30-wbm' / sheet['mask']
        output = work / 'masks' / sheet['name']
        report_path = output.with_suffix('.json')
        recipe = dict(index_sha256=source_report['sha256'], glo_sha256=sha256(source_mask),
                      code_sha256=sha256(Path(__file__)), gdal=rasterio.__gdal_version__,
                      shape=[sheet['height'], sheet['width']], transform=sheet['transform'])
        if output.exists():
            report = read_json(report_path)
            if report['recipe'] != recipe or report['sha256'] != sha256(output):
                raise ValueError('WRIS mask cache changed; use a different --label')
        else:
            report = dict(sheet=sheet['name'], recipe=recipe, **sheet_mask(index, source_mask, output))
            report['sha256'] = sha256(output)
            atomic_json(report_path, report)
        paths[sheet['name']] = output
        reports.append(report)
    recipe = dict(source_sha256=lock['sha256'], legend_url=lock['legend_url'],
                  masks={r['sheet']: r['sha256'] for r in reports})
    atomic_json(work / 'wris-masks.json', dict(recipe=recipe, sheets=reports, source=source_report))
    return paths, recipe
