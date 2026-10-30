#!/usr/bin/env python3
"""Verify pmtiles-mosaic CLI output against encoded-tile checksums.

This reads existing PMTiles; it does not write or repack any PMTiles file.
Run in the pmtiles-mosaic environment to record its exact installed versions.
"""
import argparse
from contextlib import closing
from datetime import datetime,timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import sqlite3

from pmtiles.reader import Reader,all_tiles
from pmtiles.tile import zxy_to_tileid


def sha256(path):
    with path.open('rb') as stream:
        return hashlib.file_digest(stream,'sha256').hexdigest()


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',required=True,type=Path)
    parser.add_argument('--pmtiles-dir',required=True,type=Path)
    parser.add_argument('--label',default='display')
    args=parser.parse_args()
    root=args.workdir.resolve();out=args.pmtiles_dir.resolve()
    original=json.loads((root/f'output/{args.label}/archives.json').read_text())
    report=json.loads((root/f'reports/{args.label}-tiles.json').read_text())
    expected={}
    for zoom in report['zooms']:
        db_path=root/f"tiles/{args.label}/z{zoom['zoom']}.sqlite"
        with closing(sqlite3.connect(db_path.as_uri()+'?mode=ro',uri=True)) as db:
            expected.update(db.execute('SELECT tileid,sha256 FROM tiles'))
    total=len(expected);parts=[];index={}
    for path in sorted(out.glob('*.pmtiles')):
        count=0;tile_sizes={}
        with path.open('rb') as stream:
            def read_at(offset,length):
                stream.seek(offset);return stream.read(length)
            h=Reader(read_at).header()
            if h['tile_type'].value!=4 or h['tile_compression'].value!=1:
                raise ValueError('Expected WebP without extra tile compression')
            for (z,x,y),data in all_tiles(read_at):
                if expected.pop(zxy_to_tileid(z,x,y),None)!=hashlib.sha256(data).hexdigest():
                    raise ValueError(f'Changed, duplicate or unexpected tile: {z}/{x}/{y}')
                count+=1
                tile_sizes[f'{z},{x},{y}']=len(data)
        if path.stat().st_size>=2**31:raise ValueError('Archive exceeds GitHub file size limit')
        entry=dict(file=path.name,bytes=path.stat().st_size,sha256=sha256(path),tiles=count,
                   verified_tiles=count,minzoom=h['min_zoom'],maxzoom=h['max_zoom'],
                   bounds=[h[k]/1e7 for k in ('min_lon_e7','min_lat_e7','max_lon_e7','max_lat_e7')])
        parts.append(entry)
        west,south,east,north=entry['bounds']
        index[path.name]=dict(tiles=tile_sizes,min_zoom=h['min_zoom'],max_zoom=h['max_zoom'],
                             bounds=[south,west,north,east])
        print(json.dumps(dict(event='verified_partition',**entry)),flush=True)
    if expected:raise ValueError(f'{len(expected)} encoded tiles are missing from the PMTiles')
    validation=out/'.validation';validation.mkdir(exist_ok=True)
    (validation/'partition-info.json').write_text(json.dumps(index,separators=(',',':'))+'\n')
    software={d.metadata['Name']:d.version for d in importlib.metadata.distributions()}
    (out/'requirements-mosaic.txt').write_text(''.join(f'{k}=={v}\n' for k,v in sorted(software.items())))
    record=dict(verified=datetime.now(timezone.utc).isoformat(),software=software,
                verifier_sha256=sha256(Path(__file__)),verified_tiles=total,partitions=len(parts),
                packager='pmtiles-mosaic partition CLI; thirteen --from-source MBTiles inputs')
    (out/'partition-software.json').write_text(json.dumps(record,indent=2,sort_keys=True)+'\n')
    (out/'archives.json').write_text(json.dumps({**original,'partitions':parts,
                                               'verified_tiles':total,'packaging':record},indent=2,sort_keys=True)+'\n')
    print(json.dumps(dict(event='partitions_verified',tiles=total,archives=len(parts))),flush=True)


if __name__=='__main__':main()
