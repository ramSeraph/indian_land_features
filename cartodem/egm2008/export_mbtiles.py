#!/usr/bin/env python3
"""Export completed per-zoom tile databases to standard MBTiles, without encoding.

Only coordinates and schema change: XYZ y becomes the MBTiles TMS tile_row.
Every stored image is verified against its original recorded SHA-256.
"""
import argparse
from contextlib import closing
import hashlib
import json
import math
from pathlib import Path
import sqlite3


def file_sha(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def bounds(z, xmin, xmax, ymin, ymax):
    n = 2**z
    lat = lambda y: math.degrees(math.atan(math.sinh(math.pi*(1-2*y/n))))
    return [xmin/n*360-180, lat(ymax+1), (xmax+1)/n*360-180, lat(ymin)]


def export_zoom(source, target, zoom, metadata, expected_count):
    with closing(sqlite3.connect(source.resolve().as_uri()+'?mode=ro', uri=True)) as src:
        expected = {(x, 2**zoom-1-y): digest for x,y,digest in src.execute('SELECT x,y,sha256 FROM tiles')}
        if len(expected) != expected_count:
            raise ValueError('Source tile count differs')
        extent = src.execute('SELECT MIN(x),MAX(x),MIN(y),MAX(y) FROM tiles').fetchone()
        box = bounds(zoom, *extent)
        info = {**metadata, 'minzoom': str(zoom), 'maxzoom': str(zoom),
                'bounds': ','.join(map(str,box)), 'center': f'{(box[0]+box[2])/2},{(box[1]+box[3])/2},{zoom}',
                'type': 'baselayer', 'version': '1', 'scheme': 'tms',
                'description': 'CartoDEM EGM2008; Mapbox TerrainRGB; 512 px lossless WebP; opaque 0 m for nodata.'}
        if not target.exists():
            partial = target.with_suffix('.mbtiles.part')
            if partial.exists():
                raise ValueError('Incomplete export exists: '+str(partial))
            with closing(sqlite3.connect(partial)) as dst:
                dst.execute('PRAGMA journal_mode=DELETE')
                dst.execute('PRAGMA synchronous=NORMAL')
                dst.execute('CREATE TABLE metadata (name TEXT PRIMARY KEY, value TEXT)')
                dst.execute('CREATE TABLE tiles (zoom_level INTEGER NOT NULL, tile_column INTEGER NOT NULL, '
                            'tile_row INTEGER NOT NULL, tile_data BLOB NOT NULL, '
                            'PRIMARY KEY (zoom_level,tile_column,tile_row))')
                dst.executemany('INSERT INTO metadata VALUES (?,?)', [
                    (k,v if isinstance(v,str) else json.dumps(v)) for k,v in info.items()])
                for x,y,data,digest in src.execute('SELECT x,y,data,sha256 FROM tiles ORDER BY tileid'):
                    if hashlib.sha256(data).hexdigest() != digest:
                        raise ValueError(f'Source tile changed: {zoom}/{x}/{y}')
                    dst.execute('INSERT INTO tiles VALUES (?,?,?,?)', (zoom,x,2**zoom-1-y,data))
                dst.commit()
            candidate = partial
        else:
            candidate = target
        # Check coordinates and bytes through the actual MBTiles schema.
        with closing(sqlite3.connect(candidate.resolve().as_uri()+'?mode=ro', uri=True)) as db:
            for z,x,row,data in db.execute('SELECT zoom_level,tile_column,tile_row,tile_data FROM tiles'):
                if z != zoom or expected.pop((x,row), None) != hashlib.sha256(data).hexdigest():
                    raise ValueError('MBTiles coordinates or image bytes changed')
        if expected:
            raise ValueError('MBTiles export is missing tiles')
        if candidate != target:
            candidate.replace(target)
    return dict(file=target.name, zoom=zoom, tiles=expected_count,
                bytes=target.stat().st_size, sha256=file_sha(target))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', required=True, type=Path)
    parser.add_argument('--label', default='display')
    parser.add_argument('--outdir', required=True, type=Path)
    args = parser.parse_args()
    root = args.workdir.resolve()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    report = json.loads((root/f'reports/{args.label}-tiles.json').read_text())
    metadata = json.loads((root/f'output/{args.label}/archives.json').read_text())['metadata']
    results = []
    for info in report['zooms']:
        z = info['zoom']
        result = export_zoom(root/f'tiles/{args.label}/z{z}.sqlite', out/f'z{z}.mbtiles',
                             z, metadata, info['tiles'])
        results.append(result)
        print(json.dumps(dict(event='verified_mbtiles', **result)), flush=True)
    (out/'exports.json').write_text(json.dumps(dict(sources=results, tiles=sum(r['tiles'] for r in results)),indent=2)+'\n')


if __name__ == '__main__':
    main()
