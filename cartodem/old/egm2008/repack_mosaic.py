#!/usr/bin/env python3
"""Repartition encoded tiles geographically using pmtiles-mosaic==0.2.0.

Run in an isolated environment with requirements-mosaic.txt. No raster sampling
or image encoding occurs. Every output tile is checked against the build DB.
"""
import argparse
from contextlib import closing
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
from pathlib import Path
import shutil
import sqlite3
import tempfile

from pmtiles.reader import Reader, all_tiles
from pmtiles.tile import zxy_to_tileid
from pmtiles_mosaic.logger import get_logger
from pmtiles_mosaic.partition import Partitioner
from pmtiles_mosaic.tile_sources import create_source_from_paths


def digest(path):
    with path.open('rb') as f:
        return hashlib.file_digest(f, 'sha256').hexdigest()


def write_json(path, data):
    path.write_text(json.dumps(data, indent=2, sort_keys=True)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir', type=Path, required=True)
    parser.add_argument('--label', default='display')
    parser.add_argument('--outdir', type=Path, required=True)
    parser.add_argument('--prefix', required=True)
    parser.add_argument('--size-limit', type=int, default=2_000_000_000)
    args = parser.parse_args()
    root = args.workdir.resolve()
    out = args.outdir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    if list(out.glob('*.pmtiles')):
        raise ValueError('Repacking requires a fresh output directory')
    source = root/f'output/{args.label}'
    original = json.loads((source/'archives.json').read_text())
    started = datetime.now(timezone.utc).isoformat()
    paths = []
    for entry in original['partitions']:
        path = source/entry['file']
        if path.stat().st_size != entry['bytes'] or digest(path) != entry['sha256']:
            raise ValueError('Source archive changed: '+entry['file'])
        paths.append(str(path))
    logger = get_logger('publication-repack', 'INFO')
    temp = out/'tmp'
    temp.mkdir(exist_ok=True)
    tempfile.tempdir = str(temp)
    reader = create_source_from_paths(paths, logger=logger)
    try:
        partitioner = Partitioner(reader, str(out/args.prefix), args.size_limit,
                                  should_cache=True, logger=logger, exclude_transparent=False)
        # The stacked source retains only standard keys by default. Preserve
        # the DEM datum, encoding, quantization and nodata metadata explicitly.
        partitioner.src_metadata.update(original['metadata'])
        partitioner.partition()
        partitioner.finalize()
    finally:
        reader.cleanup()
    expected = {}
    report = json.loads((root/f'reports/{args.label}-tiles.json').read_text())
    for info in report['zooms']:
        path = root/f"tiles/{args.label}/z{info['zoom']}.sqlite"
        with closing(sqlite3.connect(path.as_uri()+'?mode=ro', uri=True)) as db:
            expected.update(db.execute('SELECT tileid,sha256 FROM tiles'))
    parts = []
    verified = 0
    for path in sorted(out.glob('*.pmtiles')):
        count = 0
        with path.open('rb') as stream:
            def read_at(offset, length):
                stream.seek(offset)
                return stream.read(length)
            h = Reader(read_at).header()
            for (z, x, y), data in all_tiles(read_at):
                key = zxy_to_tileid(z, x, y)
                if expected.pop(key, None) != hashlib.sha256(data).hexdigest():
                    raise ValueError(f'Changed, duplicate or unexpected tile: {z}/{x}/{y}')
                count += 1
        if path.stat().st_size >= 2**31:
            raise ValueError('Archive exceeds GitHub size limit')
        parts.append(dict(file=path.name, bytes=path.stat().st_size, sha256=digest(path),
                          tiles=count, verified_tiles=count, minzoom=h['min_zoom'], maxzoom=h['max_zoom'],
                          bounds=[h[k]/1e7 for k in ('min_lon_e7','min_lat_e7','max_lon_e7','max_lat_e7')]))
        verified += count
        print(json.dumps(dict(event='verified_repacked_archive', **parts[-1])), flush=True)
    if expected or verified != original['verified_tiles']:
        raise ValueError('Repacked tile coverage differs')
    software = {d.metadata['Name']: d.version for d in importlib.metadata.distributions()}
    (out/'requirements-mosaic.txt').write_text(''.join(f'{k}=={v}\n' for k,v in sorted(software.items())))
    record = dict(started=started, finished=datetime.now(timezone.utc).isoformat(),
                  software=software, source_manifest_sha256=digest(source/'archives.json'),
                  script_sha256=digest(Path(__file__)), size_limit_bytes=args.size_limit,
                  verified_tiles=verified, partitions=len(parts))
    write_json(out/'repack-software.json', record)
    write_json(out/'archives.json', {**original, 'partitions': parts, 'repack': record})
    print(json.dumps(dict(event='repack_complete', tiles=verified, archives=len(parts))), flush=True)


if __name__ == '__main__':
    main()
