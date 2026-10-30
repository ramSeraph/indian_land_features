#!/usr/bin/env python3
"""Prepare and explicitly replace one existing release's CartoDEM assets.

The release object and its notes are never edited. Preparation verifies source
checksums and preserves the existing filenames. Publication is resumable and
requires --publish; it deletes conflicting assets only from the saved inventory.
"""
import argparse
from contextlib import closing
from datetime import datetime, timezone
import fcntl
import json
import math
import os
from pathlib import Path
import re
import sqlite3
import subprocess
import sys
import tarfile
import time
import traceback
from urllib.parse import quote

from build import atomic_json, read_json, sha256
from upload_release import assets


def emit(event, **fields):
    print(json.dumps(dict(event=event, **fields)), flush=True)


class GitHub:
    def __init__(self, executable, repo):
        self.executable = executable
        self.base = 'repos/' + repo
        self.last_mutation = 0

    def request(self, endpoint, method='GET', fields=None, path=None, mime=None):
        for attempt in range(7):
            if method != 'GET':
                time.sleep(max(0, 1.1 - (time.monotonic() - self.last_mutation)))
            command = [self.executable, 'api', '--include', '--method', method, endpoint]
            if path is not None:
                command += ['--input', str(path), '-H', 'Content-Type: ' + mime]
            for key, value in (fields or {}).items():
                command += ['-f', key + '=' + str(value)]
            result = subprocess.run(command, capture_output=True, text=True)
            if method != 'GET':
                self.last_mutation = time.monotonic()
            raw = result.stdout.replace('\r\n', '\n')
            head, _, body = raw.partition('\n\n')
            headers = dict(line.lower().split(': ', 1) for line in head.splitlines()[1:] if ': ' in line)
            if result.returncode == 0:
                return json.loads(body) if body.strip() else None
            message = (body or result.stderr).strip()
            match = re.search(r'HTTP/[\d.]+ (\d+)', head)
            status = int(match.group(1)) if match else None
            limited = status == 429 or (status == 403 and 'rate limit' in message.lower())
            if attempt == 6 or not (limited or status in (500, 502, 503, 504) or status is None):
                raise RuntimeError(f'{method} {endpoint}: {message}')
            # Upload retries must reconcile the name first: a lost response can
            # follow a successful upload. Let the caller perform that recovery.
            if path is not None and not limited:
                raise RuntimeError(f'Upload response failed; rerun to reconcile assets: {message}')
            delay = max(60 * 2**attempt if limited else 5 * 2**attempt,
                        float(headers.get('retry-after', '0')))
            if headers.get('x-ratelimit-remaining') == '0':
                delay = max(delay, float(headers.get('x-ratelimit-reset', '0')) - time.time() + 2)
            emit('api_retry', method=method, status=status, wait_seconds=delay)
            time.sleep(delay)

    def release(self, tag):
        return self.request(self.base + '/releases/tags/' + quote(tag, safe=''))

    def inventory(self, release_id):
        found = []
        for page in range(1, 30):
            current = self.request(f'{self.base}/releases/{release_id}/assets?per_page=100&page={page}')
            found.extend(current)
            if len(current) < 100:
                return {entry['name']: entry for entry in found}
        raise ValueError('Unexpectedly large release')


def matches(remote, local):
    return (remote is not None and remote.get('state') == 'uploaded'
            and remote.get('size') == local['bytes']
            and remote.get('digest') == 'sha256:' + local['sha256'])


def validate_files(root, entries):
    for entry in entries:
        path = root / entry['path']
        if path.stat().st_size != entry['bytes'] or sha256(path) != entry['sha256']:
            raise ValueError('Prepared asset changed: ' + entry['name'])


def publication_order(entry):
    name = entry['name']
    # Keep the existing map intact while COGs upload. Replace its PMTiles and
    # then its manifests together at the end of the job.
    rank = (0 if name.endswith('.tif') else 1 if name.endswith('.pmtiles') else
            3 if name.endswith('.mosaic.json') else 2)
    return rank, name


def tile_box(z, x, y):
    n = 2**z
    lat = lambda v: math.degrees(math.atan(math.sinh(math.pi * (1 - 2*v/n))))
    return [x/n*360-180, lat(y+1), (x+1)/n*360-180, lat(y)]


def prepare(root, label, folder, gh, repo, tag, prefix, archive_dir=None):
    from pmtiles.reader import Reader, all_tiles
    state = read_json(root / f'reports/rebuild-{label}.json')
    if state['state'] != 'complete':
        raise ValueError('Build is not complete')
    # Fully validate both products before preparing any release mutation.
    cogs = assets(root, label, 'cogs')
    out = root / f'output/{label}'
    archive_dir = archive_dir or out
    manifest = read_json(archive_dir/'archives.json')
    if archive_dir == out:
        pmtiles = assets(root, label, 'pmtiles')
    else:
        pmtiles = [dict(path=archive_dir/e['file'], name=e['file'], bytes=e['bytes'], sha256=e['sha256'])
                   for e in manifest['partitions']]
        for name in ('software.json',):
            path = out/name
            pmtiles.append(dict(path=path, name=name, bytes=path.stat().st_size, sha256=sha256(path)))
        for name in ('requirements-mosaic.txt', 'partition-software.json'):
            path = archive_dir/name
            pmtiles.append(dict(path=path, name=name, bytes=path.stat().st_size, sha256=sha256(path)))
        validate_files(root, pmtiles)
    release_path = folder / 'release-before.json'
    if not release_path.exists():
        atomic_json(release_path, gh.release(tag))
    before = read_json(release_path)
    inventory_path = folder / 'assets-before.json'
    if not inventory_path.exists():
        atomic_json(inventory_path, list(gh.inventory(before['id']).values()))
    old = {a['name']: a for a in read_json(inventory_path)}
    if {n for n in old if n.endswith('.tif')} != {f['name'] for f in cogs if f['name'].endswith('.tif')}:
        raise ValueError('Existing and replacement COG filenames differ')
    generated = folder / 'prepared'
    generated.mkdir(exist_ok=True)
    new_manifest = json.loads(json.dumps(manifest))
    renamed = {}
    for i, part in enumerate(new_manifest['partitions']):
        renamed[part['file']] = f'{prefix}-part{i:04d}.pmtiles'
        part['file'] = renamed[part['file']]
    entries = []
    def add(path, name=None):
        entries.append(dict(path=os.path.relpath(path, root), name=name or path.name,
                            bytes=path.stat().st_size, sha256=sha256(path)))
    # Sources have already been hashed by assets(); reuse those digests.
    for asset in [*cogs, *pmtiles]:
        if asset['name'] in ('archives.json', 'PMTILES_SHA256SUMS.txt', 'TILE_BUILD_NOTES.md'):
            continue
        entries.append({**asset, 'path': os.path.relpath(asset['path'], root),
                        'name': renamed.get(asset['name'], asset['name'])})
    slices = {}
    details = {}
    seen = set()
    maxzoom = max(p['maxzoom'] for p in manifest['partitions'])
    topboxes = []
    samples = []
    for original, part in zip(manifest['partitions'], new_manifest['partitions']):
        with (archive_dir / original['file']).open('rb') as stream:
            def read_at(offset, length):
                stream.seek(offset)
                return stream.read(length)
            header = Reader(read_at).header()
            fields = ('min_lon_e7', 'min_lat_e7', 'max_lon_e7', 'max_lat_e7', 'min_zoom', 'max_zoom')
            h = {key: header[key] for key in fields}
            slices[part['file']] = {'header': h}
            sizes = {}
            first = last = None
            for (z, x, y), data in all_tiles(read_at):
                key = f'{z},{x},{y}'
                if key in seen:
                    raise ValueError('Duplicate tile across archives: ' + key)
                seen.add(key)
                sizes[key] = len(data)
                sample = dict(archive=part['file'], z=z, x=x, y=y)
                if first is None:
                    first = sample
                last = sample
                if z == maxzoom:
                    topboxes.append(tile_box(z, x, y))
            if len(sizes) != original['verified_tiles']:
                raise ValueError('Archive tile count changed')
            samples.extend([first, last])
            west, south, east, north = original['bounds']
            details[part['file']] = dict(tiles=sizes, min_zoom=h['min_zoom'], max_zoom=h['max_zoom'],
                                       bounds=[south, west, north, east])
    if len(seen) != state['tiles']:
        raise ValueError('Build tile count changed')
    bounds = [min(b[0] for b in topboxes), min(b[1] for b in topboxes),
              max(b[2] for b in topboxes), max(b[3] for b in topboxes)]
    header = dict(min_zoom=min(p['minzoom'] for p in manifest['partitions']), max_zoom=maxzoom,
                  min_lon_e7=round(bounds[0]*1e7), min_lat_e7=round(bounds[1]*1e7),
                  max_lon_e7=round(bounds[2]*1e7), max_lat_e7=round(bounds[3]*1e7),
                  center_zoom=8, center_lon_e7=round((bounds[0]+bounds[2])*5e6),
                  center_lat_e7=round((bounds[1]+bounds[3])*5e6), tile_type=4, tile_compression=1)
    mosaic = dict(version=1, metadata={**manifest['metadata'], 'type': 'baselayer', 'version': '1'},
                  header=header, slices=slices)
    tool_mosaic = archive_dir/(prefix+'.mosaic.json')
    if archive_dir != out and tool_mosaic.exists():
        mosaic = read_json(tool_mosaic)
        if set(mosaic['slices']) != set(slices):
            raise ValueError('pmtiles-mosaic manifest filenames differ from the verified archives')
    for name, content in ((prefix+'.mosaic.json', mosaic), (prefix+'.partition_info.json', details),
                          ('archives.json', new_manifest)):
        atomic_json(generated / name, content)
        add(generated / name)
    atomic_json(folder/'mosaic-samples.json', samples)
    base = f'https://github.com/{repo}/releases/download/{quote(tag, safe="")}/'
    import csv
    with (generated/'listing_files.csv').open('w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['name', 'size', 'url'])
        for entry in cogs:
            if entry['name'].endswith('.tif'):
                writer.writerow([entry['name'], entry['bytes'], base+quote(entry['name'])])
    add(generated/'listing_files.csv')
    (generated/'PMTILES_SHA256SUMS.txt').write_text(''.join(
        f"{p['sha256']}  {p['file']}\n" for p in new_manifest['partitions']))
    add(generated/'PMTILES_SHA256SUMS.txt')
    notes = (out/'TILE_BUILD_NOTES.md').read_text().split('| Archive |', 1)[0]
    if archive_dir != out:
        notes = notes.replace('PMTiles partitions target at most 1.8 billion payload bytes.',
                              'PMTiles partitions were produced by pmtiles-mosaic with a 2,000,000,000-byte limit.')
        notes += ('For publication, the encoded tile databases were exported to thirteen standard\n'
                  'MBTiles files. The pmtiles-mosaic 0.2.0 partition CLI read those files through\n'
                  'thirteen --from-source arguments and created geographic PMTiles partitions.\n'
                  'Every tile was verified byte for byte against the original encoded-tile\n'
                  'database. No raster resampling or image encoding was repeated.\n'
                  'partition-software.json and requirements-mosaic.txt record this environment.\n\n')
        # Preserve the actual numeric tiler. Update the runner/instructions to
        # select MBTiles export instead of its earlier packaging step.
        source_updates = ('rebuild_tiles.py', 'TILING_SETUP.md', 'install_tiles_arm64.sh',
                          'export_mbtiles.py', 'partition_mbtiles.sh', 'verify_partitioned.py',
                          'verify_mosaic.mjs', 'replace_release.py')
        with tarfile.open(out/'tiling-source.tar.gz', 'r:gz') as src, \
                tarfile.open(generated/'tiling-source.tar.gz', 'w:gz') as dst:
            for member in src.getmembers():
                if member.name not in (*source_updates, 'requirements-mosaic.txt'):
                    dst.addfile(member, src.extractfile(member) if member.isfile() else None)
            for name in source_updates:
                dst.add(Path(__file__).parent/name, arcname=name)
            dst.add(archive_dir/'requirements-mosaic.txt', arcname='requirements-mosaic.txt')
        add(generated/'tiling-source.tar.gz')
    notes += '| Archive | Bytes | SHA-256 |\n| --- | --- | --- |\n'
    for part in new_manifest['partitions']:
        notes += f"| {part['file']} | {part['bytes']} | `{part['sha256']}` |\n"
    (generated/'TILE_BUILD_NOTES.md').write_text(notes)
    add(generated/'TILE_BUILD_NOTES.md')
    if len(entries) > 1000 or any(e['bytes'] >= 2**31 for e in entries):
        raise ValueError('GitHub asset count or size limit exceeded')
    expected = {e['name'] for e in entries}
    if len(expected) != len(entries) or set(old) - expected:
        raise ValueError('Unexpected old asset or duplicate replacement filename')
    # Data first, then metadata; the live mosaic is replaced last.
    entries.sort(key=publication_order)
    plan = dict(repository=repo, tag=tag, release_id=before['id'], prefix=prefix,
                entries=entries, assets=len(entries), bytes=sum(e['bytes'] for e in entries),
                cogs=state['sheets'], pmtiles=len(renamed), tiles=len(seen),
                release_notes_sha256=__import__('hashlib').sha256((before['body'] or '').encode()).hexdigest())
    atomic_json(folder/'plan.json', plan)
    emit('prepared', assets=len(entries), cogs=state['sheets'], pmtiles=len(renamed), tiles=len(seen))
    return plan


def publish(root, folder, gh, plan):
    validate_files(root, plan['entries'])
    before = read_json(folder/'release-before.json')
    old = {a['name']: a for a in read_json(folder/'assets-before.json')}
    release = gh.release(plan['tag'])
    protected = ('id', 'tag_name', 'name', 'body', 'draft', 'prerelease', 'target_commitish')
    if any(release.get(k) != before.get(k) for k in protected):
        raise ValueError('Release metadata changed since preparation')
    remote = gh.inventory(plan['release_id'])
    wanted = {e['name'] for e in plan['entries']}
    if set(remote) - wanted:
        raise ValueError('Unplanned assets were added to the release')
    # Reject concurrent replacements before deleting anything.
    for entry in plan['entries']:
        current = remote.get(entry['name'])
        if current and not matches(current, entry):
            saved = old.get(entry['name'])
            if not saved or (current['id'], current['size'], current.get('digest')) != (
                    saved['id'], saved['size'], saved.get('digest')):
                raise ValueError('Asset changed outside this publication: ' + entry['name'])
    progress = dict(state='publishing', pid=os.getpid(), assets=plan['assets'], complete=0,
                    started=datetime.now(timezone.utc).isoformat(), bytes_complete=0)
    status = folder/'status.json'
    atomic_json(status, progress)
    try:
        for entry in plan['entries']:
            current = remote.get(entry['name'])
            if not matches(current, entry):
                if current:
                    gh.request(gh.base+'/releases/assets/'+str(current['id']), method='DELETE')
                    emit('removed_old_asset', name=entry['name'], asset_id=current['id'])
                mime = ('image/tiff' if entry['name'].endswith('.tif') else
                        'application/json' if entry['name'].endswith('.json') else
                        'text/csv' if entry['name'].endswith('.csv') else 'application/octet-stream')
                url = release['upload_url'].split('{')[0]+'?name='+quote(entry['name'], safe='')
                current = gh.request(url, method='POST', path=root/entry['path'], mime=mime)
                if not matches(current, entry):
                    raise ValueError('Uploaded digest or size differs: ' + entry['name'])
            progress['complete'] += 1
            progress['bytes_complete'] += entry['bytes']
            progress['last_asset'] = entry['name']
            progress['updated'] = datetime.now(timezone.utc).isoformat()
            atomic_json(status, progress)
            emit('verified_asset', name=entry['name'], complete=progress['complete'], total=plan['assets'])
        final = gh.inventory(plan['release_id'])
        if set(final) != wanted or any(not matches(final.get(e['name']), e) for e in plan['entries']):
            raise ValueError('Final release inventory differs from plan')
        release = gh.release(plan['tag'])
        if any(release.get(k) != before.get(k) for k in protected):
            raise ValueError('Release metadata changed during publication')
        atomic_json(folder/'assets-after.json', list(final.values()))
        atomic_json(folder/'release-after.json', release)
        progress.update(state='complete', release_notes_unchanged=True,
                        finished=datetime.now(timezone.utc).isoformat())
        emit('publication_complete', assets=plan['assets'], bytes=plan['bytes'])
    except BaseException as error:
        progress.update(state='failed', error=str(error))
        raise
    finally:
        atomic_json(status, progress)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--workdir', required=True, type=Path)
    p.add_argument('--label', default='display')
    p.add_argument('--repo', required=True)
    p.add_argument('--tag', required=True)
    p.add_argument('--prefix', required=True)
    p.add_argument('--gh', default='gh')
    p.add_argument('--pmtiles-dir', type=Path, help='Optional verified geographic repacking output')
    p.add_argument('--publish', action='store_true')
    p.add_argument('--detach', action='store_true')
    args = p.parse_args()
    if any(not re.fullmatch(r'[A-Za-z0-9_.-]+', s) or s in ('.', '..') for s in (args.tag, args.prefix, args.label)):
        p.error('Tag, prefix and label must be plain names')
    if not re.fullmatch(r'[A-Za-z0-9_.-]+/[A-Za-z0-9_.-]+', args.repo):
        p.error('Expected OWNER/REPOSITORY')
    root = args.workdir.resolve()
    folder = root/'publication'/args.tag
    folder.mkdir(parents=True, exist_ok=True)
    if args.detach:
        if not args.publish:
            p.error('--detach requires --publish')
        command = [sys.executable, str(Path(__file__).resolve()), *[a for a in sys.argv[1:] if a != '--detach']]
        with (folder/'publication.log').open('a') as f:
            child = subprocess.Popen(command, stdin=subprocess.DEVNULL, stdout=f, stderr=f, start_new_session=True)
        emit('started', pid=child.pid)
        return
    with (folder/'publication.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        gh = GitHub(args.gh, args.repo)
        if args.publish:
            plan = read_json(folder/'plan.json')
            if (plan['repository'], plan['tag'], plan['prefix']) != (args.repo, args.tag, args.prefix):
                raise ValueError('Publication target differs from prepared plan')
            publish(root, folder, gh, plan)
        else:
            prepare(root, args.label, folder, gh, args.repo, args.tag, args.prefix,
                    args.pmtiles_dir.resolve() if args.pmtiles_dir else None)


if __name__ == '__main__':
    main()
