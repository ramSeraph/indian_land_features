#!/usr/bin/env python3
"""Upload verified COG or PMTiles assets to an existing GitHub release.

Requires gh authentication. Re-running skips assets with matching SHA-256 digests;
conflicting existing assets cause an error and are never overwritten.
"""
import argparse
import json
from pathlib import Path
import subprocess
from urllib.parse import quote

from build import read_json,sha256


def assets(root,label,kind):
    out=root/f'output/{label}'
    if kind=='cogs':
        record=read_json(root/f'reports/{label}-cog-inputs.json')
        directory=root/record['cog_directory'];entries=record['sheets']
        extras=('cogs.json','COG_SHA256SUMS.txt')
    else:
        directory=out;entries=read_json(out/'archives.json')['partitions']
        extras=('archives.json','PMTILES_SHA256SUMS.txt','TILE_BUILD_NOTES.md',
                'software.json','tiling-source.tar.gz')
    result=[]
    for entry in entries:
        if Path(entry['file']).name!=entry['file']:raise ValueError('Expected a plain asset filename')
        path=directory/entry['file']
        if path.stat().st_size!=entry['bytes'] or sha256(path)!=entry['sha256']:
            raise ValueError(f'Asset differs from the build manifest: {path.name}')
        result.append(dict(path=path,name=path.name,bytes=entry['bytes'],sha256=entry['sha256']))
    for name in extras:
        path=out/name
        result.append(dict(path=path,name=name,bytes=path.stat().st_size,sha256=sha256(path)))
    if len({e['name'] for e in result})!=len(result):raise ValueError('Duplicate asset name')
    if any(e['bytes']>=2**31 for e in result):raise ValueError('GitHub release assets must be under 2 GiB')
    return result


def gh_json(gh,*args):
    return json.loads(subprocess.check_output([gh,*args],text=True))


def upload(gh,repo,tag,files):
    release=gh_json(gh,'api',f'repos/{repo}/releases/tags/{quote(tag,safe="")}')
    endpoint=f'repos/{repo}/releases/{release["id"]}/assets?per_page=100'
    def existing():
        pages=gh_json(gh,'api','--paginate','--slurp',endpoint)
        return {a['name']:a for page in pages for a in page}
    remote=existing()
    new=[f for f in files if f['name'] not in remote]
    if len(remote)+len(new)>1000:raise ValueError('GitHub permits up to 1000 assets per release')
    # Check every collision before uploading anything, including on a resumed run.
    for entry in files:
        current=remote.get(entry['name'])
        if current and (current.get('digest')!='sha256:'+entry['sha256'] or current['size']!=entry['bytes']):
            raise ValueError(f'Existing asset differs or has no verifiable digest: {entry["name"]}')
    for entry in files:
        if entry['name'] in remote:
            print(f'Already verified: {entry["name"]}',flush=True);continue
        subprocess.run([gh,'release','upload',tag,str(entry['path']),'--repo',repo],check=True)
        current=existing().get(entry['name'])
        if not current or current.get('digest')!='sha256:'+entry['sha256'] or current['size']!=entry['bytes']:
            raise ValueError(f'Uploaded asset verification failed: {entry["name"]}')
        print(f'Uploaded and verified: {entry["name"]}',flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--label',default='display')
    parser.add_argument('--kind',choices=('cogs','pmtiles'),required=True)
    parser.add_argument('--repo',required=True,help='GitHub OWNER/REPOSITORY')
    parser.add_argument('--tag',required=True,help='Existing release tag')
    parser.add_argument('--gh',default='gh',help='GitHub CLI executable')
    parser.add_argument('--dry-run',action='store_true',help='Verify local files and print the plan; no network calls')
    args=parser.parse_args()
    if not args.label or Path(args.label).name!=args.label or args.label in ('.','..'):
        parser.error('Use a plain output label')
    if args.repo.count('/')!=1 or any(not p or p in ('.','..') for p in args.repo.split('/')):
        parser.error('Repository must be OWNER/REPOSITORY')
    files=assets(args.workdir.resolve(),args.label,args.kind)
    print(json.dumps(dict(repository=args.repo,tag=args.tag,assets=len(files),bytes=sum(f['bytes'] for f in files))))
    if args.dry_run:
        for entry in files:print(json.dumps({k:v for k,v in entry.items() if k!='path'}))
    else:upload(args.gh,args.repo,args.tag,files)


if __name__=='__main__':main()
