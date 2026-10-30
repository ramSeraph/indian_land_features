"""Add lower zooms to a completed mosaic while preserving existing archives.

python extend_lowzooms.py --workdir /path/to/build --min-zoom 0 --detach
"""
import argparse
from contextlib import ExitStack
from datetime import datetime,timezone
import fcntl
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import traceback

from build import atomic_json,emit,read_json,sha256,write_run_notes
from production_tiles import QUANTIZATION,tiles,package


def verify_archives(root,prefix,manifest):
    for entry in manifest['partitions']:
        path=root/'output'/prefix/entry['file']
        if path.stat().st_size!=entry['bytes'] or sha256(path)!=entry['sha256']:
            raise ValueError(f'Archive changed: {entry["file"]}')


def extend(root,prefix='final',min_zoom=0):
    root=Path(root)
    manifest_path=root/f'output/{prefix}/archives.json'
    current=read_json(manifest_path)
    if min(p['minzoom'] for p in current['partitions'])<=min_zoom:
        verify_archives(root,prefix,current)
        return current
    baseline_path=root/f'reports/{prefix}-lowzoom-base.json'
    if baseline_path.exists():
        baseline=read_json(baseline_path)
    else:
        baseline=dict(archives=current,tiles=read_json(root/f'reports/{prefix}-tiles.json'),
                      package=read_json(root/f'reports/{prefix}-package.json'))
        if baseline['tiles']['recipe']!=current['recipe']:
            raise ValueError('Original tile and archive recipes differ')
        atomic_json(baseline_path,baseline)
        notes=root/'reports/RUN_NOTES.md'
        if notes.exists():shutil.copy2(notes,root/f'reports/{prefix}-before-lowzooms.md')
    original=baseline['archives'];old_recipe=original['recipe']
    first=old_recipe['minzoom'];last=old_recipe['maxzoom']
    if not 0<=min_zoom<first<=5 or last>12:
        raise ValueError('This extension adds lower zooms below an existing z5-or-lower build')
    if old_recipe['round_bits'] not in ('16-z',QUANTIZATION):
        raise ValueError('Unrecognized existing quantization schedule')
    if [z['zoom'] for z in baseline['tiles']['zooms']]!=list(range(first,last+1)):
        raise ValueError('Original zoom coverage is incomplete')
    emit('lowzoom_verify_existing',archives=len(original['partitions']))
    verify_archives(root,prefix,original)
    label=f'{prefix}-z{min_zoom}-{first-1}'
    emit('lowzoom_generate',label=label,minzoom=min_zoom,maxzoom=first-1)
    tiles(root,prefix,min_zoom,first-1,bounds=old_recipe['bounds'],
          output_prefix=label,expected_inputs=old_recipe)
    package(root,label)
    added=read_json(root/f'output/{label}/archives.json')
    added_tiles=read_json(root/f'reports/{label}-tiles.json')
    # Added archives have unique names; old files and their cached byte ranges
    # remain untouched. The manifest is the final atomic publication step.
    partitions=[]
    for index,entry in enumerate(added['partitions']):
        name=f'CartoDEM_v3r1_EGM2008_z{min_zoom}_{first-1}_{index:03d}.pmtiles'
        source=root/f'output/{label}'/entry['file']
        target=root/f'output/{prefix}'/name
        if target.exists():
            if sha256(target)!=entry['sha256']:raise ValueError('Existing extension archive differs')
        else:os.link(source,target)
        partitions.append({**entry,'file':name})
    partitions+=original['partitions']
    recipe={**added['recipe'],'type':'zoom-extension','minzoom':min_zoom,'maxzoom':last,
            'base_recipe':old_recipe,'added_recipe':added['recipe'],
            'extension_code_sha256':sha256(Path(__file__)),
            'reused_zooms':list(range(first,last+1)),
            'base_record':str(baseline_path.relative_to(root))}
    zooms=[{**z,'database':f'tiles/{label}/z{z["zoom"]}.sqlite'} for z in added_tiles['zooms']]
    zooms+=baseline['tiles']['zooms']
    combined=dict(metadata=added['metadata'],partitions=partitions,recipe=recipe,
                  verified_tiles=sum(p['verified_tiles'] for p in partitions))
    if sum(z['tiles'] for z in zooms)!=combined['verified_tiles']:
        raise ValueError('Combined package and zoom counts differ')
    atomic_json(root/f'reports/{prefix}-tiles.json',dict(recipe=recipe,zooms=zooms))
    atomic_json(root/f'reports/{prefix}-package.json',dict(partitions=partitions,
                verified_tiles=combined['verified_tiles'],total_bytes=sum(p['bytes'] for p in partitions)))
    atomic_json(manifest_path,combined)
    emit('lowzoom_complete',added_tiles=sum(z['tiles'] for z in added_tiles['zooms']),
         total_tiles=combined['verified_tiles'],archives=len(partitions))
    return combined


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--prefix',default='final')
    parser.add_argument('--min-zoom',type=int,default=0)
    parser.add_argument('--detach',action='store_true')
    args=parser.parse_args();root=args.workdir.expanduser().resolve()
    if Path(args.prefix).name!=args.prefix or args.prefix in ('.','..'):
        parser.error('Use a plain output prefix')
    if not 0<=args.min_zoom<5:parser.error('The lower zoom must be between 0 and 4')
    code=Path(__file__).resolve().parent
    if args.detach:
        command=[sys.executable,str(Path(__file__).resolve()),*[v for v in sys.argv[1:] if v!='--detach']]
        with (root/f'logs/lowzoom-{args.prefix}.log').open('a') as log:
            child=subprocess.Popen(command,cwd=code,stdin=subprocess.DEVNULL,stdout=log,stderr=log,start_new_session=True)
        print(json.dumps(dict(pid=child.pid,status=f'reports/lowzoom-{args.prefix}.json')));return
    record=dict(state='running',pid=os.getpid(),started=datetime.now(timezone.utc).isoformat(),
                minzoom=args.min_zoom,source_sha256={name:sha256(code/name) for name in
                ('extend_lowzooms.py','production_tiles.py','terrain_rgb.py','build.py')})
    status=root/f'reports/lowzoom-{args.prefix}.json'
    with ExitStack() as stack:
        for name in (f'run-{args.prefix}.lock','build.lock'):
            lock=stack.enter_context((root/'logs'/name).open('a'))
            fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        atomic_json(status,record)
        try:
            result=extend(root,args.prefix,args.min_zoom)
            write_run_notes(root,args.prefix)
            record.update(state='complete',tiles=result['verified_tiles'],archives=len(result['partitions']),
                          total_bytes=sum(p['bytes'] for p in result['partitions']))
        except BaseException as error:
            record.update(state='failed',error=str(error));traceback.print_exc()
            raise
        finally:
            record['peak_rss_kib']=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/(1024 if sys.platform=='darwin' else 1))
            record['finished']=datetime.now(timezone.utc).isoformat();atomic_json(status,record)


if __name__=='__main__':main()
