"""Run restartable conversion stages sequentially with a durable status record.

The coordinator and each heavy stage run in separate processes. --detach writes
logs and survives disconnection. All processing settings live in the versioned
modules; this runner only coordinates their existing commands.
"""
import argparse
from datetime import datetime,timezone
import fcntl
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys
import traceback

STAGES=('download','prepare','geometry','routing','fit','cogs')


def write_status(path,data):
    temporary=path.with_suffix('.json.part')
    temporary.write_text(json.dumps(data,indent=2)+'\n');temporary.replace(path)


def commands(args):
    common=['--workdir',str(args.workdir)]
    selection=['--sheets',*args.sheets] if args.sheets else []
    for stage in STAGES:
        if stage=='geometry':
            command=[sys.executable,'-m','river_flattening.prepare',*common,
                     '--label',args.label,*selection]
        elif stage in ('routing','fit','cogs'):
            command=[sys.executable,'-m','river_flattening.national',
                     'apply' if stage=='cogs' else stage,*common,'--label',args.label]
        else:
            command=[sys.executable,'build.py',stage,*common,*selection]
        yield stage,command


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--label',default='final')
    parser.add_argument('--sheets',nargs='+')
    parser.add_argument('--resume-from',choices=STAGES,default='download')
    parser.add_argument('--through',choices=STAGES,default='cogs',
                        help='Stop after COGs by default; use the documented display/MBTiles pipeline for tiles')
    parser.add_argument('--detach',action='store_true')
    args=parser.parse_args();args.workdir=args.workdir.expanduser().resolve()
    if Path(args.label).name!=args.label or args.label in ('.','..'):
        parser.error('Use a plain output label')
    first,last=STAGES.index(args.resume_from),STAGES.index(args.through)
    if first>last:parser.error('The start stage follows the last stage')
    if first>STAGES.index('geometry'):
        recipe=args.workdir/'rivers'/args.label/'geometry-recipe.json'
        if not recipe.exists() or json.loads(recipe.read_text()).get('river_source')!='glo30':
            parser.error('GLO-30 river source does not match cached geometry; rebuild geometry under a fresh label')
    for directory in ('logs','reports'):(args.workdir/directory).mkdir(parents=True,exist_ok=True)
    code=Path(__file__).resolve().parent
    log_path=args.workdir/'logs'/f'build-{args.label}.log'
    if args.detach:
        child=[sys.executable,str(Path(__file__).resolve()),
               *[v for v in sys.argv[1:] if v!='--detach']]
        with log_path.open('a') as log:
            process=subprocess.Popen(child,cwd=code,stdin=subprocess.DEVNULL,
                                     stdout=log,stderr=subprocess.STDOUT,start_new_session=True)
        print(json.dumps(dict(pid=process.pid,log=f'logs/{log_path.name}',
                              status=f'reports/build-{args.label}.json')))
        return
    status_path=args.workdir/'reports'/f'build-{args.label}.json'
    record=dict(state='running',pid=os.getpid(),label=args.label,
                started=datetime.now(timezone.utc).isoformat(),completed_stages=[],
                source_sha256={str(p.relative_to(code)):hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in sorted(code.rglob('*.py')) if '.venv' not in p.parts},
                river_source='glo30',
                scope='EGM2008, unified nodata, GLO-30 ocean and rivers, candidate COGs',
                publication_requires_gradient_review=True)
    with (args.workdir/'logs'/f'run-{args.label}.lock').open('a') as lock:
        fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
        try:
            for stage,command in list(commands(args))[first:last+1]:
                record.update(stage=stage);write_status(status_path,record)
                print(json.dumps(dict(event='build_stage',stage=stage)),flush=True)
                subprocess.run(command,cwd=code,check=True)
                record['completed_stages'].append(stage)
            record.update(state='complete',finished=datetime.now(timezone.utc).isoformat())
            write_status(status_path,record)
        except BaseException as exc:
            record.update(state='failed',error=type(exc).__name__,
                          exit_code=getattr(exc,'returncode',None),
                          finished=datetime.now(timezone.utc).isoformat())
            write_status(status_path,record);traceback.print_exc();sys.exit(1)


if __name__=='__main__':main()
