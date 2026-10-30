"""Measure native mosaic tile reads with an RSS guard; keep portable trial logs."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import time
from unittest.mock import patch

import rasterio

from build import atomic_json,sha256
from production_tiles import Sampler,validate_rgba,webp_bytes
from terrain_rgb import encode_terrain_rgb
from river_flattening.whitebox import measured


def rss():
    fields=dict(line.split(':',1) for line in Path('/proc/self/status').read_text().splitlines() if ':' in line)
    return round(int(fields['VmRSS'].split()[0])/1024,2)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workdir',type=Path,required=True)
    parser.add_argument('--label',required=True)
    parser.add_argument('--pool-size',type=int,default=48)
    parser.add_argument('--zoom',type=int,choices=range(5,13),default=5)
    parser.add_argument('--worker',action='store_true')
    args=parser.parse_args()
    if Path(args.label).name!=args.label:parser.error('Use a plain label')
    folder=args.workdir/'reports'/('tile-memory-'+args.label)
    folder.mkdir(parents=True,exist_ok=True)
    if not args.worker:
        command=[sys.executable,str(Path(__file__).resolve()),*sys.argv[1:],'--worker']
        measured(command,folder/'worker.log',1200)
        print((folder/'worker.json').read_text())
        return
    original_env=rasterio.Env
    def environment(**options):
        options.update(GDAL_MAX_DATASET_POOL_SIZE=args.pool_size,VRT_NUM_THREADS='1')
        return original_env(**options)
    records=[];started=time.monotonic()
    def event(phase,x,y,**extra):
        print(json.dumps(dict(phase=phase,x=x,y=y,rss_mib=rss(),elapsed_seconds=round(time.monotonic()-started,2),**extra)),flush=True)
    with patch('production_tiles.rasterio.Env',side_effect=environment):
        with Sampler(args.workdir/'mosaics/final.vrt',args.workdir/'mosaics/final-validity.vrt',args.zoom) as sampler:
            g=sampler.grid
            for x in range(g['xmin'],g['xmax']+1):
                for y in range(g['ymin'],g['ymax']+1):
                    event('read_begin',x,y)
                    height,valid=sampler.read(x,y)
                    event('read_done',x,y,valid_pixels=int(valid.sum()))
                    rgba=encode_terrain_rgb(height,valid=valid,round_bits=16-args.zoom)
                    validate_rgba(rgba,args.zoom,valid);data=webp_bytes(rgba)
                    records.append(dict(x=x,y=y,valid=int(valid.sum()),rgba_sha256=hashlib.sha256(rgba.tobytes()).hexdigest(),webp_sha256=hashlib.sha256(data).hexdigest()))
                    event('tile_done',x,y)
    atomic_json(folder/'tiles.json',dict(zoom=args.zoom,pool_size=args.pool_size,tiles=records,
                sampler_sha256=sha256(Path(__file__).with_name('production_tiles.py')),
                elapsed_seconds=round(time.monotonic()-started,2)))


if __name__=='__main__':main()
