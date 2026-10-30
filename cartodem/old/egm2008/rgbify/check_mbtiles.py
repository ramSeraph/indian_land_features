"""Check rio-rgbify's output coverage, including silently skipped worker errors."""
import json
import math
from pathlib import Path
import sqlite3
import sys


def check(path):
    with sqlite3.connect(Path(path).resolve().as_uri()+'?mode=ro', uri=True) as db:
        meta=dict(db.execute('SELECT name,value FROM metadata'))
        if meta.get('format')!='webp': raise ValueError('Expected WebP tiles')
        raw=meta['bounds']
        west,south,east,north=json.loads(raw) if raw.startswith('[') else map(float,raw.split(','))
        minz,maxz=int(meta['minzoom']),int(meta['maxzoom'])
        count=0
        for z in range(minz,maxz+1):
            n=2**z
            def x(lon): return max(0,min(n-1,math.floor((lon+180)/360*n)))
            def y(lat):
                lat=max(-85.0511287798066,min(85.0511287798066,lat))
                return max(0,min(n-1,math.floor((1-math.asinh(math.tan(math.radians(lat)))/math.pi)/2*n)))
            xmin,xmax=x(west+1e-10),x(east-1e-10)
            ymin,ymax=y(north-1e-10),y(south+1e-10)
            expected=(xmax-xmin+1)*(ymax-ymin+1)
            actual=db.execute('SELECT COUNT(*),MIN(tile_column),MAX(tile_column),MIN(tile_row),MAX(tile_row) '
                              'FROM tiles WHERE zoom_level=?',(z,)).fetchone()
            if actual!=(expected,xmin,xmax,n-1-ymax,n-1-ymin):
                raise ValueError(f'Incomplete or unexpected coverage at z{z}: {actual}, expected {expected} tiles')
            unique=db.execute('SELECT COUNT(*) FROM (SELECT tile_column,tile_row FROM tiles '
                              'WHERE zoom_level=? GROUP BY tile_column,tile_row)',(z,)).fetchone()[0]
            if unique!=expected: raise ValueError('Duplicate tile coordinates')
            count+=expected
        if db.execute('SELECT COUNT(*) FROM tiles').fetchone()[0]!=count:
            raise ValueError('Unexpected zoom levels')
    return dict(file=Path(path).name, tiles=count, minzoom=minz, maxzoom=maxz, format='webp')


if __name__=='__main__':
    print(json.dumps(check(sys.argv[1])))
