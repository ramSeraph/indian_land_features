"""Native water contacts and applied values for complete connected river fits."""
from pathlib import Path
import numpy as np
import rasterio
from build import read_json


def adjacent_contacts(pixels,owners,heights,width,sheets=None):
    order=np.argsort(pixels);pixels=pixels[order];owners=owners[order];heights=heights[order]
    sheets=sheets[order] if sheets is not None else None
    pairs=[];values=[]
    for offset in (1,width):
        j=np.searchsorted(pixels,pixels+offset);safe=np.minimum(j,max(0,len(j)-1))
        if not len(j): continue
        good=(j<len(j))&(pixels[safe]==pixels+offset)&(owners!=owners[safe])
        if offset==1: good&=pixels%width<width-1
        good&=np.isfinite(heights)&np.isfinite(heights[safe])&(heights!=-255)&(heights[safe]!=-255)
        if sheets is not None: good&=sheets!=sheets[safe]
        i=np.flatnonzero(good);j=j[good]
        pairs.append(np.column_stack((owners[i],owners[j])))
        values.append(np.column_stack((heights[i],heights[j])))
    if not pairs: return dict(owners=np.empty((0,2),int),heights=np.empty((0,2),'float32'))
    return dict(owners=np.concatenate(pairs),heights=np.concatenate(values))


def contact_files(root,work):
    from build import atomic_npz
    records=read_json(work/'geometry.json')['sheets'];grid=read_json(work/'grid.json')
    t=rasterio.Affine(*grid['transform']);width=grid['width'];paths=[];border=[]
    folder=work/'contacts';folder.mkdir(exist_ok=True)
    for k,record in enumerate(records):
        name=record['sheet']
        with np.load(work/'geometry'/(Path(name).stem+'.npz')) as data:
            pixels,owners=data['pixels'],data['owners']
        with rasterio.open(root/'sheets/coastal'/name) as src:
            x,y=(~t)*(src.transform.c,src.transform.f);x,y=round(x),round(y)
            row=pixels//width-y;col=pixels%width-x
            heights=src.read(1)[row,col]
            edge=(row==0)|(row==src.height-1)|(col==0)|(col==src.width-1)
        path=folder/(Path(name).stem+'.npz')
        atomic_npz(path,**adjacent_contacts(pixels,owners,heights,width));paths.append(path)
        border.append((pixels[edge],owners[edge],heights[edge],np.full(edge.sum(),k,int)))
    if border:
        pixels,owners,heights,sheets=(np.concatenate([b[i] for b in border]) for i in range(4))
        path=folder/'sheet-borders.npz'
        atomic_npz(path,**adjacent_contacts(pixels,owners,heights,width,sheets));paths.append(path)
    return paths


def native_values(network):
    left,right=network['left'],network['right'];t=network['fraction']
    good=(left>=0)&(right>=0)&network['native_supported']
    values=np.full(len(left),np.nan,dtype=np.float32)
    values[good]=(1-t[good])*network['fitted'][left[good]]+t[good]*network['fitted'][right[good]]
    groups=np.full(len(left),-1,int)
    groups[good]=network['component'][left[good]]
    return values,groups


def read_contacts(table):
    if isinstance(table,(str,Path)):
        with np.load(table) as data:
            return data['owners'],data['heights']
    return table['owners'],table['heights']
