"""Make a scientific comparison figure for an existing river pilot."""
import argparse
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
import rasterio
from matplotlib.colors import LightSource


def plot(folder):
    folder=Path(folder)
    with rasterio.open(folder/'coastal_baseline.tif') as src:
        before=src.read(1); extent=[src.bounds.left,src.bounds.right,src.bounds.bottom,src.bounds.top]
    with rasterio.open(folder/'river_pilot.tif') as src:
        after=src.read(1)
    with rasterio.open(folder/'classes.tif') as src:
        classes=src.read(1)
    valid=before!=-255
    river=(classes==3)&valid
    land=np.where(valid,before,np.nan)
    shade=LightSource(315,45).hillshade(np.nan_to_num(land),vert_exag=1,dx=30,dy=30)
    fig,axes=plt.subplots(1,3,figsize=(17,8),layout='constrained',sharex=True,sharey=True)
    images=[]
    for ax,values,title in zip(axes,[before,after,after-before],
                               ['Corrected CartoDEM river heights','Candidate river surface','Change inside river mask']):
        ax.imshow(np.where(valid,shade,np.nan),extent=extent,cmap='gray',vmin=0,vmax=1,alpha=.45)
        kwargs=dict(cmap='viridis',vmin=-5,vmax=25) if len(images)<2 else dict(cmap='RdBu_r',vmin=-10,vmax=10)
        images.append(ax.imshow(np.where(river,values,np.nan),extent=extent,**kwargs))
        ax.set_title(title)
        ax.set_xlabel('Longitude')
        ax.ticklabel_format(useOffset=False)
    axes[0].set_ylabel('Latitude')
    fig.colorbar(images[0],ax=axes[:2],label='EGM2008 elevation (m)',shrink=.6)
    fig.colorbar(images[2],ax=axes[2],label='Candidate minus baseline (m)',shrink=.6)
    report=json.loads((folder/'report.json').read_text())
    fig.suptitle(f"River-method pilot: {report['treated_river_pixels']:,} of {report['valid_river_pixels']:,} valid river pixels treated\n"
                 'Land, lake class 2 and original nodata are unchanged; unsupported components retain the baseline',fontsize=12)
    path=folder/'comparison.png'
    fig.savefig(path,dpi=150)
    plt.close(fig)
    return path


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('folder',type=Path)
    print(plot(parser.parse_args().folder))
