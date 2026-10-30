"""Plot a retained native-resolution pilot crop, without downloading full sheets."""
import argparse
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource


def plot(path):
    path=Path(path)
    with np.load(path) as data:
        before,after,classes,t=(data[k] for k in ('before','after','classes','transform'))
    valid=before!=-255;river=valid&(classes==3);change=after-before
    extent=[t[2],t[2]+t[0]*before.shape[1],t[5]+t[4]*before.shape[0],t[5]]
    shade=LightSource(315,45).hillshade(np.where(valid,before,0),dx=30,dy=30)
    fig,axes=plt.subplots(1,3,figsize=(16,8),layout='constrained',sharex=True,sharey=True)
    for i,(ax,values,title) in enumerate(zip(axes,[before,after,change],
            ['EGM2008 coastal baseline','Integrated river pilot','Change; grey rivers unchanged'])):
        ax.imshow(np.where(valid,shade,np.nan),extent=extent,cmap='gray',alpha=.4)
        if i==2:
            ax.imshow(np.where(river,1,np.nan),extent=extent,cmap='gray',vmin=0,vmax=2)
            im=ax.imshow(np.where(river&(change!=0),values,np.nan),extent=extent,cmap='RdBu_r',vmin=-10,vmax=10)
        else:
            im=ax.imshow(np.where(river,values,np.nan),extent=extent,cmap='viridis',vmin=-5,vmax=25)
        ax.set_title(title);ax.set_xlabel('Longitude');ax.ticklabel_format(useOffset=False)
        fig.colorbar(im,ax=ax,shrink=.6,label='Metres')
    axes[0].set_ylabel('Latitude')
    fig.suptitle('Coastal diagnostic crop; white holes retain original nodata')
    out=path.with_suffix('.png');fig.savefig(out,dpi=150);plt.close(fig)
    return out


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('crop');print(plot(p.parse_args().crop))
