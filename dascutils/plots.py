from . import Path
from os import devnull
from numpy import arange
from pytz import UTC
from datetime import datetime
from scipy.interpolate import interp1d
from matplotlib.pyplot import draw,pause,subplots
from matplotlib.colors import LogNorm
import matplotlib.animation as anim
#
from themisasi.plots import overlayrowcol

DPI=100

def histdasc(img,wavelength,odir=None):
    """
    creates per wavelength histograms
    the entries in list img correspond to wavelength, a 1-D array
    """
    if odir: odir =Path(odir).expanduser()
    assert len(wavelength)==len(img) #works for 3-D img as well, assuming C-order.

    fg,axs = subplots(1,3,figsize=(15,5))
    for a,I,w in zip(axs,img,wavelength):
        a.hist(I.ravel(),bins=128)
        a.set_yscale('log')
        a.set_title('$\lambda={}$ nm'.format(w))
        a.set_xlabel('14-bit data numbers')

    if odir:
        ofn = odir/'DASChistogram.png'
        print('writing {}'.format(ofn))
        fg.savefig(str(ofn),bbox_inches='tight',dpi=100)

def moviedasc(img,wavelength,times,odir,cadence,rows=None,cols=None):

    if odir:
        ofn = Path(odir).expanduser()/'DASC_{}.avi'.format(times[0][0][0][0])
        write=True
        print('writing {}'.format(ofn))
    else:
        ofn = devnull
        write=False

    Writer = anim.writers['ffmpeg']
    writer = Writer(fps=5,
                    codec='ffv1')

    fg,axs = subplots(1,3,figsize=(15,5))
    hi = []; ht=[]
    for a,w,x,mm,c in zip(axs,wavelength,(0.225,0.5,0.775),
                     ((350,800),(350,9000),(350,900)),('b','g','r')):
        #a.axis('off') #this also removes xlabel,ylabel
        a.set_xticks([]); a.set_yticks([])
        a.set_xlabel('{} nm'.format(w),color=c)
        hi.append(a.imshow(img[0][0],vmin=mm[0],vmax=mm[1],
                           origin='lower',
                           norm=LogNorm(),cmap='gray'))
        ht.append(a.set_title('',color=c))
        #fg.colorbar(hi[-1],ax=a).set_label('14-bit data numbers')
        overlayrowcol(a,rows,cols)

    T = max([t[0,0] for t in times])
    Tmax = min([t[-1,0] for t in times])

    fg.tight_layout(h_pad=1.08) #get rid of big white space in between figures
#%% loop
    with writer.saving(fg, str(ofn), DPI):
        while T<=Tmax:
            for I,Hi,Ht,t in zip(img,hi,ht,times):
                ft = interp1d(t[:,0],arange(len(t)),kind='nearest')
                ind = ft(T).astype(int)
                #print(ind,end=' ')
                Hi.set_data(I[ind])
                try:
                    Ht.set_text(datetime.fromtimestamp(t[ind,0],tz=UTC))
                except OSError: #file had corrupted time
                    Ht.set_text('')

            draw(), pause(0.05) # the pause avoids random crashes
            T += cadence
            if write:
                writer.grab_frame(facecolor='k')
