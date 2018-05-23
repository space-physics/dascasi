from pathlib import Path
import xarray
import numpy as np
from datetime import timedelta,datetime
from matplotlib.pyplot import draw,pause,figure,close
from matplotlib.colors import LogNorm
#
try:
    import themisasi.plots as themisplot
except ImportError:
    themisplot = None


def histogram_dasc(imgs, odir=None):
    """
    creates per wavelength histograms
    the entries in list img correspond to wavelength, a 1-D array
    """
    if odir is not None:
        odir = Path(odir).expanduser()

    anyplot=False
    fg = figure(figsize=(15,5))
    axs = fg.subplots(1,3)
    for a,i in zip(axs,imgs.data_vars):
        if isinstance(i,str): # FIXME to reject coordinates (remove by xarray 0.11)
            continue
        else:
            anyplot=True
        a.hist(imgs[i].dropna(dim='time',how='all').values.ravel(), bins=128)
        a.set_yscale('log')
        a.set_title(f'$\lambda={i}$ nm')
        a.set_xlabel('14-bit data numbers')

    if not anyplot:
        close(fg)

    if odir:
        ofn = odir/'DASChistogram.png'
        print('writing',ofn, end='\r')
        fg.savefig(ofn, bbox_inches='tight')


def moviedasc(imgs:xarray.Dataset, odir:Path, cadence:float, rows=None, cols=None):

    if odir:
        print('writing to',odir)
        odir = Path(odir).expanduser()

    fg = figure(figsize=(15,5))

    if imgs.wavelength is not None:
        axs = np.atleast_1d(fg.subplots(1,len(np.unique(imgs.wavelength))))
    else:
        axs = [fg.gca()]

    hi = []; ht=[]
    if imgs.time.dtype == 'M8[ns]':
        time = [datetime.utcfromtimestamp(t/1e9) for t in imgs.time.values.astype(int)]
    else:
        time = imgs.time.values.astype(datetime)
# %% setup figures
    if 'white' not in imgs.data_vars:
        for ax,w,mm,c in zip(axs,
                           np.unique(imgs.wavelength),
                           ((350,800),(350,9000),(350,900)),
                           ('b','g','r')):
            #ax.axis('off') #this also removes xlabel,ylabel
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel(f'{w} nm', color=c)

            hi.append(ax.imshow(imgs[w].dropna(dim='time',how='all')[0],
                               vmin=mm[0],vmax=mm[1],
                               origin='lower',
                               norm=LogNorm(),cmap='gray'))

            ht.append(ax.set_title('', color=c))
            #fg.colorbar(hi[-1],ax=a).set_label('14-bit data numbers')
            if themisplot is not None:
                themisplot.overlayrowcol(ax, rows, cols)

        fg.tight_layout(h_pad=1.08) #get rid of big white space in between figures
    else:
        ax = axs[0]
        ax.set_xticks([])
        ax.set_yticks([])
        hi = ax.imshow(imgs.white[0],
                       vmin=(350,10000),
                       origin='lower',
                       norm=LogNorm(),cmap='gray')

        ht = ax.set_title('')
        if themisplot is not None:
            themisplot.overlayrowcol(ax, rows, cols)
#%% loop
    print('generating video until', time[-1])
    t = time[0]
    dt = timedelta(seconds=cadence)
    while t <= time[-1]:
        if 'white' not in imgs.data_vars:
            for w,Hi,Ht in zip(np.unique(imgs.wavelength),hi,ht):
                I = imgs[w].dropna(dim='time',how='all').sel(time=t, method='nearest')
                Hi.set_data(I)
                try:
                    Ht.set_text(str(I.time.values))
                except OSError: #file had corrupted time
                    Ht.set_text('')
        else:
            I = imgs.white.sel(time=t,method='nearest')
            hi.set_data(I)
            try:
                ht.set_text(str(I.time.values))
            except OSError: #file had corrupted time
                ht.set_text('')

        draw(), pause(0.05) # the pause avoids random crashes
        t += dt

        if odir:
            ofn = odir / (str(t)+'.png')
            print('saving', ofn,end='\r')
            fg.savefig(ofn, bbox_inches='tight',facecolor='k')
