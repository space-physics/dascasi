from pathlib import Path
import xarray
import numpy as np
from datetime import timedelta, datetime
from matplotlib.pyplot import draw, pause, figure
from matplotlib.colors import LogNorm
#
try:
    import themisasi.plots as themisplot
except ImportError:
    themisplot = None


def histogram_dasc(imgs: xarray.Dataset, odir=None):
    """
    creates per wavelength histograms
    the entries in list img correspond to wavelength, a 1-D array
    """
    if odir is not None:
        odir = Path(odir).expanduser()

    fg = figure(figsize=(15, 5))
    axs = fg.subplots(1, 3)
    for a, i in zip(axs, imgs.data_vars):
        a.hist(imgs[i].dropna(dim='time', how='all').values.ravel(), bins=128)
        a.set_yscale('log')
        a.set_title(f'$\lambda={i}$ nm')
        a.set_xlabel('14-bit data numbers')

    if odir:
        ofn = odir/'DASChistogram.png'
        print('writing', ofn, end='\r')
        fg.savefig(ofn, bbox_inches='tight')


def moviedasc(imgs: xarray.Dataset, odir: Path, cadence: float, rows=None, cols=None):

    if odir:
        print('writing to', odir)
        odir = Path(odir).expanduser()

    wavelen = list(imgs.data_vars)

    fg = figure(figsize=(15, 5))

    axs = np.atleast_1d(fg.subplots(1, len(np.unique(wavelen))))

    if imgs.time.dtype == 'M8[ns]':
        time = [datetime.utcfromtimestamp(t/1e9) for t in imgs.time.values.astype(int)]
    else:
        time = imgs.time.values.astype(datetime)
# %% setup figures
    if 'unknown' not in imgs.data_vars:
        Hi = []
        Ht = []
        for ax, w, mm, c in zip(axs,
                                np.unique(wavelen),
                                ((350, 800), (350, 9000), (350, 900)),
                                ('b', 'g', 'r')):
            # ax.axis('off') #this also removes xlabel,ylabel
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel(f'{w} nm', color=c)

            Hi.append(ax.imshow(imgs[w].dropna(dim='time', how='all')[0],
                                vmin=mm[0], vmax=mm[1],
                                origin='lower',
                                norm=LogNorm(), cmap='gray'))

            Ht.append(ax.set_title('', color=c))
            # fg.colorbar(hi[-1],ax=a).set_label('14-bit data numbers')
            if themisplot is not None:
                themisplot.overlayrowcol(ax, rows, cols)

        fg.tight_layout(h_pad=1.08)  # get rid of big white space in between figures
    else:
        ax = axs[0]
        ax.set_xticks([])
        ax.set_yticks([])
        hi = ax.imshow(imgs['unknown'][0],
                       vmin=(350, 10000),
                       origin='lower',
                       norm=LogNorm(), cmap='gray')

        ht = ax.set_title('')
        if themisplot is not None:
            themisplot.overlayrowcol(ax, rows, cols)
# %% loop
    print('generating video until', time[-1])
    t = time[0]
    dt = timedelta(seconds=cadence)
    while t <= time[-1]:
        if 'unknown' not in imgs.data_vars:
            for w, hi, ht in zip(np.unique(wavelen), Hi, Ht):
                im = imgs[w].dropna(dim='time', how='all').sel(time=t, method='nearest')
                hi.set_data(im)
                try:
                    ht.set_text(str(im.time.values))
                except OSError:  # file had corrupted time
                    ht.set_text('')
        else:
            im = imgs['unknown'].sel(time=t, method='nearest')
            hi.set_data(im)
            try:
                ht.set_text(str(im.time.values))
            except OSError:  # file had corrupted time
                ht.set_text('')

        draw(), pause(0.05)  # the pause avoids random crashes
        t += dt

        if odir:
            ofn = odir / (str(t)+'.png')
            print('saving', ofn, end='\r')
            fg.savefig(ofn, bbox_inches='tight', facecolor='k')
