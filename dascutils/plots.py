from pathlib import Path
import numpy as np
import typing
from datetime import timedelta, datetime
import xarray
from matplotlib.pyplot import draw, pause, figure
from matplotlib.colors import LogNorm

#
try:
    import themisasi.plots as themisplot
except ImportError:
    themisplot = None


def plot_projected_image(imgs: xarray.DataArray):
    """
    plots a projected image at an altitude (lat, lon)

    https://matplotlib.org/3.1.0/tutorials/colors/colormaps.html#sequential
    """
    if not isinstance(imgs, xarray.DataArray):
        raise TypeError("Please pass in one image wavelength at a time")

    cmap = {"0428": "Blues", "0558": "Greens", "0630": "Reds"}

    for img in imgs:
        fg = figure()
        ax = fg.gca()
        ax.pcolormesh(imgs.lon, imgs.lat, imgs[0].values, cmap=cmap.get(imgs.name, "Grays"))
        ax.set_title(f"{str(img.time.values)[:-10]}: {imgs.name} " r"$\AA$" f"at {imgs.mapping_alt_km} km altitude")
        ax.set_xlabel("geographic longitude")
        ax.set_ylabel("geographic latitude")


def histogram_dasc(imgs: typing.Dict[str, typing.Any], outdir=None):
    """
    creates per wavelength histograms
    the entries in list img correspond to wavelength, a 1-D array
    """

    fg = figure(figsize=(15, 5))
    axs = fg.subplots(1, 3)
    for a, i in zip(axs, imgs["wavelengths"]):
        a.hist(imgs[i].values.ravel(), bins=128)
        a.set_yscale("log")
        a.set_title(r"$\lambda=" + f"{i}$ nm")
        a.set_xlabel("14-bit data numbers")

    if outdir:
        outdir = Path(outdir).expanduser()
        ofn = outdir / "DASChistogram.png"
        print("writing", ofn, end="\r")
        fg.savefig(ofn, bbox_inches="tight")


def moviedasc(imgs: typing.Dict[str, typing.Any], outdir: Path, cadence: float, rows=None, cols=None):

    wavlen = imgs["wavelengths"]

    fg = figure(figsize=(15, 5))

    axs = np.atleast_1d(fg.subplots(1, len(wavlen)))

    # %% setup figures
    if "unknown" not in wavlen:
        Hi = []
        Ht = []
        for ax, w, mm, c in zip(axs, wavlen, ((350, 800), (350, 9000), (350, 900)), ("b", "g", "r")):
            # ax.axis('off') #this also removes xlabel,ylabel
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_xlabel(f"{w} nm", color=c)

            Hi.append(ax.imshow(imgs[w][0], vmin=mm[0], vmax=mm[1], origin="lower", norm=LogNorm(), cmap="gray"))

            Ht.append(ax.set_title("", color=c))
            # fg.colorbar(hi[-1],ax=a).set_label('14-bit data numbers')
            if themisplot is not None:
                themisplot.overlayrowcol(ax, rows, cols)

        fg.tight_layout(h_pad=1.08)  # get rid of big white space in between figures
    else:
        ax = axs[0]
        ax.set_xticks([])
        ax.set_yticks([])
        hi = ax.imshow(imgs["unknown"][0], vmin=(350, 10000), origin="lower", norm=LogNorm(), cmap="gray")

        ht = ax.set_title("")
        if themisplot is not None:
            themisplot.overlayrowcol(ax, rows, cols)
    # %% loop
    t = min([imgs[wl]["time"][0] for wl in wavlen]).values.astype("datetime64[us]").astype(datetime)
    t1 = max([imgs[wl]["time"][-1] for wl in wavlen]).values.astype("datetime64[us]").astype(datetime)
    dt = timedelta(seconds=cadence)

    while t <= t1:
        if "unknown" not in wavlen:
            for w, hi, ht in zip(wavlen, Hi, Ht):
                im = imgs[w].sel(time=t, method="nearest")
                _update_panel(im, hi, ht)
        else:
            im = imgs["unknown"].sel(time=t, method="nearest")
            _update_panel(im, hi, ht)

        draw(), pause(0.05)  # the pause avoids random crashes

        if outdir:
            outdir = Path(outdir).expanduser()
            print("writing to", outdir)
            ofn = outdir / (str(t) + ".png")
            print("saving", ofn, end="\r")
            fg.savefig(ofn, bbox_inches="tight", facecolor="k")

        t += dt


def _update_panel(im, hi, ht):
    try:
        hi.set_data(im)
    except AttributeError:
        hi.set_array(im.values.ravel())

    try:
        ht.set_text(str(im.time.values))
    except OSError:  # file had corrupted time
        ht.set_text("")
