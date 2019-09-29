#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 28 11:43:24 2018

@author: smrak
"""
import xarray
import typing
import numpy as np
import datetime
from scipy.spatial import Delaunay
from pymap3d import aer2geodetic
from scipy.interpolate import griddata


def interpolateCoordinate(x: np.ndarray = None, N: int = 512, method: str = "linear"):
    if x is None:
        raise ValueError("Invalid input argument. Has to be a list or np.ndarray with a length of at least 1")
    x0, y0 = np.meshgrid(np.arange(x.shape[0]), np.arange(x.shape[1]))
    mask = np.ma.masked_invalid(x)
    x0 = x0[~mask.mask]
    y0 = y0[~mask.mask]
    X = x[~mask.mask]
    x1, y1 = np.meshgrid(np.arange(N), np.arange(N))
    z = griddata((x0, y0), X.ravel(), (x1, y1), method=method)
    return z


def interp(values, vtx, wts, fill_value=np.nan):
    ret = np.einsum("nj,nj->n", np.take(values, vtx), wts)
    ret[np.any(wts < 0, axis=1)] = fill_value
    return ret


def interpWeights(x, y, provisional_image: np.ndarray, d=2, N=512):
    xgrid, ygrid = np.meshgrid(np.linspace(np.nanmin(x), np.nanmax(x), N), np.linspace(np.nanmin(y), np.nanmax(y), N))
    # Old Coordinates maked -> flattern to a 1D array
    mask = np.ma.masked_invalid(provisional_image)
    x = x[~mask.mask]
    y = y[~mask.mask]
    z = provisional_image[~mask.mask]
    # Old Coordinates -> To a single 2D array. x and y flattern to a single dimension
    xy = np.zeros([z.size, 2])
    xy[:, 0] = x
    xy[:, 1] = y
    # New Coordinates: Flattern into a 2D array. Same as for the old
    uv = np.zeros([xgrid.shape[0] * xgrid.shape[1], 2])
    uv[:, 0] = xgrid.flatten()
    uv[:, 1] = ygrid.flatten()
    # trinagulate to new grids, simplex wights for a new grid, barycentric
    # coordinates based on the simplex for the new grids are computed
    tri = Delaunay(xy)
    simplex = tri.find_simplex(uv)
    vertices = np.take(tri.simplices, simplex, axis=0)
    temp = np.take(tri.transform, simplex, axis=0)
    delta = uv - temp[:, d]
    bary = np.einsum("njk,nk->nj", temp[:, :d, :], delta)
    vtx = vertices
    wts = np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

    return vtx, wts


def interpSpeedUp(x_in: np.ndarray = None, y_in: np.ndarray = None, image: np.ndarray = None, N: int = 512, verbose: bool = True):
    """
    The speedup is based on the scipy.interpolate.griddata algorithm. Thorough
    explanation is on the stackoverflow
    https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-
    multiple-interpolations-between-two-irregular-grids
    """

    def _interpolate(values, vtx, wts, fill_value=np.nan):
        ret = np.einsum("nj,nj->n", np.take(values, vtx), wts)
        ret[np.any(wts < 0, axis=1)] = fill_value
        return ret

    def _interpWeights(xyz, uvw, d=2):
        tri = Delaunay(xyz)
        simplex = tri.find_simplex(uvw)
        vertices = np.take(tri.simplices, simplex, axis=0)
        temp = np.take(tri.transform, simplex, axis=0)
        delta = uvw - temp[:, d]
        bary = np.einsum("njk,nk->nj", temp[:, :d, :], delta)
        return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

    # Get the valid pixels, ie: np.isifinite() indexing
    if len(image.shape) > 2:
        provisional_image = image[0]
    else:
        provisional_image = image
    mask = np.ma.masked_invalid(provisional_image)
    # Make new grids constraint by input longitude/latitude boundaries and
    # resolution N
    xgrid, ygrid = np.meshgrid(np.linspace(np.nanmin(x_in), np.nanmax(x_in), N), np.linspace(np.nanmin(y_in), np.nanmax(y_in), N))
    # Old Coordinates maked -> flattern to a 1D array
    x = x_in[~mask.mask]
    y = y_in[~mask.mask]
    z = provisional_image[~mask.mask]
    # Old Coordinates -> To a single 2D array. x and y flattern to a single dimension
    xy = np.zeros([z.size, 2])
    xy[:, 0] = x
    xy[:, 1] = y
    # New Coordinates: Flattern into a 2D array. Same as for the old
    uv = np.zeros([xgrid.shape[0] * xgrid.shape[1], 2])
    uv[:, 0] = xgrid.flatten()
    uv[:, 1] = ygrid.flatten()
    # trinagulate to new grids, simplex wights for a new grid, barycentric
    # coordinates based on the simplex for the new grids are computed
    if verbose:
        print("Computing intepolation weights")
    t0 = datetime.datetime.now()
    vtx, wts = _interpWeights(xy, uv)
    # Interpolate N images
    Zim = np.copy(image) * np.nan
    if len(image.shape) > 2:
        if verbose:
            print("Interpolating {} images".format(image.shape[0]))
        for i in range(image.shape[0]):
            if verbose:
                print("Image {}/{}".format(i + 1, image.shape[0]))
            Zim[i] = _interpolate(image[i], vtx, wts).reshape(N, N)
    else:
        Zim = _interpolate(image, vtx, wts).reshape(N, N)
    if verbose:
        print("Interpolation done in {} seconds".format((datetime.datetime.now() - t0).total_seconds()))
    # Return grids and interp points with associative weights
    return xgrid, ygrid, Zim


def circular2lla(
    az: np.ndarray = None,
    el: np.ndarray = None,
    lon0: float = None,
    lat0: float = None,
    alt0: float = 0,
    mapping_altitude: float = 100.0,
):
    """
    Input mappiing_altitude parameter is in kilometers [km] !
    Azimuth and elevation units must be in degrees !
    """
    # Convert mapping altitude to meters [m]
    mapping_altitude = mapping_altitude * 1e3
    # Make sure the alt0 is a number
    alt0 = 0 if (alt0 is None) else alt0
    # Map to altitude
    r = mapping_altitude / np.sin(np.deg2rad(el))
    # Convert to WSG
    lat, lon, alt = aer2geodetic(az, el, r, lat0, lon0, alt0)
    #
    return lat, lon, alt


def getPixelBrightness(
    D: xarray.Dataset = None,
    obstimes: np.ndarray = None,
    obs_lon: np.ndarray = None,
    obs_lat: np.ndarray = None,
    coordinates: bool = False,
):
    assert isinstance(obstimes[0], datetime.datetime)
    # Coordinates
    dasc_lat = D.lat.values
    dasc_lon = D.lon.values
    # DASC obstimes
    dasc_time = np.array([datetime.datetime.utcfromtimestamp(t) for t in D.time.values])
    # Filter according to time requirement
    treq = (dasc_time >= obstimes[0]) & (dasc_time <= obstimes[-1])
    Dtimes = dasc_time[treq]

    # GET BRIGHTNESS
    # fixed position
    if isinstance(obs_lon, (int, float)) and isinstance(obs_lat, (int, float)):
        idx = abs(obs_lon - dasc_lon[0, :]).argmin()
        idy = abs(obs_lat - dasc_lat[:, 0]).argmin()
        lon = dasc_lon[0, idx]
        lat = dasc_lat[idy, 0]
        pixel_brightness = D.image[treq][:, idx, idy].values
    # non-stationary point
    else:
        lon = []
        lat = []
        pixel_brightness = []
        for i in range(Dtimes.shape[0]):
            idt = abs(Dtimes[i] - obstimes).argmin()
            idx = abs(obs_lon[idt] - dasc_lon[0, :]).argmin()
            idy = abs(obs_lat[idt] - dasc_lat[:, 0]).argmin()
            lon.append(dasc_lon[0, idx])
            lat.append(dasc_lat[idy, 0])
            pixel_brightness.append(D.image[treq][i, idx, idy].values)
        pixel_brightness = np.array(pixel_brightness)

    # Return brightness time-series
    if coordinates:
        return Dtimes, pixel_brightness, lon, lat
    else:
        return Dtimes, pixel_brightness


def getDASCimage(D: xarray.Dataset = None, ix: typing.Union[datetime.datetime, int] = None, coordinate: str = "wsg"):
    assert isinstance(ix, (datetime.datetime, int))
    # Find the closes entrance for the given timestamp
    if isinstance(ix, datetime.datetime):
        T = datetime2posix(ix)[0]
        Di = D.sel(time=T, method="nearest")
        dasc_dt = datetime.datetime.utcfromtimestamp(Di.time.values)
        if coordinate == "polar":
            img = Di.polar.values
        elif coordinate == "wsg":
            img = Di.image.values
    # Find for the given index
    if isinstance(ix, int):
        dasc_dt = datetime.datetime.utcfromtimestamp(D.time.values[ix])
        if coordinate == "polar":
            img = D.polar.values[ix]
        elif coordinate == "wsg":
            img = D.image.values[ix]

    return dasc_dt, img


def datetime2posix(dtime):
    """
    Convert an input list of datetime format timestamp to posix timestamp
    """
    if isinstance(dtime, datetime.datetime):
        dtime = [dtime]
    return np.array([i.replace(tzinfo=datetime.timezone.utc).timestamp() for i in dtime])
