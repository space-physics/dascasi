"""
Created on Fri Sep 28 11:43:24 2018

@author: smrak
"""
import numpy as np
from scipy.spatial import Delaunay
import pymap3d as pm
import datetime
from typing import List, Tuple, Sequence
from scipy.interpolate import griddata
import xarray


def project_altitude(imgs: xarray.Dataset, mapping_altitude_km: float = None) -> xarray.Dataset:
    """
    takes image az/el and maps to lat/lon at a specific altitude
    This is a common approximation used for first-order auroral / airglow analysis
    """

    if mapping_altitude_km is None:
        mapping_altitude_km = 100
    # Get rid of NaNs in the coordinates' arrays
    eli = interpolateCoordinate(imgs['el'].values, method='nearest')
    azi = interpolateCoordinate(imgs['az'].values, method='nearest')
    # Convert Coordinates to WSG84
    lat, lon, alt = circular2lla(azi, eli, imgs.lat, imgs.lon,
                                 mapping_altitude_km=mapping_altitude_km)

    imgs.coords['mapping_lat'] = (('y', 'x'), lat)
    imgs.coords['mapping_lon'] = (('y', 'x'), lon)
    imgs.attrs['mapping_alt_km'] = mapping_altitude_km

    return imgs


def interpolateCoordinate(x: np.ndarray,
                          method: str = 'linear') -> np.ndarray:

    x0, y0 = np.meshgrid(np.arange(x.shape[0]),
                         np.arange(x.shape[1]))
    good = np.isfinite(x)
    x0 = x0[good]
    y0 = y0[good]
    X = x[good]

    x1, y1 = np.meshgrid(range(x.shape[0]), range(x.shape[1]))
    z = griddata((x0, y0), X.ravel(), (x1, y1), method=method)

    return z


def interpolate(values, vtx, wts, fill_value: float=np.nan) -> np.ndarray:
    ret = np.einsum('nj,nj->n', np.take(values, vtx), wts)
    ret[np.any(wts < 0, axis=1)] = fill_value

    return ret


def interpSpeedUpParams(x_in: np.ndarray, y_in: np.ndarray,
                        provisional_image: np.ndarray,
                        N: int = 512) -> Tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    The speedup is based on the scipy.interpolate.griddata algorithm. Thorough
    explanation is on the stackoverflow
    https://stackoverflow.com/questions/20915502/speedup-scipy-griddata-for-multiple-interpolations-between-two-irregular-grids
    """

    def _interpWeights(xyz: np.ndarray, uvw: np.ndarray, d: int=2):
        tri = Delaunay(xyz)
        simplex = tri.find_simplex(uvw)
        vertices = np.take(tri.simplices, simplex, axis=0)
        temp = np.take(tri.transform, simplex, axis=0)
        delta = uvw - temp[:, d]
        bary = np.einsum('njk,nk->nj', temp[:, :d, :], delta)
        return vertices, np.hstack((bary, 1 - bary.sum(axis=1, keepdims=True)))

    # Get the valid pixels, ie: np.isifinite() indexing
    mask = np.ma.masked_invalid(provisional_image)
    # Make new grids constraint by input longitude/latitude boundaries and resolution N
    xgrid, ygrid = np.mgrid[np.nanmin(x_in).min():np.nanmax(x_in).max():N*1j,  # type: ignore
                            np.nanmin(y_in).min():np.nanmax(y_in).max():N*1j]  # type: ignore
    # Old Coordinates maked -> flattern to a 1D array
    x = x_in[~mask.mask]
    y = y_in[~mask.mask]
    z = provisional_image[~mask.mask]
    # Old Coordinates -> To a single 2D array. x and y flattern to a single dimension
    xy = np.zeros([z.shape[0], 2])
    xy[:, 0] = x
    xy[:, 1] = y
    # New Coordinates: Flattern into a 2D array. Same as for the old
    uv = np.zeros([xgrid.shape[0]*xgrid.shape[1], 2])
    uv[:, 0] = xgrid.flatten()
    uv[:, 1] = ygrid.flatten()
    # trinagulate to new grids, simplex wights for a new grid, barycentric
    # coordinates based on the simplex for the new grids are computed
    vtx, wts = _interpWeights(xy, uv)
    # Return grids and interp points with associative weights
    return xgrid, ygrid, vtx, wts


def circular2lla(az: np.ndarray, el: np.ndarray,
                 lat0: float, lon0: float, alt0: float = 0,
                 mapping_altitude_km: float = 100.0) -> Tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Input mapping_altitude parameter is in kilometers [km] !
    Azimuth and elevation units must be in degrees !
    """
    # Convert mapping altitude to meters [m]
    malt_m = mapping_altitude_km * 1e3
    # Make sure the alt0 is a number
    alt0 = 0 if alt0 is None else alt0
    # Map to altitude
    r = malt_m / np.sin(np.radians(el))
    # Convert to WSG
    lat, lon, alt = pm.aer2geodetic(az, el, r, lat0, lon0, alt0)

    return lat, lon, alt


def datetime2posix(dtime: Sequence[datetime.datetime]) -> List[float]:
    """
    Convert an input list of datetime format timestamp to posix timestamp
    https://docs.python.org/3/library/datetime.html#datetime.datetime.timestamp
    """
    if isinstance(dtime, datetime.datetime):
        dtime = [dtime]

    return np.array([(t - datetime.datetime(1970, 1, 1)) / datetime.timedelta(seconds=1) for t in dtime])
