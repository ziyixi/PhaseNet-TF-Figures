""" 
slice.py

helper functions in cuting cross-sections, make projections.
"""
from typing import List, Tuple

import numpy as np
import pyproj
import xarray as xr
from obspy.geodetics import locations2degrees
from obspy.geodetics.base import degrees2kilometers
from scipy.interpolate import RegularGridInterpolator
from scipy.spatial import KDTree

from . import generate_tmp_file


def model_interp(
    to_interp_data: xr.DataArray, lons: np.ndarray, lats: np.ndarray, deps: np.ndarray
) -> np.ndarray:
    """Give an xarray model, interp it based on the given lats, lons, deps and construct a new xarray dataset.
    mainly used to generate the vertical cross-sections

    Args:
        to_interp_data (xr.DataArray): the data array to interp
        lons (np.ndarray): the longitude array
        lats (np.ndarray): the latitude array, define a line with lons on the plane
        deps (np.ndarray): the depth array

    Returns:
        np.ndarray: the interp result
    """
    # * len(lons) should be the same as len(lats)
    profile_list = []
    for idep in range(len(deps)):
        for ilon in range(len(lons)):
            profile_list.append([lons[ilon], lats[ilon], deps[idep]])
    model_interpolating_function = RegularGridInterpolator(
        (
            to_interp_data.longitude.data,
            to_interp_data.latitude.data,
            to_interp_data.depth.data,
        ),
        to_interp_data.data,
    )
    interp_result: np.ndarray = model_interpolating_function(profile_list)
    cross_section = np.zeros((len(lons), len(deps)))

    icount = 0
    for idep in range(len(deps)):
        for ilon in range(len(lons)):
            cross_section[ilon, idep] = interp_result[icount]
            icount += 1

    return cross_section


def topo_interp(
    to_interp_data: xr.DataArray, lons: np.ndarray, lats: np.ndarray
) -> np.ndarray:
    """Give the xarray topography model, interp the elevation line along the given (lons,lats) pair.

    Args:
        to_interp_data (xr.DataArray): the input topo data array
        lons (np.ndarray): the longitude array
        lats (np.ndarray): the latitude array, define a line with lons on the plane

    Returns:
        np.ndarray: the interp topo result
    """
    profile_list = []
    for ilon in range(len(lons)):
        profile_list.append([lons[ilon], lats[ilon]])
    # the names and the transverse might be adjusted, this is the gmt format
    grd_interpolating_function = RegularGridInterpolator(
        (to_interp_data.lon.data, to_interp_data.lat.data), to_interp_data.data.T
    )

    grd_interp_result = grd_interpolating_function(profile_list)

    # * return the 1d array
    return grd_interp_result


def slab_interp(
    to_interp_data: xr.DataArray, lons: np.ndarray, lats: np.ndarray
) -> np.ndarray:
    """generate the depth of the slab interface along a (lons,lats) track

    Args:
        to_interp_data (xr.DataArray): the loaded slab interface
        lons (np.ndarray): the lons track
        lats (np.ndarray): the lats track

    Returns:
        np.ndarray: the depth track
    """
    profile_list = []
    for ilon in range(len(lons)):
        profile_list.append([lons[ilon], lats[ilon]])
    # the names and the transverse might be adjusted, this is the gmt format
    grd_interpolating_function = RegularGridInterpolator(
        (to_interp_data.x.data, to_interp_data.y.data),
        -to_interp_data.z.data.T,
        bounds_error=False,
    )

    grd_interp_result = grd_interpolating_function(profile_list)

    # * return the 1d array
    return grd_interp_result


def gmt_lat_as_dist(
    start: Tuple[float, float],
    end: Tuple[float, float],
    a_interval: float,
    g_interval: float,
    npts: int = 1001,
) -> str:
    """Generate a lebel tmp file for pxc[file name], so we can have evenly sampled ticks in great circle represented as lat/lon

    Args:
        start (Tuple[float, float]): the start position of the line (lon,lat)
        end (Tuple[float, float]): the end position of the line (lon,lat)
        a_interval (float): the interval for first level ticks
        g_interval (float): the interval for second level ticks
        npts (int, optional): control the relative error for the returned result

    Returns:
        str: a tmp file path used for frame
    """
    # * init parameters
    startlon, startlat = start
    endlon, endlat = end
    # * calculate
    g = pyproj.Geod(ellps="WGS84")
    gcarc = locations2degrees(startlat, startlon, endlat, endlon)
    if startlat > endlat:
        startlon, endlon = endlon, startlon
        startlat, endlat = endlat, startlat
    test_points = np.array(
        g.npts(startlon, startlat, endlon, endlat, (npts - 1) * 10 + 1)
    )
    tree = KDTree(test_points[:, 1].reshape(test_points.shape[0], -1))

    starta = (
        startlat
        if startlat % a_interval == 0
        else (startlat // a_interval + 1) * a_interval
    )
    enda = endlat if endlat % a_interval == 0 else (endlat // a_interval) * a_interval
    a_list = np.arange(starta, a_interval + enda, a_interval).astype(int)
    # * generate the gmt custome axis file when use lon to plot the cross-section
    num_list = []
    type_list = []
    annote_list = []
    # before first a
    for each_lat in np.arange(a_list[0] - g_interval, startlat, -1 * g_interval)[::-1]:
        num_list.append(each_lat)
        type_list.append("f")
        annote_list.append("")
    # start from first a
    for ia in range(len(a_list) - 1):
        starta = a_list[ia]
        enda = a_list[ia + 1]
        # first a
        num_list.append(starta)
        type_list.append("a")
        annote_list.append(f"{starta}")
        for each_g in np.arange(starta + g_interval, enda, g_interval):
            num_list.append(each_g)
            type_list.append("f")
            annote_list.append("")
    # last a
    num_list.append(a_list[-1])
    type_list.append("a")
    annote_list.append(f"{a_list[-1]}")
    for each_g in np.arange(a_list[-1] + g_interval, endlat, g_interval):
        num_list.append(each_g)
        type_list.append("f")
        annote_list.append("")

    # convert num_list to actual dist_list
    num_list = np.array(num_list)
    _, pos = tree.query(num_list.reshape(len(num_list), -1))
    dist_list = pos / ((npts - 1) * 10) * gcarc

    # write lists to a temp file
    tmp = generate_tmp_file()
    with open(tmp, "w") as f:
        for index in range(len(num_list)):
            f.write(f"{dist_list[index]}  {type_list[index]}  {annote_list[index]} \n")
    return tmp


def gmt_lon_as_dist(
    start: Tuple[float, float],
    end: Tuple[float, float],
    a_interval: float,
    g_interval: float,
    npts: int = 1001,
) -> str:
    """Generate a lebel tmp file for pxc[file name], so we can have evenly sampled ticks in great circle represented as lat/lon

    Args:
        start (Tuple[float, float]): the start position of the line (lon,lat)
        end (Tuple[float, float]): the end position of the line (lon,lat)
        a_interval (float): the interval for first level ticks
        g_interval (float): the interval for second level ticks
        npts (int, optional): control the relative error for the returned result

    Returns:
        str: a tmp file path used for frame
    """
    # * init parameters
    startlon, startlat = start
    endlon, endlat = end
    # * calculate
    g = pyproj.Geod(ellps="WGS84")
    gcarc = locations2degrees(startlat, startlon, endlat, endlon)
    if startlon > endlon:
        startlon, endlon = endlon, startlon
        startlat, endlat = endlat, startlat
    test_points = np.array(
        g.npts(startlon, startlat, endlon, endlat, (npts - 1) * 10 + 1)
    )
    tree = KDTree(test_points[:, 0].reshape(test_points.shape[0], -1))

    starta = (
        startlon
        if startlon % a_interval == 0
        else (startlon // a_interval + 1) * a_interval
    )
    enda = endlon if endlon % a_interval == 0 else (endlon // a_interval) * a_interval
    a_list = np.arange(starta, a_interval + enda, a_interval).astype(int)
    # * generate the gmt custome axis file when use lon to plot the cross-section
    num_list = []
    type_list = []
    annote_list = []
    # before first a
    for each_lon in np.arange(a_list[0] - g_interval, startlon, -1 * g_interval)[::-1]:
        num_list.append(each_lon)
        type_list.append("f")
        annote_list.append("")
    # start from first a
    for ia in range(len(a_list) - 1):
        starta = a_list[ia]
        enda = a_list[ia + 1]
        # first a
        num_list.append(starta)
        type_list.append("a")
        annote_list.append(f"{starta}")
        for each_g in np.arange(starta + g_interval, enda, g_interval):
            num_list.append(each_g)
            type_list.append("f")
            annote_list.append("")
    # last a
    num_list.append(a_list[-1])
    type_list.append("a")
    annote_list.append(f"{a_list[-1]}")
    for each_g in np.arange(a_list[-1] + g_interval, endlon, g_interval):
        num_list.append(each_g)
        type_list.append("f")
        annote_list.append("")

    # convert num_list to actual dist_list
    num_list = np.array(num_list)
    _, pos = tree.query(num_list.reshape(len(num_list), -1))
    dist_list = pos / ((npts - 1) * 10) * gcarc

    # write lists to a temp file
    tmp = generate_tmp_file()
    with open(tmp, "w") as f:
        for index in range(len(num_list)):
            f.write(f"{dist_list[index]}  {type_list[index]}  {annote_list[index]} \n")
    return tmp


def extend_line(
    start: Tuple[float, float], end: Tuple[float, float], length: float
) -> Tuple[float, float]:
    """Extend the current line to the specified length

    Args:
        start (Tuple[float,float]): the coordinate for the starting position
        end (Tuple[float,float]): the coordinate for the current ending position
        length (float): the expected new line length in degree

    Returns:
        Tuple[float,float]: the new ending position
    """
    startlon, startlat = start
    endlon, endlat = end
    g = pyproj.Geod(ellps="WGS84")
    az, _, _ = g.inv(startlon, startlat, endlon, endlat)
    newlon, newlat, _ = g.fwd(startlon, startlat, az, degrees2kilometers(length) * 1000)
    return newlon, newlat
