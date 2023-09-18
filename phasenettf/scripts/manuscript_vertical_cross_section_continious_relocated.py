"""
manuscript_vertical_cross_section_continious_relocated.py:

This script is used to generate the vertical cross section of the relocated catalog of the cntinious test.
"""
import json
from pathlib import Path
from string import ascii_lowercase

import numpy as np
import pandas as pd
import pygmt
import xarray as xr
from obspy.geodetics.base import degrees2kilometers
from scipy import interpolate

from phasenettf import resource, save_path
from phasenettf.utils.slice import extend_line, model_interp, slab_interp

LENGTH = 7
REGION = [-184, -172, -26, -14]
LINES = [
    # [-179, -16, -173, -19],
    [-179.5, -17, -173.5, -20],
    [-180, -18, -174, -21],
    [-180.5, -19, -174.5, -22],
    [-181, -20, -175, -23],
]


def main():
    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="14p",
        MAP_LABEL_OFFSET="12p",
        FONT_ANNOT_PRIMARY="12p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="12p",
        FONT_TITLE="14p,black",
        MAP_FRAME_PEN="1p,black",
    )

    fig.shift_origin(yshift="1.2i")

    # * plot map
    fig.basemap(
        region=REGION,
        projection="M3.5i",
        frame=["WSen", "xaf+lLongitude (degree)", "yaf+lDepth (km)"],
    )
    df_events = load_df_events()
    plot_earth_relief(fig)
    plot_lines(fig)
    plot_events(fig, df_events)
    plot_inset(fig)

    fig.shift_origin(xshift="4.5i", yshift="-0.5i")

    # * some preparation for plotting
    slab_model = xr.open_dataset(
        resource(["slab2", f"ker_slab2_depth.grd"], normal_path=True)
    )
    model = load_tomodd_model()
    ak135 = load_ak135("vp", model)
    model = (model / ak135 - 1) * 100
    pygmt.makecpt(
        cmap=resource(["cpt", "dvs_6p_nan.cpt"]),
        series=f"-3/3/1",
        continuous=True,
        background="o",
    )

    # * plot vertical cross-sections
    with fig.subplot(
        nrows=2,
        ncols=2,
        subsize=("3.85i", "1.95i"),
        margins=["0.1i", "0.1i"],
        frame=["WSen", "xaf+lDistance (km)", "yaf+lDepth (km)"],
        sharex="b",
        sharey="l",
    ):
        for ipanel in range(4):
            with fig.set_panel(panel=ipanel):
                plot_vertical_cross_section(fig, ipanel, slab_model, model, df_events)
    fig.colorbar(position="JMR+o0.8c/0c+w3i/0.15i", frame=["x+l@~d@~lnV_P(%)"])

    save_path(fig, Path(__file__).resolve().stem)


# * =========================Helper functions=========================
def plot_earth_relief(fig: pygmt.Figure):
    grd_topo = pygmt.datasets.load_earth_relief(
        resolution="02m", region=REGION, registration="gridline"
    )
    assert type(grd_topo) == xr.DataArray
    # plot 2000m contour of topography, start from -2000m to -10000m
    fig.grdcontour(
        grd_topo,
        interval=1000,
        pen="0.8p,black",
        limit="-10000/-7000",
    )
    # plot only -1000m contour
    fig.grdcontour(
        grd_topo,
        interval=1000,
        pen="0.5p,black",
        limit="-1100/-1000",
    )
    # fig.text(
    #     position="TL",
    #     text="(a)",
    #     font="18p,Helvetica-Bold,black",
    # )


def plot_inset(fig: pygmt.Figure):
    with fig.inset(position="jBL+w7c+o0.2c", margin=0):
        fig.coast(
            region="g",
            projection=f"W{(REGION[0]+REGION[1])//2}/3c",
            land="gray",
            water="white",
            frame=["wsen", "xafg", "yafg"],
        )
        fig.plot(
            x=[REGION[0], REGION[1], REGION[1], REGION[0], REGION[0]],
            y=[REGION[2], REGION[2], REGION[3], REGION[3], REGION[2]],
            pen="1p,black",
            projection=f"W{(REGION[0]+REGION[1])//2}/3c",
        )


def plot_events(fig: pygmt.Figure, df_events: pd.DataFrame):
    pygmt.makecpt(cmap="jet", series=[0, 700, 1], continuous=True, reverse=True)
    fig.plot(
        x=df_events["longitude"],
        y=df_events["latitude"],
        style="c0.05c",
        fill=df_events["depth"],
        cmap=True,
    )
    # plot the colorbar to the right of the map
    fig.colorbar(position="JBC+o0c/0.8c+w3i/0.15i", frame=["x+lDepth (km)"])


def plot_lines(fig: pygmt.Figure):
    for iline, line in enumerate(LINES):
        startlon, startlat, endlon, endlat = line
        endlon, endlat = extend_line((startlon, startlat), (endlon, endlat), LENGTH)
        fig.plot(
            x=[startlon, endlon],
            y=[startlat, endlat],
            pen="2p,magenta",
        )
        fig.text(
            x=endlon,
            y=endlat,
            text=ascii_lowercase[iline],
            font="15p,Helvetica-Bold,magenta",
            justify="LM",
        )


def load_df_events():
    relocated_df = pd.read_csv(
        resource(["catalog", "continious_relocated.csv"], normal_path=True),
        usecols=[1, 2, 3],
        names=["latitude", "longitude", "depth"],
        sep=r"\s+",
    )
    return relocated_df


def load_tomodd_model():
    with open(
        resource(name=["model", "TX2019_global.mod"], normal_path=True, check=True)
    ) as f:
        lines = f.readlines()
    _, no_lons, no_lats, no_deps = lines[0].split()
    lons = np.array(list(map(float, lines[1].split())))
    lats = np.array(list(map(float, lines[2].split())))
    deps = np.array(list(map(float, lines[3].split())))
    # load vp model only
    model = np.zeros((int(no_lons), int(no_lats), int(no_deps)))
    for idep in range(int(no_deps)):
        for ilat in range(int(no_lats)):
            lon_data = lines[idep * int(no_lats) + ilat + 4].split()
            lon_data = list(map(float, lon_data))
            model[:, ilat, idep] = lon_data
    # construct xarray dataarray
    model = xr.DataArray(
        model,
        [("longitude", lons), ("latitude", lats), ("depth", deps)],
        name="vp",
    )
    return model


def load_ak135(parameter: str, copy_model: xr.DataArray) -> xr.DataArray:
    ak135 = np.loadtxt(
        resource(["model", "AK135F_AVG.csv"], normal_path=True), delimiter=","
    )
    h = ak135[:, 0]
    if parameter == "vp":
        v = ak135[:, 2]
    else:
        v = ak135[:, 3]
    f = interpolate.interp1d(h, v, bounds_error=False, fill_value="extrapolate")
    ak135_depth = f(copy_model.depth.data)
    ak135_abs_data = copy_model.copy()
    for index in range(len(copy_model.depth.data)):
        ak135_abs_data.data[:, :, index] = ak135_depth[index]
    return ak135_abs_data


def project_catalog(
    catalog: pd.DataFrame,
    startlon: float,
    startlat: float,
    endlon: float,
    endlat: float,
) -> pd.DataFrame:
    # change column names from lat,lon,dep to y,x,z
    catalog = catalog.rename(columns={"latitude": "y", "longitude": "x", "depth": "z"})
    catalog = catalog.reindex(columns=["x", "y", "z"])
    # project the catalog to the line
    res = pygmt.project(
        data=catalog,
        center=[startlon, startlat],
        endpoint=[endlon, endlat],
        convention="pz",
        unit=True,
        sort=True,
        length=[0, degrees2kilometers(LENGTH)],
        width=[0, 100],
    )

    # change column names back
    res.columns = ["dist", "dep"]
    # convert dist from km to degree
    res["dist"] = res["dist"].apply(lambda x: x / degrees2kilometers(1))
    return res


def plot_vertical_cross_section(
    fig: pygmt.Figure,
    ipanel: int,
    slab_model: xr.Dataset,
    model: xr.DataArray,
    df_events: pd.DataFrame,
):
    fig.basemap(
        projection="X?i/-?i",
        region=[0, LENGTH * degrees2kilometers(1), 0, 800],
    )
    startlon, startlat, endlon, endlat = LINES[ipanel]
    endlon, endlat = extend_line((startlon, startlat), (endlon, endlat), LENGTH)
    points = pygmt.project(
        center=[startlon, startlat],
        endpoint=[endlon, endlat],
        generate=0.02,
    )
    lons: np.ndarray = points.r
    lons[lons < 0] += 360
    lats: np.ndarray = points.s
    slab_deps = slab_interp(slab_model, lons, lats)
    lons -= 360
    cross_section = model_interp(model, lons, lats, np.linspace(0, 800, 801))
    cross_section_xarray = xr.DataArray(
        cross_section,
        dims=("h", "v"),
        coords={
            "h": np.linspace(0, LENGTH * degrees2kilometers(1), len(lons)),
            "v": np.linspace(0, 800, 801),
        },
    )

    fig.grdimage(cross_section_xarray.T)
    # plot slab2
    fig.plot(
        x=np.linspace(0, LENGTH * degrees2kilometers(1), len(slab_deps)),
        y=slab_deps,
        pen="1.5p,magenta",
    )
    # plot catalog
    catalog_line = project_catalog(df_events, startlon, startlat, endlon, endlat)
    fig.plot(
        x=catalog_line["dist"] * degrees2kilometers(1),
        y=catalog_line["dep"],
        style="c0.1c",
        pen="0.01c,black",
    )
    fig.text(
        position="BL",
        text=f"({ascii_lowercase[ipanel]})",
        font="18p,Helvetica-Bold,black",
    )
