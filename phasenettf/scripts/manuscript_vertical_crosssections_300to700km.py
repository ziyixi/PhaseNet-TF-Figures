import pygmt
from phasenettf.utils.slice import extend_line, slab_interp
from obspy.geodetics.base import degrees2kilometers
from phasenettf import resource, save_path
from pathlib import Path
import pandas as pd
import xarray as xr
import numpy as np
from phasenettf.utils.seismicity import fit_slab_upper_interface_from_seismicity
from string import ascii_lowercase
from phasenettf.utils.gcmt import gcmt_to_psmeca

LENGTH = 7
LINES = [
    [-179, -16, -173, -19],
    [-179.5, -17, -173.5, -20],
    [-180, -18, -174, -21],
    [-180.5, -19, -174.5, -22],
    [-181, -20, -175, -23],
    [-183, -17.8, -177, -17.8],
]
DEP_START = 300
DEP_END = 700
NEW_MECA_POSITION_LIST = [[179, -15.5], [178.5, -22], [-177, -15.5], [-175, -17]]
PLOT_GCMT_VERTICAL_LISTS = [
    [],
    ["030994E", "201808190019"],
    [],
    [],
    [],
    ["200911091044A", "201809061549A", "030994E", "201808190019A"],
]


def main():
    # * prepare data
    catalog = load_catalog()
    slab_model = xr.open_dataset(
        resource(["slab2", f"ker_slab2_depth.grd"], normal_path=True)
    )
    slab1_model = xr.open_dataset(
        resource(["slab2", f"ker_slab1.0_clip.grd"], normal_path=True)
    )
    upper_interface = fit_slab_upper_interface_from_seismicity(catalog, threshold=0.05)

    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="18p",
        MAP_LABEL_OFFSET="12p",
        FONT_ANNOT_PRIMARY="16p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="12p",
        FONT_TITLE="18p,black",
        MAP_FRAME_PEN="1p,black",
    )

    fig.shift_origin(yshift="4.5i")
    with fig.subplot(
        nrows=2,
        ncols=3,
        subsize=("3i", f"{3*4/7:.2f}i"),
        sharex="b",
        sharey="l",
        margins=["0.2i", "0.2i"],
        frame=["WSen"],
        autolabel="(a)",
    ):
        for iline, line in enumerate(LINES):
            with fig.set_panel(panel=iline):
                startlon, startlat, endlon, endlat = line
                endlon, endlat = extend_line(
                    (startlon, startlat), (endlon, endlat), LENGTH
                )
                points = pygmt.project(
                    center=[startlon, startlat],
                    endpoint=[endlon, endlat],
                    generate=0.02,
                )
                lons = points.r
                lons[lons < 0] += 360
                lats = points.s

                fig.basemap(
                    projection="X?i/-?i",
                    frame=["xaf+lDistance (km)", "yaf+lDepth (km)"],
                    region=[0, LENGTH * degrees2kilometers(1), DEP_START, DEP_END],
                )

                catalog_line = project_catalog(
                    catalog, startlon, startlat, endlon, endlat
                )
                fig.plot(
                    x=catalog_line["dist"] * degrees2kilometers(1),
                    y=catalog_line["dep"],
                    style="c0.08c",
                    pen="0.01c,black",
                )

                # plot gcmt events
                gcmt_line = project_gcmt(
                    load_gcmt(), startlon, startlat, endlon, endlat
                )
                # select gcmt_line, so event_id is in PLOT_GCMT_VERTICAL_LISTS[iline]
                gcmt_line = gcmt_line[
                    gcmt_line["event_id"].isin(PLOT_GCMT_VERTICAL_LISTS[iline])
                ]
                if len(gcmt_line) > 0:
                    fig.plot(
                        x=gcmt_line["dist"] * degrees2kilometers(1),
                        y=gcmt_line["dep"],
                        style="c0.2c",
                        pen="0.01c,black",
                        fill="red",
                    )

                # plot slab2
                slab_deps = slab_interp(slab_model, lons, lats)
                fig.plot(
                    x=np.linspace(0, LENGTH * degrees2kilometers(1), len(slab_deps)),
                    y=slab_deps,
                    pen="1.5p,magenta",
                )

                # plot slab1
                slab1_deps = slab_interp(slab1_model, lons, lats)
                fig.plot(
                    x=np.linspace(0, LENGTH * degrees2kilometers(1), len(slab1_deps)),
                    y=slab1_deps,
                    pen="1.5p,red",
                )

                # plot fitted upper interface
                fitted_deps = upper_interface(lons.to_numpy(), lats.to_numpy())
                fig.plot(
                    x=np.linspace(0, LENGTH * degrees2kilometers(1), len(slab_deps)),
                    y=fitted_deps,
                    pen="1.5p,blue",
                )

    # plot map
    fig.shift_origin(yshift="-4i", xshift="6i")
    fig.basemap(
        region=[-184, -172, -26, -14],
        projection="M3i",
        frame=["WSen", "xaf+lLongitude (degree)", "yaf+lLatitude (degree)"],
    )
    plot_earth_relief(fig)
    # plot_text(fig)
    plot_lines_on_map(fig)
    plot_gcmt_events(fig)

    save_path(fig, Path(__file__).resolve().stem)


def load_catalog():
    data = pd.read_csv(
        resource(["catalog", "continious_semi.csv"], normal_path=True),
        usecols=[1, 2, 3],
    )
    return data


def project_catalog(
    catalog: pd.DataFrame,
    startlon: float,
    startlat: float,
    endlon: float,
    endlat: float,
) -> pd.DataFrame:
    # change column names from lat,lon,dep to y,x,z
    catalog = catalog.rename(columns={"latitude": "y", "longitude": "x", "depth": "z"})
    catalog = catalog.reindex(columns=["x", "y", "z", "event_id"])
    # project the catalog to the line
    res = pygmt.project(
        data=catalog,
        center=[startlon, startlat],
        endpoint=[endlon, endlat],
        convention="pz",
        unit=True,
        sort=True,
        length=[0, degrees2kilometers(LENGTH)],
        width=[0, 50],
    )

    # change column names back
    res.columns = ["dist", "dep", "event_id"]
    # convert dist from km to degree
    res["dist"] = res["dist"].apply(lambda x: x / degrees2kilometers(1))
    return res


def project_gcmt(
    catalog: pd.DataFrame,
    startlon: float,
    startlat: float,
    endlon: float,
    endlat: float,
) -> pd.DataFrame:
    # change column names from lat,lon,dep to y,x,z
    catalog = catalog.rename(columns={"latitude": "y", "longitude": "x", "depth": "z"})
    catalog = catalog.reindex(columns=["x", "y", "z", "event_name"])
    catalog["x"] = catalog["x"].apply(lambda x: x - 360 if x > 0 else x)
    # print(catalog)
    # project the catalog to the line
    try:
        # print(startlon, startlat, endlon, endlat)
        res = pygmt.project(
            data=catalog,
            center=[startlon, startlat],
            endpoint=[endlon, endlat],
            convention="pz",
            unit=True,
            sort=True,
            # length=[0, degrees2kilometers(LENGTH)],
            # width=[0, 300],
        )
    except pd.errors.EmptyDataError:
        return pd.DataFrame(columns=["dist", "dep", "event_id"])

    # change column names back
    res.columns = ["dist", "dep", "event_id"]
    # convert dist from km to degree
    res["dist"] = res["dist"].apply(lambda x: x / degrees2kilometers(1))
    res["dep"] = res["dep"].apply(lambda x: x / 1000)
    return res


def plot_earth_relief(fig: pygmt.Figure):
    grd_topo = pygmt.datasets.load_earth_relief(
        resolution="02m", region=[-184, -172, -26, -14], registration="gridline"
    )
    assert type(grd_topo) == xr.DataArray
    # plot 2000m contour of topography, start from -2000m to -10000m
    fig.grdcontour(
        grd_topo,
        interval=1000,
        pen="0.3p,black",
        limit="-10000/-7000",
    )
    # plot only -1000m contour
    fig.grdcontour(
        grd_topo,
        interval=1000,
        pen="0.1p,black",
        limit="-1100/-1000",
    )
    # plot text on the top left corner
    fig.text(
        # x=-183.4,
        # y=-14.5,
        position="TL",
        text="f",
        font="18p,Helvetica-Bold,black",
    )
    # plot back-arc spreading centers
    fig.plot(
        data=resource(["symbols", "lau_neovolcanic.xy"], normal_path=True),
        pen="2p,magenta",
    )


def plot_text(fig: pygmt.Figure):
    text_elements = [
        {
            "text": "Tonga Trench",
            "x": -173.3,
            "y": -22.5,
            "font": "12p,Helvetica-Bold,black",
            "angle": 65,
        },
        {
            "text": "Tonga Ridge",
            "x": -174.5,
            "y": -21.2,
            "font": "12p,Helvetica-Bold,black",
            "angle": 67,
        },
        {
            "text": "Lau Basin",
            "x": -176,
            "y": -18,
            "font": "12p,Helvetica-Bold,black",
            "angle": 65,
        },
        {
            "text": "Lau Ridge",
            "x": -179,
            "y": -19,
            "font": "12p,Helvetica-Bold,black",
            "angle": 65,
        },
        {
            "text": "Fiji",
            "x": 178,
            "y": -17,
            "font": "12p,Helvetica-Bold,black",
            "angle": 0,
        },
    ]

    for element in text_elements:
        fig.text(**element)


def plot_lines_on_map(fig: pygmt.Figure):
    for iline, line in enumerate(LINES):
        startlon, startlat, endlon, endlat = line
        endlon, endlat = extend_line((startlon, startlat), (endlon, endlat), LENGTH)
        fig.plot(
            x=[startlon, endlon],
            y=[startlat, endlat],
            pen="2p,blue",
        )
        # text at the start of the line
        if iline < 5:
            fig.text(
                x=endlon,
                y=endlat,
                text=f"{ascii_lowercase[iline]})",
                font="16p,Helvetica-Bold,blue",
                justify="ML",
                offset="0.1c",
                no_clip=True,
            )
        else:
            fig.text(
                x=startlon,
                y=startlat,
                text=f"{ascii_lowercase[iline]})",
                font="16p,Helvetica-Bold,blue",
                justify="MR",
                offset="0.1c",
                no_clip=True,
            )


def load_gcmt():
    large_events_dir = resource(["gcmt"], normal_path=True)
    meca_file = gcmt_to_psmeca(large_events_dir, has_text=True)
    return meca_file


def plot_gcmt_events(fig):
    meca_file = load_gcmt()
    fig.meca(
        meca_file,
        convention="mt",
        scale="16p+m+f11p,red",
        pen="0.5p,black,solid",
        compressionfill="red",
        plot_longitude=[item[0] for item in NEW_MECA_POSITION_LIST],
        plot_latitude=[item[1] for item in NEW_MECA_POSITION_LIST],
        offset="+s0.2c",
    )
