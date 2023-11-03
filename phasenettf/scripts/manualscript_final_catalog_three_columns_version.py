import pygmt
from obspy.geodetics.base import degrees2kilometers
from phasenettf import resource, save_path
from pathlib import Path
import pandas as pd
import xarray as xr
from string import ascii_lowercase, ascii_uppercase
from phasenettf.utils.slice import extend_line, slab_interp
from phasenettf.utils.gcmt import gcmt_to_psmeca
import numpy as np

LENGTH = 7
REVERSE_LENGTH = 2
LINES = [
    [-179.5, -17, -173.5, -20],
    [-180, -18, -174, -21],
    [-180.5, -19, -174.5, -22],
    [-181, -20, -175, -23],
]
DEP_START = 0
DEP_END = 700

NEW_MECA_POSITION_LIST = [[179, -15.5], [178.5, -22], [-177, -15.5], [-175, -17]]
PLOT_GCMT_VERTICAL_LISTS = [
    ["030994E", "201808190019A"],
    ["200911091044A", "201809061549A"],
    [],
    [],
]
NEW_TEXTS = [
    ["09/03/1994", "19/08/2018"],
    ["09/11/2009", "06/09/2018"],
    [],
    [],
]
PLOT_GCMT_HORIZONTAL_LISTS = [
    "200911091044A",
    "201809061549A",
    "030994E",
    "201808190019A",
]
NEW_TEXTS_ALL = ["09/11/2009", "06/09/2018", "09/03/1994", "19/08/2018"]
TITLES = [
    "Final catalogue of this study@^(11/2009-12/2010)",
    "Manually picked reference catalogue@^(11/2009-12/2010)",
    "ISC-EHB catalogure@^(1964-2020)",
]


def main():
    catalog_list = [
        load_final_catalog(),
        load_reference_catalog(),
        load_ehb_catalog(),
    ]
    slab_model = xr.open_dataset(
        resource(["slab2", f"ker_slab2_depth.grd"], normal_path=True)
    )
    slab1_model = xr.open_dataset(
        resource(["slab2", f"ker_slab1.0_clip.grd"], normal_path=True)
    )

    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="18p",
        MAP_LABEL_OFFSET="12p",
        FONT_ANNOT_PRIMARY="16p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="4p",
        FONT_TITLE="15p,black",
        MAP_FRAME_PEN="1p,black",
    )

    # plot maps
    fig.shift_origin(xshift="0.5i", yshift="15.5i")
    with fig.subplot(
        nrows=1,
        ncols=3,
        subsize="3.4i",
        sharex="b",
        sharey="l",
        margins=["0.2i", "0.2i"],
        frame=["WSen"],
        autolabel="(a)+o0.1i",
    ):
        for ipanel in range(3):
            with fig.set_panel(ipanel):
                catalog = catalog_list[ipanel]
                fig.basemap(
                    region=[-184, -172, -26, -14],
                    projection="M?",
                    frame=[
                        f"+t{TITLES[ipanel]}",
                        "xaf+lLongitude (degree)",
                        "yaf+lLatitude (degree)",
                    ],
                )
                plot_earth_relief(fig)
                # plot events
                pygmt.makecpt(
                    cmap="jet", series=[0, 700, 1], continuous=True, reverse=True
                )
                fig.plot(
                    x=catalog["longitude"],
                    y=catalog["latitude"],
                    style="c0.08c",
                    fill=catalog["depth"],
                    cmap=True,
                )
                plot_text(fig)
                plot_lines_on_map(fig)
                plot_gcmt_events(fig)

    # cross-section
    fig.shift_origin(xshift="0i", yshift="-12.5i")
    with fig.subplot(
        nrows=4,
        ncols=3,
        figsize=("11.5i", "12i"),
        sharex="b",
        sharey="l",
        margins=["0.2i", "0.2i"],
        frame=["WSen"],
        autolabel="(d)",
    ):
        for imodel in range(3):
            catalog = catalog_list[imodel]
            for iline in range(4):
                ipanel = iline * 3 + imodel
                with fig.set_panel(ipanel):
                    fig.basemap(
                        projection="X?i/-?i",
                        frame=["xaf+lDistance (km)", "yaf+lDepth (km)"],
                        region=[
                            0,
                            (LENGTH + REVERSE_LENGTH) * degrees2kilometers(1),
                            DEP_START,
                            DEP_END,
                        ],
                    )

                    line = LINES[iline]
                    startlon, startlat, endlon, endlat = line
                    endlon, endlat = extend_line(
                        (startlon, startlat), (endlon, endlat), LENGTH
                    )
                    startlon, startlat = extend_line(
                        (startlon, startlat), (endlon, endlat), -REVERSE_LENGTH
                    )
                    points = pygmt.project(
                        center=[startlon, startlat],
                        endpoint=[endlon, endlat],
                        generate=0.02,
                    )
                    lons = points.r
                    lons[lons < 0] += 360
                    lats = points.s

                    catalog_line = project_catalog(
                        catalog, startlon, startlat, endlon, endlat
                    )
                    fig.plot(
                        x=catalog_line["dist"] * degrees2kilometers(1),
                        y=catalog_line["dep"],
                        style="c0.08c",
                        pen="0.01c,black",
                    )

                    # plot slab2
                    slab_deps = slab_interp(slab_model, lons, lats)
                    fig.plot(
                        x=np.linspace(
                            0,
                            (LENGTH + REVERSE_LENGTH) * degrees2kilometers(1),
                            len(slab_deps),
                        ),
                        y=slab_deps,
                        pen="1.5p,magenta",
                        label="Slab2",
                    )

                    # plot slab1
                    slab1_deps = slab_interp(slab1_model, lons, lats)
                    fig.plot(
                        x=np.linspace(
                            0,
                            (LENGTH + REVERSE_LENGTH) * degrees2kilometers(1),
                            len(slab1_deps),
                        ),
                        y=slab1_deps,
                        pen="1.5p,red",
                        label="Slab 1.0",
                    )

                    fig.legend(
                        position="jBR+o0.2c/0.2c",
                        # box="+gwhite+p1p",
                        transparency=50,
                    )

                    plot_gcmt_vertical(
                        fig,
                        startlon,
                        startlat,
                        endlon,
                        endlat,
                        select_list=PLOT_GCMT_VERTICAL_LISTS[iline],
                        newtext_list=NEW_TEXTS[iline],
                    )

                    # plot A-A', B-B', C-C', D-D'
                    fig.text(
                        x=0,
                        y=-50,
                        text=f"{ascii_uppercase[iline]}",
                        font="14p,Helvetica-Bold,black",
                        no_clip=True,
                    )
                    fig.text(
                        x=1000,
                        y=-50,
                        text=f"{ascii_uppercase[iline]}'",
                        font="14p,Helvetica-Bold,black",
                        no_clip=True,
                    )

    pygmt.makecpt(cmap="jet", series=[0, 700, 1], continuous=True, reverse=True)
    fig.colorbar(
        position="JBC+w4i/0.8c+h+o2.8i/2.5c",
        box=False,
        frame=["a0.5f", f'"+LDepth (km)"'],
        scale=1,
    )
    # save
    save_path(fig, Path(__file__).resolve().stem)


def load_final_catalog():
    data = pd.read_csv(
        # resource(["catalog", "continuous_v3.csv"], normal_path=True),
        resource(["catalog", "continious_semi.csv"], normal_path=True),
        usecols=[1, 2, 3],
    )
    return data


def load_ehb_catalog():
    data = pd.read_csv(
        # resource(["catalog", "continuous_v3.csv"], normal_path=True),
        resource(["catalog", "isc_ehb.csv"], normal_path=True),
        usecols=[5, 6, 7],
        skiprows=1,
        names=["latitude", "longitude", "depth"],
    )
    return data


def load_reference_catalog():
    data = pd.read_csv(
        # resource(["catalog", "continuous_v3.csv"], normal_path=True),
        resource(["catalog", "tonga_catalog_updated_2023_0426.csv"], normal_path=True),
        usecols=[0, 1, 2],
        names=["latitude", "longitude", "depth"],
        skiprows=1,
    )
    return data


def plot_earth_relief(fig: pygmt.Figure):
    grd_topo = pygmt.datasets.load_earth_relief(
        resolution="02m", region=[-184, -172, -26, -14], registration="gridline"
    )
    assert type(grd_topo) == xr.DataArray
    # plot 2000m contour of topography, start from -2000m to -10000m
    fig.grdcontour(
        grd_topo,
        interval=1000,
        pen="1p,gray",
        limit="-10000/-7000",
    )
    # plot only -1000m contour
    fig.grdcontour(
        grd_topo,
        interval=1000,
        pen="1.3p,gray",
        limit="-1100/-1000",
    )
    fig.coast(land="gray")


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
            "text": "Lau Basin",
            "x": -176,
            "y": -18,
            "font": "12p,Helvetica-Bold,black",
            "angle": 65,
        },
        {
            "text": "Fiji",
            "x": 177,
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
        startlon, startlat = extend_line(
            (startlon, startlat), (endlon, endlat), -REVERSE_LENGTH
        )
        fig.plot(
            x=[startlon, endlon],
            y=[startlat, endlat],
            pen="3p,black",
        )
        # text at the start of the line
        fig.text(
            x=startlon,
            y=startlat,
            text=f"{ascii_uppercase[iline]}",
            font="10p,Helvetica-Bold,black",
            justify="MR",
            offset="0.06c",
        )
        fig.text(
            x=endlon,
            y=endlat,
            text=f"{ascii_uppercase[iline]}'",
            font="10p,Helvetica-Bold,black",
            justify="ML",
            offset="0.06c",
        )


def load_gcmt(
    vertical=False, select_list=[], newtext_list=[], plot_longitude=[], plot_latitude=[]
):
    large_events_dir = resource(["gcmt"], normal_path=True)
    meca_file = gcmt_to_psmeca(
        large_events_dir,
        has_text=True,
        use_tmp_file=vertical,
        select_list=select_list,
        newtext_list=newtext_list,
        plot_longitude=plot_longitude,
        plot_latitude=plot_latitude,
    )
    return meca_file


def plot_gcmt_events(fig):
    meca_file = load_gcmt(
        vertical=True,
        select_list=PLOT_GCMT_HORIZONTAL_LISTS,
        newtext_list=NEW_TEXTS_ALL,
        plot_longitude=[item[0] for item in NEW_MECA_POSITION_LIST],
        plot_latitude=[item[1] for item in NEW_MECA_POSITION_LIST],
    )
    # change event_name with NEW_TEXTS_ALL, meca_file is a df with event_name as one column
    fig.meca(
        meca_file,
        convention="mt",
        scale="16p+m+f9p,black",
        pen="0.5p,black,solid",
        compressionfill="red",
        offset="+s0.2c",
        # event_name=event_name,
    )


def plot_gcmt_vertical(fig, lon1, lat1, lon2, lat2, select_list, newtext_list):
    meca_file = load_gcmt(
        vertical=True, select_list=select_list, newtext_list=newtext_list
    )
    with pygmt.clib.Session() as lib:
        lib.call_module(
            module="coupe",
            args=f"{meca_file} -Sd0.6+f8p,black -Aa{lon1}/{lat1}/{lon2}/{lat2}/90/0.1/0/700+f -Gred -Q -N -M",
        )


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
        length=[0, degrees2kilometers((LENGTH + REVERSE_LENGTH))],
        width=[0, 70],
    )

    # change column names back
    res.columns = ["dist", "dep", "event_id"]
    # convert dist from km to degree
    res["dist"] = res["dist"].apply(lambda x: x / degrees2kilometers(1))
    return res
