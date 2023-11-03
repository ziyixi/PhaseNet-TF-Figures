import pygmt
from phasenettf.utils.slice import extend_line, slab_interp
from obspy.geodetics.base import degrees2kilometers
from phasenettf import resource, save_path
from pathlib import Path
import pandas as pd
import xarray as xr
import numpy as np
from string import ascii_lowercase

LENGTH = 7
LINE = [-180, -18, -174, -21]
DEP_START = 0
DEP_END = 700
LABELS = [
    "(a) Reference",
    "(b) Associated (Iteration #1)",
    "(c) Iteration #1",
    "(d) Iteration #2",
    "(e) Iteration #3",
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

    with fig.subplot(
        nrows=2,
        ncols=3,
        figsize=("9i", "6i"),
        sharex="b",
        sharey="l",
        margins=["0.2i", "0.2i"],
        frame=["WSen"],
        # autolabel="(a)",
    ):
        for idf, df in enumerate(catalog):
            with fig.set_panel(panel=idf):
                startlon, startlat, endlon, endlat = LINE
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

                catalog_line = project_catalog(df, startlon, startlat, endlon, endlat)
                fig.plot(
                    x=catalog_line["dist"] * degrees2kilometers(1),
                    y=catalog_line["dep"],
                    style="c0.08c",
                    pen="0.01c,black",
                )

                # plot slab2
                slab_deps = slab_interp(slab_model, lons, lats)
                fig.plot(
                    x=np.linspace(0, LENGTH * degrees2kilometers(1), len(slab_deps)),
                    y=slab_deps,
                    pen="1.5p,magenta",
                    label="Slab2",
                )

                # plot slab1
                slab1_deps = slab_interp(slab1_model, lons, lats)
                fig.plot(
                    x=np.linspace(0, LENGTH * degrees2kilometers(1), len(slab1_deps)),
                    y=slab1_deps,
                    pen="1.5p,red",
                    label="Slab 1.0",
                )

                fig.legend(
                    position="jBR+o0.2c/0.2c",
                    # box="+gwhite+p1p",
                    transparency=50,
                )
                # text at TR corner
                fig.text(
                    position="TL",
                    text=LABELS[idf],
                    font="12p,Helvetica,black",
                    offset="j0.1c/-0.5c",
                    no_clip=True,
                )

    # * plot map
    fig.shift_origin(xshift="6.6i", yshift="-0.3i")
    fig.basemap(
        region=[-184, -172, -26, -14],
        projection="M3i",
        frame=["WSen", "xaf+lLongitude (degree)", "yaf+lLatitude (degree)"],
    )
    plot_earth_relief(fig)
    plot_text(fig)
    plot_lines_on_map(fig)

    save_path(fig, Path(__file__).resolve().stem)


def load_catalog():
    # * load reference catalog
    reference_df = pd.read_csv(
        resource(["catalog", "tonga_catalog_updated_2023_0426.csv"], normal_path=True),
        usecols=[0, 1, 2],
    )
    # * load associated catalog
    associated_df = pd.read_csv(
        resource(["catalog", "continious_associated_catalog.csv"], normal_path=True),
        usecols=["longitude", "latitude", "z(km)", "event_index"],
        sep=r"\s+",
    )
    associated_df.rename(
        columns={
            "longitude": "longitude",
            "latitude": "latitude",
            "z(km)": "depth",
            "event_index": "event_index",
        },
        inplace=True,
    )
    associated_assignment_df = pd.read_csv(
        resource(["catalog", "continious_associated_assignment.csv"], normal_path=True),
        skiprows=1,
        names=[
            "id",
            "date",
            "time",
            "amp",
            "type",
            "prob",
            "event_index",
            "gamma_score",
        ],
        sep=r"\s+",
    )
    associated_assignment_df = associated_assignment_df[
        associated_assignment_df["gamma_score"] > 0
    ]
    # filter associated_df so event_index appears at least 10 times in associated_assignment_df
    associated_df = associated_df[
        associated_df["event_index"].isin(
            associated_assignment_df["event_index"]
            .value_counts()[associated_assignment_df["event_index"].value_counts() > 10]
            .index
        )
    ]
    # * load bootstrapped catalog
    bootstrapped_df = pd.read_csv(
        resource(["catalog", "continious_bootstrapped.csv"], normal_path=True),
        usecols=[1, 2, 3],
    )
    print(len(bootstrapped_df))
    # * load semi-supervised catalog
    semi_df = pd.read_csv(
        resource(["catalog", "continious_semi.csv"], normal_path=True),
        usecols=[1, 2, 3],
    )

    # * load semi-supervised iteration 3 catalog
    iter3_df = pd.read_csv(
        resource(["catalog", "continuous_v3.csv"], normal_path=True),
        usecols=[1, 2, 3],
    )

    catalogs = [reference_df, associated_df, bootstrapped_df, semi_df, iter3_df]
    return catalogs


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
        width=[0, 70],
    )

    # change column names back
    res.columns = ["dist", "dep", "event_id"]
    # convert dist from km to degree
    res["dist"] = res["dist"].apply(lambda x: x / degrees2kilometers(1))
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
        pen="1.6p,black",
        limit="-10000/-7000",
    )
    # plot only -1000m contour
    fig.grdcontour(
        grd_topo,
        interval=1000,
        pen="1p,black",
        limit="-1100/-1000",
    )
    # plot text on the top left corner
    fig.text(
        # x=-183.4,
        # y=-14.5,
        position="TL",
        text="(f)",
        font="18p,Helvetica,black",
        offset="j0.1c/0.1c",
    )
    # plot back-arc spreading centers
    # fig.plot(
    #     data=resource(["symbols", "lau_neovolcanic.xy"], normal_path=True),
    #     pen="2p,magenta",
    # )
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
        # {
        #     "text": "Tonga Ridge",
        #     "x": -174.5,
        #     "y": -21.2,
        #     "font": "12p,Helvetica-Bold,black",
        #     "angle": 67,
        # },
        {
            "text": "Lau Basin",
            "x": -176,
            "y": -18,
            "font": "12p,Helvetica-Bold,black",
            "angle": 65,
        },
        # {
        #     "text": "Lau Ridge",
        #     "x": -179,
        #     "y": -19,
        #     "font": "12p,Helvetica-Bold,black",
        #     "angle": 65,
        # },
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
    startlon, startlat, endlon, endlat = LINE
    endlon, endlat = extend_line((startlon, startlat), (endlon, endlat), LENGTH)
    fig.plot(
        x=[startlon, endlon],
        y=[startlat, endlat],
        pen="3p,blue",
    )
