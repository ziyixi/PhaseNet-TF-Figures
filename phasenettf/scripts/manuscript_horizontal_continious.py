"""
horizontal_top100km.py:

This script is used to show the top 100 km events distribution for continious result.
"""
from pathlib import Path

import pandas as pd
import pygmt
import xarray as xr

from phasenettf import resource, save_path

# * global settings
BASE_XSHIFTS = [0.6, 4.3, -4.3, 4.3]
BASE_YSHIFTS = [6.0, 0, -4.3, 0]
FRAMES = [
    ["Wsen", "xaf+lLongitude (degree)", "yaf+lLatitude (degree)"],
    ["wsen", "xaf+lLongitude (degree)", "yaf+lLatitude (degree)"],
    ["WSen", "xaf+lLongitude (degree)", "yaf+lLatitude (degree)"],
    ["wSen", "xaf+lLongitude (degree)", "yaf+lLatitude (degree)"],
]
LABELS = [
    "(a) Reference",
    "(b) Associated (Iteration 1)",
    "(c) Relocated (Iteration 1)",
    "(d) Semi-supervised (Iteration 2)",
]


# * main function
def main():
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

    catalogs = load_catalog()

    for ipanel in range(len(BASE_XSHIFTS)):
        fig.shift_origin(
            xshift=f"{BASE_XSHIFTS[ipanel]}i", yshift=f"{BASE_YSHIFTS[ipanel]}i"
        )
        # * plot map
        fig.basemap(
            region=[-184, -172, -26, -14],
            projection="X4ig/4ig",
            frame=FRAMES[ipanel],
        )
        plot_earth_relief(fig, ipanel)
        plot_events_map(fig, catalogs[ipanel])

    pygmt.makecpt(cmap="jet", series=[0, 700, 1], continuous=True, reverse=False)
    fig.colorbar(
        position="JBC+w5i/0.8c+h+o-2i/1.8c",
        box=False,
        frame=["a0.5f", f'"+LDepth (km)"'],
        scale=1,
    )

    save_path(fig, Path(__file__).resolve().stem)


# * utility functions
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
    catalogs = [reference_df, associated_df, bootstrapped_df, semi_df]

    return catalogs


def plot_earth_relief(fig: pygmt.Figure, ipanel: int):
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
        text=LABELS[ipanel],
        font="18p,Helvetica-Bold,black",
    )
    # plot back-arc spreading centers
    fig.plot(
        data=resource(["symbols", "lau_neovolcanic.xy"], normal_path=True),
        pen="2p,magenta",
    )


def plot_events_map(fig: pygmt.Figure, df_events: pd.DataFrame):
    pygmt.makecpt(cmap="jet", series=[0, 700, 1], continuous=True, reverse=False)
    fig.plot(
        x=df_events["longitude"],
        y=df_events["latitude"],
        style="c0.1c",
        fill=df_events["depth"],
        cmap=True,
    )
    fig.text(
        position="BL",
        text=f"{len(df_events)} events",
    )
