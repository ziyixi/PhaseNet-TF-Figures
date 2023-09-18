"""
manuscript_horizontal_benchmark_1hourwindow.py

This script is used to generate the figures for the manuscript. It is used to compare reference, associated, relocated, and bootstrapped events for 1 hour window phasenet-tf predictions.
"""
from pathlib import Path

import pandas as pd
import pygmt
import xarray as xr

from phasenettf import resource, save_path

# * global settings
BASE_XSHIFTS = [0.3, 6.3, -6.3, 6.3]
BASE_YSHIFTS = [6.3, 0, -6.0, 0]
INNER_MAP_XSHIFTS = [1.5, 0, 1.5, 0]
INNER_MAP_YSHIFTS = [0, 0, 1.5, 1.5]
INNER_LON_XSHIFTS = [1.5, 0, 1.5, 0]
INNER_LON_YSHIFTS = [4.264, 4.264, 0, 0]
INNER_LAT_XSHIFTS = [0, 4.264, 0, 4.264]
INNER_LAT_YSHIFTS = [0, 0, 1.5, 1.5]
FRAMES = [
    ["wSEn", "xaf+lLongitude (degree)", "yaf+lLatitude (degree)"],
    ["wSen", "xaf+lLongitude (degree)", "yaf"],
    ["wsEn", "xaf", "yaf"],
    ["wsen", "xaf", "yaf"],
]
LABELS = ["(a) Reference", "(b) Associated", "(c) Relocated", "(d) Bootstrapped"]


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
        fig.shift_origin(
            xshift=f"{INNER_MAP_XSHIFTS[ipanel]}i",
            yshift=f"{INNER_MAP_YSHIFTS[ipanel]}i",
        )
        fig.basemap(
            region=[-184, -172, -26, -14],
            projection="X4ig/4ig",
            frame=FRAMES[ipanel],
        )
        plot_earth_relief(fig, ipanel)
        plot_events_map(fig, catalogs[ipanel])
        fig.shift_origin(
            xshift=f"{-INNER_MAP_XSHIFTS[ipanel]}i",
            yshift=f"{-INNER_MAP_YSHIFTS[ipanel]}i",
        )

        # * plot lon
        fig.shift_origin(
            xshift=f"{INNER_LON_XSHIFTS[ipanel]}i",
            yshift=f"{INNER_LON_YSHIFTS[ipanel]}i",
        )
        if ipanel == 0:
            frame = ["Wsen", "xaf", "ya400f+lDepth (km)"]
        elif ipanel == 2:
            frame = ["Wsen", "xaf", "ya400f"]
        else:
            frame = ["wsen", "xaf", "ya400f"]
        fig.basemap(region=[-184, -172, 0, 800], projection="X4ig/-1.236i", frame=frame)
        plot_events_lon(fig, catalogs[ipanel])
        fig.shift_origin(
            xshift=f"{-INNER_LON_XSHIFTS[ipanel]}i",
            yshift=f"{-INNER_LON_YSHIFTS[ipanel]}i",
        )

        # * plot lat
        fig.shift_origin(
            xshift=f"{INNER_LAT_XSHIFTS[ipanel]}i",
            yshift=f"{INNER_LAT_YSHIFTS[ipanel]}i",
        )
        fig.basemap(
            region=[0, 800, -26, -14],
            projection="X1.236i/4ig",
            frame=["wSen", "xa400f", "yaf"]
            if ipanel < 2
            else ["wsen", "xa400f", "yaf"],
        )
        plot_events_lat(fig, catalogs[ipanel])
        fig.shift_origin(
            xshift=f"{-INNER_LAT_XSHIFTS[ipanel]}i",
            yshift=f"{-INNER_LAT_YSHIFTS[ipanel]}i",
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
        resource(["catalog", "1hour_window_associated_catalog.csv"], normal_path=True),
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
        resource(
            ["catalog", "1hour_window_associated_assignment.csv"], normal_path=True
        ),
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

    # * load relocated catalog
    relocated_df = pd.read_csv(
        resource(["catalog", "1hour_window_relocated.csv"], normal_path=True),
        usecols=[1, 2, 3],
        names=["latitude", "longitude", "depth"],
        sep=r"\s+",
    )
    # * load bootstrapped catalog
    bootstrapped_df = pd.read_csv(
        resource(["catalog", "1hour_window_refined.csv"], normal_path=True),
        usecols=[1, 2, 3],
    )
    catalogs = [reference_df, associated_df, relocated_df, bootstrapped_df]
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
    # plot text on the top left corner
    fig.text(
        # x=-183.4,
        # y=-14.5,
        position="TL",
        text=LABELS[ipanel],
        font="18p,Helvetica-Bold,black",
    )


def plot_events_map(fig: pygmt.Figure, df_events: pd.DataFrame):
    pygmt.makecpt(cmap="jet", series=[0, 800, 1], continuous=True, reverse=False)
    fig.plot(
        x=df_events["longitude"],
        y=df_events["latitude"],
        style="c0.12c",
        fill=df_events["depth"],
        cmap=True,
    )
    fig.text(
        position="BL",
        text=f"{len(df_events)} events",
    )


def plot_events_lat(fig: pygmt.Figure, df_events: pd.DataFrame):
    pygmt.makecpt(cmap="jet", series=[0, 800, 1], continuous=True, reverse=False)
    fig.plot(
        x=df_events["depth"],
        y=df_events["latitude"],
        style="c0.08c",
        fill=df_events["depth"],
        cmap=True,
    )


def plot_events_lon(fig: pygmt.Figure, df_events: pd.DataFrame):
    pygmt.makecpt(cmap="jet", series=[0, 800, 1], continuous=True, reverse=False)
    fig.plot(
        x=df_events["longitude"],
        y=df_events["depth"],
        style="c0.08c",
        fill=df_events["depth"],
        cmap=True,
    )
