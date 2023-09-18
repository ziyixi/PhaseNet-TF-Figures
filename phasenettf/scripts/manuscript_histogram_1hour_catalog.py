"""
manuscript_histogram_1hour_catalog.py:

This script generates the histogram of the 1-hour catalog used in the manuscript.
"""
from pathlib import Path

import pandas as pd
import pygmt

from phasenettf import resource, save_path
from typing import List
from obspy.geodetics.base import locations2degrees, degrees2kilometers

TITLES = [
    "Reference Catalog",
    "Associated Catalog",
    "Relocated Catalog",
    "Bootstrapped Catalog",
]
LABELS = [
    ["a", "b", "c", "d"],
    ["", "e", "f", "g"],
    ["", "h", "i", "j"],
    ["", "k", "l", "m"],
]


def main():
    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="10p",
        MAP_LABEL_OFFSET="12p",
        FONT_ANNOT_PRIMARY="12p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="12p",
        FONT_TITLE="14p,black",
        MAP_FRAME_PEN="1p,black",
    )
    catalog = load_catalog()

    with fig.subplot(
        nrows=4,
        ncols=4,
        subsize=("2i", "2i"),
        margins=["0.03i", "0.2i"],
        frame=["WSen", "xaf", "yaf"],
        sharey="l",
    ):
        for irow in range(4):
            for icol in range(4):
                with fig.set_panel(panel=(irow, icol)):
                    if irow == 0:
                        plot_depth_distribution(fig, catalog, icol)
                    elif irow == 1:
                        plot_depth_differece_distribution(fig, catalog, icol)
                    elif irow == 2:
                        plot_time_differece_distribution(fig, catalog, icol)
                    elif irow == 3:
                        plot_location_differece_distribution(fig, catalog, icol)
                    else:
                        fig.basemap(
                            projection="X?i/?i",
                            region=[0, 1, 0, 1],
                        )

        save_path(fig, Path(__file__).resolve().stem)


def plot_depth_distribution(
    fig: pygmt.Figure, catalog: List[pd.DataFrame], icatalog: int
):
    fig.basemap(
        projection="X?i/?i",
        region=[0, 700, 0, 300],
        frame=[f"+t{TITLES[icatalog]}", "xaf+lDepth (km)", "yaf+lCount"],
    )
    fig.histogram(
        data=catalog[icatalog]["depth"],
        series=20,
        pen="1p,black",
    )
    if icatalog == 0:
        fig.text(
            position="TR", text="Depth Distribution", font="10p,Helvetica-Bold,black"
        )
    fig.text(
        position="TL", text=f"({LABELS[0][icatalog]})", font="14p,Helvetica-Bold,black"
    )


def plot_depth_differece_distribution(
    fig: pygmt.Figure, catalog: List[pd.DataFrame], icatalog: int
):
    if icatalog == 0:
        return
    fig.basemap(
        projection="X?i/?i",
        region=[-100, 100, 0, 400],
        frame=["WSen", "xa80f+lDepth Difference (km)", "yaf+lCount"]
        if icatalog == 1
        else ["wSen", "xa80f+lDepth Difference (km)", "yaf"],
    )
    diff_list = generate_depth_diff(catalog[icatalog], catalog[0])
    fig.histogram(
        data=diff_list,
        series=10,
        pen="1p,black",
        center=True,
    )
    if icatalog == 1:
        fig.text(
            position="TR",
            text="Depth - Reference",
            font="10p,Helvetica-Bold,black",
        )
    fig.text(
        position="TL", text=f"({LABELS[1][icatalog]})", font="14p,Helvetica-Bold,black"
    )


def plot_time_differece_distribution(
    fig: pygmt.Figure, catalog: List[pd.DataFrame], icatalog: int
):
    if icatalog == 0:
        return
    fig.basemap(
        projection="X?i/?i",
        region=[-20, 20, 0, 500],
        frame=["WSen", "xa15f+lOrigin Time Difference (s)", "yaf+lCount"]
        if icatalog == 1
        else ["wSen", "xa15f+lOrigin Time Difference (s)", "yaf"],
    )
    diff_list = generate_time_diff(catalog[icatalog], catalog[0])
    fig.histogram(
        data=diff_list,
        series=2,
        pen="1p,black",
        center=True,
    )
    if icatalog == 1:
        fig.text(
            position="TR",
            text="Origin Time - Reference",
            font="10p,Helvetica-Bold,black",
        )
    fig.text(
        position="TL", text=f"({LABELS[2][icatalog]})", font="14p,Helvetica-Bold,black"
    )


def plot_location_differece_distribution(
    fig: pygmt.Figure, catalog: List[pd.DataFrame], icatalog: int
):
    if icatalog == 0:
        return
    fig.basemap(
        projection="X?i/?i",
        region=[0, 300, 0, 600],
        frame=["WSen", "xa200f+lGreat Circle Distance (km)", "yaf+lCount"]
        if icatalog == 1
        else ["wSen", "xa200f+lGreat Circle Distance (km)", "yaf"],
    )
    diff_list = generate_location_diff(catalog[icatalog], catalog[0])
    fig.histogram(
        data=diff_list,
        series=15,
        pen="1p,black",
    )
    if icatalog == 1:
        fig.text(
            position="TR",
            text="Distance to Reference",
            font="10p,Helvetica-Bold,black",
        )
    fig.text(
        position="TL", text=f"({LABELS[3][icatalog]})", font="14p,Helvetica-Bold,black"
    )


# *============================== helper functions ==============================*
def load_catalog():
    # * load reference catalog
    reference_df = pd.read_csv(
        resource(["catalog", "tonga_catalog_updated_2023_0426.csv"], normal_path=True),
        usecols=[0, 1, 2, 3],
    )
    reference_df["time"] = pd.to_datetime(reference_df["time"])
    # * load associated catalog
    associated_df = pd.read_csv(
        resource(["catalog", "1hour_window_associated_catalog.csv"], normal_path=True),
        usecols=["time", "longitude", "latitude", "z(km)", "event_index"],
        sep=r"\s+",
    )
    associated_df["time"] = pd.to_datetime(associated_df["time"])
    associated_df.rename(
        columns={
            "z(km)": "depth",
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
    associated_df.pop("event_index")

    # * load relocated catalog
    relocated_df = pd.read_csv(
        resource(["catalog", "1hour_window_relocated.csv"], normal_path=True),
        usecols=[1, 2, 3, 10, 11, 12, 13, 14, 15],
        names=[
            "latitude",
            "longitude",
            "depth",
            "year",
            "month",
            "day",
            "hour",
            "minute",
            "second",
        ],
        sep=r"\s+",
    )
    relocated_df["time"] = pd.to_datetime(
        relocated_df[["year", "month", "day", "hour", "minute", "second"]]
    )
    relocated_df.drop(
        columns=["year", "month", "day", "hour", "minute", "second"], inplace=True
    )
    # * load bootstrapped catalog
    bootstrapped_df = pd.read_csv(
        resource(["catalog", "1hour_window_refined.csv"], normal_path=True),
        usecols=[1, 2, 3, 4],
    )
    bootstrapped_df["time"] = pd.to_datetime(bootstrapped_df["datetime"])
    bootstrapped_df.drop(columns=["datetime"], inplace=True)
    catalogs = [reference_df, associated_df, relocated_df, bootstrapped_df]

    return catalogs


def generate_depth_diff(df1: pd.DataFrame, df2: pd.DataFrame):
    df1 = df1.sort_values("time")
    df2 = df2.sort_values("time")
    merged_df = pd.merge_asof(
        df1,
        df2,
        on="time",
        direction="nearest",
        tolerance=pd.Timedelta(seconds=30),
        suffixes=("_cur", "_ref"),
    )
    merged_df.dropna(inplace=True)
    diff_list = (merged_df["depth_cur"] - merged_df["depth_ref"]).tolist()
    return diff_list


def generate_time_diff(df1: pd.DataFrame, df2: pd.DataFrame):
    df1 = df1.copy().sort_values("time")
    df1 = df1.rename(columns={"time": "time_cur"})

    df2 = df2.copy().sort_values("time")
    df2 = df2.rename(columns={"time": "time_ref"})

    merged_df = pd.merge_asof(
        df1,
        df2,
        left_on="time_cur",
        right_on="time_ref",
        direction="nearest",
        tolerance=pd.Timedelta(seconds=30),
    )

    merged_df["time_cur"] = pd.to_datetime(merged_df["time_cur"])
    merged_df["time_ref"] = pd.to_datetime(merged_df["time_ref"])
    merged_df.dropna(subset=["time_cur", "time_ref"], inplace=True)
    merged_df["time_diff"] = (
        merged_df["time_cur"] - merged_df["time_ref"]
    ).dt.total_seconds()
    return merged_df["time_diff"].to_list()


def generate_location_diff(df1: pd.DataFrame, df2: pd.DataFrame):
    df1 = df1.sort_values("time")
    df2 = df2.sort_values("time")
    merged_df = pd.merge_asof(
        df1,
        df2,
        on="time",
        direction="nearest",
        tolerance=pd.Timedelta(seconds=30),
        suffixes=("_cur", "_ref"),
    )
    merged_df.dropna(inplace=True)
    distances = locations2degrees(
        lat1=merged_df["latitude_cur"].to_list(),
        long1=merged_df["longitude_cur"].to_list(),
        lat2=merged_df["latitude_ref"].to_list(),
        long2=merged_df["longitude_ref"].to_list(),
    )
    distance_list = [degrees2kilometers(distance) for distance in distances]
    return distance_list
