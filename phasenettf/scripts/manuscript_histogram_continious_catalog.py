"""
manuscript_histogram_continious_catalog.py:

This script generates the histogram of the continious result used in the manuscript, with the catalog time and depth distribution.
"""
from pathlib import Path

import pandas as pd
import pygmt

from phasenettf import resource, save_path
from typing import List

TITLES = [
    "Associated",
    "Relocated",
    "Bootstraped",
    "Semi-supervised",
]
LABELS = [
    ["a", "b", "c", "d"],
    ["d", "e", "f", "g"],
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
    catalogs = [
        load_associated_phases(),
        load_relocated_phases(),
        load_bootstrapped_phases(),
        load_semi_phases(),
    ]

    with fig.subplot(
        nrows=2,
        ncols=4,
        subsize=("2i", "2i"),
        margins=["0.08i", "0.2i"],
        frame=["WSen", "xaf", "yaf"],
        sharey="l",
    ):
        for irow in range(2):
            for icol in range(4):
                with fig.set_panel(panel=(irow, icol)):
                    if irow == 0:
                        plot_depth_distribution(fig, catalogs[icol], icol)
                    elif irow == 1:
                        plot_time_distribution(fig, catalogs[icol], icol)
                    else:
                        fig.basemap(
                            projection="X?i/?i",
                            region=[0, 1, 0, 1],
                        )

        save_path(fig, Path(__file__).resolve().stem)


def plot_depth_distribution(fig: pygmt.Figure, catalog: pd.DataFrame, icol: int):
    fig.basemap(
        projection="X?i/?i",
        region=[0, 700, 0, 3000],
        frame=[f"+t{TITLES[icol]}", "xaf+lDepth (km)", "yaf+lCount"],
    )
    fig.histogram(
        data=catalog["depth"],
        series=20,
        pen="1p,black",
    )


def plot_time_distribution(fig: pygmt.Figure, catalog: pd.DataFrame, icol: int):
    with pygmt.config(
        FORMAT_DATE_MAP="o",
        FORMAT_TIME_PRIMARY_MAP="Character",
    ):
        fig.basemap(
            projection="X?i/?i",
            region="2009-10-1T/2011-2-1T/0/2000",
            frame=["WSen", "pxa2Of1o", "sxa1Y", "yaf+lCount"]
            if icol == 0
            else ["wSen", "pxa2Of1o", "sxa1Y", "yaf+lCount"],
        )
        # format catalpg time to like 2009-10-1
        # catalog["time"] = catalog["time"].dt.strftime("%Y-%m-%d")
        fig.histogram(
            data=catalog["time"],
            series=f"{86400*30}",
            pen="1p,black",
        )


# *============================== helper functions ==============================*
def load_associated_phases():
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
    associated_df = pd.read_csv(
        resource(["catalog", "continious_associated_catalog.csv"], normal_path=True),
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
    associated_df = associated_df[
        associated_df["event_index"].isin(
            associated_assignment_df["event_index"]
            .value_counts()[associated_assignment_df["event_index"].value_counts() > 10]
            .index
        )
    ]

    return associated_df


def load_relocated_phases():
    relocated_df = pd.read_csv(
        resource(["catalog", "continious_relocated.csv"], normal_path=True),
        usecols=[0, 1, 2, 3, 10, 11, 12, 13, 14, 15],
        names=[
            "event_index",
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

    return relocated_df


def load_bootstrapped_phases():
    bootstrapped_df = pd.read_csv(
        resource(["catalog", "continious_bootstrapped.csv"], normal_path=True),
        usecols=[0, 1, 2, 3, 4],
    )
    bootstrapped_df.rename(
        columns={
            "id": "event_index",
        },
        inplace=True,
    )
    bootstrapped_df["time"] = pd.to_datetime(bootstrapped_df["datetime"])
    bootstrapped_df.drop(columns=["datetime"], inplace=True)

    return bootstrapped_df


def load_semi_phases():
    semi_df = pd.read_csv(
        resource(["catalog", "continious_semi.csv"], normal_path=True),
        usecols=[0, 1, 2, 3, 4],
    )
    semi_df.rename(
        columns={
            "id": "event_index",
        },
        inplace=True,
    )
    semi_df["time"] = pd.to_datetime(semi_df["datetime"])
    semi_df.drop(columns=["datetime"], inplace=True)

    return semi_df
