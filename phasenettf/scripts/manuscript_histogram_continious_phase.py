"""
manuscript_histogram_continious_phase.py:

This script generates the histogram of the continious result used in the manuscript, with the AK135 model comparision.
"""
from pathlib import Path

import pandas as pd
import pygmt

from phasenettf import resource, save_path
from typing import List
from obspy.geodetics.base import locations2degrees, degrees2kilometers

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
    associated_assignment_df = load_associated_assignment_df(is_semi=False)
    a_p, a_s = load_associated_phases(associated_assignment_df)
    r_p, r_s = load_relocated_phases(associated_assignment_df)
    b_p, b_s = load_bootstrapped_phases(associated_assignment_df)
    associated_assignment_df = load_associated_assignment_df(is_semi=True)
    s_p, s_s = load_semi_phases(associated_assignment_df)
    p_list = [a_p, r_p, b_p, s_p]
    s_list = [a_s, r_s, b_s, s_s]
    diff_lists = [p_list, s_list]

    with fig.subplot(
        nrows=2,
        ncols=4,
        subsize=("2i", "2i"),
        margins=["0.08i", "0.05i"],
        frame=["WSen", "xaf", "yaf"],
        sharex="b",
        sharey="l",
    ):
        for irow in range(2):
            for icol in range(4):
                with fig.set_panel(panel=(irow, icol)):
                    plot_time_diff(fig, diff_lists[irow][icol], irow, icol)

        save_path(fig, Path(__file__).resolve().stem)


def plot_time_diff(fig: pygmt.Figure, diff_list: List[float], irow: int, icol: int):
    fig.basemap(
        projection="X?i/?i",
        region=[-10, 10, 0, 70000] if irow == 0 else [-10, 10, 0, 15000],
        frame=[f"+t{TITLES[icol]}", "xaf+lTime Difference (s)", "yaf+lCount"]
        if irow == 0
        else ["xaf+lTime Difference (s)", "yaf+lCount"],
    )
    fig.histogram(
        data=diff_list,
        series=0.5,
        pen="1p,black",
        center=True,
    )
    phase_type = "P" if irow == 0 else "S"
    if icol == 3:
        fig.text(
            position="TR", text=f"{phase_type} Wave", font="14p,Helvetica-Bold,black"
        )
    fig.text(
        position="TL", text=f"({LABELS[irow][icol]})", font="14p,Helvetica-Bold,black"
    )


# *============================== helper functions ==============================*
def load_associated_assignment_df(is_semi=False):
    if not is_semi:
        associated_assignment_df = pd.read_csv(
            resource(
                ["catalog", "continious_associated_assignment.csv"], normal_path=True
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
    else:
        associated_assignment_df = pd.read_csv(
            resource(["catalog", "semi-v1.assignment.csv"], normal_path=True),
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
    return associated_assignment_df


def load_associated_phases(associated_assignment_df_raw: pd.DataFrame):
    associated_assignment_df = associated_assignment_df_raw.copy()
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

    # keep associated_assignment_df only with event_index in associated_df
    event_index_set = set(associated_df["event_index"].unique())

    associated_assignment_df = associated_assignment_df[
        associated_assignment_df["event_index"].isin(event_index_set)
    ]

    associated_assignment_df["sta"] = associated_assignment_df["id"].apply(
        lambda x: x.split(".")[1]
    )
    associated_assignment_df["time"] = pd.to_datetime(
        associated_assignment_df["date"] + " " + associated_assignment_df["time"]
    )
    associated_assignment_df.drop(
        columns=["date", "id", "amp", "prob", "gamma_score"], inplace=True
    )

    taup_df = pd.read_csv(
        resource(["catalog", "associated_df_taup.csv"], normal_path=True),
    )
    taup_df["p_arrival_time"] = pd.to_datetime(taup_df["p_arrival_time"])
    taup_df["s_arrival_time"] = pd.to_datetime(taup_df["s_arrival_time"])

    # generate p diff column
    associated_p = associated_assignment_df[associated_assignment_df["type"] == "p"]
    associated_p = associated_p.merge(taup_df, on=["sta", "event_index"])
    associated_p.dropna(inplace=True)
    associated_p["p_diff"] = (
        associated_p["time"] - associated_p["p_arrival_time"]
    ).dt.total_seconds()

    # generate s diff column
    associated_s = associated_assignment_df[associated_assignment_df["type"] == "s"]
    associated_s = associated_s.merge(taup_df, on=["sta", "event_index"])
    associated_s.dropna(inplace=True)
    associated_s["s_diff"] = (
        associated_s["time"] - associated_s["s_arrival_time"]
    ).dt.total_seconds()

    return associated_p["p_diff"].to_list(), associated_s["s_diff"].to_list()


def load_relocated_phases(associated_assignment_df_raw: pd.DataFrame):
    associated_assignment_df = associated_assignment_df_raw.copy()
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

    # keep associated_assignment_df only with event_index in relocated_df
    event_index_set = set(relocated_df["event_index"].unique())

    associated_assignment_df = associated_assignment_df[
        associated_assignment_df["event_index"].isin(event_index_set)
    ]

    associated_assignment_df["sta"] = associated_assignment_df["id"].apply(
        lambda x: x.split(".")[1]
    )
    associated_assignment_df["time"] = pd.to_datetime(
        associated_assignment_df["date"] + " " + associated_assignment_df["time"]
    )
    associated_assignment_df.drop(
        columns=["date", "id", "amp", "prob", "gamma_score"], inplace=True
    )

    taup_df = pd.read_csv(
        resource(["catalog", "relocated_df_taup.csv"], normal_path=True),
    )
    taup_df["p_arrival_time"] = pd.to_datetime(taup_df["p_arrival_time"])
    taup_df["s_arrival_time"] = pd.to_datetime(taup_df["s_arrival_time"])

    # generate p diff column
    associated_p = associated_assignment_df[associated_assignment_df["type"] == "p"]
    associated_p = associated_p.merge(taup_df, on=["sta", "event_index"])
    associated_p.dropna(inplace=True)
    associated_p["p_diff"] = (
        associated_p["time"] - associated_p["p_arrival_time"]
    ).dt.total_seconds()

    # generate s diff column
    associated_s = associated_assignment_df[associated_assignment_df["type"] == "s"]
    associated_s = associated_s.merge(taup_df, on=["sta", "event_index"])
    associated_s.dropna(inplace=True)
    associated_s["s_diff"] = (
        associated_s["time"] - associated_s["s_arrival_time"]
    ).dt.total_seconds()

    return associated_p["p_diff"].to_list(), associated_s["s_diff"].to_list()


def load_bootstrapped_phases(associated_assignment_df_raw: pd.DataFrame):
    associated_assignment_df = associated_assignment_df_raw.copy()
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

    # keep associated_assignment_df only with event_index in bootstrapped_df
    event_index_set = set(bootstrapped_df["event_index"].unique())

    associated_assignment_df = associated_assignment_df[
        associated_assignment_df["event_index"].isin(event_index_set)
    ]

    associated_assignment_df["sta"] = associated_assignment_df["id"].apply(
        lambda x: x.split(".")[1]
    )
    associated_assignment_df["time"] = pd.to_datetime(
        associated_assignment_df["date"] + " " + associated_assignment_df["time"]
    )
    associated_assignment_df.drop(
        columns=["date", "id", "amp", "prob", "gamma_score"], inplace=True
    )

    taup_df = pd.read_csv(
        resource(["catalog", "bootstrapped_df_taup.csv"], normal_path=True),
    )
    taup_df["p_arrival_time"] = pd.to_datetime(taup_df["p_arrival_time"])
    taup_df["s_arrival_time"] = pd.to_datetime(taup_df["s_arrival_time"])

    # generate p diff column
    associated_p = associated_assignment_df[associated_assignment_df["type"] == "p"]
    associated_p = associated_p.merge(taup_df, on=["sta", "event_index"])
    associated_p.dropna(inplace=True)
    associated_p["p_diff"] = (
        associated_p["time"] - associated_p["p_arrival_time"]
    ).dt.total_seconds()

    # generate s diff column
    associated_s = associated_assignment_df[associated_assignment_df["type"] == "s"]
    associated_s = associated_s.merge(taup_df, on=["sta", "event_index"])
    associated_s.dropna(inplace=True)
    associated_s["s_diff"] = (
        associated_s["time"] - associated_s["s_arrival_time"]
    ).dt.total_seconds()

    return associated_p["p_diff"].to_list(), associated_s["s_diff"].to_list()


def load_semi_phases(associated_assignment_df_raw: pd.DataFrame):
    associated_assignment_df = associated_assignment_df_raw.copy()
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

    # keep associated_assignment_df only with event_index in semi_df
    event_index_set = set(semi_df["event_index"].unique())

    associated_assignment_df = associated_assignment_df[
        associated_assignment_df["event_index"].isin(event_index_set)
    ]

    associated_assignment_df["sta"] = associated_assignment_df["id"].apply(
        lambda x: x.split(".")[1]
    )
    associated_assignment_df["time"] = pd.to_datetime(
        associated_assignment_df["date"] + " " + associated_assignment_df["time"]
    )
    associated_assignment_df.drop(
        columns=["date", "id", "amp", "prob", "gamma_score"], inplace=True
    )

    taup_df = pd.read_csv(
        resource(["catalog", "semi_df_taup.csv"], normal_path=True),
    )
    taup_df["p_arrival_time"] = pd.to_datetime(taup_df["p_arrival_time"])
    taup_df["s_arrival_time"] = pd.to_datetime(taup_df["s_arrival_time"])

    # generate p diff column
    associated_p = associated_assignment_df[associated_assignment_df["type"] == "p"]
    associated_p = associated_p.merge(taup_df, on=["sta", "event_index"])
    associated_p.dropna(inplace=True)
    associated_p["p_diff"] = (
        associated_p["time"] - associated_p["p_arrival_time"]
    ).dt.total_seconds()

    # generate s diff column
    associated_s = associated_assignment_df[associated_assignment_df["type"] == "s"]
    associated_s = associated_s.merge(taup_df, on=["sta", "event_index"])
    associated_s.dropna(inplace=True)
    associated_s["s_diff"] = (
        associated_s["time"] - associated_s["s_arrival_time"]
    ).dt.total_seconds()

    return associated_p["p_diff"].to_list(), associated_s["s_diff"].to_list()
