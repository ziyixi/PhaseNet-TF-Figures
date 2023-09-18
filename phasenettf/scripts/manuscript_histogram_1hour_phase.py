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
    "PhaseNet-TF Predictions",
    "Associated Phases",
    "Relocated Phases",
    "Bootstrapped Phases",
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
    catalog = load_catalog()
    cur_phases = load_cur_phases()
    cur_phases_no_filter_minus_one = load_cur_phases(filter_minus_one=False)
    cur_phases_p_list = [select_phases(df, cur_phases, "P") for df in catalog]
    cur_phases_p_list = [
        cur_phases_no_filter_minus_one[cur_phases_no_filter_minus_one["type"] == "p"]
    ] + cur_phases_p_list
    cur_phases_s_list = [select_phases(df, cur_phases, "S") for df in catalog]
    cur_phases_s_list = [
        cur_phases_no_filter_minus_one[cur_phases_no_filter_minus_one["type"] == "s"]
    ] + cur_phases_s_list

    ref_phases_p = load_ref_phases("p")
    ref_phases_s = load_ref_phases("s")
    get_phase_time_diff_list(cur_phases_p_list[0], ref_phases_p)

    with fig.subplot(
        nrows=2,
        ncols=4,
        subsize=("2i", "2i"),
        margins=["0.03i", "0.05i"],
        frame=["WSen", "xaf", "yaf"],
        sharex="b",
        sharey="l",
    ):
        for irow in range(2):
            for icol in range(4):
                with fig.set_panel(panel=(irow, icol)):
                    if irow == 0:
                        diff_list = get_phase_time_diff_list(
                            cur_phases_p_list[icol], ref_phases_p
                        )
                        plot_time_diff(fig, diff_list, irow, icol, "P")
                    elif irow == 1:
                        diff_list = get_phase_time_diff_list(
                            cur_phases_s_list[icol], ref_phases_s
                        )
                        plot_time_diff(fig, diff_list, irow, icol, "S")
                    else:
                        fig.basemap(
                            projection="X?i/?i",
                            region=[0, 1, 0, 1],
                        )

        save_path(fig, Path(__file__).resolve().stem)


def plot_time_diff(fig, diff_list: List[float], irow: int, icol: int, phase_type: str):
    fig.basemap(
        projection="X?i/?i",
        region=[-0.5, 0.5, 0, 4000],
        frame=[f"+t{TITLES[icol]}", "xa0.4f+lTime Difference (s)", "yaf+lCount"]
        if phase_type == "P"
        else ["xa0.4f+lTime Difference (s)", "yaf+lCount"],
    )
    fig.histogram(
        data=diff_list,
        series=0.02,
        pen="1p,black",
        center=True,
    )
    if icol == 3:
        fig.text(
            position="TR", text=f"{phase_type} Wave", font="14p,Helvetica-Bold,black"
        )
    fig.text(
        position="TL", text=f"({LABELS[irow][icol]})", font="14p,Helvetica-Bold,black"
    )


# *============================== helper functions ==============================*
def load_catalog():
    # * load reference catalog
    # reference_df = pd.read_csv(
    #     resource(["catalog", "tonga_catalog_updated_2023_0426.csv"], normal_path=True),
    #     usecols=[0, 1, 2, 3],
    # )
    # reference_df["time"] = pd.to_datetime(reference_df["time"])
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

    # * load relocated catalog
    relocated_df = pd.read_csv(
        resource(["catalog", "1hour_window_relocated.csv"], normal_path=True),
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
    # * load bootstrapped catalog
    bootstrapped_df = pd.read_csv(
        resource(["catalog", "1hour_window_refined.csv"], normal_path=True),
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
    catalogs = [associated_df, relocated_df, bootstrapped_df]

    return catalogs


def load_cur_phases(filter_minus_one=True):
    phases = pd.read_csv(
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
    # id column to sta, eg: Z1.PNGI..BH to PNGI
    phases["sta"] = phases["id"].apply(lambda x: x.split(".")[1])

    phases = phases[phases["gamma_score"] > 0]
    if filter_minus_one:
        phases = phases[phases["event_index"] != -1]
    phases["timestamp"] = pd.to_datetime(phases["date"] + " " + phases["time"])
    phases.drop(columns=["date", "time"], inplace=True)
    phases.rename(
        columns={
            "timestamp": "time",
        },
        inplace=True,
    )
    return phases


def load_ref_phases(phase_type: str):
    phases = pd.read_csv(
        resource(["catalog", "tonga_picks_updated_2023_0426.csv"], normal_path=True),
    )
    phases.rename(columns={"stacode": "sta"}, inplace=True)
    if phase_type.lower() == "p":
        phases["time"] = pd.to_datetime(phases["arrival_time_P"])
    elif phase_type.lower() == "s":
        phases["time"] = pd.to_datetime(phases["arrival_time_S"])
    # drop rows when time is NaT
    phases.dropna(subset=["time"], inplace=True)
    return phases


def select_phases(df: pd.DataFrame, phases: pd.DataFrame, phase_type: str):
    cur = phases[phases["type"] == phase_type.lower()]
    # keep only the phases that event_index are in the df
    event_index_set = set(df["event_index"].unique())
    cur = cur[cur["event_index"].isin(event_index_set)]
    return cur


def get_phase_time_diff_list(df_cur: pd.DataFrame, df_ref: pd.DataFrame):
    df_cur, df_ref = df_cur.copy(), df_ref.copy()
    df_cur["time"] = pd.to_datetime(df_cur["time"])
    df_cur.rename(columns={"time": "time_cur"}, inplace=True)
    df_ref["time"] = pd.to_datetime(df_ref["time"])
    df_ref.rename(columns={"time": "time_ref"}, inplace=True)

    # sort dataframes by 'sta' and 'time'
    df_cur = df_cur.sort_values(["sta", "time_cur"])
    df_ref = df_ref.sort_values(["sta", "time_ref"])

    # create an empty dataframe to store the result
    result = pd.DataFrame()

    # find the nearest 'time' in df_cur for each 'sta'
    for sta in df_ref["sta"].unique():
        subset_cur = df_cur[df_cur["sta"] == sta]
        subset_ref = df_ref[df_ref["sta"] == sta]
        merged = pd.merge_asof(
            subset_ref,
            subset_cur,
            direction="nearest",
            left_on="time_ref",
            right_on="time_cur",
        )
        result = pd.concat([result, merged])

    result.dropna(subset=["time_ref", "time_cur"], inplace=True)
    result["time_diff"] = (result["time_cur"] - result["time_ref"]).dt.total_seconds()
    return result["time_diff"].to_list()
