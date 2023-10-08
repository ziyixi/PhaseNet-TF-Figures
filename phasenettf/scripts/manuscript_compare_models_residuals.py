"""
Compare the residuals of the each iteration of the semi-supervised learning PhaseNet-TF, retrained PhaseNet, and PhaseNet
"""
from pathlib import Path

import pandas as pd
import pygmt

from phasenettf import resource, save_path
from typing import List

TITLES = [
    "PhaseNet-TF@^(iteration 1)",
    "PhaseNet-TF@^(iteration 2)",
    "PhaseNet-TF@^(iteration 3)",
    "PhaseNet@^Tonga model",
    "PhaseNet@^origional model",
]
LABELS = [
    ["a", "b", "c", "d", "e"],
    ["f", "g", "h", "i", "j"],
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

    p_wave_time_diff_list = [
        get_phase_time_diff_list(load_ptf_iter1("p"), load_ref_phases("p")),
        get_phase_time_diff_list(load_ptf_iter2("p"), load_ref_phases("p")),
        get_phase_time_diff_list(load_ptf_iter3("p"), load_ref_phases("p")),
        get_phase_time_diff_list(load_phasenet_retrained("p"), load_ref_phases("p")),
        load_phasenet_origional("p"),
    ]
    s_wave_time_diff_list = [
        get_phase_time_diff_list(load_ptf_iter1("s"), load_ref_phases("s")),
        get_phase_time_diff_list(load_ptf_iter2("s"), load_ref_phases("s")),
        get_phase_time_diff_list(load_ptf_iter3("s"), load_ref_phases("s")),
        get_phase_time_diff_list(load_phasenet_retrained("s"), load_ref_phases("s")),
        load_phasenet_origional("s"),
    ]
    load_phasenet_origional("p")

    with fig.subplot(
        nrows=2,
        ncols=5,
        subsize=("2i", "2i"),
        margins=["0.03i", "0.05i"],
        frame=["WSen", "xaf", "yaf"],
        sharex="b",
        sharey="l",
    ):
        for irow in range(2):
            for icol in range(5):
                with fig.set_panel(panel=(irow, icol)):
                    if irow == 0:
                        plot_time_diff(
                            fig, p_wave_time_diff_list[icol], irow, icol, "P"
                        )
                    elif irow == 1:
                        plot_time_diff(
                            fig, s_wave_time_diff_list[icol], irow, icol, "S"
                        )
                    else:
                        fig.basemap(
                            projection="X?i/?i",
                            region=[0, 1, 0, 1],
                        )

        save_path(fig, Path(__file__).resolve().stem)


def plot_time_diff(fig, diff_list: List[float], irow: int, icol: int, phase_type: str):
    fig.basemap(
        projection="X?i/?i",
        region=[-1.0, 1.0, 0, 10000] if phase_type == "P" else [-1.0, 1.0, 0, 2500],
        frame=[f"+t{TITLES[icol]}", "xa0.4f+lTime Difference (s)", "yaf+lCount"]
        if phase_type == "P"
        else ["xa0.4f+lTime residuals (s)", "yaf+lCount"],
    )
    fig.histogram(
        data=diff_list,
        series=0.05,
        pen="1p,black",
        center=True,
    )
    if icol == 4:
        fig.text(
            position="TR", text=f"{phase_type} Wave", font="14p,Helvetica-Bold,black"
        )
    fig.text(
        position="TL", text=f"({LABELS[irow][icol]})", font="14p,Helvetica-Bold,black"
    )


# *============================== helper functions ==============================*
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


def load_ptf_iter1(phase_type: str):
    phases = pd.read_csv(
        resource(["catalog", "phasenettf_v0_phases_1hourwindow.csv"], normal_path=True)
    )
    phases = phases[["time", "sta", "phase"]]
    phases.rename(columns={"phase": "type"}, inplace=True)
    # phase to lower case
    phases["type"] = phases["type"].str.lower()
    # remove last character of time
    phases["time"] = phases["time"].str[:-1]
    # time to datetime
    phases["time"] = pd.to_datetime(phases["time"])
    return select_phases(phases, phase_type)


def load_ptf_iter2(phase_type: str):
    phases = pd.read_csv(
        resource(["catalog", "phasenettf_v1_phases_continious.csv"], normal_path=True)
    )
    phases = phases[["time", "sta", "phase"]]
    phases.rename(columns={"phase": "type"}, inplace=True)
    # phase to lower case
    phases["type"] = phases["type"].str.lower()
    # remove last character of time
    phases["time"] = phases["time"].str[:-1]
    # time to datetime
    phases["time"] = pd.to_datetime(phases["time"])
    return select_phases(phases, phase_type)


def load_ptf_iter3(phase_type: str):
    phases = pd.read_csv(
        resource(["catalog", "phasenettf_v2_phases_continious.csv"], normal_path=True)
    )
    phases = phases[["time", "sta", "phase"]]
    phases.rename(columns={"phase": "type"}, inplace=True)
    # phase to lower case
    phases["type"] = phases["type"].str.lower()
    # remove last character of time
    phases["time"] = phases["time"].str[:-1]
    # time to datetime
    phases["time"] = pd.to_datetime(phases["time"])
    return select_phases(phases, phase_type)


def load_phasenet_origional(phase_type: str):
    if phase_type == "p":
        phases = pd.read_csv(
            resource(["catalog", "p_diff_phasenet_origional.csv"], normal_path=True),
            names=["index", "diff"],
        )
        res = phases["diff"].to_list()
    elif phase_type == "s":
        phases = pd.read_csv(
            resource(["catalog", "s_diff_phasenet_origional.csv"], normal_path=True),
            names=["index", "diff"],
        )
        res = phases["diff"].to_list()
    return res


def load_phasenet_retrained(phase_type: str):
    phases = pd.read_csv(
        resource(
            ["catalog", "phasenet_retrained_phases_1hourwindow.csv"], normal_path=True
        )
    )
    phases = phases[["time", "sta", "phase"]]
    phases.rename(columns={"phase": "type"}, inplace=True)
    # phase to lower case
    phases["type"] = phases["type"].str.lower()
    # remove last character of time
    phases["time"] = phases["time"].str[:-1]
    # time to datetime
    phases["time"] = pd.to_datetime(phases["time"])
    return select_phases(phases, phase_type)


def select_phases(phases: pd.DataFrame, phase_type: str):
    cur = phases[phases["type"] == phase_type.lower()]
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
    # only keep abs(time_diff) <= 0.5
    result = result[abs(result["time_diff"]) <= 1.0]
    return result["time_diff"].to_list()
