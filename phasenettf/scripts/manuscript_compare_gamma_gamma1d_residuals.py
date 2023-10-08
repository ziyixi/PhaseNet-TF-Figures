"""
This file is designed to plot the histogram to compare GaMMA and GaMMA-1D residuals.
"""
from pathlib import Path

import pandas as pd
import pygmt

from phasenettf import resource, save_path
from typing import List

TITLES = [
    "GaMMA-1D",
    "GaMMA",
]
LABELS = [
    ["a", "b"],
    ["c", "d"],
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
    hist_lists = [
        [
            load_gamma1d_time_diff(),
            load_gamma_time_diff(),
        ],
        [
            load_gamma1d_depth_diff(),
            load_gamma_depth_diff(),
        ],
    ]

    with fig.subplot(
        nrows=2,
        ncols=2,
        subsize=("2i", "2i"),
        margins=["0.03i", "0.2i"],
        frame=["WSen", "xaf", "yaf"],
        sharex=False,
        sharey="l",
    ):
        for irow in range(2):
            for icol in range(2):
                with fig.set_panel(panel=(irow, icol)):
                    if irow == 0:
                        plot_time_diff(fig, hist_lists[irow][icol], irow, icol)
                    elif irow == 1:
                        plot_depth_diff(fig, hist_lists[irow][icol], irow, icol)
                    else:
                        raise ValueError("Invalid irow value.")

        save_path(fig, Path(__file__).resolve().stem)


def plot_time_diff(fig, diff_list: List[float], irow: int, icol: int):
    fig.basemap(
        projection="X?i/?i",
        region=[-30.0, 30.0, 0, 400],
        frame=[f"+t{TITLES[icol]}", "xaf+lTime Residuals (s)", "yaf+lCount"],
    )
    fig.histogram(
        data=diff_list,
        series=1,
        pen="1p,black",
        center=True,
    )
    fig.text(
        position="TL", text=f"({LABELS[irow][icol]})", font="14p,Helvetica-Bold,black"
    )


def plot_depth_diff(fig, diff_list: List[float], irow: int, icol: int):
    fig.basemap(
        projection="X?i/?i",
        region=[-200.0, 200.0, 0, 250],
        frame=["x80f+lDepth Residuals (km)", "yaf+lCount"],
    )
    fig.histogram(
        data=diff_list,
        series=6,
        pen="1p,black",
        center=True,
    )
    fig.text(
        position="TL", text=f"({LABELS[irow][icol]})", font="14p,Helvetica-Bold,black"
    )


# *============================== helper functions ==============================*
def load_gamma_time_diff():
    data = pd.read_csv(
        resource(["catalog", "gamma", "gamma_origin_time_diff.csv"], normal_path=True),
    )
    return data["origin_time_diff"].to_list()


def load_gamma1d_time_diff():
    data = pd.read_csv(
        resource(
            ["catalog", "gamma", "gamma1d_origin_time_diff.csv"], normal_path=True
        ),
    )
    return data["origin_time_diff"].to_list()


def load_gamma_depth_diff():
    data = pd.read_csv(
        resource(["catalog", "gamma", "gamma_origin_depth_diff.csv"], normal_path=True),
    )
    return data["origin_depth_diff"].to_list()


def load_gamma1d_depth_diff():
    data = pd.read_csv(
        resource(
            ["catalog", "gamma", "gamma1d_origin_depth_diff.csv"], normal_path=True
        ),
    )
    return data["origin_depth_diff"].to_list()
