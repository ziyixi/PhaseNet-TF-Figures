""" 
relocated_events_horizontal.py

Plot the horizontal distribution of relocated events.
"""
from string import ascii_lowercase

import pandas as pd
import pygmt

from phasenettf import resource, save_path
from phasenettf.utils.slice import extend_line

# * configs
REGION = [-184, -172, -24, -14]
TRENCH_FILE = resource(["trench", "trench.gmt"], normal_path=True, check=True)
CATALOG_FILE = resource(
    ["catalog", "tomoDD.all_months_threshold0.reloc"], normal_path=True, check=True
)
LINES = [
    [-179, -16, -173, -19],
    [-179.5, -17, -173.5, -20],
    [-180, -18, -174, -21],
    [-180.5, -19, -174.5, -22],
    [-181, -20, -175, -23],
]
LENGTH = 7


def load_catalog():
    data = pd.read_csv(
        CATALOG_FILE, sep="\s+", usecols=[1, 2, 3], names=["lat", "lon", "depth"]
    )
    data = data[data["depth"] > 0]
    return data


def plot_map(fig: pygmt.Figure):
    fig.basemap(region=REGION, projection="M8i", frame=["WSen", "xaf", "yaf"])
    fig.coast(shorelines="0.5p,black")
    fig.plot(data=TRENCH_FILE, pen="2p,red")

    # plot catalog
    catalog = load_catalog()
    pygmt.makecpt(cmap="jet", series=[0, 800, 1], continuous=True)
    fig.plot(
        x=catalog["lon"],
        y=catalog["lat"],
        style="c0.15c",
        fill=catalog["depth"],
        cmap=True,
    )
    fig.colorbar(
        frame=[
            "x+lDepth (km)",
        ]
    )

    # plot lines
    for iline, line in enumerate(LINES):
        startlon, startlat, endlon, endlat = line
        endlon, endlat = extend_line((startlon, startlat), (endlon, endlat), LENGTH)
        fig.plot(
            x=[startlon, endlon],
            y=[startlat, endlat],
            pen="3p,blue",
        )
        fig.text(
            x=endlon,
            y=endlat,
            text=ascii_lowercase[iline],
            font="15p,Helvetica-Bold,black",
            justify="LM",
        )

    # plot slab2
    fig.grdcontour(
        resource(["slab2", f"ker_slab2_depth.grd"]),
        interval=100,
        pen="1.5p,magenta",
    )

    # plot inset
    with fig.inset(position="jTL+w5c+o0.2c", margin=0, box="+p1.5p,gold"):
        fig.coast(
            region="g",
            projection=f"G{(REGION[0]+REGION[1])/2}/{(REGION[2]+REGION[3])/2}/?",
            land="brown",
            water="lightblue",
        )
        fig.plot(
            x=[REGION[0], REGION[1], REGION[1], REGION[0], REGION[0]],
            y=[REGION[2], REGION[2], REGION[3], REGION[3], REGION[2]],
            pen="1p,black",
        )


def main():
    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL="10p", MAP_LABEL_OFFSET="5p", FONT_ANNOT_PRIMARY="10p")

    plot_map(fig)

    save_path(fig, "relocated_all_months_threshold0_horizontal", suffix="pdf")
