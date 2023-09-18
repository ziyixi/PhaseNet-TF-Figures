"""
manuscript_event_station_distribution.py:

This script is used to generate the event-station distribution figure for the manuscript.
"""
import pygmt
from pathlib import Path
from phasenettf import resource, save_path
import xarray as xr
import json
import pandas as pd

# *  =========================Main function=========================
REGION = [-184, -172, -26, -14]


def main():
    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="14p",
        MAP_LABEL_OFFSET="8p",
        FONT_ANNOT_PRIMARY="12p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="8p",
        FONT_TITLE="14p,black",
        MAP_FRAME_PEN="1p,black",
    )

    fig.basemap(
        region=REGION,
        projection="M6i",
        frame=["WSen", "xaf", "yaf"],
    )

    plot_earth_relief(fig)
    plot_text(fig)
    plot_stations(fig)
    plot_events(fig)
    plot_inset(fig)

    save_path(fig, Path(__file__).resolve().stem)


# * =========================Helper functions=========================
def plot_earth_relief(fig: pygmt.Figure):
    grd_topo = pygmt.datasets.load_earth_relief(
        resolution="02m", region=REGION, registration="gridline"
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


def plot_text(fig: pygmt.Figure):
    text_elements = [
        {
            "text": "Tonga Trench",
            "x": -173.3,
            "y": -22.5,
            "font": "18p,Helvetica-Bold,black",
            "angle": 65,
        },
        {
            "text": "Tonga Ridge",
            "x": -174.5,
            "y": -21.2,
            "font": "18p,Helvetica-Bold,black",
            "angle": 67,
        },
        {
            "text": "Lau Basin",
            "x": -176,
            "y": -18,
            "font": "18p,Helvetica-Bold,black",
            "angle": 65,
        },
        {
            "text": "Lau Ridge",
            "x": -179,
            "y": -19,
            "font": "18p,Helvetica-Bold,black",
            "angle": 65,
        },
        {
            "text": "Fiji",
            "x": 178,
            "y": -17,
            "font": "18p,Helvetica-Bold,black",
            "angle": 0,
        },
    ]

    for element in text_elements:
        fig.text(**element)


def plot_stations(fig: pygmt.Figure):
    with open(resource(["stations", "stations.json"], normal_path=True)) as f:
        stations = json.load(f)
    # json's key is station name, value is a dict, with keys longitude, latitude, and local_depth_m
    # prepare a pandas dataframe for plotting
    station_df = []
    for station in stations:
        net, sta, _, _ = station.split(".")
        station_df.append(
            [
                net,
                sta,
                stations[station]["longitude"],
                stations[station]["latitude"],
            ]
        )
    station_df = pd.DataFrame(station_df, columns=["net", "sta", "lon", "lat"])
    # plot net==YL as reverse triangle, net==Z1 as triange, II as diamond
    fig.plot(
        x=station_df[station_df["net"] == "YL"]["lon"],
        y=station_df[station_df["net"] == "YL"]["lat"],
        style="i0.2c",
        pen="2p,black",
        label="YL",
    )
    fig.plot(
        x=station_df[station_df["net"] == "Z1"]["lon"],
        y=station_df[station_df["net"] == "Z1"]["lat"],
        style="t0.3c",
        pen="2p,black",
        label="Z1",
    )
    fig.plot(
        x=station_df[station_df["net"] == "II"]["lon"],
        y=station_df[station_df["net"] == "II"]["lat"],
        style="s0.3c",
        pen="2p,black",
        label="Z1",
    )
    # plot legend in the bottom right corner
    fig.legend(
        position="JBR+jBR+o0.2c/0.2c",
        box="+gwhite+p1p,black",
    )


def plot_inset(fig: pygmt.Figure):
    with fig.inset(position="jBL+w7c+o0.2c", margin=0):
        fig.coast(
            region="g",
            projection=f"W{(REGION[0]+REGION[1])//2}/6c",
            land="gray",
            water="white",
            frame=["wsen", "xafg", "yafg"],
        )
        fig.plot(
            x=[REGION[0], REGION[1], REGION[1], REGION[0], REGION[0]],
            y=[REGION[2], REGION[2], REGION[3], REGION[3], REGION[2]],
            pen="1p,black",
            projection=f"W{(REGION[0]+REGION[1])//2}/6c",
        )


def plot_events(fig: pygmt.Figure):
    df_events = pd.read_csv(
        resource(["catalog", "tonga_catalog_updated_2023_0426.csv"], normal_path=True),
        parse_dates=["time"],
    )
    pygmt.makecpt(cmap="oslo", series=[0, 700, 1], continuous=True, reverse=True)
    fig.plot(
        x=df_events["longitude"],
        y=df_events["latitude"],
        style="c0.15c",
        fill=df_events["depth"],
        cmap=True,
    )
    # plot the colorbar to the right of the map
    fig.colorbar(position="JMR+o0.2c/0c+w5i/0.3i", frame=["x+lDepth (km)"])