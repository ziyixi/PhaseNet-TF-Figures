import pandas as pd

from phasenettf import resource, save_path
import pygmt
from pathlib import Path

TIME_START = pd.Timestamp("2010-09-15T13:43:00")
TIME_END = pd.Timestamp("2010-09-15T13:53:00")


def main():
    df, df_catalog = load_associated_phases(TIME_START, TIME_END)

    # for each unique event_index, assign a unique color in GMT
    colors_list = [
        "red3",
        "magenta3",
        "cyan3",
        "yellow3",
        "purple3",
        "green3",
        "blue3",
        "orange3",
        "gray3",
        "brown3",
        # more colors below
        "lightred",
        "lightgreen",
        "lightblue",
        "lightyellow",
        "lightpurple",
        "lightorange",
        "lightcyan",
        "lightmagenta",
        "lightgray",
        "lightbrown",
    ]
    colors = {}
    for i, event_index in enumerate(df["event_index"].unique()):
        colors[event_index] = colors_list[i]

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
    # * latitude plot
    fig.basemap(
        region=[TIME_START, TIME_END, -26, -13],
        projection="X8i/2i",
        frame=[
            "WSen",
            f'xa2Mf+l{str(TIME_START).split(" ")[0]}',
            "yaf+lLatitude (degree)",
        ],
    )
    for event_index in df["event_index"].unique():
        if event_index == 34127:
            continue
        df_event = df[df["event_index"] == event_index]
        df_event_p = df_event[df_event["type"] == "p"]
        df_event_s = df_event[df_event["type"] == "s"]
        if len(df_event_p) > 0:
            fig.plot(
                x=df_event_p.time,
                y=df_event_p.lat,
                style="c0.15c",
                fill=colors[event_index],
            )
        if len(df_event_s) > 0:
            fig.plot(
                x=df_event_s.time,
                y=df_event_s.lat,
                style="t0.2c",
                fill=colors[event_index],
            )
        # df_catalog has columns: time, longitude, latitude, depth, event_index
        # plot the catalog time and position
        df_catalog_event = df_catalog[df_catalog["event_index"] == event_index]
        fig.plot(
            x=df_catalog_event["time"],
            y=df_catalog_event["latitude"],
            style=f"+2c",
            fill=colors[event_index],
            pen=f"1p,{colors[event_index]}",
        )
        fig.plot(x="1900-01-01", y=0, style="c0.15c", fill="black", label="P")
        fig.plot(x="1900-01-01", y=0, style="t0.2c", fill="black", label="S")
        fig.plot(x="1900-01-01", y=0, style="+0.3c", fill="black", label="Event")
        fig.legend(position="JTR+jTR+o0.2c", box="+gwhite+p1p")

    # * longitude plot
    fig.shift_origin(yshift="h+0.2i")
    fig.basemap(
        region=[TIME_START, TIME_END, 171, 191],
        projection="X8i/2i",
        frame=["Wsen", "xa2Mf", "yaf+lLongitude (degree)"],
    )
    for event_index in df["event_index"].unique():
        if event_index == 34127:
            continue
        df_event = df[df["event_index"] == event_index]
        df_event_p = df_event[df_event["type"] == "p"]
        df_event_s = df_event[df_event["type"] == "s"]
        if len(df_event_p) > 0:
            fig.plot(
                x=df_event_p.time,
                y=df_event_p.lon,
                style="c0.15c",
                fill=colors[event_index],
            )
        if len(df_event_s) > 0:
            fig.plot(
                x=df_event_s.time,
                y=df_event_s.lon,
                style="t0.2c",
                fill=colors[event_index],
            )
        # df_catalog has columns: time, longitude, latitude, depth, event_index
        # plot the catalog time and position
        df_catalog_event = df_catalog[df_catalog["event_index"] == event_index]
        fig.plot(
            x=df_catalog_event["time"],
            y=df_catalog_event["longitude"],
            style=f"+2c",
            fill=colors[event_index],
            pen=f"1p,{colors[event_index]}",
        )
        fig.plot(x="1900-01-01", y=0, style="c0.15c", fill="black", label="P")
        fig.plot(x="1900-01-01", y=0, style="t0.2c", fill="black", label="S")
        fig.plot(x="1900-01-01", y=0, style="+0.3c", fill="black", label="Event")
        fig.legend(position="JTR+jTR+o0.2c", box="+gwhite+p1p")

    save_path(fig, Path(__file__).resolve().stem)


def load_associated_phases(time_start: pd.Timestamp, time_end: pd.Timestamp):
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

    # associated_df = associated_df[
    #     associated_df["event_index"].isin(
    #         associated_assignment_df["event_index"]
    #         .value_counts()[associated_assignment_df["event_index"].value_counts() > 10]
    #         .index
    #     )
    # ]
    # filter associated_assignment_df by  associated_df
    # associated_assignment_df = associated_assignment_df[
    #     associated_assignment_df["event_index"].isin(associated_df["event_index"])
    # ]

    associated_assignment_df["station"] = associated_assignment_df["id"].apply(
        lambda x: x.split(".")[1]
    )

    stations_df = pd.read_csv(
        resource(["stations", "stations.csv"], normal_path=True), sep=r"\s+"
    )
    associated_assignment_df = associated_assignment_df.merge(stations_df, on="station")
    associated_assignment_df["time"] = pd.to_datetime(
        [
            f"{d}T{t}"
            for d, t in zip(
                associated_assignment_df["date"], associated_assignment_df["time"]
            )
        ]
    )
    associated_assignment_df.pop("date")

    # Filter by time for associated_assignment_df
    associated_assignment_df = associated_assignment_df[
        (associated_assignment_df["time"] >= time_start)
        & (associated_assignment_df["time"] <= time_end)
    ]
    # convert associated_assignment_df["lon"] to the range of [0, 360]
    associated_assignment_df["lon"] = associated_assignment_df["lon"].apply(
        lambda x: x if x >= 0 else x + 360
    )
    # convert associated_df["lon"] to the range of [0, 360]
    associated_df["longitude"] = associated_df["longitude"].apply(
        lambda x: x if x >= 0 else x + 360
    )

    return associated_assignment_df, associated_df
