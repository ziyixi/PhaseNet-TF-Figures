""" 
mantle_wedge_events.py

Plot the station and time distribution of mantle wedge events.
"""
from string import ascii_lowercase

import pandas as pd
import pygmt

from phasenettf import resource, save_path
from phasenettf.utils.slice import extend_line
import json

# * configs
REGION = [-184, -172, -24, -14]

# get EVENT_TIME, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH from EVENT_ID
#     22695 -17.735397 -178.915802   306.753         -92317.4         133595.6         308927.5   1544.2   1728.6   1922.9 2010  7 23 22 59 49.140 0.0     0     0    18     4  0.000  0.050   1
# EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
#     22695,
#     -178.915802,
#     -17.735397,
#     306.753,
#     "2010-07-23T22:59:49.140000",
# )
#     14090 -18.934633 -178.654587   319.646         -65391.3           6782.0         320003.6   1906.8   2133.8   1314.7 2010  8 29  4 15 49.540 0.0     0     0    22    12  0.000  0.079   1
# EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
#     14090,
#     -178.654587,
#     -18.934633,
#     319.646,
#     "2010-08-29T04:15:49.540000",
# )
#     22012 -19.612202 -178.454498   228.193         -45900.1         -65693.6         228715.9   1569.7    760.9   1947.9 2010  4 11  9 44 51.860 0.0     0     0    42     0  0.000  0.049   1
# EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
#     22012,
#     -178.454498,
#     -19.612202,
#     228.193,
#     "2010-04-11T09:44:51.860000",
# )
#     26520 -20.337921 -178.594833   259.686         -59488.9        -142793.7         261643.9   1237.0   1919.7   1665.3 2010  7 24  1 16 26.520 0.0     0     0    30     6  0.000  0.087   1
# EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
#     26520,
#     -178.594833,
#     -20.337921,
#     259.686,
#     "2010-07-24T01:16:26.520000",
# )
#     13877 -20.440792 -178.405792   258.704         -40562.7        -153733.9         260772.4   1265.3   1974.0   1594.0 2010  6  5 19 57 11.880 0.0     0     0    30     6  0.000  0.078   1
# EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
#     13877,
#     -178.405792,
#     -20.440792,
#     258.704,
#     "2010-06-05T19:57:11.880000",
# )
#     15239 -22.057981 -176.433289    73.264         159581.8        -336672.7          84294.4   3189.1   1673.7   3024.2 2010  6 23 11  0 17.280 0.0     0     0    16     2  0.000  0.025   1
# EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
#     15239,
#     -176.433289,
#     -22.057981,
#     73.264,
#     "2010-06-23T11:00:17.280000",
# )
#      1325 -17.506905 -175.782303   702.935         209175.6         146371.5         708688.1   2532.4   2482.4   2119.4 2010  1 15  2  1 38.850 0.0     0     0    20     8  0.000  0.070   1
# EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
#     1325,
#     -175.782303,
#     -17.506905,
#     702.935,
#     "2010-01-15T02:01:38.850000",
# )
#      3611 -18.904078 -178.802872   717.857         -74939.9           9293.2         718361.1   4033.5   4190.8   1977.3 2010  9 17  6 26 27.240 0.0     0     0    22     2  0.000  0.066   1
# EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
#     3611,
#     -178.802872,
#     -18.904078,
#     717.857,
#     "2010-09-17T06:26:27.240000",
# )
#     10973 -17.837137 -178.546692   142.796         -56568.2         126309.5         144334.3   1310.4   1690.2   3317.9 2010 10  1  7 48 47.350 0.0     0     0    20    10  0.000  0.062   1
# EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
#     10973,
#     -178.546692,
#     -17.837137,
#     142.796,
#     "2010-10-01T07:48:47.350000",
# )
#      5123 -18.736052  179.945679   286.224        -206561.0          26825.2         289790.4   2068.0   2668.8   3178.4 2010  5 19  1 39 20.750 0.0     0     0   113     2  0.000  0.067   1
# EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
#     5123,
#     -180.054321,
#     -18.736052,
#     286.224,
#     "2010-05-19T01:39:20.750000",
# )
#     48872 -19.806345 -178.770203   222.565         -77758.0         -86696.6         223668.5   2204.3   3412.5   3487.4 2010  4 11  9 44 50.010 0.0     0     0   109     2  0.000  0.072   1
EVENT_ID, EVENT_LONGITUDE, EVENT_LATITUDE, EVENT_DEPTH, EVENT_TIME = (
    48872,
    -178.770203,
    -19.806345,
    222.565,
    "2010-04-11T09:44:50.010000",
)


def load_stations():
    stations = json.load(
        open(resource(["stations", "stations.json"], normal_path=True, check=True))
    )
    # eg: stations["YL.B01W..HH"]["longitude"]=-175
    # create a pd dataframe with columns: station, longitude, latitude. station should use B01W
    stations_df = []
    for station in stations:
        stations_df.append(
            {
                "station": station.split(".")[1],
                "longitude": stations[station]["longitude"],
                "latitude": stations[station]["latitude"],
            }
        )
    stations_df = pd.DataFrame(stations_df)
    return stations_df


def load_association():
    association = pd.read_csv(
        resource(
            ["association", "phasenet_tf_all_months.assignment.threshold05.csv"],
            normal_path=True,
            check=True,
        ),
        sep="\s+",
        # id	timestamp	amp	type	prob	event_index	gamma_score
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
        skiprows=1,
    )
    association["timestamp"] = association["date"] + "T" + association["time"]
    association = association[association["event_index"] == EVENT_ID]
    # get all unique station names, and the corresponding type==p arrival time
    # eg: station_name=B01W, arrival_time=2019-01-01T00:00:00.000000Z
    station_name_p = []
    arrival_time_p = []
    station_name_s = []
    arrival_time_s = []
    for _, row in association.iterrows():
        if row["type"] == "p":
            station_name_p.append(row["id"].split(".")[1])
            arrival_time_p.append(
                (
                    pd.Timestamp(row["timestamp"]) - pd.Timestamp(EVENT_TIME)
                ).total_seconds()
            )
        elif row["type"] == "s":
            station_name_s.append(row["id"].split(".")[1])
            arrival_time_s.append(
                (
                    pd.Timestamp(row["timestamp"]) - pd.Timestamp(EVENT_TIME)
                ).total_seconds()
            )
    return station_name_p, arrival_time_p, station_name_s, arrival_time_s


def plot_map(fig: pygmt.Figure):
    fig.basemap(region=REGION, projection="M8i", frame=["WSen", "xaf", "yaf"])
    fig.coast(shorelines="0.5p,black")
    fig.plot(
        data=resource(["trench", "trench.gmt"], normal_path=True, check=True),
        pen="2p,red",
    )

    # plot slab2
    fig.grdcontour(
        resource(["slab2", f"ker_slab2_depth.grd"]),
        interval=100,
        pen="1.5p,magenta",
    )

    # plot inset
    # with fig.inset(position="jTL+w5c+o0.2c", margin=0, box="+p1.5p,gold"):
    #     fig.coast(
    #         region="g",
    #         projection=f"G{(REGION[0]+REGION[1])/2}/{(REGION[2]+REGION[3])/2}/?",
    #         land="brown",
    #         water="lightblue",
    #     )
    #     fig.plot(
    #         x=[REGION[0], REGION[1], REGION[1], REGION[0], REGION[0]],
    #         y=[REGION[2], REGION[2], REGION[3], REGION[3], REGION[2]],
    #         pen="1p,black",
    #     )

    stations = load_stations()
    fig.plot(
        x=stations.longitude,
        y=stations.latitude,
        style="t0.3c",
        pen="1p,black",
        fill="black",
    )
    station_name_p, arrival_time_p, station_name_s, arrival_time_s = load_association()
    station_longitude_p = []
    station_latitude_p = []
    for station in station_name_p:
        station_longitude_p.append(
            stations[stations.station == station].iloc[0].longitude
        )
        station_latitude_p.append(
            stations[stations.station == station].iloc[0].latitude
        )
        if station_longitude_p[-1] > 0:
            station_longitude_p[-1] -= 360
    station_longitude_s = []
    station_latitude_s = []
    for station in station_name_s:
        station_longitude_s.append(
            stations[stations.station == station].iloc[0].longitude
        )
        station_latitude_s.append(
            stations[stations.station == station].iloc[0].latitude
        )
        if station_longitude_s[-1] > 0:
            station_longitude_s[-1] -= 360

    pygmt.makecpt(cmap="jet", series=[20, 100, 1], continuous=True)
    fig.plot(
        x=station_longitude_p,
        y=station_latitude_p,
        style="t0.3c",
        fill=arrival_time_p,
        cmap=True,
    )
    fig.text(
        x=station_longitude_p,
        y=station_latitude_p,
        text=station_name_p,
        font="5p,Helvetica,black",
        justify="CM",
        offset="0.2c",
    )
    # plot event location
    fig.plot(
        x=EVENT_LONGITUDE,
        y=EVENT_LATITUDE,
        style="a0.6c",
        pen="1p,black",
        fill="magenta",
    )
    fig.colorbar(
        frame=[
            "x+lSeconds (s)",
        ]
    )

    # text at top left for eventid, time, depth, magnitude
    fig.text(
        x=REGION[0] + 0.1,
        y=REGION[3] - 1,
        text=f"Event ID: {EVENT_ID}\n"
        + f"Time: {EVENT_TIME}\n"
        + f"Depth: {EVENT_DEPTH} km\n",
        font="15p,Helvetica,black",
        justify="BL",
    )


def plot_arrival(fig: pygmt.Figure):
    fig.basemap(
        region=[0, 150, -184, -172],
        projection="X6i/6i",
        frame=["WSen", "xaf+lTime (s)", "yaf+lLongitude (degree)"],
    )
    station_name_p, arrival_time_p, station_name_s, arrival_time_s = load_association()
    stations = load_stations()
    station_longitude_p = []
    for station in station_name_p:
        station_longitude_p.append(
            stations[stations.station == station].iloc[0].longitude
        )
        if station_longitude_p[-1] > 0:
            station_longitude_p[-1] -= 360
    station_longitude_s = []
    for station in station_name_s:
        station_longitude_s.append(
            stations[stations.station == station].iloc[0].longitude
        )
        if station_longitude_s[-1] > 0:
            station_longitude_s[-1] -= 360

    fig.plot(
        x=arrival_time_p,
        y=station_longitude_p,
        style="t0.3c",
        pen="1p,black",
        fill="black",
        label="P",
    )
    fig.text(
        x=arrival_time_p,
        y=station_longitude_p,
        text=station_name_p,
        font="5p,Helvetica,black",
        justify="CM",
        offset="0.3c/0.1c",
    )
    fig.plot(
        x=arrival_time_s,
        y=station_longitude_s,
        style="t0.3c",
        pen="1p,red",
        fill="red",
        label="S",
    )
    fig.text(
        x=arrival_time_s,
        y=station_longitude_s,
        text=station_name_s,
        font="5p,Helvetica,red",
        justify="CM",
        offset="0.3c/0.1c",
    )
    fig.plot(
        x=0,
        y=EVENT_LONGITUDE,
        style="a0.6c",
        pen="1p,black",
        fill="magenta",
        no_clip=True,
    )
    fig.legend(position="jTR+o0.2c", box="+gwhite+p1p")


def main():
    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL="10p", MAP_LABEL_OFFSET="5p", FONT_ANNOT_PRIMARY="10p")

    plot_map(fig)
    fig.shift_origin(xshift="r9i")
    plot_arrival(fig)

    save_path(fig, f"mantle_wedge_events_{EVENT_ID}", suffix="pdf")
