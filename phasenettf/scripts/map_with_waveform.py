""" 
map_with_waveform.py

This script is used to show the map of the station distribution of our study region, and waveform example of land/OBS stations, and their picks.
"""
import pygmt
from phasenettf import resource, save_path
import obspy
import numpy as np

# * configs
REGION = [-184, -172, -24, -14]
STATION_XMLS = resource(
    ["tongaml_continious_xml", "*xml"], normal_path=True, check=False
)
PHASES = {
    "A03": {
        "P": obspy.UTCDateTime("2009-12-02T21:03:10.404850"),
        "S": obspy.UTCDateTime("2009-12-02T21:04:21.269680"),
    },
    "C06": {
        "P": obspy.UTCDateTime("2009-12-02T21:03:06.937420"),
        "S": obspy.UTCDateTime("2009-12-02T21:04:15.712560"),
    },
    "TNGA": {
        "P": obspy.UTCDateTime("2009-12-02T21:03:09.509700"),
        "S": obspy.UTCDateTime("2009-12-02T21:04:20.529030"),
        "PS": obspy.UTCDateTime("2009-12-02T21:03:21.836071"),
    },
}


# * functions
def prepare_stations():
    invs = obspy.read_inventory(STATION_XMLS)
    lons, lats = [], []
    lons_marker, lats_marker, names_marker = [], [], []
    for inv in invs:
        lons.append(inv[0].longitude)
        lats.append(inv[0].latitude)
        if inv[0].code in PHASES:
            lons_marker.append(inv[0].longitude)
            lats_marker.append(inv[0].latitude)
            names_marker.append(inv[0].code)
    return lons, lats, lons_marker, lats_marker, names_marker


def load_waveform(sta: str):
    waveform_path = resource(
        ["waveforms", f"*{sta}*waveform.mseed"], normal_path=True, check=False
    )
    prediction_path = resource(
        ["waveforms", f"*{sta}*prediction.mseed"], normal_path=True, check=False
    )
    st = obspy.read(waveform_path)
    st_prediction = obspy.read(prediction_path)
    ptime = PHASES[sta]["P"]
    stime = PHASES[sta]["S"]
    start, end = ptime - 10, ptime + 110
    st.trim(starttime=start, endtime=end)
    st_prediction.trim(starttime=start, endtime=end)
    # detrend, demean, filter to 1HZ to 10HZ
    st.detrend("linear")
    st.detrend("demean")
    st.filter("bandpass", freqmin=1, freqmax=10, corners=4, zerophase=True)
    # get data
    res = []
    res_prediction = []
    for tr in st:
        res.append(tr.data)
    for tr in st_prediction:
        res_prediction.append(tr.data)
    x = np.linspace(0, 120, len(res[0]))
    return x, res, res_prediction, 10, stime - ptime + 10


def plot_map(fig: pygmt.Figure):
    fig.basemap(region=REGION, projection="M8i", frame=["WSen", "xaf", "yaf"])
    fig.coast(shorelines="0.5p,black")
    with fig.inset(position="jTR+w5c+o0.2c", margin=0, box="+p1.5p,gold"):
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

    sta_lons, sta_lats, lons_marker, lats_marker, names_marker = prepare_stations()
    fig.plot(
        x=sta_lons,
        y=sta_lats,
        style="t0.3c",
        pen="1p,black",
        fill="black",
    )
    # plot event
    fig.plot(
        x=-178.309,
        y=-17.7343,
        style="a0.8c",
        pen="2p,red",
    )
    # marker stations
    fig.plot(
        x=lons_marker,
        y=lats_marker,
        style="t0.3c",
        pen="1p,blue",
        fill="blue",
    )
    # arrows to show stations
    style = "v0.1i+s+e+a40+gorange+h0+p1p,orange"
    fig.plot(
        x=[lons_marker[0] - 0.1],
        y=[lats_marker[0] - 0.1],
        style=style,
        pen="0.04i,orange",
        direction=[[180], [-25.5]],
        no_clip=True,
    )
    fig.plot(
        x=[lons_marker[1] + 0.1],
        y=[lats_marker[1] + 0.1],
        style=style,
        pen="0.04i,orange",
        direction=[[-171], [-19]],
        no_clip=True,
    )
    fig.plot(
        x=[lons_marker[2] + 0.1],
        y=[lats_marker[2] - 0.1],
        style=style,
        pen="0.04i,orange",
        direction=[[-170.5], [-24.5]],
        no_clip=True,
    )
    # text event information
    fig.text(
        position="TL",
        text="ISC ID: 14226097",
        font="10p,Helvetica-Bold,black",
        offset="j0.05i/0.05i",
    )
    fig.text(
        position="TL",
        text="Mw: 4.9",
        font="10p,Helvetica-Bold,black",
        offset="j0.05i/0.18i",
    )
    fig.text(
        position="TL",
        text="Time: 2009-12-02 21:01:42",
        font="10p,Helvetica-Bold,black",
        offset="j0.05i/0.31i",
    )
    fig.text(
        position="TL",
        text="Depth: 596 km",
        font="10p,Helvetica-Bold,black",
        offset="j0.05i/0.44i",
    )


def plot_waveform(
    fig: pygmt.Figure,
    sta: str,
):
    x, data, data_prediction, poffset, soffset = load_waveform(sta)
    fig.shift_origin(yshift="r0i")
    fig.basemap(
        region=[0, 120, 0, 1],
        projection="X5i/0.8i",
        frame=["WSen", "xaf+lTime(s)", "y1.0f+lPredictions"],
    )
    fig.plot(
        x=x,
        y=data_prediction[0],
        pen="1p,red",
        label="P",
    )
    fig.plot(
        x=x,
        y=data_prediction[1],
        pen="1p,green",
        label="S",
    )
    fig.legend(position="jTR+o0.2c", box="+gwhite+p1p")
    fig.plot(x=[poffset, poffset], y=[0, 1], pen="1p,black,dashed")
    fig.plot(x=[soffset, soffset], y=[0, 1], pen="1p,black,dashed")

    d = np.array(data)
    range_abs = np.max(np.abs(d))
    range = [-range_abs, range_abs]

    fig.shift_origin(yshift="r0.8i")
    fig.basemap(
        region=[0, 120] + range, projection="X5i/0.8i", frame=["wsen", "xaf", "yaf"]
    )
    fig.plot(
        x=x,
        y=data[2],
        pen="1p,black",
    )
    fig.text(
        position="TL",
        text="Z",
        font="15p,Helvetica-Bold,black",
        offset="j0i/0i",
        no_clip=True,
    )
    fig.plot(x=[poffset, poffset], y=range, pen="1p,black,dashed")
    fig.plot(x=[soffset, soffset], y=range, pen="1p,black,dashed")

    fig.shift_origin(yshift="r0.8i")
    fig.basemap(
        region=[0, 120] + range, projection="X5i/0.8i", frame=["wsen", "xaf", "yaf"]
    )
    fig.plot(
        x=x,
        y=data[1],
        pen="1p,black",
    )
    fig.text(
        position="TL",
        text="2",
        font="15p,Helvetica-Bold,black",
        offset="j0i/0i",
        no_clip=True,
    )
    fig.plot(x=[poffset, poffset], y=range, pen="1p,black,dashed")
    fig.plot(x=[soffset, soffset], y=range, pen="1p,black,dashed")

    fig.shift_origin(yshift="r0.8i")
    fig.basemap(
        region=[0, 120] + range,
        projection="X5i/0.8i",
        frame=[f"wsen+t{sta}", "xaf", "yaf"],
    )
    fig.plot(
        x=x,
        y=data[0],
        pen="1p,black",
    )
    fig.text(
        position="TL",
        text="1",
        font="15p,Helvetica-Bold,black",
        offset="j0i/0i",
        no_clip=True,
    )
    fig.plot(x=[poffset, poffset], y=range, pen="1p,black,dashed")
    fig.plot(x=[soffset, soffset], y=range, pen="1p,black,dashed")


def main():
    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL="10p", MAP_LABEL_OFFSET="5p", FONT_ANNOT_PRIMARY="10p")

    fig.shift_origin(yshift="f5.5i")
    plot_map(fig)

    fig.shift_origin(xshift="f1i", yshift="f1i")
    plot_waveform(fig, "A03")

    fig.shift_origin(xshift="f10i", yshift="f7.5i")
    plot_waveform(fig, "C06")

    fig.shift_origin(xshift="f8i", yshift="f1.5i")
    plot_waveform(fig, "TNGA")

    save_path(fig, "map_with_waveform", suffix="jpg")
