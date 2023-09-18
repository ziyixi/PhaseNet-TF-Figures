"""
manuscript_workflow_component:

This script is used to generate the workflow component figure for the manuscript, including waveform and picks, association examples.
"""
import pygmt
from phasenettf import save_path, resource
from pathlib import Path
import h5py
import numpy as np
import obspy
from phasenettf.utils.spectrogram import GenSgram
import torch
import xarray as xr
import httpx
import pandas as pd
import json

# * =========================Global variables=========================
EVENT_ID = "4_54379"
EVENT_ID_PAIRS = ["4_54513", "4_54379"]
STATION_ID = "B03"
MAX_CLAMP = 10
STATION_LIST = ["EUAP", "A06W", "TNGA", "C15", "A09"]
OFFSETS = [-1400, 1200]
REGION = [-184, -172, -26, -14]


# * =========================Main function=========================


def plot_waveform():
    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="14p",
        MAP_LABEL_OFFSET="8p",
        FONT_ANNOT_PRIMARY="12p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="8p",
        FONT_TITLE="14p,black",
        MAP_FRAME_PEN="0.5p,white",
    )
    st, _, arrivals, arrival_types = load_waveform()

    with fig.subplot(
        nrows=3,
        ncols=1,
        subsize=(f"{3/0.618:.2f}i", "1.5i"),
        margins=["0.2i", "-0.11i"],
    ):
        colors = ["red", "magenta"]
        phases = ["P", "S"]
        panel_labels = ["1", "2", "Z"]
        for panel in range(3):
            with fig.set_panel(panel=panel):
                frame_label = ["lbrt"]
                fig.basemap(
                    region=[-10, 110, -1, 1],
                    projection="X?i/?i",
                    frame=frame_label,
                )
                fig.plot(
                    x=np.linspace(-10, 110, 4800),
                    y=st[panel].data,
                    pen="1p,black",
                    label=panel_labels[panel],
                )
                for iphase, phase in enumerate(phases):
                    for itype, type in enumerate(arrival_types):
                        if type == phase:
                            fig.plot(
                                x=[arrivals[itype] * 0.025, arrivals[itype] * 0.025],
                                y=[-1, 1],
                                pen=f"1p,{colors[iphase]},dashed",
                                label=f"{phase}" if panel == 2 else None,
                            )
                fig.legend(
                    position="JTL+jTL+o1.5c/0c",
                )

    save_path(fig, Path(__file__).resolve().stem + ".waveform", suffix=".png")


def plot_waveform_longer():
    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="14p",
        MAP_LABEL_OFFSET="8p",
        FONT_ANNOT_PRIMARY="12p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="8p",
        FONT_TITLE="14p,black",
        MAP_FRAME_PEN="0.5p,white",
    )
    st = load_waveform_longer()

    with fig.subplot(
        nrows=3,
        ncols=1,
        subsize=(f"{6/0.618:.2f}i", "1.5i"),
        margins=["0.2i", "-0.11i"],
    ):
        panel_labels = ["1", "2", "Z"]
        for panel in range(3):
            with fig.set_panel(panel=panel):
                frame_label = ["lbrt"]
                fig.basemap(
                    region=[0, 600, -1, 1],
                    projection="X?i/?i",
                    frame=frame_label,
                )
                fig.plot(
                    x=np.linspace(0, 600, len(st[panel].data)),
                    y=st[panel].data,
                    pen="1p,black",
                    label=panel_labels[panel],
                )
                fig.legend(
                    position="JTL+jTL+o1.5c/0c",
                )

    save_path(fig, Path(__file__).resolve().stem + ".waveform_longer", suffix=".png")


def plot_spectrogram():
    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="14p",
        MAP_LABEL_OFFSET="8p",
        FONT_ANNOT_PRIMARY="12p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="8p",
        FONT_TITLE="14p,black",
        MAP_FRAME_PEN="0.5p,white",
    )
    _, spec, arrivals, arrival_types = load_waveform()
    pygmt.makecpt(
        cmap="jet", series=[0, MAX_CLAMP, 0.1], continuous=True, reverse=False
    )

    with fig.subplot(
        nrows=3,
        ncols=1,
        subsize=(f"{3/0.618:.2f}i", "1.5i"),
        margins=["0.2i", "-0.11i"],
    ):
        colors = ["red", "magenta"]
        phases = ["P", "S"]
        panel_labels = ["1", "2", "Z"]
        for panel in range(3):
            with fig.set_panel(panel=panel):
                frame_label = ["lbrt"]
                fig.basemap(
                    region=[-10, 110, 0, 10],
                    projection="X?i/?i",
                    frame=frame_label,
                )
                fig.grdimage(
                    spec[panel],
                    cmap=True,
                    shading=True,
                    frame=frame_label,
                )
                for iphase, phase in enumerate(phases):
                    for itype, type in enumerate(arrival_types):
                        if type == phase:
                            fig.plot(
                                x=[arrivals[itype] * 0.025, arrivals[itype] * 0.025],
                                y=[0, 64],
                                pen=f"1p,{colors[iphase]},dashed",
                                label=f"{phase}" if panel == 2 else None,
                            )
                # text the pabel labels in the top left corner
                fig.text(
                    x=-5,
                    y=9.5,
                    text=panel_labels[panel],
                    font="14p,Helvetica-Bold,white",
                )

    save_path(fig, Path(__file__).resolve().stem + ".spectrogram", suffix=".png")


def plot_prediction():
    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="14p",
        MAP_LABEL_OFFSET="8p",
        FONT_ANNOT_PRIMARY="12p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="8p",
        FONT_TITLE="14p,black",
        MAP_FRAME_PEN="0.5p,white",
    )
    st = load_waveform_longer()
    p_prediction, s_prediction = get_phasenettf_predictions()

    fig.basemap(
        region=[0, 1600, 0, 1],
        projection=f"X{3/0.618:.2f}i/1.5i",
        frame=["wsrt", "xaf+lTime (s)"],
    )
    fig.plot(
        x=np.linspace(0, 1600, len(p_prediction)),
        y=p_prediction,
        pen="1p,red",
        label="Predicted P",
    )
    fig.plot(
        x=np.linspace(0, 1600, len(s_prediction)),
        y=s_prediction,
        pen="1p,magenta",
        label="Predicted S",
    )
    fig.legend(
        position="JTL+jTL+o1.5c/0c",
    )

    save_path(fig, Path(__file__).resolve().stem + ".prediction")


def plot_association():
    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="14p",
        MAP_LABEL_OFFSET="8p",
        FONT_ANNOT_PRIMARY="12p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="8p",
        FONT_TITLE="14p,black",
        MAP_FRAME_PEN="0.5p,white",
    )
    predictions, waves = load_association_waveforms()

    # waveform 1
    fig.shift_origin(yshift="4.5i")
    with fig.subplot(
        nrows=5,
        ncols=1,
        figsize=("4i", "4i"),
        margins=["0.2i", "-0.11i"],
    ):
        for panel in range(5):
            with fig.set_panel(panel=panel):
                fig.basemap(
                    region=[0, 150, -1, 1],
                    projection="X?i/?i",
                    frame=["lbrt"],
                )
                # only plot z component
                fig.plot(
                    x=np.linspace(0, 150, 6000),
                    y=waves[EVENT_ID_PAIRS[0]][STATION_LIST[panel]][2].data,
                    pen="1p,blue",
                )
                # plot predicted P and S in the same panel
                p = np.array(predictions[EVENT_ID_PAIRS[0]][STATION_LIST[panel]]["P"])
                s = np.array(predictions[EVENT_ID_PAIRS[0]][STATION_LIST[panel]]["S"])
                fig.plot(
                    x=np.linspace(0, 150, 6000),
                    y=p * 2 - 1,
                    pen="1p,red",
                )
                fig.plot(
                    x=np.linspace(0, 150, 6000),
                    y=s * 2 - 1,
                    pen="1p,magenta",
                )
                # text at the top left corner for station name
                fig.text(
                    x=10,
                    y=0.8,
                    text=STATION_LIST[panel],
                    font="12p,Helvetica-Bold,black",
                )

    # waveform 2
    fig.shift_origin(yshift="-4.5i")
    with fig.subplot(
        nrows=5,
        ncols=1,
        figsize=("4i", "4i"),
        margins=["0.2i", "-0.11i"],
    ):
        for panel in range(5):
            with fig.set_panel(panel=panel):
                fig.basemap(
                    region=[0, 150, -1, 1],
                    projection="X?i/?i",
                    frame=["lbrt"],
                )
                # only plot z component
                fig.plot(
                    x=np.linspace(0, 150, 6000),
                    y=waves[EVENT_ID_PAIRS[1]][STATION_LIST[panel]][2].data,
                    pen="1p,orange",
                )
                # plot predicted P and S in the same panel
                p = np.array(predictions[EVENT_ID_PAIRS[1]][STATION_LIST[panel]]["P"])
                s = np.array(predictions[EVENT_ID_PAIRS[1]][STATION_LIST[panel]]["S"])
                fig.plot(
                    x=np.linspace(0, 150, 6000),
                    y=p * 2 - 1,
                    pen="1p,red",
                )
                fig.plot(
                    x=np.linspace(0, 150, 6000),
                    y=s * 2 - 1,
                    pen="1p,magenta",
                )
                # text at the top left corner for station name
                fig.text(
                    x=10,
                    y=0.8,
                    text=STATION_LIST[panel],
                    font="12p,Helvetica-Bold,black",
                )

    # map
    with pygmt.config(MAP_FRAME_PEN="1p,black"):
        fig.shift_origin(yshift="1.5i", xshift="4.6i")
        fig.basemap(
            region=[-184, -172, -26, -14],
            projection="M6i",
            frame=["WSen", "xaf", "yaf"],
        )
        plot_earth_relief(fig)
        plot_text(fig)
        plot_stations(fig)
        plot_events(fig)
        plot_inset(fig)

    save_path(fig, Path(__file__).resolve().stem + ".association", suffix=".png")


# * =========================Helper functions=========================
def load_waveform():
    waveform_h5 = resource(["waveforms", "waveform.h5"], normal_path=True)
    waveform = h5py.File(waveform_h5, "r")
    wave = waveform[EVENT_ID][STATION_ID][:, 9200:14000]

    st = obspy.Stream()
    for i in range(len(wave)):
        st.append(obspy.Trace(data=wave[i], header={"delta": 0.025}))
    st.filter("bandpass", freqmin=1, freqmax=10, corners=4, zerophase=True)
    st.normalize(global_max=False)

    arrivals = waveform[EVENT_ID][STATION_ID].attrs["phase_index"]
    arrivals = [item - 9600 for item in arrivals]
    arrival_types = waveform[EVENT_ID][STATION_ID].attrs["phase_type"]
    sgram = wave_to_spectrogram(st)
    return st, sgram, arrivals, arrival_types


def load_waveform_longer():
    waveform_h5 = resource(["waveforms", "waveform.h5"], normal_path=True)
    waveform = h5py.File(waveform_h5, "r")
    wave = waveform[EVENT_ID_PAIRS[1]][STATION_LIST[0]][:, :]

    st = obspy.Stream()
    for i in range(len(wave)):
        st.append(obspy.Trace(data=wave[i], header={"delta": 0.025}))
    # st.filter("bandpass", freqmin=1, freqmax=10, corners=4, zerophase=True)
    st.normalize(global_max=False)

    return st


def load_association_waveforms():
    waveform_h5 = resource(["waveforms", "waveform.h5"], normal_path=True)
    waveform = h5py.File(waveform_h5, "r")
    predictions = {}
    waves = {}
    for iev, event_id in enumerate(EVENT_ID_PAIRS):
        predictions[event_id] = {}
        waves[event_id] = {}
        for station_id in STATION_LIST:
            id = f"{event_id}.{station_id}"
            offset = int(
                (
                    pd.Timestamp(waveform[event_id].attrs["event_time"])
                    - pd.Timestamp(waveform[event_id][station_id].attrs["begin_time"])
                ).total_seconds()
                * 40
            )
            wave = waveform[event_id][station_id][
                :, offset + OFFSETS[iev] : offset + OFFSETS[iev] + 6000
            ]
            timestamp = waveform[event_id][station_id].attrs["begin_time"]
            w = wave.tolist()
            sensitivity = 0.5
            return_prediction = True

            headers = {"Content-Type": "application/json"}

            request_body = {
                "id": id,
                "timestamp": timestamp,
                "waveform": w,
                "sensitivity": sensitivity,
                "return_prediction": return_prediction,
            }
            url = "http://0.0.0.0:8080/api/predict"
            response = httpx.post(url, headers=headers, json=request_body, timeout=6000)
            if response.status_code == 200:
                print("success", id)
            else:
                print("fail", id)
                raise Exception(response.text)
            p_prediction = response.json()["prediction"]["P"]
            s_prediction = response.json()["prediction"]["S"]
            predictions[event_id][station_id] = {
                "P": p_prediction[:6000],
                "S": s_prediction[:6000],
            }
            # waves
            st = obspy.Stream()
            for i in range(len(wave)):
                st.append(obspy.Trace(data=wave[i], header={"delta": 0.025}))
            st.filter("bandpass", freqmin=1, freqmax=10, corners=4, zerophase=True)
            st.normalize(global_max=False)
            waves[event_id][station_id] = st
    return predictions, waves


def wave_to_spectrogram(st):
    sgram_gen = GenSgram(max_clamp=MAX_CLAMP)
    # wrap obspy wave into torch tensor with batch==1
    wave = np.zeros((3, 4800), dtype=np.float32)
    for i in range(3):
        wave[i, :] = st[i].data
    wave = torch.from_numpy(wave).unsqueeze(0)
    sgram = sgram_gen(wave)
    sgram = sgram.squeeze(0).numpy()
    # sgram is with shape 3, 64, 4800, i.e., channel, freq, time, wrap it to xarray
    sgram = xr.DataArray(
        sgram,
        dims=["channel", "freq", "time"],
        coords={
            "channel": ["1", "2", "Z"],
            "freq": np.linspace(0, 10, 64),
            "time": np.linspace(-10, 110, 4800),
        },
    )

    return sgram


def get_phasenettf_predictions():
    # docker run --rm -p 8080:8080 --name pt ghcr.io/ziyixi/phasenet-tf:latest
    waveform_h5 = resource(["waveforms", "waveform.h5"], normal_path=True)
    waveform = h5py.File(waveform_h5, "r")
    wave = waveform[EVENT_ID_PAIRS[1]][STATION_LIST[0]][:, :]

    id = f"{EVENT_ID_PAIRS[0]}.{STATION_LIST[0]}"
    timestamp = waveform[EVENT_ID_PAIRS[0]][STATION_LIST[0]].attrs["begin_time"]
    waveform = wave.tolist()
    sensitivity = 0.5
    return_prediction = True

    headers = {"Content-Type": "application/json"}

    request_body = {
        "id": id,
        "timestamp": timestamp,
        "waveform": waveform,
        "sensitivity": sensitivity,
        "return_prediction": return_prediction,
    }

    url = "http://0.0.0.0:8080/api/predict"
    response = httpx.post(url, headers=headers, json=request_body, timeout=6000)

    # Check if the request was successful (HTTP status code 200)
    if response.status_code == 200:
        print("Request was successful!")
        # print("Response:", response.json())
    else:
        print("Request failed with status code:", response.status_code)
        raise Exception(response.text)
        # print("Response:", response.json())
    p_prediction = response.json()["prediction"]["P"]
    s_prediction = response.json()["prediction"]["S"]
    return p_prediction, s_prediction


# * =========================Map=========================
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
    selected_station_df = []
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
        if sta in STATION_LIST:
            selected_station_df.append(
                [
                    net,
                    sta,
                    stations[station]["longitude"],
                    stations[station]["latitude"],
                ]
            )
    station_df = pd.DataFrame(station_df, columns=["net", "sta", "lon", "lat"])
    selected_station_df = pd.DataFrame(
        selected_station_df, columns=["net", "sta", "lon", "lat"]
    )
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
    if len(selected_station_df[selected_station_df["net"] == "YL"]) > 0:
        fig.plot(
            x=selected_station_df[selected_station_df["net"] == "YL"]["lon"],
            y=selected_station_df[selected_station_df["net"] == "YL"]["lat"],
            style="i0.3c",
            pen="2p,red",
            label="Selected stations",
        )
    if len(selected_station_df[selected_station_df["net"] == "Z1"]) > 0:
        fig.plot(
            x=selected_station_df[selected_station_df["net"] == "Z1"]["lon"],
            y=selected_station_df[selected_station_df["net"] == "Z1"]["lat"],
            style="t0.4c",
            pen="2p,red",
        )
    if len(selected_station_df[selected_station_df["net"] == "II"]) > 0:
        fig.plot(
            x=selected_station_df[selected_station_df["net"] == "II"]["lon"],
            y=selected_station_df[selected_station_df["net"] == "II"]["lat"],
            style="s0.4c",
            pen="2p,red",
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
    # pygmt.makecpt(cmap="oslo", series=[0, 700, 1], continuous=True, reverse=True)
    # fig.plot(
    #     x=df_events["longitude"],
    #     y=df_events["latitude"],
    #     style="c0.15c",
    #     fill=df_events["depth"],
    #     cmap=True,
    # )
    # # plot the colorbar to the right of the map
    # fig.colorbar(position="JMR+o0.2c/0c+w5i/0.3i", frame=["x+lDepth (km)"])

    # plot first event
    fig.plot(
        x=df_events[df_events["originid"] == EVENT_ID_PAIRS[0]]["longitude"],
        y=df_events[df_events["originid"] == EVENT_ID_PAIRS[0]]["latitude"],
        style="a0.9c",
        fill="orange",
    )
    # plot second event
    fig.plot(
        x=df_events[df_events["originid"] == EVENT_ID_PAIRS[1]]["longitude"],
        y=df_events[df_events["originid"] == EVENT_ID_PAIRS[1]]["latitude"],
        style="a0.9c",
        fill="blue",
    )


def main():
    # plot_waveform()
    # plot_waveform_longer()
    # plot_spectrogram()
    plot_prediction()
    # plot_association()
