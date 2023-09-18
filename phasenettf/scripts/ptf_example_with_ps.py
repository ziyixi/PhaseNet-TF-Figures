"""
manuscript_phasenet_tf_example.py:

This script is used to generate PhaseNet-TF test examples to show the performance.
"""
from pathlib import Path

import h5py
import httpx
import numpy as np
import obspy
import pandas as pd
import pygmt
import torch
import xarray as xr

from phasenettf import resource, save_path
from phasenettf.utils.spectrogram import GenSgram
from string import ascii_lowercase

XOFFSET = 8.6
YOFFSET = 9.6
X_SHIFTS = [1.0, XOFFSET, -XOFFSET, XOFFSET]
Y_SHIFTS = [YOFFSET + 1.8, 0, -YOFFSET, 0]
FRAMES = {
    0: {
        "waveform": ["Wsen", "xaf", "ya0.8f0.4+lAmplitude"],
        "spectrogram": ["Wsen", "xaf", "ya4f2+lFrequency"],
        "prediction": ["Wsen", "xaf", "ya0.4f0.2+lProbability"],
    },
    1: {
        "waveform": ["wsen", "xaf", "ya0.8f0.4"],
        "spectrogram": ["wsen", "xaf", "ya4f2"],
        "prediction": ["wsen", "xaf", "ya0.4f0.2"],
    },
    2: {
        "waveform": ["Wsen", "xaf", "ya0.8f0.4+lAmplitude"],
        "spectrogram": ["Wsen", "xaf", "ya4f2+lFrequency"],
        "prediction": ["WSen", "xaf+lTime (s)", "ya0.4f0.2+lProbability"],
    },
    3: {
        "waveform": ["wsen", "xaf", "ya0.8f0.4"],
        "spectrogram": ["wsen", "xaf", "ya4f2"],
        "prediction": ["wSen", "xaf+lTime (s)", "ya0.4f0.2"],
    },
}
WAVEFORMS = [
    ["11_52111", "B01W"],
    ["11_52113", "A12W"],
    ["11_52113", "EUAP"],
    ["11_52113", "FONI"],
]
MAX_CLAMP = 3


def main():
    fig = pygmt.Figure()
    pygmt.config(
        FONT_LABEL="18p",
        MAP_LABEL_OFFSET="12p",
        FONT_ANNOT_PRIMARY="16p",
        MAP_FRAME_TYPE="plain",
        MAP_TITLE_OFFSET="12p",
        FONT_TITLE="18p,black",
        MAP_FRAME_PEN="1p,black",
    )

    for iplot in range(len(X_SHIFTS)):
        fig.shift_origin(xshift=f"{X_SHIFTS[iplot]}i", yshift=f"{Y_SHIFTS[iplot]}i")
        with fig.subplot(
            nrows=7,
            ncols=1,
            subsize=(f"{5/0.618:.2f}i", "1.3i"),
            margins=["0.2i", "-0.141i"],
        ):
            (
                st,
                sgram,
                p_prediction,
                s_prediction,
                arrivals,
                arrival_types,
                component,
            ) = load_waveform_spectrogram_prediction(
                WAVEFORMS[iplot][0], WAVEFORMS[iplot][1]
            )
            for ipanel in range(3):
                with fig.set_panel(panel=ipanel):
                    fig.basemap(
                        region=[0, 120, -1, 1],
                        projection="X?i/?i",
                        frame=FRAMES[iplot]["waveform"],
                    )
                    fig.plot(
                        x=np.linspace(0, 120, 4800),
                        y=st[ipanel].data,
                        pen="1p,black",
                    )
                    plot_manual_labels(fig, arrival_types, arrivals)
                    fig.text(
                        x=3,
                        y=0.8,
                        text=f"{component[ipanel]}",
                        font="18p,Helvetica-Bold,black",
                    )
                    if ipanel == 0:
                        fig.text(
                            x=115,
                            y=0.75,
                            text=f"({ascii_lowercase[iplot]})",
                            font="30p,Helvetica-Bold,black",
                        )
                        fig.text(
                            x=90,
                            y=0.75,
                            text=f"{WAVEFORMS[iplot][0]}.{WAVEFORMS[iplot][1]}",
                            font="18p,Helvetica-Bold,black",
                        )
            pygmt.makecpt(
                cmap="jet", series=[0, MAX_CLAMP, 0.1], continuous=True, reverse=False
            )
            for ipanel in range(3, 6):
                with fig.set_panel(panel=ipanel):
                    fig.basemap(
                        region=[0, 120, 0, 10],
                        projection="X?i/?i",
                        frame=FRAMES[iplot]["spectrogram"],
                    )
                    fig.grdimage(
                        sgram[ipanel - 3],
                        cmap=True,
                    )
                    plot_manual_labels(fig, arrival_types, arrivals)
                    fig.text(
                        x=3,
                        y=9,
                        text=f"{component[ipanel-3]}",
                        font="18p,Helvetica-Bold,white",
                    )
            with fig.set_panel(panel=6):
                fig.basemap(
                    region=[0, 120, 0, 1],
                    projection="X?i/?i",
                    frame=FRAMES[iplot]["prediction"],
                )
                fig.plot(
                    x=np.linspace(0, 120, len(p_prediction)),
                    y=p_prediction,
                    pen="1p,red",
                    label="Predicted P",
                )
                fig.plot(
                    x=np.linspace(0, 120, len(s_prediction)),
                    y=s_prediction,
                    pen="1p,magenta",
                    label="Predicted S",
                )
                plot_manual_labels(fig, arrival_types, arrivals, with_label=True)
                fig.legend(
                    position="JTR+jTR+o0.2c/0.2c",
                    box="+gwhite+p1p,black",
                )

    fig.shift_origin(xshift=f"{-XOFFSET/2+1}i")
    pygmt.makecpt(
        cmap="jet", series=[0, MAX_CLAMP, 0.1], continuous=True, reverse=False
    )
    fig.colorbar(
        position="JBC+w5i/0.8c+h+o0i/1.8c",
        box=False,
        frame=["a0.5f", f'"+LSpecrogram Amplitude"'],
        scale=1,
    )

    save_path(fig, Path(__file__).resolve().stem)


def load_waveform_spectrogram_prediction(event_id, station_id):
    waveform_h5 = resource(["waveforms", "waveform.h5"], normal_path=True)
    raw = h5py.File(waveform_h5, "r")[event_id][station_id]
    raw_wave = raw[:, 9200:14000]

    # * waveform
    st = obspy.Stream()
    for i in range(len(raw_wave)):
        st.append(obspy.Trace(data=raw_wave[i], header={"delta": 0.025}))
    st.filter("bandpass", freqmin=1, freqmax=10, corners=4, zerophase=True)
    st.normalize(global_max=False)

    arrivals = raw.attrs["phase_index"]
    arrivals = [item - 9200 for item in arrivals]
    arrival_types = raw.attrs["phase_type"]
    component = raw.attrs["component"]

    # * spectrogram
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
            "time": np.linspace(0, 120, 4800),
        },
    )

    # * prediction
    id = f"{event_id}.{station_id}"
    timestamp = "1970-01-01T00:00:00.000000"
    w = wave.squeeze(0).numpy().tolist()
    sensitivity = 0.5
    return_prediction = True
    request_body = {
        "id": id,
        "timestamp": timestamp,
        "waveform": w,
        "sensitivity": sensitivity,
        "return_prediction": return_prediction,
    }
    url = "http://0.0.0.0:8080/api/predict"
    headers = {"Content-Type": "application/json"}
    response = httpx.post(url, headers=headers, json=request_body, timeout=6000)
    if response.status_code == 200:
        print("success", id)
    else:
        print("fail", id)
        raise Exception(response.text)
    p_prediction = np.array(response.json()["prediction"]["P"][:4800])
    s_prediction = np.array(response.json()["prediction"]["S"][:4800])

    return st, sgram, p_prediction, s_prediction, arrivals, arrival_types, component


def plot_manual_labels(fig, arrival_types, arrivals, with_label=False):
    colors = ["red", "magenta", "black"]
    for iphase, phase in enumerate(["P", "S", "PS"]):
        for itype, type in enumerate(arrival_types):
            if type == phase:
                fig.plot(
                    x=[
                        arrivals[itype] * 0.025,
                        arrivals[itype] * 0.025,
                    ],
                    y=[-20, 20],
                    pen=f"2p,{colors[iphase]},.-.",
                    label=f"Manually picked {phase}" if with_label else None,
                )
