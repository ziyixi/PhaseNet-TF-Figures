from pathlib import Path

import numpy as np
import obspy
import pygmt
import torch
import xarray as xr
from obspy import Stream

from phasenettf import resource, save_path
from phasenettf.utils.spectrogram import GenSgram

WAVEFORM_LENGTH = 3600
EXAMPLE_WAVEFORM_1 = resource(
    [
        "waveforms",
        "continuous",
        "YL.B12.2010-04-07.waveform.mseed",
    ],
    normal_path=True,
)
EXAMPLE_PREDICTION_1 = resource(
    [
        "waveforms",
        "continuous",
        "YL.B12.2010-04-07.prediction.mseed",
    ],
    normal_path=True,
)
EXAMPLE_WAVEFORM_2 = resource(
    [
        "waveforms",
        "continuous",
        "Z1.VAVP.2010-04-07.waveform.mseed",
    ],
    normal_path=True,
)
EXAMPLE_PREDICTION_2 = resource(
    [
        "waveforms",
        "continuous",
        "Z1.VAVP.2010-04-07.prediction.mseed",
    ],
    normal_path=True,
)
TITLE_1 = "YL.B12 2010-04-07 0h to 1h"
TITLE_2 = "Z1.VAVP 2010-04-07 0h to 1h"
MAX_CLAMP = 0.5


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

    # * first event
    wave1 = obspy.read(EXAMPLE_WAVEFORM_1)
    wave1.normalize(global_max=False)
    pred1 = obspy.read(EXAMPLE_PREDICTION_1)
    sgram1 = gen_sgram(wave1)
    plot_one_panel(fig, wave1, pred1, sgram1, TITLE_1)

    # * second event
    fig.shift_origin(xshift="w+0.9i", yshift="6.4i")
    wave2 = obspy.read(EXAMPLE_WAVEFORM_2)
    wave2.normalize(global_max=False)
    pred2 = obspy.read(EXAMPLE_PREDICTION_2)
    sgram2 = gen_sgram(wave2)
    plot_one_panel(fig, wave2, pred2, sgram2, TITLE_2)

    # * colorbar
    fig.shift_origin(xshift="-3.45i")
    pygmt.makecpt(
        cmap="jet", series=[0, MAX_CLAMP, 0.02], continuous=True, reverse=False
    )
    fig.colorbar(
        position="JBC+w5i/0.8c+h+o0i/1.8c",
        box=False,
        frame=["a0.5f", f'"+LSpecrogram Amplitude"'],
        scale=1,
    )

    # * save
    save_path(fig, Path(__file__).resolve().stem)


def plot_one_panel(
    fig: pygmt.Figure,
    waveform: Stream,
    possibility: Stream,
    sgram: xr.DataArray,
    title: str,
):
    x = np.linspace(0, WAVEFORM_LENGTH, len(waveform[0].data))

    # * waveform
    for i in range(3):
        fig.basemap(
            region=[0, WAVEFORM_LENGTH, -1, 1],
            projection="X6i/1i",
            frame=["Wsen", "xaf", "y0.8f+lAmplitude"]
            if i != 0
            else [f"Wsen+t{title}", "xaf", "y0.8f+lAmplitude"],
        )
        fig.plot(
            x=x,
            y=waveform[i].data,
            pen="0.5p,black",
        )
        fig.text(
            x=90,
            y=0.8,
            text=f"{waveform[i].id[-1]}",
            font="18p,Helvetica-Bold,black",
        )
        fig.shift_origin(yshift="-h-0i")
    fig.shift_origin(yshift="+h0i")
    fig.shift_origin(yshift="-h-0.2i")

    # * spectrogram
    pygmt.makecpt(
        cmap="jet", series=[0, MAX_CLAMP, 0.02], continuous=True, reverse=False
    )
    for i in range(3):
        fig.basemap(
            region=[0, WAVEFORM_LENGTH, 0, 10],
            projection="X6i/1i",
            frame=["Wsen", "xaf", "y4f+lFrequency (Hz)"],
        )
        fig.grdimage(
            sgram[i],
            cmap=True,
        )
        fig.text(
            x=90,
            y=9,
            text=f"{waveform[i].id[-1]}",
            font="18p,Helvetica-Bold,white",
        )
        fig.shift_origin(yshift="-h-0i")
    fig.shift_origin(yshift="+h0i")
    fig.shift_origin(yshift="-h-0.2i")

    # * possibility
    fig.basemap(
        region=[0, WAVEFORM_LENGTH, 0, 1],
        projection="X6i/1i",
        frame=["WSen", "xaf+lTime (s)", "y0.5f+lPossibility"],
    )
    fig.plot(x=x, y=possibility[0].data, pen="1p,red", label="P Possibility")
    fig.plot(x=x, y=possibility[1].data, pen="1p,blue", label="S Possibility")
    fig.legend(
        position="JTR+jTR+o1.5i/0.2c",
        box="+gwhite+p1p,black",
    )


def gen_sgram(st: Stream):
    length = len(st[0].data)
    sgram_gen = GenSgram(max_clamp=MAX_CLAMP, width=length)
    wave = np.zeros((3, length), dtype=np.float32)
    for i in range(3):
        wave[i, :] = st[i].data
    wave = torch.from_numpy(wave).unsqueeze(0)
    sgram = sgram_gen(wave)
    sgram = sgram.squeeze(0).numpy()
    # sgram is with shape 3, 64, length, i.e., channel, freq, time, wrap it to xarray
    sgram = xr.DataArray(
        sgram,
        dims=["channel", "freq", "time"],
        coords={
            "channel": ["1", "2", "Z"],
            "freq": np.linspace(0, 10, 64),
            "time": np.linspace(0, WAVEFORM_LENGTH, length),
        },
    )
    return sgram
