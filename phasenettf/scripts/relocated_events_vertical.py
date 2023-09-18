""" 
relocated_events_vertical.py

Plot the vertical distribution of associated events along cross-sections.
"""
import numpy as np
import pandas as pd
import pygmt
import xarray as xr

from phasenettf import resource, save_path
from phasenettf.utils.slice import extend_line, slab_interp, model_interp
from obspy.geodetics.base import degrees2kilometers
from scipy import interpolate

# configs
LENGTH = 7
LINES = [
    [-179, -16, -173, -19],
    [-179.5, -17, -173.5, -20],
    [-180, -18, -174, -21],
    [-180.5, -19, -174.5, -22],
    [-181, -20, -175, -23],
]


# def load_catalog():
#     data = pd.read_csv(
#         resource(
#             ["catalog", "manual_picks", "manual_picks_association_relocated.reloc"],
#             normal_path=True,
#             check=True,
#         ),
#         sep="\s+",
#         usecols=[0, 1, 2, 3, 7, 8, 9],
#         names=["event_id", "lat", "lon", "depth", "errorx", "errory", "errorz"],
#     )
#     data = data[data["depth"] > 0]
#     # filter out abs(errors) > 5000
#     data = data[
#         (data["errorx"] < 5000)
#         & (data["errorx"] > -5000)
#         & (data["errory"] < 5000)
#         & (data["errory"] > -5000)
#         & (data["errorz"] < 5000)
#         & (data["errorz"] > -5000)
#     ]
#     return data


def load_catalog():
    data = pd.read_csv(
        resource(
            ["catalog", "manual_picks", "manual_phase_picks.association_filtered.csv"],
            normal_path=True,
            check=True,
        ),
        sep=",",
        usecols=[0, 1, 2, 3],
        names=["event_id", "lat", "lon", "depth"],
        skiprows=1,
    )
    data = data[data["depth"] > 0]
    return data


def load_tomodd_model():
    with open(
        resource(name=["model", "TX2019_global.mod"], normal_path=True, check=True)
    ) as f:
        lines = f.readlines()
    _, no_lons, no_lats, no_deps = lines[0].split()
    lons = np.array(list(map(float, lines[1].split())))
    lats = np.array(list(map(float, lines[2].split())))
    deps = np.array(list(map(float, lines[3].split())))
    # load vp model only
    model = np.zeros((int(no_lons), int(no_lats), int(no_deps)))
    for idep in range(int(no_deps)):
        for ilat in range(int(no_lats)):
            lon_data = lines[idep * int(no_lats) + ilat + 4].split()
            lon_data = list(map(float, lon_data))
            model[:, ilat, idep] = lon_data
    # construct xarray dataarray
    model = xr.DataArray(
        model,
        [("longitude", lons), ("latitude", lats), ("depth", deps)],
        name="vp",
    )
    return model


def load_ak135(parameter: str, copy_model: xr.DataArray) -> xr.DataArray:
    ak135 = np.loadtxt(
        resource(["model", "AK135F_AVG.csv"], normal_path=True), delimiter=","
    )
    h = ak135[:, 0]
    if parameter == "vp":
        v = ak135[:, 2]
    else:
        v = ak135[:, 3]
    f = interpolate.interp1d(h, v, bounds_error=False, fill_value="extrapolate")
    ak135_depth = f(copy_model.depth.data)
    ak135_abs_data = copy_model.copy()
    for index in range(len(copy_model.depth.data)):
        ak135_abs_data.data[:, :, index] = ak135_depth[index]
    return ak135_abs_data


def project_catalog(
    catalog: pd.DataFrame,
    startlon: float,
    startlat: float,
    endlon: float,
    endlat: float,
) -> pd.DataFrame:
    # change column names from lat,lon,dep to y,x,z
    catalog = catalog.rename(columns={"lat": "y", "lon": "x", "depth": "z"})
    catalog = catalog.reindex(columns=["x", "y", "z", "event_id"])
    # project the catalog to the line
    res = pygmt.project(
        data=catalog,
        center=[startlon, startlat],
        endpoint=[endlon, endlat],
        convention="pz",
        unit=True,
        sort=True,
        length=[0, degrees2kilometers(LENGTH)],
        width=[0, 100],
    )

    # change column names back
    res.columns = ["dist", "dep", "event_id"]
    # convert dist from km to degree
    res["dist"] = res["dist"].apply(lambda x: x / degrees2kilometers(1))
    return res


def main():
    fig = pygmt.Figure()
    pygmt.config(FONT_LABEL="15p", MAP_LABEL_OFFSET="10p", FONT_ANNOT_PRIMARY="13p")
    slab_model = xr.open_dataset(
        resource(["slab2", f"ker_slab2_depth.grd"], normal_path=True)
    )
    catalog = load_catalog()
    model = load_tomodd_model()
    ak135 = load_ak135("vp", model)
    model = (model / ak135 - 1) * 100
    # pygmt.makecpt(cmap="jet", series="6/12/0.3", continuous=True, background="o")
    pygmt.makecpt(
        cmap=resource(["cpt", "dvs_6p_nan.cpt"]),
        series=f"-3/3/1",
        continuous=True,
        background="o",
    )

    with fig.subplot(
        nrows=2,
        ncols=3,
        figsize=("9i", "6i"),
        sharex="b",
        sharey="l",
        margins=["0.2i", "0.12i"],
        frame=["WSen"],
        autolabel="(a)",
    ):
        for iline, line in enumerate(LINES):
            startlon, startlat, endlon, endlat = line
            endlon, endlat = extend_line((startlon, startlat), (endlon, endlat), LENGTH)
            points = pygmt.project(
                center=[startlon, startlat], endpoint=[endlon, endlat], generate=0.02
            )
            lons = points.r
            lons[lons < 0] += 360
            lats = points.s
            slab_deps = slab_interp(slab_model, lons, lats)
            lons -= 360
            cross_section = model_interp(model, lons, lats, np.linspace(0, 800, 801))
            cross_section_xarray = xr.DataArray(
                cross_section,
                dims=("h", "v"),
                coords={
                    "h": np.linspace(0, LENGTH * degrees2kilometers(1), len(lons)),
                    "v": np.linspace(0, 800, 801)[::-1],
                },
            )

            irow, icol = iline // 3, iline % 3
            fig.grdimage(cross_section_xarray.T, panel=[irow, icol])
            fig.basemap(
                projection="X?i/-?i",
                frame=["xaf+lDistance (km)", "yaf+lDepth (km)"],
                region=[0, LENGTH * degrees2kilometers(1), 0, 800],
                # panel=[irow, icol],
            )
            # plot slab2
            fig.plot(
                x=np.linspace(0, LENGTH * degrees2kilometers(1), len(slab_deps)),
                y=slab_deps,
                pen="1.5p,magenta",
                # panel=[irow, icol],
            )
            # plot catalog
            catalog_line = project_catalog(catalog, startlon, startlat, endlon, endlat)
            fig.plot(
                x=catalog_line["dist"] * degrees2kilometers(1),
                y=catalog_line["dep"],
                style="c0.1c",
                pen="0.01c,black",
                # panel=[irow, icol],
            )
            # fig.text(
            #     x=catalog_line["dist"] * degrees2kilometers(1),
            #     y=catalog_line["dep"],
            #     text=catalog_line["event_id"],
            #     font="5p,Helvetica-Bold,black",
            # )
    fig.colorbar(
        # justified inside map frame (j) at Top Center (TC)
        position="JBC+w6i/0.8c+h+o0i/0.3i",
        box=False,
        frame=["a1f"],
        scale=1,
    )

    save_path(
        fig,
        "manual_picks.manual_phase_picks.association_filtered.vertical",
        suffix="pdf",
    )


if __name__ == "__main__":
    main()
