""" 
associated_events_vertical.py

Plot the vertical distribution of associated events along cross-sections.
"""
import numpy as np
import pandas as pd
import pygmt
import xarray as xr

from phasenettf import resource, save_path
from phasenettf.utils.slice import extend_line, slab_interp
from obspy.geodetics.base import degrees2kilometers

# configs
LENGTH = 7
LINES = [
    [-179, -16, -173, -19],
    [-179.5, -17, -173.5, -20],
    [-180, -18, -174, -21],
    [-180.5, -19, -174.5, -22],
    [-181, -20, -175, -23],
]


def load_catalog():
    data = pd.read_csv(
        resource(
            ["catalog", "tomoDD.all_months_threshold0.loc"],
            normal_path=True,
            check=True,
        ),
        sep="\s+",
        usecols=[1, 2, 3],
        names=["lat", "lon", "depth"],
    )
    return data


def project_catalog(
    catalog: pd.DataFrame,
    startlon: float,
    startlat: float,
    endlon: float,
    endlat: float,
) -> pd.DataFrame:
    # change column names from lat,lon,dep to y,x,z
    catalog = catalog.rename(columns={"lat": "y", "lon": "x", "depth": "z"})
    catalog = catalog.reindex(columns=["x", "y", "z"])
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
    res.columns = ["dist", "dep"]
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

    with fig.subplot(
        nrows=2,
        ncols=3,
        figsize=("9i", "6i"),
        sharex="b",
        sharey="l",
        margins=["0.1i", "0.06i"],
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

            irow, icol = iline // 3, iline % 3
            fig.basemap(
                projection="X?i/-?i",
                frame=["xaf+lDistance (degree)", "yaf+lDepth (km)"],
                region=[0, LENGTH, 0, 800],
                panel=[irow, icol],
            )
            # plot slab2
            fig.plot(
                x=np.linspace(0, LENGTH, len(slab_deps)),
                y=slab_deps,
                pen="1.5p,magenta",
            )
            # plot catalog
            catalog_line = project_catalog(catalog, startlon, startlat, endlon, endlat)
            fig.plot(
                x=catalog_line["dist"],
                y=catalog_line["dep"],
                style="c0.1c",
                pen="0.01c,black",
            )

    save_path(fig, "associated_all_months_threshold0_vertical", suffix="pdf")


if __name__ == "__main__":
    main()
