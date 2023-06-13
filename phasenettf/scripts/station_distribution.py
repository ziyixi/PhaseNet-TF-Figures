import obspy
import pygmt

from phasenettf import resource, save_path

STATION_XMLS = resource(["tongaml_continious_xml", "*xml"],
                        normal_path=True, check=False)
REGION = [-184, -172, -24, -14]


def prepare_data():
    invs = obspy.read_inventory(STATION_XMLS)
    res = {}
    for inv in invs:
        res[inv[0].code] = (inv[0].longitude, inv[0].latitude)
    return res


def main():
    stations = prepare_data()

    fig = pygmt.Figure()
    # basemap
    with pygmt.config(MAP_FRAME_TYPE="plain", MAP_TICK_LENGTH="0p"):
        fig.basemap(region=REGION, projection="M8i",
                    frame=["WSen", "xafg", "yafg"])

    # load topo
    grd_topo = pygmt.datasets.load_earth_relief(
        resolution="01m", region=REGION)

    # * base lines
    fig.coast(water="167/194/223")
    fig.grdimage(grd_topo, cmap=resource(
        ["cpt", "land_sea.cpt"]), shading="+d")

    # fig.plot(data=resource(
    #     ["Plate_Boundaries", "nuvel1_boundaries"]), pen="2p,red")

    for sta in stations:
        fig.plot(x=stations[sta][0], y=stations[sta][1],
                 style="t0.3c", pen="1p,black", fill="black")
        fig.text(x=stations[sta][0], y=stations[sta][1] +
                 0.2, text=sta, font="5p,Helvetica-Bold,red")

    save_path(fig, "station_distribution")
