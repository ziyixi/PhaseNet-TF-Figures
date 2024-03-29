""" 
gcmt.py

handle gcmt related problems
"""
from os.path import join
from typing import List

import numpy as np
import obspy
import pandas as pd
from . import generate_tmp_file


def gcmt_to_psmeca(
    gcmt_dir: str,
    has_text: bool = False,
    use_tmp_file=False,
    select_list=[],
    newtext_list=[],
    plot_longitude=[],
    plot_latitude=[],
) -> str:
    """Generate gcmt temp file for psmeca plotting

    Args:
        gcmt_dir (str): gcmt files directory
        has_text (bool): if annotate texts

    Returns:
        str: the temp psmeca plotting path, wrapped as gmt_path
    """

    def get_tensor_dict(files: str) -> dict:
        result = {}
        for item in obspy.read_events(files):
            id = item.origins[0].resource_id.id.split("/")[2]
            tensor = item.focal_mechanisms[0].moment_tensor.tensor
            latitude = item.origins[0].latitude
            longitude = item.origins[0].longitude
            depth = item.origins[0].depth
            result[id] = (tensor, longitude, latitude, depth)
        return result

    def split_tensor_exponent(tensor):
        result = {}
        search_list = np.array(
            [
                tensor.m_rr,
                tensor.m_tt,
                tensor.m_pp,
                tensor.m_rt,
                tensor.m_rp,
                tensor.m_tp,
            ]
        )
        search_list = np.abs(search_list)
        ref = np.min(search_list)

        exp = len(str(int(ref))) - 1
        result = {
            "m_rr": tensor.m_rr / (10**exp),
            "m_tt": tensor.m_tt / (10**exp),
            "m_pp": tensor.m_pp / (10**exp),
            "m_rt": tensor.m_rt / (10**exp),
            "m_rp": tensor.m_rp / (10**exp),
            "m_tp": tensor.m_tp / (10**exp),
            "exp": exp,
        }
        return result

    tensor_dict = get_tensor_dict(join(gcmt_dir, "*"))
    if not use_tmp_file:
        res = []
        for key in tensor_dict:
            item, longitude, latitude, depth = tensor_dict[key]
            tensor = split_tensor_exponent(item)
            cur = {
                "longitude": longitude,
                "latitude": latitude,
                "depth": depth,
                "event_name": key if has_text else "",
                "mrr": tensor["m_rr"],
                "mtt": tensor["m_tt"],
                "mff": tensor["m_pp"],
                "mrt": tensor["m_rt"],
                "mrf": tensor["m_rp"],
                "mtf": tensor["m_tp"],
                "exponent": tensor["exp"],
            }
            res.append(cur)
        # convert to dataframe, with columns as the key in cur dict
        res = pd.DataFrame(res)

        # the path is always used in gmt script
        return res
    else:
        tmp_file = generate_tmp_file()
        with open(tmp_file, "w") as f:
            for key in tensor_dict:
                if key not in select_list:
                    continue
                newtext = newtext_list[select_list.index(key)]
                newlongitude = (
                    plot_longitude[select_list.index(key)] if plot_longitude else 0
                )
                newlatitude = (
                    plot_latitude[select_list.index(key)] if plot_latitude else 0
                )
                item, longitude, latitude, depth = tensor_dict[key]
                tensor = split_tensor_exponent(item)
                if has_text:
                    f.write(
                        f'{longitude} {latitude} {depth/1000:.2f} {tensor["m_rr"]:.3f} {tensor["m_tt"]:.3f} {tensor["m_pp"]:.3f} {tensor["m_rt"]:.3f} {tensor["m_rp"]:.3f} {tensor["m_tp"]:.3f} {tensor["exp"]} {newlongitude} {newlatitude} {newtext}\n'
                    )
                else:
                    f.write(
                        f'{longitude} {latitude} {depth/1000:.2f} {tensor["m_rr"]:.3f} {tensor["m_tt"]:.3f} {tensor["m_pp"]:.3f} {tensor["m_rt"]:.3f} {tensor["m_rp"]:.3f} {tensor["m_tp"]:.3f} {tensor["exp"]} {newlongitude} {newlatitude} \n'
                    )

        # the path is always used in gmt script
        return tmp_file
