import re
import tempfile

import numpy as np
from phasenettf import resource


def generate_tmp_file(content: str = "", suffix: str = "") -> str:
    """write content to a temporary file and return the file path

    Args:
        content (str): the content of the file
        suffix (str): the suffix of the file

    Returns:
        str: the tmp file path
    """
    tmp = tempfile.NamedTemporaryFile(
        delete=False, suffix=(suffix if suffix != "" else None)
    )
    with open(tmp.name, "w") as f:
        f.write(content)
    return tmp.name


def get_vol_list() -> np.ndarray:
    with open(resource(["Volcanoes", "volcanoes.tsv"], normal_path=True), "r") as f:
        data = f.readlines()
    # handle data
    pattern = re.compile(
        r""""([^"]*)"\s+"([^"]*)"\s+"[^"]*"\s+(\d+\.?\d*)\s+(\d+\.?\d*)\s+-?\d+\.?\d*"""
    )
    vol_list = []
    for each_line in data:
        thefind = re.findall(pattern, each_line)
        if len(thefind) != 0:
            result = thefind[0]
            if len(result) != 0:
                if (
                    result[1] not in ["Japan", "Philippines", "Indonesia", "Taiwan"]
                    and float(result[3]) <= 140
                ):
                    vol_list.append(result)
    vol_list = np.array(vol_list)
    return vol_list[:, 2:].astype(float)
