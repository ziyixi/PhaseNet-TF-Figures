from os.path import dirname, exists, join
from typing import List, Union

import pygmt

root_path = dirname(__file__)
resource_path = join(root_path, "data")

__version__ = "0.1.0"


def gmt_path(path: str) -> str:
    """wrap "" outside the normal path
    Args:
        path (str): normal path
    Returns:
        str: wrapped path
    """
    return f'"{path}"'


def resource(name: Union[str, List[str]], normal_path: bool = False,  check: bool = True) -> str:
    """get resource path from its name
    Args:
        name (str): the resource name, might be like [a,b] represent
        normal_path (bool): if the normal path (not gmt path), default is False
        check (bool): if check the file exists or not
    Raises:
        Exception: the input {name} doesn't have resource in {res}!
    Returns:
        str: the resouce path in the package
    """
    if type(name) == str:
        res = join(resource_path, name)
    else:
        res = join(resource_path, *name)

    if (not check) or exists(res):
        return res if normal_path else gmt_path(res)
    else:
        raise Exception(f"the input {name} doesn't have resource in {res}!")


def save_path(fig: pygmt.Figure, name: str, fig_dir_name: str = "fig", suffix: str = "pdf") -> None:
    """get the figure save path as the same level of the package
    Args:
        name (str): figure name
        name (str): saving directory as the same level of the package, default as fig
    Returns:
        str: the saving directory
    """
    base_path = dirname(root_path)
    fig_path = join(base_path, fig_dir_name, name+"."+suffix)
    fig.savefig(fig_path)