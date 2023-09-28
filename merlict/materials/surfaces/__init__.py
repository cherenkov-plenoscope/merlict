import os
import json_numpy
from ... import utils


def get_resources_dir():
    return utils.resources.path("materials", "surfaces", "resources")


def list_resources():
    return utils.resources.list(
        path=get_resources_dir(),
        glob_filename_pattern="*.json",
        only_basename=True,
        splitext=True,
    )


def init(key="perfect_absorber/rgb_0_0_0"):
    """
    Returns the surface's properties from merlict's own library-resources.
    If `key` is followed by a pattern such as `key/rgb_R_G_B`, then the
    color will be set to the integer values R,G, and B.

    Parameters
    ----------
    key : str, optional
        The key of the surface in merlict's own library. Default is
        `perfect_absorber`.
    """
    RGB = "/" in key
    basic_key = os.path.dirname(key) if RGB else key
    path = os.path.join(get_resources_dir(), basic_key + ".json")

    try:
        with open(path, "rt") as f:
            c = json_numpy.loads(f.read())
    except FileNotFoundError as e:
        print(
            "Unknown surface {:s}. Known surfaces are: {:s}".format(
                key, str(list_resources())
            )
        )
        raise e

    if RGB:
        rgb_key = os.path.basename(key)
        rgb = str.split(rgb_key, "_")
        assert rgb[0] == "rgb"
        assert len(rgb) == 4
        rgb = rgb[1:]
        c["color"] = [int(i) for i in rgb]
    return c
