import os
import json_numpy
from ... import utils


def get_resources_path():
    return utils.resources.path("materials", "surfaces", "resources")


def list_resources():
    return utils.resources.list(
        path=get_resources_path(),
        glob_filename_pattern="*.json",
        only_basename=True,
        splitext=True,
    )


def init_from_resources(key="perfect_absorber"):
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
    path = os.path.join(get_resources_path(), key + ".json")

    with open(path, "rt") as f:
        c = json_numpy.loads(f.read())

    return c
