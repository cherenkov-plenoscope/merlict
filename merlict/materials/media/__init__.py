import os
import json_numpy
from ... import utils


def get_resources_dir():
    return utils.resources.path("materials", "media", "resources")


def list_resources():
    return utils.resources.list(
        path=get_resources_dir(),
        glob_filename_pattern="*.json",
        only_basename=True,
        splitext=True,
    )


def init(key="vacuum"):
    """
    Returns the medium's properties from merlict's own library-resources.

    Parameters
    ----------
    key : str, optional
        The key of the medium in merlict's own library. Default is `vacuum`.
    """
    path = os.path.join(get_resources_dir(), key + ".json")
    try:
        with open(path, "rt") as f:
            c = json_numpy.loads(f.read())
    except FileNotFoundError as e:
        print(
            "Unknown medium {:s}. Known media are: {:s}".format(
                key, str(list_resources())
            )
        )
        raise e

    return c
