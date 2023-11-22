from .. import materials
from . import string_format
import os
import numpy as np


def init(default_medium="vacuum"):
    """
    Returns a minimal sceneryPy without any objects in it.

    Parameters
    ----------
    default_medium : str
        The key of a medium in merlict's own library of media to be used for
        the default medium in between surfaces.
    """
    sceneryPy = {
        "readme": _default_readme(),
        "materials": {
            "media": {},
            "surfaces": {},
            "boundary_layers": {},
            "default_medium": default_medium,
        },
        "geometry": {
            "objects": {},
            "relations": {"children": []},
        },
    }
    sceneryPy["materials"]["media"][default_medium] = materials.media.init(
        key=default_medium
    )
    return sceneryPy


def init_from_object(obj, default_medium="vacuum", random_seed=1):
    sceneryPy = init(default_medium=default_medium)
    prng = np.random.Generator(np.random.PCG64(random_seed))

    sceneryPy["geometry"]["objects"]["one"] = obj

    mtlkeys = list(obj["mtl"].keys())

    mtl_to_boundary_layers = {}
    for icolor, mtlkey in enumerate(mtlkeys):
        r, g, b = prng.uniform(low=64, high=192, size=3).astype(int)
        surface = materials.surfaces.init(
            "perfect_absorber/rgb_{:d}_{:d}_{:d}".format(r, g, b)
        )
        surfkey = "color_{:d}".format(icolor)
        sceneryPy["materials"]["surfaces"][surfkey] = surface

        boundkey = "bound_{:d}".format(icolor)
        sceneryPy["materials"]["boundary_layers"][boundkey] = {
            "inner": {"medium": default_medium, "surface": surfkey},
            "outer": {"medium": default_medium, "surface": surfkey},
        }
        mtl_to_boundary_layers[mtlkey] = boundkey

    sceneryPy["geometry"]["relations"]["children"].append(
        {
            "id": 1,
            "pos": [0, 0, 0],
            "rot": {"repr": "tait_bryan", "xyz_deg": [0, 0, 0]},
            "obj": "one",
            "mtl": mtl_to_boundary_layers,
        }
    )
    return sceneryPy


def write_dir(sceneryPy, path):
    sceneryStr = string_format.convert.sceneryPy_to_sceneryStr(sceneryPy)
    string_format.directory.write(sceneryStr=sceneryStr, path=path)


def write_tar(sceneryPy, path):
    sceneryStr = string_format.convert.sceneryPy_to_sceneryStr(sceneryPy)
    string_format.tapearchive.write(sceneryStr=sceneryStr, path=path)


def read_dir(path):
    sceneryStr = string_format.directory.read(path=path)
    return string_format.convert.sceneryStr_to_sceneryPy(sceneryStr)


def read_tar(path):
    sceneryStr = string_format.tapearchive.read(path=path)
    return string_format.convert.sceneryStr_to_sceneryPy(sceneryStr)


def _default_readme():
    rm = "Scenery\n"
    rm += "=======\n"
    rm += "I was written by merlict but nobody botherd to provide a README\n"
    rm += "so this is the default :sad:. Anyhow, I am a scenery of objects.\n"
    rm += "Merlict can read me and perform ray tracing on me.\n"
    return rm
