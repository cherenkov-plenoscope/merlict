import numpy as np
import tarfile
import io
import json_numpy
from .. import materials


def init(default_medium="vacuum"):
    """
    Returns a minimal sceneryPy without any objects in it.

    Parameters
    ----------
    default_medium : str
        The key of a medium in merlict's own library of media to be used for
        the default medium in between surfaces.
    """
    scenery = {
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
    scenery["materials"]["media"][default_medium] = materials.media.init(
        key=default_medium
    )
    return scenery


def convert_sceneryPy_to_sceneryDs(
    sceneryPy, indent=4, relations_indent=0, readme=None
):
    """
    Returns a sceneryDs.

    The `DS` is a dict() with all its values being of type str().
    Materials will be dumped into json-strings.
    Objects will be dumped into obj-strings.
    The key of the default_medium will be dumped into a plain string.
    These value-strings correspond to the payloads found in sceneryTar.

    Parameters
    ----------
    sceneryPy : dict
        The scenery set up by the user using dicts, lists, numpy.arrays and so
        on.
    indent : int
        Number of chars to indent in json-files. Default is 4.
    relations_indent : int
        Number of chars to indent in geometry/relations.json.
        This file can be pretty large, so the default is 0 to save space.
        One can also set `None` what avoids all linebreaks but makes reading
        the file almost impossible for humans.
    readme : str
        An optional text for the README.md
    """
    sceneryDs = {}

    sceneryDs["README.md"] = readme if readme else _default_readme()

    sceneryDs["geometry"] = {}
    sceneryDs["geometry"]["objects"] = {}
    for okey in sceneryPy["geometry"]["objects"]:
        raise NotImplementedError(
            "On my list. For this we need to break up my python-package "
            "optics_obect_wavefronts."
        )
        okey_filename = "{:s}.obj".format(okey)
        wfr = obj.init_from_mesh(mesh=sceneryPy["geometry"]["objects"][okey])
        sceneryDs["geometry"]["objects"][okey_filename] = obj.dumps(obj=wfr)

    sceneryDs["geometry"]["relations.json"] = json_numpy.dumps(
        sceneryPy["geometry"]["relations"], indent=relations_indent
    )

    sceneryDs["materials"] = {}
    sceneryDs["materials"]["media"] = {}
    for mkey in sceneryPy["materials"]["media"]:
        mkey_filename = "{:s}.json".format(mkey)
        mjsn = json_numpy.dumps(
            sceneryPy["materials"]["media"][mkey], indent=indent
        )
        sceneryDs["materials"]["media"][mkey_filename] = mjsn

    sceneryDs["materials"]["surfaces"] = {}
    for skey in sceneryPy["materials"]["surfaces"]:
        skey_filename = "{:s}.json".format(skey)
        sjsn = json_numpy.dumps(
            sceneryPy["materials"]["surfaces"][skey], indent=indent
        )
        sceneryDs["materials"]["media"][mkey_filename] = sjsn

    sceneryDs["materials"]["boundary_layers.json"] = json_numpy.dumps(
        sceneryPy["materials"]["boundary_layers"],
        indent=indent,
    )

    sceneryDs["materials"]["default_medium.txt"] = str(
        sceneryPy["materials"]["default_medium"]
    )

    return sceneryDs


def _default_readme():
    rm = "Scenery\n"
    rm += "=======\n"
    rm += "I was written by merlict but nobody botherd to provide a README\n"
    rm += "so this is the default :sad:. Anyhow, I am a scenery of objects.\n"
    rm += "Merlict can read me and perform ray tracing on me.\n"
    return rm
