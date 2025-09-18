import json_numpy as jsonp
import triangle_mesh_io as tmi
import posixpath
from . import function_csv


def sceneryPy_to_sceneryStr(sceneryPy, indent=4, relations_indent=0):
    """
    Returns a sceneryStr.

    The `Ds` is a dict() with all its values being of type str().
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
    """
    join = posixpath.join
    sceneryStr = []
    sceneryStr.append(("README.md", sceneryPy["readme"]))

    for key in sceneryPy["geometry"]["objects"]:
        _path = join("geometry", "objects", f"{key:s}.obj")
        _payload = tmi.obj.dumps(sceneryPy["geometry"]["objects"][key])
        sceneryStr.append((_path, _payload))

    _path = join("geometry", "relations.json")
    _payload = jsonp.dumps(
        sceneryPy["geometry"]["relations"], indent=relations_indent
    )
    sceneryStr.append((_path, _payload))

    for key in sceneryPy["materials"]["spectra"]:
        _path = join("materials", "spectra", f"{key:s}.csv")
        _payload = function_csv.dumps(**sceneryPy["materials"]["spectra"][key])
        sceneryStr.append((_path, _payload))

    for key in sceneryPy["materials"]["media"]:
        _path = join("materials", "media", f"{key:s}.json")
        _payload = jsonp.dumps(
            sceneryPy["materials"]["media"][key], indent=indent
        )
        sceneryStr.append((_path, _payload))

    for key in sceneryPy["materials"]["surfaces"]:
        _path = join("materials", "surfaces", f"{key:s}.json")
        _payload = jsonp.dumps(
            sceneryPy["materials"]["surfaces"][key], indent=indent
        )
        sceneryStr.append((_path, _payload))

    for key in sceneryPy["materials"]["boundary_layers"]:
        _path = join("materials", "boundary_layers", f"{key:s}.json")
        _payload = jsonp.dumps(
            sceneryPy["materials"]["boundary_layers"][key], indent=indent
        )
        sceneryStr.append((_path, _payload))

    _path = join("materials", "default_medium.txt")
    _payload = str(sceneryPy["materials"]["default_medium"])
    sceneryStr.append((_path, _payload))

    return sceneryStr


def sceneryStr_to_sceneryPy(sceneryStr):
    sceneryPy = {}
    sceneryPy["geometry"] = {}
    sceneryPy["geometry"]["objects"] = {}
    sceneryPy["materials"] = {}
    sceneryPy["materials"]["default_medium"] = None
    sceneryPy["materials"]["spectra"] = {}
    sceneryPy["materials"]["media"] = {}
    sceneryPy["materials"]["surfaces"] = {}
    sceneryPy["materials"]["boundary_layers"] = {}
    join = posixpath.join

    for item in sceneryStr:
        _path, _load = item
        if "README" in _path:
            sceneryPy["readme"] = _load

        elif join("geometry", "objects") in _path and ".obj" in _path:
            okey = _posixpath_basename_without_extension(_path)
            sceneryPy["geometry"]["objects"][okey] = tmi.obj.loads(_load)

        elif join("geometry", "relations.json") == _path:
            sceneryPy["geometry"]["relations"] = jsonp.loads(_load)

        elif join("materials", "spectra") in _path and ".csv" in _path:
            key = _posixpath_basename_without_extension(_path)
            sceneryPy["materials"]["spectra"][key] = function_csv.loads(_load)

        elif join("materials", "media") in _path and ".json" in _path:
            key = _posixpath_basename_without_extension(_path)
            sceneryPy["materials"]["media"][key] = jsonp.loads(_load)

        elif join("materials", "surfaces") in _path and ".json" in _path:
            key = _posixpath_basename_without_extension(_path)
            sceneryPy["materials"]["surfaces"][key] = jsonp.loads(_load)

        elif (
            join("materials", "boundary_layers") in _path and ".json" in _path
        ):
            key = _posixpath_basename_without_extension(_path)
            sceneryPy["materials"]["boundary_layers"][key] = jsonp.loads(_load)

        elif join("materials", "default_medium.txt") == _path:
            sceneryPy["materials"]["default_medium"] = str(_load)

    return sceneryPy


def _posixpath_basename_without_extension(path):
    path = posixpath.basename(path)
    return posixpath.splitext(path)[0]
