import json_numpy as jsonp
import triangle_mesh_io as tmi
import posixpath


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

    for okey in sceneryPy["geometry"]["objects"]:
        filepath = join("geometry", "objects", "{:s}.obj".format(okey))
        payload = tmi.obj.dumps(sceneryPy["geometry"]["objects"][okey])
        sceneryStr.append((filepath, payload))

    filepath = join("geometry", "relations.json")
    payload = jsonp.dumps(
        sceneryPy["geometry"]["relations"], indent=relations_indent
    )
    sceneryStr.append((filepath, payload))

    for mkey in sceneryPy["materials"]["media"]:
        filepath = join("materials", "media", "{:s}.json".format(mkey))
        payload = jsonp.dumps(
            sceneryPy["materials"]["media"][mkey], indent=indent
        )
        sceneryStr.append((filepath, payload))

    for skey in sceneryPy["materials"]["surfaces"]:
        filepath = join("materials", "surfaces", "{:s}.json".format(skey))
        payload = jsonp.dumps(
            sceneryPy["materials"]["surfaces"][skey], indent=indent
        )
        sceneryStr.append((filepath, payload))

    filepath = join("materials", "boundary_layers.json")
    payload = jsonp.dumps(
        sceneryPy["materials"]["boundary_layers"], indent=relations_indent
    )
    sceneryStr.append((filepath, payload))

    filepath = join("materials", "default_medium.txt")
    payload = str(sceneryPy["materials"]["default_medium"])
    sceneryStr.append((filepath, payload))

    return sceneryStr


def sceneryStr_to_sceneryPy(sceneryStr):
    sceneryPy = {}
    sceneryPy["geometry"] = {}
    sceneryPy["geometry"]["objects"] = {}
    sceneryPy["materials"] = {}
    sceneryPy["materials"]["media"] = {}
    sceneryPy["materials"]["surfaces"] = {}
    join = posixpath.join

    for item in sceneryStr:
        filepath, payload = item
        if "README" in filepath:
            sceneryPy["readme"] = payload
        elif join("geometry", "objects") in filepath and ".obj" in filepath:
            okey = _posixpath_basename_without_extension(filepath)
            sceneryPy["geometry"]["objects"][okey] = tmi.obj.loads(payload)
        elif join("geometry", "relations.json") == filepath:
            sceneryPy["geometry"]["relations"] = jsonp.loads(payload)
        elif join("materials", "media") in filepath and ".json" in filepath:
            mkey = _posixpath_basename_without_extension(filepath)
            sceneryPy["materials"]["media"][mkey] = jsonp.loads(payload)
        elif join("materials", "surfaces") in filepath and ".json" in filepath:
            mkey = _posixpath_basename_without_extension(filepath)
            sceneryPy["materials"]["surfaces"][mkey] = jsonp.loads(payload)
        elif join("materials", "default_medium.txt") == filepath:
            sceneryPy["materials"]["default_medium"] = str(payload)
    return sceneryPy


def _posixpath_basename_without_extension(path):
    path = posixpath.basename(path)
    return posixpath.splitext(path)[0]
