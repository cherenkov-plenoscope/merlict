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

    for okey in sceneryPy["geometry"]["objects"]:
        filepath = join("geometry", "objects", "{:s}.obj".format(okey))
        payload = tmi.obj.dumps(sceneryPy["geometry"]["objects"][okey])
        sceneryStr.append((filepath, payload))

    filepath = join("geometry", "relations.json")
    payload = jsonp.dumps(
        sceneryPy["geometry"]["relations"], indent=relations_indent
    )
    sceneryStr.append((filepath, payload))

    for spectra_key in sceneryPy["materials"]["spectra"]:
        spectra_dir = join("materials", "spectra", f"{spectra_key:s}")
        spectrum = sceneryPy["materials"]["spectra"][spectra_key]
        payload = function_csv.dumps(
            func=spectrum["spectrum"],
            header=("wavelength/m", "refraction/1"),
        )
        sceneryStr.append((join(media_dir, "refraction.csv"), payload))

    for media_key in sceneryPy["materials"]["media"]:
        media_dir = join("materials", "media", f"{media_key:s}")

        # refraction
        payload = function_csv.dumps(
            func=sceneryPy["materials"]["media"][media_key]["refraction"],
            header=("wavelength/m", "refraction/1"),
        )
        sceneryStr.append((join(media_dir, "refraction.csv"), payload))

        # absorbtion
        payload = function_csv.dumps(
            func=sceneryPy["materials"]["media"][media_key]["absorbtion"],
            header=("wavelength/m", "absorbtion/m^{-1}"),
        )
        sceneryStr.append((join(media_dir, "absorbtion.csv"), payload))

    for surface_key in sceneryPy["materials"]["surfaces"]:
        surface_dir = join("materials", "surfaces", f"{surface_key:s}")

        # specular_reflection
        payload = function_csv.dumps(
            func=sceneryPy["materials"]["surfaces"][surface_key][
                "specular_reflection"
            ],
            header=("wavelength/m", "reflection/1"),
        )
        sceneryStr.append(
            (join(surface_dir, "specular_reflection.csv"), payload)
        )

        # diffuse_reflection
        payload = function_csv.dumps(
            func=sceneryPy["materials"]["surfaces"][surface_key][
                "diffuse_reflection"
            ],
            header=("wavelength/m", "reflection/1"),
        )
        sceneryStr.append(
            (join(surface_dir, "diffuse_reflection.csv"), payload)
        )

        # material
        payload = sceneryPy["materials"]["surfaces"][surface_key]["material"]
        sceneryStr.append((join(surface_dir, "material.txt"), payload))

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
        elif join("materials", "boundary_layers.json") == filepath:
            sceneryPy["materials"]["boundary_layers"] = jsonp.loads(payload)
    return sceneryPy


def _posixpath_basename_without_extension(path):
    path = posixpath.basename(path)
    return posixpath.splitext(path)[0]
