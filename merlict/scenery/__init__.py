from .. import materials
from . import list_of_strings
import json_numpy
import triangle_mesh_io as tmi
import tarfile
import io
import os
import glob
import posixpath


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
    scenery["materials"]["media"][default_medium] = materials.media.init(
        key=default_medium
    )
    return scenery


def convert_sceneryPy_to_sceneryDs(sceneryPy, indent=4, relations_indent=0):
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
    """
    sceneryDs = []
    sceneryDs.append(("README.md", sceneryPy["readme"]))

    for okey in sceneryPy["geometry"]["objects"]:
        filepath = posixpath.join("geometry", "objects", "{:s}.obj".format(okey))
        payload = tmi.obj.dumps(sceneryPy["geometry"]["objects"][okey])
        sceneryDs.append((filepath, payload))

    filepath = posixpath.join("geometry", "relations.json")
    payload = json_numpy.dumps(
        sceneryPy["geometry"]["relations"], indent=relations_indent
    )
    sceneryDs.append((filepath, payload))

    for mkey in sceneryPy["materials"]["media"]:
        filepath = posixpath.join("materials", "media", "{:s}.json".format(mkey))
        payload = json_numpy.dumps(
            sceneryPy["materials"]["media"][mkey], indent=indent
        )
        sceneryDs.append((filepath, payload))

    for skey in sceneryPy["materials"]["surfaces"]:
        filepath = posixpath.join("materials", "surfaces", "{:s}.json".format(skey))
        payload = json_numpy.dumps(
            sceneryPy["materials"]["surfaces"][skey], indent=indent
        )
        sceneryDs.append((filepath, payload))

    filepath = posixpath.join("materials", "boundary_layers.json")
    payload = json_numpy.dumps(
        sceneryPy["materials"]["boundary_layers"], indent=relations_indent
    )
    sceneryDs.append((filepath, payload))

    filepath = posixpath.join("materials", "default_medium.txt")
    payload = str(sceneryPy["materials"]["default_medium"])
    sceneryDs.append((filepath, payload))

    return sceneryDs


def _default_readme():
    rm = "Scenery\n"
    rm += "=======\n"
    rm += "I was written by merlict but nobody botherd to provide a README\n"
    rm += "so this is the default :sad:. Anyhow, I am a scenery of objects.\n"
    rm += "Merlict can read me and perform ray tracing on me.\n"
    return rm


def write(sceneryPy, path, indent=4, relations_indent=0, readme=None):
    sceneryDs = convert_sceneryPy_to_sceneryDs(
        sceneryPy,
        indent=indent,
        relations_indent=relations_indent,
        readme=readme,
    )
    _write_sceneryDS_to_tar(sceneryDs=sceneryDs, path=path)


def _write_sceneryDS_to_tar(sceneryDs, path):
    with tarfile.open(name=path, mode="w|") as tar:
        for item in sceneryDs:
            filename, payload = item
            _tar_append_file(tar, filename, str.encode(payload))


def _tar_append_file(tar, filename, filebytes):
    with io.BytesIO() as buff:
        info = tarfile.TarInfo(filename)
        info.size = buff.write(filebytes)
        buff.seek(0)
        tar.addfile(info, buff)


def _write_sceneryDS_to_dir(sceneryDs, path):
    os.makedirs(path, exist_ok=True)
    for item in sceneryDs:
        relfilepath_posix, payload = item
        relfilepath = posix2os(relfilepath_posix)
        relfiledirname = posixpath.dirname(relfilepath)
        if relfiledirname:
            dirname = os.path.join(path, relfiledirname)
            os.makedirs(dirname, exist_ok=True)
        filepath = os.path.join(path, relfilepath)
        with open(filepath, "wt") as f:
            f.write(payload)


def _sceneryDS_file_order():
    return [
        "README.md",
        posixpath.join("geometry", "objects", "*.obj"),
        posixpath.join("geometry", "relations.json"),
        posixpath.join("materials", "media", "*.json"),
        posixpath.join("materials", "surfaces", "*.json"),
        posixpath.join("materials", "boundary_layers.json"),
        posixpath.join("materials", "default_medium.txt"),
    ]


def _read_sceneryDS_from_dir(path):
    sceneryDS = []
    for relfilepath_posix in _sceneryDS_file_order():
        relfilepath = posix2os(relfilepath_posix)
        filepath = os.path.join(path, relfilepath)

        if "*" in filepath:
            for ifilepath in glob.glob(filepath):
                ifilebasename = os.path.basename(ifilepath)
                ifilebasename_wo_ext = os.path.splitext(ifilebasename)[0]
                irelfilepath_posix = str.replace(
                    relfilepath_posix, "*", ifilebasename_wo_ext
                )

                payload = _read_text_file(ifilepath)
                item = (irelfilepath_posix, payload)
                sceneryDS.append(item)
        else:
            payload = _read_text_file(filepath)
            item = (relfilepath_posix, payload)
            sceneryDS.append(item)

    return sceneryDS


def convert_sceneryPs_to_sceneryDy(sceneryDs):
    sceneryPy = {}
    sceneryPy["geometry"] = {}
    sceneryPy["geometry"]["objects"] = {}
    sceneryPy["materials"] = {}
    sceneryPy["materials"]["media"] = {}
    sceneryPy["materials"]["surfaces"] = {}

    for item in sceneryDs:
        filepath, payload = item
        if "README" in filepath:
            sceneryPy["readme"] = payload
        elif "geometry/objects/" in filepath and ".obj" in filepath:
            okey = posixpath_basename_without_extension(filepath)
            sceneryPy["geometry"]["objects"][okey] = tmi.obj.loads(payload)
        elif "geometry/relations.json" == filepath:
            sceneryPy["geometry"]["relations"] = json_numpy.loads(payload)
        elif "materials/media/" in filepath and ".json" in filepath:
            mkey = posixpath_basename_without_extension(filepath)
            sceneryPy["materials"]["media"][mkey] = json_numpy.loads(payload)
        elif "materials/surfaces/" in filepath and ".json" in filepath:
            mkey = posixpath_basename_without_extension(filepath)
            sceneryPy["materials"]["surfaces"][mkey] = json_numpy.loads(payload)
        elif "materials/default_medium.txt" == filepath:
            sceneryPy["materials"]["default_medium"] = str(payload)
    return sceneryPy


def posixpath_basename_without_extension(path):
    path = posixpath.basename(path)
    return posixpath.splitext(path)[0]


def _read_text_file(path):
    with open(path, "rt") as f:
        payload = f.read()
    return payload


def os2posix(path):
    return str.replace(path, os.sep, posixpath.sep)

def posix2os(path):
    return str.replace(path, posixpath.sep, os.sep)


def is_gzipped_file(path):
    """
    Check for gzip file, see https://tools.ietf.org/html/rfc1952#page-5

    Reads in the first two bytes of a file and compares with the gzip magic
    numbers.
    """
    with open(path, 'rb') as fin:
        marker = fin.read(2)
        if len(marker) < 2:
            return False
        return marker[0] == 31 and marker[1] == 139
