import os
import glob
import posixpath
from . import fileorder


def write(sceneryStr, path):
    os.makedirs(path, exist_ok=True)
    for item in sceneryStr:
        relfilepath_posix, payload = item
        relfilepath = _posix2os(relfilepath_posix)
        relfiledirname = posixpath.dirname(relfilepath)
        if relfiledirname:
            dirname = os.path.join(path, relfiledirname)
            os.makedirs(dirname, exist_ok=True)
        filepath = os.path.join(path, relfilepath)
        with open(filepath, "wt") as f:
            f.write(payload)


def read(path):
    sceneryDS = []
    for relfilepath_posix in fileorder.list():
        relfilepath = _posix2os(relfilepath_posix)
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


def _read_text_file(path):
    with open(path, "rt") as f:
        payload = f.read()
    return payload


def _posix2os(path):
    return str.replace(path, posixpath.sep, os.sep)
