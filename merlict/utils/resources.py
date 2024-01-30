import importlib
from importlib import resources
import os
import glob


def path(*args):
    pkg_path = str(importlib.resources.files("merlict"))
    return os.path.join(pkg_path, os.path.join(*args))


def list(path, glob_filename_pattern, only_basename=True, splitext=True):
    paths = glob.glob(os.path.join(path, glob_filename_pattern))
    if only_basename:
        paths = [os.path.basename(pp) for pp in paths]
    if splitext:
        paths = [os.path.splitext(pp)[0] for pp in paths]
    return paths
