import pkg_resources
import os
import json_numpy


def get_materials_dir():
    return pkg_resources.resource_filename("merlict", "materials")


def surface(key):
    """
    Returns the surface's properties from merlict's own library of surfaces.
    """
    MATERIALS_DIR = get_materials_dir()
    SURFACES_DIR = os.path.join(MATERIALS_DIR, "surfaces")

    RGB = "/" in key

    basic_key = os.path.dirname(key) if RGB else key

    path = os.path.join(SURFACES_DIR, basic_key + ".json")
    with open(path, "rt") as f:
        c = json_numpy.loads(f.read())

    if RGB:
        rgb_key = os.path.basename(key)
        rgb = str.split(rgb_key, "_")
        assert rgb[0] == "rgb"
        assert len(rgb) == 4
        rgb = rgb[1:]
        c["color"] = [int(i) for i in rgb]
    return c


def medium(key):
    """
    Returns the medium's properties from merlict's own library of media.
    """
    MATERIALS_DIR = get_materials_dir()
    MEDIA_DIR = os.path.join(MATERIALS_DIR, "media")

    path = os.path.join(MEDIA_DIR, key + ".json")
    with open(path, "rt") as f:
        c = json_numpy.loads(f.read())
    return c
