import os
import json_numpy
from ... import utils
from .. import spectra


def get_resources_path():
    return utils.resources.path("materials", "media", "resources")


def list_resources():
    return utils.resources.list(
        path=get_resources_path(),
        glob_filename_pattern="*.json",
        only_basename=True,
        splitext=True,
    )


def init_from_resources(key):
    path = os.path.join(get_resources_path(), key + ".json")
    with open(path, "rt") as f:
        c = json_numpy.loads(f.read())
    return c


def add_to_materials_from_resources(materials, key):
    medium = init_from_resources(key=key)
    speckeys = [medium["refraction_spectrum"], medium["absorption_spectrum"]]
    for spckey in speckeys:
        if spckey not in materials["spectra"]:
            spc = spectra.init_from_resources(key=spckey)
            materials["spectra"][spckey] = spc
    materials["media"][key] = medium
