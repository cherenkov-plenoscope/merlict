import os
import json
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


def add_to_materials_from_resources(materials, key):
    path = os.path.join(get_resources_path(), key + ".json")
    with open(path, "rt") as f:
        medium = json.loads(f.read())
    speckeys = [medium["refraction_spectrum"], medium["absorption_spectrum"]]
    for spckey in speckeys:
        if spckey not in materials["spectra"]:
            spc = spectra.init_from_resources(key=spckey)
            materials["spectra"][spckey] = spc
    materials["media"][key] = medium
