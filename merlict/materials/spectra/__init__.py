import os
import numpy as np
from ... import utils
from ...scenery.string_format import function_csv


def init_from_resources(key):
    path = os.path.join(get_resources_path(), key + ".csv")

    with open(path, "rt") as f:
        x, y, x_label, y_label = function_csv.loads(f.read())
        spectrum = {"x": x, "y": y, "y_label": y_label, "x_label": x_label}

    assert_spectrum_is_valid(spectrum=spectrum)

    return spectrum


def get_resources_path():
    return utils.resources.path("materials", "spectra", "resources")


def list_resources():
    return utils.resources.list(
        path=get_resources_path(),
        glob_filename_pattern="*.csv",
        only_basename=True,
        splitext=True,
    )


def evaluate(spectrum, wavelength):
    assert wavelength >= min(spectrum["x"])
    assert wavelength <= max(spectrum["x"])
    return np.interp(x=wavelength, xp=spectrum["x"], fp=spectrum["y"])


def assert_spectrum_is_valid(spectrum, ymin=None, ymax=None):
    assert len(spectrum["x"]) >= 2, f"len(x) is {len(spectrum['x']):d}."
    assert len(spectrum["x"]) == len(spectrum["y"])

    assert utils.is_free_of_nan(spectrum["x"]), f"x :{str(spectrum['x']):s}."
    assert utils.is_all_greater_zero(spectrum["x"])
    assert utils.is_monotonically_increasing(spectrum["x"])

    assert utils.is_free_of_nan(spectrum["y"])

    if ymin is not None:
        assert np.all(spectrum["y"] >= ymin)

    if ymax is not None:
        assert np.all(spectrum["y"] <= ymax)

    assert "\n" not in spectrum["x_label"]
    assert "\n" not in spectrum["y_label"]
