import numpy as np


def Spectrum(wavelength, value, wavelength_name="wavelength", value_name="value", comment=""):
    out = {}
    out["wavelength"] = np.asarray(wavelength, dtype=float)
    out["value"] = np.asarray(value, dtype=float)
    out["value_name"] = str(value_name)
    out["wavelength_name"] = str(wavelength_name)
    out["comment"] = str(comment)
    return out
