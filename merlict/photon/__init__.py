import numpy as np
from .. import utils


def dtype():
    return [
        ("id", np.int64),
        ("ray.support.x", np.float64),
        ("ray.support.y", np.float64),
        ("ray.support.z", np.float64),
        ("ray.direction.x", np.float64),
        ("ray.direction.y", np.float64),
        ("ray.direction.z", np.float64),
        ("wavelength", np.float64),
    ]


def isphotons(r):
    return utils.recarray.isdtype(r=r, dtype=dtype())


def init(size):
    return np.core.records.recarray(shape=size, dtype=dtype())


def zeros(size):
    return utils.recarray.zeros(size=size, dtype=dtype())


def frombytes(s):
    return np.core.records.fromstring(s, dtype=dtype())


def tobytes(rays):
    return rays.tobytes(order="C")
