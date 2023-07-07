import numpy as np


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


def init(size):
    return np.core.records.recarray(shape=size, dtype=dtype(),)


def zeros(size):
    out = init(size=size)
    for key in out.dtype.names:
        out[key] = np.zeros(size, dtype=out.dtype[key])
    return out


def frombytes(s):
    return np.core.records.fromstring(s, dtype=dtype())


def tobytes(rays):
    return rays.tobytes(order="C")
