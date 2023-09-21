import numpy as np
from .. import utils


def dtype():
    return [
        ("support.x", np.float64),
        ("support.y", np.float64),
        ("support.z", np.float64),
        ("direction.x", np.float64),
        ("direction.y", np.float64),
        ("direction.z", np.float64),
    ]


def israys(r):
    return utils.recarray.isdtype(r=r, dtype=dtype())


def init(size):
    return np.core.records.recarray(
        shape=size,
        dtype=dtype(),
    )


def zeros(size):
    return utils.recarray.zeros(size=size, dtype=dtype())


def frombytes(s):
    return np.core.records.fromstring(s, dtype=dtype())


def tobytes(rays):
    return rays.tobytes(order="C")


def fromphotons(photons):
    rays = init(size=photons.shape[0])
    for ray_key in rays.dtype.names:
        photon_key = "ray.{:s}".format(ray_key)
        rays[ray_key] = photons[photon_key]
    return rays
