import numpy as np


def dtype():
    return [
        ("support.x", np.float64),
        ("support.y", np.float64),
        ("support.z", np.float64),
        ("direction.x", np.float64),
        ("direction.y", np.float64),
        ("direction.z", np.float64),
    ]


def is_rays(records):
    if not isinstance(records, np.core.records.recarray):
        return False
    for expected_name, expected_dtype in dtype():
        if expected_name not in records.dtype.names:
            return False
        if not records.dtype[expected_name] == expected_dtype:
            return False
    return True


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


def fromphotons(photons):
    rays = init(size=photons.shape[0])
    for ray_key in rays.dtype.names:
        photon_key = "ray.{:s}".format(ray_key)
        rays[ray_key] = photons[photon_key]
    return rays
