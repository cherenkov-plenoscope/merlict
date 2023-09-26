import numpy as np


def dtype():
    return [
        ("geometry_id.robj", np.uint32),
        ("geometry_id.face", np.uint32),
        ("position.x", np.float64),
        ("position.y", np.float64),
        ("position.z", np.float64),
        ("surface_normal.x", np.float64),
        ("surface_normal.y", np.float64),
        ("surface_normal.z", np.float64),
        ("position_local.x", np.float64),
        ("position_local.y", np.float64),
        ("position_local.z", np.float64),
        ("surface_normal_local.x", np.float64),
        ("surface_normal_local.y", np.float64),
        ("surface_normal_local.z", np.float64),
        ("distance_of_ray", np.float64),
        ("from_outside_to_inside", np.int64),
    ]


def init(size):
    return np.core.records.recarray(
        shape=size,
        dtype=dtype(),
    )
