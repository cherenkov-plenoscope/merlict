from .version import __version__
from . import c89
from . import materials
from . import scenery
from . import ray
from . import photon
from . import intersection
from . import intersectionSurfaceNormal
from . import viewer

import numpy as np


def photon_interactions(size):
    return np.core.records.recarray(
        shape=size,
        dtype=[
            ("on_geometry_surface", np.int32),
            ("geometry_id.robj", np.uint32),
            ("geometry_id.face", np.uint32),
            ("position.x", np.float64),
            ("position.y", np.float64),
            ("position.z", np.float64),
            ("position_local.x", np.float64),
            ("position_local.y", np.float64),
            ("position_local.z", np.float64),
            ("distance_of_ray", np.float64),
            ("medium_coming_from", np.uint64),
            ("medium_going_to", np.uint64),
            ("from_outside_to_inside", np.int32),
            ("typ", np.int32),
        ],
    )


def intersection(size):
    return np.core.records.recarray(
        shape=size,
        dtype=[
            ("geometry_id.robj", np.uint32),
            ("geometry_id.face", np.uint32),
            ("position_local.x", np.float64),
            ("position_local.y", np.float64),
            ("position_local.z", np.float64),
            ("distance_of_ray", np.float64),
        ],
    )


def IntersectionSurfaceNormal(size):
    return np.core.records.recarray(
        shape=size,
        dtype=[
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
            ("from_outside_to_inside", np.int32),
        ],
    )
