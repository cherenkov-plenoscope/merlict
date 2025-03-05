def dtype():
    return [
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
        ("type", np.int32),
    ]


def init(size):
    return np.recarray(shape=size, dtype=dtype())


def zeros(size):
    return utils.recarray.zeros(size=size, dtype=dtype())
