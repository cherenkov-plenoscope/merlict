import merlict
import numpy as np
import pkg_resources
import os


SCENERY_PATH = pkg_resources.resource_filename(
    package_or_requirement="merlict",
    resource_name=os.path.join("test", "resources", "segmented_reflector.tar"),
)


def test_intersection_size_zero():
    scenery = merlict.c89.wrapper.Server(path=SCENERY_PATH)
    rays = merlict.ray.zeros(size=0)

    assert rays.shape[0] == 0

    isecs_valid, isecs = scenery.query_intersection(rays)

    assert isecs_valid.shape[0] == 0
    assert isecs_valid.dtype == np.bool
    assert isecs.shape[0] == 0


def test_intersection_simple():
    prng = np.random.Generator(np.random.PCG64(13))
    scenery = merlict.c89.wrapper.Server(path=SCENERY_PATH)
    N = 10 * 1000
    rays = merlict.ray.zeros(size=N)

    rays["support.x"] = prng.uniform(low=-0.5, high=0.5, size=N)
    rays["support.y"] = prng.uniform(low=-0.5, high=0.5, size=N)
    rays["support.z"] = 1.0

    rays["direction.x"] = 0.0
    rays["direction.y"] = 0.0
    rays["direction.z"] = -1.0

    assert rays.shape[0] == N

    mask, isecs = scenery.query_intersection(rays)

    assert mask.shape[0] == N
    assert isecs.shape[0] == N

    assert 0.3 < np.sum(mask) / N < 0.4

    valid_isecs = isecs[mask]
    assert 0.98 < np.average(valid_isecs["distance_of_ray"]) < 1.0
