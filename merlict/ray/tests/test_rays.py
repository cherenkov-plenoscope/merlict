import merlict


def test_init():
    rays = merlict.ray.zeros(size=13)
    assert rays.shape[0] == 13


def test_bytes():
    rays = merlict.ray.zeros(size=13)
    raysb = merlict.ray.tobytes(rays)
    rays_back = merlict.ray.frombytes(raysb)

    for i in range(rays.shape[0]):
        for key in rays.dtype.names:
            assert rays_back[key][i] == rays[key][i]


def test_fromphotons():
    photons = merlict.photon.zeros(size=43)
    rays = merlict.ray.fromphotons(photons)

    for i in range(rays.shape[0]):
        for rays_key in rays.dtype.names:
            photons_key = "ray.{:s}".format(rays_key)
            assert photons[photons_key][i] == rays[rays_key][i]
