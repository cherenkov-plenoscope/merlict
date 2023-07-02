import numpy as np
cimport numpy as np
cimport cython

cdef extern from "subset.h":
    cdef struct mliVec:
        double x
        double y
        double z

    cdef struct mliArchive:
        pass

    cdef mliArchive mliArchive_init()
    cdef void mliArchive_free(mliArchive *arc)
    cdef int mliArchive_malloc_from_path(mliArchive *arc, const char *path)

    cdef struct mliScenery:
        pass

    cdef mliScenery mliScenery_init()
    cdef void mliScenery_free(mliScenery *scn)
    cdef int mliScenery_malloc_from_Archive(mliScenery *scn, const mliArchive *arc)

    cdef struct mliPrng:
        pass

    cdef mliPrng mliPrng_init_MT19937(const unsigned int seed)
    cdef mliPrng mliPrng_init_PCG32(const unsigned int seed)
    cdef unsigned int mliPrng_generate_uint32(mliPrng *prng)
    cdef double mli_random_uniform(mliPrng *prng)


def _mliVec2py(mliVec mliv):
    return np.array([mliv.x, mliv.y, mliv.z], dtype=np.float64)


def _mliVec(v):
    cdef mliVec mliv
    mliv.x = v[0]
    mliv.y = v[1]
    mliv.z = v[2]
    return mliv


cdef class Prng:
    cdef mliPrng prng

    def __init__(self, seed, engine="PCG32"):
        cdef unsigned int cseed = np.uint32(seed)
        if engine == "PCG32":
            self.prng = mliPrng_init_PCG32(cseed)
        elif engine == "MT19937":
            self.prng = mliPrng_init_MT19937(cseed)
        else:
            raise KeyError("No such engine.")

    def uint32(self):
        cdef unsigned int c = mliPrng_generate_uint32(&self.prng)
        return int(c)


cdef class Scenery:
    cdef mliScenery scenery

    def __init__(self, path):
        self.scenery = mliScenery_init()

        cdef mliArchive archive
        archive = mliArchive_init()

        cdef bytes py_bytes = path.encode()
        cdef char* cpath = py_bytes

        cdef int rc
        rc = mliArchive_malloc_from_path(&archive, cpath)
        assert rc != 0

        rc = mliScenery_malloc_from_Archive(&self.scenery, &archive)
        assert rc != 0

        mliArchive_free(&archive)

    def __exit__(self):
        mliScenery_free(&self.scenery)
