import numpy as np
cimport numpy as np
cimport cython

cdef extern from "mli_subset.h":
    cdef struct mliArchive:
        pass

    cdef struct mliVec:
        double x
        double y
        double z


def _mliVec2py(mliVec mliv):
    return np.array([mliv.x, mliv.y, mliv.z], dtype=np.float64)


def _mliVec(v):
    cdef mliVec mliv
    mliv.x = v[0]
    mliv.y = v[1]
    mliv.z = v[2]
    return mliv
