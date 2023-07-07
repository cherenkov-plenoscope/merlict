from wrapper cimport *
cimport numpy as np
cimport cython
from libc cimport stdint

import termios
import sys
import numpy as np
from .. import rays as _rays
from .. import intersections as _intersections


def _mliVec2py(mliVec mliv):
    return np.array([mliv.x, mliv.y, mliv.z], dtype=np.float64)


def _mliVec(v):
    cdef mliVec mliv
    mliv.x = v[0]
    mliv.y = v[1]
    mliv.z = v[2]
    return mliv


def _mliView(position, rotation, field_of_view):
    cdef mliView view
    view.position = _mliVec(position)
    view.rotation = _mliVec(rotation)
    view.field_of_view = field_of_view
    return view


def _mlivrConfig_from_dict(config):
    c = config
    cdef mlivrConfig _c
    _c.random_seed = c["random_seed"]
    _c.preview_num_cols = c["preview_num_cols"]
    _c.preview_num_rows = c["preview_num_rows"]
    _c.export_num_cols = c["export_num_cols"]
    _c.export_num_rows = c["export_num_rows"]
    _c.step_length = c["step_length"]
    _c.view = _mliView(
        position=c["view"]["position"],
        rotation=c["view"]["rotation"],
        field_of_view=c["view"]["field_of_view"],
    )
    _c.aperture_camera_f_stop_ratio = c["aperture_camera_f_stop_ratio"]
    _c.aperture_camera_image_sensor_width = c["aperture_camera_image_sensor_width"]
    return _c


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

"""
======================================
"""


cdef class Scenery:
    cdef mliScenery scenery

    def __cinit__(self):
        self.scenery = mliScenery_init()

    def __dealloc__(self):
        mliScenery_free(&self.scenery)

    def __init__(self, path=None, sceneryPy=None):
        if path and not sceneryPy:
            self.init_from_path(path)
        elif sceneryPy and not path:
            self.init_from_sceneryPy(sceneryPy)
        else:
            raise ValueError("Either 'path' or 'sceneryPy', but not both.")

    def view(self, config=None):
        """
        Viewer will print to stdout.
        """
        cdef mlivrConfig viewer_config
        if config:
            viewer_config = _mlivrConfig_from_dict(config)
        else:
            viewer_config = mlivrConfig_default()

        fd = sys.stdin.fileno()
        old = termios.tcgetattr(fd)
        new = termios.tcgetattr(fd)
        C_FLAG = 3

        new[C_FLAG] = new[C_FLAG] & ~termios.ICANON
        try:
            termios.tcsetattr(fd, termios.TCSADRAIN, new)

            # ----
            rc = mlivr_run_interactive_viewer(
                &self.scenery,
                viewer_config
            )
            # ----

        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old)

    def init_from_path(self, path):
        cdef int rc
        _path = str(path)
        cdef bytes _py_path = _path.encode()
        cdef char* _cpath = _py_path  # Who is responsible for this memory?

        cdef mliArchive archive = mliArchive_init()
        try:
            rc = mliArchive_malloc_from_path(&archive, _cpath)
            assert rc != 0
            rc = mliScenery_malloc_from_Archive(&self.scenery, &archive)
            assert rc != 0
        finally:
            mliArchive_free(&archive)

    def init_from_dump_in_path(self, path):
        cdef int rc
        _path = str(path)
        cdef bytes _py_path = _path.encode()
        cdef char* _cpath = _py_path
        rc = mliScenery_malloc_from_path(&self.scenery, _cpath)
        assert rc != 0

    def dump_to_path(self, path):
        cdef int rc
        _path = str(path)
        cdef bytes _py_path = _path.encode()
        cdef char* _cpath = _py_path
        rc = mliScenery_write_to_path(&self.scenery, _cpath)
        assert rc != 0

    def init_from_sceneryPy(self, sceneryPy):
        cdef int rc
        cdef mliArchive archive = mliArchive_init()
        try:
            rc = mliArchive_malloc(&archive)
            assert rc != 0

        finally:
            mliArchive_free(&archive)

    def query_intersection(self, rays):
        cdef int rc
        assert _rays.is_rays(records=rays)

        isecs = _intersections.init(size=rays.shape[0])

        cdef stdint.uint64_t num_rays = rays.shape[0]

        cdef np.ndarray[mliRay,mode="c"] crays = np.ascontiguousarray(
            rays
        )

        cdef np.ndarray[mliIntersection,mode="c"] cisecs = np.ascontiguousarray(
            isecs
        )

        cdef np.ndarray[stdint.int64_t,mode="c"] cis_valid_isecs = np.ascontiguousarray(
            np.zeros(rays.shape[0], dtype=np.int64)
        )

        rc = mliBridge_query_many_intersection(
            &self.scenery,
            num_rays,
            &crays[0],
            &cisecs[0],
            &cis_valid_isecs[0])

        assert rc == 1

        return cis_valid_isecs, isecs
