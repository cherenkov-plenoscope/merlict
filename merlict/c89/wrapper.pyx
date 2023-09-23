from .wrapper cimport *
cimport numpy as np
from libc cimport stdint

import termios
import sys
import numpy as np
from .. import ray as _ray
from .. import intersection as _intersection
from .. import intersectionSurfaceNormal as _intersectionSurfaceNormal


cdef _mliVec2py(mliVec mliv):
    return np.array([mliv.x, mliv.y, mliv.z], dtype=np.float64)


cdef _mliVec(v):
    cdef mliVec mliv
    mliv.x = v[0]
    mliv.y = v[1]
    mliv.z = v[2]
    return mliv


cdef _mliImage2py(mliImage mliimg):
    out = np.zeros(
        shape=(mliimg.num_cols, mliimg.num_rows, 3),
        dtype=np.float32)
    cdef mliColor c
    for ix in range(mliimg.num_cols):
        for iy in range(mliimg.num_rows):
            c = mliImage_at(&mliimg, ix, iy)
            out[ix, iy, 0] = c.r
            out[ix, iy, 1] = c.g
            out[ix, iy, 2] = c.b
    return out


cdef _mliView(position, rotation, field_of_view):
    cdef mliView view
    view.position = _mliVec(position)
    view.rotation = _mliVec(rotation)
    view.field_of_view = field_of_view
    return view


cdef _mlivrConfig_from_dict(config):
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
    _c.aperture_camera_image_sensor_width = c[
        "aperture_camera_image_sensor_width"
    ]
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


cdef class Server:
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
        Prints a rendered image of the scenery to stdout. Waits for user-input
        (the user pressing keys) to manipulate the view.

        Press [h] to print the help.
        """
        cdef mlivrConfig cconfig
        if config:
            cconfig = _mlivrConfig_from_dict(config)
        else:
            cconfig = mlivrConfig_default()

        fd = sys.stdin.fileno()
        old_attr = termios.tcgetattr(fd)
        new_attr = termios.tcgetattr(fd)
        C_FLAG = 3

        new_attr[C_FLAG] = new_attr[C_FLAG] & ~termios.ICANON
        try:
            termios.tcsetattr(fd, termios.TCSADRAIN, new_attr)
            _rc = mlivr_run_interactive_viewer(&self.scenery, cconfig)
        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old_attr)

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
        assert _ray.israys(rays)
        cdef stdint.uint64_t num_ray = rays.shape[0]
        isecs = _intersection.init(size=num_ray)

        cdef np.ndarray[mliRay, mode="c"] crays = np.ascontiguousarray(
            rays
        )

        cdef np.ndarray[
            mliIntersection, mode="c"
        ] cisecs = np.ascontiguousarray(
            isecs
        )

        cdef np.ndarray[
            stdint.int64_t, mode="c"
        ] cis_valid_isecs = np.ascontiguousarray(
            np.zeros(rays.shape[0], dtype=np.int64)
        )

        if num_ray:
            rc = mliBridge_query_many_intersection(
                &self.scenery,
                num_ray,
                &crays[0],
                &cisecs[0],
                &cis_valid_isecs[0])

            assert rc == 1

        isecs_mask = cis_valid_isecs.astype(dtype=np.bool_)

        return isecs_mask, isecs

    def query_intersectionSurfaceNormal(self, rays):
        cdef int rc
        assert _ray.israys(rays)
        cdef stdint.uint64_t num_ray = rays.shape[0]
        isecs = _intersectionSurfaceNormal.init(size=num_ray)

        cdef np.ndarray[mliRay, mode="c"] crays = np.ascontiguousarray(
            rays
        )

        cdef np.ndarray[
            mliIntersectionSurfaceNormal, mode="c"
        ] cisecs = np.ascontiguousarray(
            isecs
        )

        cdef np.ndarray[
            stdint.int64_t, mode="c"
        ] cis_valid_isecs = np.ascontiguousarray(
            np.zeros(rays.shape[0], dtype=np.int64)
        )

        if num_ray:
            rc = mliBridge_query_many_intersectionSurfaceNormal(
                &self.scenery,
                num_ray,
                &crays[0],
                &cisecs[0],
                &cis_valid_isecs[0])

            assert rc == 1

        isecs_mask = cis_valid_isecs.astype(dtype=np.bool_)

        return isecs_mask, isecs

    def __repr__(self):
        out = "{:s}()".format(self.__class__.__name__)
        return out
