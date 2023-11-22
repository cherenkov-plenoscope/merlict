from .wrapper cimport *

cimport numpy as cnumpy
cnumpy.import_array()

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


cdef _mliArchive_push_back_path_and_payload(
    mliArchive *archive,
    path,
    payload,
):
    cdef int rc

    _path = str(path)
    _payload = str(payload)

    cdef bytes _py_path = _path.encode()
    cdef bytes _py_payload = _payload.encode()

    cdef stdint.uint64_t path_length = np.uint64(len(_py_path))
    cdef stdint.uint64_t payload_length = len(_py_payload)

    cdef char* _cpath = _py_path
    cdef char* _cpayload = _py_payload

    rc = mliArchive_push_back_cstr(
        archive,
        _cpath,
        path_length,
        _cpayload,
        payload_length
    )
    assert rc != 0
    return


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

    def __init__(self, path=None, sceneryDs=None):
        if path and not sceneryDs:
            try:
                self.init_from_tar(path)
            except AssertionError:
                try:
                    self.init_from_dump(path)
                except AssertionError:
                    raise AssertionError(
                        "Can not read scenery from path {:s}".format(path)
                    )

        elif sceneryDs and not path:
            self.init_from_sceneryDs(sceneryDs)
        else:
            raise ValueError("Either 'path' or 'sceneryDs', but not both.")

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

    def init_from_tar(self, path):
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

    def init_from_dump(self, path):
        """
        Inits the server from a previous dump.

        Warning
        -------
        A dump is not meant to exchange sceneries with others!
        Only read dumps written on the same platform and by the same version
        of merlict. See also dump().

        Parameters
        ----------
        path : str
            Path to read the dump from.
        """
        cdef int rc
        _path = str(path)
        cdef bytes _py_path = _path.encode()
        cdef char* _cpath = _py_path
        rc = mliScenery_malloc_from_path(&self.scenery, _cpath)
        assert rc != 0

    def dump(self, path):
        """
        Dumps the compiled scenery to path.

        Warning
        -------
        A dump is not meant to exchange sceneries with others!
        A dump is specific to your platform's architecture and version of
        merlict and will probably crash on other machines.
        Think of it as a `pickle` thingy.
        The dump is meant to export and import a compiled scenery so that
        merlict does not need to compile the user's scenery again.
        So the objective of the dump is for local caching only!

        Parameters
        ----------
        path : str
            Path to write the dump to.
        """
        cdef int rc
        _path = str(path)
        cdef bytes _py_path = _path.encode()
        cdef char* _cpath = _py_path
        rc = mliScenery_write_to_path(&self.scenery, _cpath)
        assert rc != 0

    def init_from_sceneryDs(self, sceneryDs):
        cdef int rc
        cdef mliArchive tmp_archive = mliArchive_init()
        try:
            rc = mliArchive_malloc(&tmp_archive)
            assert rc != 0

            for item in sceneryDs:
                filename, payload = item
                _mliArchive_push_back_path_and_payload(
                    &tmp_archive,
                    filename,
                    payload)

            rc = mliScenery_malloc_from_Archive(&self.scenery, &tmp_archive)
            assert rc != 0
        finally:
            mliArchive_free(&tmp_archive)

    def query_intersection(self, rays):
        cdef int rc
        assert _ray.israys(rays)
        cdef stdint.uint64_t num_ray = rays.shape[0]
        isecs = _intersection.init(size=num_ray)

        cdef cnumpy.ndarray[mliRay, mode="c"] crays = np.ascontiguousarray(
            rays
        )

        cdef cnumpy.ndarray[
            mliIntersection, mode="c"
        ] cisecs = np.ascontiguousarray(
            isecs
        )

        cdef cnumpy.ndarray[
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

        cdef cnumpy.ndarray[mliRay, mode="c"] crays = np.ascontiguousarray(
            rays
        )

        cdef cnumpy.ndarray[
            mliIntersectionSurfaceNormal, mode="c"
        ] cisecs = np.ascontiguousarray(
            isecs
        )

        cdef cnumpy.ndarray[
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
