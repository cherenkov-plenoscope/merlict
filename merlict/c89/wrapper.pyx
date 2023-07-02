import termios
import sys

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
    cdef int mliArchive_malloc(mliArchive *arc)
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

    # =======
    # PHOTONS
    # =======

    # Ray
    # ---
    cdef struct mliRay:
        mliVec support
        mliVec direction

    # Photon
    # ------
    cdef struct mliPhoton:
        mliRay ray
        double wavelength
        long id

    # GeometryId
    # ----------
    cdef struct mliGeometryId:
        unsigned int robj
        unsigned int face

    # PhotonInteraction
    # -----------------
    cdef struct mliPhotonInteraction:
        int on_geometry_surface
        mliGeometryId geometry_id
        mliVec position
        mliVec position_local
        double distance_of_ray
        unsigned long medium_coming_from
        unsigned long medium_going_to
        int from_outside_to_inside
        int type

    # ======
    # VIEWER
    # ======

    cdef struct mliView:
        mliVec position
        mliVec rotation
        double field_of_view

    cdef struct mlivrConfig:
        unsigned int random_seed
        unsigned long preview_num_cols
        unsigned long preview_num_rows
        unsigned long export_num_cols
        unsigned long export_num_rows
        double step_length
        mliView view
        double aperture_camera_f_stop_ratio
        double aperture_camera_image_sensor_width

    cdef mlivrConfig mlivrConfig_default()

    cdef int mlivr_run_interactive_viewer(const mliScenery *scn, const mlivrConfig cfg)


cdef extern from "bridge.h":
    cdef int mliArchive_push_back_cstr(
        mliArchive *arc,
        const char *filename,
        const unsigned long filename_length,
        const char *payload,
        const unsigned long payload_length,
    )

    cdef int mli_bridge_propagate_photons(
        const mliScenery *scenery,
        mliPrng *prng,
        unsigned long num_photons,
        mliPhoton *photons,
    )

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

cdef class Archihe:
    cdef mliArchive archive

    def __cinit__(self):
        self.archive = mliArchive_init()

    def __dealloc__(self):
        mliArchive_free(&self.archive)

    def __init__(self, path=None):
        mliArchive_malloc(&self.archive)
        if path:
            self._from_path(path)

    def _from_path(self, path):
        cdef bytes py_bytes = path.encode()
        cdef char* cpath = py_bytes
        cdef int rc = mliArchive_malloc_from_path(&self.archive, cpath)
        assert rc != 0

    def push_back(self, filename, payload):
        filename = str(filename)
        payload = str(payload)

        cdef unsigned long cfilename_length = len(filename)
        cdef bytes py_filename = filename.encode()
        cdef char* cfilename = py_filename

        cdef unsigned long cpayload_length = len(payload)
        cdef bytes py_payload = payload.encode()
        cdef char* cpayload = py_payload

        cdef int rc = mliArchive_push_back_cstr(
            &self.archive,
            cfilename,
            cfilename_length,
            cpayload,
            cpayload_length)
        assert rc != 0

"""
cdef class Scenery:
    cdef mliScenery scenery
    cdef mlivrConfig viewer_config
    cdef mliArchive *archive

    def __cinit__(self, mliArchive *archive):
        self.scenery = mliScenery_init()
        self.viewer_config = mlivrConfig_default()

        rc = mliScenery_malloc_from_Archive(&self.scenery, archive)
        assert rc != 0

    def __dealloc__(self):
        mliScenery_free(&self.scenery)

    def __init__(self):
        pass

    def view(self, config=None):
        self.viewer_config = _mlivrConfig_from_dict(config)

        fd = sys.stdin.fileno()
        old = termios.tcgetattr(fd)
        new = termios.tcgetattr(fd)
        C_FLAG = 3

        new[C_FLAG] = new[C_FLAG] & ~termios.ICANON
        try:
            termios.tcsetattr(fd, termios.TCSADRAIN, new)

            # ----
            rc = mlivr_run_interactive_viewer(&self.scenery, self.viewer_config)
            # ----

        finally:
            termios.tcsetattr(fd, termios.TCSADRAIN, old)
"""