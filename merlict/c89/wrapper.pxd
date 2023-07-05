from libc cimport stdint

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
    cdef int mliScenery_malloc_from_Archive(
        mliScenery *scn,
        const mliArchive *arc
    )
    # this is the binary dump ```sceneryCa```
    cdef int mliScenery_malloc_from_path(
        mliScenery *scenery, const char *path
    )
    cdef int mliScenery_write_to_path(
        const mliScenery *scenery, const char *path
    )

    cdef struct mliPrng:
        pass

    cdef mliPrng mliPrng_init_MT19937(const stdint.uint32_t seed)
    cdef mliPrng mliPrng_init_PCG32(const stdint.uint32_t seed)
    cdef stdint.uint32_t mliPrng_generate_uint32(mliPrng *prng)
    cdef double mli_random_uniform(mliPrng *prng)

    # GeometryId
    # ----------
    cdef struct mliGeometryId:
        stdint.uint32_t robj
        stdint.uint32_t face

    # ===
    # RAY
    # ===

    # Ray
    # ---
    cdef struct mliRay:
        mliVec support
        mliVec direction

    # mliIntersection
    # ---------------
    cdef struct mliIntersection:
        mliGeometryId geometry_id
        mliVec position_local
        double distance_of_ray

    # mliIntersectionSurfaceNormal
    # ----------------------------
    cdef struct mliIntersectionSurfaceNormal:
        mliGeometryId geometry_id
        mliVec position
        mliVec surface_normal
        mliVec position_local
        mliVec surface_normal_local
        double distance_of_ray
        int from_outside_to_inside

    # =======
    # PHOTONS
    # =======

    # Photon
    # ------
    cdef struct mliPhoton:
        mliRay ray
        double wavelength
        long id

    # PhotonInteraction
    # -----------------
    cdef struct mliPhotonInteraction:
        int on_geometry_surface
        mliGeometryId geometry_id
        mliVec position
        mliVec position_local
        double distance_of_ray
        stdint.uint64_t medium_coming_from
        stdint.uint64_t medium_going_to
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
        stdint.uint32_t random_seed
        stdint.uint64_t preview_num_cols
        stdint.uint64_t preview_num_rows
        stdint.uint64_t export_num_cols
        stdint.uint64_t export_num_rows
        double step_length
        mliView view
        double aperture_camera_f_stop_ratio
        double aperture_camera_image_sensor_width

    cdef mlivrConfig mlivrConfig_default()

    cdef int mlivr_run_interactive_viewer(
        const mliScenery *scn,
        const mlivrConfig cfg
    )

    cdef struct mliColor:
        float r
        float g
        float b

    cdef struct mliImage:
        stdint.uint32_t num_cols
        stdint.uint32_t num_rows
        mliColor *raw


cdef extern from "bridge.h":
    cdef int mliArchive_push_back_cstr(
        mliArchive *arc,
        const char *filename,
        const stdint.uint64_t filename_length,
        const char *payload,
        const stdint.uint64_t payload_length,
    )

    cdef int mli_bridge_propagate_photons(
        const mliScenery *scenery,
        mliPrng *prng,
        stdint.uint64_t num_photons,
        mliPhoton *photons,
    )
