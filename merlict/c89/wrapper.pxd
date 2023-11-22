from libc cimport stdint


cdef extern from "merlict_c89.h":
    cdef struct mliVec:
        double x
        double y
        double z

    cdef struct mliRay:
        mliVec support
        mliVec direction

    cdef struct mliGeometryId:
        stdint.uint32_t robj
        stdint.uint32_t face

    cdef struct mliIntersection:
        mliGeometryId geometry_id
        mliVec position_local
        double distance_of_ray

    cdef struct mliIntersectionSurfaceNormal:
        mliGeometryId geometry_id
        mliVec position
        mliVec surface_normal
        mliVec position_local
        mliVec surface_normal_local
        double distance_of_ray
        stdint.int64_t from_outside_to_inside

    cdef struct mliPhoton:
        mliRay ray
        double wavelength
        stdint.int64_t id

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
        const mliArchive *arc)

    cdef int mliScenery_malloc_from_path(
        mliScenery *scenery,
        const char *path)
    cdef int mliScenery_write_to_path(
        const mliScenery *scenery,
        const char *path)

    cdef struct mliColor:
        float r
        float g
        float b

    cdef struct mliImage:
        stdint.uint32_t num_cols
        stdint.uint32_t num_rows
        mliColor *raw

    cdef mliColor mliImage_at(
        const mliImage *img,
        const stdint.uint32_t col,
        const stdint.uint32_t row)

    cdef void mliImage_set(
        const mliImage *img,
        const stdint.uint32_t col,
        const stdint.uint32_t row,
        const mliColor color)

    cdef struct mliPrng:
        pass

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
        const mlivrConfig cfg)

    cdef struct mliAtmosphere:
        double sunLatitude
        double sunHourAngle
        mliVec sunDirection
        double sunDistance
        double sunRadius
        double earthRadius
        double atmosphereRadius
        double Height_Rayleigh
        double Height_Mie
        mliColor beta_Rayleigh
        mliColor beta_Mie
        stdint.uint64_t numSamples
        stdint.uint64_t numSamplesLight
        double power
        double altitude

    cdef struct mliTracerConfig:
        stdint.uint64_t num_trails_global_light_source
        int have_atmosphere
        mliAtmosphere atmosphere
        mliColor background_color

    cdef int mliArchive_push_back_cstr(
        mliArchive *arc,
        const char *filename,
        const stdint.uint64_t filename_length,
        const char *payload,
        const stdint.uint64_t payload_length)

    cdef int mliBridge_query_many_intersection(
        const mliScenery *scenery,
        const stdint.uint64_t num_rays,
        const mliRay *rays,
        mliIntersection *isecs,
        stdint.int64_t *is_valid_isecs)

    cdef int mliBridge_query_many_intersectionSurfaceNormal(
        const mliScenery *scenery,
        const stdint.uint64_t num_rays,
        const mliRay *rays,
        mliIntersectionSurfaceNormal *isecs,
        stdint.int64_t *is_valid_isecs)
