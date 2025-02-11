from libc cimport stdint


cdef extern from "mli.h":
    cpdef enum mli_bool_states:
        MLI_FALSE,
        MLI_TRUE

    cpdef enum chk_rc_states:
        CHK_FAIL,
        CHK_SUCCESS

    cdef struct mli_Vec:
        double x
        double y
        double z

    cdef struct mli_Ray:
        mli_Vec support
        mli_Vec direction

    cdef struct mli_GeometryId:
        stdint.uint32_t robj
        stdint.uint32_t face

    cdef struct mli_Intersection:
        mli_GeometryId geometry_id
        mli_Vec position_local
        double distance_of_ray

    cdef struct mli_IntersectionSurfaceNormal:
        mli_GeometryId geometry_id
        mli_Vec position
        mli_Vec surface_normal
        mli_Vec position_local
        mli_Vec surface_normal_local
        double distance_of_ray
        stdint.int64_t from_outside_to_inside

    cdef struct mli_Photon:
        mli_Ray ray
        double wavelength
        stdint.int64_t id

    cdef struct mli_Archive:
        pass

    cdef mli_Archive mli_Archive_init()
    cdef void mli_Archive_free(mli_Archive *arc)
    cdef int mli_Archive_malloc(mli_Archive *arc)
    cdef int mli_Archive__from_path_cstr(mli_Archive *arc, const char *path)

    cdef struct mli_Scenery:
        pass

    cdef mli_Scenery mli_Scenery_init()
    cdef void mli_Scenery_free(mli_Scenery *scn)

    cdef int mli_Scenery_malloc_from_Archive(
        mli_Scenery *scn,
        const mli_Archive *arc)

    cdef int mli_Scenery_malloc_from_path(
        mli_Scenery *scenery,
        const char *path)
    cdef int mli_Scenery_write_to_path(
        const mli_Scenery *scenery,
        const char *path)

    cdef struct mli_Color:
        float r
        float g
        float b

    cdef struct mli_Image:
        pass

    cdef stdint.uint64_t mli_Image_num_cols(
        const mli_Image *img)
    cdef stdint.uint64_t mli_Image_num_rows(
        const mli_Image *img)

    cdef mli_Color mli_Image_get_by_col_row(
        const mli_Image *img,
        const stdint.uint32_t col,
        const stdint.uint32_t row)

    cdef void mli_Image_set_by_col_row(
        const mli_Image *img,
        const stdint.uint32_t col,
        const stdint.uint32_t row,
        const mli_Color color)

    cdef struct mli_Prng:
        pass

    cdef struct mli_PhotonInteraction:
        int on_geometry_surface
        mli_GeometryId geometry_id
        mli_Vec position
        mli_Vec position_local
        double distance_of_ray
        stdint.uint64_t medium_coming_from
        stdint.uint64_t medium_going_to
        int from_outside_to_inside
        int type

    cdef struct mli_View:
        mli_Vec position
        mli_Vec rotation
        double field_of_view

    cdef struct mli_viewer_Config:
        stdint.uint32_t random_seed
        stdint.uint64_t preview_num_cols
        stdint.uint64_t preview_num_rows
        stdint.uint64_t export_num_cols
        stdint.uint64_t export_num_rows
        double step_length
        mli_View view
        double gain
        double gamma

        double aperture_camera_f_stop_ratio
        double aperture_camera_image_sensor_width

    cdef mli_viewer_Config mli_viewer_Config_default()

    cdef int mli_viewer_run_interactive_viewer(
        const mli_Scenery *scn,
        const mli_viewer_Config cfg)

    cdef struct mli_ColorSpectrum:
        pass

    cdef struct mliAtmosphere:
        double sunLatitude
        double sunHourAngle
        mli_Vec sunDirection
        double sunDistance
        double sunRadius
        double earthRadius
        double atmosphereRadius
        double Height_Rayleigh
        double Height_Mie
        mli_ColorSpectrum beta_Rayleigh_spectrum
        mli_ColorSpectrum beta_Mie_spectrum
        stdint.uint64_t numSamples
        stdint.uint64_t numSamplesLight
        double power
        double altitude

    cdef int mli_Archive_push_back_cstr(
        mli_Archive *arc,
        const char *filename,
        const stdint.uint64_t filename_length,
        const char *payload,
        const stdint.uint64_t payload_length)

    cdef void mli_Bridge_query_many_intersection(
        const mli_Scenery *scenery,
        const stdint.uint64_t num_rays,
        const mli_Ray *rays,
        mli_Intersection *isecs,
        stdint.int64_t *is_valid_isecs)

    cdef void mli_Bridge_query_many_intersectionSurfaceNormal(
        const mli_Scenery *scenery,
        const stdint.uint64_t num_rays,
        const mli_Ray *rays,
        mli_IntersectionSurfaceNormal *isecs,
        stdint.int64_t *is_valid_isecs)

    # IO
    # --
    cpdef enum mli_io_type_code:
        MLI_IO_TYPE_VOID,
        MLI_IO_TYPE_FILE,
        MLI_IO_TYPE_MEMORY


    cdef struct mli_IoFile:
        pass

    cdef struct mli_IoMemory:
        unsigned char *cstr
        stdint.uint64_t capacity
        stdint.uint64_t size
        stdint.uint64_t pos

    cdef union mli_IoType:
        mli_IoFile file
        mli_IoMemory memory

    cdef struct mli_IO:
        int type
        mli_IoType data

    cdef mli_IO mli_IO_init()
    cdef int mli_IO_open_memory(mli_IO *self)
    cdef int mli_IO_close(mli_IO *self)

    cdef int mli_IoMemory__malloc_capacity(
        mli_IoMemory *self,
        const stdint.uint64_t capacity)

    cdef int mli_Scenery_to_io(const mli_Scenery *self, mli_IO *f)
    cdef int mli_Scenery_from_io(mli_Scenery *self, mli_IO *f)
