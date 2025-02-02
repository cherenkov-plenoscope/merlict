#include "bridge.h"

int mli_Archive_push_back_cstr(
        struct mli_Archive *arc,
        const char *filename,
        const uint64_t filename_length,
        const char *payload,
        const uint64_t payload_length)
{
        struct mli_String _filename = mli_String_init();
        struct mli_String _payload = mli_String_init();
        chk_msg(mli_String_malloc(&_filename, filename_length),
                "Can not malloc filename.");
        strncpy(_filename.array, filename, filename_length);
        chk_msg(mli_String_malloc(&_payload, payload_length),
                "Can not malloc payload.");
        strncpy(_payload.array, payload, payload_length);
        chk_msg(mli_Archive_push_back(arc, &_filename, &_payload),
                "Can not push back filename and payload.");
        mli_String_free(&_filename);
        mli_String_free(&_payload);
        return 1;
chk_error:
        mli_String_free(&_filename);
        mli_String_free(&_payload);
        return 0;
}

int mli_Bridge_query_many_intersection(
        const struct mli_Scenery *scenery,
        const uint64_t num_rays,
        const struct mli_Ray *rays,
        struct mli_Intersection *isecs,
        int64_t *is_valid_isecs)
{
        uint64_t i;
        for (i = 0; i < num_rays; i++) {
                struct mli_Ray ray = rays[i];
                struct mli_Intersection isec = mli_Intersection_init();
                int is_valid_isec = mli_raytracing_query_intersection(
                        scenery,
                        ray,
                        &isec
                );
                if (is_valid_isec == 1) {
                        isecs[i] = isec;
                        is_valid_isecs[i] = 1u;
                } else {
                        is_valid_isecs[i] = 0u;
                        isecs[i] = mli_Intersection_init();
                }
        }
        return 1;
}

int mli_Bridge_query_many_intersectionSurfaceNormal(
        const struct mli_Scenery *scenery,
        const uint64_t num_rays,
        const struct mli_Ray *rays,
        struct mli_IntersectionSurfaceNormal *isecs,
        int64_t *is_valid_isecs)
{
        uint64_t i;
        for (i = 0; i < num_rays; i++) {
                struct mli_Ray ray = rays[i];
                struct mli_IntersectionSurfaceNormal isec = mli_IntersectionSurfaceNormal_init();
                int is_valid_isec = mli_raytracing_query_intersection_with_surface_normal(
                        scenery,
                        ray,
                        &isec
                );
                if (is_valid_isec == 1) {
                        isecs[i] = isec;
                        is_valid_isecs[i] = 1u;
                } else {
                        is_valid_isecs[i] = 0u;
                        isecs[i] = mli_IntersectionSurfaceNormal_init();
                }
        }
        return 1;
}