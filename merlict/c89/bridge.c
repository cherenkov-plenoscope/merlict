#include "bridge.h"

int mliArchive_push_back_cstr(
        struct mliArchive *arc,
        const char *filename,
        const uint64_t filename_length,
        const char *payload,
        const uint64_t payload_length)
{
        struct mliStr str_filename = mliStr_init();
        struct mliStr str_payload = mliStr_init();
        chk_msg(mliStr_malloc(&str_filename, filename_length),
                "Can not malloc filename.");
        strncpy(str_filename.cstr, filename, filename_length);
        chk_msg(mliStr_malloc(&str_payload, payload_length),
                "Can not malloc payload.");
        strncpy(str_payload.cstr, payload, payload_length);
        chk_msg(mliArchive_push_back(arc, &str_filename, &str_payload),
                "Can not push back filename and payload.");
        mliStr_free(&str_filename);
        mliStr_free(&str_payload);
        return 1;
error:
        mliStr_free(&str_filename);
        mliStr_free(&str_payload);
        return 0;
}

int mliBridge_query_many_intersection(
        const struct mliScenery *scenery,
        const uint64_t num_rays,
        const struct mliRay *rays,
        struct mliIntersection *isecs,
        int64_t *is_valid_isecs)
{
        uint64_t i;
        for (i = 0; i < num_rays; i++) {
                struct mliRay ray = rays[i];
                struct mliIntersection isec = mliIntersection_init();
                int is_valid_isec = mli_query_intersection(
                        scenery,
                        ray,
                        &isec
                );
                if (is_valid_isec == 1) {
                        isecs[i] = isec;
                        is_valid_isecs[i] = 1u;
                } else {
                        is_valid_isecs[i] = 0u;
                        isecs[i] = mliIntersection_init();
                }
        }
        return 1;
}


int mliBridge_query_many_intersectionSurfaceNormal(
        const struct mliScenery *scenery,
        const uint64_t num_rays,
        const struct mliRay *rays,
        struct mliIntersectionSurfaceNormal *isecs,
        int64_t *is_valid_isecs)
{
        uint64_t i;
        for (i = 0; i < num_rays; i++) {
                struct mliRay ray = rays[i];
                struct mliIntersectionSurfaceNormal isec = mliIntersectionSurfaceNormal_init();
                int is_valid_isec = mli_query_intersection_with_surface_normal(
                        scenery,
                        ray,
                        &isec
                );
                if (is_valid_isec == 1) {
                        isecs[i] = isec;
                        is_valid_isecs[i] = 1u;
                } else {
                        is_valid_isecs[i] = 0u;
                        isecs[i] = mliIntersectionSurfaceNormal_init();
                }
        }
        return 1;
}