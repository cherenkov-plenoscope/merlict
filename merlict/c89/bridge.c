#include "bridge.h"


int mli_String__from_cstr_and_length(
        struct mli_String *self,
        const char *cstr,
        const uint64_t length)
{
        chk_msg(mli_String_malloc(self, length),
                "Can not malloc mli_String to length");
        strncpy(self->array, cstr, length);
        self->size = length;

        return CHK_SUCCESS;
chk_error:
        mli_String_free(self);
        return CHK_FAIL;
}


int mli_Archive_push_back_cstr(
        struct mli_Archive *arc,
        const char *filename,
        const uint64_t filename_length,
        const char *payload,
        const uint64_t payload_length)
{
        struct mli_String _filename = mli_String_init();
        struct mli_String _payload = mli_String_init();

        chk_msg(mli_String__from_cstr_and_length(
                &_filename, filename, filename_length),
                "Can not malloc filename.");
        chk_msg(mli_String__from_cstr_and_length(
                &_payload, payload, payload_length),
                "Can not malloc filename.");

        chk_msg(mli_Archive_push_back(arc, &_filename, &_payload),
                "Can not push back filename and payload.");

        mli_String_free(&_filename);
        mli_String_free(&_payload);
        return CHK_SUCCESS;
chk_error:
        mli_String_free(&_filename);
        mli_String_free(&_payload);
        return CHK_FAIL;
}

void mli_Bridge_query_many_intersection(
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
                if (is_valid_isec == MLI_TRUE) {
                        is_valid_isecs[i] = MLI_TRUE;
                        isecs[i] = isec;
                } else {
                        is_valid_isecs[i] = MLI_FALSE;
                        isecs[i] = mli_Intersection_init();
                }
        }
        return;
}

void mli_Bridge_query_many_intersectionSurfaceNormal(
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
                if (is_valid_isec == MLI_TRUE) {
                        is_valid_isecs[i] = MLI_TRUE;
                        isecs[i] = isec;
                } else {
                        is_valid_isecs[i] = MLI_FALSE;
                        isecs[i] = mli_IntersectionSurfaceNormal_init();
                }
        }
        return;
}