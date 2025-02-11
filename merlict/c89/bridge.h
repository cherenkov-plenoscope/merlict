#include "mli.h"

chk_rc mli_Archive_push_back_cstr(
        struct mli_Archive *arc,
        const char *filename,
        const uint64_t filename_length,
        const char *payload,
        const uint64_t payload_length);


void mli_Bridge_query_many_intersection(
        const struct mli_Scenery *scenery,
        const uint64_t num_rays,
        const struct mli_Ray *rays,
        struct mli_Intersection *isecs,
        int64_t *is_valid_isecs);


void mli_Bridge_query_many_intersectionSurfaceNormal(
        const struct mli_Scenery *scenery,
        const uint64_t num_rays,
        const struct mli_Ray *rays,
        struct mli_IntersectionSurfaceNormal *isecs,
        int64_t *is_valid_isecs);
