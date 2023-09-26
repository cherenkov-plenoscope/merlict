#include "chk_debug_mli_core_mli_viewer.h"

int mliArchive_push_back_cstr(
        struct mliArchive *arc,
        const char *filename,
        const uint64_t filename_length,
        const char *payload,
        const uint64_t payload_length);


int mliBridge_query_many_intersection(
        const struct mliScenery *scenery,
        const uint64_t num_rays,
        const struct mliRay *rays,
        struct mliIntersection *isecs,
        int64_t *is_valid_isecs);


int mliBridge_query_many_intersectionSurfaceNormal(
        const struct mliScenery *scenery,
        const uint64_t num_rays,
        const struct mliRay *rays,
        struct mliIntersectionSurfaceNormal *isecs,
        int64_t *is_valid_isecs);
