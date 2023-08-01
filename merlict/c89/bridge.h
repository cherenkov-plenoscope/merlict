#include "merlict_c89/merlict/mli_random_generator.h"
#include "merlict_c89/merlict/mliPhoton.h"
#include "merlict_c89/merlict/mliScenery.h"
#include "merlict_c89/merlict/mliArchive.h"
#include "merlict_c89/merlict/mliIntersection.h"
#include "merlict_c89/merlict/mliIntersectionSurfaceNormal.h"


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
