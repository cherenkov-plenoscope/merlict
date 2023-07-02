#include "merlict_c89/merlict/mli_random_generator.h"
#include "merlict_c89/merlict/mliPhoton.h"
#include "merlict_c89/merlict/mliScenery.h"
#include "merlict_c89/merlict/mliArchive.h"

int mli_bridge_propagate_photons(
        const struct mliScenery *scenery,
        struct mliPrng *prng,
        uint64_t num_photons,
        struct mliPhoton *photons);

int mliArchive_push_back_cstr(
        struct mliArchive *arc,
        const char *filename,
        const uint64_t filename_length,
        const char *payload,
        const uint64_t payload_length);