#include "merlict_c89/merlict/chk_debug.h"
#include "merlict_c89/merlict/mli_random_generator.h"
#include "merlict_c89/merlict/mliPhoton.h"
#include "merlict_c89/merlict/mliScenery.h"


int mli_bridge_propagate_photons(
        const struct mliScenery *scenery,
        struct mliPrng *prng,
        uint64_t num_photons,
        struct mliPhoton *photons);
