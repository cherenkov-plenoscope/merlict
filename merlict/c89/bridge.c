#include "bridge.h"
#include "merlict_c89/merlict/chk_debug.h"

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

int mli_bridge_propagate_photons(
        const struct mliScenery *scenery,
        struct mliPrng *prng,
        uint64_t num_photons,
        struct mliPhoton *photons)
{
        return 1;
}
