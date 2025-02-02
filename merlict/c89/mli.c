#include "mli.h"

#include <assert.h>
#include <float.h>
#include <limits.h>
#include <time.h>
#include <unistd.h>

/* Config */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_viewer_Config mli_viewer_Config_default(void)
{
        struct mli_viewer_Config cfg;
        cfg.random_seed = 0u;

        cfg.preview_num_cols = 160u;
        cfg.preview_num_rows = 90u / 2;

        cfg.export_num_cols = 1280u;
        cfg.export_num_rows = 720u;

        cfg.step_length = 1.0;
        cfg.gamma = 1.0;

        cfg.view.position.x = 0.;
        cfg.view.position.y = 0.;
        cfg.view.position.z = 0.;

        cfg.view.rotation.x = mli_math_deg2rad(90.);
        cfg.view.rotation.y = 0.;
        cfg.view.rotation.z = 0.;

        cfg.view.field_of_view = mli_math_deg2rad(80.);

        cfg.aperture_camera_f_stop_ratio = 0.95;
        cfg.aperture_camera_image_sensor_width = 64e-3;
        return cfg;
}

/* CorsikaPhotonBunch */
/* ------------------ */

/* Copyright 2016 Sebastian A. Mueller */

MLI_VECTOR_IMPLEMENTATION(
        mliDynCorsikaPhotonBunch,
        struct mli_corsika_PhotonBunch)

void mli_corsika_PhotonBunch_set_from_raw(
        struct mli_corsika_PhotonBunch *bunch,
        const float *raw)
{
        bunch->x_cm = raw[0];
        bunch->y_cm = raw[1];
        bunch->ux = raw[2];
        bunch->vy = raw[3];
        bunch->time_ns = raw[4];
        bunch->z_emission_cm = raw[5];
        bunch->weight_photons = raw[6];
        bunch->wavelength_nm = raw[7];
}

void mli_corsika_PhotonBunch_to_raw(
        const struct mli_corsika_PhotonBunch *bunch,
        float *raw)
{
        raw[0] = bunch->x_cm;
        raw[1] = bunch->y_cm;
        raw[2] = bunch->ux;
        raw[3] = bunch->vy;
        raw[4] = bunch->time_ns;
        raw[5] = bunch->z_emission_cm;
        raw[6] = bunch->weight_photons;
        raw[7] = bunch->wavelength_nm;
}

struct mli_Photon mli_corsika_PhotonBunch_to_merlict_photon(
        const struct mli_corsika_PhotonBunch bunch,
        const double production_distance_offset,
        const int64_t id)
{
        /*
        Returns an mli_Photon that will reach the observation-level in
        the same way as the corsika-photon-bunch. The weight of the
        corsika-photon-bunch is not taken into account here.

        Parameters
        ----------
        bunch :
                The corsika-photon-bunch
        production_distance_offset : double
                An arbitrary distance for the photon to travel until they
                reach the observation-level. If 0.0, the distance for a
                merlict photon is only defined by the relative arrival time
                on the observation-level.
                Ensure this offset is at least as big as your detector system
                so that photons do not start inside your detector.
        id : int64
                The photon's id.
        */

        const double vacuum_speed_of_light = MLI_PHYSICS_SPEED_OF_LIGHT_M_PER_S;
        const struct mli_Vec photon_direction_of_motion =
                mli_corsika_photon_direction_of_motion(bunch);

        const struct mli_Ray ray_running_upwards_to_production = mli_Ray_set(
                mli_corsika_photon_support_on_observation_level(bunch),
                mli_Vec_multiply(photon_direction_of_motion, -1.0));

        const double offset =
                (production_distance_offset +
                 vacuum_speed_of_light *
                         mli_corsika_photon_relative_arrival_time_on_observation_level(
                                 bunch));

        const struct mli_Vec photon_emission_position =
                mli_Ray_at(&ray_running_upwards_to_production, offset);

        struct mli_Photon photon;
        photon.ray.support = photon_emission_position;
        photon.ray.direction = photon_direction_of_motion;
        photon.wavelength = mli_corsika_photon_wavelength(bunch);
        photon.id = id;
        return photon;
}

struct mli_Vec mli_corsika_photon_direction_of_motion(
        const struct mli_corsika_PhotonBunch bunch)
{ /*
       KIT-CORSIKA coordinate-system

                         /\ z-axis
                         |
                         |\ p
                         | \ a
                         |  \ r
                         |   \ t
                         |    \ i
                         |     \ c
                         |      \ l
                         |       \ e
                         |        \
                         |  theta  \ m
                         |       ___\ o
                         |___----    \ m      ___
                         |            \ e       /| y-axis (west)
                         |             \ n    /
                         |              \ t /
                         |               \/u
                         |              / \ m
                         |            /    \
                         |          /       \
                         |        /__________\
                         |      /      ___---/
                         |    /   __---    /
                         |  /__--- phi \ /
         ________________|/--__________/______\ x-axis (north)
                        /|                    /
                      /  |
                    /    |
                  /


          Extensive Air Shower Simulation with CORSIKA, Figure 1, page 114
          (Version 7.6400 from December 27, 2017)

          Direction-cosines:

          cx = sin(theta) * cos(phi)
          cy = sin(theta) * sin(phi)

          The zenith-angle theta opens relative to the negative z-axis.

          It is the momentum of the Cherenkov-photon, which is pointing
          down towards the observation-plane.
  */
        const double z = sqrt(1.0 - bunch.ux * bunch.ux - bunch.vy * bunch.vy);
        return mli_Vec_init(bunch.ux, bunch.vy, -z);
}

struct mli_Vec mli_corsika_photon_support_on_observation_level(
        const struct mli_corsika_PhotonBunch bunch)
{
        return mli_Vec_init(
                (double)bunch.x_cm * 1e-2, (double)bunch.y_cm * 1e-2, 0.0);
}

double mli_corsika_photon_wavelength(const struct mli_corsika_PhotonBunch bunch)
{
        return fabs((double)bunch.wavelength_nm * 1e-9);
}

double mli_corsika_photon_emission_height(
        const struct mli_corsika_PhotonBunch bunch)
{
        return (double)bunch.z_emission_cm * 1e-2;
}

double mli_corsika_photon_relative_arrival_time_on_observation_level(
        const struct mli_corsika_PhotonBunch bunch)
{
        return (double)bunch.time_ns * 1e-9;
}

/* Cursor */
/* ------ */

/* Copyright 2019 Sebastian Achim Mueller                                     */

void mli_viewer_Cursor_move_up(struct mli_viewer_Cursor *cursor)
{
        if (cursor->row != 0)
                cursor->row -= 1;
}

void mli_viewer_Cursor_move_down(struct mli_viewer_Cursor *cursor)
{
        if (cursor->row + 1u < cursor->num_rows)
                cursor->row += 1;
}

void mli_viewer_Cursor_move_right(struct mli_viewer_Cursor *cursor)
{
        if (cursor->col != 0)
                cursor->col -= 1;
}

void mli_viewer_Cursor_move_left(struct mli_viewer_Cursor *cursor)
{
        if (cursor->col + 1u < cursor->num_cols)
                cursor->col += 1;
}

/* EventIo_Event */
/* ------------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */

int mliEventIo_read_telescope_offsets(
        FILE *f,
        struct mliDynEventIoTelescopeOffset *telescope_offsets,
        const uint64_t length)
{
        const int length_first_two = 4 + 4;
        int num_following_arrays;
        int32_t _narray;
        uint64_t narray;
        uint64_t n;
        float toff;

        struct mli_FloatArray xoff = mli_FloatArray_init();
        struct mli_FloatArray yoff = mli_FloatArray_init();
        struct mli_FloatArray weight = mli_FloatArray_init();

        chk_fread(&_narray, sizeof(int32_t), 1, f);
        chk_msg(_narray >= 0, "Expected num. of arrays to be positive.");
        narray = (uint64_t)_narray;

        chk_fread(&toff, sizeof(float), 1, f);

        chk(mli_FloatArray_malloc(&xoff, narray));
        chk(mli_FloatArray_malloc(&yoff, narray));
        chk(mli_FloatArray_malloc(&weight, narray));

        num_following_arrays = (int)((length - length_first_two) / narray / 4);

        chk_fread(xoff.array, sizeof(float), narray, f);
        chk_fread(yoff.array, sizeof(float), narray, f);
        memset(weight.array, 1.0, narray);

        switch (num_following_arrays) {
        case 2:
                /* do nothing, we already read the 2 arrays xoff and yoff. */
                break;
        case 3:
                chk_fread(weight.array, sizeof(float), narray, f);
                break;
        default:
                chk_bad("Expected num_following_arrays to be either 2, or 3.");
        }

        chk_msg(mliDynEventIoTelescopeOffset_malloc(telescope_offsets, narray),
                "Failed to malloc telescope_offsets.");

        for (n = 0; n < narray; n++) {
                telescope_offsets->array[n].xoff = xoff.array[n];
                telescope_offsets->array[n].yoff = yoff.array[n];
                telescope_offsets->array[n].toff = toff;
                telescope_offsets->array[n].weight = weight.array[n];
        }

        mli_FloatArray_free(&xoff);
        mli_FloatArray_free(&yoff);
        mli_FloatArray_free(&weight);

        return 1;
chk_error:
        return 0;
}

struct mliEventIoBunchHeader {
        int16_t array_id;
        int16_t telescope_id;
        float photons;
        int32_t num_bunches;
};

int mliEventIo_read_photon_bunches(
        FILE *f,
        struct mliEventIoBunchHeader *b_head,
        struct mliDynCorsikaPhotonBunch *bunches,
        const int32_t version)
{
        float tmp[8];
        int is_compact = 0;
        uint64_t row, field;

        chk_fread(&b_head->array_id, sizeof(b_head->array_id), 1, f);
        chk_fread(&b_head->telescope_id, sizeof(b_head->telescope_id), 1, f);
        chk_fread(&b_head->photons, sizeof(b_head->photons), 1, f);
        chk_fread(&b_head->num_bunches, sizeof(b_head->num_bunches), 1, f);

        is_compact = (int)(version / 1000 == 1);

        chk_msg(mliDynCorsikaPhotonBunch_malloc(bunches, b_head->num_bunches),
                "Failed to malloc bunches.");
        bunches->size = b_head->num_bunches;

        if (is_compact) {
                int16_t buff[8];
                for (row = 0; row < bunches->size; row++) {
                        chk_fread(buff, sizeof(int16_t), 8, f);

                        for (field = 0; field < 8; field++) {
                                tmp[field] = (float)buff[field];
                        }

                        bunches->array[row].x_cm = tmp[0] * 0.1;
                        bunches->array[row].y_cm = tmp[1] * 0.1;
                        bunches->array[row].ux = tmp[2] / 30000;
                        bunches->array[row].vy = tmp[3] / 30000;
                        bunches->array[row].time_ns = tmp[4] * 0.1;
                        bunches->array[row].z_emission_cm =
                                pow(10, tmp[5] * 0.001);
                        bunches->array[row].weight_photons = tmp[6] * 0.01;
                        bunches->array[row].wavelength_nm = tmp[7];
                }
        } else {
                for (row = 0; row < bunches->size; row++) {
                        chk_fread(tmp, sizeof(float), 8, f);

                        bunches->array[row].x_cm = tmp[0];
                        bunches->array[row].y_cm = tmp[1];
                        bunches->array[row].ux = tmp[2];
                        bunches->array[row].vy = tmp[3];
                        bunches->array[row].time_ns = tmp[4];
                        bunches->array[row].z_emission_cm = tmp[5];
                        bunches->array[row].weight_photons = tmp[6];
                        bunches->array[row].wavelength_nm = tmp[7];
                }
        }

        return 1;
chk_error:
        return 0;
}

struct mliEventIoEvent mliEventIoEvent_init(void)
{
        struct mliEventIoEvent event;
        memset(event.corsika_event_header, 0.0, 273 * sizeof(float));
        memset(event.corsika_event_end, 0.0, 273 * sizeof(float));
        event.telescopes = mliDynEventIoTelescope_init();
        return event;
}
void mliEventIoEvent_free(struct mliEventIoEvent *event)
{
        uint64_t i;
        for (i = 0; i < event->telescopes.size; i++) {
                mliEventIoTelescope_free(&event->telescopes.array[i]);
        }
        mliDynEventIoTelescope_free(&event->telescopes);
        (*event) = mliEventIoEvent_init();
}

int mliEventIoEvent_malloc_from_run(
        struct mliEventIoEvent *event,
        struct mliEventIoRun *run)
{
        uint64_t remaining_length = 0;
        uint64_t header_length = 0;
        uint64_t num_sub_blocks = 0;
        uint64_t i;
        struct mliDynEventIoTelescopeOffset tmp_offsets =
                mliDynEventIoTelescopeOffset_init();

        mliEventIoEvent_free(event);

        /* corsika_event_header */
        /* -------------------- */
        chk_msg(run->next_block.type == 1202, "Expected type 1202.");
        chk_msg(mliEventIoRun_read_273float32_block(
                        run->f, event->corsika_event_header),
                "Failed to read corsika_event_header 273 float block.");

        /* telescope_offsets */
        /* ----------------- */
        chk(mliEventIoRun_next_block(run, MLI_EVENTIO_TOP_LEVEL));
        chk_msg(run->next_block.type == 1203, "Expected type 1203.");
        chk_msg(mliEventIo_read_telescope_offsets(
                        run->f, &tmp_offsets, run->next_block.length),
                "Failed to read telescope_offsets.");

        /* array_header */
        /* ------------ */
        chk(mliEventIoRun_next_block(run, MLI_EVENTIO_TOP_LEVEL));
        chk_msg(run->next_block.type == 1204, "Expected type 1204.");
        chk_msg(run->next_block.only_sub_objects,
                "Expected telescope-array-block to only contain sub-blocks.");
        remaining_length = run->next_block.length;

        /* telescopes */
        /* ---------- */
        chk_msg(mliDynEventIoTelescope_malloc(
                        &event->telescopes, tmp_offsets.size),
                "Failed to malloc telescpes");
        for (i = 0; i < event->telescopes.size; i++) {
                event->telescopes.array[i] = mliEventIoTelescope_init();
                event->telescopes.array[i].offset = tmp_offsets.array[i];
        }
        mliDynEventIoTelescopeOffset_free(&tmp_offsets);

        while (remaining_length) {

                /* photon_bunches for each telescope */
                /* --------------------------------- */
                struct mliEventIoBunchHeader b_head;
                struct mliEventIoTelescope *telescope =
                        &event->telescopes.array[num_sub_blocks];

                header_length = mliEventIoHeader_read(
                        &run->next_block, run->f, MLI_EVENTIO_SUB_LEVEL);
                chk_msg(header_length,
                        "Failed to read EventIo-SUB-block-header.");
                remaining_length -= header_length;

                chk_msg(run->next_block.type == 1205,
                        "Expected subheader type 1205");

                telescope->array_id = b_head.array_id;
                telescope->telescope_id = b_head.telescope_id;
                chk_msg(mliEventIo_read_photon_bunches(
                                run->f,
                                &b_head,
                                &telescope->photon_bunches,
                                run->next_block.version),
                        "Failed to read photon_bunches.");
                remaining_length -= run->next_block.length;
                num_sub_blocks += 1;
        }

        chk_msg(remaining_length == 0, "Expected remaining_length == 0.");
        chk_msg(num_sub_blocks == event->telescopes.size,
                "Expected every telescope to have photon_bunches.");

        /* corsika_event_end */
        /* ----------------- */
        chk(mliEventIoRun_next_block(run, MLI_EVENTIO_TOP_LEVEL));
        chk_msg(run->next_block.type == 1209, "Expected type 1209.");
        chk_msg(mliEventIoRun_read_273float32_block(
                        run->f, event->corsika_event_end),
                "Failed to read corsika_event_end 273 float block.");

        /* next */
        /* ---- */
        chk(mliEventIoRun_next_block(run, MLI_EVENTIO_TOP_LEVEL));

        return 1;
chk_error:
        return 0;
}

/* EventIo_Header */
/* -------------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */

#define EXPECTED_SYNC -736130505

struct mliEventIoHeader mliEventIoHeader_init(void)
{
        struct mliEventIoHeader h;
        h.is_sync = 0;
        h.type = -1;
        h.version = -1;
        h.user = -1;
        h.extended = 0;
        h.only_sub_objects = 0;
        h.length = 0u;
        h.id = -1;
        return h;
}

struct _FirstFour {
        int32_t sync;
        int32_t type;
        int32_t id;
        int32_t length;
};
struct _FirstFour _FirstFour_zeros(void)
{
        struct _FirstFour ff;
        ff.sync = 0;
        ff.type = 0;
        ff.id = 0;
        ff.length = 0;
        return ff;
}

struct _TypeInfo {
        int32_t type;
        int32_t version;
        int user;
        int extended;
};
struct _TypeInfo _TypeInfo_zeros(void)
{
        struct _TypeInfo ti;
        ti.type = 0;
        ti.version = 0;
        ti.user = 0;
        ti.extended = 0;
        return ti;
}
struct _TypeInfo _TypeInfo_init(const int32_t _type)
{
        struct _TypeInfo info;
        info.type = _type & 0xffff;
        info.version = (_type & 0xfff00000) >> 20;
        info.user = (int)(_type & (1 << 16));
        info.extended = (int)(_type & (1 << 17));
        return info;
}

struct _LengthInfo {
        int only_sub_objects;
        int32_t length;
};
struct _LengthInfo _LengthInfo_zeros(void)
{
        struct _LengthInfo li;
        li.only_sub_objects = 0;
        li.length = 0;
        return li;
}
struct _LengthInfo _LengthInfo_init(const int32_t _length)
{
        struct _LengthInfo info;
        info.only_sub_objects = (int)(_length & 1 << 30);
        /* bit 31 of length is reserved */
        info.length = _length & 0x3fffffff;
        return info;
}

int64_t _extend_length(
        const int32_t extended,
        const struct _LengthInfo length_info)
{
        int64_t ext, len;
        ext = extended;
        ext &= 0xfff;
        len = length_info.length;
        len &= (ext << 12);
        return len;
}

int mliEventIoHeader_read(struct mliEventIoHeader *header, FILE *f, int level)
{
        int length_read = 0;
        struct _FirstFour first_four = _FirstFour_zeros();
        struct _TypeInfo type_info = _TypeInfo_zeros();
        struct _LengthInfo length_info = _LengthInfo_zeros();
        /* sub level headers do not have the 'sync' field. */

        if (level == MLI_EVENTIO_TOP_LEVEL) {
                chk_fread(&first_four.sync, sizeof(first_four.sync), 1, f);
                length_read += sizeof(first_four.sync);
                header->is_sync = EXPECTED_SYNC == first_four.sync;
        } else {
                chk_msg(level == MLI_EVENTIO_SUB_LEVEL,
                        "Level must be either 'top' or 'sub'.");
                header->is_sync = 1;
        }
        chk_fread(&first_four.type, sizeof(first_four.type), 1, f);
        length_read += sizeof(first_four.type);
        chk_fread(&first_four.id, sizeof(first_four.id), 1, f);
        length_read += sizeof(first_four.id);
        chk_fread(&first_four.length, sizeof(first_four.length), 1, f);
        length_read += sizeof(first_four.length);

        type_info = _TypeInfo_init(first_four.type);
        length_info = _LengthInfo_init(first_four.length);

        chk_msg(header->is_sync,
                "Expected EventIo-Header to be in sync, but its not.");

        header->id = first_four.id;

        header->type = type_info.type;
        header->version = type_info.version;
        header->user = type_info.user;
        header->extended = type_info.extended;
        header->only_sub_objects = length_info.only_sub_objects;

        if (!type_info.extended) {
                header->length = length_info.length;
        } else {
                int32_t extended;
                chk_fread(&extended, sizeof(int32_t), 1, f);
                length_read += sizeof(int32_t);
                header->length = _extend_length(extended, length_info);
        }

        return length_read;
chk_error:
        return 0;
}

void mliEventIoHeader_fprint(const struct mliEventIoHeader head, FILE *f)
{
        fprintf(f, "----\n");
        fprintf(f, "h.is_sync %d\n", head.is_sync);
        fprintf(f, "h.type %d\n", head.type);
        fprintf(f, "h.version %d\n", head.version);
        fprintf(f, "h.user %d\n", head.user);
        fprintf(f, "h.extended %d\n", head.extended);
        fprintf(f, "h.only_sub_objects %d\n", head.only_sub_objects);
        fprintf(f, "h.length %lu\n", head.length);
        fprintf(f, "h.id %d\n", head.id);
}

/* EventIo_Run */
/* ----------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */

struct mliEventIoRun mliEventIoRun_init(void)
{
        struct mliEventIoRun run;
        run.f = NULL;
        run.next_block = mliEventIoHeader_init();
        memset(run.corsika_run_header, 0.0, 273 * sizeof(float));
        run.corsika_input_card = mli_String_init();
        run.telescope_positions = mliDynEventIoTelescopePosition_init();
        return run;
}

int mliEventIoRun_read_273float32_block(FILE *f, float *block)
{
        int32_t block_size;
        chk_fread(&block_size, sizeof(int32_t), 1, f);
        chk_msg(block_size == 273, "Expected block-size to be 273.");
        chk_fread(block, sizeof(float), (uint64_t)block_size, f);
        return 1;
chk_error:
        return 0;
}

/*
 *      runh 1200
 *      rune 1210
 *      evth 1202
 *      evte 1209
 *      array_header 1204
 *      bunches 1205 (sub)
 */

int mliEventIoRun_read_input_card(
        FILE *f,
        struct mli_String *input_card,
        const uint64_t length)
{
        char _unknown[8];
        uint64_t input_card_length;

        chk_msg(mli_String_malloc(input_card, length + 1),
                "Failed to malloc cstr for input-card.");

        chk_fread(_unknown, sizeof(_unknown), 1, f);
        chk_msg(length >= sizeof(_unknown),
                "Expected at least 8bytes payload.");
        input_card_length = length - sizeof(_unknown);
        chk_fread(input_card->array, sizeof(char), input_card_length, f);
        return 1;
chk_error:
        return 0;
}

int mliEventIoRun_read_telescope_positions(
        FILE *f,
        struct mliDynEventIoTelescopePosition *telescope_positions,
        const uint64_t length)
{
        int32_t ntel;
        int num_following_arrays;

        chk_fread(&ntel, sizeof(int32_t), 1, f);

        num_following_arrays = (int)((length - 4) / ntel / 4);

        chk_msg(num_following_arrays == 4, "Expected exactly four arrays.")

                chk(mliDynEventIoTelescopePosition_malloc(
                        telescope_positions, ntel));
        chk_fread(
                telescope_positions->array,
                sizeof(struct mliEventIoTelescopePosition),
                (uint64_t)ntel,
                f);
        return 1;
chk_error:
        return 0;
}

int mliEventIoRun_next_block(struct mliEventIoRun *run, const int level)
{
        chk_msg(mliEventIoHeader_read(&run->next_block, run->f, level),
                "Failed to read EventIo-block-header.");
        return 1;
chk_error:
        return 0;
}

int mliEventIoRun_begin(struct mliEventIoRun *run, FILE *stream)
{
        run->f = stream;
        chk_msg(run->f, "Expected EventIo.f to be open.");

        /* corsika_run_header */
        /* ------------------ */
        chk(mliEventIoRun_next_block(run, MLI_EVENTIO_TOP_LEVEL));
        chk_msg(run->next_block.type == 1200, "Expected type 1200.");
        chk_msg(mliEventIoRun_read_273float32_block(
                        run->f, run->corsika_run_header),
                "Failed to read corsika_run_header 273 float block.");

        /* corsika_input_card */
        /* ------------------ */
        chk(mliEventIoRun_next_block(run, MLI_EVENTIO_TOP_LEVEL));
        chk_msg(run->next_block.type == 1212, "Expected type 1212.");
        chk_msg(mliEventIoRun_read_input_card(
                        run->f,
                        &run->corsika_input_card,
                        run->next_block.length),
                "Failed to read corsika-input-card.");

        /* telescope_positions */
        /* ------------------- */
        chk(mliEventIoRun_next_block(run, MLI_EVENTIO_TOP_LEVEL));
        chk_msg(run->next_block.type == 1201, "Expected type 1201.");
        chk_msg(mliEventIoRun_read_telescope_positions(
                        run->f,
                        &run->telescope_positions,
                        run->next_block.length),
                "Failed to read telescope-positions.");

        /* next */
        /* ---- */
        chk(mliEventIoRun_next_block(run, MLI_EVENTIO_TOP_LEVEL));

        return 1;
chk_error:
        return 0;
}

int mliEventIoRun_has_still_events_left(struct mliEventIoRun *run)
{
        const int32_t CORSIKA_RUN_END_TYPE = 1210;
        return run->next_block.type != CORSIKA_RUN_END_TYPE;
}

void mliEventIoRun_finalize(struct mliEventIoRun *run)
{
        mli_String_free(&run->corsika_input_card);
        mliDynEventIoTelescopePosition_free(&run->telescope_positions);
        (*run) = mliEventIoRun_init();
}

/* EventIo_Telescope */
/* ----------------- */

/* Copyright 2020 Sebastian A. Mueller*/

struct mliEventIoTelescope mliEventIoTelescope_init(void)
{
        struct mliEventIoTelescope telescope;
        telescope.array_id = -1;
        telescope.telescope_id = -1;
        telescope.offset = mliEventIoTelescopeOffset_init();
        telescope.photon_bunches = mliDynCorsikaPhotonBunch_init();
        return telescope;
}

void mliEventIoTelescope_free(struct mliEventIoTelescope *telescope)
{
        mliDynCorsikaPhotonBunch_free(&telescope->photon_bunches);
        (*telescope) = mliEventIoTelescope_init();
}

MLI_ARRAY_IMPLEMENTATION(mliDynEventIoTelescope, struct mliEventIoTelescope)

/* EventIo_TelescopeOffset */
/* ----------------------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */

struct mliEventIoTelescopeOffset mliEventIoTelescopeOffset_init(void)
{
        struct mliEventIoTelescopeOffset offset;
        offset.toff = 0.0;
        offset.xoff = 0.0;
        offset.yoff = 0.0;
        offset.weight = 0.0;
        return offset;
}

MLI_ARRAY_IMPLEMENTATION(
        mliDynEventIoTelescopeOffset,
        struct mliEventIoTelescopeOffset)

/* EventIo_TelescopePosition */
/* ------------------------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */

MLI_ARRAY_IMPLEMENTATION(
        mliDynEventIoTelescopePosition,
        struct mliEventIoTelescopePosition)

/* EventTape */
/* --------- */

/* Copyright 2020 Sebastian A. Mueller */

/* writer */
/* ====== */
struct mliEventTapeWriter mliEventTapeWriter_init(void)
{
        struct mliEventTapeWriter tio;
        tio.tar = mli_Tar_init();
        tio.flush_tar_stream_after_each_file = 1;
        tio.run_number = 0;
        tio.event_number = 0;
        tio.cherenkov_bunch_block_number = 1;
        tio.buffer = mli_FloatVector_init();
        return tio;
}

int mliEventTapeWriter_finalize(struct mliEventTapeWriter *tio)
{
        if (tio->tar.stream) {
                chk_msg(mliEventTapeWriter_flush_cherenkov_bunch_block(tio),
                        "Can't finalize cherenkov-bunch-block.");
                chk_msg(mli_Tar_write_finalize(&tio->tar),
                        "Can't finalize tar-file.");
        }
        mli_FloatVector_free(&tio->buffer);
        (*tio) = mliEventTapeWriter_init();
        return 1;
chk_error:
        return 0;
}

int mliEventTapeWriter_begin(
        struct mliEventTapeWriter *tio,
        struct mli_IO *stream,
        const uint64_t num_bunches_buffer)
{
        chk_msg(mliEventTapeWriter_finalize(tio),
                "Can't close and free previous tar-io-writer.");
        chk_msg(mli_Tar_write_begin(&tio->tar, stream), "Can't begin tar.");
        chk_msg(mli_FloatVector_malloc(&tio->buffer, 8 * num_bunches_buffer),
                "Can't malloc cherenkov-bunch-buffer.");
        return 1;
chk_error:
        return 0;
}

int mliEventTapeWriter_write_corsika_header(
        struct mliEventTapeWriter *tio,
        const char *path,
        const float *corsika_header)
{
        struct mli_TarHeader tarh = mli_TarHeader_init();
        chk_msg(mli_TarHeader_set_normal_file(
                        &tarh, path, MLI_CORSIKA_HEADER_SIZE_BYTES),
                "Can't set tar-header for corsika-header.");
        chk_msg(mli_Tar_write_header(&tio->tar, &tarh),
                "Can't write tar-header for corsika-header to tar.");
        chk_msg(mli_Tar_write_data(
                        &tio->tar,
                        corsika_header,
                        MLI_CORSIKA_HEADER_SIZE_BYTES),
                "Can't write data of corsika-header to tar.");
        if (tio->flush_tar_stream_after_each_file) {
                mli_IO_flush(tio->tar.stream);
        }
        return 1;
chk_error:
        return 0;
}

int mliEventTapeWriter_write_runh(
        struct mliEventTapeWriter *tio,
        const float *runh)
{
        char path[MLI_TAR_NAME_LENGTH] = {'\0'};
        tio->run_number = (int)(MLI_MATH_ROUND(runh[1]));
        chk_msg(tio->run_number >= 1, "Expected run_number >= 1.");
        sprintf(path, "%09d/RUNH.float32", tio->run_number);
        chk_msg(mliEventTapeWriter_write_corsika_header(tio, path, runh),
                "Can't write 'RUNH.float32' to event-tape.");
        return 1;
chk_error:
        return 0;
}

int mliEventTapeWriter_write_evth(
        struct mliEventTapeWriter *tio,
        const float *evth)
{
        char path[MLI_TAR_NAME_LENGTH] = {'\0'};
        int evth_run_number =
                (int)(MLI_MATH_ROUND(evth[MLI_CORSIKA_EVTH_RUN_NUMBER]));

        if (tio->event_number > 0) {
                chk_msg(mliEventTapeWriter_flush_cherenkov_bunch_block(tio),
                        "Can't finalize cherenkov-bunch-block.");
        }
        chk_msg(tio->run_number != 0, "Expected RUNH before EVTH.");
        chk_msg(tio->run_number == evth_run_number,
                "Expected run_number in EVTH "
                "to match run_number in last RUNH.");
        tio->event_number =
                (int)(MLI_MATH_ROUND(evth[MLI_CORSIKA_EVTH_EVENT_NUMBER]));
        chk_msg(tio->event_number > 0, "Expected event_number > 0.");
        tio->cherenkov_bunch_block_number = 1;
        sprintf(path,
                "%09d/%09d/EVTH.float32",
                tio->run_number,
                tio->event_number);
        chk_msg(mliEventTapeWriter_write_corsika_header(tio, path, evth),
                "Can't write 'EVTH.float32' to event-tape.");
        return 1;
chk_error:
        return 0;
}

int mliEventTapeWriter_flush_cherenkov_bunch_block(
        struct mliEventTapeWriter *tio)
{
        char path[MLI_TAR_NAME_LENGTH] = {'\0'};
        struct mli_TarHeader tarh = mli_TarHeader_init();
        sprintf(path,
                "%09d/%09d/%09d.cer.x8.float32",
                tio->run_number,
                tio->event_number,
                tio->cherenkov_bunch_block_number);
        chk_msg(mli_TarHeader_set_normal_file(
                        &tarh, path, tio->buffer.size * sizeof(float)),
                "Can't set cherenkov-bunch-block's tar-header.");
        chk_msg(mli_Tar_write_header(&tio->tar, &tarh),
                "Can't write tar-header for cherenkov-bunch-block to tar.");
        chk_msg(mli_Tar_write_data(
                        &tio->tar,
                        tio->buffer.array,
                        tio->buffer.size * sizeof(float)),
                "Can't write cherenkov-bunch-block to tar-file.");
        if (tio->flush_tar_stream_after_each_file) {
                mli_IO_flush(tio->tar.stream);
        }
        tio->buffer.size = 0;
        tio->cherenkov_bunch_block_number += 1;
        return 1;
chk_error:
        return 0;
}

int mliEventTapeWriter_write_cherenkov_bunch(
        struct mliEventTapeWriter *tio,
        const float *bunch)
{
        uint64_t i;
        if (tio->buffer.size == tio->buffer.capacity) {
                chk_msg(mliEventTapeWriter_flush_cherenkov_bunch_block(tio),
                        "Can't finalize cherenkov-bunch-block.");
        }
        for (i = 0; i < 8; i++) {
                tio->buffer.array[tio->buffer.size] = bunch[i];
                tio->buffer.size += 1;
        }
        return 1;
chk_error:
        return 0;
}

/* reader */
/* ====== */

struct mliEventTapeReader mliEventTapeReader_init(void)
{
        struct mliEventTapeReader tio;
        tio.tar = mli_Tar_init();
        tio.tarh = mli_TarHeader_init();
        tio.run_number = 0;
        tio.event_number = 0;
        tio.cherenkov_bunch_block_number = 0;
        tio.block_at = 0;
        tio.block_size = 0;
        tio.has_still_bunches_in_event = 0;
        return tio;
}

int mliEventTapeReader_finalize(struct mliEventTapeReader *tio)
{
        if (tio->tar.stream) {
                chk_msg(mli_Tar_read_finalize(&tio->tar),
                        "Can't finalize reading tar.");
        }
        (*tio) = mliEventTapeReader_init();
        return 1;
chk_error:
        return 0;
}

int mliEventTapeReader_begin(
        struct mliEventTapeReader *tio,
        struct mli_IO *stream)
{
        chk_msg(mliEventTapeReader_finalize(tio),
                "Can't close and free previous tar-io-reader.");
        chk_msg(mli_Tar_read_begin(&tio->tar, stream), "Can't begin tar.");
        tio->has_tarh = mli_Tar_read_header(&tio->tar, &tio->tarh);
        return 1;
chk_error:
        return 0;
}

int mliEventTapeReader_read_runh(struct mliEventTapeReader *tio, float *runh)
{
        const uint64_t NUM_DIGITS = 9;
        const uint64_t BASE = 10;
        uint64_t runh_run_number = 0;
        chk_msg(tio->has_tarh, "Expected next tar-header.");
        chk_msg(mli_cstr_match_templeate(
                        tio->tarh.name, "ddddddddd/RUNH.float32", 'd'),
                "Expected file to be 'ddddddddd/RUNH.float32.'");
        chk_msg(tio->tarh.size == MLI_CORSIKA_HEADER_SIZE_BYTES,
                "Expected RUNH to have size 273*sizeof(float)");
        chk_msg(mli_Tar_read_data(&tio->tar, (void *)runh, tio->tarh.size),
                "Can't read RUNH from tar.");
        chk_msg(runh[0] == mli_chars_to_float("RUNH"),
                "Expected RUNH[0] == 'RUNH'");
        chk_msg(mli_cstr_nto_uint64(
                        &tio->run_number, &tio->tarh.name[0], BASE, NUM_DIGITS),
                "Can't read run_number from RUNH's path.");
        runh_run_number =
                (uint64_t)(MLI_MATH_ROUND(runh[MLI_CORSIKA_RUNH_RUN_NUMBER]));
        chk_msg(tio->run_number == runh_run_number,
                "Expected run_number in RUNH's path "
                "to match run_number in RUNH.");
        tio->has_tarh = mli_Tar_read_header(&tio->tar, &tio->tarh);
        return 1;
chk_error:
        return 0;
}

int mliEventTapeReader_read_evth(struct mliEventTapeReader *tio, float *evth)
{
        const uint64_t NUM_DIGITS = 9;
        const uint64_t BASE = 10;
        uint64_t path_event_number;
        uint64_t evth_event_number;
        uint64_t path_run_number;
        uint64_t evth_run_number;
        char match[MLI_TAR_NAME_LENGTH] = "ddddddddd/ddddddddd/EVTH.float32";

        if (!tio->has_tarh) {
                return 0;
        }
        if (!mli_cstr_match_templeate(tio->tarh.name, match, 'd')) {
                return 0;
        }
        chk_msg(mli_cstr_nto_uint64(
                        &path_event_number,
                        &tio->tarh.name[10],
                        BASE,
                        NUM_DIGITS),
                "Can't parse event-number from path.");
        chk_msg(mli_cstr_nto_uint64(
                        &path_run_number, &tio->tarh.name[0], BASE, NUM_DIGITS),
                "Can't parse run-number from path.");
        chk_msg(tio->tarh.size == MLI_CORSIKA_HEADER_SIZE_BYTES,
                "Expected EVTH to have size 273*sizeof(float)");
        chk_msg(mli_Tar_read_data(&tio->tar, (void *)evth, tio->tarh.size),
                "Can't read EVTH from tar.");
        chk_msg(evth[0] == mli_chars_to_float("EVTH"),
                "Expected EVTH[0] == 'EVTH'");

        evth_event_number = (uint64_t)evth[MLI_CORSIKA_EVTH_EVENT_NUMBER];
        evth_run_number = (uint64_t)evth[MLI_CORSIKA_EVTH_RUN_NUMBER];

        chk_msg(evth_event_number == path_event_number,
                "Expected paths' event-number to match event-number in EVTH.");
        chk_msg(evth_run_number == path_run_number,
                "Expected paths' run-number to match run-number in EVTH.");

        tio->event_number = evth_event_number;
        tio->cherenkov_bunch_block_number = 1;

        /* now there must follow a cherenkov-bunch-block */
        tio->has_tarh = mli_Tar_read_header(&tio->tar, &tio->tarh);

        chk_msg(tio->has_tarh, "Expected cherenkov-bunch-block after EVTH.");
        chk_msg(mliEventTapeReader_tarh_is_valid_cherenkov_block(tio),
                "Cherenkov-bunch-block's tar-header doesn't match.");

        chk_msg(tio->tarh.size % MLI_CORSIKA_BUNCH_SIZE_BYTES == 0,
                "Expected cherenkov-bunch-block-size "
                "to be multiple of bunch-size.");
        tio->block_size = tio->tarh.size / MLI_CORSIKA_BUNCH_SIZE_BYTES;
        tio->block_at = 0;
        tio->has_still_bunches_in_event = 1;

        return 1;
chk_error:
        return 0;
}

int mliEventTapeReader_tarh_might_be_valid_cherenkov_block(
        const struct mliEventTapeReader *tio)
{
        char match[MLI_TAR_NAME_LENGTH] =
                "ddddddddd/ddddddddd/ddddddddd.cer.x8.float32";
        return mli_cstr_match_templeate(tio->tarh.name, match, 'd');
}

int mliEventTapeReader_tarh_is_valid_cherenkov_block(
        const struct mliEventTapeReader *tio)
{
        const uint64_t NUM_DIGITS = 9;
        const uint64_t BASE = 10;
        uint64_t path_run_number;
        uint64_t path_event_number;
        uint64_t path_block_number;
        chk_msg(tio->has_tarh, "Expected a next tar-header.");

        chk_msg(mliEventTapeReader_tarh_might_be_valid_cherenkov_block(tio),
                "Expected cherenkov-bunch-block-name to be valid.");

        chk_msg(mli_cstr_nto_uint64(
                        &path_run_number, &tio->tarh.name[0], BASE, NUM_DIGITS),
                "Can't parse run-number from path.");
        chk_msg(path_run_number == tio->run_number,
                "Expected consistent run-number in cherenkov-block-path.");

        chk_msg(mli_cstr_nto_uint64(
                        &path_event_number,
                        &tio->tarh.name[9 + 1],
                        BASE,
                        NUM_DIGITS),
                "Can't parse event-number from path.");
        chk_msg(path_event_number == tio->event_number,
                "Expected same event-number in cherenkov-block-path and EVTH.");

        chk_msg(mli_cstr_nto_uint64(
                        &path_block_number,
                        &tio->tarh.name[9 + 1 + 9 + 1],
                        BASE,
                        NUM_DIGITS),
                "Can't parse cherenkov-block-number from path.");
        chk_msg(path_block_number == tio->cherenkov_bunch_block_number,
                "Expected different cherenkov-bunch-block-number in "
                "cherenkov-block-path.");
        return 1;
chk_error:
        return 0;
}

int mliEventTapeReader_read_cherenkov_bunch(
        struct mliEventTapeReader *tio,
        float *bunch)
{
        if (tio->has_still_bunches_in_event == 0) {
                return 0;
        }
        if (tio->block_at == tio->block_size) {
                tio->cherenkov_bunch_block_number += 1;
                tio->has_tarh = mli_Tar_read_header(&tio->tar, &tio->tarh);
                if (!tio->has_tarh) {
                        tio->has_still_bunches_in_event = 0;
                        return 0;
                }
                if (!mliEventTapeReader_tarh_might_be_valid_cherenkov_block(
                            tio)) {
                        tio->has_still_bunches_in_event = 0;
                        return 0;
                }
                chk_msg(mliEventTapeReader_tarh_is_valid_cherenkov_block(tio),
                        "Cherenkov-bunch-block's tar-header doesn't match.");
                chk_msg(tio->tarh.size % MLI_CORSIKA_BUNCH_SIZE_BYTES == 0,
                        "Expected cherenkov-bunch-block-size "
                        "to be multiple of bunch-size.");
                tio->block_size = tio->tarh.size / MLI_CORSIKA_BUNCH_SIZE_BYTES;
                tio->block_at = 0;
        }
        chk_msg(mli_Tar_read_data(
                        &tio->tar,
                        (void *)(bunch),
                        MLI_CORSIKA_BUNCH_SIZE_BYTES),
                "Failed to read cherenkov_bunch.");

        tio->block_at += 1;
        return 1;
chk_error:
        return 0;
}

/* EventTape_testing */
/* ----------------- */

/* Copyright 2020 Sebastian A. Mueller */

void mliEventTape_testing_set_random_corsika_header(
        float *head,
        struct mli_Prng *prng)
{
        uint64_t i;
        for (i = 0; i < 273; i++) {
                head[i] = (float)mli_Prng_uniform(prng);
        }
}

void mliEventTape_testing_set_random_RUNH(
        float *runh,
        const float run_number,
        struct mli_Prng *prng)
{
        mliEventTape_testing_set_random_corsika_header(runh, prng);
        runh[0] = mli_chars_to_float("RUNH");
        runh[MLI_CORSIKA_RUNH_RUN_NUMBER] = run_number;
}

void mliEventTape_testing_set_random_EVTH(
        float *evth,
        const float event_number,
        const float run_number,
        struct mli_Prng *prng)
{
        mliEventTape_testing_set_random_corsika_header(evth, prng);
        evth[0] = mli_chars_to_float("EVTH");
        evth[MLI_CORSIKA_EVTH_EVENT_NUMBER] = event_number;
        evth[MLI_CORSIKA_EVTH_RUN_NUMBER] = run_number;
}

void mliEventTape_testing_set_random_bunch(float *bunch, struct mli_Prng *prng)
{
        bunch[0] = (float)mli_Prng_uniform(prng);
        bunch[1] = (float)mli_Prng_uniform(prng);
        bunch[2] = (float)mli_Prng_uniform(prng);
        bunch[3] = (float)mli_Prng_uniform(prng);
        bunch[4] = (float)mli_Prng_uniform(prng);
        bunch[5] = (float)mli_Prng_uniform(prng);
        bunch[6] = (float)mli_Prng_uniform(prng);
        bunch[7] = (float)mli_Prng_uniform(prng);
}

int mliEventTape_testing_bunches_are_equal(float *b1, float *b2)
{
        if (b1[0] != b2[0]) {
                fprintf(stderr, "Bunch missmatch x_cm.\n");
                return 0;
        }
        if (b1[1] != b2[1]) {
                fprintf(stderr, "Bunch missmatch y_cm.\n");
                return 0;
        }
        if (b1[2] != b2[2]) {
                fprintf(stderr, "Bunch missmatch ux.\n");
                return 0;
        }
        if (b1[3] != b2[3]) {
                fprintf(stderr, "Bunch missmatch vy.\n");
                return 0;
        }
        if (b1[4] != b2[4]) {
                fprintf(stderr, "Bunch missmatch time_ns.\n");
                return 0;
        }
        if (b1[5] != b2[5]) {
                fprintf(stderr, "Bunch missmatch z_emission_cm.\n");
                return 0;
        }
        if (b1[6] != b2[6]) {
                fprintf(stderr, "Bunch missmatch weight_photons.\n");
                return 0;
        }
        if (b1[7] != b2[7]) {
                fprintf(stderr, "Bunch missmatch wavelength_nm.\n");
                return 0;
        }
        return 1;
}

int mliEventTape_testing_corsika_headers_are_equal(
        const float *h1,
        const float *h2)
{
        int i;
        for (i = 0; i < 273; i++) {
                if (h1[i] != h2[i]) {
                        fprintf(stderr, "Corsika-header missmatch at %d.\n", i);
                        return 0;
                }
        }
        return 1;
}

int mliEventTape_testing_write_and_read(
        const char *path,
        const uint64_t num_events,
        const uint64_t buffer_size,
        const float *event_numbers,
        const uint64_t *num_bunches,
        const uint32_t random_seed)
{
        uint64_t e, b;
        float corho[273] = {0.0};
        float corhi[273] = {0.0};
        float buncho[8] = {0.0};
        float bunchi[8] = {0.0};

        struct mli_IO ostream = mli_IO_init();
        struct mli_IO istream = mli_IO_init();
        struct mliEventTapeWriter taro = mliEventTapeWriter_init();
        struct mliEventTapeReader tari = mliEventTapeReader_init();

        struct mli_Prng prng = mli_Prng_init_PCG32(random_seed);
        mli_Prng_reinit(&prng, random_seed);

        /* write RUN */
        /* ========= */
        chk(mli_IO__open_file_cstr(&ostream, path, "w"));
        chk_msg(mliEventTapeWriter_begin(&taro, &ostream, buffer_size),
                "Can't begin writer.");
        /* set RUNH */
        mliEventTape_testing_set_random_RUNH(corho, 18.0, &prng);
        chk_msg(mliEventTapeWriter_write_runh(&taro, corho),
                "Can't write RUNH.");

        for (e = 0; e < num_events; e++) {
                /* set EVTH */
                mliEventTape_testing_set_random_EVTH(
                        corho, event_numbers[e], 18.0, &prng);
                chk_msg(mliEventTapeWriter_write_evth(&taro, corho),
                        "Can't write EVTH.");
                for (b = 0; b < num_bunches[e]; b++) {
                        mliEventTape_testing_set_random_bunch(buncho, &prng);
                        chk_msg(mliEventTapeWriter_write_cherenkov_bunch(
                                        &taro, buncho),
                                "Can't write bunch.");
                }
        }

        chk_msg(mliEventTapeWriter_finalize(&taro), "Can't finalize writer.");
        mli_IO_close(&ostream);

        /* read RUN */
        /* ======== */
        mli_Prng_reinit(&prng, random_seed);

        chk(mli_IO__open_file_cstr(&istream, path, "r"));
        chk_msg(mliEventTapeReader_begin(&tari, &istream),
                "Can't begin reader.");

        /* check RUNH */
        mliEventTape_testing_set_random_RUNH(corho, 18.0, &prng);
        chk_msg(mliEventTapeReader_read_runh(&tari, corhi), "Can't read RUNH.");
        chk_msg(mliEventTape_testing_corsika_headers_are_equal(corho, corhi),
                "Expected RUNH to be equal.");

        chk_msg(tari.event_number == 0,
                "Expected reader's event-number == 0 "
                "before first EVTH is read.");

        for (e = 0; e < num_events; e++) {
                chk_msg(mliEventTapeReader_read_evth(&tari, corhi),
                        "Can't read EVTH.");
                /* check EVTH */
                mliEventTape_testing_set_random_EVTH(
                        corho, event_numbers[e], 18.0, &prng);
                chk_msg(mliEventTape_testing_corsika_headers_are_equal(
                                corho, corhi),
                        "Expected EVTH to be equal.");
                chk_msg(corhi[MLI_CORSIKA_EVTH_EVENT_NUMBER] ==
                                event_numbers[e],
                        "Expected reader's event-number to match last EVTH.");
                chk_msg(tari.event_number == event_numbers[e],
                        "Expected a different event-number.");

                for (b = 0; b < num_bunches[e]; b++) {
                        chk_msg(mliEventTapeReader_read_cherenkov_bunch(
                                        &tari, bunchi),
                                "Can't read bunch.");
                        mliEventTape_testing_set_random_bunch(buncho, &prng);
                        chk_msg(mliEventTape_testing_bunches_are_equal(
                                        bunchi, buncho),
                                "Expected bunch to be equal.");
                }
                chk_msg(!mliEventTapeReader_read_cherenkov_bunch(&tari, buncho),
                        "Did not expect another cherenkov-bunch.");
        }
        chk_msg(!mliEventTapeReader_read_evth(&tari, corho),
                "Did not expect another EVTH.");

        chk_msg(mliEventTapeReader_finalize(&tari), "Can't finalize reader.");
        mli_IO_close(&istream);
        return 1;
chk_error:
        return 0;
}

/* Histogram2d */
/* ----------- */

/* Copyright 2017 Sebastian A. Mueller */


MLI_VECTOR_IMPLEMENTATION(
        mliDynCorsikaHistogram2dBin,
        struct mli_corsika_Histogram2dBin)

struct key {
        int32_t x;
        int32_t y;
};

union i4i4_to_i8 {
        struct key i4i4;
        int64_t i8;
};

struct mli_corsika_Histogram2d mli_corsika_Histogram2d_init(void)
{
        struct mli_corsika_Histogram2d hist;
        hist.dict = mli_AvlDict_init();
        return hist;
}

void mli_corsika_Histogram2d_free(struct mli_corsika_Histogram2d *hist)
{
        mli_AvlDict_free(&hist->dict);
}

int mli_corsika_Histogram2d_malloc(
        struct mli_corsika_Histogram2d *hist,
        const uint64_t capacity)
{
        return mli_AvlDict_malloc(&hist->dict, capacity);
}

int mli_corsika_Histogram2d_assign(
        struct mli_corsika_Histogram2d *hist,
        const int32_t x,
        const int32_t y,
        const double weight)
{
        int has;
        union i4i4_to_i8 key;
        int64_t ival = 0;
        key.i4i4.x = x;
        key.i4i4.y = y;

        has = mli_AvlDict_get(&hist->dict, key.i8, &ival);
        if (has) {
                double dval = mli_math_interpret_int64_as_double(ival);
                dval += weight;
                ival = mli_math_interpret_double_as_int64(dval);

        } else {
                ival = mli_math_interpret_double_as_int64(weight);
        }
        return mli_AvlDict_set(&hist->dict, key.i8, ival);
}

uint64_t mli_corsika_Histogram2d_len(const struct mli_corsika_Histogram2d *hist)
{
        return hist->dict.len;
}

int mli_corsika_Histogram2d_flatten__(
        const struct mli_AvlNode *node,
        struct mliDynCorsikaHistogram2dBin *f)
{
        if (node == NULL) {
                return 1;
        } else {
                union i4i4_to_i8 key;
                struct mli_corsika_Histogram2dBin bin;
                key.i8 = node->key;

                bin.x = key.i4i4.x;
                bin.y = key.i4i4.y;
                bin.value = mli_math_interpret_int64_as_double(node->value);

                chk_msg(mliDynCorsikaHistogram2dBin_push_back(f, bin),
                        "Failed to push back bin-node.");

                if (node->avl.left != NULL) {
                        struct mli_AvlNode *left =
                                (struct mli_AvlNode *)(node->avl.left);
                        chk_msg(mli_corsika_Histogram2d_flatten__(left, f),
                                "Failed left");
                }
                if (node->avl.right != NULL) {
                        struct mli_AvlNode *right =
                                (struct mli_AvlNode *)(node->avl.right);
                        chk_msg(mli_corsika_Histogram2d_flatten__(right, f),
                                "Failed right");
                }
                return 1;
        }
        return 1;
chk_error:
        return 0;
}

int mli_corsika_Histogram2d_flatten(
        const struct mli_corsika_Histogram2d *hist,
        struct mliDynCorsikaHistogram2dBin *f)
{
        chk_msg(mli_corsika_Histogram2d_flatten__(
                        (const struct mli_AvlNode *)hist->dict.tree.root, f),
                "Failed to write dict.");

        return 1;
chk_error:
        return 0;
}

void mli_corsika_Histogram2d_reset(struct mli_corsika_Histogram2d *hist)
{
        mli_AvlDict_reset(&hist->dict);
}

/* aabb */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_AABB mli_AABB_set(
        const struct mli_Vec lower,
        const struct mli_Vec upper)
{
        struct mli_AABB a;
        a.lower = lower;
        a.upper = upper;
        return a;
}

struct mli_AABB mli_AABB_outermost(
        const struct mli_AABB a,
        const struct mli_AABB b)
{
        struct mli_AABB c;
        c.lower.x = MLI_MATH_MIN2(a.lower.x, b.lower.x);
        c.lower.y = MLI_MATH_MIN2(a.lower.y, b.lower.y);
        c.lower.z = MLI_MATH_MIN2(a.lower.z, b.lower.z);
        c.upper.x = MLI_MATH_MAX2(a.upper.x, b.upper.x);
        c.upper.y = MLI_MATH_MAX2(a.upper.y, b.upper.y);
        c.upper.z = MLI_MATH_MAX2(a.upper.z, b.upper.z);
        return c;
}

struct mli_Vec mli_AABB_center(const struct mli_AABB a)
{
        struct mli_Vec sum = mli_Vec_add(a.upper, a.lower);
        return mli_Vec_multiply(sum, .5);
}

int mli_AABB_valid(const struct mli_AABB a)
{
        chk_msg(!MLI_MATH_IS_NAN(a.lower.x), "aabb.lower.x is 'nan'.");
        chk_msg(!MLI_MATH_IS_NAN(a.lower.y), "aabb.lower.y is 'nan'.");
        chk_msg(!MLI_MATH_IS_NAN(a.lower.z), "aabb.lower.z is 'nan'.");

        chk_msg(!MLI_MATH_IS_NAN(a.upper.x), "aabb.upper.x is 'nan'.");
        chk_msg(!MLI_MATH_IS_NAN(a.upper.y), "aabb.upper.y is 'nan'.");
        chk_msg(!MLI_MATH_IS_NAN(a.upper.z), "aabb.upper.z is 'nan'.");

        chk_msg(a.lower.x <= a.upper.x, "Expected lower.x <= upper.x");
        chk_msg(a.lower.y <= a.upper.y, "Expected lower.y <= upper.y");
        chk_msg(a.lower.z <= a.upper.z, "Expected lower.z <= upper.z");
        return 1;
chk_error:
        return 0;
}

int mli_AABB_equal(const struct mli_AABB a, const struct mli_AABB b)
{
        chk_msg(mli_Vec_equal(a.lower, b.lower),
                "Expected 'lower'-corner to be equal.");
        chk_msg(mli_Vec_equal(a.upper, b.upper),
                "Expected 'upper'-corner to be equal.");
        return 1;
chk_error:
        return 0;
}

int mli_AABB_is_overlapping(const struct mli_AABB a, const struct mli_AABB b)
{
        const int over_x = (a.upper.x >= b.lower.x) && (b.upper.x >= a.lower.x);
        const int over_y = (a.upper.y >= b.lower.y) && (b.upper.y >= a.lower.y);
        const int over_z = (a.upper.z >= b.lower.z) && (b.upper.z >= a.lower.z);
        return (over_x && over_y) && over_z;
}

int mli_AABB_is_point_inside(
        const struct mli_AABB a,
        const struct mli_Vec point)
{
        if (a.lower.x > point.x || a.upper.x <= point.x)
                return 0;
        if (a.lower.y > point.y || a.upper.y <= point.y)
                return 0;
        if (a.lower.z > point.z || a.upper.z <= point.z)
                return 0;
        return 1;
}

/* accelerator */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Accelerator mli_Accelerator_init(void)
{
        struct mli_Accelerator out;

        out.num_objects = 0u;
        out.object_octrees = NULL;

        out.num_robjects = 0u;
        out.robject_aabbs = NULL;

        out.scenery_octree = mli_OcTree_init();

        return out;
}

void mli_Accelerator_free(struct mli_Accelerator *self)
{
        uint32_t obj;

        mli_OcTree_free(&self->scenery_octree);

        for (obj = 0; obj < self->num_objects; obj++) {
                mli_OcTree_free(&self->object_octrees[obj]);
        }
        free(self->object_octrees);

        free(self->robject_aabbs);
        (*self) = mli_Accelerator_init();
}

int mli_Accelerator_malloc(
        struct mli_Accelerator *self,
        const uint32_t num_objects,
        const uint32_t num_robjects)
{
        uint32_t obj, rob;
        mli_Accelerator_free(self);

        self->num_objects = num_objects;
        chk_malloc(self->object_octrees, struct mli_OcTree, self->num_objects);
        for (obj = 0; obj < self->num_objects; obj++) {
                self->object_octrees[obj] = mli_OcTree_init();
        }

        self->num_robjects = num_robjects;
        chk_malloc(self->robject_aabbs, struct mli_AABB, self->num_robjects);
        for (rob = 0; rob < self->num_robjects; rob++) {
                self->robject_aabbs[rob] = mli_AABB_set(
                        mli_Vec_init(0.0, 0.0, 0.0),
                        mli_Vec_init(0.0, 0.0, 0.0));
        }

        return 1;
chk_error:
        return 0;
}

int mli_Accelerator_set_robject_aabbs(
        struct mli_Accelerator *self,
        const struct mli_Geometry *geometry)
{
        uint32_t rob;
        chk_msg(self->num_robjects == geometry->num_robjects,
                "Expected num_robjects to be equal, but its not.");

        for (rob = 0; rob < self->num_robjects; rob++) {
                const uint32_t robject = geometry->robjects[rob];
                self->robject_aabbs[rob] = mli_Object_aabb(
                        &geometry->objects[robject],
                        mli_HomTraComp_from_compact(
                                geometry->robject2root[rob]));
        }
        return 1;
chk_error:
        return 0;
}

int mli_Accelerator_set_object_octrees(
        struct mli_Accelerator *self,
        const struct mli_Geometry *geometry)
{
        uint32_t obj;
        chk_msg(self->num_objects == geometry->num_objects,
                "Expected num_objects to be equal, but its not.");

        for (obj = 0; obj < self->num_objects; obj++) {
                chk_msg(mli_OcTree_malloc_from_object_wavefront(
                                &self->object_octrees[obj],
                                &geometry->objects[obj]),
                        "Failed to setup mli_OcTree for object-wavefront.");
        }

        return 1;
chk_error:
        return 0;
}

int mli_Accelerator_malloc_from_Geometry(
        struct mli_Accelerator *self,
        const struct mli_Geometry *geometry)
{
        struct mli_AABB outermost_aabb;
        struct mli_GeometryAndAccelerator accgeo;
        accgeo.accelerator = self;
        accgeo.geometry = geometry;

        chk_msg(mli_Accelerator_malloc(
                        self, geometry->num_objects, geometry->num_robjects),
                "Failed to malloc mli_Accelerator from mli_Geometry's "
                "num_robjects");

        chk_msg(mli_Accelerator_set_robject_aabbs(self, geometry),
                "Failed to set AABBs of robjects.");

        chk_msg(mli_Accelerator_set_object_octrees(self, geometry),
                "Failed to setup object octrees.");

        outermost_aabb = mli_Accelerator_outermost_aabb(self);

        chk_msg(mli_OcTree_malloc_from_Geometry(
                        &self->scenery_octree, &accgeo, outermost_aabb),
                "Failed to set up octree across all robjects in geometry.");

        return 1;
chk_error:
        return 0;
}

void mli_Accelerator_info_fprint(FILE *f, const struct mli_Accelerator *self)
{
        uint32_t rob, i;
        fprintf(f, "accelerator\n");
        fprintf(f, "-----------\n");
        fprintf(f, "\n");
        fprintf(f, "    object-references bounding-boxes (AABB)s:\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        fprintf(f, "    ");
        fprintf(f, "%5s ", "ref");
        fprintf(f, "%9s ", "-x/m");
        fprintf(f, "%9s ", "-y/m");
        fprintf(f, "%9s ", "-z/m");
        fprintf(f, "%9s ", "+x/m");
        fprintf(f, "%9s ", "+y/m");
        fprintf(f, "%9s ", "+z/m");
        fprintf(f, "\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");

        for (rob = 0; rob < self->num_robjects; rob++) {
                fprintf(f, "    ");
                fprintf(f, "%5d ", rob);
                fprintf(f, "%9.1f ", self->robject_aabbs[rob].lower.x);
                fprintf(f, "%9.1f ", self->robject_aabbs[rob].lower.y);
                fprintf(f, "%9.1f ", self->robject_aabbs[rob].lower.z);
                fprintf(f, "%9.1f ", self->robject_aabbs[rob].upper.x);
                fprintf(f, "%9.1f ", self->robject_aabbs[rob].upper.y);
                fprintf(f, "%9.1f ", self->robject_aabbs[rob].upper.z);
                fprintf(f, "\n");
        }
}

struct mli_AABB mli_Accelerator_outermost_aabb(
        const struct mli_Accelerator *self)
{
        uint32_t rob;
        struct mli_AABB aabb;
        if (self->num_robjects == 0) {
                aabb.lower =
                        mli_Vec_init(MLI_MATH_NAN, MLI_MATH_NAN, MLI_MATH_NAN);
                aabb.upper =
                        mli_Vec_init(MLI_MATH_NAN, MLI_MATH_NAN, MLI_MATH_NAN);
                return aabb;
        }
        aabb.lower = self->robject_aabbs[0].lower;
        aabb.upper = self->robject_aabbs[0].upper;
        for (rob = 0; rob < self->num_robjects; rob++) {
                aabb = mli_AABB_outermost(aabb, self->robject_aabbs[rob]);
        }
        return aabb;
}

/* accelerator_equal */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Accelerator_equal(
        const struct mli_Accelerator *a,
        const struct mli_Accelerator *b)
{
        uint32_t i = 0u;
        chk_msg(a->num_objects == b->num_objects,
                "Expected num_objects to be equal.");
        for (i = 0; i < a->num_objects; i++) {
                chk_msg(mli_OcTree_equal(
                                &a->object_octrees[i], &b->object_octrees[i]),
                        "Expected object_octrees[i] to be equal.");
        }

        chk_msg(a->num_robjects == b->num_robjects,
                "Expected num_robjects to be equal.");
        for (i = 0; i < a->num_robjects; i++) {
                chk_msg(mli_AABB_equal(
                                a->robject_aabbs[i], b->robject_aabbs[i]),
                        "Expected robject_aabbs[i] to be equal.");
        }

        chk_msg(mli_OcTree_equal(&a->scenery_octree, &b->scenery_octree),
                "Expected scenery_octree to be equal.");

        return 1;
chk_error:
        return 0;
}

/* accelerator_serialize */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Accelerator_to_io(const struct mli_Accelerator *self, struct mli_IO *f)
{
        uint64_t i = 0;

        /* magic identifier */
        struct mli_MagicId magic = mli_MagicId_init();
        chk(mli_MagicId_set(&magic, "mli_Accelerator"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        /* capacity */
        chk_IO_write(&self->num_objects, sizeof(uint32_t), 1, f);
        chk_IO_write(&self->num_robjects, sizeof(uint32_t), 1, f);

        for (i = 0; i < self->num_objects; i++) {
                chk(mli_OcTree_to_io(&self->object_octrees[i], f));
        }
        chk_IO_write(
                self->robject_aabbs,
                sizeof(struct mli_AABB),
                self->num_robjects,
                f);
        chk(mli_OcTree_to_io(&self->scenery_octree, f));

        return 1;
chk_error:
        return 0;
}

int mli_Accelerator_from_io(struct mli_Accelerator *self, struct mli_IO *f)
{
        uint64_t i = 0u;
        struct mli_MagicId magic;

        uint32_t num_robjects = 0u;
        uint32_t num_objects = 0u;

        /* magic identifier */
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Accelerator"));
        mli_MagicId_warn_version(&magic);

        /* capacity */
        chk_IO_read(&num_objects, sizeof(uint32_t), 1u, f);
        chk_IO_read(&num_robjects, sizeof(uint32_t), 1u, f);

        /* malloc */
        chk_mem(mli_Accelerator_malloc(self, num_objects, num_robjects));

        for (i = 0; i < self->num_objects; i++) {
                chk_mem(mli_OcTree_from_io(&self->object_octrees[i], f));
        }

        chk_IO_read(
                self->robject_aabbs,
                sizeof(struct mli_AABB),
                self->num_robjects,
                f);

        chk_mem(mli_OcTree_from_io(&self->scenery_octree, f));

        return 1;
chk_error:
        return 0;
}

/* accelerator_valid */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Accelerator_valid(const struct mli_Accelerator *self)
{
        uint32_t i = 0u;
        for (i = 0u; i < self->num_objects; i++) {
                chk_msg(mli_OcTree_valid(&self->object_octrees[i]),
                        "Expected object_octrees[i] to be valid.");
        }
        for (i = 0u; i < self->num_robjects; i++) {
                chk_msg(mli_AABB_valid(self->robject_aabbs[i]),
                        "Expected robject_aabbs[i] to be valid.");
        }
        chk_msg(mli_OcTree_valid(&self->scenery_octree),
                "Expected scenery_octree to be valid.");
        return 1;
chk_error:
        return 0;
}

int mli_Accelerator_valid_wrt_Geometry(
        const struct mli_Accelerator *self,
        const struct mli_Geometry *geometry)
{
        uint32_t i = 0u;
        chk_msg(self->num_objects == geometry->num_objects,
                "Expected num_objects to be equal.");
        for (i = 0u; i < self->num_objects; i++) {
                chk_msg(mli_OcTree_valid_wrt_links(
                                &self->object_octrees[i],
                                geometry->objects[i].num_faces),
                        "Expected object_octrees[i] to be valid w.r.t. "
                        "the object's num_faces.");
        }

        chk_msg(self->num_robjects == geometry->num_robjects,
                "Expected num_robjects to be equal.");
        chk_msg(mli_OcTree_valid_wrt_links(
                        &self->scenery_octree, geometry->num_robjects),
                "Expected scenery_octree to be valid w.r.t. to "
                "geometry's num_robjects.");
        return 1;
chk_error:
        return 0;
}

/* aperture */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_camera_Aperture mli_camera_Aperture_init(void)
{
        const double mm = 1e-3;
        struct mli_camera_Aperture apcam;
        apcam.focal_length = 50.0 * mm;
        apcam.aperture_radius = apcam.focal_length / 2.0;
        apcam.image_sensor_distance = apcam.focal_length;
        apcam.image_sensor_width_x = 36.0 * mm;
        apcam.image_sensor_width_y = 24.0 * mm;
        return apcam;
}

struct mli_Vec mli_camera_Aperture_pixel_center_on_image_sensor_plane(
        const double image_sensor_width_x,
        const double image_sensor_width_y,
        const double image_sensor_distance,
        const uint64_t num_pixel_x,
        const uint64_t num_pixel_y,
        const uint64_t pixel_x,
        const uint64_t pixel_y)
{
        struct mli_Vec pixel_center;
        pixel_center.x = -1.0 * mli_math_bin_center_in_linear_space(
                                        -0.5 * image_sensor_width_x,
                                        +0.5 * image_sensor_width_x,
                                        num_pixel_x,
                                        pixel_x);
        pixel_center.y = -1.0 * mli_math_bin_center_in_linear_space(
                                        -0.5 * image_sensor_width_y,
                                        +0.5 * image_sensor_width_y,
                                        num_pixel_y,
                                        pixel_y);
        pixel_center.z = -image_sensor_distance;
        return pixel_center;
}

struct mli_Vec mli_camera_Aperture_pixel_support_on_image_sensor_plane(
        const double image_sensor_width_x,
        const double image_sensor_width_y,
        const double image_sensor_distance,
        const uint64_t num_pixel_x,
        const uint64_t num_pixel_y,
        const uint64_t pixel_x,
        const uint64_t pixel_y,
        struct mli_Prng *prng)
{
        double pixel_bin_width_x = image_sensor_width_x / (double)num_pixel_x;
        double pixel_bin_width_y = image_sensor_width_y / (double)num_pixel_y;
        struct mli_Vec support =
                mli_camera_Aperture_pixel_center_on_image_sensor_plane(
                        image_sensor_width_x,
                        image_sensor_width_y,
                        image_sensor_distance,
                        num_pixel_x,
                        num_pixel_y,
                        pixel_x,
                        pixel_y);
        double rnd_x =
                (mli_Prng_uniform(prng) * pixel_bin_width_x -
                 0.5 * pixel_bin_width_x);
        double rnd_y =
                (mli_Prng_uniform(prng) * pixel_bin_width_y -
                 0.5 * pixel_bin_width_y);
        support.x = support.x + rnd_x;
        support.y = support.y + rnd_y;
        return support;
}

struct mli_Vec mli_camera_Aperture_get_object_point(
        const double focal_length,
        const struct mli_Vec pixel_support)
{
        const double object_distance =
                mli_thin_lens_get_object_given_focal_and_image(
                        focal_length, -1.0 * pixel_support.z);
        const double scaleing = object_distance / pixel_support.z;
        return mli_Vec_init(
                scaleing * pixel_support.x,
                scaleing * pixel_support.y,
                object_distance);
}

struct mli_Vec mli_camera_Aperture_ray_support_on_aperture(
        const double aperture_radius,
        struct mli_Prng *prng)
{
        /* use a perfect disc as aperture */
        return mli_Vec_random_position_on_disc(aperture_radius, prng);
}

double mli_camera_Aperture_focal_length_given_field_of_view_and_sensor_width(
        const double field_of_view,
        const double image_sensor_width)
{
        const double image_sensor_radius = 0.5 * image_sensor_width;
        const double fov_opening_angle = 0.5 * field_of_view;
        return image_sensor_radius / tan(fov_opening_angle);
}

struct mli_Ray mli_camera_Aperture_get_ray_for_pixel(
        const double focal_length,
        const double aperture_radius,
        const double image_sensor_distance,
        const double image_sensor_width_x,
        const double image_sensor_width_y,
        const uint64_t num_pixel_x,
        const uint64_t num_pixel_y,
        const uint64_t pixel_x,
        const uint64_t pixel_y,
        struct mli_Prng *prng)
{
        struct mli_Vec direction;
        struct mli_Vec aperture_support =
                mli_camera_Aperture_ray_support_on_aperture(
                        aperture_radius, prng);

        struct mli_Vec image_sensor_support =
                (mli_camera_Aperture_pixel_support_on_image_sensor_plane(
                        image_sensor_width_x,
                        image_sensor_width_y,
                        image_sensor_distance,
                        num_pixel_x,
                        num_pixel_y,
                        pixel_x,
                        pixel_y,
                        prng));

        if (fabs(1.0 - focal_length / image_sensor_distance) < 1e-6) {
                /* focus set to infinity */
                direction = mli_Vec_multiply(image_sensor_support, -1.0);
        } else {
                struct mli_Vec object_point =
                        mli_camera_Aperture_get_object_point(
                                focal_length, image_sensor_support);

                direction = mli_Vec_substract(object_point, aperture_support);
        }

        return mli_Ray_set(aperture_support, direction);
}

void mli_camera_Aperture_aquire_pixels(
        const struct mli_camera_Aperture self,
        const struct mli_Image *image,
        const struct mli_HomTraComp camera2root_comp,
        const struct mli_Shader *tracer,
        const struct mli_image_PixelVector *pixels_to_do,
        struct mli_ColorVector *colors_to_do,
        struct mli_Prng *prng)
{
        uint64_t i;
        struct mli_HomTra camera2root =
                mli_HomTraComp_from_compact(camera2root_comp);

        colors_to_do->size = 0;
        for (i = 0; i < pixels_to_do->size; i++) {
                struct mli_Ray ray_wrt_camera =
                        mli_camera_Aperture_get_ray_for_pixel(
                                self.focal_length,
                                self.aperture_radius,
                                self.image_sensor_distance,
                                self.image_sensor_width_x,
                                self.image_sensor_width_y,
                                mli_Image_num_cols(image),
                                mli_Image_num_rows(image),
                                pixels_to_do->array[i].col,
                                pixels_to_do->array[i].row,
                                prng);

                struct mli_Ray ray_wrt_root =
                        mli_HomTraComp_ray(&camera2root, ray_wrt_camera);

                struct mli_Color set_color =
                        mli_Shader_trace_ray(tracer, ray_wrt_root, prng);

                mli_ColorVector_push_back(colors_to_do, set_color);
        }

        return;
}

void mli_camera_Aperture_assign_pixel_colors_to_sum_and_exposure_image(
        const struct mli_image_PixelVector *pixels,
        const struct mli_ColorVector *colors,
        struct mli_Image *sum_image,
        struct mli_Image *exposure_image)
{
        uint64_t i;
        const struct mli_Color ONE = mli_Color_set(1.0, 1.0, 1.0);
        for (i = 0; i < pixels->size; i++) {
                const struct mli_image_Pixel pix = pixels->array[i];
                const struct mli_Color color = colors->array[i];

                const struct mli_Color _exposure =
                        mli_Image_get_by_Pixel(exposure_image, pix);
                const struct mli_Color _sum =
                        mli_Image_get_by_Pixel(sum_image, pix);

                mli_Image_set_by_Pixel(
                        sum_image, pix, mli_Color_add(_sum, color));
                mli_Image_set_by_Pixel(
                        exposure_image, pix, mli_Color_add(_exposure, ONE));
        }
}

int mli_camera_Aperture_render_image(
        const struct mli_camera_Aperture self,
        const struct mli_HomTraComp camera2root_comp,
        const struct mli_Shader *tracer,
        struct mli_Image *image,
        struct mli_Prng *prng)
{
        float noise_threshold = 0.05 * 255.0;
        uint64_t MAX_ITERATIONS = 128;
        uint64_t iteration = 0;

        struct mli_Color zero_color = mli_Color_set(0.0, 0.0, 0.0);
        const struct mli_Color HIGH_COLOR = mli_Color_set(255.0, 255.0, 255.0);
        struct mli_Image sum_image = mli_Image_init();
        struct mli_Image exposure_image = mli_Image_init();
        struct mli_Image to_do_image = mli_Image_init();
        struct mli_Image sobel_image = mli_Image_init();
        struct mli_Image previous_sobel_image = mli_Image_init();
        struct mli_Image diff_image = mli_Image_init();
        struct mli_ColorVector colors_to_do = mli_ColorVector_init();
        struct mli_image_PixelVector pixels_to_do =
                mli_image_PixelVector_init();
        chk_msg(mli_Image_malloc_same_size(&sum_image, image),
                "Failed to malloc sum_image.");
        chk_msg(mli_Image_malloc_same_size(&exposure_image, image),
                "Failed to malloc exposure_image.");
        chk_msg(mli_Image_malloc_same_size(&to_do_image, image),
                "Failed to malloc to_do_image.");
        chk_msg(mli_Image_malloc_same_size(&sobel_image, image),
                "Failed to malloc sobel_image.");
        chk_msg(mli_Image_malloc_same_size(&previous_sobel_image, image),
                "Failed to malloc previous_sobel_image.");
        chk_msg(mli_Image_malloc_same_size(&diff_image, image),
                "Failed to malloc diff_image.");
        chk_msg(mli_image_PixelVector_malloc(
                        &pixels_to_do, mli_Image_num_pixel(image)),
                "Failed to malloc pixels_to_do.");
        chk_msg(mli_ColorVector_malloc(
                        &colors_to_do, mli_Image_num_pixel(image)),
                "Failed to malloc colors_to_do.");

        mli_Image_set_all(image, zero_color);
        mli_Image_set_all(&sum_image, zero_color);
        mli_Image_set_all(&exposure_image, zero_color);
        mli_Image_set_all(&to_do_image, zero_color);
        mli_Image_set_all(&sobel_image, zero_color);
        mli_Image_set_all(&previous_sobel_image, zero_color);
        MLI_MATH_ARRAY_SET(
                colors_to_do.array, zero_color, colors_to_do.capacity);

        /*
        initial image
        =============
        */
        mli_image_PixelVector_push_back_all_from_image(&pixels_to_do, image);
        mli_camera_Aperture_aquire_pixels(
                self,
                image,
                camera2root_comp,
                tracer,
                &pixels_to_do,
                &colors_to_do,
                prng);
        mli_camera_Aperture_assign_pixel_colors_to_sum_and_exposure_image(
                &pixels_to_do, &colors_to_do, &sum_image, &exposure_image);
        mli_Image_divide_pixelwise(&sum_image, &exposure_image, image);
        chk(mli_Image_sobel(image, &sobel_image));
        mli_Image_luminance_threshold_dilatation(
                &sobel_image, 128.0, HIGH_COLOR, &to_do_image);

        fprintf(stderr, "\n");
        while (1) {
                if (iteration >= MAX_ITERATIONS)
                        break;
                mli_image_PixelVector_above_threshold(
                        &to_do_image, 0.5, &pixels_to_do);
                if (pixels_to_do.size < mli_Image_num_pixel(image) / 100.0)
                        break;
                fprintf(stderr,
                        "loop %3u / %3u, %ld,%03ld pixel left\n",
                        (uint32_t)iteration + 1,
                        (uint32_t)MAX_ITERATIONS,
                        pixels_to_do.size / 1000,
                        pixels_to_do.size % 1000);
                mli_camera_Aperture_aquire_pixels(
                        self,
                        image,
                        camera2root_comp,
                        tracer,
                        &pixels_to_do,
                        &colors_to_do,
                        prng);
                mli_camera_Aperture_assign_pixel_colors_to_sum_and_exposure_image(
                        &pixels_to_do,
                        &colors_to_do,
                        &sum_image,
                        &exposure_image);
                mli_Image_divide_pixelwise(&sum_image, &exposure_image, image);
                chk(mli_Image_copy(&sobel_image, &previous_sobel_image));
                chk(mli_Image_sobel(image, &sobel_image));
                chk(mli_Image_fabs_difference(
                        &previous_sobel_image, &sobel_image, &diff_image));
                mli_Image_set_all(&to_do_image, zero_color);
                chk(mli_Image_luminance_threshold_dilatation(
                        &diff_image,
                        noise_threshold,
                        HIGH_COLOR,
                        &to_do_image));
                iteration += 1;
        }

        mli_Image_free(&sum_image);
        mli_Image_free(&exposure_image);
        mli_Image_free(&to_do_image);
        mli_Image_free(&sobel_image);
        mli_Image_free(&previous_sobel_image);
        mli_Image_free(&diff_image);
        mli_ColorVector_free(&colors_to_do);
        mli_image_PixelVector_free(&pixels_to_do);
        return 1;
chk_error:
        mli_Image_free(&sum_image);
        mli_Image_free(&exposure_image);
        mli_Image_free(&to_do_image);
        mli_Image_free(&sobel_image);
        mli_Image_free(&previous_sobel_image);
        mli_Image_free(&diff_image);
        mli_ColorVector_free(&colors_to_do);
        mli_image_PixelVector_free(&pixels_to_do);
        return 0;
}

/* archive */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Archive mli_Archive_init(void)
{
        struct mli_Archive out;
        out.textfiles = mli_StringVector_init();
        out.filenames = mli_Map_init();
        return out;
}

void mli_Archive_free(struct mli_Archive *self)
{
        uint64_t i;
        for (i = 0; i < self->textfiles.size; i++) {
                struct mli_String *str = &self->textfiles.array[i];
                mli_String_free(str);
        }
        mli_StringVector_free(&self->textfiles);
        mli_Map_free(&self->filenames);
        (*self) = mli_Archive_init();
}

int mli_Archive_malloc(struct mli_Archive *self)
{
        mli_Archive_free(self);
        chk(mli_StringVector_malloc(&self->textfiles, 0u));
        chk(mli_Map_malloc(&self->filenames));
        return 1;
chk_error:
        return 0;
}

int mli_Archive_push_back(
        struct mli_Archive *self,
        const struct mli_String *filename,
        const struct mli_String *payload)
{
        uint64_t next;
        struct mli_String *text = NULL;
        chk_msg(filename->size < MLI_TAR_NAME_LENGTH,
                "Expected shorter filename.");
        next = mli_Map_size(&self->filenames);

        /* filename */
        /* ======== */
        chk_msg(mli_Map_insert(&self->filenames, filename, next),
                "Can not insert key.");

        /* payload */
        /* ======= */
        chk_msg(mli_StringVector_push_back(&self->textfiles, mli_String_init()),
                "Can not push back mli_String.");
        text = &self->textfiles.array[next];
        chk_msg(mli_String_copy(text, payload), "Can not copy payload.");
        return 1;
chk_error:
        return 0;
}

int mli_Archive_from_io(struct mli_Archive *self, struct mli_IO *f)
{
        struct mli_Tar tar = mli_Tar_init();
        struct mli_TarHeader tarh = mli_TarHeader_init();
        struct mli_String payload = mli_String_init();
        struct mli_String filename = mli_String_init();

        char tarh_name[MLI_TAR_NAME_LENGTH] = {'\0'};

        chk_msg(mli_Archive_malloc(self), "Can not malloc selfhive.");
        chk_msg(mli_Tar_read_begin(&tar, f), "Can't begin tar.");

        while (mli_Tar_read_header(&tar, &tarh)) {

                chk(mli_String_from_cstr(&filename, tarh.name));
                chk(mli_String_strip(&filename, &filename));
                chk(mli_path_strip_this_dir(&filename, &filename));

                chk_msg(mli_String_malloc(&payload, tarh.size),
                        "Can not allocate payload.");
                chk_msg(mli_Tar_read_data(
                                &tar, (void *)payload.array, tarh.size),
                        "Failed to read payload from tar into payload.");
                payload.size = strlen(payload.array);

                chk_msg(mli_String_convert_line_break_CRLF_CR_to_LF(
                                &payload, &payload),
                        "Failed to replace CRLF and CR linebreaks.");
                chk_msg(mli_cstr_assert_only_NUL_LF_TAB_controls(payload.array),
                        "Did not expect control codes other than "
                        "('\\n', '\\t', '\\0') in textfiles.");
                chk_msg(mli_Archive_push_back(self, &filename, &payload),
                        "Can not push back file into selfhive.");
        }

        chk_msg(mli_Tar_read_finalize(&tar), "Can't finalize reading tar.");
        mli_String_free(&payload);
        mli_String_free(&filename);
        return 1;
chk_error:
        fprintf(stderr, "tar->filename: '%s'.\n", tarh_name);
        mli_String_free(&payload);
        mli_String_free(&filename);
        mli_Archive_free(self);
        return 0;
}

int mli_Archive__from_path_cstr(struct mli_Archive *self, const char *path)
{
        struct mli_IO f = mli_IO_init();
        chk_msg(mli_IO__open_file_cstr(&f, path, "r"),
                "Cant open archive from path.");
        chk_msg(mli_Archive_from_io(self, &f),
                "Can't fread Archive from file.");
        mli_IO_close(&f);
        return 1;
chk_error:
        return 0;
}

int mli_Archive_has(
        const struct mli_Archive *self,
        const struct mli_String *filename)
{
        return mli_Map_has(&self->filenames, filename);
}

int mli_Archive_get(
        const struct mli_Archive *self,
        const struct mli_String *filename,
        struct mli_String **str)
{
        uint64_t idx;
        struct mli_String *txt = NULL;
        chk(mli_Map_find(&self->filenames, filename, &idx));
        txt = &self->textfiles.array[idx];
        (*str) = txt;
        return 1;
chk_error:
        return 0;
}

uint64_t mli_Archive_size(const struct mli_Archive *self)
{
        return mli_Map_size(&self->filenames);
}

uint64_t mli_Archive_num_filename_prefix_sufix(
        const struct mli_Archive *self,
        const char *prefix,
        const char *sufix)
{
        uint64_t i = 0;
        uint64_t match;
        uint64_t num_matches = 0;
        for (i = 0; i < self->textfiles.size; i++) {
                struct mliMapItem *item = &self->filenames.items.array[i];

                if (mli_String_starts_with_cstr(&item->key, prefix) &&
                    mli_String_ends_with_cstr(&item->key, sufix)) {
                        match = 1;
                } else {
                        match = 0;
                }

                if (match) {
                        num_matches++;
                }
        }
        return num_matches;
}

/* args */
/* ---- */

/* Copyright 2018-2024 Sebastian Achim Mueller */


int mli_StringVector_from_argc_argv(
        struct mli_StringVector *self,
        int argc,
        char *argv[])
{
        int i;
        mli_StringVector_free(self);
        chk_msg(argc >= 0, "Expected 'argc' >= 0.");
        chk_msg(mli_StringVector_malloc(self, argc), "Failed to malloc Argv.");
        for (i = 0; i < argc; i++) {
                struct mli_String *field = NULL;
                chk_msg(mli_StringVector_push_back(self, mli_String_init()),
                        "");
                field = &self->array[i];
                chk_msg(mli_String_from_cstr(field, argv[i]),
                        "Failed to malloc string in Argv.");
        }

        return 1;
chk_error:
        mli_StringVector_free(self);
        return 0;
}

/* array */
/* ----- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/* array_dummy_testing */
/* ------------------- */

/* Copyright Sebastian Achim Mueller */

MLI_ARRAY_IMPLEMENTATION(mli_ArrayTestingFloat, float)
MLI_ARRAY_IMPLEMENTATION_ZERO_TERMINATION(mli_ArrayTestingChar, char)

/* atmosphere */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
/* based on:
 * 2009-2016 Scratchapixel. Distributed under the terms of the
 * GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 */


struct mli_Atmosphere mli_Atmosphere_init(void)
{
        struct mli_Atmosphere atm;
        atm.sunLatitude = 0.0;
        atm.sunHourAngle = 12.0;
        mli_Atmosphere_set_sun_direction(
                &atm, atm.sunLatitude, atm.sunHourAngle);

        atm.sunDistance = 1.5e11;
        atm.sunRadius = 7e8;

        atm.altitude = 2300.0;
        atm.earthRadius = 6360e3;
        atm.atmosphereRadius = atm.earthRadius + 60e3;

        /* The height for the density to drop by 1 over e */
        atm.Height_Rayleigh = 7994.0;
        atm.Height_Mie = 1200.0;

        /*
        atm.beta_Rayleigh = mli_Color_set(3.8e-6, 13.5e-6, 33.1e-6);
        atm.beta_Mie = mli_Color_multiply(mli_Color_set(1.0, 1.0, 1.0), 41e-6);
        */

        mli_ColorSpectrum_set_beta_rayleigh(&atm.beta_Rayleigh_spectrum);
        mli_ColorSpectrum_set(&atm.beta_Mie_spectrum, 41e-6);
        mli_ColorSpectrum_set_radiance_of_black_body_W_per_m2_per_sr(
                &atm.sun_spectrum, 5000.0);

        atm.numSamples = 16;
        atm.numSamplesLight = 8;

        atm.power = 3000.0;

        return atm;
}

void mli_ColorSpectrum_set_beta_rayleigh(struct mli_ColorSpectrum *self)
{
        uint64_t i;
        struct mli_ColorSpectrumBinEdges edges =
                mli_ColorSpectrumBinEdges_init();
        const double beta = 1.36e-30;
        for (i = 0; i < MLI_COLORSPECTRUM_SIZE; i++) {
                const double wavelength =
                        0.5 * (edges.values[i] + edges.values[i + 1]);
                self->values[i] = beta / pow(wavelength, 4.0);
        }
}

void mli_Atmosphere_set_sun_direction(
        struct mli_Atmosphere *self,
        const double sunLatitude,
        const double sunHourAngle)
{
        self->sunHourAngle = sunHourAngle;
        self->sunLatitude = sunLatitude;

        {
                const double hours_rad =
                        MLI_MATH_PI +
                        2.0 * MLI_MATH_PI * self->sunHourAngle / 24.0;

                const struct mli_HomTraComp tc_latitude = mli_HomTraComp_set(
                        mli_Vec_init(0.0, 0.0, 0.0),
                        mli_Quaternion_set_tait_bryan(
                                self->sunLatitude, 0.0, 0.0));
                const struct mli_HomTraComp tc_hour = mli_HomTraComp_set(
                        mli_Vec_init(0.0, 0.0, 0.0),
                        mli_Quaternion_set_tait_bryan(0.0, hours_rad, 0.0));
                const struct mli_HomTraComp tc =
                        mli_HomTraComp_sequence(tc_latitude, tc_hour);
                const struct mli_HomTra t = mli_HomTraComp_from_compact(tc);
                const struct mli_Vec zenith = mli_Vec_init(0.0, 0.0, 1.0);
                self->sunDirection = mli_HomTraComp_dir(&t, zenith);
        }
}

struct mli_ColorSpectrum mli_Atmosphere_compute_depth(
        const struct mli_Atmosphere *self,
        const struct mli_Vec orig,
        const struct mli_Vec dir,
        double tmin,
        double tmax)
{
        uint64_t i, j;
        const double segmentLength = (tmax - tmin) / self->numSamples;
        double tCurrent = tmin;
        struct mli_ColorSpectrum sumR = mli_ColorSpectrum_init_zeros();
        struct mli_ColorSpectrum sumM = mli_ColorSpectrum_init_zeros();
        /* mie and rayleigh contribution */

        double opticalDepthR = 0.0;
        double opticalDepthM = 0.0;
        /*
         * mu in the paper which is the cosine of the angle between the sun
         * direction and the ray direction
         */
        const double mu = mli_Vec_dot(dir, self->sunDirection);
        const double phaseR = 3.f / (16.f * MLI_MATH_PI) * (1.0 + mu * mu);
        const double g = 0.76f;
        const double phaseM =
                3.f / (8.f * MLI_MATH_PI) *
                (((1.f - g * g) * (1.f + mu * mu)) /
                 ((2.f + g * g) * pow(1.f + g * g - 2.f * g * mu, 1.5f)));

        for (i = 0; i < self->numSamples; ++i) {
                const struct mli_Vec samplePosition = mli_Vec_add(
                        orig,
                        mli_Vec_multiply(
                                dir, (tCurrent + segmentLength * 0.5f)));
                const double height =
                        (mli_Vec_norm(samplePosition) - self->earthRadius);
                /* compute optical depth for light */
                const double hr =
                        segmentLength * exp(-height / self->Height_Rayleigh);
                const double hm =
                        segmentLength * exp(-height / self->Height_Mie);
                double t0Light = 0;
                double t1Light = 0;
                double segmentLengthLight = 0;
                double tCurrentLight = 0;
                double opticalDepthLightR = 0;
                double opticalDepthLightM = 0;

                opticalDepthR += hr;
                opticalDepthM += hm;
                /* light optical depth */
                mli_Ray_sphere_intersection(
                        samplePosition,
                        self->sunDirection,
                        self->atmosphereRadius,
                        &t0Light,
                        &t1Light);
                segmentLengthLight = t1Light / self->numSamplesLight;

                for (j = 0; j < self->numSamplesLight; ++j) {
                        const struct mli_Vec samplePositionLight = mli_Vec_add(
                                samplePosition,
                                mli_Vec_multiply(
                                        self->sunDirection,
                                        tCurrentLight +
                                                0.5 * segmentLengthLight));
                        const double heightLight =
                                mli_Vec_norm(samplePositionLight) -
                                self->earthRadius;

                        if (heightLight < 0)
                                break;

                        opticalDepthLightR +=
                                segmentLengthLight *
                                exp(-heightLight / self->Height_Rayleigh);
                        opticalDepthLightM +=
                                segmentLengthLight *
                                exp(-heightLight / self->Height_Mie);
                        tCurrentLight += segmentLengthLight;
                }

                if (j == self->numSamplesLight) {
                        const struct mli_ColorSpectrum tau =
                                mli_ColorSpectrum_add(
                                        mli_ColorSpectrum_multiply_scalar(
                                                self->beta_Rayleigh_spectrum,
                                                opticalDepthR +
                                                        opticalDepthLightR),
                                        mli_ColorSpectrum_multiply_scalar(
                                                self->beta_Mie_spectrum,
                                                1.1 * (opticalDepthM +
                                                       opticalDepthLightM)));
                        const struct mli_ColorSpectrum attenuation =
                                mli_ColorSpectrum_exp(tau, -1.0);

                        sumR = mli_ColorSpectrum_add(
                                sumR,
                                mli_ColorSpectrum_multiply_scalar(
                                        attenuation, hr));
                        sumM = mli_ColorSpectrum_add(
                                sumM,
                                mli_ColorSpectrum_multiply_scalar(
                                        attenuation, hm));
                }
                tCurrent += segmentLength;
        }

        return mli_ColorSpectrum_multiply_scalar(
                mli_ColorSpectrum_add(
                        mli_ColorSpectrum_multiply_scalar(
                                mli_ColorSpectrum_multiply(
                                        sumR, self->beta_Rayleigh_spectrum),
                                phaseR),
                        mli_ColorSpectrum_multiply_scalar(
                                mli_ColorSpectrum_multiply(
                                        sumM, self->beta_Mie_spectrum),
                                phaseM)),
                self->power);
}

struct mli_ColorSpectrum mli_Atmosphere_hit_outer_atmosphere(
        const struct mli_Atmosphere *self,
        const struct mli_Vec orig,
        const struct mli_Vec dir,
        double tmin,
        double tmax)
{
        double t_minus = -1.0;
        double t_plus = -1.0;
        int has_intersection = mli_Ray_sphere_intersection(
                orig, dir, self->atmosphereRadius, &t_minus, &t_plus);

        if (!has_intersection || t_plus < 0.0) {
                return mli_ColorSpectrum_init_zeros();
        }
        if (t_minus > tmin && t_minus > 0) {
                tmin = t_minus;
        }
        if (t_plus < tmax) {
                tmax = t_plus;
        }

        return mli_Atmosphere_compute_depth(self, orig, dir, tmin, tmax);
}

struct mli_ColorSpectrum mli_Atmosphere_hit_earth_body(
        const struct mli_Atmosphere *self,
        const struct mli_Vec orig,
        const struct mli_Vec dir)
{
        double t_minus = DBL_MAX;
        double t_plus = DBL_MAX;
        double t_max = DBL_MAX;
        int intersects_earth_body = mli_Ray_sphere_intersection(
                orig, dir, self->earthRadius, &t_minus, &t_plus);

        if (intersects_earth_body && t_minus > 0) {
                t_max = MLI_MATH_MAX2(0.0, t_minus);
        }

        return mli_Atmosphere_hit_outer_atmosphere(self, orig, dir, 0.0, t_max);
}

struct mli_ColorSpectrum mli_Atmosphere_query(
        const struct mli_Atmosphere *self,
        const struct mli_Vec orig,
        const struct mli_Vec dir)
{
        struct mli_Vec orig_up = orig;
        orig_up.z += (self->earthRadius + self->altitude);

        return mli_Atmosphere_hit_earth_body(self, orig_up, dir);
}

void mli_Atmosphere_increase_latitude(
        struct mli_Atmosphere *self,
        const double increment)
{
        if (self->sunLatitude + increment <= 0.5 * MLI_MATH_PI) {
                self->sunLatitude += increment;
                mli_Atmosphere_set_sun_direction(
                        self, self->sunLatitude, self->sunHourAngle);
        }
}

void mli_Atmosphere_decrease_latitude(
        struct mli_Atmosphere *self,
        const double increment)
{
        if (self->sunLatitude - increment >= -0.5 * MLI_MATH_PI) {
                self->sunLatitude -= increment;
                mli_Atmosphere_set_sun_direction(
                        self, self->sunLatitude, self->sunHourAngle);
        }
}

void mli_Atmosphere_increase_hours(
        struct mli_Atmosphere *self,
        const double increment)
{
        if (self->sunHourAngle + increment < 24.0) {
                self->sunHourAngle += increment;
                mli_Atmosphere_set_sun_direction(
                        self, self->sunLatitude, self->sunHourAngle);
        }
}

void mli_Atmosphere_decrease_hours(
        struct mli_Atmosphere *self,
        const double increment)
{
        if (self->sunHourAngle - increment >= 0.0) {
                self->sunHourAngle -= increment;
                mli_Atmosphere_set_sun_direction(
                        self, self->sunLatitude, self->sunHourAngle);
        }
}

void mli_Atmosphere_increase_altitude(
        struct mli_Atmosphere *self,
        const double factor)
{
        assert(factor > 1.0);
        self->altitude *= factor;
}

void mli_Atmosphere_decrease_altitude(
        struct mli_Atmosphere *self,
        const double factor)
{
        assert(factor < 1.0);
        assert(factor > 0.0);
        self->altitude *= factor;
}

/* atmosphere_json */
/* --------------- */

/* Copyright 2018-2021 Sebastian Achim Mueller */

int mli_Atmosphere_from_json_token(
        struct mli_Atmosphere *atm,
        const struct mli_Json *json,
        const uint64_t tkn)
{
        (*atm) = mli_Atmosphere_init();

        chk(mli_Json_double_by_key(
                json, tkn, &atm->sunLatitude, "sunLatitude"));
        chk(mli_Json_double_by_key(
                json, tkn, &atm->sunHourAngle, "sunHourAngle"));
        mli_Atmosphere_set_sun_direction(
                atm, atm->sunLatitude, atm->sunHourAngle);

        chk(mli_Json_double_by_key(
                json, tkn, &atm->sunDistance, "sunDistance"));
        chk_msg(atm->sunDistance > 0, "Expected atmosphere->sunDistance > 0.");
        chk(mli_Json_double_by_key(json, tkn, &atm->sunRadius, "sunRadius"));
        chk_msg(atm->sunRadius > 0, "Expected atmosphere->sunRadius > 0.");

        chk(mli_Json_double_by_key(
                json, tkn, &atm->earthRadius, "earthRadius"));
        chk_msg(atm->earthRadius > 0, "Expected atmosphere->earthRadius > 0.");
        chk(mli_Json_double_by_key(
                json, tkn, &atm->atmosphereRadius, "atmosphereRadius"));
        chk_msg(atm->atmosphereRadius > atm->earthRadius,
                "Expected atmosphere->atmosphereRadius > atm->earthRadius.");

        chk(mli_Json_double_by_key(
                json, tkn, &atm->Height_Rayleigh, "Height_Rayleigh"));
        chk(mli_Json_double_by_key(json, tkn, &atm->Height_Mie, "Height_Mie"));

        chk(mli_Json_uint64_by_key(json, tkn, &atm->numSamples, "numSamples"));
        chk_msg(atm->numSamples > 0, "Expected atmosphere->numSamples > 0.");
        chk(mli_Json_uint64_by_key(
                json, tkn, &atm->numSamplesLight, "numSamplesLight"));
        chk_msg(atm->numSamplesLight > 0,
                "Expected atmosphere->numSamplesLight > 0.");

        chk(mli_Json_double_by_key(json, tkn, &atm->power, "power"));
        chk_msg(atm->power > 0, "Expected atmosphere->power > 0.");
        chk(mli_Json_double_by_key(json, tkn, &atm->altitude, "altitude"));
        chk_msg(atm->altitude > 0, "Expected atmosphere->altitude > 0.");

        return 1;
chk_error:
        return 0;
}

/* avl_Dict */
/* -------- */


struct mli_AvlDict mli_AvlDict_init(void)
{
        struct mli_AvlDict dict;
        dict.tree.root = NULL;
        dict.tree.compare = mli_AvlNode_compare;

        dict.nodes = NULL;
        dict.capacity = 0u;
        dict.back = 0u;
        dict.len = 0u;
        return dict;
}

void mli_AvlDict_free(struct mli_AvlDict *dict)
{
        free(dict->nodes);
        (*dict) = mli_AvlDict_init();
}

int mli_AvlDict_malloc(struct mli_AvlDict *dict, const uint64_t capacity)
{
        mli_AvlDict_free(dict);
        dict->capacity = capacity;
        chk_malloc(dict->nodes, struct mli_AvlNode, dict->capacity);
        return 1;
chk_error:
        mli_AvlDict_free(dict);
        return 0;
}

struct mli_AvlNode *mli_AvlDict_find(
        struct mli_AvlDict *dict,
        const int64_t key)
{
        struct mli_AvlNode *out = NULL;
        struct mli_AvlNode probe;
        probe.key = key;
        out = (struct mli_AvlNode *)mli_AvlTree_find(
                &dict->tree, (const struct mli_Avl *)(&probe));

        return out;
}

int mli_AvlDict_update__(
        const struct mli_AvlNode *node,
        struct mli_AvlDict *out)
{
        if (node == NULL) {
                return 1;
        }
        chk_msg(mli_AvlDict_set(out, node->key, node->value),
                "Failed to insert key/value into destination dict while "
                "updating.");

        if (node->avl.left != NULL) {
                struct mli_AvlNode *left =
                        (struct mli_AvlNode *)(node->avl.left);
                chk_msg(mli_AvlDict_update__(left, out), "1");
        }
        if (node->avl.right != NULL) {
                struct mli_AvlNode *right =
                        (struct mli_AvlNode *)(node->avl.right);
                chk_msg(mli_AvlDict_update__(right, out), "2");
        }

        return 1;
chk_error:
        return 0;
}

void mli_AvlDict_swap(struct mli_AvlDict *a, struct mli_AvlDict *b)
{
        struct mli_AvlDict swap = (*a);
        (*a) = (*b);
        (*b) = swap;
}

int mli_AvlDict_grow(struct mli_AvlDict *dict)
{
        uint64_t new_capacity = (dict->capacity * 2);
        struct mli_AvlDict tmp = mli_AvlDict_init();

        chk_msg(mli_AvlDict_malloc(&tmp, new_capacity),
                "Failed to malloc bigger tmp dict in order to grow.");

        /* copy nodes */
        chk_msg(mli_AvlDict_update__(
                        (const struct mli_AvlNode *)dict->tree.root, &tmp),
                "Failed to copy nodes over to bigger tmp dict while growing.");

        mli_AvlDict_swap(dict, &tmp);
        mli_AvlDict_free(&tmp);
        return 1;
chk_error:
        return 0;
}

int mli_AvlDict_insert(
        struct mli_AvlDict *dict,
        const int64_t key,
        const int64_t value)
{
        int insert_rc;
        struct mli_AvlNode *back_node;

        if (dict->back == dict->capacity) {
                chk_msg(mli_AvlDict_grow(dict),
                        "AvlTree insertion did not work.");
        }

        dict->nodes[dict->back].key = key;
        dict->nodes[dict->back].value = value;
        back_node = &dict->nodes[dict->back];

        insert_rc =
                mli_AvlTree_insert(&dict->tree, (struct mli_Avl *)(back_node));
        chk_msg(insert_rc >= 0, "AvlTree insertion did not work.");

        dict->back += 1;
        dict->len += 1;

        return 1;
chk_error:
        return 0;
}

int mli_AvlDict_set(
        struct mli_AvlDict *dict,
        const int64_t key,
        const int64_t value)
{
        struct mli_AvlNode *nn = mli_AvlDict_find(dict, key);
        if (nn != NULL) {
                /* key is already in dict. Set it right here. */
                nn->value = value;
        } else {
                /* key is not yet in dict. Insert to grow memory. */
                chk_msg(mli_AvlDict_insert(dict, key, value),
                        "Can't insert key/value into mli_AvlDict.");
        }
        return 1;
chk_error:
        return 0;
}

int mli_AvlDict_shrink(struct mli_AvlDict *dict)
{
        uint64_t new_capacity = (dict->len * 3) / 2;
        struct mli_AvlDict tmp = mli_AvlDict_init();

        chk_msg(mli_AvlDict_malloc(&tmp, new_capacity),
                "Failed to malloc smaller tmp dict in order to shrink.");

        /* copy nodes */
        chk_msg(mli_AvlDict_update__(
                        (const struct mli_AvlNode *)dict->tree.root, &tmp),
                "Failed to copy nodes over to smaller tmp dict while "
                "shrinking.");

        mli_AvlDict_swap(dict, &tmp);
        mli_AvlDict_free(&tmp);
        return 1;
chk_error:
        return 0;
}

int mli_AvlDict_pop(struct mli_AvlDict *dict, const int64_t key)
{
        int rc_remove;
        struct mli_AvlNode *nn = mli_AvlDict_find(dict, key);
        chk_msg(nn != NULL, "key does not exist, and thus can not be removed");

        /* the key exists and can be removed */
        rc_remove = mli_AvlTree_remove(&dict->tree, (struct mli_Avl *)nn);
        chk_msg(rc_remove <= 0, "AvlTree remove did not work as expected.");

        chk_msg(dict->len > 0, "Expected len > 0.");
        dict->len -= 1;

        if (dict->len < (dict->capacity / 2)) {
                /* shrink */
                chk_msg(mli_AvlDict_shrink(dict), "Failed to shrink capacity.");
        }

        return 1;
chk_error:
        return 0;
}

int mli_AvlDict_has(struct mli_AvlDict *dict, const int64_t key)
{
        struct mli_AvlNode *nn = mli_AvlDict_find(dict, key);
        if (nn == NULL) {
                return 0;
        } else {
                return 1;
        }
}

int mli_AvlDict_get(struct mli_AvlDict *dict, const int64_t key, int64_t *value)
{
        struct mli_AvlNode *nn = mli_AvlDict_find(dict, key);
        if (nn == NULL) {
                return 0;
        } else {
                (*value) = nn->value;
                return 1;
        }
}

void mli_AvlDict_reset(struct mli_AvlDict *dict)
{
        dict->tree.root = NULL;
        dict->back = 0;
        dict->len = 0;
}

/* avl_Tree */
/* -------- */


/* Adopted from:
 *
 * ANSI C Library for maintainance of AVL Balanced Trees
 *
 * ref.:
 *  G. M. Adelson-Velskij & E. M. Landis
 *  Doklady Akad. Nauk SSSR 146 (1962), 263-266
 *
 * see also:
 *  D. E. Knuth: The Art of Computer Programming Vol.3 (Sorting and Searching)
 *
 * (C) 2000 Daniel Nagy, Budapest University of Technology and Economics
 * Released under GNU General Public License (GPL) version 2
 *
 */

/* Swing to the left
 * Warning: no balance maintainance
 */
void mli_Avl_swing_left(struct mli_Avl **root)
{
        /* no balance maintainance */
        struct mli_Avl *a = *root;
        struct mli_Avl *b = a->right;
        *root = b;
        a->right = b->left;
        b->left = a;
}

/* Swing to the right
 * Warning: no balance maintainance
 */
void mli_Avl_swing_right(struct mli_Avl **root)
{
        /* no balance maintainance */
        struct mli_Avl *a = *root;
        struct mli_Avl *b = a->left;
        *root = b;
        a->left = b->right;
        b->right = a;
}

/* Balance maintainance after especially nasty swings
 */
void mli_Avl_rebalance(struct mli_Avl *root)
{
        switch (root->balance) {
        case -1:
                root->left->balance = 0;
                root->right->balance = 1;
                break;
        case 1:
                root->left->balance = -1;
                root->right->balance = 0;
                break;
        case 0:
                root->left->balance = 0;
                root->right->balance = 0;
        }
        root->balance = 0;
}

/* Insert element a into the AVL tree t
 * returns 1 if the depth of the tree has grown
 * Warning: do not insert elements already present
 */
int mli_AvlTree_insert(struct mli_AvlTree *t, struct mli_Avl *a)
{
        /* initialize */
        a->left = 0;
        a->right = 0;
        a->balance = 0;

        /* insert into an empty tree */
        if (!t->root) {
                t->root = a;
                return 1;
        }

        if (t->compare(t->root, a) > 0) {
                /* insert into the left subtree */
                if (t->root->left) {
                        struct mli_AvlTree left_subtree;
                        left_subtree.root = t->root->left;
                        left_subtree.compare = t->compare;
                        if (mli_AvlTree_insert(&left_subtree, a)) {
                                switch (t->root->balance--) {
                                case 1:
                                        return 0;
                                case 0:
                                        return 1;
                                }
                                if (t->root->left->balance < 0) {
                                        mli_Avl_swing_right(&(t->root));
                                        t->root->balance = 0;
                                        t->root->right->balance = 0;
                                } else {
                                        mli_Avl_swing_left(&(t->root->left));
                                        mli_Avl_swing_right(&(t->root));
                                        mli_Avl_rebalance(t->root);
                                }
                        } else {
                                t->root->left = left_subtree.root;
                        }
                        return 0;
                } else {
                        t->root->left = a;
                        if (t->root->balance--) {
                                return 0;
                        }
                        return 1;
                }
        } else {
                /* insert into the right subtree */
                if (t->root->right) {
                        struct mli_AvlTree right_subtree;
                        right_subtree.root = t->root->right;
                        right_subtree.compare = t->compare;
                        if (mli_AvlTree_insert(&right_subtree, a)) {
                                switch (t->root->balance++) {
                                case -1:
                                        return 0;
                                case 0:
                                        return 1;
                                }
                                if (t->root->right->balance > 0) {
                                        mli_Avl_swing_left(&(t->root));
                                        t->root->balance = 0;
                                        t->root->left->balance = 0;
                                } else {
                                        mli_Avl_swing_right(&(t->root->right));
                                        mli_Avl_swing_left(&(t->root));
                                        mli_Avl_rebalance(t->root);
                                }
                        } else {
                                t->root->right = right_subtree.root;
                        }
                        return 0;
                } else {
                        t->root->right = a;
                        if (t->root->balance++) {
                                return 0;
                        }
                        return 1;
                }
        }
        return -1;
}

/* Remove an element a from the AVL tree t
 * returns -1 if the depth of the tree has shrunk
 * Warning: if the element is not present in the tree,
 *          returns 0 as if it had been removed succesfully.
 */
int mli_AvlTree_remove(struct mli_AvlTree *t, struct mli_Avl *a)
{
        int b;
        if (t->root == a) {
                return mli_AvlTree_removeroot(t);
        }
        b = t->compare(t->root, a);

        if (b >= 0) {
                /* remove from the left subtree */
                int ch;
                if (t->root->left) {
                        struct mli_AvlTree left_subtree;
                        left_subtree.root = t->root->left;
                        left_subtree.compare = t->compare;
                        ch = mli_AvlTree_remove(&left_subtree, a);
                        t->root->left = left_subtree.root;
                        if (ch) {
                                switch (t->root->balance++) {
                                case -1:
                                        return -1;
                                case 0:
                                        return 0;
                                }
                                switch (t->root->right->balance) {
                                case 0:
                                        mli_Avl_swing_left(&(t->root));
                                        t->root->balance = -1;
                                        t->root->left->balance = 1;
                                        return 0;
                                case 1:
                                        mli_Avl_swing_left(&(t->root));
                                        t->root->balance = 0;
                                        t->root->left->balance = 0;
                                        return -1;
                                }
                                mli_Avl_swing_right(&(t->root->right));
                                mli_Avl_swing_left(&(t->root));
                                mli_Avl_rebalance(t->root);
                                return -1;
                        }
                }
        }
        if (b <= 0) {
                /* remove from the right subtree */
                int ch;
                if (t->root->right) {
                        struct mli_AvlTree right_subtree;
                        right_subtree.root = t->root->right;
                        right_subtree.compare = t->compare;
                        ch = mli_AvlTree_remove(&right_subtree, a);
                        t->root->right = right_subtree.root;
                        if (ch) {
                                switch (t->root->balance--) {
                                case 1:
                                        return -1;
                                case 0:
                                        return 0;
                                }
                                switch (t->root->left->balance) {
                                case 0:
                                        mli_Avl_swing_right(&(t->root));
                                        t->root->balance = 1;
                                        t->root->right->balance = -1;
                                        return 0;
                                case -1:
                                        mli_Avl_swing_right(&(t->root));
                                        t->root->balance = 0;
                                        t->root->right->balance = 0;
                                        return -1;
                                }
                                mli_Avl_swing_left(&(t->root->left));
                                mli_Avl_swing_right(&(t->root));
                                mli_Avl_rebalance(t->root);
                                return -1;
                        }
                }
        }
        return 0;
}

/* Remove the root of the AVL tree t
 * Warning: dumps core if t is empty
 */
int mli_AvlTree_removeroot(struct mli_AvlTree *t)
{
        int ch;
        struct mli_Avl *a;
        if (!t->root->left) {
                if (!t->root->right) {
                        t->root = 0;
                        return -1;
                }
                t->root = t->root->right;
                return -1;
        }
        if (!t->root->right) {
                t->root = t->root->left;
                return -1;
        }
        if (t->root->balance < 0) {

                /* remove from the left subtree */
                a = t->root->left;
                while (a->right)
                        a = a->right;
        } else {
                /* remove from the right subtree */
                a = t->root->right;
                while (a->left) {
                        a = a->left;
                }
        }
        ch = mli_AvlTree_remove(t, a);
        a->left = t->root->left;
        a->right = t->root->right;
        a->balance = t->root->balance;
        t->root = a;
        if (a->balance == 0)
                return ch;
        return 0;
}

struct mli_Avl *mli_AvlTree_find(
        struct mli_AvlTree *t,
        const struct mli_Avl *probe)
{
        int64_t match;

        if (t->root == NULL) {
                return NULL;
        }

        match = t->compare(probe, t->root);

        if (match == 0) {
                return t->root;
        } else if (match < 0) {
                if (t->root->left != NULL) {
                        struct mli_AvlTree left_subtree;
                        left_subtree.root = t->root->left;
                        left_subtree.compare = t->compare;
                        return mli_AvlTree_find(&left_subtree, probe);
                } else {
                        return NULL;
                }
        } else {
                if (t->root->right != NULL) {
                        struct mli_AvlTree right_subtree;
                        right_subtree.root = t->root->right;
                        right_subtree.compare = t->compare;
                        return mli_AvlTree_find(&right_subtree, probe);
                } else {
                        return NULL;
                }
        }
        return NULL;
}

struct mli_AvlNode mli_AvlNode_init(void)
{
        struct mli_AvlNode n;
        n.avl.left = NULL;
        n.avl.right = NULL;
        n.avl.balance = 0;
        n.key = 0;
        n.value = 0;
        return n;
}

int64_t mli_AvlNode_compare(const void *a, const void *b)
{
        return ((struct mli_AvlNode *)a)->key - ((struct mli_AvlNode *)b)->key;
}

void mli_AvlNode_print(struct mli_Avl *a, int m)
{
        int n = m;
        if (a == NULL) {
                return;
        };
        if (a->right) {
                mli_AvlNode_print(a->right, m + 1);
        }
        while (n--) {
                printf("   ");
        }
        printf("%ld (%ld)\n", ((struct mli_AvlNode *)a)->key, a->balance);
        if (a->left) {
                mli_AvlNode_print(a->left, m + 1);
        }
}

/* axis_aligned_grid */
/* ----------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

/* Inspired by:
 * A Fast Voxel Traversal Algorithm for Ray Tracing
 * John Amanatides and Andrew Woo
 * Dept. of Computer Science
 * University of Toronto
 * Toronto, Ontario, Canada M5S 1A4
 */

struct mli_Idx3 mli_Idx3_set(const int64_t x, const int64_t y, const int64_t z)
{
        struct mli_Idx3 iii;
        iii.x = x;
        iii.y = y;
        iii.z = z;
        return iii;
}

struct mli_AxisAlignedGrid mli_AxisAlignedGrid_set(
        struct mli_AABB bounds,
        struct mli_Idx3 num_bins)
{
        struct mli_AxisAlignedGrid grid;
        grid.bounds = bounds;
        grid.num_bins = num_bins;
        assert(mli_AABB_valid(grid.bounds));
        assert(grid.num_bins.x > 0);
        assert(grid.num_bins.y > 0);
        assert(grid.num_bins.z > 0);
        grid.bin_width.x = (grid.bounds.upper.x - grid.bounds.lower.x) /
                           ((double)grid.num_bins.x);
        grid.bin_width.y = (grid.bounds.upper.y - grid.bounds.lower.y) /
                           ((double)grid.num_bins.y);
        grid.bin_width.z = (grid.bounds.upper.z - grid.bounds.lower.z) /
                           ((double)grid.num_bins.z);
        return grid;
}

struct mli_Idx3 mli_AxisAlignedGrid_get_voxel_idx(
        const struct mli_AxisAlignedGrid *grid,
        struct mli_Vec point)
{
        struct mli_Vec pg = mli_Vec_substract(point, grid->bounds.lower);
        struct mli_Idx3 iii;
        pg.x /= grid->bin_width.x;
        pg.y /= grid->bin_width.y;
        pg.z /= grid->bin_width.z;
        iii.x = (int64_t)floor(pg.x);
        iii.y = (int64_t)floor(pg.y);
        iii.z = (int64_t)floor(pg.z);
        return iii;
}

int mli_AxisAlignedGrid_find_voxel_of_first_interaction(
        const struct mli_AxisAlignedGrid *grid,
        const struct mli_Ray *ray,
        struct mli_Idx3 *bin)
{
        if (mli_AABB_is_point_inside(grid->bounds, ray->support)) {
                (*bin) = mli_AxisAlignedGrid_get_voxel_idx(grid, ray->support);
                return MLI_AXISALIGNEDGRID_RAY_STARTS_INSIDE_GRID;
        } else {
                double ray_parameter_near, ray_parameter_far;
                int has_intersection;
                mli_Ray_aabb_intersections(
                        (*ray),
                        grid->bounds,
                        &ray_parameter_near,
                        &ray_parameter_far);
                has_intersection =
                        mli_Ray_aabb_intersections_is_valid_given_near_and_far(
                                ray_parameter_near, ray_parameter_far);
                if (has_intersection) {
                        struct mli_Vec inner;
                        inner = mli_Ray_at(ray, ray_parameter_near);
                        (*bin) = mli_AxisAlignedGrid_get_voxel_idx(grid, inner);

                        if (bin->x >= grid->num_bins.x) {
                                bin->x -= 1;
                        }
                        if (bin->y >= grid->num_bins.y) {
                                bin->y -= 1;
                        }
                        if (bin->z >= grid->num_bins.z) {
                                bin->z -= 1;
                        }

                        return MLI_AXISALIGNEDGRID_RAY_STARTS_OUTSIDE_GRID_BUT_INTERSECTS;
                } else {
                        return MLI_AXISALIGNEDGRID_RAY_DOES_NOT_INTERSECT_GRID;
                }
        }
        return 0;
}

struct mli_Vec mli_AxisAlignedGridTraversal_first_plane(
        const struct mli_AxisAlignedGrid *grid,
        const struct mli_Idx3 voxel,
        const struct mli_Vec ray_direction)
{
        struct mli_Vec voxel_lower = mli_Vec_init(
                grid->bounds.lower.x + (double)voxel.x * grid->bin_width.x,
                grid->bounds.lower.y + (double)voxel.y * grid->bin_width.y,
                grid->bounds.lower.z + (double)voxel.z * grid->bin_width.z);
        struct mli_Vec voxel_upper = mli_Vec_init(
                grid->bounds.lower.x +
                        (double)(voxel.x + 1) * grid->bin_width.x,
                grid->bounds.lower.y +
                        (double)(voxel.y + 1) * grid->bin_width.y,
                grid->bounds.lower.z +
                        (double)(voxel.z + 1) * grid->bin_width.z);

        struct mli_Vec first;

        if (ray_direction.x >= 0.0) {
                first.x = voxel_upper.x;
        } else {
                first.x = voxel_lower.x;
        }

        if (ray_direction.y >= 0.0) {
                first.y = voxel_upper.y;
        } else {
                first.y = voxel_lower.y;
        }

        if (ray_direction.z >= 0.0) {
                first.z = voxel_upper.z;
        } else {
                first.z = voxel_lower.z;
        }

        return first;
}

double calc_t_for_x_plane(const double x_plane, const struct mli_Ray *ray)
{
        return -(ray->support.x - x_plane) / (ray->direction.x);
}

double calc_t_for_y_plane(const double y_plane, const struct mli_Ray *ray)
{
        return -(ray->support.y - y_plane) / (ray->direction.y);
}

double calc_t_for_z_plane(const double z_plane, const struct mli_Ray *ray)
{
        return -(ray->support.z - z_plane) / (ray->direction.z);
}

struct mli_AxisAlignedGridTraversal mli_AxisAlignedGridTraversal_start(
        const struct mli_AxisAlignedGrid *grid,
        const struct mli_Ray *ray)
{
        struct mli_AxisAlignedGridTraversal traversal;

        traversal.grid = grid;
        traversal.valid = mli_AxisAlignedGrid_find_voxel_of_first_interaction(
                grid, ray, &traversal.voxel);
        if (traversal.valid) {
                struct mli_Vec first_plane;
                traversal.step.x = MLI_MATH_SIGN(ray->direction.x);
                traversal.step.y = MLI_MATH_SIGN(ray->direction.y);
                traversal.step.z = MLI_MATH_SIGN(ray->direction.z);

                first_plane = mli_AxisAlignedGridTraversal_first_plane(
                        grid, traversal.voxel, ray->direction);

                traversal.tMax.x = calc_t_for_x_plane(first_plane.x, ray);
                traversal.tMax.y = calc_t_for_y_plane(first_plane.y, ray);
                traversal.tMax.z = calc_t_for_z_plane(first_plane.z, ray);

                traversal.tDelta.x = grid->bin_width.x / fabs(ray->direction.x);
                traversal.tDelta.y = grid->bin_width.y / fabs(ray->direction.y);
                traversal.tDelta.z = grid->bin_width.z / fabs(ray->direction.z);
        }
        return traversal;
}

int mli_AxisAlignedGridTraversal_next(
        struct mli_AxisAlignedGridTraversal *traversal)
{
        struct mli_AxisAlignedGridTraversal *t = traversal;
        int RAY_LEFT_GRID = 0;

        if (t->tMax.x < t->tMax.y) {
                if (t->tMax.x < t->tMax.z) {
                        t->voxel.x += t->step.x;
                        if (t->voxel.x < 0 ||
                            t->voxel.x >= t->grid->num_bins.x) {
                                traversal->valid = RAY_LEFT_GRID;
                                return traversal->valid;
                        }
                        t->tMax.x += t->tDelta.x;
                } else {
                        t->voxel.z += t->step.z;
                        if (t->voxel.z < 0 ||
                            t->voxel.z >= t->grid->num_bins.z) {
                                traversal->valid = RAY_LEFT_GRID;
                                return traversal->valid;
                        }
                        t->tMax.z += t->tDelta.z;
                }
        } else {
                if (t->tMax.y < t->tMax.z) {
                        t->voxel.y += t->step.y;
                        if (t->voxel.y < 0 ||
                            t->voxel.y >= t->grid->num_bins.y) {
                                traversal->valid = RAY_LEFT_GRID;
                                return traversal->valid;
                        }
                        t->tMax.y += t->tDelta.y;
                } else {
                        t->voxel.z += t->step.z;
                        if (t->voxel.z < 0 ||
                            t->voxel.z >= t->grid->num_bins.z) {
                                traversal->valid = RAY_LEFT_GRID;
                                return traversal->valid;
                        }
                        t->tMax.z += t->tDelta.z;
                }
        }
        return traversal->valid;
}

void mli_AxisAlignedGridTraversal_fprint(
        FILE *f,
        struct mli_AxisAlignedGridTraversal *traversal)
{
        struct mli_AxisAlignedGridTraversal *t = traversal;
        fprintf(f,
                "  grid.bounds.upper: [%f, %f, %f]\n",
                t->grid->bounds.upper.x,
                t->grid->bounds.upper.y,
                t->grid->bounds.upper.z);
        fprintf(f,
                "  grid.bounds.lower: [%f, %f, %f]\n",
                t->grid->bounds.lower.x,
                t->grid->bounds.lower.y,
                t->grid->bounds.lower.z);
        fprintf(f,
                "  grid.num_bins: [%ld, %ld, %ld]\n",
                t->grid->num_bins.x,
                t->grid->num_bins.y,
                t->grid->num_bins.z);

        fprintf(f, "  valid: %d\n", t->valid);
        fprintf(f,
                "  voxel: [%ld, %ld, %ld]\n",
                t->voxel.x,
                t->voxel.y,
                t->voxel.z);
        fprintf(f, "  step: [%f, %f, %f]\n", t->step.x, t->step.y, t->step.z);
        fprintf(f, "  tMax: [%f, %f, %f]\n", t->tMax.x, t->tMax.y, t->tMax.z);
        fprintf(f,
                "  tDelta: [%f, %f, %f]\n",
                t->tDelta.x,
                t->tDelta.y,
                t->tDelta.z);
}

void mli_Ray_fprint(FILE *f, struct mli_Ray *ray)
{
        fprintf(f,
                "[%f, %f, %f] + lam*[%f, %f, %f]",
                ray->support.x,
                ray->support.y,
                ray->support.z,
                ray->direction.x,
                ray->direction.y,
                ray->direction.z);
}

/* boundarylayer */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mli_BoundaryLayer_free(struct mli_BoundaryLayer *self)
{
        mli_String_free(&self->name);
}

struct mli_BoundaryLayer mli_BoundaryLayer_init(void)
{
        struct mli_BoundaryLayer layer;
        layer.inner.medium = 0;
        layer.inner.surface = 0;
        layer.outer.medium = 0;
        layer.outer.surface = 0;
        layer.name = mli_String_init();
        return layer;
}

int mli_BoundaryLayer_equal(
        const struct mli_BoundaryLayer *a,
        const struct mli_BoundaryLayer *b)
{
        if (!mli_String_equal(&a->name, &b->name)) {
                return 0;
        }
        if (a->inner.surface != b->inner.surface) {
                return 0;
        }
        if (a->inner.medium != b->inner.medium) {
                return 0;
        }
        if (a->outer.surface != b->outer.surface) {
                return 0;
        }
        if (a->outer.medium != b->outer.medium) {
                return 0;
        }
        return 1;
}

int mli_BoundaryLayer_to_io(
        const struct mli_BoundaryLayer *self,
        struct mli_IO *f)
{
        struct mli_MagicId magic = mli_MagicId_init();
        chk(mli_MagicId_set(&magic, "mli_BoundaryLayer"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        chk_msg(mli_String_to_io(&self->name, f),
                "Can't write boundary layer name to io.");
        chk_IO_write(&self->inner.medium, sizeof(uint64_t), 1u, f);
        chk_IO_write(&self->inner.surface, sizeof(uint64_t), 1u, f);
        chk_IO_write(&self->outer.medium, sizeof(uint64_t), 1u, f);
        chk_IO_write(&self->outer.surface, sizeof(uint64_t), 1u, f);

        return 1;
chk_error:
        return 0;
}

int mli_BoundaryLayer_from_io(struct mli_BoundaryLayer *self, struct mli_IO *f)
{
        struct mli_MagicId magic;
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_BoundaryLayer"));
        mli_MagicId_warn_version(&magic);

        chk_msg(mli_String_from_io(&self->name, f),
                "Can't read boundary layer name from io.");
        chk_IO_read(&self->inner.medium, sizeof(uint64_t), 1u, f);
        chk_IO_read(&self->inner.surface, sizeof(uint64_t), 1u, f);
        chk_IO_read(&self->outer.medium, sizeof(uint64_t), 1u, f);
        chk_IO_read(&self->outer.surface, sizeof(uint64_t), 1u, f);

        return 1;
chk_error:
        return 0;
}

/* boundarylayer_array */
/* ------------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
MLI_ARRAY_IMPLEMENTATION_FREE(
        mli_BoundaryLayerArray,
        struct mli_BoundaryLayer,
        mli_BoundaryLayer_free)

/* chk */
/* --- */

/* Copyright 2018-2021 Sebastian Achim Mueller */

int chk_eprintf(const char *format, ...)
{
        int r;
        va_list args;
        va_start(args, format);
        r = vfprintf(stderr, format, args);
        va_end(args);
        return r;
}

/* color */
/* ----- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Color mli_Color_set(const float r, const float g, const float b)
{
        struct mli_Color rgb;
        rgb.r = r;
        rgb.g = g;
        rgb.b = b;
        return rgb;
}

struct mli_Color mli_Color_mix(
        const struct mli_Color a,
        const struct mli_Color b,
        const float refl)
{
        struct mli_Color out;
        out.r = (1.f - refl) * a.r + refl * b.r;
        out.g = (1.f - refl) * a.g + refl * b.g;
        out.b = (1.f - refl) * a.b + refl * b.b;
        return out;
}

struct mli_Color mli_Color_mean(
        const struct mli_Color colors[],
        const uint32_t num_colors)
{
        struct mli_Color out = {0., 0., 0.};
        const float f_num_colors = (float)num_colors;
        uint32_t i;
        for (i = 0; i < num_colors; i++) {
                out.r = out.r + colors[i].r;
                out.g = out.g + colors[i].g;
                out.b = out.b + colors[i].b;
        }
        out.r = out.r / f_num_colors;
        out.g = out.g / f_num_colors;
        out.b = out.b / f_num_colors;
        return out;
}

struct mli_Color mli_Color_truncate(
        const struct mli_Color color,
        const float start,
        const float stop)
{
        struct mli_Color out;
        out.r = color.r;
        out.g = color.g;
        out.b = color.b;
        if (out.r > stop)
                out.r = stop;
        if (out.r < start)
                out.r = start;
        if (out.g > stop)
                out.g = stop;
        if (out.g < start)
                out.g = start;
        if (out.b > stop)
                out.b = stop;
        if (out.b < start)
                out.b = start;
        return out;
}

int mli_Color_equal(const struct mli_Color a, const struct mli_Color b)
{
        if (a.r != b.r)
                return 0;
        if (a.g != b.g)
                return 0;
        if (a.b != b.b)
                return 0;
        return 1;
}

int mli_Color_is_in_range(
        const struct mli_Color c,
        const float start,
        const float stop)
{
        if (MLI_MATH_IS_NAN(c.r))
                return 0;
        if (MLI_MATH_IS_NAN(c.g))
                return 0;
        if (MLI_MATH_IS_NAN(c.b))
                return 0;

        if (c.r < start || c.r >= stop)
                return 0;
        if (c.g < start || c.g >= stop)
                return 0;
        if (c.b < start || c.b >= stop)
                return 0;
        return 1;
}

float mli_Color_luminance(const struct mli_Color self)
{
        return self.r + self.g + self.b;
}

struct mli_Color mli_Color_add(
        const struct mli_Color u,
        const struct mli_Color v)
{
        return mli_Color_set(u.r + v.r, u.g + v.g, u.b + v.b);
}

struct mli_Color mli_Color_multiply(const struct mli_Color c, const double f)
{
        return mli_Color_set(c.r * f, c.g * f, c.b * f);
}

struct mli_Color mli_Color_multiply_elementwise(
        const struct mli_Color u,
        const struct mli_Color v)
{
        return mli_Color_set(u.r * v.r, u.g * v.g, u.b * v.b);
}

/* color_cie1931 */
/* ------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

int mli_cie1931_spectral_matching_curve_x(struct mli_Func *self)
{
        float X[][2] = {
                {2.000e-07, 0.000e+00}, {3.800e-07, 0.000e+00},
                {3.850e-07, 2.957e-03}, {3.900e-07, 5.914e-03},
                {3.950e-07, 8.871e-03}, {4.000e-07, 1.454e-02},
                {4.050e-07, 2.021e-02}, {4.100e-07, 3.462e-02},
                {4.150e-07, 5.486e-02}, {4.200e-07, 8.703e-02},
                {4.250e-07, 1.271e-01}, {4.300e-07, 1.582e-01},
                {4.350e-07, 1.844e-01}, {4.400e-07, 1.907e-01},
                {4.450e-07, 1.970e-01}, {4.500e-07, 1.875e-01},
                {4.550e-07, 1.780e-01}, {4.600e-07, 1.579e-01},
                {4.650e-07, 1.352e-01}, {4.700e-07, 1.085e-01},
                {4.750e-07, 8.082e-02}, {4.800e-07, 5.804e-02},
                {4.850e-07, 3.853e-02}, {4.900e-07, 2.404e-02},
                {4.950e-07, 1.708e-02}, {5.000e-07, 1.178e-02},
                {5.050e-07, 1.315e-02}, {5.100e-07, 1.451e-02},
                {5.150e-07, 3.042e-02}, {5.200e-07, 4.633e-02},
                {5.250e-07, 7.357e-02}, {5.300e-07, 1.036e-01},
                {5.350e-07, 1.397e-01}, {5.400e-07, 1.773e-01},
                {5.450e-07, 2.149e-01}, {5.500e-07, 2.590e-01},
                {5.550e-07, 3.076e-01}, {5.600e-07, 3.561e-01},
                {5.650e-07, 4.025e-01}, {5.700e-07, 4.489e-01},
                {5.750e-07, 4.937e-01}, {5.800e-07, 5.324e-01},
                {5.850e-07, 5.622e-01}, {5.900e-07, 5.859e-01},
                {5.950e-07, 5.949e-01}, {6.000e-07, 5.940e-01},
                {6.050e-07, 5.829e-01}, {6.100e-07, 5.566e-01},
                {6.150e-07, 5.216e-01}, {6.200e-07, 4.738e-01},
                {6.250e-07, 4.116e-01}, {6.300e-07, 3.494e-01},
                {6.350e-07, 2.991e-01}, {6.400e-07, 2.488e-01},
                {6.450e-07, 2.023e-01}, {6.500e-07, 1.608e-01},
                {6.550e-07, 1.290e-01}, {6.600e-07, 9.935e-02},
                {6.650e-07, 7.821e-02}, {6.700e-07, 5.707e-02},
                {6.750e-07, 4.652e-02}, {6.800e-07, 3.596e-02},
                {6.850e-07, 2.932e-02}, {6.900e-07, 2.365e-02},
                {6.950e-07, 2.069e-02}, {7.000e-07, 1.955e-02},
                {7.050e-07, 1.841e-02}, {7.100e-07, 1.727e-02},
                {7.150e-07, 1.625e-02}, {7.200e-07, 1.575e-02},
                {7.250e-07, 1.525e-02}, {7.300e-07, 1.475e-02},
                {7.350e-07, 1.425e-02}, {7.400e-07, 1.375e-02},
                {7.450e-07, 1.325e-02}, {7.500e-07, 1.243e-02},
                {7.550e-07, 1.036e-02}, {7.600e-07, 8.290e-03},
                {7.650e-07, 6.217e-03}, {7.700e-07, 4.145e-03},
                {7.750e-07, 2.072e-03}, {7.800e-07, 0.000e+00},
                {1.200e-06, 0.000e+00},
        };
        uint64_t i;
        chk_msg(mli_Func_malloc(self, sizeof(X) / sizeof(X[0])),
                "Can't malloc function from cie1931.X.");
        for (i = 0; i < self->num_points; i++) {
                self->x[i] = X[i][0];
                self->y[i] = X[i][1];
        }
        return 1;
chk_error:
        return 0;
}

int mli_cie1931_spectral_matching_curve_y(struct mli_Func *self)
{
        float Y[][2] = {
                {2.000e-07, 0.000e+00}, {3.800e-07, 0.000e+00},
                {3.850e-07, 1.099e-03}, {3.900e-07, 2.198e-03},
                {3.950e-07, 3.296e-03}, {4.000e-07, 4.395e-03},
                {4.050e-07, 5.494e-03}, {4.100e-07, 6.593e-03},
                {4.150e-07, 7.691e-03}, {4.200e-07, 8.676e-03},
                {4.250e-07, 9.633e-03}, {4.300e-07, 1.155e-02},
                {4.350e-07, 1.491e-02}, {4.400e-07, 1.827e-02},
                {4.450e-07, 2.162e-02}, {4.500e-07, 2.498e-02},
                {4.550e-07, 3.303e-02}, {4.600e-07, 4.109e-02},
                {4.650e-07, 4.914e-02}, {4.700e-07, 5.719e-02},
                {4.750e-07, 6.991e-02}, {4.800e-07, 8.575e-02},
                {4.850e-07, 1.047e-01}, {4.900e-07, 1.282e-01},
                {4.950e-07, 1.574e-01}, {5.000e-07, 1.950e-01},
                {5.050e-07, 2.492e-01}, {5.100e-07, 3.075e-01},
                {5.150e-07, 3.631e-01}, {5.200e-07, 4.139e-01},
                {5.250e-07, 4.602e-01}, {5.300e-07, 4.948e-01},
                {5.350e-07, 5.257e-01}, {5.400e-07, 5.422e-01},
                {5.450e-07, 5.588e-01}, {5.500e-07, 5.625e-01},
                {5.550e-07, 5.663e-01}, {5.600e-07, 5.589e-01},
                {5.650e-07, 5.487e-01}, {5.700e-07, 5.311e-01},
                {5.750e-07, 5.086e-01}, {5.800e-07, 4.825e-01},
                {5.850e-07, 4.509e-01}, {5.900e-07, 4.193e-01},
                {5.950e-07, 3.876e-01}, {6.000e-07, 3.543e-01},
                {6.050e-07, 3.198e-01}, {6.100e-07, 2.853e-01},
                {6.150e-07, 2.509e-01}, {6.200e-07, 2.195e-01},
                {6.250e-07, 1.880e-01}, {6.300e-07, 1.566e-01},
                {6.350e-07, 1.272e-01}, {6.400e-07, 1.055e-01},
                {6.450e-07, 8.373e-02}, {6.500e-07, 6.876e-02},
                {6.550e-07, 5.379e-02}, {6.600e-07, 4.365e-02},
                {6.650e-07, 3.472e-02}, {6.700e-07, 2.878e-02},
                {6.750e-07, 2.483e-02}, {6.800e-07, 2.088e-02},
                {6.850e-07, 1.693e-02}, {6.900e-07, 1.582e-02},
                {6.950e-07, 1.542e-02}, {7.000e-07, 1.503e-02},
                {7.050e-07, 1.463e-02}, {7.100e-07, 1.424e-02},
                {7.150e-07, 1.384e-02}, {7.200e-07, 1.345e-02},
                {7.250e-07, 1.305e-02}, {7.300e-07, 1.266e-02},
                {7.350e-07, 1.226e-02}, {7.400e-07, 1.186e-02},
                {7.450e-07, 1.147e-02}, {7.500e-07, 1.079e-02},
                {7.550e-07, 8.994e-03}, {7.600e-07, 7.195e-03},
                {7.650e-07, 5.396e-03}, {7.700e-07, 3.598e-03},
                {7.750e-07, 1.799e-03}, {7.800e-07, 0.000e+00},
                {1.200e-06, 0.000e+00},
        };
        uint64_t i;
        chk_msg(mli_Func_malloc(self, sizeof(Y) / sizeof(Y[0])),
                "Can't malloc function from cie1931.Y.");
        for (i = 0; i < self->num_points; i++) {
                self->x[i] = Y[i][0];
                self->y[i] = Y[i][1];
        }
        return 1;
chk_error:
        return 0;
}

int mli_cie1931_spectral_matching_curve_z(struct mli_Func *self)
{
        float Z[][2] = {
                {2.000e-07, 0.000e+00}, {3.800e-07, 0.000e+00},
                {3.850e-07, 9.873e-03}, {3.900e-07, 1.595e-02},
                {3.950e-07, 2.536e-02}, {4.000e-07, 4.352e-02},
                {4.050e-07, 7.715e-02}, {4.100e-07, 1.473e-01},
                {4.150e-07, 2.240e-01}, {4.200e-07, 4.372e-01},
                {4.250e-07, 6.496e-01}, {4.300e-07, 8.370e-01},
                {4.350e-07, 9.538e-01}, {4.400e-07, 9.983e-01},
                {4.450e-07, 1.000e+00}, {4.500e-07, 9.912e-01},
                {4.550e-07, 9.622e-01}, {4.600e-07, 9.131e-01},
                {4.650e-07, 7.944e-01}, {4.700e-07, 6.640e-01},
                {4.750e-07, 5.457e-01}, {4.800e-07, 4.268e-01},
                {4.850e-07, 3.338e-01}, {4.900e-07, 2.524e-01},
                {4.950e-07, 1.944e-01}, {5.000e-07, 1.536e-01},
                {5.050e-07, 1.215e-01}, {5.100e-07, 9.421e-02},
                {5.150e-07, 6.725e-02}, {5.200e-07, 4.686e-02},
                {5.250e-07, 3.640e-02}, {5.300e-07, 2.902e-02},
                {5.350e-07, 2.242e-02}, {5.400e-07, 1.802e-02},
                {5.450e-07, 1.509e-02}, {5.500e-07, 1.329e-02},
                {5.550e-07, 1.318e-02}, {5.600e-07, 1.307e-02},
                {5.650e-07, 1.296e-02}, {5.700e-07, 1.284e-02},
                {5.750e-07, 1.273e-02}, {5.800e-07, 1.262e-02},
                {5.850e-07, 1.260e-02}, {5.900e-07, 1.261e-02},
                {5.950e-07, 1.262e-02}, {6.000e-07, 1.263e-02},
                {6.050e-07, 1.264e-02}, {6.100e-07, 1.264e-02},
                {6.150e-07, 1.265e-02}, {6.200e-07, 1.266e-02},
                {6.250e-07, 1.267e-02}, {6.300e-07, 1.267e-02},
                {6.350e-07, 1.268e-02}, {6.400e-07, 1.269e-02},
                {6.450e-07, 1.270e-02}, {6.500e-07, 1.270e-02},
                {6.550e-07, 1.271e-02}, {6.600e-07, 1.272e-02},
                {6.650e-07, 1.273e-02}, {6.700e-07, 1.274e-02},
                {6.750e-07, 1.274e-02}, {6.800e-07, 1.275e-02},
                {6.850e-07, 1.276e-02}, {6.900e-07, 1.277e-02},
                {6.950e-07, 1.277e-02}, {7.000e-07, 1.278e-02},
                {7.050e-07, 1.279e-02}, {7.100e-07, 1.280e-02},
                {7.150e-07, 1.281e-02}, {7.200e-07, 1.281e-02},
                {7.250e-07, 1.282e-02}, {7.300e-07, 1.283e-02},
                {7.350e-07, 1.284e-02}, {7.400e-07, 1.284e-02},
                {7.450e-07, 1.285e-02}, {7.500e-07, 1.244e-02},
                {7.550e-07, 1.037e-02}, {7.600e-07, 8.295e-03},
                {7.650e-07, 6.222e-03}, {7.700e-07, 4.148e-03},
                {7.750e-07, 2.074e-03}, {7.800e-07, 0.000e+00},
                {1.200e-06, 0.000e+00},
        };
        uint64_t i;
        chk_msg(mli_Func_malloc(self, sizeof(Z) / sizeof(Z[0])),
                "Can't malloc function from cie1931.Z.");
        for (i = 0; i < self->num_points; i++) {
                self->x[i] = Z[i][0];
                self->y[i] = Z[i][1];
        }
        return 1;
chk_error:
        return 0;
}

struct mli_Mat mli_cie1931_spectral_matching_xyz_to_rgb(void)
{
        struct mli_Mat m;
        m.r00 = +3.240479;
        m.r10 = -1.537150;
        m.r20 = -0.498535;
        m.r01 = -0.969256;
        m.r11 = +1.875991;
        m.r21 = +0.041556;
        m.r02 = +0.055648;
        m.r12 = -0.204043;
        m.r22 = +1.057311;
        return m;
}

struct mli_Mat mli_cie1931_spectral_matching_rgb_to_xyz(void)
{
        struct mli_Mat m;
        m.r00 = +0.412453;
        m.r10 = +0.357580;
        m.r20 = +0.180423;
        m.r01 = +0.212671;
        m.r11 = +0.715160;
        m.r21 = +0.072169;
        m.r02 = +0.019334;
        m.r12 = +0.119193;
        m.r22 = +0.950227;
        return m;
}

int mli_cie1931_spectral_radiance_of_black_body_W_per_m2_per_sr_per_m(
        struct mli_Func *self,
        const double wavelength_start,
        const double wavelength_stop,
        const double temperature,
        const uint64_t num_points)
{
        uint64_t i;
        const double start = wavelength_start;
        const double stop = wavelength_stop;
        const double step = (stop - start) / (double)num_points;
        chk_msg(mli_Func_malloc(self, num_points),
                "Can't malloc function for black body spectrum.");
        for (i = 0; i < self->num_points; i++) {
                const double wavelength = start + step * (double)i;
                const double radiance =
                        mli_physics_plancks_spectral_radiance_law_W_per_m2_per_sr_per_m(
                                wavelength, temperature);
                self->x[i] = wavelength;
                self->y[i] = radiance;
        }
        return 1;
chk_error:
        return 0;
}

/* color_materials */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_ColorMaterials mli_ColorMaterials_init(void)
{
        struct mli_ColorMaterials out;
        out.wavelength_bin_edges = mli_ColorSpectrumBinEdges_init();

        out.observer_matching_curve_x = mli_ColorSpectrum_init_zeros();
        out.observer_matching_curve_y = mli_ColorSpectrum_init_zeros();
        out.observer_matching_curve_z = mli_ColorSpectrum_init_zeros();

        out.observer_matching_curve_xyz_to_rgb = mli_Mat_unity();

        out.spectra = mli_ColorSpectrumArray_init();
        return out;
}

void mli_ColorMaterials_free(struct mli_ColorMaterials *self)
{
        mli_ColorSpectrumArray_free(&self->spectra);
        (*self) = mli_ColorMaterials_init();
}

int mli_ColorMaterials_set_observer_cie1931(struct mli_ColorMaterials *self)
{
        struct mli_Func func = mli_Func_init();

        self->observer_matching_curve_xyz_to_rgb =
                mli_cie1931_spectral_matching_xyz_to_rgb();

        chk(mli_cie1931_spectral_matching_curve_x(&func));
        chk_msg(mli_ColorSpectrum_from_func(
                        &self->observer_matching_curve_x,
                        &self->wavelength_bin_edges,
                        &func),
                "Failed to resample cie1931.x onto ColorSpectrum.");

        chk(mli_cie1931_spectral_matching_curve_y(&func));
        chk_msg(mli_ColorSpectrum_from_func(
                        &self->observer_matching_curve_y,
                        &self->wavelength_bin_edges,
                        &func),
                "Failed to resample cie1931.y onto ColorSpectrum.");

        chk(mli_cie1931_spectral_matching_curve_z(&func));
        chk_msg(mli_ColorSpectrum_from_func(
                        &self->observer_matching_curve_z,
                        &self->wavelength_bin_edges,
                        &func),
                "Failed to resample cie1931.z onto ColorSpectrum.");

        mli_Func_free(&func);
        return 1;
chk_error:
        return 0;
}

int mli_ColorMaterials_malloc(
        struct mli_ColorMaterials *self,
        const uint64_t num_spectra)
{
        mli_ColorMaterials_free(self);
        chk_msg(mli_ColorSpectrumArray_malloc(&self->spectra, num_spectra),
                "Failed to malloc spectra");

        chk_msg(mli_ColorMaterials_set_observer_cie1931(self),
                "Failed to set observer_matching_curves.");
        return 1;
chk_error:
        mli_ColorMaterials_free(self);
        return 0;
}

int mli_ColorMaterials_malloc_from_Materials(
        struct mli_ColorMaterials *self,
        const struct mli_Materials *materials)
{
        uint64_t i;
        chk_msg(mli_ColorMaterials_malloc(self, materials->spectra.size),
                "Can't malloc ColorMaterials from Materials.");

        for (i = 0; i < materials->spectra.size; i++) {
                chk_msg(mli_ColorSpectrum_from_func(
                                &self->spectra.array[i],
                                &self->wavelength_bin_edges,
                                &materials->spectra.array[i].spectrum),
                        "Failed to resample spectrum onto ColorSpectrum.");
        }

        return 1;
chk_error:
        mli_ColorMaterials_free(self);
        return 0;
}

struct mli_Vec mli_ColorMaterials_ColorSpectrum_to_xyz(
        const struct mli_ColorMaterials *self,
        const struct mli_ColorSpectrum *spectrum)
{
        struct mli_Vec observer;
        observer.x = mli_ColorSpectrum_multiply_and_sum(
                &self->observer_matching_curve_x, spectrum);
        observer.y = mli_ColorSpectrum_multiply_and_sum(
                &self->observer_matching_curve_y, spectrum);
        observer.z = mli_ColorSpectrum_multiply_and_sum(
                &self->observer_matching_curve_z, spectrum);
        return observer;
}

/* color_spectrum */
/* -------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

struct mli_ColorSpectrumBinEdges mli_ColorSpectrumBinEdges_init(void)
{
        uint64_t i;
        const uint64_t size = MLI_COLORSPECTRUM_SIZE;
        const float start = MLI_COLORSPECTRUM_WAVELENGTH_START;
        const float stop = MLI_COLORSPECTRUM_WAVELENGTH_STOP;
        const float step = (stop - start) / (float)size;
        struct mli_ColorSpectrumBinEdges edges;
        for (i = 0; i < size + 1; i++) {
                edges.values[i] = start + (float)i * step;
        }
        return edges;
}

struct mli_ColorSpectrum mli_ColorSpectrum_init_zeros(void)
{
        struct mli_ColorSpectrum cs;
        MLI_MATH_ARRAY_SET(cs.values, 0.0, MLI_COLORSPECTRUM_SIZE);
        return cs;
}

void mli_ColorSpectrum_set(struct mli_ColorSpectrum *self, const float value)
{
        MLI_MATH_ARRAY_SET(self->values, value, MLI_COLORSPECTRUM_SIZE);
}

struct mli_ColorSpectrum mli_ColorSpectrum_add(
        const struct mli_ColorSpectrum a,
        const struct mli_ColorSpectrum b)
{
        uint64_t i;
        struct mli_ColorSpectrum c;
        for (i = 0; i < MLI_COLORSPECTRUM_SIZE; i++) {
                c.values[i] = a.values[i] + b.values[i];
        }
        return c;
}

struct mli_ColorSpectrum mli_ColorSpectrum_multiply(
        const struct mli_ColorSpectrum a,
        const struct mli_ColorSpectrum b)
{
        uint64_t i;
        struct mli_ColorSpectrum c;
        for (i = 0; i < MLI_COLORSPECTRUM_SIZE; i++) {
                c.values[i] = a.values[i] * b.values[i];
        }
        return c;
}

struct mli_ColorSpectrum mli_ColorSpectrum_multiply_scalar(
        const struct mli_ColorSpectrum a,
        const double factor)
{
        uint64_t i;
        struct mli_ColorSpectrum c;
        for (i = 0; i < MLI_COLORSPECTRUM_SIZE; i++) {
                c.values[i] = a.values[i] * factor;
        }
        return c;
}

struct mli_ColorSpectrum mli_ColorSpectrum_exp(
        const struct mli_ColorSpectrum a,
        const double factor)
{
        uint64_t i;
        struct mli_ColorSpectrum c;
        for (i = 0; i < MLI_COLORSPECTRUM_SIZE; i++) {
                c.values[i] = exp(a.values[i] * factor);
        }
        return c;
}

double mli_ColorSpectrum_multiply_and_sum(
        const struct mli_ColorSpectrum *a,
        const struct mli_ColorSpectrum *b)
{
        uint64_t i;
        double value = 0.0;
        for (i = 0; i < MLI_COLORSPECTRUM_SIZE; i++) {
                value += a->values[i] * b->values[i];
        }
        return value;
}

float mli_ColorSpectrum_sum(const struct mli_ColorSpectrum *self)
{
        uint64_t i;
        float c = 0.0;
        for (i = 0; i < MLI_COLORSPECTRUM_SIZE; i++) {
                c += self->values[i];
        }
        return c;
}

int mli_ColorSpectrum_from_func(
        struct mli_ColorSpectrum *self,
        const struct mli_ColorSpectrumBinEdges *wavelength_bin_edges,
        const struct mli_Func *func)
{
        const uint64_t NUM_SAMPLES = 100;
        uint64_t i, u;

        for (i = 0; i < MLI_COLORSPECTRUM_SIZE; i++) {
                const double start = wavelength_bin_edges->values[i];
                const double stop = wavelength_bin_edges->values[i + 1];
                const double step = (stop - start) / (double)NUM_SAMPLES;

                double val_mean = 0.0;
                for (u = 0; u < NUM_SAMPLES; u++) {
                        const double wavelength = start + (float)u * step;
                        double val = 0.0;
                        chk_msg(mli_Func_evaluate(func, wavelength, &val),
                                "Can't evaluate spectrum.");
                        val_mean += val;
                }
                val_mean = val_mean / (double)NUM_SAMPLES;
                self->values[i] = val_mean;
        }

        return 1;
chk_error:
        return 0;
}

int mli_ColorSpectrum_set_radiance_of_black_body_W_per_m2_per_sr(
        struct mli_ColorSpectrum *self,
        const double temperature)
{
        struct mli_ColorSpectrumBinEdges edges =
                mli_ColorSpectrumBinEdges_init();
        const uint64_t NUM_SAMPLES = 100;
        uint64_t i, u;

        for (i = 0; i < MLI_COLORSPECTRUM_SIZE; i++) {
                const double start = edges.values[i];
                const double stop = edges.values[i + 1];
                const double step = (stop - start) / (double)NUM_SAMPLES;

                double average_radiance = 0.0;

                for (u = 0; u < NUM_SAMPLES; u++) {
                        const double wavelength = start + (float)u * step;
                        const double radiance =
                                mli_physics_plancks_spectral_radiance_law_W_per_m2_per_sr_per_m(
                                        wavelength, temperature);
                        average_radiance += radiance;
                }
                average_radiance /= (double)NUM_SAMPLES;

                self->values[i] = average_radiance * step;
        }
        return 1;
}

/* color_spectrum_array */
/* -------------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

MLI_ARRAY_IMPLEMENTATION(mli_ColorSpectrumArray, struct mli_ColorSpectrum)

/* color_vector */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

MLI_VECTOR_IMPLEMENTATION(mli_ColorVector, struct mli_Color)

/* cstr */
/* ---- */

/* Copyright Sebastian Achim Mueller */

int mli_cstr_ends_with(const char *str, const char *sufix)
{
        uint64_t len_str, len_sufix;
        if (!str || !sufix) {
                return 0;
        }
        len_str = strlen(str);
        len_sufix = strlen(sufix);
        if (len_sufix > len_str) {
                return 0;
        }
        return strncmp(str + len_str - len_sufix, sufix, len_sufix) == 0;
}

int mli_cstr_starts_with(const char *str, const char *prefix)
{
        uint64_t len_str, len_prefix;
        if (!str || !prefix) {
                return 0;
        }
        len_str = strlen(str);
        len_prefix = strlen(prefix);
        if (len_prefix > len_str) {
                return 0;
        }
        return strncmp(str, prefix, len_prefix) == 0;
}

int mli_cstr_is_CRLF(const char *s)
{
        if (s[0] == '\0') {
                return 0;
        }
        if (s[1] == '\0') {
                return 0;
        }
        if (s[0] == '\r' && s[1] == '\n') {
                return 1;
        }
        return 0;
}

int mli_cstr_is_CR(const char *s)
{
        if (s[0] == '\0') {
                return 0;
        }
        if (s[0] == '\r') {
                return 1;
        }
        return 0;
}

int mli_cstr_assert_only_NUL_LF_TAB_controls(const char *str)
{
        return mli_cstr_assert_only_NUL_LF_TAB_controls_dbg(str, 1);
}

int mli_cstr_assert_only_NUL_LF_TAB_controls_dbg(const char *str, const int dbg)
{
        uint64_t pos = 0;
        while (str[pos] != '\0') {
                if (str[pos] >= 32 && str[pos] < 127) {
                        /* all fine */
                } else {
                        if (str[pos] == '\n') {
                                /* fine */
                        } else if (str[pos] == '\t') {
                                /* fine */
                        } else {
                                if (dbg) {
                                        fprintf(stderr,
                                                "Control code %u "
                                                "at column %ld in string.\n",
                                                (uint8_t)str[pos],
                                                pos);
                                }
                                return 0;
                        }
                }
                pos += 1;
        }
        return 1;
}

int mli_fprint_line_match(
        FILE *f,
        const int64_t line,
        const int64_t line_number)
{
        chk(fprintf(f, "% 6d", (int32_t)line));
        if (line == line_number) {
                chk(fprintf(f, "->|  "));
        } else {
                chk(fprintf(f, "  |  "));
        }
        return 1;
chk_error:
        return 0;
}

int mli_cstr_match_templeate(
        const char *s,
        const char *t,
        const char digit_wildcard)
{
        uint64_t i;
        if (strlen(s) != strlen(t)) {
                return 0;
        }
        for (i = 0; i < strlen(s); i++) {
                if (t[i] == digit_wildcard) {
                        if (!isdigit(s[i])) {
                                return 0;
                        }
                } else {
                        if (s[i] != t[i]) {
                                return 0;
                        }
                }
        }
        return 1;
}

/* cstr_numbers */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_cstr_nto_int64(
        int64_t *out,
        const char *s,
        const uint64_t base,
        const uint64_t expected_num_chars)
{
        char *end;
        uint64_t actual_num_chars = 0u;
        int64_t l;
        chk_msg(!(s[0] == '\0' || isspace(s[0])),
                "Can not convert string to int64, bad string.");
        errno = 0;
        l = strtol(s, &end, base);
        chk_msg(errno != ERANGE,
                "Can not convert string to int64, over-, under-flow.");
        chk_msg(end != NULL, "Can not convert string to int64, bad string.");
        actual_num_chars = end - s;
        chk_msg(actual_num_chars == expected_num_chars,
                "Integer has not the expected number of chars.");
        *out = l;
        return 1;
chk_error:
        return 0;
}

int mli_cstr_to_int64(int64_t *out, const char *s, const uint64_t base)
{
        chk_msg(mli_cstr_nto_int64(out, s, base, strlen(s)),
                "Can not convert string to int64.");
        return 1;
chk_error:
        return 0;
}

int mli_cstr_nto_uint64(
        uint64_t *out,
        const char *s,
        const uint64_t base,
        const uint64_t expected_num_chars)
{
        int64_t tmp;
        chk(mli_cstr_nto_int64(&tmp, s, base, expected_num_chars));
        chk_msg(tmp >= 0, "Expected a positive integer.");
        (*out) = tmp;
        return 1;
chk_error:
        return 0;
}

int mli_cstr_to_uint64(uint64_t *out, const char *s, const uint64_t base)
{
        int64_t tmp;
        chk(mli_cstr_to_int64(&tmp, s, base));
        chk_msg(tmp >= 0, "Expected a positive integer.");
        (*out) = tmp;
        return 1;
chk_error:
        return 0;
}

int mli_cstr_nto_double(
        double *out,
        const char *s,
        const uint64_t expected_num_chars)
{
        char *end;
        uint64_t actual_num_chars = 0u;
        double l;
        chk_msg(!(s[0] == '\0' || isspace(s[0])),
                "Can not convert string to float64, bad string.");
        errno = 0;
        l = strtod(s, &end);
        chk_msg(errno != ERANGE,
                "Can not convert string to float64, over-, under-flow.");
        chk_msg(end != NULL, "Can not convert string to float64.");

        actual_num_chars = end - s;
        chk_msg(actual_num_chars == expected_num_chars,
                "float64 has not the expected number of chars.");
        *out = l;
        return 1;
chk_error:
        return 0;
}

int mli_cstr_to_double(double *out, const char *s)
{
        chk_msg(mli_cstr_nto_double(out, s, strlen(s)),
                "Can not convert string to float64.");
        return 1;
chk_error:
        return 0;
}

int mli_cstr_print_uint64(
        uint64_t u,
        char *s,
        const uint64_t max_num_chars,
        const uint64_t base,
        const uint64_t min_num_digits)
{
        char literals[] = {
                '0',
                '1',
                '2',
                '3',
                '4',
                '5',
                '6',
                '7',
                '8',
                '9',
                'A',
                'B',
                'C',
                'D',
                'E',
                'F'};
        char tmp[128] = {'\0'};
        uint64_t remainder = 0u;
        uint32_t remainder32 = 0u;
        uint64_t quotient = u;
        int64_t digs = 0u;
        int64_t pos = 0;
        int64_t i = 0;
        int64_t num_leading_zeors = 0;

        chk_msg(base <= 16, "Expected base <= 16");
        chk_msg(base > 1, "Expected base > 1");
        chk_msg(max_num_chars < sizeof(tmp), "Exceeded max num. chars.");
        chk_msg(min_num_digits < max_num_chars, "Exceeded max num. chars.");

        do {
                remainder = quotient % base;
                quotient = quotient / base;
                remainder32 = (uint32_t)remainder;
                tmp[digs] = literals[remainder32];
                digs++;
                chk_msg(digs < (int64_t)sizeof(tmp),
                        "Exceeded max num. chars.");
        } while (quotient > 0u);

        num_leading_zeors = min_num_digits - digs;
        if (num_leading_zeors < 0) {
                num_leading_zeors = 0;
        }

        for (i = 0; i < num_leading_zeors; i++) {
                chk_msg(pos < (int64_t)max_num_chars,
                        "Exceeded max num. chars.");
                s[pos] = '0';
                pos++;
        }

        for (i = 0; i < digs; i++) {
                chk_msg(pos < (int64_t)max_num_chars,
                        "Exceeded max num. chars.");
                s[pos] = tmp[digs - i - 1];
                pos++;
        }

        chk_msg(pos < (int64_t)max_num_chars, "Exceeded max num. chars.");
        s[pos] = '\0';

        return 1;
chk_error:
        return 0;
}

/* cube */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Vec mli_Cube_upper(const struct mli_Cube a)
{
        return mli_Vec_add(
                a.lower,
                mli_Vec_init(a.edge_length, a.edge_length, a.edge_length));
}

struct mli_AABB mli_Cube_to_aabb(const struct mli_Cube a)
{
        struct mli_AABB out;
        out.lower = a.lower;
        out.upper = mli_Cube_upper(a);
        return out;
}

struct mli_Vec mli_Cube_center(const struct mli_Cube a)
{
        return mli_Vec_init(
                a.lower.x + a.edge_length * .5,
                a.lower.y + a.edge_length * .5,
                a.lower.z + a.edge_length * .5);
}

struct mli_Cube mli_Cube_outermost_cube(const struct mli_AABB a)
{
        struct mli_Cube cube;
        struct mli_Vec center;
        struct mli_Vec half_diagonal;
        struct mli_Vec diff;
        double max_half_length;
        center = mli_AABB_center(a);
        diff = mli_Vec_substract(a.upper, a.lower);
        max_half_length = .5 * MLI_MATH_MAX3(diff.x, diff.y, diff.z);
        half_diagonal.x = max_half_length;
        half_diagonal.y = max_half_length;
        half_diagonal.z = max_half_length;
        cube.lower = mli_Vec_substract(center, half_diagonal);
        cube.edge_length = max_half_length * 2.;
        return cube;
}

struct mli_Cube mli_Cube_octree_child(
        const struct mli_Cube cube,
        const uint32_t sx,
        const uint32_t sy,
        const uint32_t sz)
{
        struct mli_Cube child;
        struct mli_Vec length;
        struct mli_Vec center = mli_Cube_center(cube);
        length = mli_Vec_substract(center, cube.lower);
        child.lower = cube.lower;
        child.edge_length = .5 * cube.edge_length;
        if (sx) {
                child.lower.x += length.x;
        }
        if (sy) {
                child.lower.y += length.y;
        }
        if (sz) {
                child.lower.z += length.z;
        }
        return child;
}

struct mli_Cube mli_Cube_octree_child_code(
        const struct mli_Cube cube,
        const uint8_t a)
{
        struct mli_Cube child;
        struct mli_Vec length;
        struct mli_Vec center = mli_Cube_center(cube);
        length = mli_Vec_substract(center, cube.lower);
        child.lower = cube.lower;
        child.edge_length = .5 * cube.edge_length;
        if (MLI_MATH_IS_BIT(a, 2)) {
                child.lower.x += length.x;
        }
        if (MLI_MATH_IS_BIT(a, 1)) {
                child.lower.y += length.y;
        }
        if (MLI_MATH_IS_BIT(a, 0)) {
                child.lower.z += length.z;
        }
        return child;
}

int mli_Cube_equal(const struct mli_Cube a, const struct mli_Cube b)
{
        if (a.edge_length != b.edge_length)
                return 0;
        if (!mli_Vec_equal(a.lower, b.lower))
                return 0;
        return 1;
}

/* double_array */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_ARRAY_IMPLEMENTATION(mli_DoubleArray, double)

/* double_vector */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_VECTOR_IMPLEMENTATION(mli_DoubleVector, double)

/* float_array */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_ARRAY_IMPLEMENTATION(mli_FloatArray, float)

/* float_vector */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_VECTOR_IMPLEMENTATION(mli_FloatVector, float)

/* frame */
/* ----- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Frame;

struct mli_Frame mli_Frame_init(void)
{
        struct mli_Frame f;
        f.type = MLI_FRAME_TYPE_FRAME;
        f.id = 0u;
        f.frame2mother.translation = mli_Vec_init(0., 0., 0.);
        f.frame2mother.rotation = mli_Quaternion_set_tait_bryan(0., 0., 0.);
        f.frame2root = f.frame2mother;
        f.mother = NULL;
        f.children = mli_FramePtrVector_init();
        f.object = 0u;

        f.boundary_layers = mli_Uint32Vector_init();
        return f;
}

void mli_Frame_free(struct mli_Frame *f)
{
        uint64_t c;
        if (f->type == MLI_FRAME_TYPE_FRAME) {
                for (c = 0; c < f->children.size; c++) {
                        struct mli_Frame *child = f->children.array[c];
                        mli_Frame_free(child);
                }
                mli_FramePtrVector_free(&f->children);
        }
        if (f->type == MLI_FRAME_TYPE_OBJECT) {
                mli_Uint32Vector_free(&f->boundary_layers);
        }
        (*f) = mli_Frame_init();
}

int mli_Frame_malloc(struct mli_Frame *f, const uint64_t type)
{
        mli_Frame_free(f);
        f->type = type;
        if (type == MLI_FRAME_TYPE_FRAME) {
                chk_msg(mli_FramePtrVector_malloc(&f->children, 0u),
                        "Can not allocate children of frame.");
        }
        if (type == MLI_FRAME_TYPE_OBJECT) {
                chk_msg(mli_Uint32Vector_malloc(&f->boundary_layers, 0u),
                        "Failed to malloc frame's boundary_layers.");
        }
        return 1;
chk_error:
        return 0;
}

int mli_Frame_set_mother_and_child(
        struct mli_Frame *mother,
        struct mli_Frame *child)
{
        chk_msg(mother->type == MLI_FRAME_TYPE_FRAME,
                "Expected mother to be of type FRAME");
        chk_msg(mli_FramePtrVector_push_back(&mother->children, child),
                "Can not push back child-frame.");

        child->mother = (struct mli_Frame *)mother;
        return 1;
chk_error:
        return 0;
}

struct mli_Frame *mli_Frame_add(struct mli_Frame *mother, const uint64_t type)
{
        struct mli_Frame *child = NULL;
        chk_malloc(child, struct mli_Frame, 1u);
        chk_msg(mli_Frame_malloc(child, type), "Can not allocate child-frame.");
        chk_msg(mli_Frame_set_mother_and_child(mother, child),
                "Can not allocate child-pointer.");
        return child;
chk_error:
        return NULL;
}

int mli_frame_type_to_string(const uint64_t type, char *s)
{
        switch (type) {
        case MLI_FRAME_TYPE_FRAME:
                sprintf(s, "frame");
                break;
        case MLI_FRAME_TYPE_OBJECT:
                sprintf(s, "object");
                break;
        default:
                chk_bad("Type is unknown.");
                break;
        }
        return 1;
chk_error:
        return 0;
}

int mli_frame_string_to_type(const char *s, uint64_t *type)
{
        if (strcmp(s, "frame") == 0) {
                *type = MLI_FRAME_TYPE_FRAME;
        } else if (strcmp(s, "object") == 0) {
                *type = MLI_FRAME_TYPE_OBJECT;
        } else {
                chk_bad("Type is unknown.");
        }
        return 1;
chk_error:
        return 0;
}

void mli_Frame_print_walk(const struct mli_Frame *f, const uint64_t indention)
{
        uint64_t c;
        char type_string[1024];
        mli_frame_type_to_string(f->type, type_string);
        printf("%*s", (int)indention, "");
        printf(" __%s__ id:%u, at:%p\n", type_string, f->id, (void *)f);
        printf("%*s", (int)indention, "");
        printf("|-mother: id:");
        if (f->mother != NULL) {
                printf("%u,", f->mother->id);
        } else {
                printf("%p,", (void *)f->mother);
        }
        printf(" at:%p\n", (void *)f->mother);
        printf("%*s", (int)indention, "");
        printf("|-pos: (%0.1f, %0.1f, %0.1f)\n",
               f->frame2mother.translation.x,
               f->frame2mother.translation.y,
               f->frame2mother.translation.z);
        printf("%*s", (int)indention, "");
        printf("|-rotation: (%.1f| %.1f, %.1f, %.1f)\n",
               f->frame2mother.rotation.w,
               f->frame2mother.rotation.x,
               f->frame2mother.rotation.y,
               f->frame2mother.rotation.z);
        if (f->type == MLI_FRAME_TYPE_OBJECT) {
                uint32_t ii;
                printf("%*s", (int)indention, "");
                printf("|-boundary_layers [");
                for (ii = 0; ii < f->boundary_layers.size; ii++) {
                        uint32_t boundary_layer = f->boundary_layers.array[ii];
                        printf("%u,", boundary_layer);
                }
                printf("]\n");

                printf("%*s", (int)indention, "");
                printf("|-obj %u\n", f->object);
        }
        for (c = 0; c < f->children.size; c++) {
                struct mli_Frame *child = f->children.array[c];
                mli_Frame_print_walk(child, indention + 4);
        }
}

void mli_Frame_print(struct mli_Frame *f) { mli_Frame_print_walk(f, 0u); }

void mli_Frame_set_frame2root(struct mli_Frame *f)
{
        if (f->mother == NULL) {
                f->frame2root = f->frame2mother;
        } else {
                f->frame2root = mli_HomTraComp_sequence(
                        f->frame2mother, f->mother->frame2root);
        }
        if (f->type == MLI_FRAME_TYPE_FRAME) {
                uint64_t c;
                for (c = 0; c < f->children.size; c++) {
                        struct mli_Frame *child = f->children.array[c];
                        mli_Frame_set_frame2root(child);
                }
        }
}

int mli_Frame_estimate_num_robjects_and_total_num_boundary_layers_walk(
        const struct mli_Frame *frame,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers)
{
        uint64_t c;
        switch (frame->type) {
        case MLI_FRAME_TYPE_FRAME:
                for (c = 0; c < frame->children.size; c++) {
                        chk(mli_Frame_estimate_num_robjects_and_total_num_boundary_layers_walk(
                                frame->children.array[c],
                                num_robjects,
                                total_num_boundary_layers));
                }
                break;
        case MLI_FRAME_TYPE_OBJECT:
                (*num_robjects) += 1;
                (*total_num_boundary_layers) += frame->boundary_layers.size;
                break;
        default:
                chk_bad("Expected either type 'frame' or 'object'.");
                break;
        }
        return 1;
chk_error:
        return 0;
}

int mli_Frame_estimate_num_robjects_and_total_num_boundary_layers(
        const struct mli_Frame *frame,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers)
{
        (*num_robjects) = 0u;
        (*total_num_boundary_layers) = 0u;
        chk_msg(mli_Frame_estimate_num_robjects_and_total_num_boundary_layers_walk(
                        frame, num_robjects, total_num_boundary_layers),
                "Failed to walk tree of frames to estimate "
                "num_robjects and total_num_boundary_layers.");
        return 1;
chk_error:
        return 0;
}

/* frame_from_archive */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Frame_from_Archive(
        struct mli_Frame *root,
        const struct mli_Archive *archive,
        const struct mli_Map *object_names,
        const struct mli_Object *objects,
        const struct mli_Map *boundary_layer_names)
{
        uint64_t token = 0u;
        struct mli_String *fixload = NULL;
        struct mli_String fixname = mli_String_init();
        struct mli_Json tree_json = mli_Json_init();

        chk(mli_String_from_cstr(&fixname, "geometry/relations.json"));
        chk_msg(mli_Archive_get(archive, &fixname, &fixload),
                "Can not find geometry/relations.json file in archive.");
        chk_msg(mli_Json_from_string(&tree_json, fixload),
                "Failed to parse 'geometry/relations.json'.");

        chk_msg(mli_Json_token_by_key(&tree_json, 0, "children", &token),
                "Expected 'tree.json' to have key 'children'.");
        chk_msg(mli_Frame_malloc(root, MLI_FRAME_TYPE_FRAME),
                "Can not malloc root-frame.");
        chk_msg(mli_Frame_from_json(
                        root,
                        &tree_json,
                        token + 1,
                        object_names,
                        objects,
                        boundary_layer_names),
                "Failed to populate tree of Frames from 'tree.json'.");
        mli_Json_free(&tree_json);

        /* init transformations */
        mli_Frame_set_frame2root(root);

        mli_String_free(&fixname);
        return 1;
chk_error:
        mli_String_free(&fixname);
        mli_Json_free(&tree_json);
        return 0;
}

/* frame_json */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Frame_type_from_json_token(
        uint64_t *type,
        const struct mli_Json *json,
        const uint64_t token)
{
        uint64_t _t;
        int has_obj = mli_Json_token_by_key(json, token, "obj", &_t);
        int has_children = mli_Json_token_by_key(json, token, "children", &_t);

        if (has_obj && has_children) {
                chk_bad("Frame must not have both keys 'obj', and 'children'.");
        } else if (!has_obj && !has_children) {
                chk_bad("Frame must have either of keys 'obj', or 'children'.");
        } else if (has_obj && !has_children) {
                (*type) = MLI_FRAME_TYPE_OBJECT;
        } else if (!has_obj && has_children) {
                (*type) = MLI_FRAME_TYPE_FRAME;
        } else {
                chk_bad("Not expected to happen");
        }

        return 1;
chk_error:
        mli_Json_debug_token_fprint(stderr, json, token);
        return 0;
}

int mli_Frame_id_from_json_token(
        uint32_t *id,
        const struct mli_Json *json,
        const uint64_t token)
{
        uint64_t token_id;
        int64_t _id;
        chk_msg(mli_Json_token_by_key(json, token, "id", &token_id),
                "Expected Frame to have key 'id'.");
        chk_msg(mli_Json_int64_by_token(json, token_id + 1, &_id),
                "Failed to parse Frame's id.");
        chk_msg(_id >= 0, "Expected Frame's id >= 0.");
        (*id) = _id;

        return 1;
chk_error:
        mli_Json_debug_token_fprint(stderr, json, token);
        return 0;
}

int mli_Frame_pos_rot_from_json_token(
        struct mli_HomTraComp *frame2mother,
        const struct mli_Json *json,
        const uint64_t token)
{
        uint64_t token_pos, token_rot;
        /* pos */
        chk_msg(mli_Json_token_by_key(json, token, "pos", &token_pos),
                "Expected Frame to have key 'pos'.");
        chk_msg(mli_Vec_from_json_token(
                        &frame2mother->translation, json, token_pos + 1),
                "Failed to parse Frame's 'pos' from json.");

        /* rot */
        chk_msg(mli_Json_token_by_key(json, token, "rot", &token_rot),
                "Expected Frame to have key 'rot'.");
        chk_msg(mli_Quaternion_from_json(
                        &frame2mother->rotation, json, token_rot + 1),
                "Failed to parse Frame's 'rot' from json.");
        return 1;
chk_error:
        mli_Json_debug_token_fprint(stderr, json, token);
        return 0;
}

int mli_Frame_boundary_layers_form_json_token(
        struct mli_Uint32Vector *boundary_layers,
        const uint32_t object_idx,
        const struct mli_Object *objects,
        const struct mli_Map *boundary_layer_names,
        const struct mli_Json *json,
        const uint64_t token)
{
        uint64_t token_mtl_key, token_mtl;
        uint64_t material_idx;
        chk_msg(mli_Json_token_by_key(json, token, "mtl", &token_mtl_key),
                "Expected 'mtl' in Frame.");
        token_mtl = token_mtl_key + 1;
        chk_msg(json->tokens[token_mtl].type == JSMN_OBJECT,
                "Expected 'mtl' to be a json-object {}.");

        for (material_idx = 0u;
             material_idx < objects[object_idx].num_materials;
             material_idx++) {
                const char *material_key_in_object =
                        objects[object_idx].material_names[material_idx].array;

                uint64_t token_material_key = 0u;
                uint32_t boundary_layer_idx = 0u;
                chk_msg(mli_Json_token_by_key(
                                json,
                                token_mtl,
                                material_key_in_object,
                                &token_material_key),
                        "Expected object's material-key to be in "
                        "object-reference's mtls in tree.json.");

                chk_msg(mli_Map_get_value_for_string_from_json(
                                boundary_layer_names,
                                json,
                                token_material_key + 1,
                                &boundary_layer_idx),
                        "Expected boundary-layer to exist in materials.");

                chk_msg(mli_Uint32Vector_push_back(
                                boundary_layers, boundary_layer_idx),
                        "Failed to push-back boundary_layer_idx into "
                        "frame's boundary_layers.");
        }

        return 1;
chk_error:
        mli_Json_debug_token_fprint(stderr, json, token);
        return 0;
}

int mli_Frame_object_reference_form_json_token(
        uint32_t *object_reference,
        const struct mli_Json *json,
        const uint64_t token,
        const struct mli_Map *object_names)
{
        uint64_t token_obj_key;
        chk_msg(mli_Json_token_by_key(json, token, "obj", &token_obj_key),
                "Expected object to have key 'obj'.");
        chk_msg(mli_Map_get_value_for_string_from_json(
                        object_names,
                        json,
                        token_obj_key + 1,
                        object_reference),
                "Failed to get object-reference 'obj' from map");
        return 1;
chk_error:
        mli_Json_debug_token_fprint(stderr, json, token);
        return 0;
}

int mli_Frame_from_json(
        struct mli_Frame *mother,
        const struct mli_Json *json,
        const uint64_t token_children,
        const struct mli_Map *object_names,
        const struct mli_Object *objects,
        const struct mli_Map *boundary_layer_names)
{
        uint64_t num_children;
        uint64_t c;
        chk_msg(json->tokens[token_children].type == JSMN_ARRAY,
                "Expected Frame's children to be a json-array '[]'.");
        num_children = json->tokens[token_children].size;
        for (c = 0; c < num_children; c++) {
                uint64_t token_child = mli_Json__token_by_index_unsafe(
                        json, token_children, c);
                struct mli_Frame *child = NULL;
                uint64_t type;
                uint64_t token_grandchildren;

                chk_msg(mli_Frame_type_from_json_token(
                                &type, json, token_child),
                        "Failed to read type of Frame.");

                child = mli_Frame_add(mother, type);
                chk_msg(child, "Failed to add child to frame.");

                chk_msg(mli_Frame_pos_rot_from_json_token(
                                &child->frame2mother, json, token_child),
                        "Failed to set pos, and rot of Frame from json.");

                chk_msg(mli_Frame_id_from_json_token(
                                &child->id, json, token_child),
                        "Failed to set id of Frame from json.");

                switch (type) {
                case MLI_FRAME_TYPE_FRAME:
                        chk_msg(mli_Json_token_by_key(
                                        json,
                                        token_child,
                                        "children",
                                        &token_grandchildren),
                                "Expected child of type Frame to have "
                                "key 'children'.");
                        chk_msg(mli_Frame_from_json(
                                        child,
                                        json,
                                        token_grandchildren + 1,
                                        object_names,
                                        objects,
                                        boundary_layer_names),
                                "Failed to populate grandchildren "
                                "Frames from json.");
                        break;
                case MLI_FRAME_TYPE_OBJECT:
                        chk_msg(mli_Frame_object_reference_form_json_token(
                                        &child->object,
                                        json,
                                        token_child,
                                        object_names),
                                "Failed to parse object-reference "
                                "from json.");
                        chk_msg(mli_Frame_boundary_layers_form_json_token(
                                        &child->boundary_layers,
                                        child->object,
                                        objects,
                                        boundary_layer_names,
                                        json,
                                        token_child),
                                "Failed to set boundary_layers of Frame "
                                "from json.");
                        break;
                default:
                        chk_bad("Unknown type of frame.");
                        break;
                }
        }
        return 1;
chk_error:
        return 0;
}

/* frame_ptr_vector */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_VECTOR_IMPLEMENTATION(mli_FramePtrVector, struct mli_Frame *)

/* fresnel */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Fresnel mli_Fresnel_init(
        const struct mli_Vec incident,
        const struct mli_Vec normal,
        const double n_from,
        const double n_to)
{
        struct mli_Fresnel fresnel;
        fresnel.incident = incident;
        fresnel.normal = normal;
        fresnel.n_from = n_from;
        fresnel.n_to = n_to;

        fresnel._cosI = -1.0 * mli_Vec_dot(normal, incident);
        fresnel._n_from_over_n_to = n_from / n_to;
        fresnel._sinT2 =
                (fresnel._n_from_over_n_to * fresnel._n_from_over_n_to) *
                (1.0 - (fresnel._cosI * fresnel._cosI));
        fresnel._cosT = sqrt(1.0 - fresnel._sinT2);
        return fresnel;
}

double mli_Fresnel_reflection_propability(const struct mli_Fresnel fresnel)
{
        if (fresnel._sinT2 > 1.0) {
                /* total internal reflection */
                return 1.0;
        } else {
                const struct mli_Fresnel f = fresnel;
                const double nFromCosI = f.n_from * f._cosI;
                const double nFromCosT = f.n_from * f._cosT;
                const double nToCosI = f.n_to * f._cosI;
                const double nToCosT = f.n_to * f._cosT;
                const double rOrth =
                        (nFromCosI - nToCosT) / (nFromCosI + nToCosT);
                const double rPara =
                        (nToCosI - nFromCosT) / (nToCosI + nFromCosT);
                return (rOrth * rOrth + rPara * rPara) / 2.0;
        }
}

struct mli_Vec mli_Fresnel_reflection_direction(
        const struct mli_Fresnel fresnel)
{
        return mli_Vec_add(
                fresnel.incident,
                mli_Vec_multiply(fresnel.normal, fresnel._cosI * 2.0));
}

struct mli_Vec mli_Fresnel_refraction_direction(
        const struct mli_Fresnel fresnel)
{
        return mli_Vec_add(
                mli_Vec_multiply(fresnel.incident, fresnel._n_from_over_n_to),
                mli_Vec_multiply(
                        fresnel.normal,
                        fresnel._n_from_over_n_to * fresnel._cosI -
                                fresnel._cosT));
}

/* from_outside_to_inside */
/* ---------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_raytracing_from_outside_to_inside(
        const struct mli_Vec ray_direction_local,
        const struct mli_Vec surface_normal_local)
{
        const double proj =
                mli_Vec_dot(surface_normal_local, ray_direction_local);
        if (proj < 0.)
                return 1;
        else
                return 0;
}

/* func */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Func mli_Func_init(void)
{
        struct mli_Func f;
        f.num_points = 0u;
        f.x = NULL;
        f.y = NULL;
        return f;
}

void mli_Func_free(struct mli_Func *f)
{
        free(f->x);
        free(f->y);
        (*f) = mli_Func_init();
}

int mli_Func_malloc(struct mli_Func *f, const uint64_t num_points)
{
        mli_Func_free(f);
        f->num_points = num_points;
        chk_malloc(f->x, double, f->num_points);
        chk_malloc(f->y, double, f->num_points);
        return 1;
chk_error:
        mli_Func_free(f);
        return 0;
}

int mli_Func_x_is_strictly_increasing(const struct mli_Func *f)
{
        uint64_t i;
        for (i = 1; i < f->num_points; i++) {
                if (f->x[i] <= f->x[i - 1]) {
                        return 0;
                }
        }
        return 1;
}

int mli_Func_evaluate(const struct mli_Func *f, const double xarg, double *out)
{
        double y1, y0, x1, x0;
        uint64_t idx = mli_math_upper_compare_double(f->x, f->num_points, xarg);
        if (idx == 0) {
                chk_bad("mli_Func argument below lower bound.");
        } else if (idx == f->num_points) {
                chk_bad("mli_Func argument above upper bound.");
        } else {
                y1 = f->y[idx];
                y0 = f->y[idx - 1u];
                x1 = f->x[idx];
                x0 = f->x[idx - 1u];
                (*out) = mli_math_linear_interpolate_2d(xarg, x0, y0, x1, y1);
        }
        return 1;
chk_error:
        return 0;
}

int mli_Func_in_range(const struct mli_Func *f, const double xarg)
{
        if (f->num_points < 2) {
                return 0;
        }
        if (xarg >= f->x[0]) {
                if (xarg < f->x[f->num_points - 1]) {
                        return 1;
                }
        }
        return 0;
}

double mli_Func_evaluate_with_default_when_out_of_range(
        const struct mli_Func *f,
        const double xarg,
        const double default_value)
{
        double y1, y0, x1, x0;
        uint64_t idx = mli_math_upper_compare_double(f->x, f->num_points, xarg);
        if (idx == 0) {
                /* mli_Func argument below lower bound */
                return default_value;
        } else if (idx == f->num_points) {
                /* mli_Func argument above upper bound */
                return default_value;
        } else {
                y1 = f->y[idx];
                y0 = f->y[idx - 1u];
                x1 = f->x[idx];
                x0 = f->x[idx - 1u];
                return mli_math_linear_interpolate_2d(xarg, x0, y0, x1, y1);
        }
}

double mli_Func_evaluate_with_default_closest(
        const struct mli_Func *f,
        const double xarg)
{
        double y1, y0, x1, x0;
        uint64_t idx = mli_math_upper_compare_double(f->x, f->num_points, xarg);
        if (idx == 0) {
                /* mli_Func argument below lower bound */
                return f->y[0];
        } else if (idx == f->num_points) {
                /* mli_Func argument above upper bound */
                return f->y[f->num_points - 1];
        } else {
                y1 = f->y[idx];
                y0 = f->y[idx - 1u];
                x1 = f->x[idx];
                x0 = f->x[idx - 1u];
                return mli_math_linear_interpolate_2d(xarg, x0, y0, x1, y1);
        }
}

int mli_Func_fold_numeric(
        const struct mli_Func *a,
        const struct mli_Func *b,
        double *fold)
{
        uint64_t i;
        const uint64_t NUM_STEPS = 1024 * 8;
        const double xmin = a->x[0];
        const double xmax = a->x[a->num_points - 1];
        const double step_size = (xmax - xmin) / (double)NUM_STEPS;
        chk_msg(a->num_points >= 2u, "Expect a->num_points >= 2.");
        chk_msg(b->num_points >= 2u, "Expect b->num_points >= 2.");
        chk_msg(a->x[0] == b->x[0], "Expect a->x[0] == b->x[0].");
        chk_msg(a->x[a->num_points - 1] == b->x[b->num_points - 1],
                "Expect a->x[:-1] == b->x[:-1].");
        (*fold) = 0.0;
        for (i = 0; i < NUM_STEPS; i++) {
                double ra = MLI_MATH_NAN;
                double rb = MLI_MATH_NAN;
                double x = xmin + (double)i * step_size;
                chk(mli_Func_evaluate(a, x, &ra));
                chk(mli_Func_evaluate(b, x, &rb));
                (*fold) += (ra * rb) * step_size;
        }
        return 1;
chk_error:
        return 0;
}

int mli_Func_fold_numeric_default_closest(
        const struct mli_Func *a,
        const struct mli_Func *b,
        double *fold)
{
        double x_start, x_stop, x_step, x_range, x_weight;
        uint64_t i;
        const uint64_t NUM_STEPS = 1024 * 8;

        chk_msg(a->num_points >= 2u, "Expect a->num_points >= 2.");
        chk_msg(b->num_points >= 2u, "Expect b->num_points >= 2.");

        chk_msg(mli_Func_x_is_strictly_increasing(a),
                "Expected function a to be strictly_increasing.");
        chk_msg(mli_Func_x_is_strictly_increasing(b),
                "Expected function b to be strictly_increasing.");

        x_start = MLI_MATH_MAX2(a->x[0], b->x[0]);
        x_stop =
                MLI_MATH_MIN2(a->x[a->num_points - 1], b->x[b->num_points - 1]);
        x_range = x_stop - x_start;
        x_step = (x_range) / (double)NUM_STEPS;
        x_weight = x_step / x_range;

        (*fold) = 0.0;
        if (x_start < x_stop) {
                for (i = 0; i < NUM_STEPS; i++) {
                        double ra = MLI_MATH_NAN;
                        double rb = MLI_MATH_NAN;
                        double x = x_start + (double)i * x_step;
                        ra = mli_Func_evaluate_with_default_closest(a, x);
                        rb = mli_Func_evaluate_with_default_closest(b, x);
                        (*fold) += (ra * rb) * x_weight;
                }
        }

        return 1;
chk_error:
        return 0;
}

int mli_Func_equal(const struct mli_Func a, const struct mli_Func b)
{
        uint64_t i;
        if (a.num_points != b.num_points)
                return 0;
        for (i = 0; i < a.num_points; i++) {
                if (a.x[i] != b.x[i])
                        return 0;
                if (a.y[i] != b.y[i])
                        return 0;
        }
        return 1;
}

int mli_Func_is_valid(const struct mli_Func *func)
{
        uint64_t i;
        chk_msg(func->num_points >= 2,
                "Expected function to have at least two points. "
                "Evaluation is not possible when there is no valid range "
                "between two points.");

        for (i = 0; i < func->num_points; i++) {
                chk_msg(!MLI_MATH_IS_NAN(func->x[i]),
                        "Expected x-argument to be a real number, "
                        "but it is 'nan'.");
                chk_msg(!MLI_MATH_IS_NAN(func->y[i]),
                        "Expected y-value to be a real number, "
                        "but it is 'nan'.");
        }

        chk_msg(mli_Func_x_is_strictly_increasing(func),
                "Expected x-arguments to be strictly increasing, "
                "but they do not.");

        return 1;
chk_error:
        return 0;
}

int mli_Func_malloc_constant(
        struct mli_Func *self,
        const double start,
        const double stop,
        const double value)
{
        chk(mli_Func_malloc(self, 2));
        self->x[0] = start;
        self->y[0] = value;
        self->x[1] = stop;
        self->y[1] = value;

        return 1;
chk_error:
        return 0;
}

/* func_csv */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Func_from_csv_split_line(
        struct mli_String *line,
        struct mli_String *sx,
        struct mli_String *sy)
{
        size_t i;
        int num_commas = 0;
        chk(mli_String_malloc(sx, 16));
        chk(mli_String_malloc(sy, 16));
        for (i = 0; i < line->size; i++) {
                char c = line->array[i];
                if (c == ',') {
                        if (num_commas > 0) {
                                return 1;
                        }
                        num_commas += 1;
                } else {
                        if (num_commas > 0) {
                                chk(mli_String_push_back(sy, c));
                        } else {
                                chk(mli_String_push_back(sx, c));
                        }
                }
        }
        chk(mli_String_strip(sx, sx));
        chk(mli_String_strip(sy, sy));
        return 1;
chk_error:
        return 0;
}

int mli_Func_from_csv(
        struct mli_Func *func,
        struct mli_String *xname,
        struct mli_String *yname,
        struct mli_IO *io)
{
        uint64_t line_number = 0u;
        double _x, _y;
        struct mli_String line = mli_String_init();
        struct mli_String sx = mli_String_init();
        struct mli_String sy = mli_String_init();
        struct mli_DoubleVector x = mli_DoubleVector_init();
        struct mli_DoubleVector y = mli_DoubleVector_init();

        chk(mli_DoubleVector_malloc(&x, 0u));
        chk(mli_DoubleVector_malloc(&y, 0u));

        mli_Func_free(func);
        while (1) {
                line_number += 1;
                chk_msg(line_number < 1000 * 1000 * 1000,
                        "Expected less than 1e9 lines in csv file. "
                        "Something went wrong.");

                if (mli_IO_eof(io)) {
                        break;
                }
                chk_msg(mli_IO_text_read_line(io, &line, '\n'),
                        "Can not read line.");
                chk_msg(line.size < 4096, "Expected line.size < 4096.");

                if (line.size > 0) {
                        chk_msg(mli_Func_from_csv_split_line(&line, &sx, &sy),
                                "Can not split line into tokens.");
                        if (line_number == 1) {
                                chk(mli_String_copy(xname, &sx));
                                chk(mli_String_copy(yname, &sy));
                        } else {
                                chk(mli_String_to_double(&_x, &sx));
                                chk(mli_String_to_double(&_y, &sy));
                                chk(mli_DoubleVector_push_back(&x, _x));
                                chk(mli_DoubleVector_push_back(&y, _y));
                        }
                }
        }

        chk_msg(x.size == y.size, "Expected same number x, y values.");
        chk_msg(mli_Func_malloc(func, x.size),
                "Failed to malloc mli_Func from file.");

        MLI_MATH_NCPY(x.array, func->x, x.size);
        MLI_MATH_NCPY(y.array, func->y, y.size);

        mli_String_free(&line);
        mli_String_free(&sx);
        mli_String_free(&sy);
        mli_DoubleVector_free(&x);
        mli_DoubleVector_free(&y);
        return 1;
chk_error:
        mli_String_free(&line);
        mli_String_free(&sx);
        mli_String_free(&sy);
        mli_DoubleVector_free(&x);
        mli_DoubleVector_free(&y);

        mli_String_free(xname);
        mli_String_free(yname);
        mli_Func_free(func);
        return 0;
}

/* func_fprint */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Func_fprint(
        FILE *f,
        const struct mli_Func *func,
        struct mli_Func_fprint_Config cfg)
{
        int iy, ix;
        const double x_range = cfg.x_stop - cfg.x_start;
        const double y_range = cfg.y_stop - cfg.y_start;
        const double x_step = x_range / cfg.x_num;
        const double y_step = y_range / cfg.y_num;

        const int NUM_Y_SUB_STEPS = 3;
        chk(mli_Func_is_valid(func));

        for (iy = cfg.y_num - 1; iy >= 0; iy--) {
                if (iy == cfg.y_num - 1) {
                        fprintf(f, "  %-10.2e -|", cfg.y_stop);
                } else {
                        fprintf(f, "    %-10s|", "");
                }

                for (ix = 0; ix < cfg.x_num; ix++) {
                        double x, y;
                        double jy;
                        int jy_sub;
                        x = cfg.x_start + x_step * ix;
                        if (mli_Func_in_range(func, x)) {
                                chk(mli_Func_evaluate(func, x, &y));
                                jy = (y - cfg.y_start) / y_step;
                                jy_sub =
                                        (int)(NUM_Y_SUB_STEPS *
                                              (jy - floor(jy)));
                                if (iy == (int)(jy)) {
                                        switch (jy_sub) {
                                        case 0:
                                                fprintf(f, ",");
                                                break;
                                        case 1:
                                                fprintf(f, "-");
                                                break;
                                        case 2:
                                                fprintf(f, "'");
                                                break;
                                        default:
                                                fprintf(f, "?");
                                                break;
                                        }
                                } else {
                                        fprintf(f, " ");
                                }
                        }
                }
                fprintf(f, "\n");
        }
        fprintf(f, "  %-10.2e  +", cfg.y_start);
        for (ix = 0; ix < cfg.x_num; ix++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");

        fprintf(f, "          ");
        fprintf(f, "%-8.2e ", cfg.x_start);

        for (ix = 0; ix < cfg.x_num - 8; ix++) {
                fprintf(f, " ");
        }
        fprintf(f, "%-8.2e ", cfg.x_stop);
        fprintf(f, "\n");
        return 1;
chk_error:
        return 0;
}

/* func_info */
/* --------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

struct mli_FuncInfo mli_FuncInfo_init(void)
{
        struct mli_FuncInfo out;
        out.x = mli_String_init();
        out.y = mli_String_init();
        out.comment = mli_String_init();
        return out;
}

void mli_FuncInfo_free(struct mli_FuncInfo *self)
{
        mli_String_free(&self->x);
        mli_String_free(&self->y);
        mli_String_free(&self->comment);
        (*self) = mli_FuncInfo_init();
}

int mli_FuncInfo_malloc(struct mli_FuncInfo *self)
{
        chk(mli_String_malloc(&self->x, 0u));
        chk(mli_String_malloc(&self->y, 0u));
        chk(mli_String_malloc(&self->comment, 0u));
        return 1;
chk_error:
        return 0;
}

int mli_FuncInfo_equal(
        const struct mli_FuncInfo *a,
        const struct mli_FuncInfo *b)
{
        if (!mli_String_equal(&a->x, &b->x))
                return 0;
        if (!mli_String_equal(&a->y, &b->y))
                return 0;
        if (!mli_String_equal(&a->comment, &b->comment))
                return 0;
        return 1;
}

int mli_FuncInfo_to_io(const struct mli_FuncInfo *self, struct mli_IO *f)
{
        struct mli_MagicId magic = mli_MagicId_init();
        chk(mli_MagicId_set(&magic, "mli_FuncInfo"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_String_to_io(&self->x, f));
        chk(mli_String_to_io(&self->y, f));
        chk(mli_String_to_io(&self->comment, f));
        return 1;
chk_error:
        return 0;
}

int mli_FuncInfo_from_io(struct mli_FuncInfo *self, struct mli_IO *f)
{
        struct mli_MagicId magic = mli_MagicId_init();
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_FuncInfo"));
        mli_MagicId_warn_version(&magic);
        chk(mli_String_from_io(&self->x, f));
        chk(mli_String_from_io(&self->y, f));
        chk(mli_String_from_io(&self->comment, f));
        return 1;
chk_error:
        return 0;
}

/* func_serialize */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Func_to_io(const struct mli_Func *func, struct mli_IO *f)
{
        struct mli_MagicId magic;
        chk(mli_MagicId_set(&magic, "mli_Func"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        chk_IO_write(&func->num_points, sizeof(uint64_t), 1u, f);
        chk_IO_write(func->x, sizeof(double), func->num_points, f);
        chk_IO_write(func->y, sizeof(double), func->num_points, f);
        return 1;
chk_error:
        return 0;
}

int mli_Func_from_io(struct mli_Func *func, struct mli_IO *f)
{
        uint64_t num_points;
        struct mli_MagicId magic;
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Func"));
        mli_MagicId_warn_version(&magic);

        chk_IO_read(&num_points, sizeof(uint64_t), 1u, f);
        chk(mli_Func_malloc(func, num_points));
        chk_IO_read(func->x, sizeof(double), func->num_points, f);
        chk_IO_read(func->y, sizeof(double), func->num_points, f);
        chk_msg(mli_Func_is_valid(func), "Expected function to be valid.");
        return 1;
chk_error:
        mli_Func_free(func);
        return 0;
}

/* geometry */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mli_Geometry_init_objects(struct mli_Geometry *self)
{
        self->num_objects = 0u;
        self->objects = NULL;
        self->object_names = NULL;
}

void mli_Geometry_init_references(struct mli_Geometry *self)
{
        self->num_robjects = 0u;
        self->robjects = NULL;
        self->robject_ids = NULL;
        self->robject2root = NULL;
}

struct mli_Geometry mli_Geometry_init(void)
{
        struct mli_Geometry out;
        mli_Geometry_init_objects(&out);
        mli_Geometry_init_references(&out);
        return out;
}

void mli_Geometry_free_objects(struct mli_Geometry *self)
{
        uint32_t i;
        for (i = 0; i < self->num_objects; i++) {
                mli_Object_free(&(self->objects[i]));
                mli_String_free(&(self->object_names[i]));
        }
        free(self->objects);
        free(self->object_names);
        mli_Geometry_init_objects(self);
}

void mli_Geometry_free_references(struct mli_Geometry *self)
{
        free(self->robjects);
        free(self->robject_ids);
        free(self->robject2root);
        mli_Geometry_init_references(self);
}

void mli_Geometry_free(struct mli_Geometry *self)
{
        mli_Geometry_free_objects(self);
        mli_Geometry_free_references(self);
        (*self) = mli_Geometry_init();
}

int mli_Geometry_malloc_objects(
        struct mli_Geometry *self,
        const uint32_t num_objects)
{
        uint32_t i;
        self->num_objects = num_objects;
        chk_malloc(self->objects, struct mli_Object, self->num_objects);
        chk_malloc(self->object_names, struct mli_String, self->num_objects);
        for (i = 0; i < self->num_objects; i++) {
                self->objects[i] = mli_Object_init();
                self->object_names[i] = mli_String_init();
        }
        return 1;
chk_error:
        mli_Geometry_free_objects(self);
        return 0;
}

int mli_Geometry_malloc_references(
        struct mli_Geometry *self,
        const uint32_t num_robjects)
{
        self->num_robjects = num_robjects;
        chk_malloc(self->robjects, uint32_t, self->num_robjects);
        chk_malloc(self->robject_ids, uint32_t, self->num_robjects);
        chk_malloc(
                self->robject2root, struct mli_HomTraComp, self->num_robjects);
        return 1;
chk_error:
        mli_Geometry_free_references(self);
        return 0;
}

int mli_Geometry_malloc(
        struct mli_Geometry *self,
        const uint32_t num_objects,
        const uint32_t num_robjects)
{
        mli_Geometry_free(self);
        chk(mli_Geometry_malloc_objects(self, num_objects));
        chk(mli_Geometry_malloc_references(self, num_robjects));
        return 1;
chk_error:
        mli_Geometry_free(self);
        return 0;
}

void mli_Geometry_info_fprint(FILE *f, const struct mli_Geometry *self)
{
        uint32_t rob, i;
        fprintf(f, "geometry\n");
        fprintf(f, "--------\n");
        fprintf(f, "\n");
        fprintf(f, "    objects\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        fprintf(f, "    ");
        fprintf(f, "%5s ", "obj");
        fprintf(f, "%24s ", "name");
        fprintf(f, "%8s ", "#v");
        fprintf(f, "%8s ", "#vn");
        fprintf(f, "%8s ", "#f");
        fprintf(f, "%8s ", "#mtl");
        fprintf(f, "\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");

        for (i = 0; i < self->num_objects; i++) {
                fprintf(f, "    ");
                fprintf(f, "%5d ", i);
                fprintf(f, "%24s ", self->object_names[i].array);
                fprintf(f, "%8d ", self->objects[i].num_vertices);
                fprintf(f, "%8d ", self->objects[i].num_vertex_normals);
                fprintf(f, "%8d ", self->objects[i].num_faces);
                fprintf(f, "%8d ", self->objects[i].num_materials);
                fprintf(f, "\n");
        }
        fprintf(f, "\n");

        fprintf(f, "\n");
        fprintf(f, "    object-references\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        fprintf(f, "    ");
        fprintf(f, "%6s ", "ref");
        fprintf(f, "%6s ", "obj");
        fprintf(f, "%6s ", "id");
        fprintf(f, "%20s ", "translation(xyz)/m");
        fprintf(f, "%28s ", "quarternion(xyz;w)");
        fprintf(f, "\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");

        for (rob = 0; rob < self->num_robjects; rob++) {
                fprintf(f, "    ");
                fprintf(f, "%6d ", rob);
                fprintf(f, "%6d ", self->robjects[rob]);
                fprintf(f, "%6d ", self->robject_ids[rob]);
                fprintf(f, "% 6.1f,", self->robject2root[rob].translation.x);
                fprintf(f, "% 6.1f,", self->robject2root[rob].translation.y);
                fprintf(f, "% 6.1f ", self->robject2root[rob].translation.z);
                fprintf(f, " ");
                fprintf(f, "% 6.1f,", self->robject2root[rob].rotation.x);
                fprintf(f, "% 6.1f,", self->robject2root[rob].rotation.y);
                fprintf(f, "% 6.1f;", self->robject2root[rob].rotation.z);
                fprintf(f, "% 6.1f ", self->robject2root[rob].rotation.w);
                fprintf(f, "\n");
        }
}

int mli_Geometry_warn_objects(const struct mli_Geometry *self)
{
        uint32_t o;
        for (o = 0; o < self->num_objects; o++) {
                uint32_t v, vn, mtl;
                chk_msg(mli_Object_num_unused(
                                &(self->objects)[o], &v, &vn, &mtl),
                        "Can't estimate num unused v/vn/mtl in object.");
                if (v > 0 || vn > 0 || mtl > 0) {
                        fprintf(stderr,
                                "[WARNING] Object '%s' at [%u] ",
                                self->object_names[o].array,
                                o);
                }
                if (v > 0) {
                        fprintf(stderr, "has %u unused vertices (v).\n", v);
                }
                if (vn > 0) {
                        fprintf(stderr,
                                "has %u unused vertex-normals (vn).\n",
                                vn);
                }
                if (mtl > 0) {
                        fprintf(stderr,
                                "has %u unused materials (mtl).\n",
                                mtl);
                }
        }

        return 1;
chk_error:
        return 0;
}

/* geometry_aabb */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Geometry_robject_has_overlap_aabb(
        const struct mli_GeometryAndAccelerator *accgeo,
        const uint32_t robject_idx,
        const struct mli_AABB aabb)
{
        const struct mli_AABB robject_aabb =
                accgeo->accelerator->robject_aabbs[robject_idx];
        if (mli_AABB_is_overlapping(aabb, robject_aabb)) {
                /* Only test the object's faces when its own aabb has an
                 * overlap with the aabb to test for.
                 */
                const uint32_t obj_idx =
                        accgeo->geometry->robjects[robject_idx];
                const struct mli_Object *obj_ptr =
                        &accgeo->geometry->objects[obj_idx];
                const struct mli_HomTra robj2root = mli_HomTraComp_from_compact(
                        accgeo->geometry->robject2root[robject_idx]);
                return mli_Object_has_overlap_aabb(obj_ptr, robj2root, aabb);
        } else {
                return 0;
        }
}

int mli_Geometry_robject_has_overlap_aabb_void(
        const void *accgeo,
        const uint32_t robject_idx,
        const struct mli_AABB aabb)
{
        return mli_Geometry_robject_has_overlap_aabb(
                (const struct mli_GeometryAndAccelerator *)accgeo,
                robject_idx,
                aabb);
}

/* geometry_and_accelerator */
/* ------------------------ */



/* geometry_equal */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Geometry_objects_equal(
        const struct mli_Geometry *a,
        const struct mli_Geometry *b)
{
        uint32_t i = 0u;
        chk_msg(a->num_objects == b->num_objects,
                "Expected num_objects to be equal.");

        for (i = 0; i < a->num_objects; i++) {
                chk_msg(mli_Object_equal(&a->objects[i], &b->objects[i]),
                        "Expected object to be equal.");
                chk_msg(mli_String_equal(
                                &a->object_names[i], &b->object_names[i]),
                        "Expected object_name to be equal.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In geometry.object[%u]\n", i);
        return 0;
}

int mli_Geometry_object_references_equal(
        const struct mli_Geometry *a,
        const struct mli_Geometry *b)
{
        uint64_t rob = 0u;
        chk_msg(a->num_robjects == b->num_robjects,
                "Expected num_robjects to be equal.");

        for (rob = 0; rob < a->num_robjects; rob++) {
                chk_msg(a->robjects[rob] == b->robjects[rob],
                        "Expected object-references to be equal.");
                chk_msg(a->robject_ids[rob] == b->robject_ids[rob],
                        "Expected the users object-ids to be equal.");

                chk_msg(a->robject_ids[rob] == b->robject_ids[rob],
                        "Expected the users object-ids to be equal.");

                chk_msg(mli_HomTraComp_equal(
                                a->robject2root[rob], b->robject2root[rob]),
                        "Expected homogenous transformation of "
                        "object-references to be equal");
        }
        return 1;
chk_error:
        fprintf(stderr, "In geometry.object_reference[%lu]\n", rob);
        return 0;
}

int mli_Geometry_equal(
        const struct mli_Geometry *a,
        const struct mli_Geometry *b)
{
        chk_msg(mli_Geometry_objects_equal(a, b),
                "Expected objects to be equal.");
        chk_msg(mli_Geometry_object_references_equal(a, b),
                "Expected object-references to be equal.");
        return 1;
chk_error:
        return 0;
}

/* geometry_from_archive */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Geometry_from_archive(
        struct mli_Geometry *geometry,
        struct mli_Map *object_names,
        const struct mli_Archive *archive)
{
        uint64_t arc_idx = 0u;
        uint64_t obj_idx = 0u;
        struct mli_String key = mli_String_init();
        struct mli_String basename = mli_String_init();
        struct mli_String extension = mli_String_init();

        chk_dbg;
        /* objects */
        obj_idx = 0u;
        for (arc_idx = 0u; arc_idx < mli_Archive_size(archive); arc_idx++) {
                struct mli_String *filename =
                        &archive->filenames.items.array[arc_idx].key;

                chk_dbg;
                if (mli_String_starts_with_cstr(
                            filename, "geometry/objects/") &&
                    mli_String_ends_with_cstr(filename, ".obj")) {
                        struct mli_String *payload =
                                &archive->textfiles.array[arc_idx];
                        struct mli_IO buff = mli_IO_init();
                        chk_dbg;
                        chk(mli_IO_open_memory(&buff));
                        chk(mli_IO_text_write_String(&buff, payload));
                        mli_IO_rewind(&buff);

                        chk_dbg;
                        chk_msg(obj_idx < geometry->num_objects,
                                "Expected less objects in archive.");
                        chk_dbg;
                        chk(mli_path_basename(filename, &basename));
                        chk(mli_path_splitext(&basename, &key, &extension));

                        chk_msg(mli_Map_insert(object_names, &key, obj_idx),
                                "Failed to insert object-filename into map.");
                        chk_dbg;
                        chk_msg(mli_Object_malloc_from_wavefront(
                                        &geometry->objects[obj_idx], &buff),
                                "Failed to parse wave-front-object.");
                        chk_dbg;
                        chk(mli_String_copy(
                                &geometry->object_names[obj_idx], &key));
                        obj_idx += 1u;

                        mli_IO_close(&buff);
                }
        }
        chk_dbg;

        mli_String_free(&extension);
        mli_String_free(&basename);
        mli_String_free(&key);
        return 1;
chk_error:
        mli_String_free(&extension);
        mli_String_free(&basename);
        mli_String_free(&key);
        return 0;
}

/* geometry_id */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_GeometryId mli_GeometryId_init(void)
{
        struct mli_GeometryId id;
        id.robj = 0u;
        id.face = 0u;
        return id;
}

/* geometry_serialize */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Geometry_from_io(struct mli_Geometry *geometry, struct mli_IO *f)
{
        uint32_t i;
        uint32_t num_objects = 0u;
        uint32_t num_robjects = 0u;
        struct mli_MagicId magic;

        /* magic identifier */
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Geometry"));
        mli_MagicId_warn_version(&magic);

        /* payload */
        chk_IO_read(&num_objects, sizeof(uint32_t), 1u, f);
        chk_IO_read(&num_robjects, sizeof(uint32_t), 1u, f);

        chk_msg(mli_Geometry_malloc(geometry, num_objects, num_robjects),
                "Failed to malloc robjects in mli_Geometry.");

        for (i = 0; i < geometry->num_objects; i++) {
                chk_msg(mli_Object_from_io(&geometry->objects[i], f),
                        "Failed to read object into geometry.");
                chk_msg(mli_String_from_io(&geometry->object_names[i], f),
                        "Failed to read object name into geometry.");
        }

        chk_IO_read(
                geometry->robjects,
                sizeof(uint32_t),
                geometry->num_robjects,
                f);
        chk_IO_read(
                geometry->robject_ids,
                sizeof(uint32_t),
                geometry->num_robjects,
                f);
        chk_IO_read(
                geometry->robject2root,
                sizeof(struct mli_HomTraComp),
                geometry->num_robjects,
                f);

        return 1;
chk_error:
        return 0;
}

int mli_Geometry_to_io(const struct mli_Geometry *geometry, struct mli_IO *f)
{
        uint32_t i;
        struct mli_MagicId magic;

        /* magic identifier */
        chk(mli_MagicId_set(&magic, "mli_Geometry"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        /* payload */
        chk_IO_write(&geometry->num_objects, sizeof(uint32_t), 1, f);
        chk_IO_write(&geometry->num_robjects, sizeof(uint32_t), 1, f);

        for (i = 0; i < geometry->num_objects; i++) {
                chk_msg(mli_Object_to_io(&geometry->objects[i], f),
                        "Failed to write objects.");
                chk_msg(mli_String_to_io(&geometry->object_names[i], f),
                        "Failed to write object name.");
        }

        chk_IO_write(
                geometry->robjects,
                sizeof(uint32_t),
                geometry->num_robjects,
                f);
        chk_IO_write(
                geometry->robject_ids,
                sizeof(uint32_t),
                geometry->num_robjects,
                f);
        chk_IO_write(
                geometry->robject2root,
                sizeof(struct mli_HomTraComp),
                geometry->num_robjects,
                f);

        return 1;
chk_error:
        return 0;
}

/* geometry_set_from_frame */
/* ----------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Geometry__set_robjects_and_material_map_from_frame_walk(
        const struct mli_Frame *frame,
        struct mli_Geometry *geometry,
        struct mli_GeometryToMaterialMap *geomap,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers)
{
        uint64_t c;
        uint64_t material_idx;
        uint64_t robject_idx;
        switch (frame->type) {
        case MLI_FRAME_TYPE_FRAME:
                for (c = 0; c < frame->children.size; c++) {
                        chk(mli_Geometry__set_robjects_and_material_map_from_frame_walk(
                                frame->children.array[c],
                                geometry,
                                geomap,
                                num_robjects,
                                total_num_boundary_layers));
                }
                break;
        case MLI_FRAME_TYPE_OBJECT:
                robject_idx = (*num_robjects);

                chk_msg(frame->object < geometry->num_objects,
                        "Expected frame->object < num_objects.");
                /* geometry */
                geometry->robjects[robject_idx] = frame->object;
                geometry->robject2root[robject_idx] = frame->frame2root;
                geometry->robject_ids[robject_idx] = frame->id;

                /* materials map */
                chk_msg(frame->boundary_layers.size ==
                                geometry->objects[frame->object].num_materials,
                        "Expected Frame to have same "
                        "num boundary_layers as object.");

                geomap->first_boundary_layer_in_robject[robject_idx] =
                        (*total_num_boundary_layers);
                for (material_idx = 0;
                     material_idx < frame->boundary_layers.size;
                     material_idx++) {
                        mli_GeometryToMaterialMap_set(
                                geomap,
                                robject_idx,
                                material_idx,
                                frame->boundary_layers.array[material_idx]);
                        (*total_num_boundary_layers) += 1;
                }

                (*num_robjects) += 1;
                break;
        default:
                chk_bad("Expected either type 'frame' or 'object'.");
                break;
        }
        return 1;
chk_error:
        return 0;
}

int mli_Geometry_set_robjects_and_material_map_from_frame(
        const struct mli_Frame *frame,
        struct mli_Geometry *geometry,
        struct mli_GeometryToMaterialMap *geomap)
{
        uint64_t num_robjects = 0u;
        uint64_t total_num_boundary_layers = 0u;
        chk_msg(mli_Geometry__set_robjects_and_material_map_from_frame_walk(
                        frame,
                        geometry,
                        geomap,
                        &num_robjects,
                        &total_num_boundary_layers),
                "Failed to walk tree of frames to set "
                "robjects and material map.");
        return 1;
chk_error:
        return 0;
}

/* geometry_valid */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Geometry_valid_objects(const struct mli_Geometry *geometry)
{
        uint32_t i = 0;
        for (i = 0; i < geometry->num_objects; i++) {

                struct mli_String *name = &geometry->object_names[i];
                int64_t size = mli_String__discover_size(name);
                chk_msg(size > 0,
                        "Expected object_names to have at least "
                        "size '1' and to be '\\0' terminated.");
                chk_msg(size == (int64_t)name->size,
                        "Expected object_names size to "
                        "match zero termination.");

                chk_msg(mli_Object_is_valid(&geometry->objects[i]),
                        "Expected object to be valid.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In geometry.objects[%u]\n", i);
        return 0;
}

int mli_Geometry_valid_robjects_HomTras(const struct mli_Geometry *geometry)
{
        uint32_t i;
        for (i = 0; i < geometry->num_robjects; i++) {
                const struct mli_Vec t = geometry->robject2root[i].translation;
                const struct mli_Quaternion q =
                        geometry->robject2root[i].rotation;
                chk_msg(!MLI_MATH_IS_NAN(t.x), "translation.x is 'nan'.");
                chk_msg(!MLI_MATH_IS_NAN(t.y), "translation.y is 'nan'.");
                chk_msg(!MLI_MATH_IS_NAN(t.z), "translation.z is 'nan'.");

                chk_msg(!MLI_MATH_IS_NAN(q.w), "quaternion.w is 'nan'.");
                chk_msg(!MLI_MATH_IS_NAN(q.x), "quaternion.x is 'nan'.");
                chk_msg(!MLI_MATH_IS_NAN(q.y), "quaternion.y is 'nan'.");
                chk_msg(!MLI_MATH_IS_NAN(q.z), "quaternion.z is 'nan'.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In geometry.robject2root[%u]\n", i);
        return 0;
}

int mli_Geometry_valid_object_references(const struct mli_Geometry *geometry)
{
        uint32_t i;
        for (i = 0; i < geometry->num_robjects; i++) {
                chk_msg(geometry->robjects[i] < geometry->num_objects,
                        "Expected robjects to refer to valid objects.");
                /*
                 *       robject_ids are set by the user and can be whatever
                 *       the user wants.
                 */
        }
        return 1;
        fprintf(stderr, "In geometry.robject[%u]\n", i);
chk_error:
        return 0;
}

int mli_Geometry_valid(const struct mli_Geometry *geometry)
{
        chk_msg(mli_Geometry_valid_objects(geometry),
                "Expected objects to be valid.");
        chk_msg(mli_Geometry_valid_robjects_HomTras(geometry),
                "Expected robject transformations to be free of 'nan'.");
        chk_msg(mli_Geometry_valid_object_references(geometry),
                "Expected object-references to be valid.");
        return 1;
chk_error:
        return 0;
}

/* geometrytomaterialmap */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_GeometryToMaterialMap mli_GeometryToMaterialMap_init(void)
{
        struct mli_GeometryToMaterialMap map;
        map.num_robjects = 0u;
        map.total_num_boundary_layers = 0u;
        map.boundary_layers = NULL;
        map.first_boundary_layer_in_robject = NULL;
        return map;
}

void mli_GeometryToMaterialMap_free(struct mli_GeometryToMaterialMap *map)
{
        free(map->first_boundary_layer_in_robject);
        free(map->boundary_layers);
        (*map) = mli_GeometryToMaterialMap_init();
}

int mli_GeometryToMaterialMap_malloc(
        struct mli_GeometryToMaterialMap *map,
        const uint32_t num_robjects,
        const uint32_t total_num_boundary_layers)
{
        mli_GeometryToMaterialMap_free(map);
        map->num_robjects = num_robjects;
        map->total_num_boundary_layers = total_num_boundary_layers;

        chk_malloc(
                map->boundary_layers, uint32_t, map->total_num_boundary_layers);
        chk_malloc(
                map->first_boundary_layer_in_robject,
                uint32_t,
                map->num_robjects);
        MLI_MATH_ARRAY_SET(
                map->first_boundary_layer_in_robject, 0, map->num_robjects);
        MLI_MATH_ARRAY_SET(
                map->boundary_layers, 0, map->total_num_boundary_layers);
        return 1;
chk_error:
        return 0;
}

uint32_t mli_GeometryToMaterialMap_resolve_idx(
        const struct mli_GeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx)
{
        uint32_t idx = map->first_boundary_layer_in_robject[robject_idx];
        idx += material_idx;
        return idx;
}

uint32_t mli_GeometryToMaterialMap_get(
        const struct mli_GeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx)
{
        return map->boundary_layers[mli_GeometryToMaterialMap_resolve_idx(
                map, robject_idx, material_idx)];
}

void mli_GeometryToMaterialMap_set(
        const struct mli_GeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx,
        const uint32_t boundary_layer_idx)
{
        map->boundary_layers[mli_GeometryToMaterialMap_resolve_idx(
                map, robject_idx, material_idx)] = boundary_layer_idx;
}

uint32_t mli_GeometryToMaterialMap_num_boundary_layers_in_robject(
        const struct mli_GeometryToMaterialMap *map,
        const uint32_t robject_idx)
{
        const uint32_t start =
                mli_GeometryToMaterialMap_resolve_idx(map, robject_idx, 0u);
        uint32_t end = start;
        if (robject_idx + 1 < map->num_robjects) {
                end = mli_GeometryToMaterialMap_resolve_idx(
                        map, robject_idx + 1, 0u);
        } else {
                end = map->total_num_boundary_layers;
        }
        return end - start;
}

void mli_GeometryToMaterialMap_info_fprint(
        FILE *f,
        const struct mli_GeometryToMaterialMap *map)
{
        uint32_t robj, bdl, bdl_start, num_bdls, i;
        fprintf(f, "    geometry to material map\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        fprintf(f, "    ");
        fprintf(f, "%5s ", "ref");
        fprintf(f, "%24s ", "boundary-layers");
        fprintf(f, "\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");

        for (robj = 0u; robj < map->num_robjects; robj++) {
                bdl_start =
                        mli_GeometryToMaterialMap_resolve_idx(map, robj, 0u);
                num_bdls =
                        mli_GeometryToMaterialMap_num_boundary_layers_in_robject(
                                map, robj);
                fprintf(f, "    ");
                fprintf(f, "% 5d ", robj);
                fprintf(f, "  [");
                for (bdl = 0u; bdl < num_bdls; bdl++) {
                        fprintf(f,
                                "%d,",
                                map->boundary_layers[bdl_start + bdl]);
                }
                fprintf(f, "]\n");
        }
}

/* geometrytomaterialmap_equal */
/* --------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_GeometryToMaterialMap_equal(
        const struct mli_GeometryToMaterialMap *a,
        const struct mli_GeometryToMaterialMap *b)
{
        uint32_t i = 0u;
        chk_msg(a->num_robjects == b->num_robjects,
                "Expected num_robjects to be equal.");
        chk_msg(a->total_num_boundary_layers == b->total_num_boundary_layers,
                "Expected total_num_boundary_layers to be equal.");
        for (i = 0u; i < a->total_num_boundary_layers; i++) {
                chk_msg(a->boundary_layers[i] == b->boundary_layers[i],
                        "Expected all boundary_layers to be equal.");
        }
        for (i = 0u; i < a->num_robjects; i++) {
                chk_msg(a->first_boundary_layer_in_robject[i] ==
                                b->first_boundary_layer_in_robject[i],
                        "Expected all first_boundary_layer_in_robject "
                        "to be equal.");
        }
        return 1;
chk_error:
        return 0;
}

/* geometrytomaterialmap_serialize */
/* ------------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_GeometryToMaterialMap_from_io(
        struct mli_GeometryToMaterialMap *geomap,
        struct mli_IO *f)
{
        uint32_t num_robjects = 0u;
        uint32_t total_num_boundary_layers = 0u;
        struct mli_MagicId magic;

        /* magic identifier */
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_GeometryToMaterialMap"));
        mli_MagicId_warn_version(&magic);

        /* payload */
        chk_IO_read(&num_robjects, sizeof(uint32_t), 1u, f);
        chk_IO_read(&total_num_boundary_layers, sizeof(uint32_t), 1u, f);
        chk_msg(mli_GeometryToMaterialMap_malloc(
                        geomap, num_robjects, total_num_boundary_layers),
                "Failed to malloc mli_GeometryToMaterialMap.");
        chk_IO_read(
                geomap->boundary_layers,
                sizeof(uint32_t),
                geomap->total_num_boundary_layers,
                f);
        chk_IO_read(
                geomap->first_boundary_layer_in_robject,
                sizeof(uint32_t),
                geomap->num_robjects,
                f);
        return 1;
chk_error:
        return 0;
}

int mli_GeometryToMaterialMap_to_io(
        const struct mli_GeometryToMaterialMap *geomap,
        struct mli_IO *f)
{
        struct mli_MagicId magic;

        /* magic identifier */
        chk(mli_MagicId_set(&magic, "mli_GeometryToMaterialMap"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        /* payload */
        chk_IO_write(&geomap->num_robjects, sizeof(uint32_t), 1, f);
        chk_IO_write(
                &geomap->total_num_boundary_layers, sizeof(uint32_t), 1, f);
        chk_IO_write(
                geomap->boundary_layers,
                sizeof(uint32_t),
                geomap->total_num_boundary_layers,
                f);
        chk_IO_write(
                geomap->first_boundary_layer_in_robject,
                sizeof(uint32_t),
                geomap->num_robjects,
                f);
        return 1;
chk_error:
        return 0;
}

/* geometrytomaterialmap_valid */
/* --------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_GeometryToMaterialMap_valid(
        const struct mli_GeometryToMaterialMap *geomap)
{
        uint32_t i = 0u;
        chk_msg(geomap->total_num_boundary_layers >= geomap->num_robjects,
                "Expected at least as many boundary-layer-references as "
                "object-references.");
        for (i = 0u; i < geomap->num_robjects; i++) {
                chk_msg(geomap->first_boundary_layer_in_robject[i] <
                                geomap->total_num_boundary_layers,
                        "Expected all position of first_boundary_layer[i] < "
                        "total_num_boundary_layers");
        }
        return 1;
chk_error:
        return 0;
}

int mli_GeometryToMaterialMap_valid_wrt_Geometry(
        const struct mli_GeometryToMaterialMap *geomap,
        const struct mli_Geometry *geometry)
{
        uint32_t robj = 0u;
        uint32_t total_num_boundary_layers = 0u;
        chk_msg(geomap->num_robjects == geometry->num_robjects,
                "Expected num_robjects to be equal in Geometry and GeoMtlMap.");

        for (robj = 0u; robj < geomap->num_robjects; robj++) {
                const uint32_t obj = geometry->robjects[robj];
                const uint32_t obj_num_materials =
                        geometry->objects[obj].num_materials;
                chk_msg(mli_GeometryToMaterialMap_num_boundary_layers_in_robject(
                                geomap, robj) == obj_num_materials,
                        "Expected robject to have same num boundary-layers.");
                total_num_boundary_layers += obj_num_materials;
        }
        chk_msg(total_num_boundary_layers == geomap->total_num_boundary_layers,
                "Expected total_num_boundary_layers to match the Geometry.");

        return 1;
chk_error:
        return 0;
}

int mli_GeometryToMaterialMap_valid_wrt_Materials(
        const struct mli_GeometryToMaterialMap *geomap,
        const struct mli_Materials *materials)
{
        uint32_t i = 0u;
        for (i = 0u; i < geomap->total_num_boundary_layers; i++) {
                chk_msg(geomap->boundary_layers[i] <
                                materials->boundary_layers.size,
                        "Expected geomap's boundary_layers[i] to refer to "
                        "a valid boundary_layer in Materials.");
        }
        return 1;
chk_error:
        return 0;
}

/* homtra */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_HomTraComp mli_HomTraComp_set(
        const struct mli_Vec translation,
        const struct mli_Quaternion rotation)
{
        struct mli_HomTraComp comp;
        comp.translation = translation;
        comp.rotation = rotation;
        return comp;
}

struct mli_HomTra mli_HomTraComp_from_compact(const struct mli_HomTraComp trafo)
{
        struct mli_HomTra t;
        t.translation = trafo.translation;
        t.rotation = mli_Quaternion_to_matrix(trafo.rotation);
        return t;
}

struct mli_Vec mli_transform_orientation(
        const struct mli_Mat *rotation,
        const struct mli_Vec ori)
{
        struct mli_Vec out;
        out.x =
                (ori.x * rotation->r00 + ori.y * rotation->r01 +
                 ori.z * rotation->r02);
        out.y =
                (ori.x * rotation->r10 + ori.y * rotation->r11 +
                 ori.z * rotation->r12);
        out.z =
                (ori.x * rotation->r20 + ori.y * rotation->r21 +
                 ori.z * rotation->r22);
        return out;
}

struct mli_Vec mli_transform_orientation_inverse(
        const struct mli_Mat *rotation,
        const struct mli_Vec ori)
{
        struct mli_Vec out;
        out.x =
                (ori.x * rotation->r00 + ori.y * rotation->r10 +
                 ori.z * rotation->r20);
        out.y =
                (ori.x * rotation->r01 + ori.y * rotation->r11 +
                 ori.z * rotation->r21);
        out.z =
                (ori.x * rotation->r02 + ori.y * rotation->r12 +
                 ori.z * rotation->r22);
        return out;
}

struct mli_Vec mli_transform_position(
        const struct mli_Mat *rotation,
        const struct mli_Vec translation,
        const struct mli_Vec pos)
{
        struct mli_Vec out;
        out.x = (pos.x * rotation->r00 + pos.y * rotation->r01 +
                 pos.z * rotation->r02) +
                translation.x;
        out.y = (pos.x * rotation->r10 + pos.y * rotation->r11 +
                 pos.z * rotation->r12) +
                translation.y;
        out.z = (pos.x * rotation->r20 + pos.y * rotation->r21 +
                 pos.z * rotation->r22) +
                translation.z;
        return out;
}

struct mli_Vec mli_transform_position_inverse(
        const struct mli_Mat *rotation,
        const struct mli_Vec translation,
        const struct mli_Vec pos)
{
        struct mli_Vec out;
        out.x = (pos.x * rotation->r00 + pos.y * rotation->r10 +
                 pos.z * rotation->r20) -
                (rotation->r00 * translation.x + rotation->r10 * translation.y +
                 rotation->r20 * translation.z);
        out.y = (pos.x * rotation->r01 + pos.y * rotation->r11 +
                 pos.z * rotation->r21) -
                (rotation->r01 * translation.x + rotation->r11 * translation.y +
                 rotation->r21 * translation.z);
        out.z = (pos.x * rotation->r02 + pos.y * rotation->r12 +
                 pos.z * rotation->r22) -
                (rotation->r02 * translation.x + rotation->r12 * translation.y +
                 rotation->r22 * translation.z);
        return out;
}

struct mli_Ray mli_transform_ray(
        const struct mli_Mat *rotation,
        const struct mli_Vec translation,
        const struct mli_Ray in)
{
        struct mli_Ray out;
        out.support = mli_transform_position(rotation, translation, in.support);
        out.direction = mli_transform_orientation(rotation, in.direction);
        return out;
}

struct mli_Ray mli_transform_ray_inverse(
        const struct mli_Mat *rotation,
        const struct mli_Vec translation,
        const struct mli_Ray in)
{
        struct mli_Ray out;
        out.support = mli_transform_position_inverse(
                rotation, translation, in.support);
        out.direction =
                mli_transform_orientation_inverse(rotation, in.direction);
        return out;
}

struct mli_Ray mli_HomTraComp_ray(
        const struct mli_HomTra *t,
        const struct mli_Ray in)
{
        return mli_transform_ray(&t->rotation, t->translation, in);
}

struct mli_Ray mli_HomTraComp_ray_inverse(
        const struct mli_HomTra *t,
        const struct mli_Ray in)
{
        return mli_transform_ray_inverse(&t->rotation, t->translation, in);
}

struct mli_Vec mli_HomTraComp_pos(
        const struct mli_HomTra *t,
        const struct mli_Vec in)
{
        return mli_transform_position(&t->rotation, t->translation, in);
}

struct mli_Vec mli_HomTraComp_pos_inverse(
        const struct mli_HomTra *t,
        const struct mli_Vec in)
{
        return mli_transform_position_inverse(&t->rotation, t->translation, in);
}

struct mli_Vec mli_HomTraComp_dir(
        const struct mli_HomTra *t,
        const struct mli_Vec in)
{
        return mli_transform_orientation(&t->rotation, in);
}

struct mli_Vec mli_HomTraComp_dir_inverse(
        const struct mli_HomTra *t,
        const struct mli_Vec in)
{
        return mli_transform_orientation_inverse(&t->rotation, in);
}

int mli_HomTraComp_equal(
        const struct mli_HomTraComp a,
        const struct mli_HomTraComp b)
{
        if (!mli_Vec_equal(a.translation, b.translation))
                return 0;
        if (!mli_Quaternion_equal(a.rotation, b.rotation))
                return 0;
        return 1;
}

struct mli_HomTraComp mli_HomTraComp_sequence(
        const struct mli_HomTraComp a,
        const struct mli_HomTraComp b)
{
        struct mli_HomTra b_;
        struct mli_HomTraComp s;
        b_ = mli_HomTraComp_from_compact(b);
        s.translation = mli_HomTraComp_pos(&b_, a.translation);
        s.rotation = mli_Quaternion_product(b.rotation, a.rotation);
        return s;
}

/* image */
/* ----- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

#define MLI_IMAGE_CHUNK_EDGE_SIZE (2u * 3u * 5u)

struct mli_Image mli_Image_init(void)
{
        struct mli_Image out;
        out.geometry = mli_image_ChunkGeometry_set(0, 0, 1);
        out.chunks = NULL;
        return out;
}

void mli_Image_free(struct mli_Image *self)
{
        uint64_t crow, ccol;
        for (crow = 0; crow < self->geometry.num_chunks_row; crow++) {
                for (ccol = 0; ccol < self->geometry.num_chunks_col; ccol++) {
                        mli_image_Chunk_free(&self->chunks[crow][ccol]);
                }
                free(self->chunks[crow]);
        }
        free(self->chunks);
        (*self) = mli_Image_init();
}

int mli_Image__malloc(
        struct mli_Image *self,
        const uint32_t num_cols,
        const uint32_t num_rows)
{
        uint64_t crow, ccol;
        mli_Image_free(self);

        self->geometry = mli_image_ChunkGeometry_set(
                num_cols, num_rows, MLI_IMAGE_CHUNK_EDGE_SIZE);

        chk_malloc(
                self->chunks,
                struct mli_image_Chunk *,
                self->geometry.num_chunks_row);

        for (crow = 0; crow < self->geometry.num_chunks_row; crow++) {
                chk_malloc(
                        self->chunks[crow],
                        struct mli_image_Chunk,
                        self->geometry.num_chunks_col);
                for (ccol = 0; ccol < self->geometry.num_chunks_col; ccol++) {
                        chk(mli_image_Chunk_malloc(
                                &self->chunks[crow][ccol],
                                self->geometry.chunk_edge_size));
                }
        }
        return 1;
chk_error:
        mli_Image_free(self);
        return 0;
}

int mli_Image_malloc(
        struct mli_Image *self,
        const uint32_t num_cols,
        const uint32_t num_rows)
{
        if (self->geometry.num_cols == num_cols &&
            self->geometry.num_rows == num_rows) {
                return 1;
        } else {
                return mli_Image__malloc(self, num_cols, num_rows);
        }
}

int mli_Image_malloc_same_size(
        struct mli_Image *self,
        const struct mli_Image *other)
{
        return mli_Image_malloc(
                self, mli_Image_num_cols(other), mli_Image_num_rows(other));
}

uint64_t mli_Image_num_pixel(const struct mli_Image *self)
{
        return self->geometry.num_cols * self->geometry.num_rows;
}

uint64_t mli_Image_num_cols(const struct mli_Image *self)
{
        return self->geometry.num_cols;
}

uint64_t mli_Image_num_rows(const struct mli_Image *self)
{
        return self->geometry.num_rows;
}

/* SET GET PIXELWALK */

void mli_Image__set_by_PixelWalk(
        const struct mli_Image *self,
        const struct mli_image_PixelWalk walk,
        const struct mli_Color color)
{
        mli_image_Chunk_set(
                &self->chunks[walk.chunk_row][walk.chunk_col],
                walk.sub_col,
                walk.sub_row,
                color);
}

struct mli_Color mli_Image__get_by_PixelWalk(
        const struct mli_Image *self,
        const struct mli_image_PixelWalk walk)
{
        return mli_image_Chunk_get(
                &self->chunks[walk.chunk_row][walk.chunk_col],
                walk.sub_col,
                walk.sub_row);
}

struct mli_Color *mli_Image__get_ptr_by_PixelWalk(
        const struct mli_Image *self,
        const struct mli_image_PixelWalk walk)
{
        return mli_image_Chunk_get_ptr(
                &self->chunks[walk.chunk_row][walk.chunk_col],
                walk.sub_col,
                walk.sub_row);
}

/* SET GET PIXEL */

struct mli_Color mli_Image_get_by_Pixel(
        const struct mli_Image *self,
        const struct mli_image_Pixel px)
{
        struct mli_image_PixelWalk w =
                mli_image_PixelWalk_from_pixel(&self->geometry, px);
        return mli_Image__get_by_PixelWalk(self, w);
}

void mli_Image_set_by_Pixel(
        struct mli_Image *self,
        const struct mli_image_Pixel px,
        const struct mli_Color color)
{
        struct mli_image_PixelWalk w =
                mli_image_PixelWalk_from_pixel(&self->geometry, px);
        mli_Image__set_by_PixelWalk(self, w, color);
}

/* SET GET ROW COL */

void mli_Image_set_by_col_row(
        struct mli_Image *self,
        const uint32_t col,
        const uint32_t row,
        const struct mli_Color color)
{
        struct mli_image_Pixel px;
        px.col = col;
        px.row = row;
        mli_Image_set_by_Pixel(self, px, color);
}

struct mli_Color mli_Image_get_by_col_row(
        const struct mli_Image *self,
        const uint32_t col,
        const uint32_t row)
{
        struct mli_image_Pixel px;
        px.col = col;
        px.row = row;
        return mli_Image_get_by_Pixel(self, px);
}

void mli_Image_set_all(
        const struct mli_Image *self,
        const struct mli_Color color)
{
        uint64_t i;
        const uint64_t NUM = mli_Image_num_pixel(self);
        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();
        for (i = 0; i < NUM; i++) {
                mli_Image__set_by_PixelWalk(self, w, color);
                mli_image_PixelWalk_walk(&w, &self->geometry);
        }
}

int mli_Image_scale_down_twice(
        const struct mli_Image *src,
        struct mli_Image *dst)
{
        uint64_t sr, sc, i;
        uint64_t DST_NUM = mli_Image_num_pixel(dst);
        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();

        chk_msg(mli_Image_num_cols(src) % 2 == 0,
                "Expected num_cols in src image to be exact multiple of 2.");
        chk_msg(mli_Image_num_rows(src) % 2 == 0,
                "Expected num_rows in src image to be exact multiple of 2.");
        chk(mli_Image_malloc(
                dst, mli_Image_num_cols(src) / 2, mli_Image_num_rows(src) / 2));

        for (i = 0; i < DST_NUM; i++) {
                struct mli_Color mix[4];
                const struct mli_image_Pixel px =
                        mli_image_PixelWalk_get_Pixel(&w, &dst->geometry);
                sr = px.row * 2u;
                sc = px.col * 2u;

                mix[0] = mli_Image_get_by_col_row(src, sc + 0, sr + 0);
                mix[1] = mli_Image_get_by_col_row(src, sc + 0, sr + 1);
                mix[2] = mli_Image_get_by_col_row(src, sc + 1, sr + 0);
                mix[3] = mli_Image_get_by_col_row(src, sc + 1, sr + 1);

                mli_Image__set_by_PixelWalk(dst, w, mli_Color_mean(mix, 4));
                mli_image_PixelWalk_walk(&w, &dst->geometry);
        }

        return 1;
chk_error:
        return 0;
}

int mli_Image_sobel(const struct mli_Image *src, struct mli_Image *dst)
{
        uint64_t i;
        uint64_t NUM = mli_Image_num_pixel(src);

        uint16_t row, col;
        struct mli_image_Pixel c_0_r_0;
        struct mli_image_Pixel cm1_rp1;
        struct mli_image_Pixel cm1_r_0;
        struct mli_image_Pixel cm1_rm1;
        struct mli_image_Pixel cp1_rp1;
        struct mli_image_Pixel cp1_r_0;
        struct mli_image_Pixel cp1_rm1;
        struct mli_image_Pixel c_0_rp1;
        struct mli_image_Pixel c_0_rm1;

        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();

        struct mli_Color x, y, color_result;

        struct mli_Color C_cm1_rp1;
        struct mli_Color C_cm1_r_0;
        struct mli_Color C_cm1_rm1;

        struct mli_Color C_c_0_rp1;
        /*               C_c_0_r_0 */
        struct mli_Color C_c_0_rm1;

        struct mli_Color C_cp1_rp1;
        struct mli_Color C_cp1_r_0;
        struct mli_Color C_cp1_rm1;

        /*
         *       relative coordinates of folding kernel
         *       'c' is column and 'r' is row.
         *       ______ ______ ______
         *      | c:-1 | c: 0 | c:+1 |
         *      | r:+1 | r:+1 | r:+1 |
         *      |______|______|______|
         *      | c:-1 | c: 0 | c:+1 |
         *      | r: 0 | r: 0 | r: 0 |
         *      |______|______|______|
         *      | c:-1 | c: 0 | c:+1 |
         *      | r:-1 | r:-1 | r:-1 |
         *      |______|______|______|
         */

        /*
         *       weights for x-gradient-kernel
         *       ______ ______ ______
         *      | w:-1 | w: 0 | w:+1 |
         *      |      |      |      |
         *      |______|______|______|
         *      | w:-2 | w: 0 | w:+2 |
         *      |      |      |      |
         *      |______|______|______|
         *      | w:-1 | w: 0 | w:+1 |
         *      |      |      |      |
         *      |______|______|______|
         */

        /*
         *       weights for y-gradient-kernel
         *       ______ ______ ______
         *      | w:+1 | w:+2 | w:+1 |
         *      |      |      |      |
         *      |______|______|______|
         *      | w: 0 | w: 0 | w: 0 |
         *      |      |      |      |
         *      |______|______|______|
         *      | w:-1 | w:-2 | w:-1 |
         *      |      |      |      |
         *      |______|______|______|
         */

        chk(mli_Image_malloc_same_size(dst, src));
        for (i = 0; i < NUM; i++) {
                int is_on_edge = 0;
                c_0_r_0 = mli_image_PixelWalk_get_Pixel(&w, &src->geometry);
                row = c_0_r_0.row;
                col = c_0_r_0.col;
                if (row == 0 || row == mli_Image_num_rows(src) - 1)
                        is_on_edge = 1;
                if (col == 0 || col == mli_Image_num_cols(src) - 1)
                        is_on_edge = 1;

                if (is_on_edge) {
                        color_result = mli_Color_set(0.0, 0.0, 0.0);
                } else {
                        cm1_rp1 = mli_image_Pixel_set_col_row(col - 1, row + 1);
                        cm1_r_0 = mli_image_Pixel_set_col_row(col - 1, row + 0);
                        cm1_rm1 = mli_image_Pixel_set_col_row(col - 1, row - 1);

                        c_0_rp1 = mli_image_Pixel_set_col_row(col, row + 1);
                        /* c_0_r_0 */
                        c_0_rm1 = mli_image_Pixel_set_col_row(col, row - 1);

                        cp1_rp1 = mli_image_Pixel_set_col_row(col + 1, row + 1);
                        cp1_r_0 = mli_image_Pixel_set_col_row(col + 1, row + 0);
                        cp1_rm1 = mli_image_Pixel_set_col_row(col + 1, row - 1);

                        x = mli_Color_set(0., 0., 0.);
                        y = mli_Color_set(0., 0., 0.);

                        C_cm1_rp1 = mli_Image_get_by_Pixel(src, cm1_rp1);
                        C_cm1_r_0 = mli_Image_get_by_Pixel(src, cm1_r_0);
                        C_cm1_rm1 = mli_Image_get_by_Pixel(src, cm1_rm1);

                        C_c_0_rp1 = mli_Image_get_by_Pixel(src, c_0_rp1);

                        C_c_0_rm1 = mli_Image_get_by_Pixel(src, c_0_rm1);

                        C_cp1_rp1 = mli_Image_get_by_Pixel(src, cp1_rp1);
                        C_cp1_r_0 = mli_Image_get_by_Pixel(src, cp1_r_0);
                        C_cp1_rm1 = mli_Image_get_by_Pixel(src, cp1_rm1);

                        x = mli_Color_add(
                                x, mli_Color_multiply(C_cm1_rp1, -1.0));
                        x = mli_Color_add(
                                x, mli_Color_multiply(C_cm1_r_0, -2.0));
                        x = mli_Color_add(
                                x, mli_Color_multiply(C_cm1_rm1, -1.0));

                        x = mli_Color_add(
                                x, mli_Color_multiply(C_cp1_rp1, +1.0));
                        x = mli_Color_add(
                                x, mli_Color_multiply(C_cp1_r_0, +2.0));
                        x = mli_Color_add(
                                x, mli_Color_multiply(C_cp1_rm1, +1.0));

                        y = mli_Color_add(
                                y, mli_Color_multiply(C_cm1_rp1, -1.0));
                        y = mli_Color_add(
                                y, mli_Color_multiply(C_c_0_rp1, -2.0));
                        y = mli_Color_add(
                                y, mli_Color_multiply(C_cp1_rp1, -1.0));

                        y = mli_Color_add(
                                y, mli_Color_multiply(C_cm1_rm1, +1.0));
                        y = mli_Color_add(
                                y, mli_Color_multiply(C_c_0_rm1, +2.0));
                        y = mli_Color_add(
                                y, mli_Color_multiply(C_cp1_rm1, +1.0));

                        color_result = mli_Color_set(
                                mli_math_hypot(x.r, y.r),
                                mli_math_hypot(x.g, y.g),
                                mli_math_hypot(x.b, y.b));
                }
                mli_Image__set_by_PixelWalk(dst, w, color_result);
                mli_image_PixelWalk_walk(&w, &src->geometry);
        }
        return 1;
chk_error:
        return 0;
}

int mli_Image_luminance_threshold_dilatation(
        const struct mli_Image *self,
        const float threshold,
        const struct mli_Color marker,
        struct mli_Image *out)
{
        uint64_t i;
        const uint64_t NUM = mli_Image_num_pixel(self);
        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();
        const int32_t num_rows = mli_Image_num_rows(self);
        const int32_t num_cols = mli_Image_num_cols(self);

        chk(mli_Image_malloc_same_size(out, self));

        for (i = 0u; i < NUM; i++) {
                const struct mli_Color c = mli_Image__get_by_PixelWalk(self, w);
                const float luminance = mli_Color_luminance(c);
                struct mli_image_Pixel px =
                        mli_image_PixelWalk_get_Pixel(&w, &self->geometry);
                if (luminance > threshold) {
                        int32_t orow, ocol;
                        for (orow = -1; orow < 2; orow++) {
                                for (ocol = -1; ocol < 2; ocol++) {
                                        if (px.row + orow >= 0 &&
                                            px.col + ocol >= 0 &&
                                            px.row + orow < num_rows &&
                                            px.col + ocol < num_cols) {
                                                mli_Image_set_by_col_row(
                                                        out,
                                                        px.col + ocol,
                                                        px.row + orow,
                                                        marker);
                                        }
                                }
                        }
                }
        }
        return 1;
chk_error:
        return 0;
}

void mli_image_PixelVector_push_back_all_from_image(
        struct mli_image_PixelVector *pixels,
        const struct mli_Image *image)
{
        uint64_t i;
        const uint64_t NUM = mli_Image_num_pixel(image);
        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();

        pixels->size = 0;
        for (i = 0u; i < NUM; i++) {
                struct mli_image_Pixel px =
                        mli_image_PixelWalk_get_Pixel(&w, &image->geometry);
                mli_image_PixelVector_push_back(pixels, px);
                mli_image_PixelWalk_walk(&w, &image->geometry);
        }
}

int mli_image_PixelVector_above_threshold(
        const struct mli_Image *self,
        const float threshold,
        struct mli_image_PixelVector *pixels)
{
        uint64_t i;
        const uint64_t NUM = mli_Image_num_pixel(self);
        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();

        pixels->size = 0;
        for (i = 0u; i < NUM; i++) {
                struct mli_image_Pixel px =
                        mli_image_PixelWalk_get_Pixel(&w, &self->geometry);
                double lum = 0.0;
                struct mli_Color c = mli_Image__get_by_PixelWalk(self, w);
                lum = c.r + c.g + c.b;
                if (lum > threshold) {
                        chk(mli_image_PixelVector_push_back(pixels, px));
                }
                mli_image_PixelWalk_walk(&w, &self->geometry);
        }
        return 1;
chk_error:
        return 0;
}

int mli_Image_copy(const struct mli_Image *src, struct mli_Image *dst)
{
        uint64_t i;
        const uint64_t NUM = mli_Image_num_pixel(src);
        struct mli_image_PixelWalk walk = mli_image_PixelWalk_init();
        struct mli_Color color;

        chk(mli_Image_malloc_same_size(dst, src));

        for (i = 0u; i < NUM; i++) {
                color = mli_Image__get_by_PixelWalk(src, walk);
                mli_Image__set_by_PixelWalk(dst, walk, color);
                mli_image_PixelWalk_walk(&walk, &src->geometry);
        }

        return 1;
chk_error:
        return 0;
}

int mli_Image_fabs_difference(
        const struct mli_Image *a,
        const struct mli_Image *b,
        struct mli_Image *out)
{
        uint64_t i;
        uint64_t NUM = mli_Image_num_pixel(a);
        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();

        chk_msg(mli_image_ChunkGeometry_equal(a->geometry, b->geometry),
                "Expected images 'a' and 'b' to have the same geometry.");
        chk(mli_Image_malloc_same_size(out, a));

        for (i = 0; i < NUM; i++) {
                const struct mli_Color _a = mli_Image__get_by_PixelWalk(a, w);
                const struct mli_Color _b = mli_Image__get_by_PixelWalk(b, w);
                struct mli_Color _o;
                _o.r = fabs(_a.r - _b.r);
                _o.g = fabs(_a.g - _b.g);
                _o.b = fabs(_a.b - _b.b);
                mli_Image__set_by_PixelWalk(out, w, _o);
                mli_image_PixelWalk_walk(&w, &out->geometry);
        }
        return 1;
chk_error:
        return 0;
}

void mli_Image_histogram(
        struct mli_Image *self,
        const double *col_bin_edges,
        const double *row_bin_edges,
        const double col_val,
        const double row_val,
        const struct mli_Color weight)
{
        uint64_t col_upper_idx, row_upper_idx;
        int valid_col, valid_row;

        col_upper_idx = mli_math_upper_compare_double(
                col_bin_edges, mli_Image_num_cols(self) + 1, col_val);

        row_upper_idx = mli_math_upper_compare_double(
                row_bin_edges, mli_Image_num_rows(self) + 1, row_val);

        valid_col = col_upper_idx > 0 &&
                    col_upper_idx < mli_Image_num_cols(self) + 1;
        valid_row = row_upper_idx > 0 &&
                    row_upper_idx < mli_Image_num_rows(self) + 1;

        if (valid_col && valid_row) {
                const uint32_t col = col_upper_idx - 1;
                const uint32_t row = row_upper_idx - 1;
                const struct mli_image_PixelWalk w =
                        mli_image_PixelWalk_from_pixel(
                                &self->geometry,
                                mli_image_Pixel_set_col_row(col, row));
                struct mli_Color *c = mli_Image__get_ptr_by_PixelWalk(self, w);
                (*c) = mli_Color_add((*c), weight);
        }
}

struct mli_Color mli_Image_max(const struct mli_Image *self)
{
        struct mli_Color max = mli_Color_set(-FLT_MAX, -FLT_MAX, -FLT_MAX);
        uint64_t i;
        uint64_t NUM = mli_Image_num_pixel(self);
        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();

        for (i = 0; i < NUM; i++) {
                struct mli_Color c = mli_Image__get_by_PixelWalk(self, w);
                if (c.r > max.r) {
                        max.r = c.r;
                }
                if (c.g > max.g) {
                        max.g = c.g;
                }
                if (c.b > max.b) {
                        max.b = c.b;
                }
                mli_image_PixelWalk_walk(&w, &self->geometry);
        }
        return max;
}

void mli_Image_multiply(struct mli_Image *self, const struct mli_Color color)
{
        uint64_t i;
        uint64_t NUM = mli_Image_num_pixel(self);
        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();

        for (i = 0; i < NUM; i++) {
                struct mli_Color *_self =
                        mli_Image__get_ptr_by_PixelWalk(self, w);
                _self->r *= color.r;
                _self->g *= color.g;
                _self->b *= color.b;
                mli_image_PixelWalk_walk(&w, &self->geometry);
        }
}

void mli_Image_power(struct mli_Image *self, const struct mli_Color power)
{
        uint64_t i;
        uint64_t NUM = mli_Image_num_pixel(self);
        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();

        for (i = 0; i < NUM; i++) {
                struct mli_Color *_self =
                        mli_Image__get_ptr_by_PixelWalk(self, w);
                _self->r = pow(_self->r, power.r);
                _self->g = pow(_self->g, power.g);
                _self->b = pow(_self->b, power.b);
                mli_image_PixelWalk_walk(&w, &self->geometry);
        }
}

int mli_Image_divide_pixelwise(
        const struct mli_Image *numerator,
        const struct mli_Image *denominator,
        struct mli_Image *out)
{
        uint64_t i;
        uint64_t NUM = mli_Image_num_pixel(out);
        struct mli_image_PixelWalk w = mli_image_PixelWalk_init();

        chk_msg(mli_image_ChunkGeometry_equal(
                        numerator->geometry, denominator->geometry),
                "Expected images 'numerator' and 'denominator' to have "
                "the same geometry.");
        chk(mli_Image_malloc_same_size(out, numerator));

        for (i = 0; i < NUM; i++) {
                const struct mli_Color _numerator =
                        mli_Image__get_by_PixelWalk(numerator, w);
                const struct mli_Color _denominator =
                        mli_Image__get_by_PixelWalk(denominator, w);
                struct mli_Color _out = mli_Color_set(
                        _numerator.r / _denominator.r,
                        _numerator.g / _denominator.g,
                        _numerator.b / _denominator.b);
                mli_Image__set_by_PixelWalk(out, w, _out);
                mli_image_PixelWalk_walk(&w, &out->geometry);
        }
        return 1;
chk_error:
        return 0;
}

/* image_Pixel */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_image_Pixel mli_image_Pixel_set_col_row(
        const uint16_t col,
        const uint16_t row)
{
        struct mli_image_Pixel out;
        out.col = col;
        out.row = row;
        return out;
}

void mli_image_Pixel_fprint(FILE *f, const struct mli_image_Pixel *self)
{
        fprintf(f, "(col: %d, row: %d)", self->col, self->row);
}

/* image_PixelVector */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_VECTOR_IMPLEMENTATION(mli_image_PixelVector, struct mli_image_Pixel)

/* image_PixelWalk */
/* --------------- */

/* Copyright 2020-2021 Sebastian Achim Mueller */

struct mli_image_PixelWalk mli_image_PixelWalk_init(void)
{
        struct mli_image_PixelWalk walk;
        walk.i = 0;

        walk.chunk_row = 0u;
        walk.sub_row = 0u;

        walk.chunk_col = 0u;
        walk.sub_col = 0u;
        return walk;
}

struct mli_image_Pixel mli_image_PixelWalk_get_Pixel(
        const struct mli_image_PixelWalk *self,
        const struct mli_image_ChunkGeometry *chunk_geometry)
{
        struct mli_image_Pixel px;
        const struct mli_image_ChunkGeometry *cg = chunk_geometry;
        px.row = self->chunk_row * cg->chunk_edge_size + self->sub_row;
        px.col = self->chunk_col * cg->chunk_edge_size + self->sub_col;
        return px;
}

void mli_image_PixelWalk_walk(
        struct mli_image_PixelWalk *self,
        const struct mli_image_ChunkGeometry *chunk_geometry)
{
        struct mli_image_Pixel px;
        const struct mli_image_ChunkGeometry *cg = chunk_geometry;

        self->sub_row += 1u;
        px = mli_image_PixelWalk_get_Pixel(self, chunk_geometry);
        if (self->sub_row < cg->chunk_edge_size && px.row < cg->num_rows) {
                return;
        }

        self->sub_row = 0u;
        self->sub_col += 1u;
        px = mli_image_PixelWalk_get_Pixel(self, chunk_geometry);
        if (self->sub_col < cg->chunk_edge_size && px.col < cg->num_cols) {
                return;
        }

        self->sub_col = 0u;
        self->chunk_row += 1u;
        if (self->chunk_row < cg->num_chunks_row) {
                return;
        }

        self->chunk_row = 0u;
        self->chunk_col += 1u;
        if (self->chunk_col < cg->num_chunks_col) {
                return;
        }

        self->chunk_col = 0;
        return;
}

struct mli_image_PixelWalk mli_image_PixelWalk_from_pixel(
        const struct mli_image_ChunkGeometry *geometry,
        struct mli_image_Pixel pixel)
{
        struct mli_image_PixelWalk address;
        address.chunk_col = pixel.col / geometry->chunk_edge_size;
        address.sub_col = pixel.col % geometry->chunk_edge_size;
        address.chunk_row = pixel.row / geometry->chunk_edge_size;
        address.sub_row = pixel.row % geometry->chunk_edge_size;
        return address;
}

void mli_image_PixelWalk_fprint(FILE *f, const struct mli_image_PixelWalk *self)
{
        fprintf(f,
                "(col: %d/%d, row: %d/%d)",
                self->chunk_col,
                self->sub_col,
                self->chunk_row,
                self->sub_row);
}

/* image_chunk */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_image_Chunk mli_image_Chunk_init(void)
{
        struct mli_image_Chunk out;
        out.edge_size = 0;
        out.array = NULL;
        return out;
}

void mli_image_Chunk_free(struct mli_image_Chunk *self)
{
        free(self->array);
        (*self) = mli_image_Chunk_init();
}

int mli_image_Chunk_malloc(
        struct mli_image_Chunk *self,
        const uint64_t edge_size)
{
        const uint64_t num_pixel = edge_size * edge_size;
        self->edge_size = edge_size;
        chk_malloc(self->array, struct mli_Color, num_pixel);
        return 1;
chk_error:
        return 0;
}

void mli_image_Chunk_set(
        struct mli_image_Chunk *self,
        const uint64_t col,
        const uint64_t row,
        const struct mli_Color color)
{
        self->array[mli_image_Chunk__idx(self, col, row)] = color;
}

struct mli_Color mli_image_Chunk_get(
        const struct mli_image_Chunk *self,
        const uint64_t col,
        const uint64_t row)
{
        return self->array[mli_image_Chunk__idx(self, col, row)];
}

struct mli_Color *mli_image_Chunk_get_ptr(
        const struct mli_image_Chunk *self,
        const uint64_t col,
        const uint64_t row)
{
        return &self->array[mli_image_Chunk__idx(self, col, row)];
}

uint64_t mli_image_Chunk__idx(
        const struct mli_image_Chunk *self,
        const uint64_t col,
        const uint64_t row)
{
        return col * self->edge_size + row;
}

struct mli_image_ChunkGeometry mli_image_ChunkGeometry_set(
        const uint64_t num_cols,
        const uint64_t num_rows,
        const uint64_t chunk_edge_size)
{
        struct mli_image_ChunkGeometry out;
        uint64_t full_row, full_col, rest_row, rest_col;

        assert(chunk_edge_size > 0);

        out.num_cols = num_cols;
        out.num_rows = num_rows;
        out.chunk_edge_size = chunk_edge_size;

        full_row = num_rows / chunk_edge_size;
        full_col = num_cols / chunk_edge_size;
        rest_row = num_rows % chunk_edge_size;
        rest_col = num_cols % chunk_edge_size;

        out.num_chunks_row = rest_row ? full_row + 1 : full_row;
        out.num_chunks_col = rest_col ? full_col + 1 : full_col;
        return out;
}

int mli_image_ChunkGeometry_equal(
        const struct mli_image_ChunkGeometry a,
        const struct mli_image_ChunkGeometry b)
{
        if (a.num_cols != b.num_cols)
                return 0;
        if (a.num_rows != b.num_rows)
                return 0;
        if (a.chunk_edge_size != b.chunk_edge_size)
                return 0;
        if (a.num_chunks_row != b.num_chunks_row)
                return 0;
        if (a.num_chunks_col != b.num_chunks_col)
                return 0;
        return 1;
}

/* image_ppm */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Image_from_io(struct mli_Image *img, struct mli_IO *f)
{
        struct mli_String line = mli_String_init();
        uint64_t num_comment_lines = 0;
        uint64_t num_cols;
        uint64_t num_rows;
        uint64_t col;
        uint64_t row;
        chk_msg(mli_IO_text_read_line(f, &line, '\n'),
                "Can't read header-line.");
        chk_msg(mli_String_equal_cstr(&line, "P6"),
                "Expected starts with 'P6'.");

        while (1) {
                chk_msg(num_comment_lines < 1024,
                        "Expected < 1024 comment lines.");
                chk_msg(mli_IO_text_read_line(f, &line, '\n'),
                        "Can't read header-line.");
                if (mli_String_starts_with_cstr(&line, "#")) {
                        num_comment_lines += 1u;
                } else {
                        break;
                }
        }
        chk_msg(mli_String_to_uint64(&num_cols, &line, 10),
                "Can't parse num. columns.");
        chk_msg(mli_IO_text_read_line(f, &line, '\n'),
                "Can't read header-line.");
        chk_msg(mli_String_to_uint64(&num_rows, &line, 10),
                "Can't parse num. rows.");
        chk_msg(mli_IO_text_read_line(f, &line, '\n'),
                "Can't read header-line.");
        chk_msg(mli_String_equal_cstr(&line, "255"),
                "Expected 8bit range '255'.");

        chk_mem(mli_Image_malloc(img, num_cols, num_rows));
        for (row = 0; row < mli_Image_num_rows(img); row++) {
                for (col = 0; col < mli_Image_num_cols(img); col++) {
                        uint8_t r, g, b;
                        struct mli_Color color;
                        chk_IO_read(&r, sizeof(uint8_t), 1u, f);
                        chk_IO_read(&g, sizeof(uint8_t), 1u, f);
                        chk_IO_read(&b, sizeof(uint8_t), 1u, f);
                        color.r = (float)r;
                        color.g = (float)g;
                        color.b = (float)b;
                        mli_Image_set_by_col_row(img, col, row, color);
                }
        }
        mli_String_free(&line);
        return 1;
chk_error:
        mli_String_free(&line);
        mli_Image_free(img);
        return 0;
}

int mli_Image_to_io(const struct mli_Image *img, struct mli_IO *f)
{

        uint32_t col;
        uint32_t row;
        chk(mli_IO_text_write_cstr(f, "P6\n"));
        chk(mli_IO_text_write_cstr(f, "# merlict_c89\n"));
        chk(mli_IO_text_write_cstr_format(
                f,
                "# MLI_VERSION %d.%d.%d\n",
                MLI_VERSION_MAYOR,
                MLI_VERSION_MINOR,
                MLI_VERSION_PATCH));
        chk(mli_IO_text_write_cstr_format(f, "%ld\n", mli_Image_num_cols(img)));
        chk(mli_IO_text_write_cstr_format(f, "%ld\n", mli_Image_num_rows(img)));
        chk(mli_IO_text_write_cstr(f, "255\n"));
        for (row = 0; row < mli_Image_num_rows(img); row++) {
                for (col = 0; col < mli_Image_num_cols(img); col++) {
                        struct mli_Color color =
                                mli_Image_get_by_col_row(img, col, row);
                        struct mli_Color out =
                                mli_Color_truncate(color, 0., 255.);
                        uint8_t r = (uint8_t)out.r;
                        uint8_t g = (uint8_t)out.g;
                        uint8_t b = (uint8_t)out.b;
                        chk_IO_write(&r, sizeof(uint8_t), 1u, f);
                        chk_IO_write(&g, sizeof(uint8_t), 1u, f);
                        chk_IO_write(&b, sizeof(uint8_t), 1u, f);
                }
        }
        return 1;
chk_error:
        return 0;
}

/* image_print */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mli_Image_print(const struct mli_Image *img, const uint64_t print_mode)
{
        const uint64_t num_symbols = 0;
        mli_Image_print_chars(img, NULL, NULL, NULL, num_symbols, print_mode);
}

void mli_Image_print_chars(
        const struct mli_Image *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols,
        const uint64_t print_mode)
{
        if (print_mode == MLI_ANSI_ESCAPE_COLOR) {
                mli_Image_print_ansi_escape_chars(
                        img, symbols, rows, cols, num_symbols);
        } else {
                mli_Image_print_ascii_chars(
                        img, symbols, rows, cols, num_symbols);
                if (print_mode != MLI_ASCII_MONOCHROME) {
                        fprintf(stderr,
                                "Do not know print_mode %u\n",
                                (uint32_t)print_mode);
                }
        }
        return;
}

void mli_Image_print_ansi_escape_chars(
        const struct mli_Image *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols)
{
        uint32_t col, row, sym;
        char symbol;
        for (row = 0; row < mli_Image_num_rows(img); row++) {
                for (col = 0; col < mli_Image_num_cols(img); col++) {
                        struct mli_Color color =
                                mli_Image_get_by_col_row(img, col, row);
                        struct mli_Color out =
                                mli_Color_truncate(color, 0., 255.);
                        uint8_t r = (uint8_t)out.r;
                        uint8_t g = (uint8_t)out.g;
                        uint8_t b = (uint8_t)out.b;
                        symbol = ' ';
                        for (sym = 0; sym < num_symbols; sym++) {
                                if (rows[sym] == row && cols[sym] == col) {
                                        symbol = symbols[sym];
                                        break;
                                }
                        }
                        printf("\033[48;2;%u;%u;%um%c\033[0m", r, g, b, symbol);
                }
                putchar('\n');
        }
        fflush(stdout);
}

void mli_Image_print_ascii_chars(
        const struct mli_Image *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols)
{
        uint32_t col, row, sym;
        char symbol;
        char chars_with_ascending_fill[16] = {
                ' ',
                '.',
                ',',
                '-',
                '~',
                '^',
                '+',
                '=',
                'i',
                'j',
                '$',
                'l',
                'L',
                'Y',
                '8',
                '#'};
        for (row = 0; row < mli_Image_num_rows(img); row++) {
                for (col = 0; col < mli_Image_num_cols(img); col++) {
                        struct mli_Color color =
                                mli_Image_get_by_col_row(img, col, row);
                        struct mli_Color out =
                                mli_Color_truncate(color, 0., 255.);
                        float lum = 1.0 / 3.0 * (out.r + out.g + out.b);
                        int64_t l = lum / 16.0;
                        if (l < 0) {
                                l = 0;
                        }
                        if (l >= 16) {
                                l = 15;
                        }
                        symbol = chars_with_ascending_fill[l];
                        for (sym = 0; sym < num_symbols; sym++) {
                                if (rows[sym] == row && cols[sym] == col) {
                                        symbol = symbols[sym];
                                        break;
                                }
                        }
                        putchar(symbol);
                }
                putchar('\n');
        }
        fflush(stdout);
}

/* intersection */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Intersection mli_Intersection_init(void)
{
        struct mli_Intersection psec;
        psec.geometry_id = mli_GeometryId_init();
        psec.position_local = mli_Vec_init(0.0, 0.0, 0.0);
        psec.distance_of_ray = DBL_MAX;
        return psec;
}

/* intersection_and_scenery */
/* ------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_IntersectionLayerSide mli_IntersectionLayerSide_init(void)
{
        struct mli_IntersectionLayerSide side;
        side.surface = NULL;
        side.surface_idx = 0;
        side.medium = NULL;
        side.medium_idx = 0;
        return side;
}

struct mli_IntersectionLayer mli_IntersectionLayer_init(void)
{
        struct mli_IntersectionLayer ilay;
        ilay.side_coming_from = mli_IntersectionLayerSide_init();
        ilay.side_going_to = mli_IntersectionLayerSide_init();
        return ilay;
}

struct mli_IntersectionLayer mli_raytracing_get_intersection_layer(
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec)
{
        const uint64_t idx = mli_Scenery_resolve_boundary_layer_idx(
                scenery, isec->geometry_id);
        const struct mli_BoundaryLayer layer =
                scenery->materials.boundary_layers.array[idx];
        struct mli_IntersectionLayer ilay = mli_IntersectionLayer_init();

        const struct mli_Surface *inner_surface =
                &scenery->materials.surfaces.array[layer.inner.surface];
        const struct mli_Medium *inner_medium =
                &scenery->materials.media.array[layer.inner.medium];
        const struct mli_Surface *outer_surface =
                &scenery->materials.surfaces.array[layer.outer.surface];
        const struct mli_Medium *outer_medium =
                &scenery->materials.media.array[layer.outer.medium];

        if (isec->from_outside_to_inside) {
                ilay.side_coming_from.surface = outer_surface;
                ilay.side_coming_from.medium = outer_medium;
                ilay.side_coming_from.surface_idx = layer.outer.surface;
                ilay.side_coming_from.medium_idx = layer.outer.medium;

                ilay.side_going_to.surface = inner_surface;
                ilay.side_going_to.medium = inner_medium;
                ilay.side_going_to.surface_idx = layer.inner.surface;
                ilay.side_going_to.medium_idx = layer.inner.medium;
        } else {
                ilay.side_coming_from.surface = inner_surface;
                ilay.side_coming_from.medium = inner_medium;
                ilay.side_coming_from.surface_idx = layer.inner.surface;
                ilay.side_coming_from.medium_idx = layer.inner.medium;

                ilay.side_going_to.surface = outer_surface;
                ilay.side_going_to.medium = outer_medium;
                ilay.side_going_to.surface_idx = layer.outer.surface;
                ilay.side_going_to.medium_idx = layer.outer.medium;
        }

        return ilay;
}

struct mli_BoundaryLayer_Side mli_raytracing_get_side_coming_from(
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec)
{
        const uint64_t idx = mli_Scenery_resolve_boundary_layer_idx(
                scenery, isec->geometry_id);
        struct mli_BoundaryLayer layer =
                scenery->materials.boundary_layers.array[idx];
        if (isec->from_outside_to_inside)
                return layer.outer;
        else
                return layer.inner;
}

struct mli_BoundaryLayer_Side mli_raytracing_get_side_going_to(
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec)
{
        const uint64_t idx = mli_Scenery_resolve_boundary_layer_idx(
                scenery, isec->geometry_id);
        struct mli_BoundaryLayer layer =
                scenery->materials.boundary_layers.array[idx];
        if (isec->from_outside_to_inside)
                return layer.inner;
        else
                return layer.outer;
}

const struct mli_Func *mli_raytracing_get_refractive_index_going_to(
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec)
{
        const struct mli_BoundaryLayer_Side going_to =
                mli_raytracing_get_side_going_to(scenery, isec);
        const uint64_t idx = scenery->materials.media.array[going_to.medium]
                                     .refraction_spectrum;
        return &scenery->materials.spectra.array[idx].spectrum;
}

const struct mli_Func *mli_raytracing_get_refractive_index_coming_from(
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec)
{
        const struct mli_BoundaryLayer_Side coming_from =
                mli_raytracing_get_side_coming_from(scenery, isec);
        const uint64_t idx = scenery->materials.media.array[coming_from.medium]
                                     .refraction_spectrum;
        return &scenery->materials.spectra.array[idx].spectrum;
}

/* intersection_surface_normal */
/* --------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_IntersectionSurfaceNormal mli_IntersectionSurfaceNormal_init(void)
{
        struct mli_IntersectionSurfaceNormal isec;
        isec.geometry_id = mli_GeometryId_init();
        isec.position = mli_Vec_init(0.0, 0.0, 0.0);
        isec.surface_normal = mli_Vec_init(0.0, 0.0, 1.0);
        isec.position_local = mli_Vec_init(0.0, 0.0, 0.0);
        isec.surface_normal_local = mli_Vec_init(0.0, 0.0, 1.0);
        isec.distance_of_ray = DBL_MAX;
        isec.from_outside_to_inside = 1;
        return isec;
}

/* io */
/* -- */

/* Copyright 2018-2023 Sebastian Achim Mueller */

struct mli_IO mli_IO_init(void)
{
        struct mli_IO out;
        out.type = MLI_IO_TYPE_VOID;
        return out;
}

int mli_IO_open_memory(struct mli_IO *self)
{
        mli_IO_close(self);
        self->type = MLI_IO_TYPE_MEMORY;
        self->data.memory = mli_IoMemory_init();
        return mli_IoMemory_open(&self->data.memory);
}

int mli_IO_open_file(
        struct mli_IO *self,
        const struct mli_String *filename,
        const struct mli_String *mode)
{
        mli_IO_close(self);
        self->type = MLI_IO_TYPE_FILE;
        self->data.file = mli_IoFile_init();
        return mli_IoFile_open(&self->data.file, filename, mode);
}

int mli_IO_adopt_file(struct mli_IO *self, FILE *cfile)
{
        mli_IO_close(self);
        self->type = MLI_IO_TYPE_FILE;
        self->data.file = mli_IoFile_init();
        return mli_IoFile_adopt_cfile(&self->data.file, cfile);
}

int mli_IO__open_file_cstr(
        struct mli_IO *self,
        const char *filename,
        const char *mode)
{
        struct mli_String _filename = mli_String_init();
        struct mli_String _mode = mli_String_init();

        chk(mli_String_from_cstr(&_filename, filename));
        chk(mli_String_from_cstr(&_mode, mode));
        chk(mli_IO_open_file(self, &_filename, &_mode));

        mli_String_free(&_filename);
        mli_String_free(&_mode);
        return 1;
chk_error:
        mli_String_free(&_filename);
        mli_String_free(&_mode);
        mli_IO_close(self);
        return 0;
}

int mli_IO_close(struct mli_IO *self)
{
        int rc = 0;
        switch (self->type) {
        case MLI_IO_TYPE_MEMORY:
                rc = mli_IoMemory_close(&self->data.memory);
                break;
        case MLI_IO_TYPE_FILE:
                rc = mli_IoFile_close(&self->data.file);
                break;
        case MLI_IO_TYPE_VOID:
                rc = 1;
                break;
        }
        (*self) = mli_IO_init();
        return rc;
}

size_t mli_IO_write(
        const void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IO *self)
{
        size_t rc = 0;
        switch (self->type) {
        case MLI_IO_TYPE_MEMORY:
                rc = mli_IoMemory_write(ptr, size, count, &self->data.memory);
                break;
        case MLI_IO_TYPE_FILE:
                rc = mli_IoFile_write(ptr, size, count, &self->data.file);
                break;
        default:
                rc = 0;
        }
        return rc;
}

size_t mli_IO_read(
        void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IO *self)
{
        size_t rc;
        switch (self->type) {
        case MLI_IO_TYPE_MEMORY:
                rc = mli_IoMemory_read(ptr, size, count, &self->data.memory);
                break;
        case MLI_IO_TYPE_FILE:
                rc = mli_IoFile_read(ptr, size, count, &self->data.file);
                break;
        default:
                rc = 0;
        }
        return rc;
}

void mli_IO_rewind(struct mli_IO *self)
{
        switch (self->type) {
        case MLI_IO_TYPE_MEMORY:
                mli_IoMemory_rewind(&self->data.memory);
                break;
        case MLI_IO_TYPE_FILE:
                rewind(self->data.file.cfile);
                break;
        }
}

int64_t mli_IO_tell(struct mli_IO *self)
{
        int64_t rc = -1;
        switch (self->type) {
        case MLI_IO_TYPE_MEMORY:
                rc = mli_IoMemory_tell(&self->data.memory);
                break;
        case MLI_IO_TYPE_FILE:
                rc = mli_IoFile_tell(&self->data.file);
                break;
        }
        return rc;
}

int64_t mli_IO_seek(
        struct mli_IO *self,
        const int64_t offset,
        const int64_t origin)
{
        int64_t rc = -1;
        switch (self->type) {
        case MLI_IO_TYPE_MEMORY:
                rc = mli_IoMemory_seek(&self->data.memory, offset, origin);
                break;
        case MLI_IO_TYPE_FILE:
                rc = mli_IoFile_seek(&self->data.file, offset, origin);
                break;
        }
        return rc;
}

int mli_IO_eof(const struct mli_IO *self)
{
        int64_t rc = EOF;
        switch (self->type) {
        case MLI_IO_TYPE_MEMORY:
                rc = mli_IoMemory_eof(&self->data.memory);
                break;
        case MLI_IO_TYPE_FILE:
                rc = mli_IoFile_eof(&self->data.file);
                break;
        }
        return rc;
}

int mli_IO_flush(struct mli_IO *self)
{
        int64_t rc = -1;
        switch (self->type) {
        case MLI_IO_TYPE_MEMORY:
                rc = 1;
                break;
        case MLI_IO_TYPE_FILE:
                rc = mli_IoFile_flush(&self->data.file);
                break;
        }
        return rc;
}

/* io_file */
/* ------- */

/* Copyright 2018-2023 Sebastian Achim Mueller */

struct mli_IoFile mli_IoFile_init(void)
{
        struct mli_IoFile out;
        out.cfile = NULL;
        return out;
}

int mli_IoFile_close(struct mli_IoFile *self)
{
        int rc = 1;
        if (self->cfile != NULL) {
                if (!mli_IoFile__cfile_is_stdin_or_stdout_stderr(self)) {
                        int fclose_rc = fclose(self->cfile);
                        if (fclose_rc == EOF) {
                                rc = EOF;
                        }
                }
        }

        (*self) = mli_IoFile_init();
        return rc;
}

int mli_IoFile_open(
        struct mli_IoFile *self,
        const struct mli_String *filename,
        const struct mli_String *mode)
{
        mli_IoFile_close(self);
        self->cfile = fopen(filename->array, mode->array);
        chk_msg(self->cfile != NULL, "Failed to open file.");
        return 1;
chk_error:
        mli_IoFile_close(self);
        return 0;
}

int mli_IoFile_adopt_cfile(struct mli_IoFile *self, FILE *cfile)
{
        mli_IoFile_close(self);
        self->cfile = cfile;
        return 1;
}

size_t mli_IoFile_write(
        const void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IoFile *self)
{
        const size_t actual_count = fwrite(ptr, size, count, self->cfile);
        chk_msg(actual_count == count, "Can not write to file.");

chk_error:
        return actual_count;
}

size_t mli_IoFile_read(
        void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IoFile *self)
{
        const size_t actual_count = fread(ptr, size, count, self->cfile);
        chk_msg(actual_count == count, "Can not read from file.");

chk_error:
        return actual_count;
}

void mli_IoFile_rewind(struct mli_IoFile *self) { rewind(self->cfile); }

int64_t mli_IoFile_tell(struct mli_IoFile *self) { return ftell(self->cfile); }

int64_t mli_IoFile_seek(
        struct mli_IoFile *self,
        const int64_t offset,
        const int64_t origin)
{
        return fseek(self->cfile, offset, origin);
}

int mli_IoFile_eof(const struct mli_IoFile *self) { return feof(self->cfile); }

int mli_IoFile_flush(struct mli_IoFile *self) { return fflush(self->cfile); }

int mli_IoFile__cfile_is_stdin_or_stdout_stderr(const struct mli_IoFile *self)
{
        if (self->cfile == stdin) {
                return 1;
        }
        if (self->cfile == stdout) {
                return 1;
        }
        if (self->cfile == stderr) {
                return 1;
        }
        return 0;
}

/* io_memory */
/* --------- */

/* Copyright 2018-2023 Sebastian Achim Mueller */

struct mli_IoMemory mli_IoMemory_init(void)
{
        struct mli_IoMemory byt;
        byt.cstr = NULL;
        byt.capacity = 0u;
        byt.size = 0u;
        byt.pos = 0u;
        return byt;
}

int mli_IoMemory_close(struct mli_IoMemory *self)
{
        const int RC_IS_ALWAYS_1 = 1;
        free(self->cstr);
        (*self) = mli_IoMemory_init();
        return RC_IS_ALWAYS_1;
}

int mli_IoMemory__malloc_capacity(
        struct mli_IoMemory *self,
        const uint64_t capacity)
{
        mli_IoMemory_close(self);
        self->capacity = MLI_MATH_MAX2(2u, capacity);
        self->size = 0u;
        chk_malloc(self->cstr, unsigned char, self->capacity + 1u);
        memset(self->cstr, '\0', self->capacity + 1u);
        return 1;
chk_error:
        return 0;
}

int mli_IoMemory__malloc(struct mli_IoMemory *self)
{
        chk(mli_IoMemory__malloc_capacity(self, 0u));
        return 1;
chk_error:
        return 0;
}

int mli_IoMemory__realloc_capacity(
        struct mli_IoMemory *self,
        const uint64_t new_capacity)
{
        uint64_t numcpy;
        struct mli_IoMemory tmp = mli_IoMemory_init();

        chk_msg(mli_IoMemory__malloc_capacity(&tmp, self->capacity),
                "Failed to allocate temporary swap.");
        memcpy(tmp.cstr, self->cstr, self->capacity);
        tmp.pos = self->pos;
        tmp.size = self->size;

        chk_msg(mli_IoMemory__malloc_capacity(self, new_capacity),
                "Faild to allocate new capacity.");

        if (new_capacity >= tmp.capacity) {
                /* growing or same */
                numcpy = tmp.capacity;
        } else {
                /* shrinking */
                numcpy = new_capacity;
                chk_msg(tmp.pos <= new_capacity,
                        "Expected cursor 'pos' to be within new capacity.");
                chk_msg(tmp.size <= new_capacity,
                        "Expected 'size' to be within new capacity.");
        }

        memcpy(self->cstr, tmp.cstr, numcpy);
        self->pos = tmp.pos;
        self->size = tmp.size;

        mli_IoMemory_close(&tmp);
        return 1;
chk_error:
        return 0;
}

int mli_IoMemory_open(struct mli_IoMemory *self)
{
        if (self->cstr == NULL) {
                chk(mli_IoMemory__malloc(self));
        } else {
                chk_msg(self->capacity >= 2u,
                        "Expected minimal capacity of 2.");
                self->pos = 0u;
                self->size = 0u;
                memset(self->cstr, '\0', self->capacity + 1u);
        }
        return 1;
chk_error:
        return 0;
}

int mli_IoMemory__shrink_to_fit(struct mli_IoMemory *self)
{
        chk_msg(mli_IoMemory__realloc_capacity(self, self->size),
                "Failed to reallocate to size.");
        return 1;
chk_error:
        return 0;
}

int mli_IoMemory__write_unsigned_char(
        struct mli_IoMemory *self,
        const unsigned char *c)
{
        /* 'puts' a single selfe (char) to the selfesIo-buffer */
        uint64_t new_size;

        if (self->pos == self->size) {
                new_size = self->size + 1u;
        } else {
                new_size = self->size;
        }

        if (self->cstr == NULL) {
                chk(mli_IoMemory__malloc(self));
        }

        if (new_size >= self->capacity) {
                const uint64_t min_new_capacity =
                        MLI_MATH_MAX2(new_size, 2 * self->capacity);
                chk_msg(mli_IoMemory__realloc_capacity(self, min_new_capacity),
                        "Failed to reallocate.");
        }
        self->cstr[self->pos] = (*c);
        self->size = new_size;
        self->pos += 1;

        return 1;
chk_error:
        return 0;
}

int mli_IoMemory__read_unsigned_char(
        struct mli_IoMemory *self,
        unsigned char *c)
{
        if (self->pos >= self->size) {
                return EOF;
        }
        (*c) = self->cstr[self->pos];
        self->pos += 1;
        return 1;
}

size_t mli_IoMemory_write(
        const void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IoMemory *self)
{
        const size_t block_size = size * count;
        unsigned char *block = (unsigned char *)ptr;

        size_t i;
        for (i = 0u; i < block_size; i++) {
                chk_msg(mli_IoMemory__write_unsigned_char(self, &block[i]), "");
        }
chk_error:
        return (i + 1u) / size;
}

size_t mli_IoMemory_read(
        void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IoMemory *self)
{
        const size_t block_size = size * count;
        unsigned char *block = (unsigned char *)ptr;

        size_t i;
        for (i = 0u; i < block_size; i++) {
                unsigned char c;
                const int rc = mli_IoMemory__read_unsigned_char(self, &c);
                if (rc == EOF) {
                        return i / size;
                        ;
                } else {
                        block[i] = c;
                }
        }
        return (i + 1u) / size;
}

void mli_IoMemory_rewind(struct mli_IoMemory *self) { self->pos = 0u; }

int64_t mli_IoMemory_tell(struct mli_IoMemory *self)
{
        return (int64_t)self->pos;
}

int64_t mli_IoMemory_seek(
        struct mli_IoMemory *self,
        const int64_t offset,
        const int64_t origin)
{
        const int64_t SUCCESS = 0;
        const int64_t FAIL = -1;
        /* If successful, the function returns zero.
         * Otherwise, it returns non-zero value.
         */
        if (origin == SEEK_CUR) {
                const int64_t new_pos = self->pos + offset;
                return mli_IoMemory_seek(self, new_pos, SEEK_SET);
        } else if (origin == SEEK_SET) {
                if (offset < 0 || offset >= (int64_t)self->size) {
                        return FAIL;
                } else {
                        self->pos = offset;
                        return SUCCESS;
                }
        } else {
                return FAIL;
        }
}

int mli_IoMemory_eof(const struct mli_IoMemory *self)
{
        if (self->pos < self->size) {
                return 0;
        } else {
                return EOF;
        }
}

int mli_IoMemory__write_cstr(struct mli_IoMemory *self, const char *cstr)
{
        /* For testing purposes */
        size_t i = 0;
        while (cstr[i] != '\0') {
                chk(mli_IoMemory_write(
                        (unsigned char *)(&cstr[i]), sizeof(char), 1, self));
                i += 1;
        }

        return 1;
chk_error:
        return 0;
}

/* io_text */
/* ------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

int mli_IO_text_getc(struct mli_IO *self)
{
        char c = EOF;
        size_t rc = mli_IO_read((void *)(&c), sizeof(char), 1, self);
        if (rc == 0) {
                return EOF;
        } else {
                return (int)c;
        }
}

int mli_IO_text_putc(struct mli_IO *self, const char c)
{
        chk_msg(mli_IO_write((void *)(&c), sizeof(char), 1, self),
                "Can not write char to IO.");
        return 1;
chk_error:
        return 0;
}

int mli_IO_text_write_cstr(struct mli_IO *self, const char *cstr)
{
        struct mli_String tmp = mli_String_init();
        chk_msg(mli_String_from_cstr(&tmp, cstr),
                "Can't malloc mli_String from cstr");
        chk(mli_IO_text_write_String(self, &tmp));
        mli_String_free(&tmp);
        return 1;
chk_error:
        mli_String_free(&tmp);
        return 0;
}

int mli_IO_text_write_cstr_format(struct mli_IO *self, const char *format, ...)
{
        struct mli_String tmp = mli_String_init();
        va_list args;
        va_start(args, format);
        chk(mli_String_from_vargs(&tmp, format, args));
        chk(mli_IO_text_write_String(self, &tmp));
        va_end(args);
        mli_String_free(&tmp);
        return 1;
chk_error:
        va_end(args);
        mli_String_free(&tmp);
        return 0;
}

int mli_IO_text_read_string(struct mli_IO *self, struct mli_String *str)
{
        chk(mli_IO_text_read_line(self, str, '\0'));
        return 1;
chk_error:
        return 0;
}

int mli_IO_text_read_line(
        struct mli_IO *self,
        struct mli_String *line,
        const char delimiter)
{
        chk(mli_String_malloc(line, 0));

        if (mli_IO_eof(self)) {
                return EOF;
        }

        while (!mli_IO_eof(self)) {
                const int c = mli_IO_text_getc(self);
                if (c == delimiter) {
                        break;
                } else if (c == '\0') {
                        break;
                } else if (c == EOF) {
                        break;
                } else {
                        chk(mli_String_push_back(line, c));
                }
        }

        return 1;
chk_error:
        return 0;
}

int mli_IO_text_write_String(struct mli_IO *self, const struct mli_String *str)
{
        uint64_t i;
        for (i = 0; i < str->size; i++) {
                char c = str->array[i];
                if (c == '\0') {
                        break;
                }
                chk(mli_IO_text_putc(self, c));
        }
        return 1;
chk_error:
        return 0;
}

int mli_IO_text_write_multi_line_debug_view_line_match(
        struct mli_IO *self,
        const int64_t line_number,
        const int64_t line_number_of_interest)
{
        chk(mli_IO_text_write_cstr_format(self, "% 6d", (int32_t)line_number));
        if (line_number == line_number_of_interest) {
                chk(mli_IO_text_write_cstr_format(self, "->|  "));
        } else {
                chk(mli_IO_text_write_cstr_format(self, "  |  "));
        }
        return 1;
chk_error:
        return 0;
}

int mli_IO_text_write_multi_line_debug_view(
        struct mli_IO *self,
        const struct mli_String *text,
        const uint64_t line_number,
        const uint64_t line_radius)
{
        int64_t _line_number = (int64_t)line_number;
        int64_t _line_radius = (int64_t)line_radius;
        int64_t line_start = MLI_MATH_MAX2(_line_number - _line_radius, 1);
        int64_t line_stop = line_number + line_radius;
        int64_t line = 1;
        uint64_t i = 0;

        chk_msg(line_radius > 1, "Expected line_radius > 1.");

        chk(mli_IO_text_write_cstr_format(self, "  line     text\n"));
        chk(mli_IO_text_write_cstr_format(self, "        |\n"));

        while (i < text->size && text->array[i]) {
                int prefix = (line + 1 >= line_start) && (line < line_stop);
                int valid = (line >= line_start) && (line <= line_stop);
                if (text->array[i] == '\n') {
                        line++;
                }
                if (prefix && i == 0) {
                        chk(mli_IO_text_write_multi_line_debug_view_line_match(
                                self, line, _line_number));
                }
                if (valid) {
                        chk(mli_IO_text_putc(self, text->array[i]));
                }
                if (prefix && text->array[i] == '\n') {
                        chk(mli_IO_text_write_multi_line_debug_view_line_match(
                                self, line, _line_number));
                }
                i++;
        }
        chk(mli_IO_text_putc(self, '\n'));

        return 1;
chk_error:
        return 0;
}

/* json */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Json mli_Json_init(void)
{
        struct mli_Json j;
        j.raw = mli_String_init();
        j.num_tokens = 0u;
        j.tokens = NULL;
        return j;
}

void mli_Json_free(struct mli_Json *json)
{
        mli_String_free(&json->raw);
        free(json->tokens);
        (*json) = mli_Json_init();
}

int mli_Json__malloc_tokens(struct mli_Json *json)
{
        struct jsmntok_t default_token = {JSMN_UNDEFINED, 0, 0, 0};
        chk_msg(&json->raw.array != NULL, "Expected raw cstr to be malloced.");
        json->num_tokens = json->raw.size / 2;
        chk_malloc(json->tokens, struct jsmntok_t, json->num_tokens);
        MLI_MATH_ARRAY_SET(json->tokens, default_token, json->num_tokens);
        return 1;
chk_error:
        return 0;
}

int mli_Json__parse_tokens(struct mli_Json *json)
{
        int64_t num_tokens_parsed;
        struct jsmn_parser parser;

        chk_msg(&json->tokens != NULL, "Expected tokens to be malloced.");
        jsmn_init(&parser);
        num_tokens_parsed = jsmn_parse(
                &parser,
                json->raw.array,
                json->raw.size,
                json->tokens,
                json->num_tokens);
        chk_msgf(
                num_tokens_parsed != JSMN_ERROR_NOMEM,
                ("Not enough tokens. Only got " PRIu64, json->num_tokens));
        chk_msg(num_tokens_parsed != JSMN_ERROR_INVAL,
                "Invalid character inside JSON string.");
        chk_msg(num_tokens_parsed != JSMN_ERROR_PART,
                "The string is not a full JSON packet, more "
                "bytes expected.");
        chk_msg(num_tokens_parsed >= 0, "Can't parse Json-string");
        json->num_tokens = num_tokens_parsed;
        return 1;
chk_error:
        return 0;
}

int mli_Json_from_string(struct mli_Json *self, const struct mli_String *str)
{
        mli_Json_free(self);
        chk_msg(mli_String_copy(&self->raw, str),
                "Failed to copy string into json->raw.");
        chk_msg(mli_Json__malloc_tokens(self), "Can't malloc Json's tokens.");
        chk_msg(mli_Json__parse_tokens(self), "Can't parse Json into tokens.");
        return 1;
chk_error:
        mli_Json_free(self);
        return 0;
}

int mli_Json_from_io(struct mli_Json *self, struct mli_IO *io)
{
        mli_Json_free(self);
        chk_msg(mli_IO_text_read_string(io, &self->raw),
                "Failed to read file into String.");
        chk_msg(mli_Json__malloc_tokens(self), "Can't malloc Json's tokens.");
        chk_msg(mli_Json__parse_tokens(self), "Can't parse Json into tokens.");
        return 1;
chk_error:
        mli_Json_free(self);
        return 0;
}

int mli_Json_cstr_by_token(
        const struct mli_Json *json,
        const uint64_t token,
        char *return_string,
        const uint64_t return_string_size)
{
        const struct jsmntok_t t = json->tokens[token];
        const uint64_t actual_length = t.end - t.start;
        chk_msg(actual_length < return_string_size,
                "Expected return_string_size to be sufficiently large for "
                "json-string, but it is not.");
        memcpy(return_string, json->raw.array + t.start, actual_length);
        return_string[actual_length] = '\0';
        return 1;
chk_error:
        return 0;
}

int mli_Json_string_by_token(
        const struct mli_Json *json,
        const uint64_t token,
        struct mli_String *return_string)
{
        const struct jsmntok_t t = json->tokens[token];
        const uint64_t size = t.end - t.start;
        chk(mli_String_malloc(return_string, size));
        memcpy(return_string->array, json->raw.array + t.start, size);
        return_string->size = size;
        return 1;
chk_error:
        return 0;
}

int mli_Json_int64_by_token(
        const struct mli_Json *json,
        const uint64_t token,
        int64_t *return_int64)
{
        const struct jsmntok_t t = json->tokens[token];
        const uint64_t token_length = t.end - t.start;
        chk_msg(t.type == JSMN_PRIMITIVE,
                "Json int64 expected json-token-to be JSMN_PRIMITIVE.");
        chk_msg(mli_cstr_nto_int64(
                        return_int64,
                        (char *)&json->raw.array[t.start],
                        10,
                        token_length),
                "Can't parse int64.");
        return 1;
chk_error:
        return 0;
}

int mli_Json_uint64_by_token(
        const struct mli_Json *json,
        const uint64_t token,
        uint64_t *val)
{
        int64_t tmp;
        chk(mli_Json_int64_by_token(json, token, &tmp));
        chk_msg(tmp >= 0, "Expected value to be unsigned.");
        (*val) = (uint64_t)tmp;
        return 1;
chk_error:
        return 0;
}

int mli_Json_int64_by_key(
        const struct mli_Json *json,
        const uint64_t token,
        int64_t *val,
        const char *key)
{
        uint64_t token_n;
        chk(mli_Json_token_by_key_eprint(json, token, key, &token_n));
        chk_msgf(
                mli_Json_int64_by_token(json, token_n + 1, val),
                ("Can't parse value of '%s' into int64.", key));
        return 1;
chk_error:
        return 0;
}

int mli_Json_uint64_by_key(
        const struct mli_Json *json,
        const uint64_t token,
        uint64_t *val,
        const char *key)
{
        int64_t tmp;
        chk(mli_Json_int64_by_key(json, token, &tmp, key));
        chk_msg(tmp >= 0, "Expected value to be unsigned.");
        (*val) = (uint64_t)tmp;
        return 1;
chk_error:
        return 0;
}

int mli_Json_double_by_token(
        const struct mli_Json *json,
        const uint64_t token,
        double *val)
{
        const struct jsmntok_t t = json->tokens[token];
        const uint64_t token_length = t.end - t.start;
        chk_msg(t.type == JSMN_PRIMITIVE,
                "Json float64 expected json-token-to be JSMN_PRIMITIVE.");
        chk_msg(mli_cstr_nto_double(
                        val, (char *)&json->raw.array[t.start], token_length),
                "Can't parse double.");
        return 1;
chk_error:
        return 0;
}

int mli_Json_double_by_key(
        const struct mli_Json *json,
        const uint64_t token,
        double *val,
        const char *key)
{
        uint64_t token_n;
        chk(mli_Json_token_by_key_eprint(json, token, key, &token_n));
        chk_msgf(
                mli_Json_double_by_token(json, token_n + 1, val),
                ("Can't parse value of '%s' into double.", key));

        return 1;
chk_error:
        return 0;
}

int mli_Json_cstrcmp(
        const struct mli_Json *json,
        const uint64_t token,
        const char *str)
{
        uint64_t i;
        const struct jsmntok_t t = json->tokens[token];
        const uint64_t token_length = t.end - t.start;
        const uint64_t str_length = strlen(str);

        if (token_length != str_length) {
                return 0;
        }
        for (i = 0; i < token_length; i++) {
                const char token_char = (char)json->raw.array[t.start + i];
                const char str_char = str[i];
                if (token_char != str_char) {
                        return 0;
                }
        }
        return 1;
}

int mli_Json_token_by_key(
        const struct mli_Json *json,
        const uint64_t token,
        const char *key,
        uint64_t *key_token)
{
        int64_t found = 0;
        int64_t child = 0;
        int64_t subchild_balance = 0;
        int64_t idx = token + 1;

        while (child < json->tokens[token].size) {
                if (mli_Json_cstrcmp(json, idx, key)) {
                        (*key_token) = idx;
                        found += 1;
                }
                subchild_balance += json->tokens[idx].size;
                while (subchild_balance > 0) {
                        idx += 1;
                        subchild_balance += json->tokens[idx].size;
                        subchild_balance -= 1;
                }
                idx += 1;
                child += 1;
        }
        return found;
}

int mli_Json_token_by_key_eprint(
        const struct mli_Json *json,
        const uint64_t token,
        const char *key,
        uint64_t *key_token)
{
        chk_msgf(
                mli_Json_token_by_key(json, token, key, key_token),
                ("Expected key '%s' in json.", key));
        return 1;
chk_error:
        return 0;
}

uint64_t mli_Json__token_by_index_unsafe(
        const struct mli_Json *json,
        const uint64_t token,
        const uint64_t index)
{
        uint64_t child = 0;
        int64_t subchild_balance = 0;
        uint64_t idx = token + 1;

        while (child < index) {
                subchild_balance += json->tokens[idx].size;
                while (subchild_balance > 0) {
                        idx += 1;
                        subchild_balance += json->tokens[idx].size;
                        subchild_balance -= 1;
                }
                idx += 1;
                child += 1;
        }
        return idx;
}

int mli_Json_token_by_idx(
        const struct mli_Json *json,
        const uint64_t token,
        const uint64_t idx,
        uint64_t *idx_token)
{
        chk_msgf(
                (int64_t)idx < json->tokens[token].size,
                ("Array idx '%lu' is out of range for size '%lu'.",
                 token,
                 json->tokens[token].size));
        (*idx_token) = mli_Json__token_by_index_unsafe(json, token, idx);
        return 1;
chk_error:
        return 0;
}

int mli_Json_debug_token_fprint(
        FILE *f,
        const struct mli_Json *self,
        const uint64_t token)
{
        struct mli_IO io = mli_IO_init();
        mli_IO_adopt_file(&io, f);
        chk(mli_Json_debug_token_to_io(self, token, &io));
        return 1;
chk_error:
        return 0;
}

int mli_Json_debug_token_to_io(
        const struct mli_Json *self,
        const uint64_t token,
        struct mli_IO *io)
{
        uint64_t i = 0u;
        struct jsmntok_t t = self->tokens[token];
        uint32_t token_size = t.end - t.start;
        uint64_t line_number =
                1u + mli_String_countn(&self->raw, '\n', t.start);

        chk(mli_IO_text_write_cstr_format(
                io, "line: %u, ", (uint32_t)line_number));
        chk(mli_IO_text_write_cstr_format(io, "token: %u, ", (uint32_t)token));
        chk(mli_IO_text_write_cstr_format(io, "type: %d, ", t.type));
        chk(mli_IO_text_write_cstr_format(io, "children: %d, ", t.size));
        chk(mli_IO_text_write_cstr_format(
                io, "chars: (%d -> %d, %d)\n", t.start, t.end, token_size));

        for (i = 0; i < token_size; i++) {
                chk(mli_IO_text_putc(io, self->raw.array[t.start + i]));
        }
        chk(mli_IO_text_putc(io, '\n'));
        return 1;
chk_error:
        return 0;
}

int mli_Json_debug_to_io(const struct mli_Json *self, struct mli_IO *io)
{
        uint64_t i;
        for (i = 0; i < self->num_tokens; i++) {
                chk_msg(mli_Json_debug_token_to_io(self, i, io),
                        "Failed to write json-token debug-info to io.");
        }
        return 1;
chk_error:
        return 0;
}

/* json_jsmn */
/* --------- */

/*
 * MIT License
 *
 * Copyright (c) 2010 Serge Zaitsev
 *               2018-2020 Sebastian Achim Mueller
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 *
 * - removed JSMN_PARENT_LINKS
 * - always JSMN_STRICT
 */


void jsmn_init(struct jsmn_parser *parser);

/**
 * Run JSON parser. It parses a JSON data string into and array of tokens, each
 * describing
 * a single JSON object.
 */
int jsmn_parse(
        struct jsmn_parser *parser,
        const char *js,
        const size_t len,
        struct jsmntok_t *tokens,
        const unsigned int num_tokens);

/**
 * Allocates a fresh unused token from the token pool.
 */
static struct jsmntok_t *jsmn_alloc_token(
        struct jsmn_parser *parser,
        struct jsmntok_t *tokens,
        const size_t num_tokens)
{
        struct jsmntok_t *tok;
        if (parser->toknext >= num_tokens) {
                return NULL;
        }
        tok = &tokens[parser->toknext++];
        tok->start = tok->end = -1;
        tok->size = 0;
        return tok;
}

/**
 * Fills token type and boundaries.
 */
static void jsmn_fill_token(
        struct jsmntok_t *token,
        const enum jsmntype_t type,
        const int start,
        const int end)
{
        token->type = type;
        token->start = start;
        token->end = end;
        token->size = 0;
}

/**
 * Fills next available token with JSON primitive.
 */
static int jsmn_parse_primitive(
        struct jsmn_parser *parser,
        const char *js,
        const size_t len,
        struct jsmntok_t *tokens,
        const size_t num_tokens)
{
        struct jsmntok_t *token;
        int start;

        start = parser->pos;

        for (; parser->pos < len && js[parser->pos] != '\0'; parser->pos++) {
                switch (js[parser->pos]) {
                case '\t':
                case '\r':
                case '\n':
                case ' ':
                case ',':
                case ']':
                case '}':
                        goto found;
                }
                if (js[parser->pos] < 32 || js[parser->pos] >= 127) {
                        parser->pos = start;
                        return JSMN_ERROR_INVAL;
                }
        }

        /* In strict mode primitive must be followed by a comma/object/array */
        parser->pos = start;
        return JSMN_ERROR_PART;

found:
        if (tokens == NULL) {
                parser->pos--;
                return 0;
        }
        token = jsmn_alloc_token(parser, tokens, num_tokens);
        if (token == NULL) {
                parser->pos = start;
                return JSMN_ERROR_NOMEM;
        }
        jsmn_fill_token(token, JSMN_PRIMITIVE, start, parser->pos);
        parser->pos--;
        return 0;
}

/**
 * Fills next token with JSON string.
 */
static int jsmn_parse_string(
        struct jsmn_parser *parser,
        const char *js,
        const size_t len,
        struct jsmntok_t *tokens,
        const size_t num_tokens)
{
        struct jsmntok_t *token;

        int start = parser->pos;

        parser->pos++;

        /* Skip starting quote */
        for (; parser->pos < len && js[parser->pos] != '\0'; parser->pos++) {
                char c = js[parser->pos];

                /* Quote: end of string */
                if (c == '\"') {
                        if (tokens == NULL) {
                                return 0;
                        }
                        token = jsmn_alloc_token(parser, tokens, num_tokens);
                        if (token == NULL) {
                                parser->pos = start;
                                return JSMN_ERROR_NOMEM;
                        }
                        jsmn_fill_token(
                                token, JSMN_STRING, start + 1, parser->pos);
                        return 0;
                }

                /* Backslash: Quoted symbol expected */
                if (c == '\\' && parser->pos + 1 < len) {
                        int i;
                        parser->pos++;
                        switch (js[parser->pos]) {
                        /* Allowed escaped symbols */
                        case '\"':
                        case '/':
                        case '\\':
                        case 'b':
                        case 'f':
                        case 'r':
                        case 'n':
                        case 't':
                                break;
                        /* Allows escaped symbol \uXXXX */
                        case 'u':
                                parser->pos++;
                                for (i = 0; i < 4 && parser->pos < len &&
                                            js[parser->pos] != '\0';
                                     i++) {
                                        /* If it isn't a hex character we have
                                         * an error */
                                        if (!((js[parser->pos] >= 48 &&
                                               js[parser->pos] <=
                                                       57) || /* 0-9 */
                                              (js[parser->pos] >= 65 &&
                                               js[parser->pos] <=
                                                       70) || /* A-F */
                                              (js[parser->pos] >= 97 &&
                                               js[parser->pos] <=
                                                       102))) { /* a-f */
                                                parser->pos = start;
                                                return JSMN_ERROR_INVAL;
                                        }
                                        parser->pos++;
                                }
                                parser->pos--;
                                break;
                        /* Unexpected symbol */
                        default:
                                parser->pos = start;
                                return JSMN_ERROR_INVAL;
                        }
                }
        }
        parser->pos = start;
        return JSMN_ERROR_PART;
}

/**
 * Parse JSON string and fill tokens.
 */
int jsmn_parse(
        struct jsmn_parser *parser,
        const char *js,
        const size_t len,
        struct jsmntok_t *tokens,
        const unsigned int num_tokens)
{
        int r;
        int i;
        struct jsmntok_t *token;
        int count = parser->toknext;

        for (; parser->pos < len && js[parser->pos] != '\0'; parser->pos++) {
                char c;
                enum jsmntype_t type;

                c = js[parser->pos];
                switch (c) {
                case '{':
                case '[':
                        count++;
                        if (tokens == NULL) {
                                break;
                        }
                        token = jsmn_alloc_token(parser, tokens, num_tokens);
                        if (token == NULL) {
                                return JSMN_ERROR_NOMEM;
                        }
                        if (parser->toksuper != -1) {
                                struct jsmntok_t *t = &tokens[parser->toksuper];

                                /* In strict mode an object or array can't
                                 * become a key */
                                if (t->type == JSMN_OBJECT) {
                                        return JSMN_ERROR_INVAL;
                                }

                                t->size++;
                        }
                        token->type = (c == '{' ? JSMN_OBJECT : JSMN_ARRAY);
                        token->start = parser->pos;
                        parser->toksuper = parser->toknext - 1;
                        break;
                case '}':
                case ']':
                        if (tokens == NULL) {
                                break;
                        }
                        type = (c == '}' ? JSMN_OBJECT : JSMN_ARRAY);
                        for (i = parser->toknext - 1; i >= 0; i--) {
                                token = &tokens[i];
                                if (token->start != -1 && token->end == -1) {
                                        if (token->type != type) {
                                                return JSMN_ERROR_INVAL;
                                        }
                                        parser->toksuper = -1;
                                        token->end = parser->pos + 1;
                                        break;
                                }
                        }
                        /* Error if unmatched closing bracket */
                        if (i == -1) {
                                return JSMN_ERROR_INVAL;
                        }
                        for (; i >= 0; i--) {
                                token = &tokens[i];
                                if (token->start != -1 && token->end == -1) {
                                        parser->toksuper = i;
                                        break;
                                }
                        }
                        break;
                case '\"':
                        r = jsmn_parse_string(
                                parser, js, len, tokens, num_tokens);
                        if (r < 0) {
                                return r;
                        }
                        count++;
                        if (parser->toksuper != -1 && tokens != NULL) {
                                tokens[parser->toksuper].size++;
                        }
                        break;
                case '\t':
                case '\r':
                case '\n':
                case ' ':
                        break;
                case ':':
                        parser->toksuper = parser->toknext - 1;
                        break;
                case ',':
                        if (tokens != NULL && parser->toksuper != -1 &&
                            tokens[parser->toksuper].type != JSMN_ARRAY &&
                            tokens[parser->toksuper].type != JSMN_OBJECT) {
                                for (i = parser->toknext - 1; i >= 0; i--) {
                                        if (tokens[i].type == JSMN_ARRAY ||
                                            tokens[i].type == JSMN_OBJECT) {
                                                if (tokens[i].start != -1 &&
                                                    tokens[i].end == -1) {
                                                        parser->toksuper = i;
                                                        break;
                                                }
                                        }
                                }
                        }
                        break;

                /* In strict mode primitives are: numbers and booleans */
                case '-':
                case '0':
                case '1':
                case '2':
                case '3':
                case '4':
                case '5':
                case '6':
                case '7':
                case '8':
                case '9':
                case 't':
                case 'f':
                case 'n':
                        /* And they must not be keys of the object */
                        if (tokens != NULL && parser->toksuper != -1) {
                                const struct jsmntok_t *t =
                                        &tokens[parser->toksuper];
                                if (t->type == JSMN_OBJECT ||
                                    (t->type == JSMN_STRING && t->size != 0)) {
                                        return JSMN_ERROR_INVAL;
                                }
                        }

                        r = jsmn_parse_primitive(
                                parser, js, len, tokens, num_tokens);
                        if (r < 0) {
                                return r;
                        }
                        count++;
                        if (parser->toksuper != -1 && tokens != NULL) {
                                tokens[parser->toksuper].size++;
                        }
                        break;

                /* Unexpected char in strict mode */
                default:
                        return JSMN_ERROR_INVAL;
                }
        }

        if (tokens != NULL) {
                for (i = parser->toknext - 1; i >= 0; i--) {
                        /* Unmatched opened object or array */
                        if (tokens[i].start != -1 && tokens[i].end == -1) {
                                return JSMN_ERROR_PART;
                        }
                }
        }

        return count;
}

/**
 * Creates a new parser based over a given buffer with an array of tokens
 * available.
 */
void jsmn_init(struct jsmn_parser *parser)
{
        parser->pos = 0;
        parser->toknext = 0;
        parser->toksuper = -1;
}

/* json_walk */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_JsonWalk mli_JsonWalk_init(void)
{
        struct mli_JsonWalk out;
        out.json = NULL;
        out.token = 0;
        return out;
}

struct mli_JsonWalk mli_JsonWalk_set(const struct mli_Json *json)
{
        struct mli_JsonWalk out;
        out.json = json;
        out.token = 0;
        return out;
}

struct mli_JsonWalk mli_JsonWalk_copy(const struct mli_JsonWalk *self)
{
        struct mli_JsonWalk out;
        out.json = self->json;
        out.token = self->token;
        return out;
}

int mli_JsonWalk__type(const struct mli_JsonWalk *self)
{
        return self->json->tokens[self->token].type;
}

int mli_JsonWalk_to_key(struct mli_JsonWalk *self, const char *key)
{
        uint64_t key_token;
        chk_msg(mli_JsonWalk__type(self) == JSMN_OBJECT,
                "Can only go to key-value in json-object.");
        chk_msgf(
                mli_Json_token_by_key(self->json, self->token, key, &key_token),
                ("No key '%s'", key));
        self->token = key_token + 1;
        return 1;
chk_error:
        return 0;
}

int mli_JsonWalk_to_idx(struct mli_JsonWalk *self, const uint64_t idx)
{
        uint64_t idx_token;
        chk_msg(mli_JsonWalk__type(self) == JSMN_ARRAY,
                "Can only go to idx-value in json-array.");
        chk_msgf(
                mli_Json_token_by_idx(self->json, self->token, idx, &idx_token),
                ("No index '%lu'", idx));
        self->token = idx_token;
        return 1;
chk_error:
        return 0;
}

void mli_JsonWalk_to_root(struct mli_JsonWalk *self) { self->token = 0; }

int mli_JsonWalk_get_array_size(const struct mli_JsonWalk *self, uint64_t *size)
{
        chk_msg(mli_JsonWalk__type(self) == JSMN_ARRAY,
                "Can only get size of json-array.");
        (*size) = self->json->tokens[self->token].size;
        return 1;
chk_error:
        return 0;
}

int mli_JsonWalk_get_string(
        const struct mli_JsonWalk *self,
        struct mli_String *val)
{
        return mli_Json_string_by_token(self->json, self->token, val);
}

int mli_JsonWalk_get_int64(const struct mli_JsonWalk *self, int64_t *val)
{
        return mli_Json_int64_by_token(self->json, self->token, val);
}

int mli_JsonWalk_get_double(const struct mli_JsonWalk *self, double *val)
{
        return mli_Json_double_by_token(self->json, self->token, val);
}

/* lambertian */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Vec mli_lambertian_cosine_law_draw_direction_wrt_z(
        struct mli_Prng *prng)
{
        double azimuth;
        double sin_theta, cos_theta;
        azimuth = MLI_MATH_2PI * mli_Prng_uniform(prng);
        sin_theta = mli_Prng_uniform(prng);
        cos_theta = sqrt(1.0 - sin_theta * sin_theta);
        return mli_Vec_init(
                sin_theta * cos(azimuth), sin_theta * sin(azimuth), cos_theta);
}

struct mli_Vec mli_lambertian_cosine_law_draw_direction_wrt_surface_normal(
        struct mli_Prng *prng,
        const struct mli_Vec surface_normal)
{
        const struct mli_Vec z = mli_Vec_init(0, 0, 1);
        const struct mli_Vec lambertian_wrt_z =
                mli_lambertian_cosine_law_draw_direction_wrt_z(prng);
        const double rho = mli_Vec_angle_between(z, surface_normal);
        if (rho > 0.0) {
                const struct mli_Mat rot = mli_Mat_init_axis_angle(
                        mli_Vec_cross(z, surface_normal), -1.0 * rho);
                return mli_transform_orientation(&rot, lambertian_wrt_z);
        } else {
                return lambertian_wrt_z;
        }
}

/* magicid */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_MagicId mli_MagicId_init(void)
{
        struct mli_MagicId magic;
        memset(magic.word, '\0', sizeof(magic.word));
        magic.mayor = MLI_VERSION_MAYOR;
        magic.minor = MLI_VERSION_MINOR;
        magic.patch = MLI_VERSION_PATCH;
        return magic;
}

int mli_MagicId_set(struct mli_MagicId *magic, const char *word)
{
        uint64_t i, len;
        (*magic) = mli_MagicId_init();
        chk_msg(strlen(word) < sizeof(magic->word),
                "Expected magic word to be shorter.");

        len = MLI_MATH_MIN2(strlen(word), sizeof(magic->word));

        for (i = 0; i < len; i++) {
                magic->word[i] = word[i];
        }
        while (i < sizeof(magic->word)) {
                magic->word[i] = '\0';
                i += 1;
        }
        return 1;
chk_error:
        return 0;
}

int mli_MagicId_has_word(const struct mli_MagicId *magic, const char *word)
{
        uint64_t i, len;
        chk_msg(strlen(word) < sizeof(magic->word),
                "Expected magic word to be shorter.");

        len = MLI_MATH_MIN2(strlen(word), sizeof(magic->word));

        for (i = 0; i < len; i++) {
                if (magic->word[i] != word[i]) {
                        return 0;
                }
        }
        while (i < sizeof(magic->word)) {
                if (magic->word[i] != '\0') {
                        return 0;
                }
                i += 1;
        }
        return 1;
chk_error:
        return 0;
}

void mli_MagicId_warn_version(const struct mli_MagicId *magic)
{
        if ((magic->mayor != MLI_VERSION_MAYOR) ||
            (magic->minor != MLI_VERSION_MINOR) ||
            (magic->patch != MLI_VERSION_PATCH)) {
                fprintf(stderr,
                        "[WARNING] Expected magic version to be %d.%d.%d, "
                        "but it is %d.%d.%d.\n",
                        MLI_VERSION_MAYOR,
                        MLI_VERSION_MINOR,
                        MLI_VERSION_PATCH,
                        magic->mayor,
                        magic->minor,
                        magic->patch);
        }
}

/* map */
/* --- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

MLI_VECTOR_IMPLEMENTATION(mli_MapItemVector, struct mliMapItem)

struct mli_Map mli_Map_init(void)
{
        struct mli_Map map;
        map.items = mli_MapItemVector_init();
        return map;
}

void mli_Map_free(struct mli_Map *map)
{
        size_t i;
        for (i = 0; i < map->items.size; i++) {
                struct mliMapItem *item = &map->items.array[i];
                mli_String_free(&item->key);
        }
        mli_MapItemVector_free(&map->items);
}

int mli_Map_malloc(struct mli_Map *map)
{
        chk_msg(mli_MapItemVector_malloc(&map->items, 0u),
                "Failed to malloc map items.");
        return 1;
chk_error:
        return 0;
}

uint64_t mli_Map_size(const struct mli_Map *map) { return map->items.size; }

int mli_Map_find(
        const struct mli_Map *map,
        const struct mli_String *key,
        uint64_t *idx)
{
        uint64_t i;
        for (i = 0; i < map->items.size; i++) {
                struct mliMapItem *item = &map->items.array[i];
                if (mli_String_compare(&item->key, key) == 0) {
                        (*idx) = i;
                        return 1;
                }
        }
        return 0;
}

int mli_Map_has(const struct mli_Map *map, const struct mli_String *key)
{
        uint64_t idx;
        return mli_Map_find(map, key, &idx);
}

int mli_Map_insert(
        struct mli_Map *map,
        const struct mli_String *key,
        uint64_t value)
{
        struct mliMapItem item;
        if (map->items.size == 0u) {
                chk_msg(mli_Map_malloc(map),
                        "Failed to initially malloc dyn-map.");
        }
        chk_msg(0 == mli_Map_has(map, key), "Key already in use.");
        item.key = mli_String_init();
        chk_msg(mli_String_copy(&item.key, key),
                "Failed to malloc map item key.") item.value = value;
        chk_msg(mli_MapItemVector_push_back(&map->items, item),
                "Failed to mmaloc item.");
        return 1;
chk_error:
        return 0;
}

int mli_Map_get(
        const struct mli_Map *map,
        const struct mli_String *key,
        uint64_t *value)
{
        uint64_t idx;
        struct mliMapItem *item = NULL;
        chk_msg(mli_Map_find(map, key, &idx), "Key does not exist.");
        item = &map->items.array[idx];
        (*value) = item->value;
        return 1;
chk_error:
        return 0;
}

/* map_json */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Map_insert_key_from_json(
        struct mli_Map *map,
        const struct mli_Json *json,
        const uint64_t token,
        const uint64_t value)
{
        struct mli_String buff = mli_String_init();
        const uint64_t name_strlen =
                (json->tokens[token].end - json->tokens[token].start);
        chk_msg(json->tokens[token].type == JSMN_STRING,
                "Expected key to be of type string.");
        chk_msg(mli_String_malloc(&buff, name_strlen),
                "Can not malloc String.");
        buff.size = buff.capacity;
        chk_msg(mli_Json_cstr_by_token(
                        json, token, buff.array, buff.capacity + 1),
                "Failed to extract string from json.");
        chk_msg(mli_Map_insert(map, &buff, value),
                "Failed to insert name and value into map.");
        mli_String_free(&buff);
        return 1;
chk_error:
        mli_String_free(&buff);
        mli_Json_debug_token_fprint(stderr, json, token);
        return 0;
}

int mli_Map_get_value_for_string_from_json(
        const struct mli_Map *map,
        const struct mli_Json *json,
        const uint64_t token,
        uint32_t *out_value)
{
        struct mli_String buff = mli_String_init();
        uint64_t value;
        uint64_t name_strlen =
                (json->tokens[token].end - json->tokens[token].start);
        chk_msg(json->tokens[token].type == JSMN_STRING,
                "Expected token to be of type string to be given to mliMap.");
        chk_msg(mli_String_malloc(&buff, name_strlen),
                "Can not malloc String.");
        chk_msg(mli_Json_cstr_by_token(
                        json, token, buff.array, buff.capacity + 1),
                "Failed to extract string from json.");
        buff.size = buff.capacity;
        chk_msg(mli_Map_get(map, &buff, &value),
                "Failed to get value for json-string-key from map.");
        (*out_value) = (uint32_t)value;

        mli_String_free(&buff);
        return 1;
chk_error:
        mli_String_free(&buff);
        mli_Json_debug_token_fprint(stderr, json, token);
        return 0;
}

/* mat */
/* --- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mli_Mat_set(struct mli_Mat *a, uint64_t col, uint64_t row, const double v)
{
        switch (col) {
        case 0:
                switch (row) {
                case 0:
                        a->r00 = v;
                        break;
                case 1:
                        a->r01 = v;
                        break;
                case 2:
                        a->r02 = v;
                        break;
                default:
                        assert(0);
                        break;
                }
                break;
        case 1:
                switch (row) {
                case 0:
                        a->r10 = v;
                        break;
                case 1:
                        a->r11 = v;
                        break;
                case 2:
                        a->r12 = v;
                        break;
                default:
                        assert(0);
                        break;
                }
                break;
        case 2:
                switch (row) {
                case 0:
                        a->r20 = v;
                        break;
                case 1:
                        a->r21 = v;
                        break;
                case 2:
                        a->r22 = v;
                        break;
                default:
                        assert(0);
                        break;
                }
                break;
        default:
                assert(0);
                break;
        }
}

double mli_Mat_get(const struct mli_Mat *a, uint64_t col, uint64_t row)
{
        double o = 0;
        switch (col) {
        case 0:
                switch (row) {
                case 0:
                        o = a->r00;
                        break;
                case 1:
                        o = a->r01;
                        break;
                case 2:
                        o = a->r02;
                        break;
                default:
                        assert(0);
                        break;
                }
                break;
        case 1:
                switch (row) {
                case 0:
                        o = a->r10;
                        break;
                case 1:
                        o = a->r11;
                        break;
                case 2:
                        o = a->r12;
                        break;
                default:
                        assert(0);
                        break;
                }
                break;
        case 2:
                switch (row) {
                case 0:
                        o = a->r20;
                        break;
                case 1:
                        o = a->r21;
                        break;
                case 2:
                        o = a->r22;
                        break;
                default:
                        assert(0);
                        break;
                }
                break;
        default:
                assert(0);
                break;
        }
        return o;
}

struct mli_Mat mli_Mat_unity(void)
{
        struct mli_Mat u;
        u.r00 = 1.0;
        u.r01 = 0.0;
        u.r02 = 0.0;
        u.r10 = 0.0;
        u.r11 = 1.0;
        u.r12 = 0.0;
        u.r20 = 0.0;
        u.r21 = 0.0;
        u.r22 = 1.0;
        return u;
}

int mli_Mat_equal_margin(
        const struct mli_Mat a,
        const struct mli_Mat b,
        const double margin)
{
        int i, j;
        for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                        const double diff = fabs(
                                mli_Mat_get(&a, i, j) - mli_Mat_get(&b, i, j));
                        if (diff > margin) {
                                return 0;
                        }
                }
        }
        return 1;
}

struct mli_Mat mli_Mat_init_tait_bryan(
        const double rx,
        const double ry,
        const double rz)
{
        struct mli_Mat rot;
        const double cosRx = cos(rx);
        const double cosRy = cos(ry);
        const double cosRz = cos(rz);
        const double sinRx = sin(rx);
        const double sinRy = sin(ry);
        const double sinRz = sin(rz);
        rot.r00 = cosRy * cosRz;
        rot.r01 = cosRx * sinRz + sinRx * sinRy * cosRz;
        rot.r02 = sinRx * sinRz - cosRx * sinRy * cosRz;
        rot.r10 = -cosRy * sinRz;
        rot.r11 = cosRx * cosRz - sinRx * sinRy * sinRz;
        rot.r12 = sinRx * cosRz + cosRx * sinRy * sinRz;
        rot.r20 = sinRy;
        rot.r21 = -sinRx * cosRy;
        rot.r22 = cosRx * cosRy;
        return rot;
}

struct mli_Mat mli_Mat_init_axis_angle(
        const struct mli_Vec axis,
        const double angle)
{
        struct mli_Mat rot;
        const double norm = mli_Vec_norm(axis);
        const double rx = axis.x / norm;
        const double ry = axis.y / norm;
        const double rz = axis.z / norm;
        const double sinR = sin(-angle);
        const double cosR = cos(-angle);
        rot.r00 = cosR + rx * rx * (1. - cosR);
        rot.r01 = rx * ry * (1. - cosR) - rz * sinR;
        rot.r02 = rx * rz * (1. - cosR) + ry * sinR;
        rot.r10 = ry * rx * (1. - cosR) + rz * sinR;
        rot.r11 = cosR + ry * ry * (1. - cosR);
        rot.r12 = ry * rz * (1. - cosR) - rx * sinR;
        rot.r20 = rz * rx * (1. - cosR) - ry * sinR;
        rot.r21 = rz * ry * (1. - cosR) + rx * sinR;
        rot.r22 = cosR + rz * rz * (1. - cosR);
        return rot;
}

struct mli_Mat mli_Mat_init_columns(
        const struct mli_Vec c0,
        const struct mli_Vec c1,
        const struct mli_Vec c2)
{
        struct mli_Mat m;
        m.r00 = c0.x;
        m.r01 = c1.x;
        m.r02 = c2.x;
        m.r10 = c0.y;
        m.r11 = c1.y;
        m.r12 = c2.y;
        m.r20 = c0.z;
        m.r21 = c1.z;
        m.r22 = c2.z;
        return m;
}

struct mli_Mat mli_Mat_covariance(
        const struct mli_Vec *vecs,
        const uint64_t num_vecs,
        const struct mli_Vec vecs_mean)
{
        /*
        Estimates the is the 'sample covariance'.

        Parameters
        ----------
        vecs
                The sample of vectors to estimate the covariance of.
        num_vecs
                Number of vectors in sample.
        vecs_mean
                The mean of the vectors.
        */
        uint64_t i;

        double sum_xx = 0.0;
        double sum_xy = 0.0;
        double sum_xz = 0.0;

        double sum_yx = 0.0;
        double sum_yy = 0.0;
        double sum_yz = 0.0;

        double sum_zx = 0.0;
        double sum_zy = 0.0;
        double sum_zz = 0.0;

        double sum_x = 0.0;
        double sum_y = 0.0;
        double sum_z = 0.0;

        double inum = 1.0 / (((double)num_vecs) - 1.0);
        struct mli_Vec v;
        struct mli_Mat cov;

        for (i = 0; i < num_vecs; i++) {
                v = mli_Vec_substract(vecs[i], vecs_mean);
                sum_xx += v.x * v.x;
                sum_xy += v.x * v.y;
                sum_xz += v.x * v.z;

                /* sum_yx += v.y * v.x; */
                sum_yy += v.y * v.y;
                sum_yz += v.y * v.z;

                /* sum_zx += v.z * v.x; */
                /* sum_zy += v.z * v.y; */
                sum_zz += v.z * v.z;

                sum_x += v.x;
                sum_y += v.y;
                sum_z += v.z;
        }

        sum_yx = sum_xy;
        sum_zx = sum_xz;
        sum_zy = sum_yz;

        cov.r00 = (sum_xx * inum) - (sum_x * inum) * (sum_x * inum);
        cov.r01 = (sum_xy * inum) - (sum_x * inum) * (sum_y * inum);
        cov.r02 = (sum_xz * inum) - (sum_x * inum) * (sum_z * inum);
        cov.r10 = (sum_yx * inum) - (sum_y * inum) * (sum_x * inum);
        cov.r11 = (sum_yy * inum) - (sum_y * inum) * (sum_y * inum);
        cov.r12 = (sum_yz * inum) - (sum_y * inum) * (sum_z * inum);
        cov.r20 = (sum_zx * inum) - (sum_z * inum) * (sum_x * inum);
        cov.r21 = (sum_zy * inum) - (sum_z * inum) * (sum_y * inum);
        cov.r22 = (sum_zz * inum) - (sum_z * inum) * (sum_z * inum);
        return cov;
}

struct mli_Mat mli_Mat_transpose(const struct mli_Mat m)
{
        struct mli_Mat t;
        t.r00 = m.r00;
        t.r01 = m.r10;
        t.r02 = m.r20;

        t.r10 = m.r01;
        t.r11 = m.r11;
        t.r12 = m.r21;

        t.r20 = m.r02;
        t.r21 = m.r12;
        t.r22 = m.r22;
        return t;
}

struct mli_Mat mli_Mat_multiply(const struct mli_Mat x, const struct mli_Mat y)
{
        struct mli_Mat p;
        p.r00 = x.r00 * y.r00 + x.r01 * y.r10 + x.r02 * y.r20;
        p.r10 = x.r10 * y.r00 + x.r11 * y.r10 + x.r12 * y.r20;
        p.r20 = x.r20 * y.r00 + x.r21 * y.r10 + x.r22 * y.r20;
        p.r01 = x.r00 * y.r01 + x.r01 * y.r11 + x.r02 * y.r21;
        p.r11 = x.r10 * y.r01 + x.r11 * y.r11 + x.r12 * y.r21;
        p.r21 = x.r20 * y.r01 + x.r21 * y.r11 + x.r22 * y.r21;
        p.r02 = x.r00 * y.r02 + x.r01 * y.r12 + x.r02 * y.r22;
        p.r12 = x.r10 * y.r02 + x.r11 * y.r12 + x.r12 * y.r22;
        p.r22 = x.r20 * y.r02 + x.r21 * y.r12 + x.r22 * y.r22;
        return p;
}

struct mli_Mat mli_Mat_minor(const struct mli_Mat x, const int d)
{
        struct mli_Mat m = x;
        switch (d) {
        case 0:
                break;
        case 1:
                m.r00 = 1.0;

                m.r01 = 0.0;
                m.r02 = 0.0;

                m.r10 = 0.0;
                m.r20 = 0.0;
                break;
        case 2:
                m.r00 = 1.0;
                m.r11 = 1.0;

                m.r01 = 0.0;
                m.r02 = 0.0;

                m.r10 = 0.0;
                m.r20 = 0.0;

                m.r12 = 0.0;
                m.r21 = 0.0;
                break;
        default:
                assert(0);
                break;
        }
        return m;
}

struct mli_Mat mli_Mat_vector_outer_product(const struct mli_Vec v)
{
        struct mli_Mat x;
        x.r00 = -2.0 * v.x * v.x + 1;
        x.r01 = -2.0 * v.x * v.y;
        x.r02 = -2.0 * v.x * v.z;
        x.r10 = -2.0 * v.y * v.x;
        x.r11 = -2.0 * v.y * v.y + 1;
        x.r12 = -2.0 * v.y * v.z;
        x.r20 = -2.0 * v.z * v.x;
        x.r21 = -2.0 * v.z * v.y;
        x.r22 = -2.0 * v.z * v.z + 1;
        return x;
}

void mli_Mat_qr_decompose(
        const struct mli_Mat m,
        struct mli_Mat *q,
        struct mli_Mat *r)
{
        /* housholder */
        /* ========== */

        struct mli_Mat q0, q1;
        struct mli_Mat z = m;

        struct mli_Vec e;
        struct mli_Vec x;
        double xnorm;

        /* column 0 */
        /* -------- */
        z = mli_Mat_minor(z, 0);
        x = mli_Vec_init(z.r00, z.r10, z.r20);
        xnorm = mli_Vec_norm(x);
        if (m.r00 > 0) {
                xnorm = -xnorm;
        };

        e = mli_Vec_init(xnorm, 0.0, 0.0);
        e = mli_Vec_add(x, e);
        e = mli_Vec_normalized(e);

        q0 = mli_Mat_vector_outer_product(e);
        z = mli_Mat_multiply(q0, z);

        /* column 1 */
        /* -------- */
        z = mli_Mat_minor(z, 1);
        x = mli_Vec_init(z.r01, z.r11, z.r21);
        xnorm = mli_Vec_norm(x);
        if (m.r11 > 0) {
                xnorm = -xnorm;
        };

        e = mli_Vec_init(0.0, xnorm, 0.0);
        e = mli_Vec_add(x, e);
        e = mli_Vec_normalized(e);

        q1 = mli_Mat_vector_outer_product(e);
        z = mli_Mat_multiply(q1, z);

        /* finalize */
        /* ======== */
        (*q) = q0;
        (*r) = mli_Mat_multiply(q0, m);

        (*q) = mli_Mat_multiply(q1, *q);

        (*r) = mli_Mat_multiply(*q, m);
        (*q) = mli_Mat_transpose(*q);
}

int mli_Mat_has_shurform(const struct mli_Mat m, const double margin)
{
        if (fabs(m.r10) > margin)
                return 0;
        if (fabs(m.r20) > margin)
                return 0;
        if (fabs(m.r21) > margin)
                return 0;
        return 1;
}

void mli_Mat_find_eigenvalues(
        const struct mli_Mat a,
        double *e0,
        double *e1,
        double *e2,
        const double margin,
        const uint64_t max_num_iterations)
{
        struct mli_Mat ak = a;
        struct mli_Mat qq = mli_Mat_unity();
        uint64_t k = 0;

        while (k < max_num_iterations) {
                struct mli_Mat q, r;
                k += 1;
                mli_Mat_qr_decompose(ak, &q, &r);
                ak = mli_Mat_multiply(r, q);
                qq = mli_Mat_multiply(qq, q);

                if (mli_Mat_has_shurform(ak, margin)) {
                        break;
                }
        }
        (*e0) = ak.r00;
        (*e1) = ak.r11;
        (*e2) = ak.r22;
}

int mli_Mat_find_eigenvector_for_eigenvalue(
        struct mli_Mat A,
        const double eigen_value,
        struct mli_Vec *eigen_vector,
        const double tolerance)
{
        int pivots[4] = {0, 0, 0, 0};
        struct mli_Vec right_hand_side = {1.0, 1.0, 1.0};

        chk_msg(tolerance > 0.0, "Expected tolerance > 0.0");

        A.r00 -= eigen_value;
        A.r11 -= eigen_value;
        A.r22 -= eigen_value;

        chk_msg(mli_Mat_lup_decompose(&A, pivots, tolerance),
                "Can not decompose LU-matices");
        mli_Mat_lup_solve(&A, pivots, &right_hand_side, eigen_vector);
        (*eigen_vector) = mli_Vec_normalized(*eigen_vector);

        return 1;
chk_error:
        return 0;
}

int mli_Mat_lup_decompose(struct mli_Mat *A, int *P, const double tolerance)
{
        /*
        Parameters
        ----------
        A
                Matrix to decompose.
        tolerance
                small tolerance number to detect failure when the
                matrix is near degenerate
        Output
        ------
        A
                Matrix A is changed, it contains a copy of both matrices
                L-E and U as A=(L-E)+U such that P*A=L*U.
        P : int[N+1]
                The permutation matrix is not stored as a matrix,
                but in an integer vector P of size N+1 containing column
                indexes where the permutation matrix has "1".
                The last element P[N]=S+N,
                where S is the number of row exchanges needed for determinant
                computation, det(P)=(-1)^S

        Inspired by https://en.wikipedia.org/wiki/LU_decomposition
        */
        int i, j, k, imax;
        double maxA, absA;

        for (i = 0; i <= 3; i++) {
                /* Unit permutation matrix, P[N] initialized with N */
                P[i] = i;
        }

        for (i = 0; i < 3; i++) {
                maxA = 0.0;
                imax = i;

                for (k = i; k < 3; k++) {
                        absA = fabs(mli_Mat_get(A, k, i));
                        if (absA > maxA) {
                                maxA = absA;
                                imax = k;
                        }
                }

                if (maxA < tolerance) {
                        return 0; /* failure, matrix is degenerate */
                }

                if (imax != i) {
                        /* pivoting */
                        double tmp[3];
                        j = P[i];
                        P[i] = P[imax];
                        P[imax] = j;

                        /* pivoting rows of A */
                        tmp[0] = mli_Mat_get(A, i, 0);
                        tmp[1] = mli_Mat_get(A, i, 1);
                        tmp[2] = mli_Mat_get(A, i, 2);

                        mli_Mat_set(A, i, 0, mli_Mat_get(A, imax, 0));
                        mli_Mat_set(A, i, 1, mli_Mat_get(A, imax, 1));
                        mli_Mat_set(A, i, 2, mli_Mat_get(A, imax, 2));

                        mli_Mat_set(A, imax, 0, tmp[0]);
                        mli_Mat_set(A, imax, 1, tmp[1]);
                        mli_Mat_set(A, imax, 2, tmp[2]);

                        /* counting pivots starting from 3 (for determinant) */
                        P[3]++;
                }

                for (j = i + 1; j < 3; j++) {
                        double tmp =
                                mli_Mat_get(A, j, i) / mli_Mat_get(A, i, i);
                        mli_Mat_set(A, j, i, tmp);

                        for (k = i + 1; k < 3; k++) {
                                double tmp =
                                        (mli_Mat_get(A, j, k) -
                                         mli_Mat_get(A, j, i) *
                                                 mli_Mat_get(A, i, k));
                                mli_Mat_set(A, j, k, tmp);
                        }
                }
        }

        return 1; /* decomposition done */
}

void mli_Mat_lup_solve(
        const struct mli_Mat *A,
        const int *P,
        const struct mli_Vec *b,
        struct mli_Vec *x)
{
        /*
        Solve linear equation A*x = b for x.

        Parameters
        ----------
        A and P
                Filled in mli_Mat_lup_decompose.
        b
                Right hand side vector b in A*x = b.

        Output
        ------
        x
                The solution vector x of A*x = b
        */
        const uint64_t N = 3;
        uint64_t i, k, idown;
        for (i = 0; i < N; i++) {
                mli_Vec_set(x, i, mli_Vec_get(b, P[i]));
                for (k = 0; k < i; k++) {
                        const double tmp =
                                (mli_Vec_get(x, i) -
                                 mli_Mat_get(A, i, k) * mli_Vec_get(x, k));
                        mli_Vec_set(x, i, tmp);
                }
        }

        i = N - 1;
        for (idown = 0; idown < N; idown++) {
                double tmp;
                for (k = i + 1; k < N; k++) {
                        const double tmp =
                                (mli_Vec_get(x, i) -
                                 mli_Mat_get(A, i, k) * mli_Vec_get(x, k));
                        mli_Vec_set(x, i, tmp);
                }
                tmp = mli_Vec_get(x, i) / mli_Mat_get(A, i, i);
                mli_Vec_set(x, i, tmp);
                i--;
        }
}

struct mli_Vec mli_Mat_dot_product(
        const struct mli_Mat *m,
        const struct mli_Vec v)
{
        struct mli_Vec o;
        o.x = m->r00 * v.x + m->r10 * v.y + m->r20 * v.z;
        o.y = m->r01 * v.x + m->r11 * v.y + m->r21 * v.z;
        o.z = m->r02 * v.x + m->r12 * v.y + m->r22 * v.z;
        return o;
}

/* materials */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_MaterialsCapacity mli_MaterialsCapacity_init(void)
{
        struct mli_MaterialsCapacity cap;
        cap.num_spectra = 0;
        cap.num_surfaces = 0;
        cap.num_media = 0;
        cap.num_boundary_layers = 0;
        return cap;
}

struct mli_Materials mli_Materials_init(void)
{
        struct mli_Materials res;
        res.spectra = mli_SpectrumArray_init();
        res.surfaces = mli_SurfaceArray_init();
        res.media = mli_MediumArray_init();
        res.boundary_layers = mli_BoundaryLayerArray_init();

        res.default_medium = 0u;
        return res;
}

void mli_Materials_free(struct mli_Materials *self)
{
        mli_SpectrumArray_free(&self->spectra);
        mli_SurfaceArray_free(&self->surfaces);
        mli_MediumArray_free(&self->media);
        mli_BoundaryLayerArray_free(&self->boundary_layers);
        (*self) = mli_Materials_init();
}

int mli_Materials_malloc(
        struct mli_Materials *self,
        const struct mli_MaterialsCapacity rescap)
{
        uint64_t i;
        mli_Materials_free(self);

        chk(mli_SpectrumArray_malloc(&self->spectra, rescap.num_spectra));
        for (i = 0; i < self->spectra.size; i++) {
                self->spectra.array[i] = mli_Spectrum_init();
        }

        chk(mli_SurfaceArray_malloc(&self->surfaces, rescap.num_surfaces));
        for (i = 0; i < self->surfaces.size; i++) {
                self->surfaces.array[i] = mli_Surface_init();
        }

        chk(mli_MediumArray_malloc(&self->media, rescap.num_media));
        for (i = 0; i < self->media.size; i++) {
                self->media.array[i] = mli_Medium_init();
        }

        chk(mli_BoundaryLayerArray_malloc(
                &self->boundary_layers, rescap.num_boundary_layers));
        for (i = 0; i < self->boundary_layers.size; i++) {
                self->boundary_layers.array[i] = mli_BoundaryLayer_init();
        }

        return 1;
chk_error:
        mli_Materials_free(self);
        return 0;
}

int mli_Materials_info_fprint(FILE *f, const struct mli_Materials *self)
{
        uint32_t i = 0;
        struct mli_String tmp = mli_String_init();

        fprintf(f, "materials\n");
        fprintf(f, "---------\n");
        fprintf(f, "\n");

        fprintf(f, "    media\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        fprintf(f, "    ");
        fprintf(f, "%3s ", "#");
        fprintf(f, "%24s ", "name");
        fprintf(f, "%12s ", "default");
        fprintf(f, "\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        for (i = 0; i < self->media.size; i++) {
                struct mli_Medium *medium = &self->media.array[i];
                fprintf(f, "    ");
                fprintf(f, "% 3d ", i);
                fprintf(f, "%24s ", medium->name.array);

                if (i == self->default_medium) {
                        fprintf(f, "%12s", "True");
                }
                fprintf(f, "\n");
        }
        fprintf(f, "\n");

        fprintf(f, "    surfaces\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        for (i = 0; i < self->surfaces.size; i++) {
                struct mli_Surface *surface = &self->surfaces.array[i];
                struct mli_String tmp = mli_String_init();
                chk(mli_Surface_type_to_string(surface->type, &tmp));
                fprintf(f, "    ");
                fprintf(f, "% 3d ", i);
                fprintf(f, "%24s ", surface->name.array);
                fprintf(f, "%24s ", tmp.array);
                fprintf(f, "\n");
        }
        fprintf(f, "\n");

        fprintf(f, "    boundary layers\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        for (i = 0; i < self->boundary_layers.size; i++) {
                struct mli_BoundaryLayer *layer =
                        &self->boundary_layers.array[i];
                fprintf(f, "    ");
                fprintf(f, "% 3d ", i);
                fprintf(f, "%24s ", layer->name.array);

                fprintf(f, "%8ld ", layer->inner.medium);
                fprintf(f, "%8ld ", layer->inner.surface);

                fprintf(f, "%8ld ", layer->outer.medium);
                fprintf(f, "%8ld ", layer->outer.surface);

                fprintf(f, "\n");
        }

        mli_String_free(&tmp);
        return 1;
chk_error:
        mli_String_free(&tmp);
        return 0;
}

int mli_Materials__has_surface_name_cstr(
        const struct mli_Materials *self,
        const char *name)
{
        uint64_t i;
        for (i = 0; i < self->surfaces.size; i++) {
                if (mli_String_equal_cstr(
                            &self->surfaces.array[i].name, name)) {
                        return 1;
                }
        }
        return 0;
}

int mli_Materials__has_medium_name_cstr(
        const struct mli_Materials *self,
        const char *name)
{
        uint64_t i;
        for (i = 0; i < self->media.size; i++) {
                if (mli_String_equal_cstr(&self->media.array[i].name, name)) {
                        return 1;
                }
        }
        return 0;
}

/* materials_equal */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Materials_media_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b)
{
        uint64_t i = 0u;
        chk_msg(a->media.size == b->media.size, "Different number of media.");
        for (i = 0; i < a->media.size; i++) {
                chk_msg(mli_Medium_equal(
                                &a->media.array[i], &b->media.array[i]),
                        "Medium is different.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In materials.media[%lu].\n", i);
        return 0;
}

int mli_Materials_surfaces_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b)
{
        uint64_t i = 0u;
        chk_msg(a->surfaces.size == b->surfaces.size,
                "Different number of surfaces.");
        for (i = 0; i < a->surfaces.size; i++) {
                chk_msg(mli_Surface_equal(
                                &a->surfaces.array[i], &b->surfaces.array[i]),
                        "Surface is different.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In materials.surfaces[%lu].\n", i);
        return 0;
}

int mli_Materials_boundary_layers_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b)
{
        uint64_t i = 0u;
        chk_msg(a->boundary_layers.size == b->boundary_layers.size,
                "Different number of boundary_layers.");
        for (i = 0; i < a->boundary_layers.size; i++) {
                chk_msg(mli_BoundaryLayer_equal(
                                &a->boundary_layers.array[i],
                                &b->boundary_layers.array[i]),
                        "Boundary layer is different.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In materials.boundary_layers[%lu].\n", i);
        return 0;
}

int mli_Materials_spectra_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b)
{
        uint64_t i = 0u;
        chk_msg(a->spectra.size == b->spectra.size,
                "Different number of spectra.");
        for (i = 0; i < a->spectra.size; i++) {
                chk_msg(mli_Spectrum_equal(
                                &a->spectra.array[i], &b->spectra.array[i]),
                        "Spectrum is not equal.");
        }
        return 1;
chk_error:
        return 0;
}

int mli_Materials_default_medium_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b)
{
        chk_msg(a->default_medium == b->default_medium,
                "Different default_medium.");
        return 1;
chk_error:
        return 0;
}

int mli_Materials_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b)
{
        chk(mli_Materials_spectra_equal(a, b));
        chk(mli_Materials_media_equal(a, b));
        chk(mli_Materials_surfaces_equal(a, b));
        chk(mli_Materials_boundary_layers_equal(a, b));
        chk(mli_Materials_default_medium_equal(a, b));
        return 1;
chk_error:
        return 0;
}

/* materials_from_archive */
/* ---------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Materials__key_from_filename(
        struct mli_String *key,
        const struct mli_String *filename)
{
        struct mli_String basename = mli_String_init();
        struct mli_String extension = mli_String_init();
        chk(mli_path_basename(filename, &basename));
        chk(mli_path_splitext(&basename, key, &extension));
        mli_String_free(&extension);
        mli_String_free(&basename);
        return 1;
chk_error:
        return 0;
}

int mli_Materials_from_Archive__set_spectra(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive)
{
        uint64_t spc_idx = 0;
        uint64_t arc_idx = 0;

        struct mli_String *filename = NULL;
        struct mli_String key = mli_String_init();

        spc_idx = 0u;
        for (arc_idx = 0u; arc_idx < mli_Archive_size(archive); arc_idx++) {
                filename = &archive->filenames.items.array[arc_idx].key;

                if (mli_String_starts_with_cstr(
                            filename, "materials/spectra/") &&
                    mli_String_ends_with_cstr(filename, ".csv")) {
                        struct mli_String *payload =
                                &archive->textfiles.array[arc_idx];
                        struct mli_IO buff = mli_IO_init();

                        chk_msg(spc_idx < materials->spectra.size,
                                "Expected sufficient capacity for spectra.");

                        chk(mli_FuncInfo_malloc(
                                &materials->spectra.array[spc_idx].info));
                        chk(mli_IO_open_memory(&buff));
                        chk(mli_IO_text_write_String(&buff, payload));
                        mli_IO_rewind(&buff);
                        chk_msg(mli_Func_from_csv(
                                        &materials->spectra.array[spc_idx]
                                                 .spectrum,
                                        &materials->spectra.array[spc_idx]
                                                 .info.x,
                                        &materials->spectra.array[spc_idx]
                                                 .info.y,
                                        &buff),
                                "Failed to parse spectral function from "
                                "archive.");
                        mli_IO_close(&buff);

                        chk(mli_Materials__key_from_filename(&key, filename));

                        chk_msg(mli_Map_insert(&names->spectra, &key, spc_idx),
                                "Failed to insert spectrum-name into map.");

                        chk(mli_String_copy(
                                &materials->spectra.array[spc_idx].name,
                                &names->spectra.items.array[spc_idx].key));
                        spc_idx += 1u;
                }
        }

        mli_String_free(&key);

        return 1;
chk_error:
        return 0;
}

int mli_Materials_from_Archive__set_media(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive)
{
        uint64_t med_idx = 0;
        uint64_t arc_idx = 0;
        struct mli_String *filename = NULL;
        struct mli_String key = mli_String_init();

        med_idx = 0u;
        for (arc_idx = 0u; arc_idx < mli_Archive_size(archive); arc_idx++) {
                filename = &archive->filenames.items.array[arc_idx].key;

                if (mli_String_starts_with_cstr(filename, "materials/media/") &&
                    mli_String_ends_with_cstr(filename, ".json")) {
                        struct mli_String *payload =
                                &archive->textfiles.array[arc_idx];

                        chk_msg(med_idx < materials->media.size,
                                "Expected sufficient capacity for media.");

                        chk(mli_Materials__key_from_filename(&key, filename));

                        chk_msg(mli_Medium_from_json_string_and_name(
                                        &materials->media.array[med_idx],
                                        &names->spectra,
                                        payload,
                                        &key),
                                "Can't parse medium from json string.");

                        chk_msg(mli_Map_insert(&names->media, &key, med_idx),
                                "Failed to insert media-name into map.");

                        med_idx += 1u;
                }
        }

        chk_msgf(
                med_idx == materials->media.size,
                ("Expected to parse %lu media but only found %lu.",
                 materials->media.size,
                 med_idx));

        mli_String_free(&key);

        return 1;
chk_error:
        return 0;
}

int mli_Materials_from_Archive__set_surfaces(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive)
{
        uint64_t srf_idx = 0;
        uint64_t arc_idx = 0;
        struct mli_String *filename = NULL;
        struct mli_String key = mli_String_init();

        srf_idx = 0u;
        for (arc_idx = 0u; arc_idx < mli_Archive_size(archive); arc_idx++) {
                filename = &archive->filenames.items.array[arc_idx].key;
                if (mli_String_starts_with_cstr(
                            filename, "materials/surfaces/") &&
                    mli_String_ends_with_cstr(filename, ".json")) {
                        struct mli_String *payload =
                                &archive->textfiles.array[arc_idx];

                        chk_msg(srf_idx < materials->surfaces.size,
                                "Expected sufficient capacity for surfaces.");

                        chk(mli_Materials__key_from_filename(&key, filename));

                        chk_msg(mli_Surface_from_json_string_and_name(
                                        &materials->surfaces.array[srf_idx],
                                        &names->spectra,
                                        payload,
                                        &key),
                                "Can't parse surface from json string.");

                        chk_msg(mli_Map_insert(&names->surfaces, &key, srf_idx),
                                "Failed to insert surface-name into map.");

                        srf_idx += 1u;
                }
        }

        chk_msgf(
                srf_idx == materials->surfaces.size,
                ("Expected to parse %lu surfaces but only found %lu.",
                 materials->surfaces.size,
                 srf_idx));

        mli_String_free(&key);

        return 1;
chk_error:
        return 0;
}

int mli_BoundaryLayer_from_json_string_and_name(
        struct mli_BoundaryLayer *self,
        const struct mli_Map *surface_names,
        const struct mli_Map *media_names,
        const struct mli_String *json_string,
        const struct mli_String *name)
{
        struct mli_String key = mli_String_init();
        struct mli_Json json = mli_Json_init();
        struct mli_JsonWalk walk = mli_JsonWalk_init();

        chk_msg(mli_String_copy(&self->name, name),
                "Failed to copy boundary layer name.");
        chk_msg(mli_Json_from_string(&json, json_string),
                "Failed to parse boundary layer json string.");

        walk = mli_JsonWalk_set(&json);
        chk(mli_JsonWalk_to_key(&walk, "inner"));
        chk(mli_JsonWalk_to_key(&walk, "medium"));
        chk(mli_JsonWalk_get_string(&walk, &key));
        chk_msgf(
                mli_Map_get(media_names, &key, &self->inner.medium),
                ("Expected 'inner->medium':'%s' to be in media_names.",
                 key.array));

        walk = mli_JsonWalk_set(&json);
        chk(mli_JsonWalk_to_key(&walk, "inner"));
        chk(mli_JsonWalk_to_key(&walk, "surface"));
        chk(mli_JsonWalk_get_string(&walk, &key));
        chk_msgf(
                mli_Map_get(surface_names, &key, &self->inner.surface),
                ("Expected 'inner->surface':'%s' to be in surface_names.",
                 key.array));

        walk = mli_JsonWalk_set(&json);
        chk(mli_JsonWalk_to_key(&walk, "outer"));
        chk(mli_JsonWalk_to_key(&walk, "medium"));
        chk(mli_JsonWalk_get_string(&walk, &key));
        chk_msgf(
                mli_Map_get(media_names, &key, &self->outer.medium),
                ("Expected 'outer->medium':'%s' to be in media_names.",
                 key.array));

        walk = mli_JsonWalk_set(&json);
        chk(mli_JsonWalk_to_key(&walk, "outer"));
        chk(mli_JsonWalk_to_key(&walk, "surface"));
        chk(mli_JsonWalk_get_string(&walk, &key));
        chk_msgf(
                mli_Map_get(surface_names, &key, &self->outer.surface),
                ("Expected 'outer->surface':'%s' to be in surface_names.",
                 key.array));

        mli_Json_free(&json);
        mli_String_free(&key);
        return 1;
chk_error:
        return 0;
}

int mli_Materials_from_Archive__set_boundary_layers(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive)
{
        uint64_t bdl_idx = 0;
        uint64_t arc_idx = 0;
        struct mli_String *filename = NULL;
        struct mli_String key = mli_String_init();

        for (arc_idx = 0u; arc_idx < mli_Archive_size(archive); arc_idx++) {
                filename = &archive->filenames.items.array[arc_idx].key;
                if (mli_String_starts_with_cstr(
                            filename, "materials/boundary_layers/") &&
                    mli_String_ends_with_cstr(filename, ".json")) {
                        struct mli_String *payload =
                                &archive->textfiles.array[arc_idx];

                        chk_msg(bdl_idx < materials->boundary_layers.size,
                                "Expected sufficient capacity for boundary "
                                "layers.");

                        chk(mli_Materials__key_from_filename(&key, filename));

                        chk_msg(mli_BoundaryLayer_from_json_string_and_name(
                                        &materials->boundary_layers
                                                 .array[bdl_idx],
                                        &names->surfaces,
                                        &names->media,
                                        payload,
                                        &key),
                                "Can't set boundary layer from json string.");

                        chk_msg(mli_Map_insert(
                                        &names->boundary_layers, &key, bdl_idx),
                                "Failed to insert boundary layer name into "
                                "map.");

                        bdl_idx += 1u;
                }
        }

        chk_msgf(
                bdl_idx == materials->boundary_layers.size,
                ("Expected to parse %lu boundary layers but only found %lu.",
                 materials->boundary_layers.size,
                 bdl_idx));

        mli_String_free(&key);

        return 1;
chk_error:
        return 0;
}

int mli_Materials_from_Archive__set_default_medium(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive)
{

        struct mli_String filename = mli_String_init();
        struct mli_String key = mli_String_init();
        struct mli_String *payload = NULL;

        chk(mli_String_from_cstr(&filename, "materials/default_medium.txt"));
        chk_msg(mli_Archive_get(archive, &filename, &payload),
                "Can not find 'materials/default_medium.txt' in scenery.");

        chk(mli_String_strip(payload, &key));
        chk_msg(mli_Map_get(&names->media, &key, &materials->default_medium),
                "Failed to assign the 'default_medium'.");

        mli_String_free(&filename);
        mli_String_free(&key);

        return 1;
chk_error:
        return 0;
}

int mli_Materials_from_Archive(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive)
{
        uint64_t iii = 0;
        struct mli_IO ioerr = mli_IO_init();
        struct mli_MaterialsCapacity capacity = mli_MaterialsCapacity_init();

        /* free */
        /* ---- */
        mli_Materials_free(materials);
        mli_materials_Names_free(names);
        chk(mli_materials_Names_malloc(names));

        /* estimate capacity */
        /* ----------------- */
        capacity.num_boundary_layers = mli_Archive_num_filename_prefix_sufix(
                archive, "materials/boundary_layers/", ".json");
        capacity.num_media = mli_Archive_num_filename_prefix_sufix(
                archive, "materials/media/", ".json");
        capacity.num_surfaces = mli_Archive_num_filename_prefix_sufix(
                archive, "materials/surfaces/", ".json");
        capacity.num_spectra = mli_Archive_num_filename_prefix_sufix(
                archive, "materials/spectra/", ".csv");
        chk_msg(mli_Materials_malloc(materials, capacity),
                "Can't malloc materials.");

        /* set fields */
        /* ---------- */
        chk_msg(mli_Materials_from_Archive__set_spectra(
                        materials, names, archive),
                "Can't set spectra from archive.");
        chk_msg(mli_Materials_from_Archive__set_media(
                        materials, names, archive),
                "Can't set media from archive.");
        chk_msg(mli_Materials_from_Archive__set_surfaces(
                        materials, names, archive),
                "Can't set surfaces from archive.");
        chk_msg(mli_Materials_from_Archive__set_boundary_layers(
                        materials, names, archive),
                "Can't set boundary layers from archive.");
        chk_msg(mli_Materials_from_Archive__set_default_medium(
                        materials, names, archive),
                "Can't set default_medium from archive.");

        chk(mli_IO_adopt_file(&ioerr, stderr));
        chk(mli_IO_text_write_cstr_format(&ioerr, "materials->spectra\n"));
        chk(mli_IO_text_write_cstr_format(&ioerr, "------------------\n"));
        for (iii = 0; iii < materials->spectra.size; iii++) {
                struct mli_Spectrum *spectrum = &materials->spectra.array[iii];
                chk(mli_Spectrum_print_to_io(spectrum, &ioerr));
        }

        return 1;
chk_error:
        mli_Materials_free(materials);
        mli_materials_Names_free(names);
        return 0;
}

/* materials_names */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_materials_Names mli_materials_Names_init(void)
{
        struct mli_materials_Names nm;
        nm.spectra = mli_Map_init();
        nm.media = mli_Map_init();
        nm.surfaces = mli_Map_init();
        nm.boundary_layers = mli_Map_init();
        return nm;
}

int mli_materials_Names_malloc(struct mli_materials_Names *namemap)
{
        mli_materials_Names_free(namemap);
        chk_mem(mli_Map_malloc(&namemap->spectra));
        chk_mem(mli_Map_malloc(&namemap->media));
        chk_mem(mli_Map_malloc(&namemap->surfaces));
        chk_mem(mli_Map_malloc(&namemap->boundary_layers));
        return 1;
chk_error:
        return 0;
}

void mli_materials_Names_free(struct mli_materials_Names *namemap)
{
        mli_Map_free(&namemap->spectra);
        mli_Map_free(&namemap->media);
        mli_Map_free(&namemap->surfaces);
        mli_Map_free(&namemap->boundary_layers);
}

/* materials_serialize */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Materials_to_io(const struct mli_Materials *self, struct mli_IO *f)
{
        uint64_t i;
        struct mli_MagicId magic = mli_MagicId_init();
        chk(mli_MagicId_set(&magic, "mli_Materials"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        chk_IO_write(&self->spectra.size, sizeof(uint64_t), 1u, f);
        chk_IO_write(&self->surfaces.size, sizeof(uint64_t), 1u, f);
        chk_IO_write(&self->media.size, sizeof(uint64_t), 1u, f);
        chk_IO_write(&self->boundary_layers.size, sizeof(uint64_t), 1u, f);

        for (i = 0; i < self->spectra.size; i++) {
                chk(mli_Spectrum_to_io(&self->spectra.array[i], f));
        }
        for (i = 0; i < self->surfaces.size; i++) {
                chk(mli_Surface_to_io(&self->surfaces.array[i], f));
        }
        for (i = 0; i < self->media.size; i++) {
                chk(mli_Medium_to_io(&self->media.array[i], f));
        }
        for (i = 0; i < self->boundary_layers.size; i++) {
                chk(mli_BoundaryLayer_to_io(
                        &self->boundary_layers.array[i], f));
        }

        chk_IO_write(&self->default_medium, sizeof(uint64_t), 1, f);

        return 1;
chk_error:
        return 0;
}

int mli_Materials_from_io(struct mli_Materials *self, struct mli_IO *f)
{
        uint64_t i;
        struct mli_MagicId magic;
        struct mli_MaterialsCapacity cap;

        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Materials"));
        mli_MagicId_warn_version(&magic);

        chk_IO_read(&cap.num_spectra, sizeof(uint64_t), 1u, f);
        chk_IO_read(&cap.num_surfaces, sizeof(uint64_t), 1u, f);
        chk_IO_read(&cap.num_media, sizeof(uint64_t), 1u, f);
        chk_IO_read(&cap.num_boundary_layers, sizeof(uint64_t), 1u, f);

        chk(mli_Materials_malloc(self, cap));

        for (i = 0; i < self->spectra.size; i++) {
                chk(mli_Spectrum_from_io(&self->spectra.array[i], f));
        }
        for (i = 0; i < self->surfaces.size; i++) {
                chk(mli_Surface_from_io(&self->surfaces.array[i], f));
        }
        for (i = 0; i < self->media.size; i++) {
                chk(mli_Medium_from_io(&self->media.array[i], f));
        }
        for (i = 0; i < self->boundary_layers.size; i++) {
                chk(mli_BoundaryLayer_from_io(
                        &self->boundary_layers.array[i], f));
        }
        chk_IO_read(&self->default_medium, sizeof(uint64_t), 1, f);

        return 1;
chk_error:
        return 0;
}

/* materials_valid */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Materials_valid_media(const struct mli_Materials *self)
{
        uint64_t i = 0u;
        for (i = 0; i < self->media.size; i++) {
                chk_msg(mli_Medium_valid_wrt_materials(
                                &self->media.array[i], self),
                        "Medium is not valid.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In materials.media[%lu]\n", i);
        return 0;
}

int mli_Materials_valid_spectra(const struct mli_Materials *self)
{
        uint64_t i = 0u;
        for (i = 0; i < self->spectra.size; i++) {
                chk_msg(mli_Func_is_valid(&self->spectra.array[i].spectrum),
                        "Expected spectrum function to be valid.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In materials.spectra[%lu]\n", i);
        return 0;
}

int mli_Materials_valid_surfaces(const struct mli_Materials *self)
{
        uint64_t i = 0u;
        for (i = 0; i < self->surfaces.size; i++) {
                chk(self);
                chk_warning("IMPLEMENT ME!!!");
        }
        return 1;
chk_error:
        fprintf(stderr, "In materials.surface[%lu]\n", i);
        return 0;
}

int mli_Materials_valid_boundary_layers(const struct mli_Materials *self)
{
        uint64_t i = 0u;
        for (i = 0; i < self->boundary_layers.size; i++) {
                struct mli_BoundaryLayer *layer =
                        &self->boundary_layers.array[i];
                chk_msg(layer->inner.surface < self->surfaces.size,
                        "inner.surface is invalid.");
                chk_msg(layer->outer.surface < self->surfaces.size,
                        "outer.surface is invalid.");
                chk_msg(layer->inner.medium < self->media.size,
                        "inner.medium is invalid.");
                chk_msg(layer->outer.medium < self->media.size,
                        "outer.medium is invalid.");

                chk_msg(mli_String_valid(&layer->name, 1),
                        "boundary_layer.name is invalid.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In materials.boundary_layers[%lu]\n", i);
        return 0;
}

int mli_Materials_valid(const struct mli_Materials *self)
{
        chk_msg(self->default_medium <= self->media.size,
                "Expected default-medium to reference a valid medium.");
        chk(mli_Materials_valid_spectra(self));
        chk(mli_Materials_valid_media(self));
        chk(mli_Materials_valid_surfaces(self));
        chk(mli_Materials_valid_boundary_layers(self));
        return 1;
chk_error:
        return 0;
}

/* math */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

double mli_math_rad2deg(const double angle_in_rad)
{
        return 180. * angle_in_rad / MLI_MATH_PI;
}

double mli_math_deg2rad(const double angle_in_deg)
{
        return angle_in_deg * (1. / 180.) * MLI_MATH_PI;
}

double mli_math_hypot(const double a, const double b)
{
        return sqrt(a * a + b * b);
}

double mli_math_square(const double a) { return a * a; }

/*
 *  parameters
 *  ----------
 *      points          Sorted array in ascending order.
 *      num_points      Number of points.
 *      point_arg       The point to find the upper-bound for.
 */
uint64_t mli_math_upper_compare_double(
        const double *points,
        const uint64_t num_points,
        const double point_arg)
{
        uint64_t upper_index = 0;
        MLI_MATH_UPPER_COMPARE(points, num_points, point_arg, upper_index);
        return upper_index;
}

void mli_math_histogram(
        const double *bin_edges,
        const uint64_t num_bin_edges,
        uint64_t *underflow_bin,
        uint64_t *bins,
        uint64_t *overflow_bin,
        const double point)
{
        uint64_t idx_upper =
                mli_math_upper_compare_double(bin_edges, num_bin_edges, point);
        if (idx_upper == 0) {
                (*underflow_bin) += 1u;
        } else if (idx_upper == num_bin_edges) {
                (*overflow_bin) += 1u;
        } else {
                bins[idx_upper - 1] += 1u;
        }
}

void mli_math_linspace(
        const double start,
        const double stop,
        double *points,
        const uint64_t num_points)
{
        uint64_t i;
        const double range = stop - start;
        const double step = range / (double)(num_points - 1u);
        for (i = 0; i < num_points; i++) {
                points[i] = (double)i * step + start;
        }
}

double mli_math_mean(const double vals[], const uint64_t size)
{
        uint64_t i;
        double sum = 0;
        for (i = 0; i < size; i++) {
                sum = sum + vals[i];
        }
        return sum / (double)size;
}

double mli_math_std(
        const double vals[],
        const uint64_t size,
        const double vals_mean)
{
        uint64_t i;
        double s = 0.;
        for (i = 0; i < size; i++) {
                s = s + (vals[i] - vals_mean) * (vals[i] - vals_mean);
        }
        return sqrt(s / (double)size);
}

double mli_math_bin_center_in_linear_space(
        const double start,
        const double stop,
        const uint64_t num_bins,
        const uint64_t bin)
{
        const double width = stop - start;
        const double bin_width = width / (double)num_bins;
        return start + bin * bin_width + 0.5 * bin_width;
}

double mli_math_linear_interpolate_1d(
        const double weight,
        const double start,
        const double end)
{
        return start + weight * (end - start);
}

double mli_math_linear_interpolate_2d(
        const double xarg,
        const double x0,
        const double y0,
        const double x1,
        const double y1)
{
        /*
         *      |
         *  y1 -|            o
         *      |
         *  y0 -|    o
         *      |       xarg
         *      +----|---|---|----
         *          x0       x1
         *
         *  f(x) = m*x + b
         *  m = (y1 - y0)/(x1 - x0)
         *  y0 = m*x0 + b
         *  b = y0 - m*x0
         */
        const double m = (y1 - y0) / (x1 - x0);
        const double b = y0 - m * x0;
        return m * xarg + b;
}

double mli_math_relative_ratio(const double a, const double b)
{
        return fabs(a - b) / (0.5 * (a + b));
}

double mli_math_interpret_int64_as_double(int64_t i)
{
        double f;
        memcpy(&f, &i, sizeof(double));
        return f;
}

int64_t mli_math_interpret_double_as_int64(double d)
{
        int64_t i;
        memcpy(&i, &d, sizeof(int64_t));
        return i;
}

/* math_quadratic_equation */
/* ----------------------- */

/* Copyright 2019 Sebastian A. Mueller */

int mli_math_quadratic_equation(
        const double p,
        const double q,
        double *minus_solution,
        double *plus_solution)
{
        /*
         *  y = a*x^2 + b*x + c
         *  p = b/a
         *  q = c/a
         *  x_m = -p/2 - sqrt((-p/2)^2 - q)
         *  x_p = -p/2 + sqrt((-p/2)^2 - q)
         */
        const double p_over_2 = 0.5 * p;
        const double inner_part_of_squareroot = p_over_2 * p_over_2 - q;
        double squareroot;
        if (inner_part_of_squareroot >= 0.0) {
                squareroot = sqrt(inner_part_of_squareroot);
                (*minus_solution) = -p_over_2 - squareroot;
                (*plus_solution) = -p_over_2 + squareroot;
                return 1;
        } else {
                return 0;
        }
}

/* medium */
/* ------ */

/* Copyright 2018-2024 Sebastian Achim Mueller */

struct mli_Medium mli_Medium_init(void)
{
        struct mli_Medium out;
        out.name = mli_String_init();
        out.refraction_spectrum = 0;
        out.absorbtion_spectrum = 0;
        return out;
}

void mli_Medium_free(struct mli_Medium *self)
{
        mli_String_free(&self->name);
        (*self) = mli_Medium_init();
}

int mli_Medium_valid_wrt_materials(
        const struct mli_Medium *self,
        const struct mli_Materials *materials)
{
        chk_msg(mli_String_valid(&self->name, 1), "name is invalid.");

        chk_msg(self->refraction_spectrum < materials->spectra.size,
                "refraction_spectrum index is not in materials.");
        chk_msg(self->absorbtion_spectrum < materials->spectra.size,
                "absorbtion_spectrum index is not in materials.");

        return 1;
chk_error:
        return 0;
}

int mli_Medium_equal(const struct mli_Medium *a, const struct mli_Medium *b)
{
        chk_msg(mli_String_equal(&a->name, &b->name),
                "Different names of medium models.");
        chk_msg(a->refraction_spectrum == b->refraction_spectrum,
                "Different refraction_spectrum.");
        chk_msg(a->absorbtion_spectrum == b->absorbtion_spectrum,
                "Different absorbtion_spectrum.");

        return 1;
chk_error:
        return 0;
}

int mli_Medium_to_io(const struct mli_Medium *self, struct mli_IO *f)
{
        struct mli_MagicId magic = mli_MagicId_init();
        chk(mli_MagicId_set(&magic, "mli_Medium"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        chk_msg(mli_String_to_io(&self->name, f),
                "Can't write medium.name to io.");
        chk_IO_write(&self->refraction_spectrum, sizeof(int64_t), 1u, f);
        chk_IO_write(&self->absorbtion_spectrum, sizeof(int64_t), 1u, f);

        return 1;
chk_error:
        return 0;
}

int mli_Medium_from_io(struct mli_Medium *self, struct mli_IO *f)
{
        struct mli_MagicId magic;
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Medium"));
        mli_MagicId_warn_version(&magic);

        chk_msg(mli_String_from_io(&self->name, f),
                "Can't read medium.name from io.");
        chk_IO_read(&self->refraction_spectrum, sizeof(int64_t), 1u, f);
        chk_IO_read(&self->absorbtion_spectrum, sizeof(int64_t), 1u, f);

        return 1;
chk_error:
        return 0;
}

int mli_Medium_from_json_string_and_name(
        struct mli_Medium *self,
        const struct mli_Map *spectra_names,
        const struct mli_String *json_string,
        const struct mli_String *name)
{
        struct mli_Json json = mli_Json_init();
        struct mli_JsonWalk walk = mli_JsonWalk_init();
        struct mli_String key = mli_String_init();

        mli_Medium_free(self);

        chk_msg(mli_String_copy(&self->name, name), "Can't copy medium name.");

        chk_msg(mli_Json_from_string(&json, json_string),
                "Can't parse medium from json string.");

        walk = mli_JsonWalk_set(&json);
        chk_msg(mli_JsonWalk_to_key(&walk, "refraction_spectrum"),
                "Expected field 'refraction_spectrum' in medium json string.");
        chk_msg(mli_JsonWalk_get_string(&walk, &key),
                "Expected 'refraction_spectrum' to hold a string.");
        chk_msg(mli_Map_get(spectra_names, &key, &self->refraction_spectrum),
                "Expected 'refraction_spectrum' to be in spectra_names.");

        walk = mli_JsonWalk_set(&json);
        chk_msg(mli_JsonWalk_to_key(&walk, "absorbtion_spectrum"),
                "Expected field 'absorbtion_spectrum' in medium json string.");
        chk_msg(mli_JsonWalk_get_string(&walk, &key),
                "Expected 'absorbtion_spectrum' to hold a string.");
        chk_msg(mli_Map_get(spectra_names, &key, &self->absorbtion_spectrum),
                "Expected 'absorbtion_spectrum' to be in spectra_names.");

        mli_String_free(&key);
        mli_Json_free(&json);
        return 1;
chk_error:
        return 0;
}

/* medium_array */
/* ------------ */

/* Copyright 2018-2024 Sebastian Achim Mueller */
MLI_ARRAY_IMPLEMENTATION_FREE(
        mli_MediumArray,
        struct mli_Medium,
        mli_Medium_free)

/* object */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Object mli_Object_init(void)
{
        struct mli_Object obj;
        obj.num_vertices = 0;
        obj.vertices = NULL;

        obj.num_vertex_normals = 0;
        obj.vertex_normals = NULL;

        obj.num_faces = 0;
        obj.faces_vertices = NULL;
        obj.faces_vertex_normals = NULL;
        obj.faces_materials = NULL;

        obj.num_materials = 0;
        obj.material_names = NULL;
        return obj;
}

void mli_Object_free(struct mli_Object *obj)
{
        size_t i;
        free(obj->vertices);
        free(obj->vertex_normals);
        free(obj->faces_vertices);
        free(obj->faces_vertex_normals);
        free(obj->faces_materials);

        for (i = 0; i < obj->num_materials; i++) {
                mli_String_free(&obj->material_names[i]);
        }
        free(obj->material_names);
        *obj = mli_Object_init();
}

int mli_Object_malloc(
        struct mli_Object *obj,
        const uint64_t num_vertices,
        const uint64_t num_vertex_normals,
        const uint64_t num_faces,
        const uint64_t num_materials)
{
        uint32_t i;
        chk_msg(num_vertices < UINT32_MAX, "Expected num_vertices < uint32");
        chk_msg(num_vertex_normals < UINT32_MAX,
                "Expected num_vertex_normals < uint32");
        chk_msg(num_faces < UINT32_MAX, "Expected num_faces < uint32");
        chk_msg(num_materials < UINT32_MAX, "Expected num_materials < uint32");
        chk_msg(num_materials > 0, "Expected num_materials > 0");
        chk_msg(num_materials <= num_faces,
                "Expected num_materials <= num_faces");

        mli_Object_free(obj);
        obj->num_vertices = num_vertices;
        obj->num_vertex_normals = num_vertex_normals;
        obj->num_faces = num_faces;
        obj->num_materials = num_materials;
        chk_malloc(obj->vertices, struct mli_Vec, obj->num_vertices);
        chk_malloc(
                obj->vertex_normals, struct mli_Vec, obj->num_vertex_normals);
        chk_malloc(obj->faces_vertices, struct mli_object_Face, obj->num_faces);
        chk_malloc(
                obj->faces_vertex_normals,
                struct mli_object_Face,
                obj->num_faces);
        chk_malloc(obj->faces_materials, uint16_t, obj->num_faces);
        chk_malloc(obj->material_names, struct mli_String, obj->num_materials);
        for (i = 0; i < obj->num_materials; i++) {
                obj->material_names[i] = mli_String_init();
        }
        return 1;
chk_error:
        return 0;
}

int mli_Object_equal(const struct mli_Object *a, const struct mli_Object *b)
{
        uint64_t i;
        chk(a->num_vertices == b->num_vertices);
        chk(a->num_vertex_normals == b->num_vertex_normals);
        chk(a->num_faces == b->num_faces);
        chk(a->num_materials == b->num_materials);

        for (i = 0; i < a->num_vertices; i++) {
                chk(mli_Vec_equal(a->vertices[i], b->vertices[i]));
        }
        for (i = 0; i < a->num_vertex_normals; i++) {
                chk(mli_Vec_equal(a->vertex_normals[i], b->vertex_normals[i]));
        }
        for (i = 0; i < a->num_faces; i++) {
                chk(mli_object_Face_equal(
                        a->faces_vertices[i], b->faces_vertices[i]));
                chk(mli_object_Face_equal(
                        a->faces_vertex_normals[i],
                        b->faces_vertex_normals[i]));
                chk(a->faces_materials[i] == b->faces_materials[i]);
        }
        for (i = 0; i < a->num_materials; i++) {
                chk(mli_String_equal(
                        &a->material_names[i], &b->material_names[i]));
        }
        return 1;
chk_error:
        return 0;
}

uint32_t mli_Object_resolve_material_idx(
        const struct mli_Object *obj,
        const uint32_t face_idx)
{
        assert(face_idx < obj->num_faces);
        return obj->faces_materials[face_idx];
}

/* object_AABB */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Object_face_in_local_frame_has_overlap_aabb(
        const struct mli_Object *obj,
        const uint64_t face_idx,
        const struct mli_AABB aabb)
{
        struct mli_object_Face face;
        if (face_idx >= obj->num_faces) {
                return 0;
        }
        face = obj->faces_vertices[face_idx];
        if (mli_Triangle_has_overlap_aabb(
                    obj->vertices[face.a],
                    obj->vertices[face.b],
                    obj->vertices[face.c],
                    aabb)) {
                return 1;
        }
        return 0;
}

int mli_Object_has_overlap_aabb(
        const struct mli_Object *obj,
        const struct mli_HomTra local2root,
        const struct mli_AABB aabb)
{
        uint64_t face_idx;
        for (face_idx = 0; face_idx < obj->num_faces; face_idx++) {
                const struct mli_object_Face face =
                        obj->faces_vertices[face_idx];

                const struct mli_Vec a_local = obj->vertices[face.a];
                const struct mli_Vec b_local = obj->vertices[face.b];
                const struct mli_Vec c_local = obj->vertices[face.c];

                const struct mli_Vec a_root =
                        mli_HomTraComp_pos(&local2root, a_local);
                const struct mli_Vec b_root =
                        mli_HomTraComp_pos(&local2root, b_local);
                const struct mli_Vec c_root =
                        mli_HomTraComp_pos(&local2root, c_local);

                if (mli_Triangle_has_overlap_aabb(
                            a_root, b_root, c_root, aabb)) {
                        return 1;
                }
        }
        return 0;
}

int mli_Object_face_in_local_frame_has_overlap_aabb_void(
        const void *obj,
        const uint32_t face_idx,
        const struct mli_AABB aabb)
{
        return mli_Object_face_in_local_frame_has_overlap_aabb(
                (const struct mli_Object *)obj, (uint64_t)face_idx, aabb);
}

struct mli_AABB mli_Object_aabb(
        const struct mli_Object *obj,
        const struct mli_HomTra local2root)
{
        struct mli_AABB aabb;
        struct mli_Vec seed_vertex_local;
        struct mli_Vec seed_vertex_root;
        struct mli_object_Face seed_face;
        uint64_t face_idx;

        if (obj->num_faces == 0) {
                aabb.lower =
                        mli_Vec_init(MLI_MATH_NAN, MLI_MATH_NAN, MLI_MATH_NAN);
                aabb.upper =
                        mli_Vec_init(MLI_MATH_NAN, MLI_MATH_NAN, MLI_MATH_NAN);
                return aabb;
        }

        seed_face = obj->faces_vertices[0];
        seed_vertex_local = obj->vertices[seed_face.a];
        seed_vertex_root = mli_HomTraComp_pos(&local2root, seed_vertex_local);

        aabb.lower = seed_vertex_root;
        aabb.upper = seed_vertex_root;
        for (face_idx = 0; face_idx < obj->num_faces; face_idx++) {
                const struct mli_object_Face face =
                        obj->faces_vertices[face_idx];
                const struct mli_Vec a_local = obj->vertices[face.a];
                const struct mli_Vec b_local = obj->vertices[face.b];
                const struct mli_Vec c_local = obj->vertices[face.c];
                const struct mli_Vec a_root =
                        mli_HomTraComp_pos(&local2root, a_local);
                const struct mli_Vec b_root =
                        mli_HomTraComp_pos(&local2root, b_local);
                const struct mli_Vec c_root =
                        mli_HomTraComp_pos(&local2root, c_local);

                aabb.lower.x = MLI_MATH_MIN2(aabb.lower.x, a_root.x);
                aabb.lower.y = MLI_MATH_MIN2(aabb.lower.y, a_root.y);
                aabb.lower.z = MLI_MATH_MIN2(aabb.lower.z, a_root.z);

                aabb.upper.x = MLI_MATH_MAX2(aabb.upper.x, a_root.x);
                aabb.upper.y = MLI_MATH_MAX2(aabb.upper.y, a_root.y);
                aabb.upper.z = MLI_MATH_MAX2(aabb.upper.z, a_root.z);

                aabb.lower.x = MLI_MATH_MIN2(aabb.lower.x, b_root.x);
                aabb.lower.y = MLI_MATH_MIN2(aabb.lower.y, b_root.y);
                aabb.lower.z = MLI_MATH_MIN2(aabb.lower.z, b_root.z);

                aabb.upper.x = MLI_MATH_MAX2(aabb.upper.x, b_root.x);
                aabb.upper.y = MLI_MATH_MAX2(aabb.upper.y, b_root.y);
                aabb.upper.z = MLI_MATH_MAX2(aabb.upper.z, b_root.z);

                aabb.lower.x = MLI_MATH_MIN2(aabb.lower.x, c_root.x);
                aabb.lower.y = MLI_MATH_MIN2(aabb.lower.y, c_root.y);
                aabb.lower.z = MLI_MATH_MIN2(aabb.lower.z, c_root.z);

                aabb.upper.x = MLI_MATH_MAX2(aabb.upper.x, c_root.x);
                aabb.upper.y = MLI_MATH_MAX2(aabb.upper.y, c_root.y);
                aabb.upper.z = MLI_MATH_MAX2(aabb.upper.z, c_root.z);
        }
        return aabb;
}

struct mli_AABB mli_Object_aabb_in_local_frame(const struct mli_Object *obj)
{
        struct mli_HomTra unity;
        unity.translation = mli_Vec_init(0.0, 0.0, 0.0);
        unity.rotation = mli_Mat_unity();
        return mli_Object_aabb(obj, unity);
}

/* object_face */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_object_Face mli_object_Face_set(
        const uint32_t a,
        const uint32_t b,
        const uint32_t c)
{
        struct mli_object_Face face;
        face.a = a;
        face.b = b;
        face.c = c;
        return face;
}

int mli_object_Face_equal(
        const struct mli_object_Face a,
        const struct mli_object_Face b)
{
        if (a.a != b.a)
                return 0;
        if (a.b != b.b)
                return 0;
        if (a.c != b.c)
                return 0;
        return 1;
}

/* object_face_vector */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_VECTOR_IMPLEMENTATION(mli_object_FaceVector, struct mli_object_Face)

/* object_serialize */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Object_to_io(const struct mli_Object *obj, struct mli_IO *f)
{
        uint64_t i;
        struct mli_MagicId magic;
        chk(mli_MagicId_set(&magic, "mli_Object"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        chk_IO_write(&obj->num_vertices, sizeof(uint32_t), 1u, f);
        chk_IO_write(&obj->num_vertex_normals, sizeof(uint32_t), 1u, f);
        chk_IO_write(&obj->num_faces, sizeof(uint32_t), 1u, f);
        chk_IO_write(&obj->num_materials, sizeof(uint32_t), 1u, f);

        chk_IO_write(
                obj->vertices, sizeof(struct mli_Vec), obj->num_vertices, f);

        chk_IO_write(
                obj->vertex_normals,
                sizeof(struct mli_Vec),
                obj->num_vertex_normals,
                f);

        chk_IO_write(
                obj->faces_vertices,
                sizeof(struct mli_object_Face),
                obj->num_faces,
                f);
        chk_IO_write(
                obj->faces_vertex_normals,
                sizeof(struct mli_object_Face),
                obj->num_faces,
                f);
        chk_IO_write(obj->faces_materials, sizeof(uint16_t), obj->num_faces, f);

        for (i = 0; i < obj->num_materials; i++) {
                chk(mli_String_to_io(&obj->material_names[i], f));
        }

        return 1;
chk_error:
        return 0;
}

int mli_Object_from_io(struct mli_Object *obj, struct mli_IO *f)
{
        uint64_t i;
        uint32_t num_vertices;
        uint32_t num_vertex_normals;
        uint32_t num_faces;
        uint32_t num_materials;
        struct mli_MagicId magic;
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Object"));
        mli_MagicId_warn_version(&magic);

        chk_IO_read(&num_vertices, sizeof(uint32_t), 1u, f);
        chk_IO_read(&num_vertex_normals, sizeof(uint32_t), 1u, f);
        chk_IO_read(&num_faces, sizeof(uint32_t), 1u, f);
        chk_IO_read(&num_materials, sizeof(uint32_t), 1u, f);

        chk(mli_Object_malloc(
                obj,
                num_vertices,
                num_vertex_normals,
                num_faces,
                num_materials));

        chk_IO_read(
                obj->vertices, sizeof(struct mli_Vec), obj->num_vertices, f);
        chk_IO_read(
                obj->vertex_normals,
                sizeof(struct mli_Vec),
                obj->num_vertex_normals,
                f);

        chk_IO_read(
                obj->faces_vertices,
                sizeof(struct mli_object_Face),
                obj->num_faces,
                f);
        chk_IO_read(
                obj->faces_vertex_normals,
                sizeof(struct mli_object_Face),
                obj->num_faces,
                f);
        chk_IO_read(obj->faces_materials, sizeof(uint16_t), obj->num_faces, f);

        for (i = 0; i < obj->num_materials; i++) {
                chk(mli_String_from_io(&obj->material_names[i], f));
        }

        chk_msg(mli_Object_has_valid_faces(obj),
                "A face refers to a not existing vertex/vertex_normal.");
        return 1;
chk_error:
        mli_Object_free(obj);
        return 0;
}

/* object_valid */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Object_is_valid(const struct mli_Object *obj)
{
        chk_msg(mli_Object_has_valid_vertices(obj), "Bad vertex.");
        chk_msg(mli_Object_has_valid_faces(obj), "Expected faces to be valid.");
        chk_msg(mli_Object_has_valid_normals(obj, MLI_MATH_EPSILON),
                "Bad vertex-normal.");
        chk_msg(mli_Object_has_valid_materials(obj), "Bad material.");
        return 1;
chk_error:
        return 0;
}

int mli_Object_has_valid_faces(const struct mli_Object *obj)
{
        uint32_t i = 0;
        for (i = 0; i < obj->num_faces; i++) {
                chk_msg(obj->faces_vertices[i].a <= obj->num_vertices,
                        "Expected face.a <= num_vertices");
                chk_msg(obj->faces_vertices[i].b <= obj->num_vertices,
                        "Expected face.b <= num_vertices");
                chk_msg(obj->faces_vertices[i].c <= obj->num_vertices,
                        "Expected face.c <= num_vertices");

                chk_msg(obj->faces_vertex_normals[i].a <=
                                obj->num_vertex_normals,
                        "Expected faces_vertex_normals.a <= "
                        "num_vertex_normals");
                chk_msg(obj->faces_vertex_normals[i].b <=
                                obj->num_vertex_normals,
                        "Expected faces_vertex_normals.b <= "
                        "num_vertex_normals");
                chk_msg(obj->faces_vertex_normals[i].c <=
                                obj->num_vertex_normals,
                        "Expected faces_vertex_normals.c <= "
                        "num_vertex_normals");
                chk_msg(obj->faces_materials[i] < obj->num_materials,
                        "Expected faces_materials < "
                        "num_materials");
        }
        return 1;
chk_error:
        fprintf(stderr, "In obj.faces[%u]\n", i);
        return 0;
}

int mli_Object_has_valid_vertices(const struct mli_Object *obj)
{
        uint32_t i = 0;
        for (i = 0; i < obj->num_vertices; i++) {
                chk_msg(!MLI_MATH_IS_NAN(obj->vertices[i].x), "X is 'nan'.");
                chk_msg(!MLI_MATH_IS_NAN(obj->vertices[i].y), "Y is 'nan'.");
                chk_msg(!MLI_MATH_IS_NAN(obj->vertices[i].z), "Z is 'nan'.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In obj.vertices[%u]\n", i);
        return 0;
}

int mli_Object_has_valid_normals(
        const struct mli_Object *obj,
        const double epsilon)
{
        uint32_t i = 0;
        for (i = 0; i < obj->num_vertex_normals; i++) {
                double norm;
                chk_msg(!MLI_MATH_IS_NAN(obj->vertex_normals[i].x),
                        "X is 'nan'.");
                chk_msg(!MLI_MATH_IS_NAN(obj->vertex_normals[i].y),
                        "Y is 'nan'.");
                chk_msg(!MLI_MATH_IS_NAN(obj->vertex_normals[i].z),
                        "Z is 'nan'.");

                norm = mli_Vec_norm(obj->vertex_normals[i]);
                chk_msg(fabs(norm - 1.0) <= epsilon,
                        "Expected vertex_normals to be normalized.");
        }
        return 1;
chk_error:
        fprintf(stderr, "In obj.vertex_normals[%u]\n", i);
        return 0;
}

int mli_Object_has_valid_materials(const struct mli_Object *obj)
{
        uint32_t i = 0;
        for (i = 0; i < obj->num_materials; i++) {
                chk(mli_String_valid(&obj->material_names[i], 1));
        }
        return 1;
chk_error:
        fprintf(stderr, "In obj.material_names[%u]\n", i);
        return 0;
}

int mli_Object_num_unused(
        const struct mli_Object *obj,
        uint32_t *num_unused_vertices,
        uint32_t *num_unused_vertex_normals,
        uint32_t *num_unused_materials)
{
        uint64_t *v = NULL;
        uint64_t *vn = NULL;
        uint64_t *mtl = NULL;
        uint64_t i, f;

        (*num_unused_vertices) = 0;
        (*num_unused_vertex_normals) = 0;
        (*num_unused_materials) = 0;

        chk_malloc(v, uint64_t, obj->num_vertices);
        MLI_MATH_ARRAY_SET(v, 0, obj->num_vertices);

        chk_malloc(vn, uint64_t, obj->num_vertex_normals);
        MLI_MATH_ARRAY_SET(vn, 0, obj->num_vertex_normals);

        chk_malloc(mtl, uint64_t, obj->num_materials);
        MLI_MATH_ARRAY_SET(mtl, 0, obj->num_materials);

        for (f = 0; f < obj->num_faces; f++) {
                v[obj->faces_vertices[f].a] += 1;
                v[obj->faces_vertices[f].b] += 1;
                v[obj->faces_vertices[f].c] += 1;

                vn[obj->faces_vertex_normals[f].a] += 1;
                vn[obj->faces_vertex_normals[f].b] += 1;
                vn[obj->faces_vertex_normals[f].c] += 1;

                mtl[obj->faces_materials[f]] += 1;
        }

        for (i = 0; i < obj->num_vertices; i++) {
                if (v[i] == 0) {
                        (*num_unused_vertices) += 1;
                }
        }

        for (i = 0; i < obj->num_vertex_normals; i++) {
                if (vn[i] == 0) {
                        (*num_unused_vertex_normals) += 1;
                }
        }

        for (i = 0; i < obj->num_materials; i++) {
                if (mtl[i] == 0) {
                        (*num_unused_materials) += 1;
                }
        }

        free(v);
        free(vn);
        free(mtl);
        return 1;
chk_error:
        return 0;
}

/* object_wavefront */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

#define MLI_WAVEFRONT_FACE_LINE_V 7
#define MLI_WAVEFRONT_FACE_LINE_V_VN 37
#define MLI_WAVEFRONT_FACE_LINE_V_VT_VN 25

#define MLI_WAVEFRONT_LINE_BUFF_LENGTH 64

int mli_Object_is_face_line_toggle(const int state)
{
        switch (state) {
        case 3:
        case 5:
        case 7:

        case 27:
        case 29:
        case 32:
        case 34:
        case 37:

        case 11:
        case 13:
        case 15:
        case 17:
        case 19:
        case 21:
        case 23:
        case 25:
                return 1;
                break;
        }
        return 0;
}

int mli_String_to_uint32(uint32_t *out, const struct mli_String *str)
{
        uint64_t u = 0;
        chk(mli_String_to_uint64(&u, str, 10));
        (*out) = (uint32_t)u;
        return 1;
chk_error:
        return 0;
}

int mli_Object_parse_face_line(
        const struct mli_String *line,
        struct mli_object_Face *faces_vertices,
        struct mli_object_Face *faces_texture_points,
        struct mli_object_Face *faces_vertex_normals,
        int *line_mode)
{
        struct mli_object_Face *v = faces_vertices;
        struct mli_object_Face *vt = faces_texture_points;
        struct mli_object_Face *vn = faces_vertex_normals;
        /*
        statemachine
        ============

        ( 0)
         |
         | 'f'
         V
        ( 1)
         |
         | ws
         V  ___
        ( 2)<__] ws
         |
         | dig
         V  ___
        ( 3)<__] dig
         |\________________________
         | ws                      | '/'
         V  ___                    V
        ( 4)<__] ws               (10)
         |                         |\___________________________
         | dig                     | dig                        | '/'
         V  ___                    V  ___                       V
        ( 5)<__] dig              (11)<__] dig                 (26)
         |                         |                            |
         | ws                      | '/'                        | dig
         V  ___                    V                            V  ___
        ( 6)<__] ws               (12)                         (27)<__] dig
         |                         |                            |
         | dig                     | dig                        | ws
         V  ___                    V  ___                       V  ___
        ( 7)<__] dig              (13)<__] dig                 (28)<__] ws
         |                         |                            |
         | !dig OR ws              | ws                         | dig
         V                         V  ___                       V  ___
        (END)                     (14)<__] ws                  (29)<__] dig
                                   |                            |
        f v v v                    | dig                        | '/'
                                   V  ___                       V
                                  (15)<__] dig                 (30)
                                   |                            |
                                   | '/'                        | '/'
                                   V                            V
                                  (16)                         (31)
                                   |                            |
                                   | dig                        | dig
                                   V  ___                       V  ___
                                  (17)<__] dig                 (32)<__] dig
                                   |                            |
                                   | '/'                        | ws
                                   V                            V  ___
                                  (18)                         (33)<__] ws
                                   |                            |
                                   | dig                        | dig
                                   V  ___                       V  ___
                                  (19)<__] dig                 (34)<__] dig
                                   |                            |
                                   | ws                         | '/'
                                   V  ___                       V
                                  (20)<__] ws                  (35)
                                   |                            |
                                   | dig                        | '/'
                                   V  ___                       V
                                  (21)<__] dig                 (36)
                                   |                            |
                                   | '/'                        | dig
                                   V                            V  ___
                                  (22)                         (37)<__] dig
                                   |                            |
                                   | dig                        | !dig OR ws
                                   V  ___                       V
                                  (23)<__] dig                 (END)
                                   |
                                   | '/'                 f v//vn v//vn v//vn
                                   V
                                  (24)
                                   |
                                   | dig
                                   V  ___
                                  (25)<__] dig
                                   |
                                   | !dig OR ws
                                   V
                                  (END)

                                f v/vt/vn v/vt/vn v/vt/vn

        */
        const int final_state = 99;
        const int error_state = -1;
        int statemachine[][5] = {
                /*       'f'  ws  dig !dig '/' */
                {01, -1, -1, -1, -1}, /* 00 */

                {-1, 02, -1, -1, -1}, /* 01 */
                {-1, 02, 03, -1, -1}, /* 02 */
                {-1, 04, 03, -1, 10}, /* 03 */
                {-1, 04, 05, -1, -1}, /* 04 */
                {-1, 06, 05, -1, -1}, /* 05 */
                {-1, 06, 07, -1, -1}, /* 06 */
                {-1, 99, 07, 99, -1}, /* 07 */

                {-1, -1, -1, -1, -1}, /* 08 */
                {-1, -1, -1, -1, -1}, /* 09 */

                {-1, -1, 11, -1, 26}, /* 10 */
                {-1, -1, 11, -1, 12}, /* 11 */
                {-1, -1, 13, -1, -1}, /* 12 */
                {-1, 14, 13, -1, -1}, /* 13 */
                {-1, 14, 15, -1, -1}, /* 14 */
                {-1, -1, 15, -1, 16}, /* 15 */
                {-1, -1, 17, -1, -1}, /* 16 */
                {-1, -1, 17, -1, 18}, /* 17 */
                {-1, -1, 19, -1, -1}, /* 18 */
                {-1, 20, 19, -1, -1}, /* 19 */
                {-1, 20, 21, -1, -1}, /* 20 */
                {-1, -1, 21, -1, 22}, /* 21 */
                {-1, -1, 23, -1, -1}, /* 22 */
                {-1, -1, 23, -1, 24}, /* 23 */
                {-1, -1, 25, -1, -1}, /* 24 */
                {-1, 99, 25, 99, -1}, /* 25 */

                {-1, -1, 27, -1, -1}, /* 26 */
                {-1, 28, 27, -1, -1}, /* 27 */
                {-1, 28, 29, -1, -1}, /* 28 */
                {-1, -1, 29, -1, 30}, /* 29 */
                {-1, -1, -1, -1, 31}, /* 30 */
                {-1, -1, 32, -1, -1}, /* 31 */
                {-1, 33, 32, -1, -1}, /* 32 */
                {-1, 33, 34, -1, -1}, /* 33 */
                {-1, -1, 34, -1, 35}, /* 34 */
                {-1, -1, -1, -1, 36}, /* 35 */
                {-1, -1, 37, -1, -1}, /* 36 */
                {-1, 99, 37, 99, -1}, /* 37 */
        };

        int state = 0;
        int old_state = state;
        const uint64_t MAX_NUM_CHARS = 256;
        uint64_t i = 0;
        char c;

        struct mli_String wuff = mli_String_init();
        chk(mli_String_malloc(&wuff, MAX_NUM_CHARS));
        wuff.size = 0;

        while (state != final_state) {
                chk_msg(i <= MAX_NUM_CHARS, "Expected less chars in line.");
                if (i < line->size) {
                        chk_msg(mli_String_get(line, i, &c), "Beyond 'line'.");
                } else {
                        c = '\0';
                }

                /* printf("i: %d, c:%c state: %d \n", i, c, state); */

                if (state == error_state) {
                        *line_mode = -1;
                        chk_bad("Can not parse line.");
                }

                /* next state */
                /* ========== */
                if (c == 'f') {
                        state = statemachine[state][0];
                } else if (c == ' ') {
                        state = statemachine[state][1];
                } else if (c == '/') {
                        state = statemachine[state][4];
                } else if (isdigit(c)) {
                        state = statemachine[state][2];
                } else if (!isdigit(c)) {
                        state = statemachine[state][3];
                } else {
                        state = error_state;
                }

                /* line mode */
                /* ========= */
                if (state == MLI_WAVEFRONT_FACE_LINE_V) {
                        *line_mode = MLI_WAVEFRONT_FACE_LINE_V;
                }
                if (state == MLI_WAVEFRONT_FACE_LINE_V_VN) {
                        *line_mode = MLI_WAVEFRONT_FACE_LINE_V_VN;
                }
                if (state == MLI_WAVEFRONT_FACE_LINE_V_VT_VN) {
                        *line_mode = MLI_WAVEFRONT_FACE_LINE_V_VT_VN;
                }

                if (mli_Object_is_face_line_toggle(state)) {
                        wuff.array[wuff.size] = c;
                        wuff.size++;
                } else if (mli_Object_is_face_line_toggle(old_state)) {
                        uint64_t r;
                        for (r = wuff.size; r < MAX_NUM_CHARS; r++) {
                                wuff.array[r] = '\0';
                        }

                        switch (old_state) {
                        /* MLI_WAVEFRONT_FACE_LINE_V                          */
                        case 3:
                                chk(mli_String_to_uint32(&v->a, &wuff));
                                break;
                        case 5:
                                chk(mli_String_to_uint32(&v->b, &wuff));
                                break;
                        case 7:
                                chk(mli_String_to_uint32(&v->c, &wuff));
                                break;

                        /* MLI_WAVEFRONT_FACE_LINE_V_VN                       */
                        /*    3                        v->a                   */
                        case 27:
                                chk(mli_String_to_uint32(&vn->a, &wuff));
                                break;
                        case 29:
                                chk(mli_String_to_uint32(&v->b, &wuff));
                                break;
                        case 32:
                                chk(mli_String_to_uint32(&vn->b, &wuff));
                                break;
                        case 34:
                                chk(mli_String_to_uint32(&v->c, &wuff));
                                break;
                        case 37:
                                chk(mli_String_to_uint32(&vn->c, &wuff));
                                break;

                        /* MLI_WAVEFRONT_FACE_LINE_V_VT_VN                    */
                        /*    3                        v->a                   */
                        case 11:
                                chk(mli_String_to_uint32(&vt->a, &wuff));
                                break;
                        case 13:
                                chk(mli_String_to_uint32(&vn->a, &wuff));
                                break;
                        case 15:
                                chk(mli_String_to_uint32(&v->b, &wuff));
                                break;
                        case 17:
                                chk(mli_String_to_uint32(&vt->b, &wuff));
                                break;
                        case 19:
                                chk(mli_String_to_uint32(&vn->b, &wuff));
                                break;
                        case 21:
                                chk(mli_String_to_uint32(&v->c, &wuff));
                                break;
                        case 23:
                                chk(mli_String_to_uint32(&vt->c, &wuff));
                                break;
                        case 25:
                                chk(mli_String_to_uint32(&vn->c, &wuff));
                                break;

                        default:
                                break;
                        }
                        wuff.size = 0;
                }

                old_state = state;
                i++;
        }
        mli_String_free(&wuff);
        return 1;
chk_error:
        mli_String_free(&wuff);
        return 0;
}

int mli_Object_is_vert_line_toggle(const int state)
{
        switch (state) {
        case 2:
        case 4:
        case 6:
                return 1;
                break;
        }
        return 0;
}

int mli_Object_parse_three_float_line(
        const struct mli_String *line,
        struct mli_Vec *v)
{
        /*
        statemachine
        ============

        ( 0)
         |
         | ws
         V  ___
        ( 1)<__] ws
         |
         | print
         V  ___
        ( 2)<__] print
         |
         | ws
         V  ___
        ( 3)<__] ws
         |
         | print
         V  ___
        ( 4)<__] print
         |
         | ws
         V  ___
        ( 5)<__] ws
         |
         | print
         V  ___
        ( 6)<__] print
         |
         | ws OR eol
         V
        (END)

        */
        const int final_state = 99;
        const int error_state = -1;
        int statemachine[][3] = {
                /* ws print eol */
                {01, -1, -1}, /* 0 */
                {01, 02, -1}, /* 1 */
                {03, 02, -1}, /* 2 */
                {03, 04, -1}, /* 3 */
                {05, 04, -1}, /* 4 */
                {05, 06, -1}, /* 5 */
                {99, 06, 99}, /* 6 */
        };

        int state = 0;
        int old_state = state;
        const uint64_t MAX_NUM_CHARS = 256;
        uint64_t i = 0;
        char c;

        struct mli_String wuff = mli_String_init();
        chk(mli_String_malloc(&wuff, MAX_NUM_CHARS));
        wuff.size = 0;

        while (state != final_state) {
                chk_msg(i <= MAX_NUM_CHARS, "Expected less chars in line.");
                if (i < line->size) {
                        chk_msg(mli_String_get(line, i, &c), "Beyond 'line'.");
                } else {
                        c = '\0';
                }

                if (state == error_state) {
                        fprintf(stderr,
                                "[ERROR] Can not parse line '%s'\n",
                                line->array);
                }
                chk_msg(state != error_state,
                        "Can not parse three float line.");

                /* next state */
                /* ========== */
                if (c == ' ') {
                        state = statemachine[state][0];
                } else if (isprint(c)) {
                        state = statemachine[state][1];
                } else if (c == '\0') {
                        state = statemachine[state][2];
                } else {
                        state = error_state;
                }

                if (mli_Object_is_vert_line_toggle(state)) {
                        wuff.array[wuff.size] = c;
                        wuff.size++;
                } else if (mli_Object_is_vert_line_toggle(old_state)) {
                        uint64_t r;
                        for (r = wuff.size; r < MAX_NUM_CHARS; r++) {
                                wuff.array[r] = '\0';
                        }

                        switch (old_state) {
                        case 2:
                                chk(mli_String_to_double(&v->x, &wuff));
                                break;
                        case 4:
                                chk(mli_String_to_double(&v->y, &wuff));
                                break;
                        case 6:
                                chk(mli_String_to_double(&v->z, &wuff));
                                break;
                        default:
                                break;
                        }
                        wuff.size = 0;
                }

                old_state = state;
                i++;
        }
        return 1;
chk_error:
        return 0;
}

int mli_Object_parse_face_vertices_and_normals(
        const struct mli_String *line,
        struct mli_object_Face *fv,
        struct mli_object_Face *fvn)
{
        int line_mode = -1;
        struct mli_object_Face tmp_fvt;

        chk_msg(mli_Object_parse_face_line(line, fv, &tmp_fvt, fvn, &line_mode),
                "Can not parse face-line.");

        chk_msg((line_mode == MLI_WAVEFRONT_FACE_LINE_V_VT_VN) ||
                        (line_mode == MLI_WAVEFRONT_FACE_LINE_V_VN),
                "Expected faces to have vertex-normals.");

        chk_msg(fv->a >= 1, "Expected fv.a >= 1");
        chk_msg(fv->b >= 1, "Expected fv.b >= 1");
        chk_msg(fv->c >= 1, "Expected fv.c >= 1");
        fv->a -= 1;
        fv->b -= 1;
        fv->c -= 1;

        chk_msg(fvn->a >= 1, "Expected fvn.a >= 1");
        chk_msg(fvn->b >= 1, "Expected fvn.b >= 1");
        chk_msg(fvn->c >= 1, "Expected fvn.c >= 1");
        fvn->a -= 1;
        fvn->b -= 1;
        fvn->c -= 1;
        return 1;
chk_error:
        return 0;
}

int mli_Object_malloc_from_wavefront(struct mli_Object *obj, struct mli_IO *io)
{
        uint64_t i = 0u;
        uint64_t line_number = 0u;
        uint64_t mtl = 0u;

        struct mli_String line = mli_String_init();
        struct mli_String tmp = mli_String_init();
        struct mli_String svn = mli_String_init();
        struct mli_String sv = mli_String_init();
        struct mli_String sf = mli_String_init();
        struct mli_String susemtl = mli_String_init();

        /* init dyn */
        struct mli_VecVector v = mli_VecVector_init();
        struct mli_VecVector vn = mli_VecVector_init();

        struct mli_object_FaceVector fv = mli_object_FaceVector_init();
        struct mli_object_FaceVector fvn = mli_object_FaceVector_init();
        struct mli_Uint32Vector fm = mli_Uint32Vector_init();

        struct mli_Map material_names = mli_Map_init();

        /* malloc dyn */
        chk(mli_VecVector_malloc(&v, 0u));
        chk(mli_VecVector_malloc(&vn, 0u));

        chk(mli_object_FaceVector_malloc(&fv, 0u));
        chk(mli_object_FaceVector_malloc(&fvn, 0u));
        chk(mli_Uint32Vector_malloc(&fm, 0u));

        chk(mli_Map_malloc(&material_names));

        chk(mli_String_from_cstr(&svn, "vn "));
        chk(mli_String_from_cstr(&sv, "v "));
        chk(mli_String_from_cstr(&sf, "f "));
        chk(mli_String_from_cstr(&susemtl, "usemtl "));

        chk_dbg;
        /* parse wavefront into dyn */
        while (1) {
                line_number += 1;
                chk_msg(line_number < 1000 * 1000 * 1000,
                        "Expected less than 1e9 lines in wavefront obj file. "
                        "Something went wrong.");

                if (mli_IO_eof(io)) {
                        break;
                }
                chk_msg(mli_IO_text_read_line(io, &line, '\n'),
                        "Can not read line.");
                chk_msg(line.size < 1024, "Expected line.size < 1024.");

                if (line.size > 0) {
                        if (mli_String_starts_with(&line, &svn)) {
                                struct mli_Vec tmp_vn;
                                chk(mli_String_copyn(
                                        &tmp, &line, 2, line.size - 2));
                                chk_msg(mli_Object_parse_three_float_line(
                                                &tmp, &tmp_vn),
                                        "Can not parse vertex-normal-line.");
                                chk_msg(mli_Vec_dot(tmp_vn, tmp_vn) > 0.0,
                                        "vn can not be normalized.") tmp_vn =
                                        mli_Vec_normalized(tmp_vn);

                                chk(mli_VecVector_push_back(&vn, tmp_vn));
                        } else if (mli_String_starts_with(&line, &sv)) {
                                struct mli_Vec tmp_v;
                                chk(mli_String_copyn(
                                        &tmp, &line, 1, line.size - 1));
                                chk_msg(mli_Object_parse_three_float_line(
                                                &tmp, &tmp_v),
                                        "Can not parse vertex-line.");
                                chk(mli_VecVector_push_back(&v, tmp_v));
                        } else if (mli_String_starts_with(&line, &sf)) {
                                struct mli_object_Face tmp_fv;
                                struct mli_object_Face tmp_fvn;
                                chk_msg(mli_Map_size(&material_names) > 0,
                                        "Expected 'usemtl' before first "
                                        "face 'f'.");
                                chk_msg(mli_Object_parse_face_vertices_and_normals(
                                                &line, &tmp_fv, &tmp_fvn),
                                        "Failed to parse face-line.");
                                chk(mli_object_FaceVector_push_back(
                                        &fv, tmp_fv));
                                chk(mli_object_FaceVector_push_back(
                                        &fvn, tmp_fvn));
                                chk(mli_Uint32Vector_push_back(&fm, mtl));
                        } else if (mli_String_starts_with(&line, &susemtl)) {
                                chk(mli_String_copyn(
                                        &tmp, &line, 7, line.size - 7));
                                if (!mli_Map_has(&material_names, &tmp)) {
                                        chk(mli_Map_insert(
                                                &material_names, &tmp, 0));
                                }
                                chk(mli_Map_find(&material_names, &tmp, &mtl));
                        } else {
                                /*fprintf(stderr, "no match: '%s'\n",
                                 * line.array);*/
                        }
                } /* line_length > 0 */
        }
        chk_dbg

                /* copy dyn into static mli_Object */
                chk_msg(fv.size == fvn.size,
                        "Expected num. vertex-indices == num. "
                        "vertex-normal-indices.");
        chk_msg(mli_Object_malloc(
                        obj,
                        v.size,
                        vn.size,
                        fv.size,
                        mli_Map_size(&material_names)),
                "Failed to malloc mli_Object from file.");
        chk_dbg MLI_MATH_NCPY(v.array, obj->vertices, v.size);
        MLI_MATH_NCPY(vn.array, obj->vertex_normals, vn.size);
        MLI_MATH_NCPY(fv.array, obj->faces_vertices, fv.size);
        MLI_MATH_NCPY(fvn.array, obj->faces_vertex_normals, fvn.size);
        MLI_MATH_NCPY(fm.array, obj->faces_materials, fm.size);
        for (i = 0; i < mli_Map_size(&material_names); i++) {
                chk(mli_String_copy(
                        &obj->material_names[i],
                        &material_names.items.array[i].key));
        }

        chk_msg(mli_Object_is_valid(obj), "Expected object to be valid.");

        /* free dyn */
        mli_VecVector_free(&v);
        mli_VecVector_free(&vn);

        mli_object_FaceVector_free(&fv);
        mli_object_FaceVector_free(&fvn);
        mli_Uint32Vector_free(&fm);

        mli_Map_free(&material_names);

        mli_String_free(&line);
        mli_String_free(&tmp);
        mli_String_free(&svn);
        mli_String_free(&sv);
        mli_String_free(&sf);
        mli_String_free(&susemtl);

        return 1;
chk_error:
        mli_Object_free(obj);

        /* free dyn */
        mli_VecVector_free(&v);
        mli_VecVector_free(&vn);

        mli_object_FaceVector_free(&fv);
        mli_object_FaceVector_free(&fvn);
        mli_Uint32Vector_free(&fm);

        mli_Map_free(&material_names);

        mli_String_free(&line);
        mli_String_free(&tmp);
        mli_String_free(&svn);
        mli_String_free(&sv);
        mli_String_free(&sf);
        mli_String_free(&susemtl);

        return 0;
}

int mli_Object_fprint_to_wavefront(
        struct mli_IO *f,
        const struct mli_Object *obj)
{
        uint32_t i, mtl, face;
        chk(mli_IO_text_write_cstr_format(f, "# vertices\n"));
        for (i = 0; i < obj->num_vertices; i++) {
                chk(mli_IO_text_write_cstr_format(
                        f,
                        "v %.6f %.6f %.6f\n",
                        obj->vertices[i].x,
                        obj->vertices[i].y,
                        obj->vertices[i].z));
        }

        chk(mli_IO_text_write_cstr_format(f, "# vertex normals\n"));
        for (i = 0; i < obj->num_vertex_normals; i++) {
                chk(mli_IO_text_write_cstr_format(
                        f,
                        "vn %.6f %.6f %.6f\n",
                        obj->vertex_normals[i].x,
                        obj->vertex_normals[i].y,
                        obj->vertex_normals[i].z));
        }

        chk(mli_IO_text_write_cstr_format(f, "# faces\n"));
        for (face = 0; face < obj->num_faces; face++) {
                if ((face == 0) || (mtl != obj->faces_materials[face])) {
                        mtl = obj->faces_materials[face];
                        chk(mli_IO_text_write_cstr_format(
                                f,
                                "usemtl %s\n",
                                obj->material_names[mtl].array));
                }

                chk(mli_IO_text_write_cstr_format(
                        f,
                        "f %d//%d %d//%d %d//%d\n",
                        obj->faces_vertices[face].a + 1,
                        obj->faces_vertex_normals[face].a + 1,
                        obj->faces_vertices[face].b + 1,
                        obj->faces_vertex_normals[face].b + 1,
                        obj->faces_vertices[face].c + 1,
                        obj->faces_vertex_normals[face].c + 1));
        }

        return 1;
chk_error:
        return 0;
}

/* octree */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/*
 * mli_OcTree
 * =========
 *
 * A cache-aware mli_OcTree to accelerate the traversal.
 * mli_OcTree is created from the dynamic, mli_octree_TmpOcTree.
 *
 * Structure of mli_octree_TmpOcTree
 * -------------------------
 *
 *                               mli_octree_TmpNode(0,-)
 *                                    |-objects[a,b,c,d,e,f,g,h]
 *                     _______________|________________
 *                    |                                |
 *                mli_octree_TmpNode(1,-) mli_octree_TmpNode(2,-)
 *                    |objects[a,b,c,d]                |objects[e,f,g,h]
 *             _______|_______                 ________|________
 *            |               |               |                 |
 *        mli_octree_TmpNode(3,0) mli_octree_TmpNode(4,1)
 * mli_octree_TmpNode(5,2)   mli_octree_TmpNode(6,3) |objects[a,b,c]
 * |objects[c,d]  |objects[e,f]     |objects[f,g,h]
 *
 *      mli_octree_TmpNode link to each other via pointer.
 *      Each mli_octree_TmpNode is allocated separately.
 *
 *      advantages:
 *              easy to grow and fill
 *      disadvantage:
 *              memory is scattered, difficult to serialize
 *
 * Structure of mli_OcTree
 * ----------------------
 *
 *                object_links w.r.t. mli_Geometry ->
 *          leafs +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
 *            |---|  a  |  b  |  c  |  c  |  d  |  e  |  f  |  f  |  g  |  h  |
 *            |   +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
 * mli_OcTree  |   0     1     2     3     4     5     6     7     8     9 10 |
 * |   ^                 ^           ^           ^ |         |
 *  |---------|   addresses
 *  |         |   +-----+-----+-----+-----+
 *  |         |---|  0  |  3  |  4  |  6  |
 *  |             +-----+-----+-----+-----+ first w.r.t object_links
 *  |             |  3  |  2  |  2  |  3  |
 *  |             +-----+-----+-----+-----+ size w.r.t object_links
 *  |             0     1     2     3
 *  |
 *  |         nodes ->
 *  |         +-----+-----+-----+
 *  |---------|node |leaf |leaf |
 *            +-----+-----+-----+ ^
 *            |node |leaf |leaf | |
 *            +-----+-----+-----+ types   type == node
 *            |  1  |  0  |  2  |                  w.r.t. mli_OcTree.nodes
 *            +-----+-----+-----+ ^       type == leaf
 *            |  2  |  1  |  3  | |               w.r.t. leafs.addresses
 *            +-----+-----+-----+ children
 *            0     1     2
 *
 *      advantages:
 *              small memory footprint, contoinous memory blocks,
 *              easy to serialze
 *      disadvantage:
 *              can not grow
 */

struct mli_octree_LeafAddress mli_octree_LeafAddress_init(void)
{
        struct mli_octree_LeafAddress address;
        address.first_object_link = 0u;
        address.num_object_links = 0u;
        return address;
}

struct mli_octree_LeafArray mli_octree_LeafArray_init(void)
{
        struct mli_octree_LeafArray leafs;
        leafs.num_leafs = 0;
        leafs.adresses = NULL;
        leafs.num_object_links = 0u;
        leafs.object_links = NULL;
        return leafs;
}

void mli_octree_LeafArray_free(struct mli_octree_LeafArray *leafs)
{
        free(leafs->object_links);
        free(leafs->adresses);
        *leafs = mli_octree_LeafArray_init();
}

int mli_octree_LeafArray_malloc(
        struct mli_octree_LeafArray *leafs,
        const uint64_t num_leafs,
        const uint64_t num_object_links)
{
        uint64_t i;
        mli_octree_LeafArray_free(leafs);
        leafs->num_leafs = num_leafs;
        chk_malloc(
                leafs->adresses,
                struct mli_octree_LeafAddress,
                leafs->num_leafs);
        for (i = 0; i < leafs->num_leafs; i++) {
                leafs->adresses[i] = mli_octree_LeafAddress_init();
        }
        leafs->num_object_links = num_object_links;
        chk_malloc(leafs->object_links, uint32_t, leafs->num_object_links);
        for (i = 0; i < leafs->num_object_links; i++) {
                leafs->object_links[i] = 0u;
        }
        return 1;
chk_error:
        return 0;
}

struct mli_octree_Node mli_octree_Node_init(void)
{
        uint64_t c = 0;
        struct mli_octree_Node node;
        for (c = 0; c < 8u; c++) {
                node.children[c] = 0u;
                node.types[c] = MLI_OCTREE_TYPE_NONE;
        }
        return node;
}

struct mli_OcTree mli_OcTree_init(void)
{
        struct mli_OcTree tree;
        tree.cube.lower = mli_Vec_init(0., 0., 0.);
        tree.cube.edge_length = 0.;
        tree.num_nodes = 0u;
        tree.nodes = NULL;
        tree.leafs = mli_octree_LeafArray_init();
        tree.root_type = MLI_OCTREE_TYPE_NONE;
        return tree;
}

void mli_OcTree_free(struct mli_OcTree *tree)
{
        free(tree->nodes);
        mli_octree_LeafArray_free(&tree->leafs);
        *tree = mli_OcTree_init();
}

int mli_OcTree_malloc(
        struct mli_OcTree *tree,
        const uint64_t num_nodes,
        const uint64_t num_leafs,
        const uint64_t num_object_links)
{
        uint64_t i;
        mli_OcTree_free(tree);
        tree->num_nodes = num_nodes;
        chk_malloc(tree->nodes, struct mli_octree_Node, tree->num_nodes);
        for (i = 0; i < tree->num_nodes; i++)
                tree->nodes[i] = mli_octree_Node_init();
        chk(mli_octree_LeafArray_malloc(
                &tree->leafs, num_leafs, num_object_links));
        return 1;
chk_error:
        return 0;
}

void mli_OcTree_set_node(
        struct mli_OcTree *tree,
        const struct mli_octree_TmpNode *dynnode)
{
        uint64_t c;
        uint64_t i = dynnode->node_index;
        assert(i < tree->num_nodes);
        for (c = 0; c < 8u; c++) {
                if (dynnode->children[c] != NULL) {
                        if (dynnode->children[c]->node_index >= 0) {
                                tree->nodes[i].children[c] =
                                        dynnode->children[c]->node_index;
                                tree->nodes[i].types[c] = MLI_OCTREE_TYPE_NODE;
                        } else if (dynnode->children[c]->leaf_index >= 0) {
                                tree->nodes[i].children[c] =
                                        dynnode->children[c]->leaf_index;
                                tree->nodes[i].types[c] = MLI_OCTREE_TYPE_LEAF;
                        } else {
                                tree->nodes[i].children[c] = 0;
                                tree->nodes[i].types[c] = MLI_OCTREE_TYPE_NONE;
                        }
                } else {
                        tree->nodes[i].children[c] = 0;
                        tree->nodes[i].types[c] = MLI_OCTREE_TYPE_NONE;
                }
        }
}

void mli_OcTree_set_leaf(
        struct mli_OcTree *tree,
        const struct mli_octree_TmpNode *dynnode,
        uint64_t *object_link_size)
{
        uint64_t o;
        uint64_t i = dynnode->leaf_index;
        assert(i < tree->leafs.num_leafs);
        tree->leafs.adresses[i].first_object_link = *object_link_size;
        tree->leafs.adresses[i].num_object_links = dynnode->num_objects;
        *object_link_size += dynnode->num_objects;
        for (o = 0; o < dynnode->num_objects; o++) {
                uint64_t l = tree->leafs.adresses[i].first_object_link + o;
                assert(l < tree->leafs.num_object_links);
                tree->leafs.object_links[l] = dynnode->objects[o];
        }
}

void mli_OcTree_set_walk(
        struct mli_OcTree *tree,
        const struct mli_octree_TmpNode *dynnode,
        uint64_t *object_link_size)
{
        uint64_t c;
        if (dynnode->node_index >= 0) {
                mli_OcTree_set_node(tree, dynnode);
        } else if (dynnode->leaf_index >= 0) {
                mli_OcTree_set_leaf(tree, dynnode, object_link_size);
        }

        for (c = 0; c < 8u; c++) {
                if (dynnode->children[c] != NULL) {
                        mli_OcTree_set_walk(
                                tree, dynnode->children[c], object_link_size);
                }
        }
}

void mli_OcTree_set(
        struct mli_OcTree *tree,
        const struct mli_octree_TmpOcTree *dyntree)
{
        uint64_t object_link_size = 0u;
        tree->cube = dyntree->cube;
        if (mli_octree_TmpNode_num_children(&dyntree->root) > 0) {
                tree->root_type = MLI_OCTREE_TYPE_NODE;
        } else {
                tree->root_type = MLI_OCTREE_TYPE_LEAF;
        }

        mli_OcTree_set_walk(tree, &dyntree->root, &object_link_size);
}

uint64_t mli_OcTree_node_num_children(
        const struct mli_OcTree *tree,
        const uint64_t node_idx)
{
        uint64_t num = 0u;
        uint64_t c;
        for (c = 0; c < 8u; c++) {
                if (tree->nodes[node_idx].types[c] > MLI_OCTREE_TYPE_NONE) {
                        num++;
                }
        }
        return num;
}

uint64_t mli_OcTree_leaf_num_objects(
        const struct mli_OcTree *tree,
        const uint64_t leaf)
{
        return tree->leafs.adresses[leaf].num_object_links;
}

uint32_t mli_OcTree_leaf_object_link(
        const struct mli_OcTree *tree,
        const uint64_t leaf,
        const uint64_t object_link)
{
        uint64_t i = tree->leafs.adresses[leaf].first_object_link + object_link;
        return tree->leafs.object_links[i];
}

int mli_OcTree_equal_payload_walk(
        const struct mli_OcTree *tree,
        const int32_t node_idx,
        const int32_t node_type,
        const struct mli_octree_TmpNode *tmp_node)
{
        if (node_type == MLI_OCTREE_TYPE_LEAF) {
                /* leaf */
                uint64_t leaf_idx;
                uint64_t obj;
                chk_msg(mli_octree_TmpNode_num_children(tmp_node) == 0,
                        "Expect tmp_node to have 0 children when "
                        "node_type == LEAF.");
                leaf_idx = node_idx;
                chk_msg(leaf_idx < tree->leafs.num_leafs,
                        "The leaf_idx is out of range.");
                chk_msg(tree->leafs.adresses[leaf_idx].num_object_links ==
                                tmp_node->num_objects,
                        "Expected leafs to have equal num_object_links.");
                for (obj = 0; obj < tmp_node->num_objects; obj++) {
                        uint64_t l = obj + tree->leafs.adresses[leaf_idx]
                                                   .first_object_link;
                        chk_msg(tree->leafs.object_links[l] ==
                                        tmp_node->objects[obj],
                                "Expected object_links in leaf to be equal.");
                }
                chk_msg(tree->leafs.adresses[leaf_idx].num_object_links ==
                                tmp_node->num_objects,
                        "Expected leafs to have equal num_object_links.");
        } else if (node_type == MLI_OCTREE_TYPE_NODE) {
                /* node */
                uint64_t c;
                chk_msg(node_idx >= 0,
                        "This node_idx is expected to point to a node.");
                for (c = 0; c < 8u; c++) {
                        if (tmp_node->children[c] == NULL) {
                                chk_msg(tree->nodes[node_idx].children[c] == 0,
                                        "Expected node's child == "
                                        "0 when tmp_node's child == NULL");
                        } else {
                                if (tmp_node->children[c]->num_objects == 0) {
                                        chk_msg(tree->nodes[node_idx]
                                                                .children[c] ==
                                                        0,
                                                "Expected node's child != "
                                                "0 when tmp_node's child != "
                                                "NULL");
                                }
                        }
                }

                for (c = 0; c < 8u; c++) {
                        if (tmp_node->children[c] != NULL) {
                                uint64_t child_node_idx =
                                        tree->nodes[node_idx].children[c];
                                int32_t child_node_type =
                                        tree->nodes[node_idx].types[c];
                                chk_msg(mli_OcTree_equal_payload_walk(
                                                tree,
                                                child_node_idx,
                                                child_node_type,
                                                tmp_node->children[c]),
                                        "Expected tree to be euqal further "
                                        "down");
                        }
                }
        } else if (node_type == MLI_OCTREE_TYPE_NONE) {

        } else {
                chk_bad("node_idx must be either node, leaf or none");
        }

        return 1;
chk_error:
        return 0;
}

int mli_OcTree_equal_payload(
        const struct mli_OcTree *tree,
        const struct mli_octree_TmpOcTree *tmp_octree)
{
        int32_t root_node_idx = 0;
        int32_t root_node_type = MLI_OCTREE_TYPE_NODE;
        chk_msg(mli_Cube_equal(tree->cube, tmp_octree->cube),
                "Cubes are not equal");
        chk_msg(mli_OcTree_equal_payload_walk(
                        tree, root_node_idx, root_node_type, &tmp_octree->root),
                "Tree is not equal");
        return 1;
chk_error:
        return 0;
}

void mli_OcTree_print_walk(
        const struct mli_OcTree *tree,
        const int32_t node_idx,
        const uint8_t node_type,
        const uint32_t indent,
        const uint32_t child)
{
        uint32_t i;
        uint32_t c;
        for (i = 0u; i < indent; i++)
                printf(" ");
        if (node_type == MLI_OCTREE_TYPE_NONE) {
                printf("|-Leaf[%d, %d] %u: %u []", -1, -1, child, 0);
        } else if (node_type == MLI_OCTREE_TYPE_LEAF) {
                int32_t leaf_idx = node_idx;
                uint32_t j;
                assert(leaf_idx < (int32_t)tree->leafs.num_leafs);
                printf("|-Leaf[%d, %d] %u: %u [",
                       -1,
                       leaf_idx,
                       child,
                       tree->leafs.adresses[leaf_idx].num_object_links);
                for (j = 0; j < tree->leafs.adresses[leaf_idx].num_object_links;
                     j++) {
                        int32_t l = j + tree->leafs.adresses[leaf_idx]
                                                .first_object_link;
                        printf("%u, ", tree->leafs.object_links[l]);
                }
                printf("]");
        } else if (node_type == MLI_OCTREE_TYPE_NODE) {
                printf("Node[%d, %d]: %u", node_idx, -1, child);
        }
        printf("\n");
        for (c = 0u; c < 8u; c++) {
                if (node_type == MLI_OCTREE_TYPE_NODE) {
                        int32_t child_node_idx;
                        int32_t child_node_type;
                        assert(node_idx < (int32_t)tree->num_nodes);
                        child_node_idx = tree->nodes[node_idx].children[c];
                        child_node_type = tree->nodes[node_idx].types[c];
                        if (child_node_type != MLI_OCTREE_TYPE_NONE) {
                                mli_OcTree_print_walk(
                                        tree,
                                        child_node_idx,
                                        child_node_type,
                                        indent + 2,
                                        c);
                        }
                }
        }
}

void mli_OcTree_print(const struct mli_OcTree *tree)
{
        printf("__ mli_OcTree __\n");
        printf("- num_nodes %u\n", (uint32_t)tree->num_nodes);
        printf("- num_leafs %u\n", (uint32_t)tree->leafs.num_leafs);
        printf("- root_type %d\n", tree->root_type);
        printf("- cube.lower %.1f, %.1f, %.1f\n",
               tree->cube.lower.x,
               tree->cube.lower.y,
               tree->cube.lower.z);
        printf("- cube.edge_length %.1f\n", tree->cube.edge_length);
        if (tree->num_nodes > 0) {
                int32_t root_idx = 0;
                mli_OcTree_print_walk(tree, root_idx, tree->root_type, 0u, 0u);
        }
}

int mli_OcTree_malloc_from_object_wavefront(
        struct mli_OcTree *octree,
        const struct mli_Object *object)
{
        uint64_t num_nodes = 0;
        uint64_t num_leafs = 0;
        uint64_t num_object_links = 0;
        struct mli_octree_TmpOcTree tmp_octree = mli_octree_TmpOcTree_init();
        chk_msg(mli_octree_TmpOcTree_malloc_from_bundle(
                        &tmp_octree,
                        (const void *)object,
                        object->num_faces,
                        mli_Object_face_in_local_frame_has_overlap_aabb_void,
                        mli_Object_aabb_in_local_frame(object)),
                "Failed to create dynamic, and temporary TmpOcTree "
                "from mli_Object");
        mli_octree_TmpNode_set_flat_index(&tmp_octree.root);
        mli_octree_TmpNode_num_nodes_leafs_objects(
                &tmp_octree.root, &num_nodes, &num_leafs, &num_object_links);

        chk_msg(mli_OcTree_malloc(
                        octree, num_nodes, num_leafs, num_object_links),
                "Failed to allocate cache-aware octree from dynamic octree.");
        mli_OcTree_set(octree, &tmp_octree);
        mli_octree_TmpOcTree_free(&tmp_octree);

        return 1;
chk_error:
        return 0;
}

int mli_OcTree_malloc_from_Geometry(
        struct mli_OcTree *octree,
        const struct mli_GeometryAndAccelerator *accgeo,
        const struct mli_AABB outermost_aabb)
{
        uint64_t num_nodes = 0;
        uint64_t num_leafs = 0;
        uint64_t num_object_links = 0;
        struct mli_octree_TmpOcTree tmp_octree = mli_octree_TmpOcTree_init();
        chk_msg(mli_octree_TmpOcTree_malloc_from_bundle(
                        &tmp_octree,
                        (const void *)accgeo,
                        accgeo->geometry->num_robjects,
                        mli_Geometry_robject_has_overlap_aabb_void,
                        outermost_aabb),
                "Failed to create dynamic, and temporary TmpOcTree "
                "from scenery(Geometry, Accelerator)");
        mli_octree_TmpNode_set_flat_index(&tmp_octree.root);
        mli_octree_TmpNode_num_nodes_leafs_objects(
                &tmp_octree.root, &num_nodes, &num_leafs, &num_object_links);

        chk_msg(mli_OcTree_malloc(
                        octree, num_nodes, num_leafs, num_object_links),
                "Failed to allocate cache-aware octree from dynamic octree.");
        mli_OcTree_set(octree, &tmp_octree);
        mli_octree_TmpOcTree_free(&tmp_octree);

        return 1;
chk_error:
        return 0;
}

/* octree_equal */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_OcTree_equal(const struct mli_OcTree *a, const struct mli_OcTree *b)
{
        uint64_t i, j;
        chk_msg(a->num_nodes == b->num_nodes, "num_nodes not equal.");
        chk_msg(a->leafs.num_leafs == b->leafs.num_leafs,
                "leafs.num_leafs not equal.");
        chk_msg(a->leafs.num_object_links == b->leafs.num_object_links,
                "num_object_links not equal.");

        chk_msg(a->cube.lower.x == b->cube.lower.x, "cube.lower.x not eaual");
        chk_msg(a->cube.lower.y == b->cube.lower.y, "cube.lower.y not eaual");
        chk_msg(a->cube.lower.z == b->cube.lower.z, "cube.lower.z not eaual");
        chk_msg(a->cube.edge_length == b->cube.edge_length,
                "cube.edge_length not eaual");

        chk_msg(a->root_type == b->root_type, "cube.root_type not eaual");

        for (i = 0; i < a->num_nodes; i++) {
                for (j = 0; j < 8; j++) {
                        chk_msg(a->nodes[i].children[j] ==
                                        b->nodes[i].children[j],
                                "octree.nodes[i].children[j] not equal");
                        chk_msg(a->nodes[i].types[j] == b->nodes[i].types[j],
                                "octree.nodes[i].types[j] not equal");
                }
        }

        for (i = 0; i < a->leafs.num_leafs; i++) {
                chk_msg(a->leafs.adresses[i].first_object_link ==
                                b->leafs.adresses[i].first_object_link,
                        "octree.leafs.adresses[i].first_object_link not equal");
                chk_msg(a->leafs.adresses[i].num_object_links ==
                                b->leafs.adresses[i].num_object_links,
                        "octree.leafs.adresses[i].num_object_links not equal");
        }
        for (i = 0; i < a->leafs.num_object_links; i++) {
                chk_msg(a->leafs.object_links[i] == b->leafs.object_links[i],
                        "octree.leafs.object_links[i] not equal");
        }

        return 1;
chk_error:
        return 0;
}

/* octree_overlaps */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/* octree_serialize */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_OcTree_to_io(const struct mli_OcTree *octree, struct mli_IO *f)
{
        struct mli_MagicId magic;

        /* magic identifier */
        chk(mli_MagicId_set(&magic, "mli_OcTree"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        /* capacity */
        chk_IO_write(&octree->num_nodes, sizeof(uint64_t), 1u, f);
        chk_IO_write(&octree->leafs.num_leafs, sizeof(uint64_t), 1u, f);
        chk_IO_write(&octree->leafs.num_object_links, sizeof(uint64_t), 1u, f);

        /* nodes */
        chk_IO_write(
                octree->nodes,
                sizeof(struct mli_octree_Node),
                octree->num_nodes,
                f);

        /* leaf addresses */
        chk_IO_write(
                octree->leafs.adresses,
                sizeof(struct mli_octree_LeafAddress),
                octree->leafs.num_leafs,
                f);

        /* leaf links */
        chk_IO_write(
                octree->leafs.object_links,
                sizeof(uint32_t),
                octree->leafs.num_object_links,
                f);

        /* mli_Cube */
        chk_IO_write(&octree->cube.lower.x, sizeof(double), 1u, f);
        chk_IO_write(&octree->cube.lower.y, sizeof(double), 1u, f);
        chk_IO_write(&octree->cube.lower.z, sizeof(double), 1u, f);
        chk_IO_write(&octree->cube.edge_length, sizeof(double), 1u, f);

        /* root type */
        chk_IO_write(&octree->root_type, sizeof(uint8_t), 1u, f);

        return 1;
chk_error:
        return 0;
}

int mli_OcTree_from_io(struct mli_OcTree *octree, struct mli_IO *f)
{
        uint64_t num_nodes;
        uint64_t num_leafs;
        uint64_t num_object_links;
        struct mli_MagicId magic;

        /* magic identifier */
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_OcTree"));
        mli_MagicId_warn_version(&magic);

        /* capacity */
        chk_IO_read(&num_nodes, sizeof(uint64_t), 1u, f);
        chk_IO_read(&num_leafs, sizeof(uint64_t), 1u, f);
        chk_IO_read(&num_object_links, sizeof(uint64_t), 1u, f);

        chk_msg(mli_OcTree_malloc(
                        octree, num_nodes, num_leafs, num_object_links),
                "Can not malloc octree from file.");

        chk_IO_read(
                octree->nodes,
                sizeof(struct mli_octree_Node),
                octree->num_nodes,
                f);
        chk_IO_read(
                octree->leafs.adresses,
                sizeof(struct mli_octree_LeafAddress),
                octree->leafs.num_leafs,
                f);
        chk_IO_read(
                octree->leafs.object_links,
                sizeof(uint32_t),
                octree->leafs.num_object_links,
                f);

        /* mli_Cube */
        chk_IO_read(&octree->cube.lower.x, sizeof(double), 1u, f);
        chk_IO_read(&octree->cube.lower.y, sizeof(double), 1u, f);
        chk_IO_read(&octree->cube.lower.z, sizeof(double), 1u, f);
        chk_IO_read(&octree->cube.edge_length, sizeof(double), 1u, f);

        /* root type */
        chk_IO_read(&octree->root_type, sizeof(uint8_t), 1u, f);

        return 1;
chk_error:
        return 0;
}

/* octree_tmp */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

uint64_t mli_octree_guess_depth_based_on_num_objects(const uint64_t num_objects)
{
        return 3u + (uint64_t)ceil(log((double)num_objects) / log(8.0));
}

/*
 * The dynamic node
 * ================
 */

struct mli_octree_TmpNode mli_octree_TmpNode_init(void)
{
        struct mli_octree_TmpNode n;
        uint64_t c;
        for (c = 0; c < 8u; c++) {
                n.children[c] = NULL;
        }
        n.num_objects = 0u;
        n.objects = NULL;
        n.flat_index = MLI_OCTREE_TMPNODE_FLAT_INDEX_NONE;
        n.node_index = MLI_OCTREE_TMPNODE_FLAT_INDEX_NONE;
        n.leaf_index = MLI_OCTREE_TMPNODE_FLAT_INDEX_NONE;
        return n;
}

void mli_octree_TmpNode_free(struct mli_octree_TmpNode *n)
{
        uint32_t c;
        for (c = 0; c < 8u; c++)
                if (n->children[c] != NULL)
                        mli_octree_TmpNode_free(n->children[c]);
        free(n->objects);
}

int mli_octree_TmpNode_malloc(
        struct mli_octree_TmpNode *n,
        const uint32_t num_objects)
{
        mli_octree_TmpNode_free(n);
        n->num_objects = num_objects;
        chk_malloc(n->objects, uint32_t, n->num_objects);
        return 1;
chk_error:
        return 0;
}

uint32_t mli_octree_TmpNode_signs_to_child(
        const uint32_t sx,
        const uint32_t sy,
        const uint32_t sz)
{
        return 4 * sx + 2 * sy + 1 * sz;
}

int mli_octree_TmpNode_add_children(
        struct mli_octree_TmpNode *node,
        const void *bundle,
        int (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mli_AABB),
        const struct mli_Cube cube,
        const uint64_t depth,
        const uint64_t max_depth)
{
        uint32_t c;
        uint32_t sx, sy, sz, obj;
        struct mli_Cube child_cubes[8];
        struct mli_octree_OverlapVector overlap[8];

        if (node->num_objects <= 32u) {
                return 1;
        }

        if (depth == max_depth) {
                return 1;
        }

        for (c = 0u; c < 8u; c++) {
                overlap[c] = mli_octree_OverlapVector_init();
                chk(mli_octree_OverlapVector_malloc(
                        &overlap[c], node->num_objects));
        }

        /* sense possible children */
        for (sx = 0u; sx < 2u; sx++) {
                for (sy = 0u; sy < 2u; sy++) {
                        for (sz = 0u; sz < 2u; sz++) {
                                const uint32_t child =
                                        mli_octree_TmpNode_signs_to_child(
                                                sx, sy, sz);
                                child_cubes[child] =
                                        mli_Cube_octree_child(cube, sx, sy, sz);
                                for (obj = 0u; obj < node->num_objects; obj++) {
                                        const uint32_t object_idx =
                                                node->objects[obj];
                                        if (item_in_bundle_has_overlap_aabb(
                                                    bundle,
                                                    object_idx,
                                                    mli_Cube_to_aabb(
                                                            child_cubes
                                                                    [child]))) {
                                                chk(mli_octree_OverlapVector_push_back(
                                                        &overlap[child],
                                                        object_idx));
                                        }
                                }
                        }
                }
        }

        for (c = 0; c < 8u; c++) {
                chk_malloc(node->children[c], struct mli_octree_TmpNode, 1u);
                (*node->children[c]) = mli_octree_TmpNode_init();
                chk(mli_octree_TmpNode_malloc(
                        node->children[c], overlap[c].size));
                MLI_MATH_NCPY(
                        overlap[c].array,
                        node->children[c]->objects,
                        overlap[c].size);
        }

        for (c = 0u; c < 8u; c++) {
                mli_octree_OverlapVector_free(&overlap[c]);
        }

        for (c = 0; c < 8u; c++) {
                mli_octree_TmpNode_add_children(
                        node->children[c],
                        bundle,
                        item_in_bundle_has_overlap_aabb,
                        child_cubes[c],
                        depth + 1u,
                        max_depth);
        }

        return 1;
chk_error:
        return 0;
}

int mli_octree_TmpNode_malloc_tree_from_bundle(
        struct mli_octree_TmpNode *root_node,
        const void *bundle,
        const uint32_t num_items_in_bundle,
        int (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mli_AABB),
        const struct mli_Cube bundle_cube)
{
        uint32_t idx, start_depth, max_depth;
        start_depth = 0u;
        max_depth = mli_octree_guess_depth_based_on_num_objects(
                num_items_in_bundle);

        chk_msg(mli_octree_TmpNode_malloc(root_node, num_items_in_bundle),
                "Failed to allocate root-node in dynamic octree.");

        for (idx = 0; idx < root_node->num_objects; idx++) {
                root_node->objects[idx] = idx;
        }
        mli_octree_TmpNode_add_children(
                root_node,
                bundle,
                item_in_bundle_has_overlap_aabb,
                bundle_cube,
                start_depth,
                max_depth);
        return 1;
chk_error:
        return 0;
}

int mli_octree_TmpNode_num_children(const struct mli_octree_TmpNode *node)
{
        uint32_t c, num = 0;
        for (c = 0u; c < 8u; c++)
                if (node->children[c] != NULL)
                        num++;
        return num;
}

void mli_octree_TmpNode_print(
        const struct mli_octree_TmpNode *node,
        const uint32_t indent,
        const uint32_t child)
{
        uint32_t i;
        uint32_t c;
        uint32_t num_children = mli_octree_TmpNode_num_children(node);
        for (i = 0u; i < indent; i++)
                printf(" ");
        if (num_children == 0) {
                uint32_t j;
                printf("|-Leaf[%d, %d] %u: %u [",
                       node->node_index,
                       node->leaf_index,
                       child,
                       node->num_objects);
                for (j = 0; j < node->num_objects; j++) {
                        printf("%u, ", node->objects[j]);
                }
                printf("]");
        } else {
                printf("Node[%d, %d]: %u",
                       node->node_index,
                       node->leaf_index,
                       child);
        }
        printf("\n");
        for (c = 0u; c < 8u; c++) {
                if (node->children[c] != NULL) {
                        mli_octree_TmpNode_print(
                                node->children[c], indent + 2, c);
                }
        }
}

/*
 * Assign a flat_index to every node that carries objects.
 * The ordering is:
 *      in tree:
 *                          3
 *                      1--{
 *                     /    4
 *                 0--{
 *                    \    5
 *                     2--{
 *                         6
 *
 *      in flat list:
 *      0, 1, 2, 3, 4, 5, 6, ...
 */

int mli_octree_TmpNode_exists_and_has_objects(
        const struct mli_octree_TmpNode *node)
{
        if (node != NULL) {
                if (node->num_objects > 0u) {
                        return 1;
                }
        }
        return 0;
}

void mli_octree_TmpNode_set_flat_index_walk(
        struct mli_octree_TmpNode *node,
        int32_t *flat_index,
        int32_t *node_index,
        int32_t *leaf_index)
{
        uint64_t c;
        for (c = 0u; c < 8u; c++) {
                if (mli_octree_TmpNode_exists_and_has_objects(
                            node->children[c])) {
                        (*flat_index)++;
                        node->children[c]->flat_index = *flat_index;

                        if (mli_octree_TmpNode_num_children(
                                    node->children[c]) == 0) {
                                node->children[c]->leaf_index = *leaf_index;
                                (*leaf_index)++;
                        } else {
                                (*node_index)++;
                                node->children[c]->node_index = *node_index;
                        }
                }
        }
        for (c = 0u; c < 8u; c++) {
                if (mli_octree_TmpNode_exists_and_has_objects(
                            node->children[c])) {
                        mli_octree_TmpNode_set_flat_index_walk(
                                node->children[c],
                                flat_index,
                                node_index,
                                leaf_index);
                }
        }
}

void mli_octree_TmpNode_set_flat_index(struct mli_octree_TmpNode *root_node)
{
        int32_t flat_index = 0;
        int32_t node_index = 0;
        int32_t leaf_index = 0;
        root_node->flat_index = flat_index;

        if (mli_octree_TmpNode_num_children(root_node) == 0) {
                root_node->leaf_index = leaf_index;
        } else {
                root_node->node_index = node_index;
        }

        mli_octree_TmpNode_set_flat_index_walk(
                root_node, &flat_index, &node_index, &leaf_index);
}

/*
 * Find the number of valid nodes in dynamic tree
 */

void mli_octree_TmpNode_num_nodes_leafs_objects_walk(
        const struct mli_octree_TmpNode *node,
        uint64_t *num_nodes,
        uint64_t *num_leafs,
        uint64_t *num_object_links)
{
        uint64_t c;
        if (node->node_index != MLI_OCTREE_TMPNODE_FLAT_INDEX_NONE) {
                (*num_nodes)++;
        }
        if (node->leaf_index != MLI_OCTREE_TMPNODE_FLAT_INDEX_NONE) {
                (*num_leafs)++;
                (*num_object_links) += node->num_objects;
        }
        for (c = 0; c < 8u; c++) {
                if (node->children[c] != NULL) {
                        mli_octree_TmpNode_num_nodes_leafs_objects_walk(
                                node->children[c],
                                num_nodes,
                                num_leafs,
                                num_object_links);
                }
        }
}

void mli_octree_TmpNode_num_nodes_leafs_objects(
        const struct mli_octree_TmpNode *root_node,
        uint64_t *num_nodes,
        uint64_t *num_leafs,
        uint64_t *num_object_links)
{
        *num_nodes = 0;
        *num_leafs = 0;
        *num_object_links = 0;
        mli_octree_TmpNode_num_nodes_leafs_objects_walk(
                root_node, num_nodes, num_leafs, num_object_links);
}

/*
 * The dynamic octree
 * ==================
 */

struct mli_octree_TmpOcTree mli_octree_TmpOcTree_init(void)
{
        struct mli_octree_TmpOcTree octree;
        octree.cube.lower = mli_Vec_init(0., 0., 0.);
        octree.cube.edge_length = 0.;
        octree.root = mli_octree_TmpNode_init();
        return octree;
}

void mli_octree_TmpOcTree_free(struct mli_octree_TmpOcTree *octree)
{
        mli_octree_TmpNode_free(&octree->root);
}

int mli_octree_TmpOcTree_malloc_from_bundle(
        struct mli_octree_TmpOcTree *octree,
        const void *bundle,
        const uint32_t num_items_in_bundle,
        int (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mli_AABB),
        struct mli_AABB bundle_aabb)
{
        mli_octree_TmpOcTree_free(octree);
        octree->cube = mli_Cube_outermost_cube(bundle_aabb);
        chk_msg(mli_octree_TmpNode_malloc_tree_from_bundle(
                        &octree->root,
                        bundle,
                        num_items_in_bundle,
                        item_in_bundle_has_overlap_aabb,
                        octree->cube),
                "Failed to allocate dynamic octree from bundle.");
        return 1;
chk_error:
        return 0;
}

void mli_octree_TmpOcTree_print(const struct mli_octree_TmpOcTree *octree)
{
        mli_octree_TmpNode_print(&octree->root, 0u, 0u);
}

/* octree_valid */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_OcTree_valid(const struct mli_OcTree *octree)
{
        uint32_t n, c;
        chk_msg(!MLI_MATH_IS_NAN(octree->cube.lower.x),
                "cube.lower.x is 'nan'.");
        chk_msg(!MLI_MATH_IS_NAN(octree->cube.lower.y),
                "cube.lower.y is 'nan'.");
        chk_msg(!MLI_MATH_IS_NAN(octree->cube.lower.z),
                "cube.lower.z is 'nan'.");
        chk_msg(!MLI_MATH_IS_NAN(octree->cube.edge_length),
                "cube.edge_length is 'nan'.");
        chk_msg(octree->cube.edge_length >= 0.0,
                "Expected cube.edge_length >= 0.0.");

        for (n = 0u; n < octree->num_nodes; n++) {
                for (c = 0u; c < 8u; c++) {
                        if (octree->nodes[n].types[c] == MLI_OCTREE_TYPE_NONE) {
                                chk_msg(octree->nodes[n].children[c] == 0u,
                                        "Expected the address of a 'NONE' "
                                        "child to be '0'.");
                        } else if (
                                octree->nodes[n].types[c] ==
                                MLI_OCTREE_TYPE_NODE) {
                                chk_msg(octree->nodes[n].children[c] <
                                                octree->num_nodes,
                                        "Expected the address of a 'NODE' "
                                        "child to be < num_nodes.");
                        } else if (
                                octree->nodes[n].types[c] ==
                                MLI_OCTREE_TYPE_LEAF) {
                                chk_msg(octree->nodes[n].children[c] <
                                                octree->leafs.num_leafs,
                                        "Expected the address of a 'LEAF' "
                                        "child to be < leafs.num_leafs.");
                        } else {
                                chk_bad("Expected octree->nodes[n].type[c] "
                                        "to be either NONE, NODE, or LEAF.");
                        }
                }
        }
        return 1;
chk_error:
        return 0;
}

int mli_OcTree_valid_wrt_links(
        const struct mli_OcTree *octree,
        const uint32_t num_links)
{
        uint32_t n;
        for (n = 0u; n < octree->leafs.num_object_links; n++) {
                chk_msg(octree->leafs.object_links[n] < num_links,
                        "Expected object_links[n] <  num_links.");
        }
        return 1;
chk_error:
        return 0;
}

/* path */
/* ---- */


int mli_path_strip_this_dir(
        const struct mli_String *src,
        struct mli_String *dst)
{
        uint64_t start = 0;
        uint64_t length = 0;
        struct mli_String cpysrc = mli_String_init();
        chk_msg(src->array, "Expected src-string to be allocated.");
        chk_msg(mli_String_copy(&cpysrc, src), "Can not copy input.");
        mli_String_free(dst);

        while (start + 1 < cpysrc.size) {
                if (cpysrc.array[start] == '.' &&
                    cpysrc.array[start + 1] == '/') {
                        start += 2;
                } else {
                        break;
                }
        }
        length = cpysrc.size - start;
        chk(mli_String_copyn(dst, &cpysrc, start, length));
        mli_String_free(&cpysrc);
        return 1;
chk_error:
        mli_String_free(&cpysrc);
        mli_String_free(dst);
        return 0;
}

int mli_path_basename(const struct mli_String *src, struct mli_String *dst)
{
        int64_t pos_last_del = -1;
        mli_String_free(dst);
        chk_msg(src->array != NULL, "Expected src-path to be allocated");

        pos_last_del = mli_String_rfind(src, '/');

        if (pos_last_del < 0) {
                chk(mli_String_from_cstr_fromat(dst, src->array));
        } else {
                chk(mli_String_from_cstr_fromat(
                        dst, &src->array[pos_last_del + 1]));
        }
        return 1;
chk_error:
        mli_String_free(dst);
        return 0;
}

int mli_path_splitext(
        const struct mli_String *src,
        struct mli_String *dst,
        struct mli_String *ext)
{
        int64_t p = -1;
        int64_t d = -1;
        struct mli_String tmp = mli_String_init();
        chk_msg(src->array != NULL, "Expected src-path to be allocated");
        chk(mli_String_copy(&tmp, src));

        mli_String_free(dst);
        mli_String_free(ext);

        p = mli_String_rfind(&tmp, '.');
        d = mli_String_rfind(&tmp, '/');

        if (p <= 0 || d > p || ((d + 1 == p) && (p + 1 < (int64_t)tmp.size))) {
                chk(mli_String_copy(dst, &tmp));
                chk(mli_String_from_cstr_fromat(ext, ""));
        } else {
                chk(mli_String_copyn(dst, &tmp, 0, p));
                chk(mli_String_copyn(ext, &tmp, p + 1, tmp.size - p - 1));
        }

        mli_String_free(&tmp);
        return 1;
chk_error:
        mli_String_free(&tmp);
        mli_String_free(dst);
        mli_String_free(ext);
        return 0;
}

/* photon */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/* photon_interaction */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_photon_interaction_type_to_string(const int32_t type, char *s)
{
        switch (type) {
        case MLI_PHOTON_CREATION:
                sprintf(s, "creation");
                break;
        case MLI_PHOTON_ABSORBTION:
                sprintf(s, "absorbtion");
                break;
        case MLI_PHOTON_ABSORBTION_MEDIUM:
                sprintf(s, "absorbtion in medium");
                break;
        case MLI_PHOTON_FRESNEL_REFLECTION:
                sprintf(s, "Fresnel reflection");
                break;
        case MLI_PHOTON_REFRACTION:
                sprintf(s, "refraction");
                break;
        case MLI_PHOTON_SPECULAR_REFLECTION:
                sprintf(s, "specular reflection");
                break;
        case MLI_PHOTON_DIFFUSE_REFLECTION:
                sprintf(s, "diffuse reflection");
                break;
        default:
                chk_bad("PhotonInteraction.type is unknown.");
                break;
        }
        return 1;
chk_error:
        return 0;
}

int mli_photon_time_of_flight(
        const struct mli_Materials *materials,
        const struct mli_PhotonInteraction *phisec,
        const double wavelength,
        double *time_of_flight)
{
        double refractive_index;
        const uint64_t spc_idx =
                materials->media.array[phisec->medium_coming_from]
                        .refraction_spectrum;
        const struct mli_Func *spc =
                &materials->spectra.array[spc_idx].spectrum;
        chk_msg(mli_Func_evaluate(spc, wavelength, &refractive_index),
                "Failed to eval. refraction for wavelength.");

        (*time_of_flight) = (refractive_index * phisec->distance_of_ray) /
                            MLI_PHYSICS_SPEED_OF_LIGHT_M_PER_S;
        return 1;
chk_error:
        return 0;
}

/* photon_interaction_vector */
/* ------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

MLI_VECTOR_IMPLEMENTATION(
        mli_PhotonInteractionVector,
        struct mli_PhotonInteraction)

int mli_PhotonInteractionVector_time_of_flight(
        const struct mli_PhotonInteractionVector *history,
        const struct mli_Scenery *scenery,
        const double wavelength,
        double *total_time_of_flight)
{
        uint64_t i;
        (*total_time_of_flight) = 0.0;
        for (i = 0; i < history->size; i++) {
                double time_of_flight = 0.0;
                chk_msg(mli_photon_time_of_flight(
                                &scenery->materials,
                                &history->array[i],
                                wavelength,
                                &time_of_flight),
                        "Can't estimate time-of-flight.");
                (*total_time_of_flight) += time_of_flight;
        }
        return 1;
chk_error:
        return 0;
}

void mli_PhotonInteractionVector_print(
        const struct mli_PhotonInteractionVector *history,
        const struct mli_Scenery *scenery)
{
        uint64_t i;
        char type_string[1024];
        char out_in[] = "out->in";
        char in_out[] = "in->out";
        printf("History(%usize, %ucapacity)\n",
               (uint32_t)history->size,
               (uint32_t)history->capacity);
        printf("==================\n");

        printf("   #  "
               "(   id; robj,  obj, face)  "
               "[     x/m,     y/m,     z/m]  "
               "type         "
               "refr.    "
               "abs.    "
               "dist/m"
               "\n");
        printf("-----------------------------------------------------"
               "-----------------------------------------------------\n");

        for (i = 0; i < history->size; i++) {
                struct mli_PhotonInteraction phisec = history->array[i];
                printf(" % 3d  ", (int32_t)i);

                if (phisec.on_geometry_surface == 1) {
                        printf("(% 5d;% 5d,% 5d,% 5d)  ",
                               scenery->geometry
                                       .robject_ids[phisec.geometry_id.robj],
                               phisec.geometry_id.robj,
                               scenery->geometry
                                       .robjects[phisec.geometry_id.robj],
                               phisec.geometry_id.face);
                } else {
                        printf("            n/a            ");
                }

                printf("[% -.1e,% -.1e,% -.1e]  ",
                       phisec.position.x,
                       phisec.position.y,
                       phisec.position.z);

                mli_photon_interaction_type_to_string(phisec.type, type_string);

                printf("%-24s ", type_string);

                printf("{%-12s,%-12s}  ",
                       scenery->materials.media.array[phisec.medium_coming_from]
                               .name.array,
                       scenery->materials.media.array[phisec.medium_going_to]
                               .name.array);

                if (phisec.type == MLI_PHOTON_CREATION) {
                        printf(" n/a  ");
                } else {
                        printf(" %1.1e  ", phisec.distance_of_ray);
                }

                if (phisec.on_geometry_surface == 1) {
                        if (phisec.from_outside_to_inside) {
                                printf("%s", out_in);
                        } else {
                                printf("%s", in_out);
                        }
                } else {
                        printf("n/a");
                }
                printf("\n");
        }
}

/* photon_propagation */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_PhotonInteraction mliPhotonInteraction_from_Intersection(
        const int64_t type,
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec)
{
        struct mli_PhotonInteraction phia;

        struct mli_BoundaryLayer_Side side_coming_from, side_going_to;

        side_coming_from = mli_raytracing_get_side_coming_from(scenery, isec);
        side_going_to = mli_raytracing_get_side_going_to(scenery, isec);

        phia.type = type;
        phia.position = isec->position;
        phia.position_local = isec->position_local;
        phia.medium_coming_from = side_coming_from.medium;
        phia.medium_going_to = side_going_to.medium;
        phia.distance_of_ray = isec->distance_of_ray;
        phia.on_geometry_surface = 1;
        phia.geometry_id = isec->geometry_id;
        phia.from_outside_to_inside = isec->from_outside_to_inside;
        return phia;
}

int mli_propagate_photon_cooktorrance(
        struct mli_PhotonPropagation *env,
        const struct mli_IntersectionSurfaceNormal *isec)
{
        double spectral_reflection_propability;
        double diffuse;
        double specular;

        double rnd;
        struct mli_BoundaryLayer_Side side_coming_from =
                mli_raytracing_get_side_coming_from(env->scenery, isec);

        const struct mli_Surface_CookTorrance *cook =
                &env->scenery->materials.surfaces
                         .array[side_coming_from.surface]
                         .data.cooktorrance;
        const struct mli_Func *reflection_spectrum =
                &env->scenery->materials.spectra
                         .array[cook->reflection_spectrum]
                         .spectrum;

        chk_msg(mli_Func_evaluate(
                        reflection_spectrum,
                        env->photon->wavelength,
                        &spectral_reflection_propability),
                "Failed to eval. spectral_reflection_propability for "
                "wavelength.");

        diffuse = cook->diffuse_weight * spectral_reflection_propability;
        specular = cook->specular_weight * spectral_reflection_propability;

        rnd = mli_Prng_uniform(env->prng);
        /*
                                                      absorbtion
                  diffuse      specular        (1.0 - diffuse - specular)
                  __/\____  _____/\__________  ___________/\____________
                 /        \/                 \/                         \

                |.........|..................|..........................|
                |         |                  |                          |
               0.0      diffuse      diffuse + specular                1.0
        */
        if (rnd < diffuse) {
                chk(mli_PhotonInteractionVector_push_back(
                        env->history,
                        mliPhotonInteraction_from_Intersection(
                                MLI_PHOTON_DIFFUSE_REFLECTION,
                                env->scenery,
                                isec)));
                env->photon->ray = mli_Ray_set(
                        isec->position,
                        mli_lambertian_cosine_law_draw_direction_wrt_surface_normal(
                                env->prng, isec->surface_normal));
                chk_msg(mli_propagate_photon_env(env),
                        "Failed to continue after diffuse reflection "
                        "cooktorrance.");
        } else if (rnd < (specular + diffuse)) {
                chk(mli_PhotonInteractionVector_push_back(
                        env->history,
                        mliPhotonInteraction_from_Intersection(
                                MLI_PHOTON_SPECULAR_REFLECTION,
                                env->scenery,
                                isec)));
                env->photon->ray = mli_Ray_set(
                        isec->position,
                        mli_Vec_mirror(
                                env->photon->ray.direction,
                                isec->surface_normal));
                chk_msg(mli_propagate_photon_env(env),
                        "Failed to continue after specular reflection "
                        "cooktorrance.");
        } else {
                chk(mli_PhotonInteractionVector_push_back(
                        env->history,
                        mliPhotonInteraction_from_Intersection(
                                MLI_PHOTON_ABSORBTION, env->scenery, isec)));
        }
        return 1;
chk_error:
        return 0;
}

int mli_propagate_photon_pass_boundary_layer(
        struct mli_PhotonPropagation *env,
        const struct mli_IntersectionSurfaceNormal *isec,
        const struct mli_Fresnel fresnel)
{
        chk(mli_PhotonInteractionVector_push_back(
                env->history,
                mliPhotonInteraction_from_Intersection(
                        MLI_PHOTON_REFRACTION, env->scenery, isec)));
        env->photon->ray = mli_Ray_set(
                isec->position, mli_Fresnel_refraction_direction(fresnel));
        chk_msg(mli_propagate_photon_env(env),
                "Failed to continue after passing boundary layer");
        return 1;
chk_error:
        return 0;
}

int mli_propagate_photon_probability_passing_medium_coming_from(
        const struct mli_Scenery *scenery,
        const struct mli_Photon *photon,
        const struct mli_IntersectionSurfaceNormal *isec,
        double *probability_passing)
{
        const struct mli_BoundaryLayer_Side side_coming_from =
                mli_raytracing_get_side_coming_from(scenery, isec);
        const struct mli_Medium *medium_coming_from =
                &scenery->materials.media.array[side_coming_from.medium];

        const struct mli_Func *absorbtion_spectrum =
                &scenery->materials.spectra
                         .array[medium_coming_from->absorbtion_spectrum]
                         .spectrum;
        double one_over_e_way;

        chk_msg(mli_Func_evaluate(
                        absorbtion_spectrum,
                        photon->wavelength,
                        &one_over_e_way),
                "Photon's wavelength is out of range to "
                "evaluate absorbtion in medium coming from");
        (*probability_passing) = exp(-isec->distance_of_ray / one_over_e_way);

        return 1;
chk_error:
        return 0;
}

int mli_propagate_photon_fresnel_refraction_and_reflection(
        struct mli_PhotonPropagation *env,
        const struct mli_IntersectionSurfaceNormal *isec)
{
        struct mli_Fresnel fresnel;
        double n_going_to;
        double n_coming_from;
        double reflection_propability;
        struct mli_Vec facing_surface_normal;
        chk_msg(mli_Func_evaluate(
                        mli_raytracing_get_refractive_index_going_to(
                                env->scenery, isec),
                        env->photon->wavelength,
                        &n_going_to),
                "Failed to eval. refraction going to for wavelength.");
        chk_msg(mli_Func_evaluate(
                        mli_raytracing_get_refractive_index_coming_from(
                                env->scenery, isec),
                        env->photon->wavelength,
                        &n_coming_from),
                "Failed to eval. refraction coming from for wavelength.");
        facing_surface_normal =
                isec->from_outside_to_inside
                        ? isec->surface_normal
                        : mli_Vec_multiply(isec->surface_normal, -1.0);
        fresnel = mli_Fresnel_init(
                env->photon->ray.direction,
                facing_surface_normal,
                n_coming_from,
                n_going_to);
        reflection_propability = mli_Fresnel_reflection_propability(fresnel);
        if (reflection_propability > mli_Prng_uniform(env->prng)) {
                chk(mli_PhotonInteractionVector_push_back(
                        env->history,
                        mliPhotonInteraction_from_Intersection(
                                MLI_PHOTON_FRESNEL_REFLECTION,
                                env->scenery,
                                isec)));
                env->photon->ray = mli_Ray_set(
                        isec->position,
                        mli_Fresnel_reflection_direction(fresnel));
                chk_msg(mli_propagate_photon_env(env),
                        "Failed to continue after reflection");
        } else {
                chk_msg(mli_propagate_photon_pass_boundary_layer(
                                env, isec, fresnel),
                        "Failed to pass boundary");
        }
        return 1;
chk_error:
        return 0;
}

int mli_propagate_photon_interact_with_object(
        struct mli_PhotonPropagation *env,
        const struct mli_IntersectionSurfaceNormal *isec)
{
        const struct mli_BoundaryLayer_Side side_coming_from =
                mli_raytracing_get_side_coming_from(env->scenery, isec);
        const struct mli_Surface *surface_coming_from =
                &env->scenery->materials.surfaces
                         .array[side_coming_from.surface];

        switch (surface_coming_from->type) {
        case MLI_SURFACE_TYPE_TRANSPARENT:
                chk_msg(mli_propagate_photon_fresnel_refraction_and_reflection(
                                env, isec),
                        "Failed Fresnel.");
                break;
        case MLI_SURFACE_TYPE_COOKTORRANCE:
                chk_msg(mli_propagate_photon_cooktorrance(env, isec),
                        "Failed cook-torrance.");
                break;
        default:
                chk_bad("Unkown type of surface.");
                break;
        }
        return 1;
chk_error:
        return 0;
}

int mli_propagate_photon_distance_until_absorbtion(
        const struct mli_Func *absorbtion_in_medium_passing_through,
        const double wavelength,
        struct mli_Prng *prng,
        double *distance_until_absorbtion)
{
        double one_over_e_way;
        chk_msg(mli_Func_evaluate(
                        absorbtion_in_medium_passing_through,
                        wavelength,
                        &one_over_e_way),
                "Failed to eval. absorbtion for wavelength.");
        (*distance_until_absorbtion) =
                mli_Prng_expovariate(prng, 1. / one_over_e_way);
        return 1;
chk_error:
        return 0;
}

int mli_propagate_photon_work_on_causal_intersection(
        struct mli_PhotonPropagation *env)
{
        int ray_does_intersect_surface = 0;
        double distance_until_absorbtion = 0.0;
        struct mli_IntersectionSurfaceNormal next_intersection;
        struct mli_Func *absorbtion_in_medium_passing_through;
        struct mli_Medium *medium_passing_through;
        struct mli_PhotonInteraction phia;

        ray_does_intersect_surface =
                mli_raytracing_query_intersection_with_surface_normal(
                        env->scenery, env->photon->ray, &next_intersection);

        if (ray_does_intersect_surface) {
                int photon_is_absorbed_before_reaching_surface;
                struct mli_IntersectionLayer layer;

                layer = mli_raytracing_get_intersection_layer(
                        env->scenery, &next_intersection);

                absorbtion_in_medium_passing_through =
                        &env->scenery->materials.spectra
                                 .array[layer.side_coming_from.medium
                                                ->absorbtion_spectrum]
                                 .spectrum;

                chk(mli_propagate_photon_distance_until_absorbtion(
                        absorbtion_in_medium_passing_through,
                        env->photon->wavelength,
                        env->prng,
                        &distance_until_absorbtion));

                photon_is_absorbed_before_reaching_surface =
                        distance_until_absorbtion <
                        next_intersection.distance_of_ray;

                if (env->history->size == 0) {
                        /* creation */
                        phia.type = MLI_PHOTON_CREATION;
                        phia.position = env->photon->ray.support;
                        phia.position_local = phia.position;
                        phia.distance_of_ray = 0.0;
                        phia.on_geometry_surface = 0;
                        phia.geometry_id = mli_GeometryId_init();
                        phia.from_outside_to_inside = 1;

                        phia.medium_coming_from =
                                layer.side_coming_from.medium_idx;
                        phia.medium_going_to =
                                layer.side_coming_from.medium_idx;

                        chk(mli_PhotonInteractionVector_push_back(
                                env->history, phia));
                }

                if (photon_is_absorbed_before_reaching_surface) {
                        /* absorbtion in medium */
                        phia.type = MLI_PHOTON_ABSORBTION_MEDIUM;
                        phia.position = mli_Ray_at(
                                &env->photon->ray, distance_until_absorbtion);
                        ;
                        phia.position_local = phia.position;
                        phia.distance_of_ray = distance_until_absorbtion;
                        phia.on_geometry_surface = 0;
                        phia.geometry_id = mli_GeometryId_init();
                        phia.from_outside_to_inside = 1;

                        phia.medium_coming_from =
                                layer.side_coming_from.medium_idx;
                        phia.medium_going_to =
                                layer.side_coming_from.medium_idx;

                        chk(mli_PhotonInteractionVector_push_back(
                                env->history, phia));
                } else {
                        chk_msg(mli_propagate_photon_interact_with_object(
                                        env, &next_intersection),
                                "Failed to interact photon with object "
                                "surface.");
                }
        } else {
                const uint64_t default_medium =
                        env->scenery->materials.default_medium;

                medium_passing_through =
                        &env->scenery->materials.media.array[default_medium];
                absorbtion_in_medium_passing_through =
                        &env->scenery->materials.spectra
                                 .array[medium_passing_through
                                                ->absorbtion_spectrum]
                                 .spectrum;

                chk(mli_propagate_photon_distance_until_absorbtion(
                        absorbtion_in_medium_passing_through,
                        env->photon->wavelength,
                        env->prng,
                        &distance_until_absorbtion));

                if (env->history->size == 0) {
                        /* creation */
                        phia.type = MLI_PHOTON_CREATION;
                        phia.position = env->photon->ray.support;
                        phia.position_local = phia.position;
                        phia.distance_of_ray = 0.0;
                        phia.on_geometry_surface = 0;
                        phia.geometry_id = mli_GeometryId_init();
                        phia.from_outside_to_inside = 1;

                        phia.medium_coming_from = default_medium;
                        phia.medium_going_to = default_medium;

                        chk(mli_PhotonInteractionVector_push_back(
                                env->history, phia));
                }

                /* absorbtion in medium */
                phia.type = MLI_PHOTON_ABSORBTION_MEDIUM;
                phia.position = mli_Ray_at(
                        &env->photon->ray, distance_until_absorbtion);
                phia.position_local = phia.position;
                phia.distance_of_ray = distance_until_absorbtion;
                phia.on_geometry_surface = 0;
                phia.geometry_id = mli_GeometryId_init();
                phia.from_outside_to_inside = 1;

                phia.medium_coming_from = default_medium;
                phia.medium_going_to = default_medium;

                chk(mli_PhotonInteractionVector_push_back(env->history, phia));
        }

        return 1;
chk_error:
        return 0;
}

int mli_propagate_photon_env(struct mli_PhotonPropagation *env)
{
        if (env->max_interactions > env->history->size) {
                chk_msg(mli_propagate_photon_work_on_causal_intersection(env),
                        "Failed to work on intersection.");
        }
        return 1;
chk_error:
        return 0;
}

int mli_propagate_photon(
        const struct mli_Scenery *scenery,
        struct mli_PhotonInteractionVector *history,
        struct mli_Photon *photon,
        struct mli_Prng *prng,
        const uint64_t max_interactions)
{
        struct mli_PhotonPropagation env;
        env.scenery = scenery;
        env.history = history;
        env.photon = photon;
        env.prng = prng;
        env.max_interactions = max_interactions;
        chk(mli_propagate_photon_env(&env));
        return 1;
chk_error:
        return 0;
}

/* photon_source */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_photon_source_parallel_towards_z_from_xy_disc(
        struct mli_PhotonVector *out_photons,
        const double wavelength,
        const double radius,
        const uint64_t num_photons,
        struct mli_Prng *prng)
{
        uint64_t i;
        const struct mli_Vec direction = mli_Vec_init(0., 0., 1.);
        for (i = 0; i < num_photons; i++) {
                struct mli_Photon ph;
                ph.ray.support = mli_Vec_random_position_on_disc(radius, prng);
                ph.ray.direction = direction;
                ph.wavelength = wavelength;
                ph.id = i;
                chk(mli_PhotonVector_push_back(out_photons, ph));
        }
        return 1;
chk_error:
        return 0;
}

int mli_photon_source_point_like_opening_cone_towards_z(
        struct mli_PhotonVector *out_photons,
        const double wavelength,
        const double opening_angle,
        const uint64_t num_photons,
        struct mli_Prng *prng)
{
        uint64_t i;
        struct mli_prng_UniformRange azimuth =
                mli_prng_UniformRange_set(0.0, 2.0 * MLI_MATH_PI);
        struct mli_prng_ZenithRange zenith =
                mli_prng_ZenithRange_set(0.0, opening_angle);
        for (i = 0; i < num_photons; i++) {
                struct mli_Vec direction =
                        mli_Vec_random_draw_direction_in_zenith_azimuth_range(
                                zenith, azimuth, prng);
                struct mli_Photon ph;
                ph.ray.support = mli_Vec_init(0., 0., 0.);
                ph.ray.direction = direction;
                ph.wavelength = wavelength;
                ph.id = i;
                chk(mli_PhotonVector_push_back(out_photons, ph));
        }
        return 1;
chk_error:
        return 0;
}

/* photon_vector */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_VECTOR_IMPLEMENTATION(mli_PhotonVector, struct mli_Photon)

/* physics */
/* ------- */


double mli_physics_plancks_spectral_radiance_law_W_per_m2_per_sr_per_m(
        const double wavelength_m,
        const double temperature_K)
{
        const double c = MLI_PHYSICS_SPEED_OF_LIGHT_M_PER_S;
        const double k = MLI_PHYSICS_BOLTZMANN_J_PER_K;
        const double h = MLI_PHYSICS_PLANCK_KG_M2_PER_S;
        const double T = temperature_K;
        const double l = wavelength_m;
        const double A = (2.0 * h * pow(c, 2.0)) / pow(l, 5.0);
        const double _denominator = exp((h * c) / (l * k * T)) - 1.0;
        const double B = (1.0 / _denominator);
        return A * B;
}

/* pinhole */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_camera_PinHole mli_camera_PinHole_set(
        const double field_of_view,
        const struct mli_Image *image,
        const double row_over_column_pixel_ratio)
{
        struct mli_camera_PinHole out;

        assert(field_of_view > 0.);
        assert(field_of_view < mli_math_deg2rad(180.));

        out.optical_axis = mli_Vec_init(0.0, 0.0, 1.0);
        out.col_axis = mli_Vec_init(1.0, 0.0, 0.0);
        out.row_axis = mli_Vec_init(0.0, row_over_column_pixel_ratio, 0.0);

        out.distance_to_principal_point =
                ((0.5 * mli_Image_num_cols(image)) / tan(0.5 * field_of_view));
        out.principal_point = mli_Vec_multiply(
                out.optical_axis, out.distance_to_principal_point);
        return out;
}

struct mli_Ray mli_camera_PinHole_ray_at_row_col(
        const struct mli_camera_PinHole *self,
        const struct mli_Image *image,
        const uint32_t row,
        const uint32_t col)
{
        const struct mli_Vec pin_hole_position = mli_Vec_init(0.0, 0.0, 0.0);
        int row_idx_on_sensor = row - mli_Image_num_rows(image) / 2;
        int col_idx_on_sensor = col - mli_Image_num_cols(image) / 2;
        struct mli_Vec s_row =
                mli_Vec_multiply(self->row_axis, row_idx_on_sensor);
        struct mli_Vec s_col =
                mli_Vec_multiply(self->col_axis, col_idx_on_sensor);
        struct mli_Vec sensor_intersection = self->principal_point;
        sensor_intersection = mli_Vec_add(sensor_intersection, s_row);
        sensor_intersection = mli_Vec_add(sensor_intersection, s_col);
        return mli_Ray_set(pin_hole_position, sensor_intersection);
}

void mli_camera_PinHole_render_image(
        struct mli_camera_PinHole self,
        const struct mli_HomTraComp camera2root_comp,
        const struct mli_Shader *shader,
        struct mli_Image *image,
        struct mli_Prng *prng)
{
        struct mli_image_PixelWalk walk = mli_image_PixelWalk_init();
        struct mli_HomTra camera2root =
                mli_HomTraComp_from_compact(camera2root_comp);
        uint32_t i;
        const uint32_t num_pixel = mli_Image_num_pixel(image);

        for (i = 0; i < num_pixel; i++) {
                struct mli_image_Pixel px =
                        mli_image_PixelWalk_get_Pixel(&walk, &image->geometry);
                struct mli_Ray ray_wrt_camera =
                        mli_camera_PinHole_ray_at_row_col(
                                &self, image, px.row, px.col);

                struct mli_Ray ray_wrt_root =
                        mli_HomTraComp_ray(&camera2root, ray_wrt_camera);

                struct mli_Color color =
                        mli_Shader_trace_ray(shader, ray_wrt_root, prng);
                mli_Image_set_by_col_row(image, px.col, px.row, color);
                mli_image_PixelWalk_walk(&walk, &image->geometry);
        }
}

void mli_camera_PinHole_render_image_with_view(
        const struct mli_View view,
        const struct mli_Shader *shader,
        struct mli_Image *image,
        const double row_over_column_pixel_ratio,
        struct mli_Prng *prng)
{
        struct mli_camera_PinHole camera = mli_camera_PinHole_set(
                view.field_of_view, image, row_over_column_pixel_ratio);
        struct mli_HomTraComp camera2root_comp = mli_View_to_HomTraComp(view);
        mli_camera_PinHole_render_image(
                camera, camera2root_comp, shader, image, prng);
}

/* prng */
/* ---- */


uint32_t mli_Prng_generate_uint32(struct mli_Prng *prng)
{
        return prng->generate_uint32(&prng->_storage);
}

void mli_Prng_reinit(struct mli_Prng *prng, const uint32_t seed)
{
        prng->reinit(&prng->_storage, seed);
}

/**
 *      Mersenne Twister 19937
 *      ----------------------
 */

struct mli_Prng mli_Prng_init_MT19937(const uint32_t seed)
{
        struct mli_Prng prng;
        prng._storage.mt19937 = mli_prng_MT19937_init(seed);
        prng.generate_uint32 = mli_Prng_MT19937_generate_uint32;
        prng.reinit = mli_Prng_MT19937_reinit;
        return prng;
}

uint32_t mli_Prng_MT19937_generate_uint32(void *mt)
{
        return mli_prng_MT19937_generate_uint32((struct mli_prng_MT19937 *)mt);
}

void mli_Prng_MT19937_reinit(void *mt, const uint32_t seed)
{
        mli_prng_MT19937_reinit((struct mli_prng_MT19937 *)mt, seed);
}

/**
 *      PCG32
 *      -----
 */

struct mli_Prng mli_Prng_init_PCG32(const uint32_t seed)
{
        struct mli_Prng prng;
        prng._storage.pcg32 = mli_prng_PCG32_init(seed);
        prng.generate_uint32 = mli_Prng_PCG32_generate_uint32;
        prng.reinit = mli_Prng_PCG32_reinit;
        return prng;
}

uint32_t mli_Prng_PCG32_generate_uint32(void *pcg)
{
        return mli_prng_PCG32_generate_uint32((struct mli_prng_PCG32 *)pcg);
}

void mli_Prng_PCG32_reinit(void *pcg, const uint32_t seed)
{
        mli_prng_PCG32_reinit((struct mli_prng_PCG32 *)pcg, seed);
}

/**
 *      API
 *      ---
 *
 */

double mli_Prng_uniform(struct mli_Prng *prng)
{
        uint32_t rn_int = mli_Prng_generate_uint32(prng);
        const double rn = (double)rn_int;
        const double max_uint32 = (double)UINT32_MAX;
        return rn / max_uint32;
}

double mli_Prng_expovariate(struct mli_Prng *prng, const double rate)
{
        /*      Sampling from a poisson distribution */
        return -log(mli_Prng_uniform(prng)) / rate;
}

double mli_Prng_normal_Irwin_Hall_approximation(
        struct mli_Prng *prng,
        const double mean,
        const double std)
{
        uint64_t i;
        double sum_of_12 = 0.;
        double std1 = 0.;
        for (i = 0; i < 12; i++) {
                sum_of_12 += mli_Prng_uniform(prng);
        }
        std1 = sum_of_12 - 6.;
        return mean + std1 * std;
}

/*
        uniform linear range
        ====================
        Draw a uniform distribution within a limited range.
*/
struct mli_prng_UniformRange mli_prng_UniformRange_set(
        double start,
        double stop)
{
        struct mli_prng_UniformRange p;
        p.start = start;
        p.range = stop - start;
        assert(p.range >= 0.0);
        return p;
}

double mli_prng_draw_uniform(
        const struct mli_prng_UniformRange uniform_range,
        struct mli_Prng *prng)
{
        return uniform_range.range * mli_Prng_uniform(prng) +
               uniform_range.start;
}

/*
        zenith range
        ============
        Draw zenith-distances for an even distribution of points on a sphere.
*/
struct mli_prng_ZenithRange mli_prng_ZenithRange_set(
        const double min_zenith_distance,
        const double max_zenith_distance)
{
        struct mli_prng_ZenithRange zp;
        zp.z_min = (cos(min_zenith_distance) + 1.0) / 2.0;
        zp.z_range = (cos(max_zenith_distance) + 1.0) / 2.0 - zp.z_min;
        return zp;
}

double mli_prng_draw_zenith(
        const struct mli_prng_ZenithRange range,
        struct mli_Prng *prng)
{
        const double z = (range.z_range * mli_Prng_uniform(prng)) + range.z_min;
        return acos(2.0 * z - 1.0);
}

/* prng_MT19937 */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/*
 *      Adopted from https://en.wikipedia.org/wiki/Mersenne_Twister
 */
void mli_prng_MT19937_set_constants(struct mli_prng_MT19937 *mt)
{
        /*
         *      Define MT19937 constants (32-bit RNG)
         *      Assumes W = 32 (omitting this)
         */
        mt->N = 624, mt->M = 397;
        mt->R = 31;
        mt->A = 0x9908B0DF;
        mt->F = 1812433253;
        mt->U = 11;
        /*       Assumes D = 0xFFFFFFFF (omitting this) */
        mt->S = 7;
        mt->B = 0x9D2C5680;
        mt->T = 15;
        mt->C = 0xEFC60000;
        mt->L = 18;
        mt->MASK_LOWER = (1u << mt->R) - 1;
        mt->MASK_UPPER = (1u << mt->R);
}

void mli_prng_MT19937_reinit(struct mli_prng_MT19937 *mt, const uint32_t seed)
{
        uint32_t i;
        mli_prng_MT19937_set_constants(mt);
        mt->mt[0] = seed;
        for (i = 1; i < mt->N; i++) {
                mt->mt[i] =
                        (mt->F * (mt->mt[i - 1] ^ (mt->mt[i - 1] >> 30)) + i);
        }
        mt->index = mt->N;
}

struct mli_prng_MT19937 mli_prng_MT19937_init(const uint32_t seed)
{
        struct mli_prng_MT19937 mt;
        mli_prng_MT19937_reinit(&mt, seed);
        return mt;
}

void mli_prng_MT19937_twist(struct mli_prng_MT19937 *mt)
{
        uint32_t i, x, xA;
        for (i = 0; i < mt->N; i++) {
                x = (mt->mt[i] & mt->MASK_UPPER) +
                    (mt->mt[(i + 1) % mt->N] & mt->MASK_LOWER);
                xA = x >> 1;
                if (x & 0x1) {
                        xA ^= mt->A;
                }
                mt->mt[i] = mt->mt[(i + mt->M) % mt->N] ^ xA;
        }
        mt->index = 0;
}

uint32_t mli_prng_MT19937_generate_uint32(struct mli_prng_MT19937 *mt)
{
        uint32_t y;
        int i = mt->index;
        if (mt->index >= mt->N) {
                mli_prng_MT19937_twist(mt);
                i = mt->index;
        }
        y = mt->mt[i];
        mt->index = i + 1;
        y ^= (mt->mt[i] >> mt->U);
        y ^= (y << mt->S) & mt->B;
        y ^= (y << mt->T) & mt->C;
        y ^= (y >> mt->L);
        return y;
}

/* prng_PCG32 */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_prng_PCG32 mli_prng_PCG32_init(const uint32_t seed)
{
        struct mli_prng_PCG32 pcg32;
        mli_prng_pcg_setseq_64_srandom_r(&pcg32.state_setseq_64, seed, 0u);
        return pcg32;
}

uint32_t mli_prng_PCG32_generate_uint32(struct mli_prng_PCG32 *pcg32)
{
        return mli_prng_pcg_setseq_64_xsh_rr_32_random_r(
                &pcg32->state_setseq_64);
}

void mli_prng_PCG32_reinit(struct mli_prng_PCG32 *pcg32, const uint32_t seed)
{
        mli_prng_pcg_setseq_64_srandom_r(&pcg32->state_setseq_64, seed, 0u);
}

/* prng_pcg_variants_32bit_subset */
/* ------------------------------ */

/**
 * PCG Random Number Generation for C.
 *
 * Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *     http://www.pcg-random.org
 */


#define MLI_PRNG_PCG_DEFAULT_MULTIPLIER_64 6364136223846793005U
#define MLI_PRNG_PCG_DEFAULT_INCREMENT_64 1442695040888963407U

/**     Rotate helper functions.
 */

uint32_t mli_prng_pcg_rotr_32(uint32_t value, unsigned int rot)
{
        return (value >> rot) | (value << ((-rot) & 31));
}

/**     Output functions.  These are the core of the PCG generation scheme.
 *
 *      XSH RR
 */

uint32_t mli_prng_pcg_output_xsh_rr_64_32(uint64_t state)
{
        return mli_prng_pcg_rotr_32(
                ((state >> 18u) ^ state) >> 27u, state >> 59u);
}

/**     Functions to advance the underlying LCG.
 *      These functions are considered semi-private.
 *      There is rarely a good reason to call them directly.
 */

void mli_prng_pcg_setseq_64_step_r(struct mli_prng_pcg_state_setseq_64 *rng)
{
        rng->state = rng->state * MLI_PRNG_PCG_DEFAULT_MULTIPLIER_64 + rng->inc;
}

/**     Seed the RNG state.
 */

void mli_prng_pcg_setseq_64_srandom_r(
        struct mli_prng_pcg_state_setseq_64 *rng,
        uint64_t initstate,
        uint64_t initseq)
{
        rng->state = 0U;
        rng->inc = (initseq << 1u) | 1u;
        mli_prng_pcg_setseq_64_step_r(rng);
        rng->state += initstate;
        mli_prng_pcg_setseq_64_step_r(rng);
}

/**     Now, finally we provide
 *      a random_r function that provides a random number of the appropriate
 *      type (using the full range of the type).
 *
 *      Generation functions for XSH RR
 */

uint32_t mli_prng_pcg_setseq_64_xsh_rr_32_random_r(
        struct mli_prng_pcg_state_setseq_64 *rng)
{
        uint64_t oldstate = rng->state;
        mli_prng_pcg_setseq_64_step_r(rng);
        return mli_prng_pcg_output_xsh_rr_64_32(oldstate);
}

/* quaternion */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Quaternion mli_Quaternion_set(
        const double w,
        const double x,
        const double y,
        const double z)
{
        struct mli_Quaternion o;
        o.w = w;
        o.x = x;
        o.y = y;
        o.z = z;
        return o;
}

struct mli_Quaternion mli_Quaternion_set_unit_xyz(
        const double x,
        const double y,
        const double z)
{
        /*
         *       Recover 4th element: q.w.
         *       Expect unit-quaternion:
         *       1.0 != q.w**2 + q.x**2 + q.y**2 + q.z**2
         *       thus:
         *       q.w**2 = 1.0 - q.x**2 - q.y**2 - q.z**2
         *       q.w = sqrt(1.0 - q.x**2 - q.y**2 - q.z**2)
         */
        const double w = sqrt(1.0 - x * x - y * y - z * z);
        return mli_Quaternion_set(w, x, y, z);
}

int mli_Quaternion_equal(
        const struct mli_Quaternion a,
        const struct mli_Quaternion b)
{
        if (fabs(a.w - b.w) > DBL_EPSILON)
                return 0;
        if (fabs(a.x - b.x) > DBL_EPSILON)
                return 0;
        if (fabs(a.y - b.y) > DBL_EPSILON)
                return 0;
        if (fabs(a.z - b.z) > DBL_EPSILON)
                return 0;
        return 1;
}

int mli_Quaternion_equal_margin(
        const struct mli_Quaternion a,
        const struct mli_Quaternion b,
        const double margin)
{
        if (fabs(a.w - b.w) >= margin) {
                return 0;
        }
        if (fabs(a.x - b.x) >= margin) {
                return 0;
        }
        if (fabs(a.y - b.y) >= margin) {
                return 0;
        }
        if (fabs(a.z - b.z) >= margin) {
                return 0;
        }
        return 1;
}

struct mli_Quaternion mli_Quaternion_complex_conjugate(
        const struct mli_Quaternion q)
{
        struct mli_Quaternion c;
        c.w = q.w;
        c.x = -q.x;
        c.y = -q.y;
        c.z = -q.z;
        return c;
}

struct mli_Quaternion mli_Quaternion_product(
        const struct mli_Quaternion p,
        const struct mli_Quaternion q)
{
        struct mli_Quaternion pq;
        const struct mli_Vec P = mli_Vec_init(p.x, p.y, p.z);
        const struct mli_Vec Q = mli_Vec_init(q.x, q.y, q.z);
        const struct mli_Vec P_cross_Q = mli_Vec_cross(P, Q);
        pq.w = p.w * q.w - mli_Vec_dot(P, Q);
        pq.x = p.w * Q.x + q.w * P.x + P_cross_Q.x;
        pq.y = p.w * Q.y + q.w * P.y + P_cross_Q.y;
        pq.z = p.w * Q.z + q.w * P.z + P_cross_Q.z;
        return pq;
}

double mli_Quaternion_product_complex_conjugate(const struct mli_Quaternion p)
{
        return p.w * p.w + p.x * p.x + p.y * p.y + p.z * p.z;
}

double mli_Quaternion_norm(const struct mli_Quaternion q)
{
        return sqrt(mli_Quaternion_product_complex_conjugate(q));
}

struct mli_Quaternion mli_Quaternion_set_rotaxis_and_angle(
        const struct mli_Vec rot_axis,
        const double angle)
{
        const struct mli_Vec normed_rot_axis = mli_Vec_normalized(rot_axis);
        struct mli_Quaternion quat;
        const double angle_half = .5 * angle;
        const double sin_angle_half = sin(angle_half);
        quat.w = cos(angle_half);
        quat.x = normed_rot_axis.x * sin_angle_half;
        quat.y = normed_rot_axis.y * sin_angle_half;
        quat.z = normed_rot_axis.z * sin_angle_half;
        return quat;
}

struct mli_Mat mli_Quaternion_to_matrix(const struct mli_Quaternion quat)
{
        struct mli_Mat o;
        const double w2 = quat.w * quat.w;
        const double x2 = quat.x * quat.x;
        const double y2 = quat.y * quat.y;
        const double z2 = quat.z * quat.z;
        const double _2wx = 2. * quat.w * quat.x;
        const double _2wy = 2. * quat.w * quat.y;
        const double _2wz = 2. * quat.w * quat.z;
        /* const double 2xw */
        const double _2xy = 2. * quat.x * quat.y;
        const double _2xz = 2. * quat.x * quat.z;
        /* const double 2yw */
        /* const double 2yx */
        const double _2yz = 2. * quat.y * quat.z;
        o.r00 = w2 + x2 - y2 - z2;
        o.r01 = _2xy - _2wz;
        o.r02 = _2xz + _2wy;
        o.r10 = _2xy + _2wz;
        o.r11 = w2 - x2 + y2 - z2;
        o.r12 = _2yz - _2wx;
        o.r20 = _2xz - _2wy;
        o.r21 = _2yz + _2wx;
        o.r22 = w2 - x2 - y2 + z2;
        return o;
}

struct mli_Quaternion mli_Quaternion_set_tait_bryan(
        const double rx,
        const double ry,
        const double rz)
{
        const struct mli_Quaternion qz = mli_Quaternion_set_rotaxis_and_angle(
                mli_Vec_init(0, 0, 1), -rz);
        const struct mli_Quaternion qy = mli_Quaternion_set_rotaxis_and_angle(
                mli_Vec_init(0, 1, 0), -ry);
        const struct mli_Quaternion qx = mli_Quaternion_set_rotaxis_and_angle(
                mli_Vec_init(1, 0, 0), -rx);
        const struct mli_Quaternion qz_qy = mli_Quaternion_product(qz, qy);
        return mli_Quaternion_product(qz_qy, qx);
}

/* quaternion_json */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Quaternion_tait_bryan_from_json(
        struct mli_Quaternion *quat,
        const struct mli_Json *json,
        const uint64_t token)
{
        uint64_t token_xyz;
        struct mli_Vec xyz_deg;
        chk_msg(mli_Json_token_by_key(json, token, "xyz_deg", &token_xyz),
                "Expected tait_bryan to have key 'xyz_deg'.");
        chk_msg(mli_Vec_from_json_token(&xyz_deg, json, token_xyz + 1),
                "Failed to parse tait_bryan's 'xyz_deg' from json.");
        *quat = mli_Quaternion_set_tait_bryan(
                mli_math_deg2rad(xyz_deg.x),
                mli_math_deg2rad(xyz_deg.y),
                mli_math_deg2rad(xyz_deg.z));
        return 1;
chk_error:
        return 0;
}

int mli_Quaternion_axis_angle_from_json(
        struct mli_Quaternion *quat,
        const struct mli_Json *json,
        const uint64_t token)
{
        uint64_t token_axis, token_angle;
        double angle_deg;
        struct mli_Vec axis;
        chk_msg(mli_Json_token_by_key(json, token, "axis", &token_axis),
                "Expected axis_angle to have key 'axis'.");
        chk_msg(mli_Vec_from_json_token(&axis, json, token_axis + 1),
                "Failed to parse axis_angle's 'axis' from json.");
        chk_msg(mli_Json_token_by_key(json, token, "angle_deg", &token_angle),
                "Expected axis_angle to have key 'angle_deg'.");
        chk_msg(mli_Json_double_by_token(json, token_angle + 1, &angle_deg),
                "Failed to parse axis_angle's 'angle_deg' from json.");
        *quat = mli_Quaternion_set_rotaxis_and_angle(
                axis, mli_math_deg2rad(angle_deg));
        return 1;
chk_error:
        return 0;
}

int mli_Quaternion_quaternion_from_json(
        struct mli_Quaternion *quat,
        const struct mli_Json *json,
        const uint64_t token)
{
        uint64_t token_xyz;
        struct mli_Vec q;
        chk_msg(mli_Json_token_by_key(json, token, "xyz", &token_xyz),
                "Expected quaternion to have key 'xyz'.");
        chk_msg(mli_Vec_from_json_token(&q, json, token_xyz + 1),
                "Failed to parse quaternion's 'xyz' from json.");
        *quat = mli_Quaternion_set_unit_xyz(q.x, q.y, q.z);
        chk_msg(fabs(mli_Quaternion_norm(*quat) - 1.) < 1e-6,
                "Expected norm(quaternion) < 1e-6. Expected unit-quaternion.");
        return 1;
chk_error:
        return 0;
}

int mli_Quaternion_from_json(
        struct mli_Quaternion *quat,
        const struct mli_Json *json,
        const uint64_t token)
{
        uint64_t token_repr = 0u;
        uint64_t token_repr_value = 0u;
        chk_msg(mli_Json_token_by_key(json, token, "repr", &token_repr),
                "Expected 'rot' to have key 'repr'.");
        token_repr_value = token_repr + 1;

        if (mli_Json_cstrcmp(json, token_repr_value, "tait_bryan")) {
                chk_msg(mli_Quaternion_tait_bryan_from_json(quat, json, token),
                        "Failed to parse tait_bryan rotation.");
        } else if (mli_Json_cstrcmp(json, token_repr_value, "axis_angle")) {
                chk_msg(mli_Quaternion_axis_angle_from_json(quat, json, token),
                        "Failed to parse axis_angle rotation.");
        } else if (mli_Json_cstrcmp(json, token_repr_value, "quaternion")) {
                chk_msg(mli_Quaternion_quaternion_from_json(quat, json, token),
                        "Failed to parse quaternion rotation.");
        } else {
                chk_bad("Unknown representation 'repr' in rotation.");
        }
        return 1;
chk_error:
        return 0;
}

/* ray */
/* --- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Ray mli_Ray_set(
        const struct mli_Vec support,
        const struct mli_Vec direction)
{
        struct mli_Ray ray;
        ray.support = support;
        ray.direction =
                mli_Vec_multiply(direction, 1. / mli_Vec_norm(direction));
        return ray;
}

struct mli_Vec mli_Ray_at(const struct mli_Ray *ray, const double t)
{
        struct mli_Vec out;
        out.x = ray->support.x + t * ray->direction.x;
        out.y = ray->support.y + t * ray->direction.y;
        out.z = ray->support.z + t * ray->direction.z;
        return out;
}

int mli_Ray_sphere_intersection(
        const struct mli_Vec support,
        const struct mli_Vec direction,
        const double radius,
        double *minus_solution,
        double *plus_solution)
{
        const double sup_times_dir = mli_Vec_dot(support, direction);
        const double dir_times_dir = mli_Vec_dot(direction, direction);
        const double sup_times_sup = mli_Vec_dot(support, support);
        const double radius_square = radius * radius;

        const double p = 2.0 * (sup_times_dir / dir_times_dir);
        const double q = sup_times_sup / dir_times_dir - radius_square;

        return mli_math_quadratic_equation(p, q, minus_solution, plus_solution);
}

/* ray_AABB */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mli_Ray_aabb_intersections(
        const struct mli_Ray ray,
        const struct mli_AABB aabb,
        double *t_near,
        double *t_far)
{
        struct mli_Vec frac;
        struct mli_Vec lower, upper;
        struct mli_Vec t1, t2;

        frac.x = 1. / ray.direction.x;
        frac.y = 1. / ray.direction.y;
        frac.z = 1. / ray.direction.z;

        lower.x = (aabb.lower.x - ray.support.x) * frac.x;
        lower.y = (aabb.lower.y - ray.support.y) * frac.y;
        lower.z = (aabb.lower.z - ray.support.z) * frac.z;

        upper.x = (aabb.upper.x - ray.support.x) * frac.x;
        upper.y = (aabb.upper.y - ray.support.y) * frac.y;
        upper.z = (aabb.upper.z - ray.support.z) * frac.z;

        t1.x = MLI_MATH_MIN2(lower.x, upper.x);
        t1.y = MLI_MATH_MIN2(lower.y, upper.y);
        t1.z = MLI_MATH_MIN2(lower.z, upper.z);

        t2.x = MLI_MATH_MAX2(lower.x, upper.x);
        t2.y = MLI_MATH_MAX2(lower.y, upper.y);
        t2.z = MLI_MATH_MAX2(lower.z, upper.z);

        (*t_near) = MLI_MATH_MAX3(t1.x, t1.y, t1.z);
        (*t_far) = MLI_MATH_MIN3(t2.x, t2.y, t2.z);
}

int mli_Ray_aabb_intersections_is_valid_given_near_and_far(
        const double t_near,
        const double t_far)
{
        if (t_far < 0) {
                return 0;
        }
        if (t_near > t_far) {
                return 0;
        }
        return 1;
}

int mli_Ray_has_overlap_aabb(
        const struct mli_Ray ray,
        const struct mli_AABB aabb)
{
        double t_near, t_far;
        mli_Ray_aabb_intersections(ray, aabb, &t_near, &t_far);
        return mli_Ray_aabb_intersections_is_valid_given_near_and_far(
                t_near, t_far);
}

/* ray_octree_traversal */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/*
 *  Based on
 *
 *  @article{revelles2000efficient,
 *      title={An efficient parametric algorithm for octree traversal},
 *      author={Revelles, Jorge and Urena, Carlos and Lastra, Miguel},
 *      year={2000},
 *      publisher={V{\'a}clav Skala-UNION Agency}
 *  }
 *
 *  Profiting from additional comments by Jeroen Baert.
 */

/*
 *                    __ Z
 *        Y            /|
 *         /\         /
 *         |   -----------------
 *         |  /   3   /   7   /|
 *         | /---------------/ |
 *         |/   2   /   6   /|7|
 *         |---------------- | |
 *         |       |       |6|/|
 *         |   2   |   6   | / |
 *         |       |       |/|5|
 *         |-------|-------| | |
 *         |       |       |4|/
 *         |   0   |   4   | /
 *         |       |       |/
 *         L--------------------> X
 *
 *      Figure 1: Labeled Octree (the hidden one has label 1).
 *      cite{revelles2000efficient}
 */

/*  Current  |     Exit-plane
 *  sub-node |
 *  (state)  |  YZ      XZ      XY
 * ----------|----------------------
 *      0    |  4       2       1
 *      1    |  5       3       END
 *      2    |  6       END     3
 *      3    |  7       END     END
 *      4    |  END     6       5
 *      5    |  END     7       END
 *      6    |  END     END     7
 *      7    |  END     END     END
 *
 *      Table 3: State transition.
 *      cite{revelles2000efficient}
 */

int mli_raytracing_ray_octree_traversal_first_octree_node(
        const struct mli_Vec t0,
        const struct mli_Vec tm)
{
        uint8_t child = 0;
        if (t0.x > t0.y) {
                if (t0.x > t0.z) {
                        /* Y-Z-plane */
                        if (tm.y < t0.x)
                                child |= 2;
                        if (tm.z < t0.x)
                                child |= 1;
                        return (int)child;
                }
        } else {
                if (t0.y > t0.z) {
                        /* X-Z-plane */
                        if (tm.x < t0.y)
                                child |= 4;
                        if (tm.z < t0.y)
                                child |= 1;
                        return (int)child;
                }
        }
        /* X-Y-plane */
        if (tm.x < t0.z)
                child |= 4;
        if (tm.y < t0.z)
                child |= 2;
        return (int)child;
}

int mli_raytracing_ray_octree_traversal_next_octree_node(
        const struct mli_Vec tm,
        int x,
        int y,
        int z)
{
        if (tm.x < tm.y) {
                if (tm.x < tm.z) {
                        /* Y-Z-plane */
                        return x;
                }
        } else {
                if (tm.y < tm.z) {
                        /* X-Z-plane */
                        return y;
                }
        }
        return z; /* X-Y-plane */
}

void mli_raytracing_ray_octree_traversal_sub(
        struct mli_Vec t0,
        struct mli_Vec t1,
        const struct mli_OcTree *octree,
        const int32_t node_idx,
        const int32_t node_type,
        uint8_t permutation,
        void *work,
        void (*work_on_leaf_node)(
                void *,
                const struct mli_OcTree *,
                const uint32_t))
{
        int32_t proc_node;
        struct mli_Vec tm, nt0, nt1;

        if (t1.x < 0 || t1.y < 0 || t1.z < 0) {
                return;
        }

        if (node_type == MLI_OCTREE_TYPE_NONE) {
                return;
        }

        if (node_type == MLI_OCTREE_TYPE_LEAF) {
                /* callback */
                work_on_leaf_node(work, octree, node_idx);
                return;
        }

        tm.x = 0.5 * (t0.x + t1.x);
        tm.y = 0.5 * (t0.y + t1.y);
        tm.z = 0.5 * (t0.z + t1.z);

        proc_node =
                mli_raytracing_ray_octree_traversal_first_octree_node(t0, tm);

        do {
                switch (proc_node) {
                case 0: {
                        nt0 = mli_Vec_init(t0.x, t0.y, t0.z);
                        nt1 = mli_Vec_init(tm.x, tm.y, tm.z);
                        mli_raytracing_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx].children[permutation],
                                octree->nodes[node_idx].types[permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node =
                                mli_raytracing_ray_octree_traversal_next_octree_node(
                                        nt1, 4, 2, 1);
                        break;
                }
                case 1: {
                        nt0 = mli_Vec_init(t0.x, t0.y, tm.z);
                        nt1 = mli_Vec_init(tm.x, tm.y, t1.z);
                        mli_raytracing_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[1 ^ permutation],
                                octree->nodes[node_idx].types[1 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node =
                                mli_raytracing_ray_octree_traversal_next_octree_node(
                                        nt1, 5, 3, 8);
                        break;
                }
                case 2: {
                        nt0 = mli_Vec_init(t0.x, tm.y, t0.z);
                        nt1 = mli_Vec_init(tm.x, t1.y, tm.z);
                        mli_raytracing_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[2 ^ permutation],
                                octree->nodes[node_idx].types[2 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node =
                                mli_raytracing_ray_octree_traversal_next_octree_node(
                                        nt1, 6, 8, 3);
                        break;
                }
                case 3: {
                        nt0 = mli_Vec_init(t0.x, tm.y, tm.z);
                        nt1 = mli_Vec_init(tm.x, t1.y, t1.z);
                        mli_raytracing_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[3 ^ permutation],
                                octree->nodes[node_idx].types[3 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node =
                                mli_raytracing_ray_octree_traversal_next_octree_node(
                                        nt1, 7, 8, 8);
                        break;
                }
                case 4: {
                        nt0 = mli_Vec_init(tm.x, t0.y, t0.z);
                        nt1 = mli_Vec_init(t1.x, tm.y, tm.z);
                        mli_raytracing_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[4 ^ permutation],
                                octree->nodes[node_idx].types[4 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node =
                                mli_raytracing_ray_octree_traversal_next_octree_node(
                                        nt1, 8, 6, 5);
                        break;
                }
                case 5: {
                        nt0 = mli_Vec_init(tm.x, t0.y, tm.z);
                        nt1 = mli_Vec_init(t1.x, tm.y, t1.z);
                        mli_raytracing_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[5 ^ permutation],
                                octree->nodes[node_idx].types[5 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node =
                                mli_raytracing_ray_octree_traversal_next_octree_node(
                                        nt1, 8, 7, 8);
                        break;
                }
                case 6: {
                        nt0 = mli_Vec_init(tm.x, tm.y, t0.z);
                        nt1 = mli_Vec_init(t1.x, t1.y, tm.z);
                        mli_raytracing_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[6 ^ permutation],
                                octree->nodes[node_idx].types[6 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node =
                                mli_raytracing_ray_octree_traversal_next_octree_node(
                                        nt1, 8, 8, 7);
                        break;
                }
                case 7: {
                        nt0 = mli_Vec_init(tm.x, tm.y, tm.z);
                        nt1 = mli_Vec_init(t1.x, t1.y, t1.z);
                        mli_raytracing_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[7 ^ permutation],
                                octree->nodes[node_idx].types[7 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node = 8;
                        break;
                }
                } /* end switch */
        } while (proc_node < 8);
}

void mli_raytracing_ray_octree_traversal(
        const struct mli_OcTree *octree,
        const struct mli_Ray ray,
        void *work,
        void (*work_on_leaf_node)(
                void *,
                const struct mli_OcTree *,
                const uint32_t))
{
        struct mli_Vec t0, t1;
        struct mli_Vec div;
        struct mli_Ray ray_wrt_octree;
        struct mli_Vec cube_upper, cube_size;
        struct mli_Cube cube;
        int32_t octree_root_node, octree_root_type;
        uint8_t permutation = 0;
        cube = octree->cube;
        cube_upper = mli_Cube_upper(cube);
        octree_root_node = 0u;
        octree_root_type = octree->root_type;
        cube_size = mli_Vec_add(cube.lower, cube_upper);

        ray_wrt_octree = ray;

        if (ray_wrt_octree.direction.x < 0) {
                ray_wrt_octree.support.x =
                        -ray_wrt_octree.support.x + cube_size.x;
                ray_wrt_octree.direction.x = -ray_wrt_octree.direction.x;
                permutation |= 4;
        } else if (ray_wrt_octree.direction.x == 0.0) {
                ray_wrt_octree.direction.x =
                        MLI_RAYTRACING_RAY_OCTREE_TRAVERSAL_EPSILON;
        }

        if (ray_wrt_octree.direction.y < 0) {
                ray_wrt_octree.support.y =
                        -ray_wrt_octree.support.y + cube_size.y;
                ray_wrt_octree.direction.y = -ray_wrt_octree.direction.y;
                permutation |= 2;
        } else if (ray_wrt_octree.direction.y == 0.0) {
                ray_wrt_octree.direction.y =
                        MLI_RAYTRACING_RAY_OCTREE_TRAVERSAL_EPSILON;
        }

        if (ray_wrt_octree.direction.z < 0) {
                ray_wrt_octree.support.z =
                        -ray_wrt_octree.support.z + cube_size.z;
                ray_wrt_octree.direction.z = -ray_wrt_octree.direction.z;
                permutation |= 1;
        } else if (ray_wrt_octree.direction.z == 0.0) {
                ray_wrt_octree.direction.z =
                        MLI_RAYTRACING_RAY_OCTREE_TRAVERSAL_EPSILON;
        }

        div.x = 1.0 / ray_wrt_octree.direction.x;
        div.y = 1.0 / ray_wrt_octree.direction.y;
        div.z = 1.0 / ray_wrt_octree.direction.z;

        t0.x = (cube.lower.x - ray_wrt_octree.support.x) * div.x;
        t1.x = (cube_upper.x - ray_wrt_octree.support.x) * div.x;

        t0.y = (cube.lower.y - ray_wrt_octree.support.y) * div.y;
        t1.y = (cube_upper.y - ray_wrt_octree.support.y) * div.y;

        t0.z = (cube.lower.z - ray_wrt_octree.support.z) * div.z;
        t1.z = (cube_upper.z - ray_wrt_octree.support.z) * div.z;

        if (MLI_MATH_MAX3(t0.x, t0.y, t0.z) < MLI_MATH_MIN3(t1.x, t1.y, t1.z)) {
                mli_raytracing_ray_octree_traversal_sub(
                        t0,
                        t1,
                        octree,
                        octree_root_node,
                        octree_root_type,
                        permutation,
                        work,
                        work_on_leaf_node);
        }
}

/* ray_scenery_query */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mli_raytracing_inner_object_traversal(
        void *_inner,
        const struct mli_OcTree *object_octree,
        const uint32_t object_octree_leaf_idx)
{
        /* traverse faces in an object-wavefront */
        struct mli_raytracing_QueryInnerWork *inner =
                (struct mli_raytracing_QueryInnerWork *)_inner;

        uint32_t f;
        const uint32_t num_faces_in_object_leaf = mli_OcTree_leaf_num_objects(
                object_octree, object_octree_leaf_idx);

        struct mli_Intersection tmp_isec = mli_Intersection_init();

        for (f = 0; f < num_faces_in_object_leaf; f++) {

                uint32_t face_idx = mli_OcTree_leaf_object_link(
                        object_octree, object_octree_leaf_idx, f);

                struct mli_object_Face fv =
                        inner->object->faces_vertices[face_idx];

                int32_t hit = mli_Ray_intersects_triangle(
                        inner->ray_object,
                        inner->object->vertices[fv.a],
                        inner->object->vertices[fv.b],
                        inner->object->vertices[fv.c],
                        &tmp_isec.distance_of_ray);

                tmp_isec.geometry_id.face = face_idx;

                if (hit) {
                        tmp_isec.position_local = mli_Ray_at(
                                &inner->ray_object, tmp_isec.distance_of_ray);

                        inner->has_intersection = 1;
                        if (tmp_isec.distance_of_ray <
                            inner->intersection->distance_of_ray) {
                                (*inner->intersection) = tmp_isec;
                        }
                }
        }
        return;
}

int mli_raytracing_query_object_reference(
        const struct mli_Object *object,
        const struct mli_OcTree *object_octree,
        const struct mli_HomTraComp robject2root_comp,
        const struct mli_Ray ray_root,
        struct mli_Intersection *isec)
{
        struct mli_HomTra robject2root =
                mli_HomTraComp_from_compact(robject2root_comp);

        struct mli_raytracing_QueryInnerWork inner;
        inner.has_intersection = 0;
        inner.intersection = isec;
        inner.ray_object = mli_HomTraComp_ray_inverse(&robject2root, ray_root);
        inner.object = object;

        mli_raytracing_ray_octree_traversal(
                object_octree,
                inner.ray_object,
                (void *)&inner,
                mli_raytracing_inner_object_traversal);

        return inner.has_intersection;
}

void mli_raytracing_outer_scenery_traversal(
        void *_outer,
        const struct mli_OcTree *scenery_octree,
        const uint32_t scenery_octree_leaf_idx)
{
        /* traverse object-wavefronts in a scenery */
        struct mli_raytracing_QueryOuterWork *outer =
                (struct mli_raytracing_QueryOuterWork *)_outer;

        uint32_t ro;
        const uint32_t num_robjects_in_scenery_leaf =
                mli_OcTree_leaf_num_objects(
                        scenery_octree, scenery_octree_leaf_idx);

        struct mli_Intersection tmp_isec = mli_Intersection_init();

        for (ro = 0; ro < num_robjects_in_scenery_leaf; ro++) {

                uint32_t robject_idx = mli_OcTree_leaf_object_link(
                        scenery_octree, scenery_octree_leaf_idx, ro);
                uint32_t object_idx = outer->geometry->robjects[robject_idx];

                int32_t hit = mli_raytracing_query_object_reference(
                        &outer->geometry->objects[object_idx],
                        &outer->accelerator->object_octrees[object_idx],
                        outer->geometry->robject2root[robject_idx],
                        outer->ray_root,
                        &tmp_isec);

                tmp_isec.geometry_id.robj = robject_idx;

                if (hit) {
                        if (tmp_isec.distance_of_ray <
                            outer->intersection->distance_of_ray) {
                                (*outer->intersection) = tmp_isec;
                        }
                }
        }
        return;
}

int mli_raytracing_query_intersection(
        const struct mli_Scenery *scenery,
        const struct mli_Ray ray_root,
        struct mli_Intersection *isec)
{
        struct mli_raytracing_QueryOuterWork outer;

        (*isec) = mli_Intersection_init();

        outer.intersection = isec;
        outer.geometry = &scenery->geometry;
        outer.accelerator = &scenery->accelerator;
        outer.ray_root = ray_root;

        mli_raytracing_ray_octree_traversal(
                &scenery->accelerator.scenery_octree,
                ray_root,
                (void *)&outer,
                mli_raytracing_outer_scenery_traversal);

        if (isec->distance_of_ray == DBL_MAX) {
                return 0;
        } else {
                return 1;
        }
}

int mli_raytracing_query_intersection_with_surface_normal(
        const struct mli_Scenery *scenery,
        const struct mli_Ray ray_root,
        struct mli_IntersectionSurfaceNormal *isecsrf)
{
        struct mli_Intersection isec = mli_Intersection_init();

        const int has_intersection =
                mli_raytracing_query_intersection(scenery, ray_root, &isec);

        if (has_intersection) {
                uint32_t robject_idx = isec.geometry_id.robj;
                uint32_t object_idx =
                        scenery->geometry.robjects[isec.geometry_id.robj];
                uint32_t face_idx = isec.geometry_id.face;

                struct mli_HomTra robject2root = mli_HomTraComp_from_compact(
                        scenery->geometry.robject2root[robject_idx]);
                struct mli_Ray ray_object =
                        mli_HomTraComp_ray_inverse(&robject2root, ray_root);

                struct mli_Object *obj = &scenery->geometry.objects[object_idx];

                struct mli_object_Face fv = obj->faces_vertices[face_idx];
                struct mli_object_Face fvn =
                        obj->faces_vertex_normals[face_idx];

                (*isecsrf) = mli_IntersectionSurfaceNormal_init();
                isecsrf->distance_of_ray = isec.distance_of_ray;
                isecsrf->geometry_id = isec.geometry_id;
                isecsrf->position = mli_Ray_at(&ray_root, isec.distance_of_ray);
                isecsrf->position_local = isec.position_local;

                /* find surface-normal */
                isecsrf->surface_normal_local = mli_Triangle_surface_normal(
                        obj->vertices[fv.a],
                        obj->vertices[fv.b],
                        obj->vertices[fv.c],
                        obj->vertex_normals[fvn.a],
                        obj->vertex_normals[fvn.b],
                        obj->vertex_normals[fvn.c],
                        isecsrf->position_local);

                isecsrf->surface_normal = mli_HomTraComp_dir(
                        &robject2root, isecsrf->surface_normal_local);

                isecsrf->from_outside_to_inside =
                        mli_raytracing_from_outside_to_inside(
                                ray_object.direction,
                                isecsrf->surface_normal_local);

                return 1;
        } else {
                return 0;
        }
}

/* ray_voxel_overlap */
/* ----------------- */

/* Copyright 2017 Sebastian A. Mueller */

void intersection_plane(
        double *support,
        double *direction,
        double off,
        double *intersection,
        unsigned int dim)
{
        double a = (off - support[dim]) / direction[dim];
        intersection[0] = support[0] + direction[0] * a;
        intersection[1] = support[1] + direction[1] * a;
        intersection[2] = support[2] + direction[2] * a;
}

double distance(double *v, double *u)
{
        double dx = v[0] - u[0];
        double dy = v[1] - u[1];
        double dz = v[2] - u[2];
        return sqrt(dx * dx + dy * dy + dz * dz);
}

double c_ray_box_overlap(
        double *support,
        double *direction,
        double xl,
        double xu,
        double yl,
        double yu,
        double zl,
        double zu)
{
        double hits_l[3][3];
        unsigned int nl = 0;

        double hits_u[3][3];
        unsigned int nu = 0;

        /* X plane */
        /* ------- */
        if (direction[0] != 0.0) {
                double ixl[3];
                double ixu[3];
                intersection_plane(support, direction, xl, ixl, 0);
                if ((ixl[1] >= yl && ixl[1] < yu) &&
                    (ixl[2] >= zl && ixl[2] < zu)) {
                        nl = nl + 1;
                        hits_l[nl - 1][0] = ixl[0];
                        hits_l[nl - 1][1] = ixl[1];
                        hits_l[nl - 1][2] = ixl[2];
                }

                intersection_plane(support, direction, xu, ixu, 0);
                if ((ixu[1] >= yl && ixu[1] <= yu) &&
                    (ixu[2] >= zl && ixu[2] <= zu)) {
                        nu = nu + 1;
                        hits_u[nu - 1][0] = ixu[0];
                        hits_u[nu - 1][1] = ixu[1];
                        hits_u[nu - 1][2] = ixu[2];
                }
        }

        /* Y plane */
        /* ------- */
        if (direction[1] != 0.0) {
                double iyl[3];
                double iyu[3];
                intersection_plane(support, direction, yl, iyl, 1);
                if ((iyl[0] >= xl && iyl[0] < xu) &&
                    (iyl[2] >= zl && iyl[2] < zu)) {
                        nl = nl + 1;
                        hits_l[nl - 1][0] = iyl[0];
                        hits_l[nl - 1][1] = iyl[1];
                        hits_l[nl - 1][2] = iyl[2];
                }

                intersection_plane(support, direction, yu, iyu, 1);
                if ((iyu[0] >= xl && iyu[0] <= xu) &&
                    (iyu[2] >= zl && iyu[2] <= zu)) {
                        nu = nu + 1;
                        hits_u[nu - 1][0] = iyu[0];
                        hits_u[nu - 1][1] = iyu[1];
                        hits_u[nu - 1][2] = iyu[2];
                }
        }

        /* Z plane */
        /* ------- */
        if (direction[2] != 0.0) {
                double izl[3];
                double izu[3];
                intersection_plane(support, direction, zl, izl, 2);
                if ((izl[0] >= xl && izl[0] < xu) &&
                    (izl[1] >= yl && izl[1] < yu)) {
                        nl = nl + 1;
                        hits_l[nl - 1][0] = izl[0];
                        hits_l[nl - 1][1] = izl[1];
                        hits_l[nl - 1][2] = izl[2];
                }

                intersection_plane(support, direction, zu, izu, 2);
                if ((izu[0] >= xl && izu[0] <= xu) &&
                    (izu[1] >= yl && izu[1] <= yu)) {
                        nu = nu + 1;
                        hits_u[nu - 1][0] = izu[0];
                        hits_u[nu - 1][1] = izu[1];
                        hits_u[nu - 1][2] = izu[2];
                }
        }

        if (nl == 2 && nu == 0) {
                return distance(hits_l[0], hits_l[1]);

        } else if (nl == 0 && nu == 2) {
                return distance(hits_u[0], hits_u[1]);

        } else if (nl == 1 && nu == 1) {
                return distance(hits_u[0], hits_l[0]);

        } else if (nl == 3 && nu == 3) {
                return distance(hits_u[0], hits_l[0]);

        } else if (nl == 2 && nu == 2) {
                return distance(hits_u[0], hits_l[0]);

        } else if (nl == 3 && nu == 0) {
                return sqrt(
                        (xu - xl) * (xu - xl) + (yu - yl) * (yu - yl) +
                        (zu - zl) * (zu - zl));

        } else if (nl > 0 && nu > 0) {
                return distance(hits_u[0], hits_l[0]);

        } else if (nl == 0 && nu == 3) {
                return sqrt(
                        (xu - xl) * (xu - xl) + (yu - yl) * (yu - yl) +
                        (zu - zl) * (zu - zl));

        } else {
                return 0.0;
        }
}

void c_next_space_partitions(
        unsigned int *dim_range,
        unsigned int dim_partitions[2][2],
        unsigned int *n)
{
        if (dim_range[1] - dim_range[0] < 1) {
                *n = *n + 1;
                dim_partitions[*n - 1][0] = dim_range[0];
                dim_partitions[*n - 1][1] = dim_range[1];
                return;
        } else {
                unsigned int cut = (dim_range[1] - dim_range[0]) / 2;
                *n = *n + 1;
                dim_partitions[*n - 1][0] = dim_range[0];
                dim_partitions[*n - 1][1] = dim_range[0] + cut;
                *n = *n + 1;
                dim_partitions[*n - 1][0] = dim_range[0] + cut;
                dim_partitions[*n - 1][1] = dim_range[1];
                return;
        }
}

void c_overlap_of_ray_with_voxels(
        double *support,
        double *direction,
        double *x_bin_edges,
        double *y_bin_edges,
        double *z_bin_edges,
        unsigned int *x_range,
        unsigned int *y_range,
        unsigned int *z_range,

        unsigned int *number_overlaps,
        unsigned int *x_idxs,
        unsigned int *y_idxs,
        unsigned int *z_idxs,
        double *overlaps)
{
        unsigned int x_partitions[2][2];
        unsigned int nxp = 0;

        unsigned int y_partitions[2][2];
        unsigned int nyp = 0;

        unsigned int z_partitions[2][2];
        unsigned int nzp = 0;

        unsigned int xp;
        unsigned int yp;
        unsigned int zp;

        c_next_space_partitions(x_range, x_partitions, &nxp);
        c_next_space_partitions(y_range, y_partitions, &nyp);
        c_next_space_partitions(z_range, z_partitions, &nzp);

        for (xp = 0; xp < nxp; xp = xp + 1) {
                for (yp = 0; yp < nxp; yp = yp + 1) {
                        for (zp = 0; zp < nxp; zp = zp + 1) {
                                double overlap = c_ray_box_overlap(
                                        support,
                                        direction,
                                        x_bin_edges[x_partitions[xp][0]],
                                        x_bin_edges[x_partitions[xp][1]],
                                        y_bin_edges[y_partitions[yp][0]],
                                        y_bin_edges[y_partitions[yp][1]],
                                        z_bin_edges[z_partitions[zp][0]],
                                        z_bin_edges[z_partitions[zp][1]]);

                                if (x_partitions[xp][1] - x_partitions[xp][0] ==
                                            1 &&
                                    y_partitions[yp][1] - y_partitions[yp][0] ==
                                            1 &&
                                    z_partitions[zp][1] - z_partitions[zp][0] ==
                                            1 &&
                                    overlap > 0.0) {
                                        *number_overlaps = *number_overlaps + 1;
                                        x_idxs[*number_overlaps - 1] =
                                                x_partitions[xp][0];
                                        y_idxs[*number_overlaps - 1] =
                                                y_partitions[yp][0];
                                        z_idxs[*number_overlaps - 1] =
                                                z_partitions[zp][0];
                                        overlaps[*number_overlaps - 1] =
                                                overlap;
                                } else if (overlap > 0.0) {
                                        c_overlap_of_ray_with_voxels(
                                                support,
                                                direction,
                                                x_bin_edges,
                                                y_bin_edges,
                                                z_bin_edges,
                                                x_partitions[xp],
                                                y_partitions[yp],
                                                z_partitions[zp],
                                                number_overlaps,
                                                x_idxs,
                                                y_idxs,
                                                z_idxs,
                                                overlaps);
                                }
                        }
                }
        }
        return;
}

void mli_corsika_overlap_of_ray_with_voxels(
        const struct mli_corsika_PhotonBunch *bunch,
        const struct mli_DoubleVector *x_bin_edges,
        const struct mli_DoubleVector *y_bin_edges,
        const struct mli_DoubleVector *z_bin_edges,
        struct mli_Uint32Vector *x_idxs,
        struct mli_Uint32Vector *y_idxs,
        struct mli_Uint32Vector *z_idxs,
        struct mli_DoubleVector *overlaps)
{
        double support[3];
        double direction[3];
        struct mli_Vec _direction = mli_Vec_init(0.0, 0.0, 0.0);

        unsigned int number_overlaps = 0u;
        unsigned int x_range[2];
        unsigned int y_range[2];
        unsigned int z_range[2];

        support[0] = bunch->x_cm;
        support[1] = bunch->y_cm;
        support[2] = 0.0;

        _direction = mli_corsika_photon_direction_of_motion(*bunch);
        direction[0] = _direction.x;
        direction[1] = _direction.y;
        direction[2] = _direction.z;

        x_range[0] = 0u;
        x_range[1] = x_bin_edges->size - 1u;

        y_range[0] = 0u;
        y_range[1] = y_bin_edges->size - 1u;

        z_range[0] = 0u;
        z_range[1] = z_bin_edges->size - 1u;

        x_idxs->size = 0;
        y_idxs->size = 0;
        z_idxs->size = 0;
        overlaps->size = 0;

        c_overlap_of_ray_with_voxels(
                support,
                direction,
                x_bin_edges->array,
                y_bin_edges->array,
                z_bin_edges->array,
                x_range,
                y_range,
                z_range,
                &number_overlaps,
                x_idxs->array,
                y_idxs->array,
                z_idxs->array,
                overlaps->array);

        assert(number_overlaps <= x_idxs->capacity);

        x_idxs->size = number_overlaps;
        y_idxs->size = number_overlaps;
        z_idxs->size = number_overlaps;
        overlaps->size = number_overlaps;
}

/* scenery */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Scenery mli_Scenery_init(void)
{
        struct mli_Scenery out;
        out.geometry = mli_Geometry_init();
        out.accelerator = mli_Accelerator_init();
        out.materials = mli_Materials_init();
        out.geomap = mli_GeometryToMaterialMap_init();
        return out;
}

void mli_Scenery_free(struct mli_Scenery *self)
{
        mli_Geometry_free(&self->geometry);
        mli_Accelerator_free(&self->accelerator);
        mli_Materials_free(&self->materials);
        mli_GeometryToMaterialMap_free(&self->geomap);
}

void mli_Scenery_info_fprint(FILE *f, const struct mli_Scenery *self)
{
        mli_Geometry_info_fprint(f, &self->geometry);
        fprintf(f, "\n");
        mli_Accelerator_info_fprint(f, &self->accelerator);
        fprintf(f, "\n");
        mli_Materials_info_fprint(f, &self->materials);
        fprintf(f, "\n");
        mli_GeometryToMaterialMap_info_fprint(f, &self->geomap);
        fprintf(f, "\n");
}

uint32_t mli_Scenery_resolve_boundary_layer_idx(
        const struct mli_Scenery *scenery,
        const struct mli_GeometryId geometry_id)
{
        const uint32_t robject_idx = geometry_id.robj;
        const uint32_t object_idx = scenery->geometry.robjects[robject_idx];
        const uint32_t face_idx = geometry_id.face;
        const uint32_t obj_mtl_idx = mli_Object_resolve_material_idx(
                &scenery->geometry.objects[object_idx], face_idx);
        const uint32_t boundary_layer_idx = mli_GeometryToMaterialMap_get(
                &scenery->geomap, robject_idx, obj_mtl_idx);
        return boundary_layer_idx;
}

/* scenery_equal */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Scenery_equal(const struct mli_Scenery *a, const struct mli_Scenery *b)
{
        chk_msg(mli_Geometry_equal(&a->geometry, &b->geometry),
                "Expected geometry to be valid.");
        chk_msg(mli_Materials_equal(&a->materials, &b->materials),
                "Expected materials to be valid.");
        chk_msg(mli_Accelerator_equal(&a->accelerator, &b->accelerator),
                "Expected accelerator to be valid.");
        chk_msg(mli_GeometryToMaterialMap_equal(&a->geomap, &b->geomap),
                "Expected geomap to be valid.");
        return 1;
chk_error:
        return 0;
}

/* scenery_minimal_object */
/* ---------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Func_malloc_color_spectrum(
        struct mli_Func *self,
        const struct mli_Color color)
{
        /* r: 600nm
         * g: 550nm
         * b: 450nm
         */

        chk(mli_Func_malloc(self, 11));
        self->x[0] = 200e-9;
        self->y[0] = 0.0;

        /* blue peak */
        self->x[1] = 400e-9;
        self->y[1] = 0.0;
        self->x[2] = 450e-9;
        self->y[2] = color.b;
        self->x[3] = 500e-9;
        self->y[3] = 0.0;

        /* green peak */
        self->x[4] = 501e-9;
        self->y[4] = 0.0;
        self->x[5] = 550e-9;
        self->y[5] = color.g;
        self->x[6] = 575e-9;
        self->y[6] = 0.0;

        /* red peak */
        self->x[7] = 576e-9;
        self->y[7] = 0.0;
        self->x[8] = 600e-9;
        self->y[8] = color.r;
        self->x[9] = 650e-9;
        self->y[9] = 0.0;

        self->x[10] = 1200e-9;
        self->y[10] = 0.0;

        return 1;
chk_error:
        return 0;
}

struct mli_Color mli_Color_random_uniform(struct mli_Prng *prng)
{
        struct mli_prng_UniformRange uniform_range;
        struct mli_Color out;
        uniform_range.start = 0.0;
        uniform_range.range = 1.0;

        out = mli_Color_set(
                mli_prng_draw_uniform(uniform_range, prng),
                mli_prng_draw_uniform(uniform_range, prng),
                mli_prng_draw_uniform(uniform_range, prng));
        return out;
}

int mli_Scenery_malloc_minimal_from_wavefront(
        struct mli_Scenery *self,
        const char *path)
{
        uint32_t i, total_num_boundary_layers;
        uint64_t spec;
        struct mli_Prng prng = mli_Prng_init_MT19937(1u);
        struct mli_IO str = mli_IO_init();
        struct mli_MaterialsCapacity mtlcap = mli_MaterialsCapacity_init();
        struct mli_String _path = mli_String_init();
        struct mli_String _mode = mli_String_init();
        struct mli_Spectrum *spectrum = NULL;
        struct mli_Surface *surface = NULL;
        struct mli_Medium *medium = NULL;
        struct mli_BoundaryLayer *layer = NULL;

        mli_Scenery_free(self);

        chk_msg(mli_Geometry_malloc(&self->geometry, 1u, 1u),
                "Failed to malloc geometry.");

        /* set object */
        chk(mli_String_from_cstr(&_path, path));
        chk(mli_String_from_cstr(&_mode, "r"));
        chk_msg(mli_IO_open_file(&str, &_path, &_mode), "Failed to open file.");
        chk_msg(mli_Object_malloc_from_wavefront(
                        &self->geometry.objects[0], &str),
                "Failed to malloc wavefront-object from string.");
        mli_IO_close(&str);

        mli_String_free(&_path);
        mli_String_free(&_mode);

        chk(mli_String_from_cstr(
                &self->geometry.object_names[0], "default-object"));

        /* set reference */
        self->geometry.robjects[0] = 0u;
        self->geometry.robject_ids[0] = 0u;
        self->geometry.robject2root[0] = mli_HomTraComp_set(
                mli_Vec_init(0.0, 0.0, 0.0),
                mli_Quaternion_set_tait_bryan(0.0, 0.0, 0.0));

        /* materials */
        total_num_boundary_layers = self->geometry.objects[0].num_materials;

        mtlcap.num_media = 1u;
        mtlcap.num_boundary_layers = total_num_boundary_layers;
        mtlcap.num_surfaces = total_num_boundary_layers;
        mtlcap.num_spectra = 2u + total_num_boundary_layers;

        chk_msg(mli_Materials_malloc(&self->materials, mtlcap),
                "Failed to malloc materials.");

        spec = 0;

        spectrum = &self->materials.spectra.array[spec];
        chk(mli_String_from_cstr(&spectrum->name, "vacuum_refraction"));
        chk(mli_Func_malloc_constant(
                &spectrum->spectrum, 200e-9, 1200e-9, 1.0));
        spec += 1;

        spectrum = &self->materials.spectra.array[spec];
        chk(mli_String_from_cstr(&spectrum->name, "vacuum_absorbtion"));
        chk(mli_Func_malloc_constant(
                &spectrum->spectrum, 200e-9, 1200e-9, 0.0));
        spec += 1;

        medium = &self->materials.media.array[0];
        mli_Medium_free(medium);
        chk(mli_String_from_cstr(&medium->name, "vacuum"));
        medium->refraction_spectrum = 0;
        medium->absorbtion_spectrum = 1;

        for (i = 0u; i < total_num_boundary_layers; i++) {
                spectrum = &self->materials.spectra.array[spec];
                mli_Spectrum_free(spectrum);
                chk(mli_String_from_cstr_fromat(
                        &spectrum->name, "reflection_%06u", i));
                chk(mli_Func_malloc_color_spectrum(
                        &spectrum->spectrum, mli_Color_random_uniform(&prng)));

                surface = &self->materials.surfaces.array[i];
                mli_Surface_free(surface);
                chk(mli_String_from_cstr_fromat(
                        &surface->name, "surface_%06u", i));
                surface->type = MLI_SURFACE_TYPE_COOKTORRANCE;
                surface->data.cooktorrance.reflection_spectrum = spec;
                surface->data.cooktorrance.diffuse_weight = 1.0;
                surface->data.cooktorrance.specular_weight = 0.0;
                surface->data.cooktorrance.roughness = 0.2;

                layer = &self->materials.boundary_layers.array[i];
                mli_BoundaryLayer_free(layer);
                chk(mli_String_from_cstr_fromat(&layer->name, "layer_%06u", i));
                layer->inner.medium = 0;
                layer->inner.surface = i;
                layer->outer.medium = 0;
                layer->outer.surface = i;

                spec += 1;
        }

        chk_msg(mli_GeometryToMaterialMap_malloc(
                        &self->geomap,
                        self->geometry.num_robjects,
                        total_num_boundary_layers),
                "Failed to malloc geometry to materials map.");

        /* set map */
        for (i = 0u; i < total_num_boundary_layers; i++) {
                mli_GeometryToMaterialMap_set(&self->geomap, 0u, i, i);
        }

        chk_msg(mli_Accelerator_malloc_from_Geometry(
                        &self->accelerator, &self->geometry),
                "Failed to malloc accelerator from geometry.");

        chk_msg(mli_Scenery_valid(self), "Expected scenery to be valid.");
        return 1;
chk_error:
        return 0;
}

/* scenery_serialize */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Scenery_to_io(const struct mli_Scenery *self, struct mli_IO *f)
{
        struct mli_MagicId magic;
        chk(mli_MagicId_set(&magic, "mli_Scenery"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        chk(mli_Geometry_to_io(&self->geometry, f));
        chk(mli_Accelerator_to_io(&self->accelerator, f));
        chk(mli_Materials_to_io(&self->materials, f));
        chk(mli_GeometryToMaterialMap_to_io(&self->geomap, f));
        return 1;
chk_error:
        return 0;
}

int mli_Scenery_from_io(struct mli_Scenery *self, struct mli_IO *f)
{
        struct mli_MagicId magic;

        mli_Scenery_free(self);

        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Scenery"));
        mli_MagicId_warn_version(&magic);

        chk(mli_Geometry_from_io(&self->geometry, f));
        chk(mli_Accelerator_from_io(&self->accelerator, f));
        chk(mli_Materials_from_io(&self->materials, f));
        chk(mli_GeometryToMaterialMap_from_io(&self->geomap, f));
        return 1;
chk_error:
        mli_Scenery_free(self);
        return 0;
}

int mli_Scenery_malloc_from_path(struct mli_Scenery *self, const char *path)
{
        struct mli_IO f = mli_IO_init();
        chk_msg(mli_IO__open_file_cstr(&f, path, "r"), "Can't open file.");
        chk_msg(mli_Scenery_from_io(self, &f), "Can't read from file.");
        mli_IO_close(&f);
        return 1;
chk_error:
        mli_IO_close(&f);
        mli_Scenery_free(self);
        return 0;
}

int mli_Scenery_write_to_path(const struct mli_Scenery *self, const char *path)
{
        struct mli_IO f = mli_IO_init();
        chk_msg(mli_IO__open_file_cstr(&f, path, "w"), "Can't write to file.");
        chk_msg(mli_Scenery_to_io(self, &f), "Can't write to file.");
        mli_IO_close(&f);
        return 1;
chk_error:
        mli_IO_close(&f);
        return 0;
}

/* scenery_tar */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Scenery_from_io_tar(struct mli_Scenery *self, struct mli_IO *f)
{
        struct mli_Archive archive = mli_Archive_init();
        chk_dbg;
        chk_msg(mli_Archive_from_io(&archive, f),
                "Can't read archive from file.");
        chk_dbg;
        chk_msg(mli_Scenery_malloc_from_Archive(self, &archive),
                "Can't malloc Scenery from Archive.");
        chk_dbg;
        mli_Archive_free(&archive);
        return 1;
chk_error:
        return 0;
}

int mli_Scenery__from_path_cstr(struct mli_Scenery *self, const char *path)
{
        struct mli_IO f = mli_IO_init();
        chk_msgf(
                mli_IO__open_file_cstr(&f, path, "r"),
                ("Can't open path '%s'.", path));
        chk_msgf(
                mli_Scenery_from_io_tar(self, &f),
                ("Can't fread Scenery from path '%s'.", path));
        mli_IO_close(&f);
        return 1;
chk_error:
        return 0;
}

int mli_Scenery_malloc_from_Archive(
        struct mli_Scenery *self,
        const struct mli_Archive *archive)
{
        uint64_t num_robjects = 0u;
        uint64_t num_objects = 0u;
        uint64_t total_num_boundary_layers = 0u;

        struct mli_materials_Names material_names = mli_materials_Names_init();
        struct mli_Map object_names = mli_Map_init();
        struct mli_Frame root = mli_Frame_init();

        chk_dbg;
        chk_msg(mli_Materials_from_Archive(
                        &self->materials, &material_names, archive),
                "Failed to malloc materials.");

        chk_dbg;
        chk(mli_Map_malloc(&object_names));

        chk_dbg;
        num_objects = mli_Archive_num_filename_prefix_sufix(
                archive, "geometry/objects/", ".obj");

        chk_dbg;
        chk_msg(mli_Geometry_malloc_objects(&self->geometry, num_objects),
                "Failed to malloc geometry.objects.");

        chk_dbg;
        chk_msg(mli_Geometry_from_archive(
                        &self->geometry, &object_names, archive),
                "Failed to malloc geometry.objects.");

        chk_dbg;
        chk_msg(mli_Frame_from_Archive(
                        &root,
                        archive,
                        &object_names,
                        self->geometry.objects,
                        &material_names.boundary_layers),
                "Failed to malloc and populate tree of frames.");

        chk_msg(mli_Frame_estimate_num_robjects_and_total_num_boundary_layers(
                        &root, &num_robjects, &total_num_boundary_layers),
                "Can not estimate num_robjects from tree of frames.");

        chk_msg(mli_GeometryToMaterialMap_malloc(
                        &self->geomap, num_robjects, total_num_boundary_layers),
                "Failed to malloc geometry to materials map.");

        chk_dbg;
        chk_msg(mli_Geometry_malloc_references(&self->geometry, num_robjects),
                "Failed to malloc geometry.references.");

        chk_msg(mli_Geometry_set_robjects_and_material_map_from_frame(
                        &root, &self->geometry, &self->geomap),
                "Can not set robjects.");

        mli_materials_Names_free(&material_names);
        mli_Map_free(&object_names);
        mli_Frame_free(&root);

        chk_msg(mli_Accelerator_malloc_from_Geometry(
                        &self->accelerator, &self->geometry),
                "Failed to malloc accelerator from geometry.");

        chk_msg(mli_Scenery_valid(self), "Expected scenery to be valid.");
        chk_msg(mli_Geometry_warn_objects(&self->geometry),
                "Failed to warn about objects.");

        chk_dbg;
        return 1;
chk_error:
        return 0;
}

/* scenery_valid */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Scenery_valid(const struct mli_Scenery *self)
{
        /* check in itself */
        chk_msg(mli_Materials_valid(&self->materials),
                "Expected materials to be valid.");
        chk_msg(mli_Geometry_valid(&self->geometry),
                "Expected geometry to be valid.");
        chk_msg(mli_Accelerator_valid(&self->accelerator),
                "Expected accelerator to be valid");
        chk_msg(mli_GeometryToMaterialMap_valid(&self->geomap),
                "Expected geometry-to-materials-map to be valid.");

        /* check interplay */
        chk_msg(mli_Accelerator_valid_wrt_Geometry(
                        &self->accelerator, &self->geometry),
                "Expected accelerator to be valid w.r.t. geometry.");
        chk_msg(mli_GeometryToMaterialMap_valid_wrt_Geometry(
                        &self->geomap, &self->geometry),
                "Expected geomap to be valid w.r.t. geometry.");
        chk_msg(mli_GeometryToMaterialMap_valid_wrt_Materials(
                        &self->geomap, &self->materials),
                "Expected geomap to be valid w.r.t. materials.");
        return 1;
chk_error:
        return 0;
}

/* shader */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Shader mli_Shader_init(void)
{
        struct mli_Shader tracer;
        tracer.scenery = NULL;
        tracer.scenery_color_materials = NULL;
        tracer.config = NULL;
        return tracer;
}

double mli_Shader_estimate_sun_visibility_weight(
        const struct mli_Shader *tracer,
        const struct mli_Vec position,
        struct mli_Prng *prng)
{
        return (1.0 - mli_Shader_estimate_sun_obstruction_weight(
                              tracer, position, prng));
}

double mli_Shader_estimate_sun_obstruction_weight(
        const struct mli_Shader *tracer,
        const struct mli_Vec position,
        struct mli_Prng *prng)
{
        uint64_t i;
        double num_obstructions = 0.0;

        for (i = 0; i < tracer->config->num_trails_global_light_source; i++) {
                struct mli_Vec pos_in_source = mli_Vec_add(
                        mli_Vec_multiply(
                                tracer->config->atmosphere.sunDirection,
                                tracer->config->atmosphere.sunDistance),
                        mli_Vec_multiply(
                                mli_Vec_random_position_inside_unit_sphere(
                                        prng),
                                tracer->config->atmosphere.sunRadius));

                struct mli_Ray line_of_sight_to_source = mli_Ray_set(
                        position, mli_Vec_substract(pos_in_source, position));

                struct mli_Intersection isec;

                const int has_intersection = mli_raytracing_query_intersection(
                        tracer->scenery, line_of_sight_to_source, &isec);

                if (has_intersection) {
                        num_obstructions += 1.0;
                }
        }

        return num_obstructions /
               tracer->config->num_trails_global_light_source;
}

struct mli_Color mli_Shader_trace_ray(
        const struct mli_Shader *tracer,
        const struct mli_Ray ray,
        struct mli_Prng *prng)
{
        struct mli_ColorSpectrum spectrum;
        struct mli_Vec xyz, rgb;

        spectrum =
                mli_Shader_trace_path_to_next_intersection(tracer, ray, prng);

        xyz = mli_ColorMaterials_ColorSpectrum_to_xyz(
                tracer->scenery_color_materials, &spectrum);

        rgb = mli_Mat_dot_product(
                &tracer->scenery_color_materials
                         ->observer_matching_curve_xyz_to_rgb,
                xyz);

        return mli_Color_set(rgb.x, rgb.y, rgb.z);
}

struct mli_ColorSpectrum mli_Shader_trace_path_to_next_intersection(
        const struct mli_Shader *tracer,
        const struct mli_Ray ray,
        struct mli_Prng *prng)
{
        struct mli_IntersectionSurfaceNormal intersection =
                mli_IntersectionSurfaceNormal_init();
        struct mli_ColorSpectrum out;
        int has_intersection =
                mli_raytracing_query_intersection_with_surface_normal(
                        tracer->scenery, ray, &intersection);

        if (has_intersection) {
                out = mli_Shader_trace_next_intersection(
                        tracer, &intersection, prng);
        } else {
                out = mli_Shader_trace_ambient_background(tracer, ray);
        }
        return out;
}

struct mli_ColorSpectrum mli_Shader_trace_ambient_background(
        const struct mli_Shader *tracer,
        const struct mli_Ray ray)
{
        if (tracer->config->have_atmosphere) {
                return mli_Shader_trace_ambient_background_atmosphere(
                        tracer, ray);
        } else {
                return mli_Shader_trace_ambient_background_whitebox(tracer);
        }
}

struct mli_ColorSpectrum mli_Shader_trace_ambient_background_atmosphere(
        const struct mli_Shader *tracer,
        const struct mli_Ray ray)
{
        return mli_Atmosphere_query(
                &tracer->config->atmosphere, ray.support, ray.direction);
}

struct mli_ColorSpectrum mli_Shader_trace_ambient_background_whitebox(
        const struct mli_Shader *tracer)
{
        return mli_ColorSpectrum_multiply_scalar(
                tracer->config->ambient_radiance_W_per_m2_per_sr, 0.5);
}

struct mli_ColorSpectrum mli_Shader_trace_ambient_sun(
        const struct mli_Shader *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        struct mli_Prng *prng)
{
        if (tracer->config->have_atmosphere) {
                return mli_Shader_trace_ambient_sun_atmosphere(
                        tracer, intersection, prng);
        } else {
                return mli_Shader_trace_ambient_sun_whitebox(
                        tracer, intersection, prng);
        }
}

struct mli_ColorSpectrum mli_Shader_trace_ambient_sun_atmosphere(
        const struct mli_Shader *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        struct mli_Prng *prng)
{
        struct mli_ColorSpectrum tone;

        const double sun_visibility = mli_Shader_estimate_sun_visibility_weight(
                tracer, intersection->position, prng);

        if (sun_visibility > 0.0) {
                tone = mli_raytracing_color_tone_of_sun(
                        tracer->config, intersection->position);
                tone = mli_ColorSpectrum_multiply_scalar(tone, sun_visibility);
        } else {
                tone = mli_raytracing_color_tone_of_diffuse_sky(
                        tracer, intersection, prng);
        }
        return tone;
}

struct mli_ColorSpectrum mli_Shader_trace_ambient_sun_whitebox(
        const struct mli_Shader *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        struct mli_Prng *prng)
{
        const double sun_visibility = mli_Shader_estimate_sun_visibility_weight(
                tracer, intersection->position, prng);

        return mli_ColorSpectrum_multiply_scalar(
                tracer->config->ambient_radiance_W_per_m2_per_sr,
                sun_visibility);
}

struct mli_ColorSpectrum mli_Shader_trace_next_intersection(
        const struct mli_Shader *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        struct mli_Prng *prng)
{
        struct mli_ColorSpectrum out;
        struct mli_IntersectionLayer intersection_layer;

        intersection_layer = mli_raytracing_get_intersection_layer(
                tracer->scenery, intersection);

        switch (intersection_layer.side_coming_from.surface->type) {
        case MLI_SURFACE_TYPE_COOKTORRANCE:
                out = mli_Shader_trace_intersection_cooktorrance(
                        tracer, intersection, &intersection_layer, prng);
                break;
        case MLI_SURFACE_TYPE_TRANSPARENT:
                out = mli_Shader_trace_intersection_transparent(
                        tracer, intersection, &intersection_layer, prng);
                break;
        default:
                chk_warning("surface type is not implemented.");
                out = mli_ColorSpectrum_init_zeros();
        }
        return out;
}

struct mli_ColorSpectrum mli_Shader_trace_intersection_cooktorrance(
        const struct mli_Shader *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        const struct mli_IntersectionLayer *intersection_layer,
        struct mli_Prng *prng)
{
        struct mli_ColorSpectrum incoming;
        struct mli_ColorSpectrum outgoing;
        struct mli_ColorSpectrum reflection;
        const struct mli_Surface_CookTorrance *cook =
                &intersection_layer->side_coming_from.surface->data
                         .cooktorrance;
        double theta;
        double lambert_factor;
        double factor;

        assert(intersection_layer->side_coming_from.surface->type ==
               MLI_SURFACE_TYPE_COOKTORRANCE);

        incoming = mli_Shader_trace_ambient_sun(tracer, intersection, prng);

        reflection = tracer->scenery_color_materials->spectra
                             .array[cook->reflection_spectrum];

        theta = mli_Vec_angle_between(
                tracer->config->atmosphere.sunDirection,
                intersection->surface_normal);

        lambert_factor = fabs(cos(theta));

        factor = lambert_factor * cook->diffuse_weight;

        outgoing = mli_ColorSpectrum_multiply(incoming, reflection);
        outgoing = mli_ColorSpectrum_multiply_scalar(outgoing, factor);
        return outgoing;
}

struct mli_ColorSpectrum mli_Shader_trace_intersection_transparent(
        const struct mli_Shader *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        const struct mli_IntersectionLayer *intersection_layer,
        struct mli_Prng *prng)
{
        const struct mli_Surface_Transparent *transparent =
                &intersection_layer->side_coming_from.surface->data.transparent;
        assert(intersection_layer->side_coming_from.surface->type ==
               MLI_SURFACE_TYPE_TRANSPARENT);

        return mli_ColorSpectrum_init_zeros();
}

/* shader_atmosphere */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_ColorSpectrum mli_raytracing_color_tone_of_sun(
        const struct mli_shader_Config *config,
        const struct mli_Vec support)
{
        struct mli_ColorSpectrum sun_spectrum = config->atmosphere.sun_spectrum;
        double width_atmosphere = config->atmosphere.atmosphereRadius -
                                  config->atmosphere.earthRadius;

        if (config->atmosphere.altitude < width_atmosphere) {
                struct mli_ColorSpectrum color_close_to_sun =
                        mli_Atmosphere_query(
                                &config->atmosphere,
                                support,
                                config->atmosphere.sunDirection);

                double f = config->atmosphere.sunDirection.z;
                uint64_t argmax = 0;
                double vmax = 1.0;
                MLI_MATH_ARRAY_ARGMAX(
                        color_close_to_sun.values,
                        MLI_COLORSPECTRUM_SIZE,
                        argmax);
                vmax = color_close_to_sun.values[argmax];
                color_close_to_sun = mli_ColorSpectrum_multiply_scalar(
                        color_close_to_sun, (1.0 / vmax));

                return mli_ColorSpectrum_add(
                        mli_ColorSpectrum_multiply_scalar(sun_spectrum, f),
                        mli_ColorSpectrum_multiply_scalar(
                                color_close_to_sun, (1.0 - f)));
        } else {
                return sun_spectrum;
        }
}

struct mli_ColorSpectrum mli_raytracing_color_tone_of_diffuse_sky(
        const struct mli_Shader *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        struct mli_Prng *prng)
{
        int i;
        struct mli_ColorSpectrum sky = mli_ColorSpectrum_init_zeros();
        struct mli_Ray obstruction_ray;
        struct mli_Vec facing_surface_normal;
        struct mli_Intersection isec;
        int has_direct_view_to_sky = 0;
        int num_samples = 5;

        facing_surface_normal =
                intersection->from_outside_to_inside
                        ? intersection->surface_normal
                        : mli_Vec_multiply(intersection->surface_normal, -1.0);

        for (i = 0; i < num_samples; i++) {
                struct mli_Vec rnd_dir = mli_Vec_random_direction_in_hemisphere(
                        prng, facing_surface_normal);

                obstruction_ray.support = intersection->position;
                obstruction_ray.direction = rnd_dir;

                has_direct_view_to_sky = !mli_raytracing_query_intersection(
                        tracer->scenery, obstruction_ray, &isec);

                if (has_direct_view_to_sky) {
                        struct mli_ColorSpectrum sample = mli_Atmosphere_query(
                                &tracer->config->atmosphere,
                                intersection->position,
                                rnd_dir);

                        double theta = mli_Vec_angle_between(
                                rnd_dir, facing_surface_normal);
                        double lambert_factor = fabs(cos(theta));

                        sky = mli_ColorSpectrum_add(
                                sky,
                                mli_ColorSpectrum_multiply_scalar(
                                        sample, lambert_factor));
                }
        }

        return mli_ColorSpectrum_multiply_scalar(
                sky, 1.0 / (255.0 * num_samples));
}

/* shader_config */
/* ------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

struct mli_shader_Config mli_shader_Config_init(void)
{
        struct mli_shader_Config config;
        mli_ColorSpectrum_set_radiance_of_black_body_W_per_m2_per_sr(
                &config.ambient_radiance_W_per_m2_per_sr, 5000.0);

        config.num_trails_global_light_source = 3;
        config.have_atmosphere = 0;
        config.atmosphere = mli_Atmosphere_init();
        return config;
}

/* shader_config_json */
/* ------------------ */

/* Copyright 2018-2021 Sebastian Achim Mueller */

int mli_shader_Config_from_json_token(
        struct mli_shader_Config *tc,
        const struct mli_Json *json,
        const uint64_t tkn)
{
        uint64_t atmtkn;
        uint64_t have_atmosphere;
        chk(mli_Json_uint64_by_key(
                json,
                tkn,
                &tc->num_trails_global_light_source,
                "num_trails_global_light_source"));
        chk_msg(tc->num_trails_global_light_source > 0,
                "Expected num_trails_global_light_source > 0.");

        chk(mli_Json_uint64_by_key(
                json, tkn, &have_atmosphere, "have_atmosphere"));
        tc->have_atmosphere = (int)have_atmosphere;

        chk(mli_Json_token_by_key(json, tkn, "atmosphere", &atmtkn));
        chk(mli_Atmosphere_from_json_token(&tc->atmosphere, json, atmtkn + 1));

        return 1;
chk_error:
        return 0;
}

/* spectrum */
/* -------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

void mli_Spectrum_free(struct mli_Spectrum *self)
{
        mli_Func_free(&self->spectrum);
        mli_FuncInfo_free(&self->info);
        mli_String_free(&self->name);
}

struct mli_Spectrum mli_Spectrum_init(void)
{
        struct mli_Spectrum out;
        out.spectrum = mli_Func_init();
        out.info = mli_FuncInfo_init();
        out.name = mli_String_init();
        return out;
}

int mli_Spectrum_equal(
        const struct mli_Spectrum *a,
        const struct mli_Spectrum *b)
{
        chk_msg(mli_Func_equal(a->spectrum, b->spectrum),
                "Different spectral function.");
        chk_msg(mli_FuncInfo_equal(&a->info, &b->info),
                "Different spectrum name.");
        chk_msg(mli_String_equal(&a->name, &b->name),
                "Different spectrum name.");
        return 1;
chk_error:
        return 0;
}


int mli_Spectrum_to_io(const struct mli_Spectrum *self, struct mli_IO *f)
{
        struct mli_MagicId magic = mli_MagicId_init();
        chk(mli_MagicId_set(&magic, "mli_Spectrum"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        chk(mli_String_to_io(&self->name, f));
        chk(mli_FuncInfo_to_io(&self->info, f));
        chk(mli_Func_to_io(&self->spectrum, f));

        return 1;
chk_error:
        return 0;
}

int mli_Spectrum_from_io(struct mli_Spectrum *self, struct mli_IO *f)
{
        struct mli_MagicId magic;
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Spectrum"));
        mli_MagicId_warn_version(&magic);

        mli_Spectrum_free(self);
        chk_msg(mli_String_from_io(&self->name, f),
                "Failed to read spectrum name.");
        chk_msg(mli_FuncInfo_from_io(&self->info, f),
                "Failed to read spectrum info.");
        chk_msg(mli_Func_from_io(&self->spectrum, f),
                "Failed to read spectrum.");
        return 1;
chk_error:
        return 0;
}

int mli_Spectrum_print_to_io(const struct mli_Spectrum *self, struct mli_IO *f)
{
        uint64_t xamin, xamax, yamin, yamax = 0;
        MLI_MATH_ARRAY_ARGMIN(
                self->spectrum.x, self->spectrum.num_points, xamin);
        MLI_MATH_ARRAY_ARGMAX(
                self->spectrum.x, self->spectrum.num_points, xamax);
        MLI_MATH_ARRAY_ARGMIN(
                self->spectrum.y, self->spectrum.num_points, yamin);
        MLI_MATH_ARRAY_ARGMAX(
                self->spectrum.y, self->spectrum.num_points, yamax);

        chk(mli_IO_text_write_cstr_format(f, "name: '%s', ", self->name.array));
        chk(mli_IO_text_write_cstr_format(
                f, "comment: '%s', ", self->info.comment.array));
        chk(mli_IO_text_write_cstr_format(
                f, "num. points: %lu\n", self->spectrum.num_points));
        chk(mli_IO_text_write_cstr_format(
                f,
                "x: [%e, %e), '%s'\n",
                self->spectrum.x[xamin],
                self->spectrum.x[xamax],
                self->info.x.array));
        chk(mli_IO_text_write_cstr_format(
                f,
                "y: [%e, %e), '%s'\n",
                self->spectrum.y[yamin],
                self->spectrum.y[yamax],
                self->info.y.array));
        return 1;
chk_error:
        return 0;
}

/* spectrum_array */
/* -------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
MLI_ARRAY_IMPLEMENTATION_FREE(
        mli_SpectrumArray,
        struct mli_Spectrum,
        mli_Spectrum_free)

/* string */
/* ------ */

/* Copyright Sebastian Achim Mueller */

MLI_VECTOR_IMPLEMENTATION_ZERO_TERMINATION(mli_String, char)

int mli_String_from_cstr_fromat(
        struct mli_String *self,
        const char *format,
        ...)
{
        va_list args;
        va_start(args, format);
        chk_msg(mli_String_from_vargs(self, format, args),
                "Failed to malloc String from variadic args.");
        va_end(args);
        return 1;
chk_error:
        va_end(args);
        return 0;
}

int mli_String_from_vargs(
        struct mli_String *self,
        const char *format,
        va_list args)
{
        struct mli_String tmp = mli_String_init();
        const int64_t tmp_capacity = 10 * strlen(format);
        chk(mli_String_malloc(&tmp, tmp_capacity));

        vsprintf(tmp.array, format, args);
        chk_msg(mli_String__discover_size(&tmp) <= tmp_capacity,
                "Probably 'vsprintf' caused a buffer overflow.");
        tmp.size = strlen(tmp.array);
        chk(mli_String_copy(self, &tmp));

        mli_String_free(&tmp);
        return 1;
chk_error:
        mli_String_free(&tmp);
        return 0;
}

int mli_String_from_cstr(struct mli_String *self, const char *s)
{
        size_t length = strlen(s);
        chk(mli_String_malloc(self, length));
        self->size = length;
        strncpy(self->array, s, self->size);
        return 1;
chk_error:
        return 0;
}

int mli_String_ends_with(
        const struct mli_String *self,
        const struct mli_String *suffix)
{
        if (!self->array || !suffix->array) {
                return 0;
        }
        if (suffix->size > self->size) {
                return 0;
        }
        return strncmp(self->array + self->size - suffix->size,
                       suffix->array,
                       suffix->size) == 0;
}

int mli_String_ends_with_cstr(const struct mli_String *self, const char *cstr)
{
        return mli_cstr_ends_with(self->array, cstr);
}

int mli_String_starts_with(
        const struct mli_String *self,
        const struct mli_String *prefix)
{
        if (!self->array || !prefix->array) {
                return 0;
        }
        if (prefix->size > self->size) {
                return 0;
        }
        return strncmp(self->array, prefix->array, prefix->size) == 0;
}

int mli_String_starts_with_cstr(const struct mli_String *self, const char *cstr)
{
        return mli_cstr_starts_with(self->array, cstr);
}

int mli_String_has_prefix_suffix(
        const struct mli_String *self,
        const struct mli_String *prefix,
        const struct mli_String *suffix)
{
        uint64_t has_pre = 1;
        uint64_t has_suf = 1;
        if (prefix->array != NULL) {
                has_pre = mli_String_starts_with(self, prefix);
        }

        if (suffix->array != NULL) {
                has_suf = mli_String_ends_with(self, suffix);
        }

        if (has_pre == 1 && has_suf == 1) {
                return 1;
        } else {
                return 0;
        }
}

int64_t mli_String_rfind(const struct mli_String *self, const char c)
{
        int64_t i;
        for (i = self->size - 1; i >= 0; i--) {
                if (self->array[i] == '\0') {
                        return -1;
                } else if (self->array[i] == c) {
                        return i;
                }
        }
        return -1;
}

int64_t mli_String_find(const struct mli_String *self, const char c)
{
        int64_t i;
        for (i = 0; i < (int64_t)self->size; i++) {
                if (self->array[i] == '\0') {
                        return -1;
                } else if (self->array[i] == c) {
                        return i;
                }
        }
        return -1;
}

int mli_String_match_templeate(
        const struct mli_String *s,
        const struct mli_String *t,
        const char digit_wildcard)
{
        uint64_t i;
        if (s->size != t->size) {
                return 0;
        }
        for (i = 0; i < s->size; i++) {
                if (t->array[i] == digit_wildcard) {
                        if (!isdigit(s->array[i])) {
                                return 0;
                        }
                } else {
                        if (s->array[i] != t->array[i]) {
                                return 0;
                        }
                }
        }
        return 1;
}

int mli_String_strip(const struct mli_String *src, struct mli_String *dst)
{
        int64_t start = 0;
        int64_t stop = 0;
        int64_t len = -1;
        struct mli_String cpysrc = mli_String_init();

        chk_msg(src->array, "Expected src-string to be allocated.");
        chk_msg(mli_String_copy(&cpysrc, src), "Can not copy input.");
        mli_String_free(dst);

        while (start < (int64_t)cpysrc.size && isspace(cpysrc.array[start])) {
                start += 1;
        }

        stop = cpysrc.size - 1;
        while (stop >= 0 && isspace(cpysrc.array[stop])) {
                stop -= 1;
        }

        len = stop - start;

        if (len < 0) {
                chk(mli_String_from_cstr_fromat(dst, ""));
        } else {
                chk(mli_String_copyn(dst, &cpysrc, start, len + 1));
        }
        mli_String_free(&cpysrc);
        return 1;
chk_error:
        mli_String_free(&cpysrc);
        mli_String_free(dst);
        return 0;
}

uint64_t mli_String_countn(
        const struct mli_String *self,
        const char match,
        const uint64_t num_chars_to_scan)
{
        uint64_t i = 0;
        uint64_t count = 0u;
        while (i < self->size && i < num_chars_to_scan) {
                if (self->array[i] == match) {
                        count++;
                }
                i++;
        }
        return count;
}

int mli_String_equal_cstr(const struct mli_String *self, const char *cstr)
{
        if (self->array == NULL) {
                return 0;
        }
        if (strcmp(self->array, cstr) == 0) {
                return 1;
        } else {
                return 0;
        }
}

int64_t mli_String_compare(
        const struct mli_String *s1,
        const struct mli_String *s2)
{
        return strcmp(s1->array, s2->array);
}

int mli_String_convert_line_break_CRLF_CR_to_LF(
        struct mli_String *dst,
        const struct mli_String *src)
{
        uint64_t i = 0;
        struct mli_String cpysrc = mli_String_init();
        chk(mli_String_malloc(&cpysrc, src->size));

        while (i < src->size) {
                if (mli_cstr_is_CRLF((char *)&src->array[i])) {
                        chk(mli_String_push_back(&cpysrc, '\n'));
                        i += 2;
                } else if (mli_cstr_is_CR((char *)&src->array[i])) {
                        chk(mli_String_push_back(&cpysrc, '\n'));
                        i += 1;
                } else {
                        chk(mli_String_push_back(&cpysrc, src->array[i]));
                        i += 1;
                }
        }
        chk(mli_String_copy(dst, &cpysrc));
        mli_String_free(&cpysrc);
        return 1;
chk_error:
        mli_String_free(&cpysrc);
        return 0;
}

int64_t mli_String__discover_size(const struct mli_String *self)
{
        int64_t i;
        for (i = 0; i < (int64_t)self->capacity; i++) {
                if (self->array[i] == '\0') {
                        break;
                }
        }

        if (i == (int64_t)self->capacity - 1) {
                if (self->array[i] != '\0') {
                        i = -1;
                }
        }
        return i;
}

int mli_String_equal(
        const struct mli_String *self,
        const struct mli_String *other)
{
        if (self->array == NULL || other->array == NULL) {
                return 0;
        }
        if (self->size == other->size) {
                size_t i;
                for (i = 0; i < self->size; i++) {
                        if (self->array[i] != other->array[i]) {
                                return 0;
                        }
                }
                return 1;
        } else {
                return 0;
        }
}

int mli_String_valid(const struct mli_String *self, const size_t min_size)
{
        const int64_t size = mli_String__discover_size(self);
        chk_msg(size >= (int64_t)min_size,
                "Expected string to have min_size "
                "and to be '\\0' terminated.");
        chk_msg(size == (int64_t)self->size,
                "Expected string.size to "
                "match zero termination.");
        return 1;
chk_error:
        return 0;
}

int mli_String__find_idx_with_cstr(
        const struct mli_String *names,
        const uint64_t num_names,
        const char *key,
        uint64_t *idx)
{
        uint64_t i;
        (*idx) = 0u;
        for (i = 0; i < num_names; i++) {
                if (mli_String_equal_cstr(&names[i], key)) {
                        (*idx) = i;
                        return 1;
                }
        }
        return 0;
}

/* string_numbers */
/* -------------- */

/* Copyright 2018-2023 Sebastian Achim Mueller */

int mli_String_nto_double(
        double *out,
        const struct mli_String *str,
        const uint64_t expected_num_chars)
{
        char *end;
        uint64_t actual_num_chars = 0u;
        double l;
        chk_msg(str->array != NULL, "Expected str to be allocated.");
        chk_msg(!(str->array[0] == '\0' || isspace(str->array[0])),
                "Can not convert string to double, bad string.");
        errno = 0;
        l = strtod(str->array, &end);
        chk_msg(errno != ERANGE,
                "Can not convert string to double, over-, under-flow.");
        chk_msg(end != NULL, "Can not convert string to double.");

        actual_num_chars = end - str->array;
        chk_msg(actual_num_chars == expected_num_chars,
                "double has not the expected number of chars.");
        *out = l;
        return 1;
chk_error:
        return 0;
}

int mli_String_to_double(double *out, const struct mli_String *str)
{
        chk_msg(mli_String_nto_double(out, str, str->size),
                "Can not convert mli_String to double.");
        return 1;
chk_error:
        return 0;
}

int mli_String_nto_int64(
        int64_t *out,
        const struct mli_String *str,
        const uint64_t base,
        const uint64_t expected_num_chars)
{
        char *end;
        uint64_t actual_num_chars = 0u;
        int64_t l;
        chk_msg(str->array != NULL, "Expected str to be allocated.");
        chk_msg(!(str->array[0] == '\0' || isspace(str->array[0])),
                "Can not convert string to int64, bad string.");
        errno = 0;
        l = strtol(str->array, &end, base);
        chk_msg(errno != ERANGE,
                "Can not convert string to int64, over-, under-flow.");
        chk_msg(end != NULL, "Can not convert string to int64, bad string.");
        actual_num_chars = end - str->array;
        chk_msg(actual_num_chars == expected_num_chars,
                "Integer has not the expected number of chars.");
        *out = l;
        return 1;
chk_error:
        return 0;
}

int mli_String_to_int64(
        int64_t *out,
        const struct mli_String *str,
        const uint64_t base)
{
        chk_msg(mli_String_nto_int64(out, str, base, str->size),
                "Can not convert string to int64.");
        return 1;
chk_error:
        return 0;
}

int mli_String_nto_uint64(
        uint64_t *out,
        const struct mli_String *str,
        const uint64_t base,
        const uint64_t expected_num_chars)
{
        int64_t tmp;
        chk(mli_String_nto_int64(&tmp, str, base, expected_num_chars));
        chk_msg(tmp >= 0, "Expected a positive integer.");
        (*out) = tmp;
        return 1;
chk_error:
        return 0;
}

int mli_String_to_uint64(
        uint64_t *out,
        const struct mli_String *str,
        const uint64_t base)
{
        int64_t tmp;
        chk(mli_String_to_int64(&tmp, str, base));
        chk_msg(tmp >= 0, "Expected a positive integer.");
        (*out) = tmp;
        return 1;
chk_error:
        return 0;
}

int mli_String_reverse_print_uint64(
        const uint64_t u,
        struct mli_String *str,
        const uint64_t base)
{
        char literals[] = {
                '0',
                '1',
                '2',
                '3',
                '4',
                '5',
                '6',
                '7',
                '8',
                '9',
                'A',
                'B',
                'C',
                'D',
                'E',
                'F'};
        uint64_t remainder = 0u;
        uint32_t remainder32 = 0u;
        uint64_t quotient = u;

        chk_msg(base <= 16, "Expected base <= 16");
        chk_msg(base > 1, "Expected base > 1");
        mli_String_malloc(str, 16);

        do {
                remainder = quotient % base;
                quotient = quotient / base;
                remainder32 = (uint32_t)remainder;
                chk(mli_String_push_back(str, literals[remainder32]));
                chk_msg(str->size < 127, "Exceeded max_num_chars.");
        } while (quotient > 0u);

        return 1;
chk_error:
        mli_String_free(str);
        return 0;
}

int mli_String_print_uint64(
        const uint64_t u,
        struct mli_String *str,
        const uint64_t base,
        const uint64_t min_num_digits,
        const char leading_char)
{
        struct mli_String tmp = mli_String_init();
        int64_t i = 0;
        int64_t length = 0;
        int64_t num_leading = 0;
        int64_t MAX_NUM_CHARS = 128;

        chk_msg(base <= 16, "Expected base <= 16");
        chk_msg(base > 1, "Expected base > 1");

        chk(mli_String_reverse_print_uint64(u, &tmp, base));

        num_leading = min_num_digits - tmp.size;
        if (num_leading < 0) {
                num_leading = 0;
        }
        length = num_leading + tmp.size;
        chk(mli_String_malloc(str, length));
        chk_msg(length < MAX_NUM_CHARS, "Exceeded max_num_chars.");

        for (i = 0; i < num_leading; i++) {
                chk(mli_String_push_back(str, leading_char));
        }

        for (i = 0; i < (int64_t)tmp.size; i++) {
                char cc;
                size_t at = (int64_t)tmp.size - i - 1;
                chk(mli_String_get(&tmp, at, &cc));
                chk(mli_String_push_back(str, cc));
        }

        mli_String_free(&tmp);
        return 1;
chk_error:
        /*mli_String_free(str);*/
        mli_String_free(&tmp);
        return 0;
}

/* string_serialize */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_String_to_io(const struct mli_String *self, struct mli_IO *f)
{
        struct mli_MagicId magic = mli_MagicId_init();
        chk(mli_MagicId_set(&magic, "mli_String"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk_IO_write(&self->size, sizeof(uint64_t), 1u, f);
        chk_IO_write(self->array, sizeof(char), self->size, f);
        return 1;
chk_error:
        return 0;
}

int mli_String_from_io(struct mli_String *self, struct mli_IO *f)
{
        uint64_t size;
        struct mli_MagicId magic;
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_String"));
        mli_MagicId_warn_version(&magic);
        chk_IO_read(&size, sizeof(uint64_t), 1u, f);
        chk(mli_String_malloc(self, size));
        self->size = size;
        chk_IO_read(self->array, sizeof(char), self->size, f);
        return 1;
chk_error:
        return 0;
}

/* string_vector */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_VECTOR_IMPLEMENTATION(mli_StringVector, struct mli_String)

/* surface */
/* ------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

struct mli_Surface mli_Surface_init(void)
{
        struct mli_Surface out;
        out.name = mli_String_init();
        out.type = MLI_SURFACE_TYPE_NONE;
        return out;
}

void mli_Surface_free(struct mli_Surface *self)
{
        mli_String_free(&self->name);
        (*self) = mli_Surface_init();
}

int mli_Surface_equal(const struct mli_Surface *a, const struct mli_Surface *b)
{
        chk_msg(a->type == b->type, "Different types of surface models.");
        chk_msg(mli_String_equal(&a->name, &b->name),
                "Different names of surface models.");

        switch (a->type) {
        case MLI_SURFACE_TYPE_TRANSPARENT:
                chk_msg(mli_Surface_Transparent_equal(
                                &a->data.transparent, &b->data.transparent),
                        "'transparent' surfaces are not equal.");
                break;
        case MLI_SURFACE_TYPE_COOKTORRANCE:
                chk_msg(mli_Surface_CookTorrance_equal(
                                &a->data.cooktorrance, &b->data.cooktorrance),
                        "'cook-torrance' surfaces are not equal.");
                break;
        default:
                chk_badf(("surface-type-id '%lu' is unknown.", a->type));
        }

        return 1;
chk_error:
        return 0;
}

int mli_Surface_type_to_string(const uint64_t type, struct mli_String *s)
{
        switch (type) {
        case MLI_SURFACE_TYPE_TRANSPARENT:
                chk(mli_String_from_cstr(s, "transparent"));
                break;
        case MLI_SURFACE_TYPE_COOKTORRANCE:
                chk(mli_String_from_cstr(s, "cook-torrance"));
                break;
        default:
                chk_badf(("surface-type-id '%lu' is unknown.", type));
        }
        return 1;
chk_error:
        return 0;
}

int mli_Surface_type_from_string(const struct mli_String *s, uint64_t *id)
{
        if (mli_String_equal_cstr(s, "transparent")) {
                (*id) = MLI_SURFACE_TYPE_TRANSPARENT;
                return 1;
        } else if (mli_String_equal_cstr(s, "cook-torrance")) {
                (*id) = MLI_SURFACE_TYPE_COOKTORRANCE;
                return 1;
        } else {
                chk_badf(("surface-type-string '%s' is unknown.", s->array));
        }
        return 1;
chk_error:
        return 0;
}

int mli_Surface_to_io(const struct mli_Surface *self, struct mli_IO *f)
{
        struct mli_MagicId magic = mli_MagicId_init();
        chk(mli_MagicId_set(&magic, "mli_Surface"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        chk_msg(mli_String_to_io(&self->name, f),
                "Can't write surface.name to io.");
        chk_IO_write(&self->type, sizeof(uint64_t), 1u, f);

        switch (self->type) {
        case MLI_SURFACE_TYPE_TRANSPARENT:
                chk_msg(mli_Surface_Transparent_to_io(
                                &self->data.transparent, f),
                        "Can't write 'transparent' surface to io.");
                break;
        case MLI_SURFACE_TYPE_COOKTORRANCE:
                chk_msg(mli_Surface_CookTorrance_to_io(
                                &self->data.cooktorrance, f),
                        "Can't write 'cook-torrance' surface to io.");
                break;
        default:
                chk_badf(("surface-type-id '%lu' is unknown.", self->type));
        }
        return 1;
chk_error:
        return 0;
}

int mli_Surface_from_io(struct mli_Surface *self, struct mli_IO *f)
{
        struct mli_MagicId magic;
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Surface"));
        mli_MagicId_warn_version(&magic);

        chk_msg(mli_String_from_io(&self->name, f),
                "Can't read surface.name from io.");
        chk_IO_read(&self->type, sizeof(uint64_t), 1u, f);

        switch (self->type) {
        case MLI_SURFACE_TYPE_TRANSPARENT:
                chk_msg(mli_Surface_Transparent_from_io(
                                &self->data.transparent, f),
                        "Can't read 'transparent' surface from io.");
                break;
        case MLI_SURFACE_TYPE_COOKTORRANCE:
                chk_msg(mli_Surface_CookTorrance_from_io(
                                &self->data.cooktorrance, f),
                        "Can't read 'cook-torrance' surface from io.");
                break;
        default:
                chk_badf(("surface-type-id '%lu' is unknown.", self->type));
        }
        return 1;
chk_error:
        return 0;
}

int mli_Surface_from_json_string_and_name(
        struct mli_Surface *self,
        const struct mli_Map *spectra_names,
        const struct mli_String *json_string,
        const struct mli_String *name)
{
        struct mli_Json json = mli_Json_init();
        struct mli_JsonWalk walk = mli_JsonWalk_init();
        struct mli_String key = mli_String_init();

        chk_msg(mli_Json_from_string(&json, json_string),
                "Can't parse surface from json string.");
        walk = mli_JsonWalk_set(&json);

        chk_msg(mli_JsonWalk_to_key(&walk, "type"),
                "Expected field 'type' in surface json string.");

        chk_msg(mli_JsonWalk_get_string(&walk, &key),
                "Expected field 'type' to hold a string.");

        chk_msgf(
                mli_Surface_type_from_string(&key, &self->type),
                ("Can't map surface 'type':'%s' from json string.", key.array));

        chk_msg(mli_String_copy(&self->name, name), "Can't copy surface name.");

        switch (self->type) {
        case MLI_SURFACE_TYPE_TRANSPARENT:
                chk_msg(mli_Surface_Transparent_from_json_string(
                                &self->data.transparent,
                                spectra_names,
                                json_string),
                        "Can't parse 'transparent' surface from json.");
                break;
        case MLI_SURFACE_TYPE_COOKTORRANCE:
                chk_msg(mli_Surface_CookTorrance_from_json_string(
                                &self->data.cooktorrance,
                                spectra_names,
                                json_string),
                        "Can't parse 'cook-torrance' surface from json.");
                break;
        default:
                chk_badf(("surface-type-id '%lu' is unknown.", self->type));
        }

        mli_String_free(&key);
        mli_Json_free(&json);
        return 1;
chk_error:
        return 0;
}

/* surface_array */
/* ------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
MLI_ARRAY_IMPLEMENTATION_FREE(
        mli_SurfaceArray,
        struct mli_Surface,
        mli_Surface_free)

/* surface_cooktorrance */
/* -------------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

int mli_Surface_CookTorrance_equal(
        const struct mli_Surface_CookTorrance *a,
        const struct mli_Surface_CookTorrance *b)
{
        if (a->reflection_spectrum != b->reflection_spectrum) {
                return 0;
        }
        if (a->diffuse_weight != b->diffuse_weight) {
                return 0;
        }
        if (a->specular_weight != b->specular_weight) {
                return 0;
        }
        if (a->roughness != b->roughness) {
                return 0;
        }
        return 1;
}

int mli_Surface_CookTorrance_to_io(
        const struct mli_Surface_CookTorrance *self,
        struct mli_IO *f)
{
        struct mli_MagicId magic = mli_MagicId_init();
        chk(mli_MagicId_set(&magic, "mli_Surface_CookTorrance"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);

        chk_IO_write(&self->reflection_spectrum, sizeof(uint64_t), 1u, f);
        chk_IO_write(&self->diffuse_weight, sizeof(double), 1u, f);
        chk_IO_write(&self->specular_weight, sizeof(double), 1u, f);
        chk_IO_write(&self->roughness, sizeof(double), 1u, f);

        return 1;
chk_error:
        return 0;
}

int mli_Surface_CookTorrance_from_io(
        struct mli_Surface_CookTorrance *self,
        struct mli_IO *f)
{
        struct mli_MagicId magic;
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Surface_CookTorrance"));
        mli_MagicId_warn_version(&magic);

        chk_IO_read(&self->reflection_spectrum, sizeof(uint64_t), 1u, f);
        chk_IO_read(&self->diffuse_weight, sizeof(double), 1u, f);
        chk_IO_read(&self->specular_weight, sizeof(double), 1u, f);
        chk_IO_read(&self->roughness, sizeof(double), 1u, f);

        return 1;
chk_error:
        return 0;
}

int mli_Surface_CookTorrance_from_json_string(
        struct mli_Surface_CookTorrance *self,
        const struct mli_Map *spectra_names,
        const struct mli_String *json_string)
{
        struct mli_Json json = mli_Json_init();
        struct mli_JsonWalk walk = mli_JsonWalk_init();
        struct mli_String key = mli_String_init();

        chk_msg(mli_Json_from_string(&json, json_string),
                "Can't parse cook-torrance surface from json string.");
        walk = mli_JsonWalk_set(&json);

        chk_msg(mli_JsonWalk_to_key(&walk, "reflection_spectrum"),
                "Expected field 'reflection_spectrum' in cook-torrance "
                "surface json string.");
        chk_msg(mli_JsonWalk_get_string(&walk, &key),
                "Expected 'reflection_spectrum' to hold a string.");
        chk_msg(mli_Map_get(spectra_names, &key, &self->reflection_spectrum),
                "Expected 'reflection_spectrum' to be in spectra_names.");

        mli_JsonWalk_to_root(&walk);
        chk_msg(mli_JsonWalk_to_key(&walk, "diffuse_weight"),
                "Expected field 'diffuse_weight' in cook-torrance.");
        chk_msg(mli_JsonWalk_get_double(&walk, &self->diffuse_weight),
                "Expected 'diffuse_weight' to hold a double.");

        mli_JsonWalk_to_root(&walk);
        chk_msg(mli_JsonWalk_to_key(&walk, "specular_weight"),
                "Expected field 'specular_weight' in cook-torrance.");
        chk_msg(mli_JsonWalk_get_double(&walk, &self->specular_weight),
                "Expected 'specular_weight' to hold a double.");

        mli_JsonWalk_to_root(&walk);
        chk_msg(mli_JsonWalk_to_key(&walk, "roughness"),
                "Expected field 'roughness' in cook-torrance.");
        chk_msg(mli_JsonWalk_get_double(&walk, &self->roughness),
                "Expected 'roughness' to hold a double.");

        mli_String_free(&key);
        mli_Json_free(&json);

        return 1;
chk_error:
        return 0;
}

/* surface_transparent */
/* ------------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

int mli_Surface_Transparent_equal(
        const struct mli_Surface_Transparent *a,
        const struct mli_Surface_Transparent *b)
{
        if (a->nothing != b->nothing) {
                return 0;
        }
        return 1;
}

int mli_Surface_Transparent_to_io(
        const struct mli_Surface_Transparent *self,
        struct mli_IO *f)
{
        struct mli_MagicId magic = mli_MagicId_init();
        chk(mli_MagicId_set(&magic, "mli_Surface_Transparent"));
        chk_IO_write(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk_IO_write(&self->nothing, sizeof(uint64_t), 1u, f);
        return 1;
chk_error:
        return 0;
}

int mli_Surface_Transparent_from_io(
        struct mli_Surface_Transparent *self,
        struct mli_IO *f)
{
        struct mli_MagicId magic;
        chk_IO_read(&magic, sizeof(struct mli_MagicId), 1u, f);
        chk(mli_MagicId_has_word(&magic, "mli_Surface_Transparent"));
        mli_MagicId_warn_version(&magic);
        chk_IO_read(&self->nothing, sizeof(uint64_t), 1u, f);
        return 1;
chk_error:
        return 0;
}

int mli_Surface_Transparent_from_json_string(
        struct mli_Surface_Transparent *self,
        const struct mli_Map *spectra_names,
        const struct mli_String *json_string)
{
        chk(spectra_names != NULL);
        chk(json_string != NULL);

        self->nothing = 0;

        return 1;
chk_error:
        return 0;
}

/* tar */
/* --- */

/**
 * Copyright (c) 2017 rxi
 * Copyright (c) 2019 Sebastian A. Mueller
 *                    Max-Planck-Institute for nuclear-physics, Heidelberg
 */


/*                             basics                                         */
/* ========================================================================== */

uint64_t mli_Tar_round_up(uint64_t n, uint64_t incr)
{
        return n + (incr - n % incr) % incr;
}

int mli_Tar_field_to_uint(
        uint64_t *out,
        const char *field,
        const uint64_t fieldsize)
{
        char buff[MLI_TAR_NAME_LENGTH] = {'\0'};
        chk(fieldsize < MLI_TAR_NAME_LENGTH);
        memcpy(buff, field, fieldsize);

        /* Take care of historic 'space' (32 decimal) termination */
        /* Convert all 'space' terminations to '\0' terminations. */

        if (buff[fieldsize - 1] == 32) {
                buff[fieldsize - 1] = 0;
        }
        if (buff[fieldsize - 2] == 32) {
                buff[fieldsize - 2] = 0;
        }

        chk(mli_cstr_to_uint64(out, buff, MLI_TAR_OCTAL));
        return 1;
chk_error:
        return 0;
}

int mli_Tar_uint_to_field(
        const uint64_t val,
        char *field,
        const uint64_t fieldsize)
{
        chk(mli_cstr_print_uint64(
                val, field, fieldsize, MLI_TAR_OCTAL, fieldsize - 1));
        return 1;
chk_error:
        return 0;
}

int mli_Tar_uint64_to_field12_2001star_base256(uint64_t val, char *field)
{
        uint8_t tmp[12];
        int64_t i = 0;
        for (i = 11; i > 0; i--) {
                tmp[i] = (uint8_t)(val % 256u);
                val = val / 256u;
        }

        chk_msg(val == 0u, "Expected value to be less than 256**11.");
        /* set highest bit in leftmost byte to 1 */
        tmp[0] = (uint8_t)128u;

        memcpy(field, tmp, 12);
        return 1;
chk_error:
        return 0;
}

int mli_Tar_field12_to_uint64_2001star_base256(const char *field, uint64_t *val)
{
        uint8_t tmp[12];
        uint64_t i = 0u;
        const uint64_t powers[] = {
                0x100000000000000,
                0x1000000000000,
                0x10000000000,
                0x100000000,
                0x1000000,
                0x10000,
                0x100,
                0x1,
        };

        memcpy(tmp, field, 12);
        chk_msg(tmp[0] == 128u,
                "Expected field[0] == 128, indicating 256-base, 2001star.");
        chk_msg(tmp[1] == 0u,
                "Expected field[1] == 0, 256**10 exceeds uint64.");
        chk_msg(tmp[2] == 0u,
                "Expected field[2] == 0, 256**09 exceeds uint64.");
        chk_msg(tmp[3] == 0u,
                "Expected field[3] == 0, 256**08 exceeds uint64.");

        (*val) = 0u;
        for (i = 4; i < 12; i++) {
                (*val) = (*val) + powers[i - 4] * (uint64_t)tmp[i];
        }
        return 1;
chk_error:
        return 0;
}

/*                               raw header                                   */
/* ========================================================================== */

uint64_t mli_TarRawHeader_checksum(const struct mli_TarRawHeader *rh)
{
        uint64_t i;
        unsigned char *p = (unsigned char *)rh;
        uint64_t res = 256;
        for (i = 0; i < offsetof(struct mli_TarRawHeader, checksum); i++) {
                res += p[i];
        }
        for (i = offsetof(struct mli_TarRawHeader, type); i < sizeof(*rh);
             i++) {
                res += p[i];
        }
        return res;
}

int mli_TarRawHeader_is_null(const struct mli_TarRawHeader *rh)
{
        uint64_t i = 0u;
        unsigned char *p = (unsigned char *)rh;
        for (i = 0; i < sizeof(struct mli_TarRawHeader); i++) {
                if (p[i] != '\0') {
                        return 0;
                }
        }
        return 1;
}

int mli_TarRawHeader_from_header(
        struct mli_TarRawHeader *rh,
        const struct mli_TarHeader *h)
{
        uint64_t chksum;

        /* Load header into raw header */
        memset(rh, 0, sizeof(*rh));
        chk_msg(mli_Tar_uint_to_field(h->mode, rh->mode, sizeof(rh->mode)),
                "bad mode");
        chk_msg(mli_Tar_uint_to_field(h->owner, rh->owner, sizeof(rh->owner)),
                "bad owner");
        if (h->size >= MLI_TAR_MAX_FILESIZE_OCTAL) {
                chk_msg(mli_Tar_uint64_to_field12_2001star_base256(
                                h->size, rh->size),
                        "bad size, mode: base-256");
        } else {
                chk_msg(mli_Tar_uint_to_field(
                                h->size, rh->size, sizeof(rh->size)),
                        "bad size, mode: base-octal");
        }
        chk_msg(mli_Tar_uint_to_field(h->mtime, rh->mtime, sizeof(rh->mtime)),
                "bad mtime");
        rh->type = h->type ? h->type : MLI_TAR_NORMAL_FILE;
        memcpy(rh->name, h->name, sizeof(rh->name));
        memcpy(rh->linkname, h->linkname, sizeof(rh->linkname));

        /* Calculate and write checksum */
        chksum = mli_TarRawHeader_checksum(rh);
        chk_msg(mli_cstr_print_uint64(
                        chksum,
                        rh->checksum,
                        sizeof(rh->checksum),
                        MLI_TAR_OCTAL,
                        sizeof(rh->checksum) - 2),
                "bad checksum");

        rh->checksum[sizeof(rh->checksum) - 1] = 32;

        chk_msg(rh->checksum[sizeof(rh->checksum) - 2] == 0,
                "Second last char in checksum must be '\\0', i.e. 0(decimal).");
        chk_msg(rh->checksum[sizeof(rh->checksum) - 1] == 32,
                "Last char in checksum must be ' ', i.e. 32(decimal).");

        return 1;
chk_error:
        return 0;
}

/*                                  header                                    */
/* ========================================================================== */

struct mli_TarHeader mli_TarHeader_init(void)
{
        struct mli_TarHeader h;
        h.mode = 0;
        h.owner = 0;
        h.size = 0;
        h.mtime = 0;
        h.type = 0;
        memset(h.name, '\0', sizeof(h.name));
        memset(h.linkname, '\0', sizeof(h.linkname));
        return h;
}

int mli_TarHeader_set_directory(struct mli_TarHeader *h, const char *name)
{
        (*h) = mli_TarHeader_init();
        chk_msg(strlen(name) < sizeof(h->name), "Dirname is too long.");
        memcpy(h->name, name, strlen(name));
        h->type = MLI_TAR_DIRECTORY;
        h->mode = 0775;
        return 1;
chk_error:
        return 0;
}

int mli_TarHeader_set_normal_file(
        struct mli_TarHeader *h,
        const char *name,
        const uint64_t size)
{
        (*h) = mli_TarHeader_init();
        chk_msg(strlen(name) < sizeof(h->name), "Filename is too long.");
        memcpy(h->name, name, strlen(name));
        h->size = size;
        h->type = MLI_TAR_NORMAL_FILE;
        h->mode = 0664;
        return 1;
chk_error:
        return 0;
}

int mli_TarHeader_from_raw(
        struct mli_TarHeader *h,
        const struct mli_TarRawHeader *rh)
{
        uint64_t chksum_actual, chksum_expected;
        chksum_actual = mli_TarRawHeader_checksum(rh);

        /* Build and compare checksum */
        chk_msg(mli_Tar_field_to_uint(
                        &chksum_expected, rh->checksum, sizeof(rh->checksum)),
                "bad checksum string.");
        chk_msg(chksum_actual == chksum_expected, "bad checksum.");

        /* Load raw header into header */
        chk_msg(mli_Tar_field_to_uint(&h->mode, rh->mode, sizeof(rh->mode)),
                "bad mode");
        chk_msg(mli_Tar_field_to_uint(&h->owner, rh->owner, sizeof(rh->owner)),
                "bad owner");
        if (rh->size[0] == -128) {
                chk_msg(mli_Tar_field12_to_uint64_2001star_base256(
                                rh->size, &h->size),
                        "bad size, mode: base-256");
        } else {
                chk_msg(mli_Tar_field_to_uint(
                                &h->size, rh->size, sizeof(rh->size)),
                        "bad size, mode: base-octal");
        }
        chk_msg(mli_Tar_field_to_uint(&h->mtime, rh->mtime, sizeof(rh->mtime)),
                "bad mtime");
        h->type = rh->type;
        memcpy(h->name, rh->name, sizeof(h->name));
        memcpy(h->linkname, rh->linkname, sizeof(h->linkname));

        return 1;
chk_error:
        return 0;
}

/* tar */
/* === */

struct mli_Tar mli_Tar_init(void)
{
        struct mli_Tar out;
        out.stream = NULL;
        out.pos = 0u;
        out.remaining_data = 0u;
        return out;
}

/*                                 read                                       */
/* ========================================================================== */

int mli_Tar_read_begin(struct mli_Tar *tar, struct mli_IO *stream)
{
        chk_msg(tar->stream == NULL,
                "Can't begin reading tar. "
                "tar is either still open or not initialized.");
        (*tar) = mli_Tar_init();
        tar->stream = stream;
        chk_msg(tar->stream, "Can't begin reading tar. Tar->stream is NULL.");
        return 1;
chk_error:
        return 0;
}

int mli_Tar_tread(struct mli_Tar *tar, void *data, const uint64_t size)
{
        int64_t res = mli_IO_read(data, 1, size, tar->stream);
        chk_msg(res >= 0, "Failed reading from tar.");
        chk_msg((uint64_t)res == size, "Failed reading from tar.");
        tar->pos += size;
        return 1;
chk_error:
        return 0;
}

int mli_Tar_read_header(struct mli_Tar *tar, struct mli_TarHeader *h)
{
        struct mli_TarRawHeader rh;

        chk_msg(mli_Tar_tread(tar, &rh, sizeof(rh)),
                "Failed to read raw header");

        if (mli_TarRawHeader_is_null(&rh)) {
                (*h) = mli_TarHeader_init();
                return 0;
        }

        chk_msg(mli_TarHeader_from_raw(h, &rh), "Failed to parse raw header.");
        tar->remaining_data = h->size;
        return 1;
chk_error:
        return 0;
}

int mli_Tar_read_data(struct mli_Tar *tar, void *ptr, uint64_t size)
{
        chk_msg(tar->remaining_data >= size,
                "Expect size to be read >= remaining_data");
        chk_msg(mli_Tar_tread(tar, ptr, size), "Failed to read payload-data.");
        tar->remaining_data -= size;

        if (tar->remaining_data == 0) {
                uint64_t i;
                const uint64_t next_record = mli_Tar_round_up(tar->pos, 512);
                const uint64_t padding_size = next_record - tar->pos;
                char padding;

                for (i = 0; i < padding_size; i++) {
                        chk_msg(mli_Tar_tread(tar, &padding, 1),
                                "Failed to read padding-block "
                                "to reach next record.");
                }
        }

        return 1;
chk_error:
        return 0;
}

int mli_Tar_read_finalize(struct mli_Tar *tar)
{
        struct mli_TarHeader h = mli_TarHeader_init();
        chk_msg(mli_Tar_read_header(tar, &h) == 0,
                "Failed to read the 2nd final block of zeros.");
        chk(h.mode == 0);
        chk(h.owner == 0);
        chk(h.size == 0);
        chk(h.mtime == 0);
        chk(h.type == 0);
        chk(h.name[0] == '\0');
        chk(h.linkname[0] == '\0');
        return 1;
chk_error:
        return 0;
}

/*                                  write                                     */
/* ========================================================================== */

int mli_Tar_write_begin(struct mli_Tar *tar, struct mli_IO *stream)
{
        chk_msg(tar->stream == NULL,
                "Can't begin writing tar. "
                "tar is either still open or not initialized.");
        (*tar) = mli_Tar_init();
        tar->stream = stream;
        chk_msg(tar->stream, "Can't begin writing tar. Tar->stream is NULL.");
        return 1;
chk_error:
        return 0;
}

int mli_Tar_twrite(struct mli_Tar *tar, const void *data, const uint64_t size)
{
        int64_t res = mli_IO_write(data, 1, size, tar->stream);
        chk_msg(res >= 0, "Failed writing to tar.");
        chk_msg((uint64_t)res == size, "Failed writing to tar.");
        tar->pos += size;
        return 1;
chk_error:
        return 0;
}

int mli_Tar_write_header(struct mli_Tar *tar, const struct mli_TarHeader *h)
{
        struct mli_TarRawHeader rh;
        chk_msg(mli_TarRawHeader_from_header(&rh, h),
                "Failed to make raw-header");
        tar->remaining_data = h->size;
        chk_msg(mli_Tar_twrite(tar, &rh, sizeof(rh)),
                "Failed to write header.");
        return 1;
chk_error:
        return 0;
}

int mli_Tar_write_null_bytes(struct mli_Tar *tar, uint64_t n)
{
        uint64_t i;
        char nul = '\0';
        for (i = 0; i < n; i++) {
                chk_msg(mli_Tar_twrite(tar, &nul, 1), "Failed to write nulls");
        }
        return 1;
chk_error:
        return 0;
}

int mli_Tar_write_data(struct mli_Tar *tar, const void *data, uint64_t size)
{
        chk_msg(tar->remaining_data >= size,
                "Expect tar->remaining_data >= size to be written.");
        chk_msg(mli_Tar_twrite(tar, data, size),
                "Failed to write payload-data.");
        tar->remaining_data -= size;

        if (tar->remaining_data == 0) {
                const uint64_t next_record = mli_Tar_round_up(tar->pos, 512);
                const uint64_t padding_size = next_record - tar->pos;
                chk_msg(mli_Tar_write_null_bytes(tar, padding_size),
                        "Failed to write padding zeros.");
        }
        return 1;
chk_error:
        return 0;
}

int mli_Tar_write_finalize(struct mli_Tar *tar)
{
        chk_msg(mli_Tar_write_null_bytes(
                        tar, sizeof(struct mli_TarRawHeader) * 2),
                "Failed to write two final null records.");
        return 1;
chk_error:
        return 0;
}

/* tar_io */
/* ------ */

/* Copyright 2018-2024 Sebastian Achim Mueller */

int mli_Tar_read_data_to_IO(
        struct mli_Tar *tar,
        struct mli_IO *buff,
        const uint64_t size)
{
        uint64_t i;
        for (i = 0; i < size; i++) {
                unsigned char c;
                chk(mli_Tar_read_data(tar, (void *)(&c), 1));
                chk(mli_IO_write((void *)(&c), sizeof(unsigned char), 1, buff));
        }

        return 1;
chk_error:
        return 0;
}

int mli_Tar_write_data_from_IO(
        struct mli_Tar *tar,
        struct mli_IO *buff,
        const uint64_t size)
{
        uint64_t i;
        chk_msg(tar->stream, "tar is not open.");
        for (i = 0; i < size; i++) {
                unsigned char c;
                chk(mli_IO_read((void *)(&c), sizeof(unsigned char), 1, buff));
                chk(mli_Tar_write_data(tar, (void *)(&c), 1));
        }

        return 1;
chk_error:
        return 0;
}

/* thin_lens */
/* --------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

double mli_thin_lens_get_object_given_focal_and_image(
        const double focal_length,
        const double image_distance)
{
        /* 1/f = 1/g + 1/b */
        /* 1/g = 1/f - 1/b */
        /* g = 1/(1/f - 1/b) */
        return 1.0 / (1.0 / focal_length - 1.0 / image_distance);
}

double mli_thin_lens_get_image_given_focal_and_object(
        const double focal_length,
        const double object_distance)
{
        /* 1/f = 1/g + 1/b */
        /* 1/b = 1/f - 1/g */
        /* b = 1/(1/f - 1/g) */
        return 1.0 / (1.0 / focal_length - 1.0 / object_distance);
}

/* toggle_stdin */
/* ------------ */


#ifdef __unix__


struct termios mli_viewer_non_canonical_stdin(void)
{
        struct termios old_terminal;
        struct termios new_terminal;
        tcgetattr(STDIN_FILENO, &old_terminal);
        new_terminal = old_terminal;
        new_terminal.c_lflag &= ~(ICANON);
        tcsetattr(STDIN_FILENO, TCSANOW, &new_terminal);
        return old_terminal;
}

void mli_viewer_restore_stdin(struct termios *old_terminal)
{
        tcsetattr(STDIN_FILENO, TCSANOW, old_terminal);
}

#endif /* __unix__ */

/* triangle */
/* -------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */

/* triangle_aabb */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

#define MLI_TRIANGLE_INSIDE 0
#define MLI_TRIANGLE_OUTSIDE 1

/* Voorhies, Douglas,
 * Triangle-Cube Intersection,
 * Graphics Gems III, p. 236-239, code: p. 521-526 */

/* Which of the six face-plane(s) is point P outside of? */

int64_t mli_triangle_aabb_face_plane(struct mli_Vec p)
{
        int64_t outcode;
        outcode = 0;
        if (p.x > .5)
                outcode |= 0x01;
        if (p.x < -.5)
                outcode |= 0x02;
        if (p.y > .5)
                outcode |= 0x04;
        if (p.y < -.5)
                outcode |= 0x08;
        if (p.z > .5)
                outcode |= 0x10;
        if (p.z < -.5)
                outcode |= 0x20;
        return (outcode);
}

/* Which of the twelve edge plane(s) is point P outside of? */

int64_t mli_triangle_aabb_bevel_2d(struct mli_Vec p)
{
        int64_t outcode;
        outcode = 0;
        if (p.x + p.y > 1.0)
                outcode |= 0x001;
        if (p.x - p.y > 1.0)
                outcode |= 0x002;
        if (-p.x + p.y > 1.0)
                outcode |= 0x004;
        if (-p.x - p.y > 1.0)
                outcode |= 0x008;
        if (p.x + p.z > 1.0)
                outcode |= 0x010;
        if (p.x - p.z > 1.0)
                outcode |= 0x020;
        if (-p.x + p.z > 1.0)
                outcode |= 0x040;
        if (-p.x - p.z > 1.0)
                outcode |= 0x080;
        if (p.y + p.z > 1.0)
                outcode |= 0x100;
        if (p.y - p.z > 1.0)
                outcode |= 0x200;
        if (-p.y + p.z > 1.0)
                outcode |= 0x400;
        if (-p.y - p.z > 1.0)
                outcode |= 0x800;
        return (outcode);
}

/* Which of the eight corner plane(s) is point P outside of? */

int64_t mli_triangle_aabb_bevel_3d(struct mli_Vec p)
{
        int64_t outcode;
        outcode = 0;
        if ((p.x + p.y + p.z) > 1.5)
                outcode |= 0x01;
        if ((p.x + p.y - p.z) > 1.5)
                outcode |= 0x02;
        if ((p.x - p.y + p.z) > 1.5)
                outcode |= 0x04;
        if ((p.x - p.y - p.z) > 1.5)
                outcode |= 0x08;
        if ((-p.x + p.y + p.z) > 1.5)
                outcode |= 0x10;
        if ((-p.x + p.y - p.z) > 1.5)
                outcode |= 0x20;
        if ((-p.x - p.y + p.z) > 1.5)
                outcode |= 0x40;
        if ((-p.x - p.y - p.z) > 1.5)
                outcode |= 0x80;
        return (outcode);
}

/* Test the point "alpha" of the way from P1 to P2 */
/* See if it is on a face of the cube              */
/* Consider only faces in "mask"                   */

int64_t mli_triangle_aabb_check_point(
        struct mli_Vec p1,
        struct mli_Vec p2,
        double alpha,
        int64_t mask)
{
        struct mli_Vec plane_point;
        plane_point.x = mli_math_linear_interpolate_1d(alpha, p1.x, p2.x);
        plane_point.y = mli_math_linear_interpolate_1d(alpha, p1.y, p2.y);
        plane_point.z = mli_math_linear_interpolate_1d(alpha, p1.z, p2.z);
        return (mli_triangle_aabb_face_plane(plane_point) & mask);
}

/* Compute intersection of P1 --> P2 line segment with face planes */
/* Then test intersection point to see if it is on cube face       */
/* Consider only face planes in "outcode_diff"                     */
/* Note: Zero bits in "outcode_diff" means face line is outside of */

int64_t mli_triangle_aabb_check_line(
        struct mli_Vec p1,
        struct mli_Vec p2,
        int64_t outcode_diff)
{
        if ((0x01 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (0.5f - p1.x) / (p2.x - p1.x), 0x3e) ==
                    MLI_TRIANGLE_INSIDE)
                        return (MLI_TRIANGLE_INSIDE);
        if ((0x02 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (-0.5f - p1.x) / (p2.x - p1.x), 0x3d) ==
                    MLI_TRIANGLE_INSIDE)
                        return (MLI_TRIANGLE_INSIDE);
        if ((0x04 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (0.5f - p1.y) / (p2.y - p1.y), 0x3b) ==
                    MLI_TRIANGLE_INSIDE)
                        return (MLI_TRIANGLE_INSIDE);
        if ((0x08 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (-0.5f - p1.y) / (p2.y - p1.y), 0x37) ==
                    MLI_TRIANGLE_INSIDE)
                        return (MLI_TRIANGLE_INSIDE);
        if ((0x10 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (0.5f - p1.z) / (p2.z - p1.z), 0x2f) ==
                    MLI_TRIANGLE_INSIDE)
                        return (MLI_TRIANGLE_INSIDE);
        if ((0x20 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (-0.5f - p1.z) / (p2.z - p1.z), 0x1f) ==
                    MLI_TRIANGLE_INSIDE)
                        return (MLI_TRIANGLE_INSIDE);
        return (MLI_TRIANGLE_OUTSIDE);
}

/* Test if 3D point is inside 3D triangle */

int64_t mli_Triangle_intersects_point(struct mli_Triangle t, struct mli_Vec p)
{
        int64_t sign12_bitmask, sign23_bitmask, sign31_bitmask;
        struct mli_Vec vect12, vect23, vect31, vect1h, vect2h, vect3h;
        struct mli_Vec cross12_1p, cross23_2p, cross31_3p;

        /* First, a quick bounding-box test:                               */
        /* If P is outside triangle bbox, there cannot be an intersection. */

        if (p.x > MLI_MATH_MAX3(t.v1.x, t.v2.x, t.v3.x))
                return (MLI_TRIANGLE_OUTSIDE);
        if (p.y > MLI_MATH_MAX3(t.v1.y, t.v2.y, t.v3.y))
                return (MLI_TRIANGLE_OUTSIDE);
        if (p.z > MLI_MATH_MAX3(t.v1.z, t.v2.z, t.v3.z))
                return (MLI_TRIANGLE_OUTSIDE);
        if (p.x < MLI_MATH_MIN3(t.v1.x, t.v2.x, t.v3.x))
                return (MLI_TRIANGLE_OUTSIDE);
        if (p.y < MLI_MATH_MIN3(t.v1.y, t.v2.y, t.v3.y))
                return (MLI_TRIANGLE_OUTSIDE);
        if (p.z < MLI_MATH_MIN3(t.v1.z, t.v2.z, t.v3.z))
                return (MLI_TRIANGLE_OUTSIDE);

        /* For each triangle side, make a vector out of it by subtracting
         * vertexes; */
        /* make another vector from one vertex to point P. */
        /* The crossproduct of these two vectors is orthogonal to both and the
         */
        /* signs of its X,Y,Z components indicate whether P was to the inside or
         */
        /* to the outside of this triangle side. */

        vect12 = mli_Vec_substract(t.v1, t.v2);
        vect1h = mli_Vec_substract(t.v1, p);
        cross12_1p = mli_Vec_cross(vect12, vect1h);
        /*sign12_bitmask = MLI_MATH_SIGN3(cross12_1p);*/
        sign12_bitmask = mli_Vec_sign3_bitmask(cross12_1p, MLI_MATH_EPSILON);
        /* Extract X,Y,Z signs as 0..7 or 0...63 integer */

        vect23 = mli_Vec_substract(t.v2, t.v3);
        vect2h = mli_Vec_substract(t.v2, p);
        cross23_2p = mli_Vec_cross(vect23, vect2h);
        sign23_bitmask = mli_Vec_sign3_bitmask(cross23_2p, MLI_MATH_EPSILON);

        vect31 = mli_Vec_substract(t.v3, t.v1);
        vect3h = mli_Vec_substract(t.v3, p);
        cross31_3p = mli_Vec_cross(vect31, vect3h);
        sign31_bitmask = mli_Vec_sign3_bitmask(cross31_3p, MLI_MATH_EPSILON);

        /* If all three crossproduct vectors agree in their component signs,  */
        /* then the point must be inside all three.                           */
        /* P cannot be MLI_TRIANGLE_OUTSIDE all three sides simultaneously. */

        /* this is the old test; with the revised MLI_MATH_SIGN3() macro, the
         * test needs to be revised. */
        return ((sign12_bitmask & sign23_bitmask & sign31_bitmask) == 0)
                       ? MLI_TRIANGLE_OUTSIDE
                       : MLI_TRIANGLE_INSIDE;
}

/**********************************************/
/* This is the main algorithm procedure.      */
/* Triangle t is compared with a unit cube,   */
/* centered on the origin.                    */
/* It returns MLI_TRIANGLE_INSIDE (0) or MLI_TRIANGLE_OUTSIDE(1) if t   */
/* intersects or does not intersect the cube. */
/**********************************************/

int64_t mli_Triangle_intersects_norm_aabb(struct mli_Triangle t)
{
        int64_t v1_test, v2_test, v3_test;
        double d, denom;
        struct mli_Vec vect12, vect13, norm;
        struct mli_Vec hitpp, hitpn, hitnp, hitnn;

        /* First compare all three vertexes with all six face-planes */
        /* If any vertex is inside the cube, return immediately!     */

        if ((v1_test = mli_triangle_aabb_face_plane(t.v1)) ==
            MLI_TRIANGLE_INSIDE)
                return (MLI_TRIANGLE_INSIDE);
        if ((v2_test = mli_triangle_aabb_face_plane(t.v2)) ==
            MLI_TRIANGLE_INSIDE)
                return (MLI_TRIANGLE_INSIDE);
        if ((v3_test = mli_triangle_aabb_face_plane(t.v3)) ==
            MLI_TRIANGLE_INSIDE)
                return (MLI_TRIANGLE_INSIDE);

        /* If all three vertexes were outside of one or more face-planes, */
        /* return immediately with a trivial rejection!                   */

        if ((v1_test & v2_test & v3_test) != 0)
                return (MLI_TRIANGLE_OUTSIDE);

        /* Now do the same trivial rejection test for the 12 edge planes */

        v1_test |= mli_triangle_aabb_bevel_2d(t.v1) << 8;
        v2_test |= mli_triangle_aabb_bevel_2d(t.v2) << 8;
        v3_test |= mli_triangle_aabb_bevel_2d(t.v3) << 8;
        if ((v1_test & v2_test & v3_test) != 0)
                return (MLI_TRIANGLE_OUTSIDE);

        /* Now do the same trivial rejection test for the 8 corner planes */

        v1_test |= mli_triangle_aabb_bevel_3d(t.v1) << 24;
        v2_test |= mli_triangle_aabb_bevel_3d(t.v2) << 24;
        v3_test |= mli_triangle_aabb_bevel_3d(t.v3) << 24;
        if ((v1_test & v2_test & v3_test) != 0)
                return (MLI_TRIANGLE_OUTSIDE);

        /* If vertex 1 and 2, as a pair, cannot be trivially rejected */
        /* by the above tests, then see if the v1-->v2 triangle edge  */
        /* intersects the cube.  Do the same for v1-->v3 and v2-->v3. */
        /* Pass to the intersection algorithm the "OR" of the outcode */
        /* bits, so that only those cube faces which are spanned by   */
        /* each triangle edge need be tested.                         */

        if ((v1_test & v2_test) == 0)
                if (mli_triangle_aabb_check_line(
                            t.v1, t.v2, v1_test | v2_test) ==
                    MLI_TRIANGLE_INSIDE)
                        return (MLI_TRIANGLE_INSIDE);
        if ((v1_test & v3_test) == 0)
                if (mli_triangle_aabb_check_line(
                            t.v1, t.v3, v1_test | v3_test) ==
                    MLI_TRIANGLE_INSIDE)
                        return (MLI_TRIANGLE_INSIDE);
        if ((v2_test & v3_test) == 0)
                if (mli_triangle_aabb_check_line(
                            t.v2, t.v3, v2_test | v3_test) ==
                    MLI_TRIANGLE_INSIDE)
                        return (MLI_TRIANGLE_INSIDE);

        /* By now, we know that the triangle is not off to any side,     */
        /* and that its sides do not penetrate the cube.  We must now    */
        /* test for the cube intersecting the interior of the triangle.  */
        /* We do this by looking for intersections between the cube      */
        /* diagonals and the triangle...first finding the intersection   */
        /* of the four diagonals with the plane of the triangle, and     */
        /* then if that intersection is inside the cube, pursuing        */
        /* whether the intersection point is inside the triangle itself. */

        /* To find plane of the triangle, first perform crossproduct on  */
        /* two triangle side vectors to compute the normal vector.       */

        vect12 = mli_Vec_substract(t.v1, t.v2);
        vect13 = mli_Vec_substract(t.v1, t.v3);
        norm = mli_Vec_cross(vect12, vect13);

        /* The normal vector "norm" X,Y,Z components are the coefficients */
        /* of the triangles AX + BY + CZ + D = 0 plane equation.  If we   */
        /* solve the plane equation for X=Y=Z (a diagonal), we get        */
        /* -D/(A+B+C) as a metric of the distance from cube center to the */
        /* diagonal/plane intersection.  If this is between -0.5 and 0.5, */
        /* the intersection is inside the cube.  If so, we continue by    */
        /* doing a point/triangle intersection.                           */
        /* Do this for all four diagonals.                                */

        d = norm.x * t.v1.x + norm.y * t.v1.y + norm.z * t.v1.z;

        /* if one of the diagonals is parallel to the plane,
         * the other will intersect the plane */
        if (fabs(denom = (norm.x + norm.y + norm.z)) > MLI_MATH_EPSILON)
        /* skip parallel diagonals to the plane; division by 0 can occur */
        {
                hitpp.x = hitpp.y = hitpp.z = d / denom;
                if (fabs(hitpp.x) <= 0.5)
                        if (mli_Triangle_intersects_point(t, hitpp) ==
                            MLI_TRIANGLE_INSIDE)
                                return (MLI_TRIANGLE_INSIDE);
        }
        if (fabs(denom = (norm.x + norm.y - norm.z)) > MLI_MATH_EPSILON) {
                hitpn.z = -(hitpn.x = hitpn.y = d / denom);
                if (fabs(hitpn.x) <= 0.5)
                        if (mli_Triangle_intersects_point(t, hitpn) ==
                            MLI_TRIANGLE_INSIDE)
                                return (MLI_TRIANGLE_INSIDE);
        }
        if (fabs(denom = (norm.x - norm.y + norm.z)) > MLI_MATH_EPSILON) {
                hitnp.y = -(hitnp.x = hitnp.z = d / denom);
                if (fabs(hitnp.x) <= 0.5)
                        if (mli_Triangle_intersects_point(t, hitnp) ==
                            MLI_TRIANGLE_INSIDE)
                                return (MLI_TRIANGLE_INSIDE);
        }
        if (fabs(denom = (norm.x - norm.y - norm.z)) > MLI_MATH_EPSILON) {
                hitnn.y = hitnn.z = -(hitnn.x = d / denom);
                if (fabs(hitnn.x) <= 0.5)
                        if (mli_Triangle_intersects_point(t, hitnn) ==
                            MLI_TRIANGLE_INSIDE)
                                return (MLI_TRIANGLE_INSIDE);
        }

        /* No edge touched the cube; no cube diagonal touched the triangle. */
        /* We're done...there was no intersection.                          */

        return (MLI_TRIANGLE_OUTSIDE);
}

struct mli_Triangle mli_Triangle_set_in_norm_aabb(
        const struct mli_Vec a,
        const struct mli_Vec b,
        const struct mli_Vec c,
        const struct mli_AABB aabb)
{
        struct mli_Triangle tri;
        struct mli_Vec aabb_center = mli_AABB_center(aabb);
        const double inv_scale_x = 1.0 / (aabb.upper.x - aabb.lower.x);
        const double inv_scale_y = 1.0 / (aabb.upper.y - aabb.lower.y);
        const double inv_scale_z = 1.0 / (aabb.upper.z - aabb.lower.z);
        /* translate */
        tri.v1 = mli_Vec_substract(a, aabb_center);
        tri.v2 = mli_Vec_substract(b, aabb_center);
        tri.v3 = mli_Vec_substract(c, aabb_center);
        /* scale */
        tri.v1.x *= inv_scale_x;
        tri.v2.x *= inv_scale_x;
        tri.v3.x *= inv_scale_x;

        tri.v1.y *= inv_scale_y;
        tri.v2.y *= inv_scale_y;
        tri.v3.y *= inv_scale_y;

        tri.v1.z *= inv_scale_z;
        tri.v2.z *= inv_scale_z;
        tri.v3.z *= inv_scale_z;
        return tri;
}

int mli_Triangle_has_overlap_aabb(
        const struct mli_Vec a,
        const struct mli_Vec b,
        const struct mli_Vec c,
        const struct mli_AABB aabb)
{
        struct mli_Triangle tri = mli_Triangle_set_in_norm_aabb(a, b, c, aabb);
        if (mli_Triangle_intersects_norm_aabb(tri) == MLI_TRIANGLE_INSIDE)
                return 1;
        else
                return 0;
}

struct mli_AABB mli_Triangle_aabb(
        const struct mli_Vec a,
        const struct mli_Vec b,
        const struct mli_Vec c)
{
        struct mli_AABB aabb;
        aabb.lower.x = MLI_MATH_MIN3(a.x, b.x, c.x);
        aabb.lower.y = MLI_MATH_MIN3(a.y, b.y, c.y);
        aabb.lower.z = MLI_MATH_MIN3(a.z, b.z, c.z);
        aabb.upper.x = MLI_MATH_MAX3(a.x, b.x, c.x);
        aabb.upper.y = MLI_MATH_MAX3(a.y, b.y, c.y);
        aabb.upper.z = MLI_MATH_MAX3(a.z, b.z, c.z);
        return aabb;
}

/* triangle_barycentric */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_triangle_BarycentrigWeights mli_triangle_barycentric_weights(
        const struct mli_Vec a,
        const struct mli_Vec b,
        const struct mli_Vec c,
        const struct mli_Vec t)
{
        /*
                                         c
                                        /|\
                                       / | \
                                      /  |  \
                                     /   |   \
                                    /    |    \
                                   /     |     \
                                  /      |      \
                                 /       |       \
                                /        |        \
                               /         |         \
                              /          |          \
                             /           |           \
                            /            |            \
                           /             |             \
                         ac              |              \
                         /               |               \
                        /                |                \
                       /      atc        |                 \
                      /                  |                  \
                     /                   t                   \
                    /                   /  \                  \
                   /                 /        \                \
                  /               /              \              \
                 /             /                    \            \
                /           at                         \          \
               /         /                                \        \
              /       /                                      \      \
             /     /                    abt                     \    \
            /   /                                                  \  \
           / /                                                        \\
          a------------------------------ab-----------------------------b

        */
        struct mli_triangle_BarycentrigWeights weights;

        const struct mli_Vec ab = mli_Vec_substract(b, a);
        const struct mli_Vec ac = mli_Vec_substract(c, a);

        const double ab_ab = mli_Vec_dot(ab, ab);
        const double ab_ac = mli_Vec_dot(ab, ac);
        const double ac_ac = mli_Vec_dot(ac, ac);

        const double abc_area_pow2 = ab_ab * ac_ac - ab_ac * ab_ac;

        const struct mli_Vec at = mli_Vec_substract(t, a);
        const double at_ab = mli_Vec_dot(at, ab);
        const double at_ac = mli_Vec_dot(at, ac);

        const double atc_area_pow2 = ab_ab * at_ac - ab_ac * at_ab;
        const double abt_area_pow2 = ac_ac * at_ab - ab_ac * at_ac;

        weights.c = atc_area_pow2 / abc_area_pow2;
        weights.b = abt_area_pow2 / abc_area_pow2;
        weights.a = 1.0 - weights.b - weights.c;

        return weights;
}

/* triangle_intersection */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Ray_intersects_triangle(
        const struct mli_Ray ray,
        const struct mli_Vec vertex_a,
        const struct mli_Vec vertex_b,
        const struct mli_Vec vertex_c,
        double *intersection_ray_parameter)
{
        /* Moeller-Trumbore-intersection-algorithm */
        struct mli_Vec edge1;
        struct mli_Vec edge2;
        struct mli_Vec h, s, q;
        double a, f, u, v, t;
        edge1 = mli_Vec_substract(vertex_b, vertex_a);
        edge2 = mli_Vec_substract(vertex_c, vertex_a);
        h = mli_Vec_cross(ray.direction, edge2);
        a = mli_Vec_dot(edge1, h);

        if (a > -MLI_MATH_EPSILON && a < MLI_MATH_EPSILON)
                return 0; /* This ray is parallel to this triangle. */
        f = 1.0 / a;
        s = mli_Vec_substract(ray.support, vertex_a);
        u = f * mli_Vec_dot(s, h);
        if (u < 0.0 || u > 1.0)
                return 0;
        q = mli_Vec_cross(s, edge1);
        v = f * mli_Vec_dot(ray.direction, q);
        if (v < 0.0 || u + v > 1.0)
                return 0;
        /* At this stage we can compute t to find out where the intersection */
        /* point is on the line. */
        t = f * mli_Vec_dot(edge2, q);
        if (t > MLI_MATH_EPSILON) {
                (*intersection_ray_parameter) = t;
                return 1;
        } else {
                /* This means that there is a line intersection but not a */
                /* ray intersection. */
                return 0;
        }
}

struct mli_Vec mli_Triangle_interpolate_surface_normal(
        const struct mli_Vec vertex_normal_a,
        const struct mli_Vec vertex_normal_b,
        const struct mli_Vec vertex_normal_c,
        const struct mli_triangle_BarycentrigWeights weights)
{
        return mli_Vec_init(
                vertex_normal_a.x * weights.a + vertex_normal_b.x * weights.b +
                        vertex_normal_c.x * weights.c,

                vertex_normal_a.y * weights.a + vertex_normal_b.y * weights.b +
                        vertex_normal_c.y * weights.c,

                vertex_normal_a.z * weights.a + vertex_normal_b.z * weights.b +
                        vertex_normal_c.z * weights.c);
}

struct mli_Vec mli_Triangle_surface_normal(
        const struct mli_Vec vertex_a,
        const struct mli_Vec vertex_b,
        const struct mli_Vec vertex_c,
        const struct mli_Vec vertex_normal_a,
        const struct mli_Vec vertex_normal_b,
        const struct mli_Vec vertex_normal_c,
        const struct mli_Vec intersection_position)
{
        struct mli_Vec surface_normal;
        struct mli_triangle_BarycentrigWeights normal_weights =
                mli_triangle_barycentric_weights(
                        vertex_a, vertex_b, vertex_c, intersection_position);

        surface_normal = mli_Triangle_interpolate_surface_normal(
                vertex_normal_a,
                vertex_normal_b,
                vertex_normal_c,
                normal_weights);

        surface_normal = mli_Vec_normalized(surface_normal);

        return surface_normal;
}

/* uint32_vector */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_VECTOR_IMPLEMENTATION(mli_Uint32Vector, uint32_t)

/* utils */
/* ----- */

/* Copyright 2020 Sebastian A. Mueller*/

float mli_chars_to_float(const char *four_char_word)
{
        float f;
        memcpy(&f, four_char_word, sizeof(float));
        return f;
}

double mli_corsika_ux_to_cx(const double ux) { return -ux; }
double mli_corsika_vy_to_cy(const double vy) { return -vy; }
double mli_corsika_wz_to_cz(const double wz) { return -wz; }

double mli_corsika_cx_to_ux(const double cx) { return -cx; }
double mli_corsika_cy_to_vy(const double cy) { return -cy; }
double mli_corsika_cz_to_wz(const double cz) { return -cz; }

double mli_corsika_restore_direction_z_component(const double x, const double y)
{
        return sqrt(1.0 - x * x - y * y);
}

/* vec */
/* --- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Vec mli_Vec_init(const double x, const double y, const double z)
{
        struct mli_Vec out;
        out.x = x;
        out.y = y;
        out.z = z;
        return out;
}

struct mli_Vec mli_Vec_add(const struct mli_Vec a, const struct mli_Vec b)
{
        struct mli_Vec out;
        out.x = a.x + b.x;
        out.y = a.y + b.y;
        out.z = a.z + b.z;
        return out;
}

struct mli_Vec mli_Vec_substract(const struct mli_Vec a, const struct mli_Vec b)
{
        struct mli_Vec out;
        out.x = a.x - b.x;
        out.y = a.y - b.y;
        out.z = a.z - b.z;
        return out;
}

struct mli_Vec mli_Vec_cross(const struct mli_Vec a, const struct mli_Vec b)
{
        struct mli_Vec out;
        out.x = (a.y * b.z - a.z * b.y);
        out.y = (a.z * b.x - a.x * b.z);
        out.z = (a.x * b.y - a.y * b.x);
        return out;
}

double mli_Vec_dot(const struct mli_Vec a, const struct mli_Vec b)
{
        return a.x * b.x + a.y * b.y + a.z * b.z;
}

struct mli_Vec mli_Vec_multiply(const struct mli_Vec v, const double a)
{
        struct mli_Vec out;
        out.x = v.x * a;
        out.y = v.y * a;
        out.z = v.z * a;
        return out;
}

double mli_Vec_norm(const struct mli_Vec a) { return sqrt(mli_Vec_dot(a, a)); }

struct mli_Vec mli_Vec_normalized(struct mli_Vec a)
{
        return mli_Vec_multiply(a, 1. / mli_Vec_norm(a));
}

double mli_Vec_angle_between(const struct mli_Vec a, const struct mli_Vec b)
{
        struct mli_Vec a_normalized = mli_Vec_multiply(a, 1. / mli_Vec_norm(a));
        struct mli_Vec b_normalized = mli_Vec_multiply(b, 1. / mli_Vec_norm(b));
        return acos(mli_Vec_dot(a_normalized, b_normalized));
}

double mli_Vec_norm_between(const struct mli_Vec a, const struct mli_Vec b)
{
        return mli_Vec_norm(mli_Vec_substract(a, b));
}

struct mli_Vec mli_Vec_mirror(
        const struct mli_Vec in,
        const struct mli_Vec normal)
{
        /*
         *      This is taken from
         *      (OPTI 421/521 – Introductory Optomechanical Engineering)
         *      J.H. Bruge
         *      University of Arizona
         *
         *                     k1    n     k2
         *                      \    /\   /
         *                       \   |   /
         *                        \  |  /
         *                         \ | /
         *      ____________________\|/______________________
         *      mirror-surface
         *
         *      k1: incidate ray
         *      k2: reflected ray
         *      n:  surface normal
         *
         *      n = [nx,ny,nz]^T
         *
         *      It can be written:
         *
         *      k2 = M*k1
         *
         *      M = EYE - 2*n*n^T
         *
         *      using EYE =  [1 0 0]
         *                   [0 1 0]
         *                   [0 0 1]
         */
        struct mli_Vec out;
        out.x = (1. - 2. * normal.x * normal.x) * in.x +
                -2. * normal.x * normal.y * in.y +
                -2. * normal.x * normal.z * in.z;

        out.y = -2. * normal.x * normal.y * in.x +
                (1. - 2. * normal.y * normal.y) * in.y +
                -2. * normal.y * normal.z * in.z;

        out.z = -2. * normal.x * normal.z * in.x +
                -2. * normal.y * normal.z * in.y +
                (1. - 2. * normal.z * normal.z) * in.z;
        return out;
}

int mli_Vec_equal_margin(
        const struct mli_Vec a,
        const struct mli_Vec b,
        const double distance_margin)
{
        struct mli_Vec diff;
        double distance_squared;
        diff = mli_Vec_substract(a, b);
        distance_squared = mli_Vec_dot(diff, diff);
        return distance_squared <= distance_margin * distance_margin;
}

int mli_Vec_equal(const struct mli_Vec a, const struct mli_Vec b)
{
        if (fabs(a.x - b.x) > DBL_EPSILON)
                return 0;
        if (fabs(a.y - b.y) > DBL_EPSILON)
                return 0;
        if (fabs(a.z - b.z) > DBL_EPSILON)
                return 0;
        return 1;
}

uint32_t mli_Vec_octant(const struct mli_Vec a)
{
        /*
         *  encodes the octant sectors where the vector is pointing to
         *      x y z sector
         *      - - -   0
         *      - - +   1
         *      - + -   2
         *      - + +   3
         *      + - -   4
         *      + - +   5
         *      + + -   6
         *      + + +   7
         */
        const uint32_t sx = a.x >= 0.;
        const uint32_t sy = a.y >= 0.;
        const uint32_t sz = a.z >= 0.;
        return 4 * sx + 2 * sy + 1 * sz;
}

int mli_Vec_sign3_bitmask(const struct mli_Vec a, const double epsilon)
{
        /* bits: 7  6  5  4  3  2  1  0  */
        /*             xp yp zp xn yn zn */

        const int xn = a.x < epsilon ? 4 : 0;   /* 2**2 */
        const int xp = a.x > -epsilon ? 32 : 0; /* 2**5 */

        const int yn = a.y < epsilon ? 2 : 0;   /* 2**1 */
        const int yp = a.y > -epsilon ? 16 : 0; /* 2**4 */

        const int zn = a.z < epsilon ? 1 : 0;  /* 2**0 */
        const int zp = a.z > -epsilon ? 8 : 0; /* 2**3 */

        return (xn | xp | yn | yp | zn | zp);
}

struct mli_Vec mli_Vec_mean(const struct mli_Vec *vecs, const uint64_t num_vecs)
{
        uint64_t i;
        struct mli_Vec mean = mli_Vec_init(0.0, 0.0, 0.0);
        for (i = 0; i < num_vecs; i++) {
                mean = mli_Vec_add(mean, vecs[i]);
        }
        return mli_Vec_multiply(mean, (1.0 / num_vecs));
}

void mli_Vec_set(struct mli_Vec *a, const uint64_t dim, const double v)
{
        switch (dim) {
        case 0:
                a->x = v;
                break;
        case 1:
                a->y = v;
                break;
        case 2:
                a->z = v;
                break;
        default:
                assert(0);
                break;
        }
}

double mli_Vec_get(const struct mli_Vec *a, const uint64_t dim)
{
        double o = 0.0;
        switch (dim) {
        case 0:
                o = a->x;
                break;
        case 1:
                o = a->y;
                break;
        case 2:
                o = a->z;
                break;
        default:
                assert(0);
                break;
        }
        return o;
}

/* vec_AABB */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Vec_overlap_aabb(
        const struct mli_Vec a,
        const struct mli_Vec aabb_lower,
        const struct mli_Vec aabb_upper)
{
        if (a.x >= aabb_lower.x && a.x <= aabb_upper.x && a.y >= aabb_lower.y &&
            a.y <= aabb_upper.y && a.z >= aabb_lower.z && a.z <= aabb_upper.z) {
                return 1;
        } else {
                return 0;
        }
}

/* vec_json */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_Vec_from_json_token(
        struct mli_Vec *v,
        const struct mli_Json *json,
        const uint64_t token)
{
        chk_msg(json->tokens[token].type == JSMN_ARRAY,
                "Expected vec-token to be a json-array.");
        chk_msg(json->tokens[token].size == 3,
                "Expected vec-token to contain exactly 3 tokens.");
        chk_msg(mli_Json_double_by_token(json, token + 1, &v->x),
                "Can not parse mli_Vec-x-value.");
        chk_msg(mli_Json_double_by_token(json, token + 2, &v->y),
                "Can not parse mli_Vec y-value.");
        chk_msg(mli_Json_double_by_token(json, token + 3, &v->z),
                "Can not parse mli_Vec z-value.");
        return 1;
chk_error:
        return 0;
}

/* vec_random */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/*
        direction
        ==========
*/
struct mli_Vec mli_Vec_random_draw_direction_in_zenith_azimuth_range(
        const struct mli_prng_ZenithRange zenith,
        const struct mli_prng_UniformRange azimuth,
        struct mli_Prng *prng)
{
        const double az = mli_prng_draw_uniform(azimuth, prng);
        const double zd = mli_prng_draw_zenith(zenith, prng);
        const double sin_zd = sin(zd);
        return mli_Vec_init(sin_zd * cos(az), sin_zd * sin(az), cos(zd));
}

struct mli_Vec mli_Vec_random_position_on_disc(
        const double radius,
        struct mli_Prng *prng)
{
        const double r = sqrt(mli_Prng_uniform(prng)) * radius;
        const double azimuth = mli_Prng_uniform(prng) * MLI_MATH_2PI;
        return mli_Vec_init(r * cos(azimuth), r * sin(azimuth), 0.0);
}

struct mli_Vec mli_Vec_random_position_inside_unit_sphere(struct mli_Prng *prng)
{
        /* rejection sampling */
        struct mli_Vec pos;
        do {
                pos.x = -1.0 + 2.0 * mli_Prng_uniform(prng);
                pos.y = -1.0 + 2.0 * mli_Prng_uniform(prng);
                pos.z = -1.0 + 2.0 * mli_Prng_uniform(prng);
        } while ((pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) > 1.0);
        return pos;
}

struct mli_Vec mli_Vec_random_direction_in_hemisphere(
        struct mli_Prng *prng,
        struct mli_Vec normal)
{
        struct mli_Vec rnd_dir;
        do {
                rnd_dir = mli_Vec_random_position_inside_unit_sphere(prng);
        } while (mli_Vec_dot(normal, rnd_dir) <= 0.0);
        return mli_Vec_normalized(rnd_dir);
}

/* vec_vector */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLI_VECTOR_IMPLEMENTATION(mli_VecVector, struct mli_Vec)

/* vector */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/* vector_dummy_testing */
/* -------------------- */

/* Copyright Sebastian Achim Mueller */

MLI_VECTOR_IMPLEMENTATION(mtlDynDummy, struct mtlDummy)
MLI_VECTOR_IMPLEMENTATION(mtlDynDummyPtr, struct mtlDummy *)

MTL_VEC_TESTING_IMPLEMENTATION(mtlDynDummy, struct mtlDummy)
MTL_VEC_TESTING_IMPLEMENTATION(mtlDynDummyPtr, struct mtlDummy *)

MLI_VECTOR_IMPLEMENTATION_ZERO_TERMINATION(mtl_VectorChar, char)

/* version */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mli_version_logo_fprint(FILE *f)
{
        fprintf(f,
                "\n  "
                "                                        _/  _/              "
                "_/\n  "
                "   _/_/_/  _/_/      _/_/    _/  _/_/  _/        _/_/_/  "
                "_/_/_/_/\n  "
                "  _/    _/    _/  _/_/_/_/  _/_/      _/  _/  _/          "
                "_/\n  "
                " _/    _/    _/  _/        _/        _/  _/  _/          _/\n "
                " "
                "_/    _/    _/    _/_/_/  _/        _/  _/    _/_/_/      "
                "_/_/\n  "
                "\n");
}

void mli_version_authors_and_affiliations_fprint(FILE *f)
{
        fprintf(f,
                "  Sebastian Achim Mueller (1,2*,3^)\n"
                "\n"
                "  [1] Max-Planck-Institute for Nuclear Physics, \n"
                "      Saupfercheckweg 1, 69117 Heidelberg, Germany\n"
                "\n"
                "  [2] Institute for Particle Physics and Astrophysics,\n"
                "      ETH-Zurich, Otto-Stern-Weg 5, 8093 Zurich, Switzerland\n"
                "\n"
                "  [3] Experimental Physics Vb, Astroparticle Physics,\n"
                "      TU-Dortmund, Otto-Hahn-Str. 4a, 44227 Dortmund, "
                "Germany\n"
                "\n"
                "   *  (2015 - 2019)\n"
                "   ^  (2013 - 2015)\n");
}

/* view */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mli_Vec mli_View_optical_axis(const struct mli_View cam)
{
        struct mli_Mat rotation = mli_Mat_init_tait_bryan(
                cam.rotation.x, cam.rotation.y, cam.rotation.z);
        return mli_transform_orientation(&rotation, mli_Vec_init(0., 0., 1.));
}

struct mli_Vec mli_View_direction_right(const struct mli_View cam)
{
        struct mli_Mat rotation = mli_Mat_init_tait_bryan(
                cam.rotation.x, cam.rotation.y, cam.rotation.z);
        return mli_transform_orientation(&rotation, mli_Vec_init(1., 0., 0.));
}

struct mli_Vec mli_View_direction_up(const struct mli_View cam)
{
        struct mli_Mat rotation = mli_Mat_init_tait_bryan(
                cam.rotation.x, cam.rotation.y, cam.rotation.z);
        return mli_transform_orientation(&rotation, mli_Vec_init(0., 1., 0.));
}

struct mli_View mli_View_move_forward(
        const struct mli_View camin,
        const double rate)
{
        struct mli_View camout = camin;
        struct mli_Vec optical_axis = mli_View_optical_axis(camin);
        camout.position = mli_Vec_add(
                camout.position, mli_Vec_multiply(optical_axis, rate));
        return camout;
}

struct mli_View mli_View_move_right(
        const struct mli_View camin,
        const double rate)
{
        struct mli_View camout = camin;
        struct mli_Vec direction_right = mli_View_direction_right(camout);
        camout.position = mli_Vec_add(
                camout.position, mli_Vec_multiply(direction_right, rate));
        return camout;
}

struct mli_View mli_View_move_up(const struct mli_View camin, const double rate)
{
        struct mli_View camout = camin;
        camout.position.z += rate;
        return camout;
}

struct mli_View mli_View_look_right(
        const struct mli_View camin,
        const double rate)
{
        struct mli_View camout = camin;
        const double diff = camin.field_of_view * rate;
        camout.rotation.z = fmod(camout.rotation.z - diff, (2. * MLI_MATH_PI));
        return camout;
}

struct mli_View mli_View_look_down_when_possible(
        const struct mli_View camin,
        const double rate)
{
        struct mli_View camout = camin;
        const double diff = camin.field_of_view * rate;
        const double next_rotation_x = camout.rotation.x + diff;
        const int fals_forward_over = next_rotation_x > MLI_MATH_PI;
        if (fals_forward_over) {
                camout.rotation.x = MLI_MATH_PI;
        } else {
                camout.rotation.x = next_rotation_x;
        }
        return camout;
}

struct mli_View mli_View_increase_fov(
        const struct mli_View camin,
        const double rate)
{
        struct mli_View camout = camin;
        if (camout.field_of_view * rate > mli_math_deg2rad(170)) {
                camout.field_of_view = mli_math_deg2rad(170);
        } else {
                camout.field_of_view *= rate;
        }
        return camout;
}

struct mli_View mli_View_decrease_fov(
        const struct mli_View camin,
        const double rate)
{
        struct mli_View camout = camin;
        if (camout.field_of_view / rate < mli_math_deg2rad(.1)) {
                camout.field_of_view = mli_math_deg2rad(.1);
        } else {
                camout.field_of_view /= rate;
        }
        return camout;
}

struct mli_View mli_View_look_up_when_possible(
        const struct mli_View camin,
        const double rate)
{
        struct mli_View camout = camin;
        const double diff = -1.0 * camin.field_of_view * rate;
        const double next_rotation_x = camout.rotation.x + diff;
        const int fals_backwards_over = next_rotation_x < 0.0;
        if (fals_backwards_over) {
                camout.rotation.x = 0.0;
        } else {
                camout.rotation.x = next_rotation_x;
        }
        return camout;
}

struct mli_HomTraComp mli_View_to_HomTraComp(const struct mli_View view)
{
        struct mli_HomTraComp view2root_comp;
        view2root_comp.translation = view.position;
        view2root_comp.rotation = mli_Quaternion_set_tait_bryan(
                view.rotation.x, view.rotation.y, view.rotation.z);
        return view2root_comp;
}

/* viewer */
/* ------ */

/* Copyright 2019 Sebastian Achim Mueller                                     */

void mli_viewer_clear_screen(void)
{
        uint64_t n = 20;
        while (n) {
                putchar('\n');
                n--;
        }
}

void mli_viewer_print_help(void)
{
        mli_viewer_clear_screen();
        mli_version_logo_fprint(stdout);
        printf("  v%d.%d.%d\n",
               MLI_VERSION_MAYOR,
               MLI_VERSION_MINOR,
               MLI_VERSION_PATCH);
        printf("\n");
        printf("  General                       Inspect\n");
        printf("    exit              [ esc ]     focus/cursor      [  c  ]\n");
        printf("    this help         [  h  ]     take picture      [space]\n");
        printf("    scenery-info      [  p  ]\n");
        printf("\n");
        printf("  Look                          Move\n");
        printf("    up                [  i  ]     forward           [  w  ]\n");
        printf("    down              [  k  ]     backward          [  s  ]\n");
        printf("    left              [  j  ]     left              [  a  ]\n");
        printf("    right             [  l  ]     right             [  d  ]\n");
        printf("                                  up                [  q  ]\n");
        printf("                                  down              [  e  ]\n");
        printf("\n");
        printf("  Field-of-view                 Quality\n");
        printf("    increase          [  m  ]     super sampling    [  b  ]\n");
        printf("    decrease          [  n  ]     color/monochrome  [  g  ]\n");
        printf("\n");
        printf("  Atmosphere                    Sun\n");
        printf("    on/off            [  0  ]     later daytime     [  9  ]\n");
        printf("    - altitude        [  4  ]     earlier daytime   [  8  ]\n");
        printf("    + altitude        [  5  ]     + latitude        [  7  ]\n");
        printf("                                  - latitude        [  6  ]\n");
        printf("  Gamma                                                    \n");
        printf("    + increase        [  x  ]                              \n");
        printf("    - decrease        [  z  ]                              \n");
        printf("\n");
        mli_version_authors_and_affiliations_fprint(stdout);
}

void mli_viewer_print_info_line(
        const struct mli_View view,
        const struct mli_viewer_Cursor cursor,
        const struct mli_shader_Config tracer_config,
        const double gamma)
{
        printf("Help 'h', "
               "Cam: "
               "pos[% -.2e, % -.2e, % -.2e]m, "
               "rot[% -.1f, % -.1f, % -.1f]deg, "
               "fov %.2fdeg, ",
               view.position.x,
               view.position.y,
               view.position.z,
               mli_math_rad2deg(view.rotation.x),
               mli_math_rad2deg(view.rotation.y),
               mli_math_rad2deg(view.rotation.z),
               mli_math_rad2deg(view.field_of_view));
        printf("gamma %.2f, ", gamma);
        printf("Sun: lat % 3.0fdeg, %02d:%02dh, alt % 3.1fkm",
               mli_math_rad2deg(tracer_config.atmosphere.sunLatitude),
               (int)(tracer_config.atmosphere.sunHourAngle),
               (int)(tracer_config.atmosphere.sunHourAngle * 60) % 60,
               tracer_config.atmosphere.altitude * 1e-3);
        if (cursor.active) {
                printf(", Cursor[%3ld, %3ld]pix", cursor.col, cursor.row);
        }
        printf(".\n");
}

void mli_viewer_timestamp_now_19chars(char *buffer)
{
        time_t now = time(0);
        struct tm *nowtm;
        nowtm = localtime(&now);
        sprintf(buffer,
                "%04u-%02u-%02u_%02u-%02u-%02u",
                nowtm->tm_year + 1900,
                nowtm->tm_mon + 1,
                nowtm->tm_mday,
                nowtm->tm_hour,
                nowtm->tm_min,
                nowtm->tm_sec);
}

int mli_viewer_get_key(void)
{
        /* Waits for keystroke and returns ascii-code.
         */
        const int key = getchar();
        if (key == EOF) {
                /* In case of EOF, return EOF */
                return key;
        } else {
                /* On some systems we have to truncate. I dont know why. */
                return key & 255;
        }
}

int mli_viewer_image_to_path(const struct mli_Image *img, const char *path)
{
        struct mli_IO f = mli_IO_init();
        chk_msg(mli_IO__open_file_cstr(&f, path, "w"),
                "Can't open path to write image.");
        chk_msg(mli_Image_to_io(img, &f), "Can't write image to file.");
        mli_IO_close(&f);
        return 1;
chk_error:
        mli_IO_close(&f);
        return 0;
}

int mli_viewer_export_image(
        const struct mli_Shader *tracer,
        const struct mli_viewer_Config config,
        const struct mli_View view,
        struct mli_Prng *prng,
        const double object_distance,
        const double gamma,
        const char *path)
{
        struct mli_Image full = mli_Image_init();
        struct mli_HomTraComp camera2root_comp;
        struct mli_camera_Aperture apcam = mli_camera_Aperture_init();

        const double image_ratio =
                ((double)config.export_num_cols /
                 (double)config.export_num_rows);
        chk_mem(mli_Image_malloc(
                &full, config.export_num_cols, config.export_num_rows));
        camera2root_comp = mli_View_to_HomTraComp(view);
        apcam.focal_length =
                mli_camera_Aperture_focal_length_given_field_of_view_and_sensor_width(
                        view.field_of_view,
                        config.aperture_camera_image_sensor_width);
        apcam.aperture_radius = 0.5 * (apcam.focal_length /
                                       config.aperture_camera_f_stop_ratio);
        apcam.image_sensor_distance =
                mli_thin_lens_get_image_given_focal_and_object(
                        apcam.focal_length, object_distance);
        apcam.image_sensor_width_x = config.aperture_camera_image_sensor_width;
        apcam.image_sensor_width_y = apcam.image_sensor_width_x / image_ratio;
        mli_camera_Aperture_render_image(
                apcam, camera2root_comp, tracer, &full, prng);
        mli_Image_power(&full, mli_Color_set(gamma, gamma, gamma));
        chk_msg(mli_viewer_image_to_path(&full, path), "Failed to write ppm.");
        mli_Image_free(&full);
        return 1;
chk_error:
        return 0;
}

int mli_viewer_run_interactive_viewer_try_non_canonical_stdin(
        const struct mli_Scenery *scenery,
        const struct mli_viewer_Config config)
{
#ifdef HAVE_TERMIOS_H
        struct termios old_terminal = mli_viewer_non_canonical_stdin();
#endif
        int rc = mli_viewer_run_interactive_viewer(scenery, config);

#ifdef HAVE_TERMIOS_H
        mli_viewer_restore_stdin(&old_terminal);
#endif
        return rc;
}

int mli_viewer_run_interactive_viewer(
        const struct mli_Scenery *scenery,
        const struct mli_viewer_Config config)
{
        struct mli_Prng prng = mli_Prng_init_MT19937(config.random_seed);
        struct mli_shader_Config tracer_config = mli_shader_Config_init();
        struct mli_Shader tracer = mli_Shader_init();
        struct mli_ColorMaterials color_materials = mli_ColorMaterials_init();
        char path[1024];
        int key;
        int super_resolution = 0;
        struct mli_viewer_Cursor cursor;
        uint64_t num_screenshots = 0;
        uint64_t print_mode = MLI_ASCII_MONOCHROME;
        char timestamp[20];
        struct mli_View view = config.view;
        struct mli_Image img = mli_Image_init();
        struct mli_Image img_gamma = mli_Image_init();
        struct mli_Image img2 = mli_Image_init();
        const double row_over_column_pixel_ratio = 2.0;
        int update_image = 1;
        int print_help = 0;
        int print_scenery_info = 0;
        int has_probing_intersection = 0;
        struct mli_IntersectionSurfaceNormal probing_intersection;
        double gamma = config.gamma;

        chk_msg(mli_ColorMaterials_malloc_from_Materials(
                        &color_materials, &scenery->materials),
                "Can't malloc color materials from scenery materials.");

        tracer.scenery = scenery;
        tracer.config = &tracer_config;
        tracer.scenery_color_materials = &color_materials;

        mli_viewer_timestamp_now_19chars(timestamp);
        chk_mem(mli_Image_malloc(
                &img, config.preview_num_cols, config.preview_num_rows));
        chk_mem(mli_Image_malloc(
                &img2,
                config.preview_num_cols * 2u,
                config.preview_num_rows * 2u));

        cursor.active = 0;
        cursor.col = config.preview_num_cols / 2;
        cursor.row = config.preview_num_rows / 2;
        cursor.num_cols = config.preview_num_cols;
        cursor.num_rows = config.preview_num_rows;
        goto show_image;

        while ((key = mli_viewer_get_key()) != MLI_VIEWER_ESCAPE_KEY) {
                update_image = 1;
                print_help = 0;
                print_scenery_info = 0;
                if (cursor.active) {
                        update_image = 0;
                        super_resolution = 0;
                        switch (key) {
                        case 'i':
                                mli_viewer_Cursor_move_up(&cursor);
                                break;
                        case 'k':
                                mli_viewer_Cursor_move_down(&cursor);
                                break;
                        case 'l':
                                mli_viewer_Cursor_move_left(&cursor);
                                break;
                        case 'j':
                                mli_viewer_Cursor_move_right(&cursor);
                                break;
                        case 'c':
                                cursor.active = !cursor.active;
                                break;
                        case 'h':
                                print_help = 1;
                                break;
                        case MLI_VIEWER_SPACE_KEY:
                                sprintf(path,
                                        "%s_%06lu.ppm",
                                        timestamp,
                                        num_screenshots);
                                num_screenshots++;
                                chk(mli_viewer_export_image(
                                        &tracer,
                                        config,
                                        view,
                                        &prng,
                                        probing_intersection.distance_of_ray,
                                        gamma,
                                        path));
                                update_image = 0;
                                break;
                        default:
                                printf("Key Press unknown: %d\n", key);
                                break;
                        }
                } else {
                        switch (key) {
                        case 'w':
                                view = mli_View_move_forward(
                                        view, config.step_length);
                                break;
                        case 's':
                                view = mli_View_move_forward(
                                        view, -config.step_length);
                                break;
                        case 'a':
                                view = mli_View_move_right(
                                        view, -config.step_length);
                                break;
                        case 'd':
                                view = mli_View_move_right(
                                        view, config.step_length);
                                break;
                        case 'q':
                                view = mli_View_move_up(
                                        view, config.step_length);
                                break;
                        case 'e':
                                view = mli_View_move_up(
                                        view, -config.step_length);
                                break;
                        case 'i':
                                view = mli_View_look_up_when_possible(
                                        view, .05);
                                break;
                        case 'k':
                                view = mli_View_look_down_when_possible(
                                        view, .05);
                                break;
                        case 'l':
                                view = mli_View_look_right(view, -.05);
                                break;
                        case 'j':
                                view = mli_View_look_right(view, .05);
                                break;
                        case 'n':
                                view = mli_View_decrease_fov(view, 1.05);
                                break;
                        case 'm':
                                view = mli_View_increase_fov(view, 1.05);
                                break;
                        case 'h':
                                print_help = 1;
                                update_image = 0;
                                break;
                        case 'b':
                                super_resolution = !super_resolution;
                                break;
                        case 'c':
                                cursor.active = !cursor.active;
                                update_image = 0;
                                break;
                        case MLI_VIEWER_SPACE_KEY:
                                printf("Go into cursor-mode first.\n");
                                break;
                        case 'g':
                                if (print_mode == MLI_ASCII_MONOCHROME) {
                                        print_mode = MLI_ANSI_ESCAPE_COLOR;
                                } else if (
                                        print_mode == MLI_ANSI_ESCAPE_COLOR) {
                                        print_mode = MLI_ASCII_MONOCHROME;
                                } else {
                                        print_mode = MLI_ASCII_MONOCHROME;
                                }
                                break;
                        case 'p':
                                print_scenery_info = 1;
                                update_image = 0;
                                break;

                        case '4':
                                mli_Atmosphere_decrease_altitude(
                                        &tracer_config.atmosphere, 0.9);
                                break;
                        case '5':
                                mli_Atmosphere_increase_altitude(
                                        &tracer_config.atmosphere, 1.1);
                                break;
                        case '6':
                                mli_Atmosphere_decrease_latitude(
                                        &tracer_config.atmosphere,
                                        mli_math_deg2rad(2.0));
                                break;
                        case '7':
                                mli_Atmosphere_increase_latitude(
                                        &tracer_config.atmosphere,
                                        mli_math_deg2rad(2.0));
                                break;
                        case '8':
                                mli_Atmosphere_decrease_hours(
                                        &tracer_config.atmosphere, 0.1);
                                break;
                        case '9':
                                mli_Atmosphere_increase_hours(
                                        &tracer_config.atmosphere, 0.1);
                                break;
                        case '0':
                                tracer_config.have_atmosphere =
                                        !tracer_config.have_atmosphere;
                                break;
                        case 'x':
                                gamma *= 1.05;
                                break;
                        case 'z':
                                gamma *= 0.95;
                                break;
                        default:
                                printf("Key Press unknown: %d\n", key);
                                update_image = 0;
                                break;
                        }
                }
        show_image:
                if (update_image) {
                        if (super_resolution) {
                                mli_camera_PinHole_render_image_with_view(
                                        view,
                                        &tracer,
                                        &img2,
                                        row_over_column_pixel_ratio,
                                        &prng);
                                mli_Image_scale_down_twice(&img2, &img);
                        } else {
                                mli_camera_PinHole_render_image_with_view(
                                        view,
                                        &tracer,
                                        &img,
                                        row_over_column_pixel_ratio,
                                        &prng);
                        }
                        chk(mli_Image_copy(&img, &img_gamma));
                        mli_Image_power(
                                &img_gamma, mli_Color_set(gamma, gamma, gamma));
                }
                mli_viewer_clear_screen();
                if (cursor.active) {
                        char symbols[1];
                        uint64_t rows[1];
                        uint64_t cols[1];
                        const uint64_t num_symbols = 1u;
                        symbols[0] = 'X';
                        rows[0] = cursor.row;
                        cols[0] = cursor.col;
                        mli_Image_print_chars(
                                &img_gamma,
                                symbols,
                                rows,
                                cols,
                                num_symbols,
                                print_mode);
                        {
                                struct mli_camera_PinHole pin_hole_camera =
                                        mli_camera_PinHole_set(
                                                view.field_of_view,
                                                &img,
                                                row_over_column_pixel_ratio);

                                struct mli_HomTraComp camera2root_comp =
                                        mli_View_to_HomTraComp(view);
                                struct mli_HomTra camera2root =
                                        mli_HomTraComp_from_compact(
                                                camera2root_comp);
                                struct mli_Ray probing_ray_wrt_camera;
                                struct mli_Ray probing_ray_wrt_root;

                                probing_ray_wrt_camera =
                                        mli_camera_PinHole_ray_at_row_col(
                                                &pin_hole_camera,
                                                &img,
                                                cursor.row,
                                                cursor.col);
                                probing_ray_wrt_root = mli_HomTraComp_ray(
                                        &camera2root, probing_ray_wrt_camera);

                                probing_intersection =
                                        mli_IntersectionSurfaceNormal_init();

                                has_probing_intersection =
                                        mli_raytracing_query_intersection_with_surface_normal(
                                                scenery,
                                                probing_ray_wrt_root,
                                                &probing_intersection);
                        }
                } else {
                        mli_Image_print(&img_gamma, print_mode);
                }
                mli_viewer_print_info_line(view, cursor, tracer_config, gamma);
                if (cursor.active) {
                        printf("Intersection: ");
                        if (has_probing_intersection) {

                                printf("(% 5d;% 5d,% 5d,% 5d)"
                                       "/"
                                       "(id;ref,obj,face), "
                                       "dist % 6.2fm, "
                                       "pos [% -.2e, % -.2e, % -.2e]m, "
                                       "normal [% -.3f, % -.3f, % -.3f], ",
                                       scenery->geometry.robject_ids
                                               [probing_intersection.geometry_id
                                                        .robj],
                                       probing_intersection.geometry_id.robj,
                                       scenery->geometry.robjects
                                               [probing_intersection.geometry_id
                                                        .robj],
                                       probing_intersection.geometry_id.face,
                                       probing_intersection.distance_of_ray,
                                       probing_intersection.position.x,
                                       probing_intersection.position.y,
                                       probing_intersection.position.z,
                                       probing_intersection.surface_normal.x,
                                       probing_intersection.surface_normal.y,
                                       probing_intersection.surface_normal.z);

                                if (probing_intersection
                                            .from_outside_to_inside) {
                                        printf("outside");
                                } else {
                                        printf(" inside");
                                }
                                printf("\n");
                        } else {
                                printf("None\n");
                        }
                } else {
                        printf("\n");
                }
                if (print_help) {
                        mli_viewer_print_help();
                }
                if (print_scenery_info) {
                        mli_viewer_clear_screen();
                        mli_Scenery_info_fprint(stdout, scenery);
                }
        }

        mli_ColorMaterials_free(&color_materials);
        mli_Image_free(&img);
        mli_Image_free(&img2);
        mli_Image_free(&img_gamma);
        return 1;
chk_error:
        mli_ColorMaterials_free(&color_materials);
        mli_Image_free(&img);
        mli_Image_free(&img2);
        mli_Image_free(&img_gamma);
        return 0;
}

