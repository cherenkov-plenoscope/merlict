#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <string.h>
#include <termios.h>

/* chk_debug */
/* --------- */

/* Copyright 2018-2021 Sebastian Achim Mueller */
#ifndef CHK_DEBUG_H_
#define CHK_DEBUG_H_


/*
 *  Based on Zed Shawn's awesome Debug Macros from his book:
 *  Learn C the hard way
 */

int chk_eprintf(const char *format, ...);

#define chk_clean_errno() (errno == 0 ? "None" : strerror(errno))

#define chk_eprint_head()                                                      \
        chk_eprintf(                                                           \
                "[ERROR] (%s:%d: errno: %s) ",                                 \
                __FILE__,                                                      \
                __LINE__,                                                      \
                chk_clean_errno())

#define chk_eprint_line(MSG)                                                   \
        {                                                                      \
                chk_eprint_head();                                             \
                chk_eprintf("%s", MSG);                                        \
                chk_eprintf("\n");                                             \
        }

#define chk_msg(C, MSG)                                                        \
        if (!(C)) {                                                            \
                chk_eprint_line(MSG);                                          \
                errno = 0;                                                     \
                goto error;                                                    \
        }

#define chk_msgf(C, MSGFMT)                                                    \
        if (!(C)) {                                                            \
                chk_eprint_head();                                             \
                chk_eprintf MSGFMT;                                            \
                chk_eprintf("\n");                                             \
                errno = 0;                                                     \
                goto error;                                                    \
        }

#define chk_bad(MSG)                                                           \
        {                                                                      \
                chk_eprint_line(MSG);                                          \
                errno = 0;                                                     \
                goto error;                                                    \
        }

#define chk(C) chk_msg(C, "Not expected.")

#define chk_mem(C) chk_msg((C), "Out of memory.")

#define chk_malloc(PTR, TYPE, NUM)                                             \
        {                                                                      \
                PTR = (TYPE *)malloc(NUM * sizeof(TYPE));                      \
                chk_mem(PTR);                                                  \
        }

#define chk_fwrite(PTR, SIZE_OF_TYPE, NUM, F)                                  \
        {                                                                      \
                const uint64_t num_written =                                   \
                        fwrite(PTR, SIZE_OF_TYPE, NUM, F);                     \
                chk_msg(num_written == NUM, "Can not write to file.");         \
        }

#define chk_fread(PTR, SIZE_OF_TYPE, NUM, F)                                   \
        {                                                                      \
                const uint64_t num_read = fread(PTR, SIZE_OF_TYPE, NUM, F);    \
                chk_msg(num_read == NUM, "Can not read from file.");           \
        }

#endif

/* mliBoundaryLayer */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIBOUNDARYLAYER_H_
#define MLIBOUNDARYLAYER_H_


struct mliSide {
        uint32_t surface;
        uint32_t medium;
};

struct mliBoundaryLayer {
        struct mliSide inner;
        struct mliSide outer;
};

int mliBoundaryLayer_equal(
        const struct mliBoundaryLayer a,
        const struct mliBoundaryLayer b);

void mliBoundaryLayer_print(const struct mliBoundaryLayer a);
#endif

/* mliColor */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLICOLOR_H_
#define MLICOLOR_H_


struct mliColor {
        float r;
        float g;
        float b;
};

int mliColor_equal(const struct mliColor a, const struct mliColor b);
struct mliColor mliColor_truncate_to_uint8(const struct mliColor color);
struct mliColor mliColor_mean(
        const struct mliColor colors[],
        const uint32_t num_colors);
struct mliColor mliColor_mix(
        const struct mliColor a,
        const struct mliColor b,
        const float refl);
struct mliColor mliColor_set(const float r, const float g, const float b);
int mliColor_is_valid_8bit_range(const struct mliColor c);

struct mliColor mliColor_add(const struct mliColor a, const struct mliColor b);
struct mliColor mliColor_multiply(const struct mliColor c, const double f);
struct mliColor mliColor_multiply_elementwise(
        const struct mliColor a,
        const struct mliColor b);
#endif

/* mliDynArray */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIDYNARRAY_H_
#define MLIDYNARRAY_H_


#define MLIDYNARRAY_DEFINITON(LIB, NAME, PAYLOAD_TYPE)                         \
                                                                               \
        struct LIB##Dyn##NAME {                                                \
                uint64_t capacity;                                             \
                uint64_t size;                                                 \
                PAYLOAD_TYPE *array;                                           \
        };                                                                     \
                                                                               \
        struct LIB##Dyn##NAME LIB##Dyn##NAME##_init(void);                     \
                                                                               \
        void LIB##Dyn##NAME##_free(struct LIB##Dyn##NAME *dh);                 \
                                                                               \
        int LIB##Dyn##NAME##_malloc(                                           \
                struct LIB##Dyn##NAME *dh, const uint64_t size);               \
                                                                               \
        int LIB##Dyn##NAME##_malloc_set_size(                                  \
                struct LIB##Dyn##NAME *dh, const uint64_t size);               \
                                                                               \
        int LIB##Dyn##NAME##_push_back(                                        \
                struct LIB##Dyn##NAME *dh, PAYLOAD_TYPE item);

#define MLIDYNARRAY_IMPLEMENTATION(LIB, NAME, PAYLOAD_TYPE)                    \
                                                                               \
        struct LIB##Dyn##NAME LIB##Dyn##NAME##_init(void)                      \
        {                                                                      \
                struct LIB##Dyn##NAME dh;                                      \
                dh.capacity = 0u;                                              \
                dh.size = 0u;                                                  \
                dh.array = NULL;                                               \
                return dh;                                                     \
        }                                                                      \
                                                                               \
        void LIB##Dyn##NAME##_free(struct LIB##Dyn##NAME *dh)                  \
        {                                                                      \
                free(dh->array);                                               \
                (*dh) = LIB##Dyn##NAME##_init();                               \
        }                                                                      \
                                                                               \
        int LIB##Dyn##NAME##_malloc(                                           \
                struct LIB##Dyn##NAME *dh, const uint64_t size)                \
        {                                                                      \
                LIB##Dyn##NAME##_free(dh);                                     \
                dh->capacity = MLI_MAX2(2, size);                              \
                dh->size = 0;                                                  \
                chk_malloc(dh->array, PAYLOAD_TYPE, dh->capacity);             \
                return 1;                                                      \
        error:                                                                 \
                return 0;                                                      \
        }                                                                      \
                                                                               \
        int LIB##Dyn##NAME##_malloc_set_size(                                  \
                struct LIB##Dyn##NAME *dh, const uint64_t size)                \
        {                                                                      \
                chk(LIB##Dyn##NAME##_malloc(dh, size));                        \
                dh->size = size;                                               \
                return 1;                                                      \
        error:                                                                 \
                return 0;                                                      \
        }                                                                      \
                                                                               \
        int LIB##Dyn##NAME##_push_back(                                        \
                struct LIB##Dyn##NAME *dh, PAYLOAD_TYPE item)                  \
        {                                                                      \
                if (dh->size == dh->capacity) {                                \
                        dh->capacity = dh->capacity * 2;                       \
                        dh->array = (PAYLOAD_TYPE *)realloc(                   \
                                (void *)dh->array,                             \
                                dh->capacity * sizeof(PAYLOAD_TYPE));          \
                        chk_mem(dh->array);                                    \
                }                                                              \
                                                                               \
                dh->array[dh->size] = item;                                    \
                dh->size += 1;                                                 \
                                                                               \
                return 1;                                                      \
        error:                                                                 \
                return 0;                                                      \
        }

#endif

/* mliDynColor */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIDYNARRAY_COLOR_H_
#define MLIDYNARRAY_COLOR_H_

MLIDYNARRAY_DEFINITON(mli, Color, struct mliColor)
MLIDYNARRAY_DEFINITON(mli, ColorPtr, struct mliColor *)
#endif

/* mliDynDouble */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIDYNARRAY_DOUBLE_H_
#define MLIDYNARRAY_DOUBLE_H_

MLIDYNARRAY_DEFINITON(mli, Double, double)
#endif

/* mliDynFloat */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIDYNARRAY_FLOAT_H_
#define MLIDYNARRAY_FLOAT_H_

MLIDYNARRAY_DEFINITON(mli, Float, float)
#endif

/* mliDynUint32 */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIDYNARRAY_UINT32_H_
#define MLIDYNARRAY_UINT32_H_

MLIDYNARRAY_DEFINITON(mli, Uint32, uint32_t)
#endif

/* mliFace */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIFACE_H_
#define MLIFACE_H_


struct mliFace {
        uint32_t a;
        uint32_t b;
        uint32_t c;
};

int mliFace_equal(const struct mliFace a, const struct mliFace b);
struct mliFace mliFace_set(
        const uint32_t a,
        const uint32_t b,
        const uint32_t c);
#endif

/* mliFunc */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIFUNC_H_
#define MLIFUNC_H_


struct mliFunc {
        uint32_t num_points;
        double *x;
        double *y;
};

int mliFunc_equal(const struct mliFunc a, const struct mliFunc b);
int mliFunc_fold_numeric(
        const struct mliFunc *a,
        const struct mliFunc *b,
        double *fold);
int mliFunc_evaluate(const struct mliFunc *f, const double xarg, double *out);
int mliFunc_x_is_strictly_increasing(const struct mliFunc *f);
int mliFunc_malloc(struct mliFunc *f, const uint32_t num_points);
void mliFunc_free(struct mliFunc *f);
struct mliFunc mliFunc_init(void);
int mliFunc_is_valid(const struct mliFunc *func);
#endif

/* mliFunc_serialize */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIFUNC_SERIALIZE_H_
#define MLIFUNC_SERIALIZE_H_


int mliFunc_malloc_fread(struct mliFunc *func, FILE *f);
int mliFunc_fwrite(const struct mliFunc *func, FILE *f);
#endif

/* mliGeometryId */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRYID_H_
#define MLIGEOMETRYID_H_


struct mliGeometryId {
        uint32_t robj;
        uint32_t face;
};

struct mliGeometryId mliGeometryId_init(void);
#endif

/* mliMagicId */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIMAGICID_H_
#define MLIMAGICID_H_


#define MLI_MAGICID_WORD_CAPACITY 52
#define MLI_MAGICID_SIZE MLI_MAGICID_WORD_CAPACITY + 12

struct mliMagicId {
        char word[MLI_MAGICID_WORD_CAPACITY];
        uint32_t mayor;
        uint32_t minor;
        uint32_t patch;
};

struct mliMagicId mliMagicId_init(void);
int mliMagicId_set(struct mliMagicId *magic, const char *word);
int mliMagicId_has_word(const struct mliMagicId *magic, const char *word);
void mliMagicId_warn_version(const struct mliMagicId *magic);
#endif

/* mliMedium */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIMEDIUM_H_
#define MLIMEDIUM_H_


struct mliMedium {
        struct mliFunc refraction;
        struct mliFunc absorbtion;
};
struct mliMedium mliMedium_init(void);
void mliMedium_free(struct mliMedium *medium);

int mliMedium_equal(const struct mliMedium *a, const struct mliMedium *b);

#endif

/* mliMedium_serialize */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIMEDIUM_SERIALIZE_H_
#define MLIMEDIUM_SERIALIZE_H_


int mliMedium_fwrite(const struct mliMedium *med, FILE *f);
int mliMedium_malloc_fread(struct mliMedium *med, FILE *f);
#endif

/* mliName */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLINAME_H_
#define MLINAME_H_


#define MLI_NAME_CAPACITY 128

struct mliName {
        char cstr[MLI_NAME_CAPACITY];
};

struct mliName mliName_init(void);
int mliName_valid(const struct mliName *name);
int mliName_equal(const struct mliName *a, const struct mliName *b);
int mliName_find_idx(
        const struct mliName *names,
        const uint64_t num_names,
        const char *key,
        uint64_t *idx);
#endif

/* mliOctOverlaps */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIOCTOVERLAPS_H_
#define MLIOCTOVERLAPS_H_


#define mliOctOverlap mliDynUint32
#define mliOctOverlap_init mliDynUint32_init
#define mliOctOverlap_malloc mliDynUint32_malloc
#define mliOctOverlap_free mliDynUint32_free
#define mliOctOverlap_push_back mliDynUint32_push_back

#endif

/* mliPixels */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIPIXELS_H_
#define MLIPIXELS_H_


struct mliPixel {
        uint16_t row;
        uint16_t col;
};

struct mliPixels {
        uint32_t num_pixels_to_do;
        uint32_t num_pixels;
        struct mliPixel *pixels;
};

int mliPixels_malloc(struct mliPixels *pix, const uint32_t num_pixels);
void mliPixels_free(struct mliPixels *pix);
struct mliPixels mliPixels_init(void);
#endif

/* mliStr */
/* ------ */

/* Copyright 2018-2023 Sebastian Achim Mueller */
#ifndef mliStr_H_
#define mliStr_H_

struct mliStr {
        uint64_t length;
        char *cstr;
};

struct mliStr mliStr_init(void);
int mliStr_malloc_copy(struct mliStr *str, const struct mliStr *src);
int mliStr_malloc_copyn(
        struct mliStr *str,
        const struct mliStr *src,
        const uint64_t start,
        const uint64_t length);
int mliStr_malloc(struct mliStr *str, const uint64_t length);
void mliStr_free(struct mliStr *str);
int mliStr_mallocf(struct mliStr *str, const char *format, ...);
int mliStr_malloc_cstr(struct mliStr *str, const char *s);

int mliStr_ends_with(const struct mliStr *str, const struct mliStr *suffix);
int mliStr_starts_with(const struct mliStr *str, const struct mliStr *prefix);
int mliStr_has_prefix_suffix(
        const struct mliStr *str,
        const struct mliStr *prefix,
        const struct mliStr *suffix);

int64_t mliStr_rfind(const struct mliStr *str, const char c);
int64_t mliStr_find(const struct mliStr *str, const char c);
int mliStr_strip(const struct mliStr *src, struct mliStr *dst);
uint64_t mliStr_countn(
        const struct mliStr *str,
        const char c,
        const uint64_t num_chars_to_scan);

#endif

/* mliStr_numbers */
/* -------------- */

/* Copyright 2018-2023 Sebastian Achim Mueller */
#ifndef mliStr_numbers_H_
#define mliStr_numbers_H_

int mliStr_nto_double(
        double *out,
        const struct mliStr *str,
        const uint64_t expected_num_chars);
int mliStr_to_double(double *out, const struct mliStr *str);
int mliStr_nto_int64(
        int64_t *out,
        const struct mliStr *str,
        const uint64_t base,
        const uint64_t expected_num_chars);
int mliStr_to_int64(
        int64_t *out,
        const struct mliStr *str,
        const uint64_t base);
int mliStr_nto_uint64(
        uint64_t *out,
        const struct mliStr *str,
        const uint64_t base,
        const uint64_t expected_num_chars);
int mliStr_to_uint64(
        uint64_t *out,
        const struct mliStr *str,
        const uint64_t base);

int mliStr_print_uint64(
        const uint64_t u,
        struct mliStr *str,
        const uint64_t base,
        const uint64_t min_num_digits,
        const char leading_char);

#endif

/* mliTar */
/* ------ */

/**
 * Copyright (c) 2017 rxi
 * Copyright (c) 2019 Sebastian A. Mueller
 *                    Max-Planck-Institute for nuclear-physics, Heidelberg
 *
 * This library is free software; you can redistribute it and/or modify it
 * under the terms of the MIT license.
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to
 * deal in the Software without restriction, including without limitation the
 * rights to use, copy, modify, merge, publish, distribute, sublicense, and/or
 * sell copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS
 * IN THE SOFTWARE.
 */

#ifndef MLITAR_H_
#define MLITAR_H_


#define MLI_TAR_VERSION_MAYOR 1
#define MLI_TAR_VERSION_MINOR 0
#define MLI_TAR_VERSION_PATCH 0

#define MLI_TAR_NORMAL_FILE '0'
#define MLI_TAR_HARD_LINK '1'
#define MLI_TAR_SYMBOLIC_LINK '2'
#define MLI_TAR_CHARACTER_SPECIAL '3'
#define MLI_TAR_BLOCK_SPECIAL '4'
#define MLI_TAR_DIRECTORY '5'
#define MLI_TAR_FIFO '6'
#define MLI_TAR_NAME_LENGTH 100

#define MLI_TAR_OCTAL 8u
#define MLI_TAR_MAX_FILESIZE_OCTAL 8589934592lu /* 8^11 */

/* basics */
/* ====== */
uint64_t mliTar_round_up(uint64_t n, uint64_t incr);
int mliTar_field_to_uint(
        uint64_t *out,
        const char *field,
        const uint64_t fieldsize);
int mliTar_uint_to_field(
        const uint64_t value,
        char *field,
        const uint64_t fieldsize);
int mliTar_uint64_to_field12_2001star_base256(uint64_t val, char *field);
int mliTar_field12_to_uint64_2001star_base256(const char *field, uint64_t *val);

/* header and raw header */
/* ===================== */
struct mliTarRawHeader {
        char name[MLI_TAR_NAME_LENGTH];
        char mode[8];
        char owner[8];
        char group[8];
        char size[12];
        char mtime[12];
        char checksum[8];
        char type;
        char linkname[MLI_TAR_NAME_LENGTH];
        char _padding[255];
};

struct mliTarHeader {
        uint64_t mode;
        uint64_t owner;
        uint64_t size;
        uint64_t mtime;
        uint64_t type;
        char name[MLI_TAR_NAME_LENGTH];
        char linkname[MLI_TAR_NAME_LENGTH];
};

uint64_t mliTarRawHeader_checksum(const struct mliTarRawHeader *rh);
int mliTarRawHeader_is_null(const struct mliTarRawHeader *rh);
int mliTarRawHeader_from_header(
        struct mliTarRawHeader *rh,
        const struct mliTarHeader *h);

struct mliTarHeader mliTarHeader_init(void);
int mliTarHeader_set_directory(struct mliTarHeader *h, const char *name);
int mliTarHeader_set_normal_file(
        struct mliTarHeader *h,
        const char *name,
        const uint64_t size);
int mliTarHeader_from_raw(
        struct mliTarHeader *h,
        const struct mliTarRawHeader *rh);

/* tar */
/* === */
struct mliTar {
        FILE *stream;
        uint64_t pos;
        uint64_t remaining_data;
};

struct mliTar mliTar_init(void);

int mliTar_read_begin(struct mliTar *tar, FILE *file);
int mliTar_read_header(struct mliTar *tar, struct mliTarHeader *h);
int mliTar_read_data(struct mliTar *tar, void *ptr, uint64_t size);
int mliTar_read_finalize(struct mliTar *tar);

int mliTar_write_begin(struct mliTar *tar, FILE *file);
int mliTar_write_header(struct mliTar *tar, const struct mliTarHeader *h);
int mliTar_write_data(struct mliTar *tar, const void *data, uint64_t size);
int mliTar_write_finalize(struct mliTar *tar);

#endif

/* mliVec */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIVEC_H_
#define MLIVEC_H_


struct mliVec {
        double x;
        double y;
        double z;
};

void mliVec_print(const struct mliVec v);
uint32_t mliVec_octant(const struct mliVec a);
int mliVec_equal(const struct mliVec a, const struct mliVec b);
int mliVec_equal_margin(
        const struct mliVec a,
        const struct mliVec b,
        const double distance_margin);
struct mliVec mliVec_mirror(const struct mliVec in, const struct mliVec normal);
double mliVec_norm_between(const struct mliVec a, const struct mliVec b);
double mliVec_angle_between(const struct mliVec a, const struct mliVec b);
struct mliVec mliVec_normalized(struct mliVec a);
double mliVec_norm(const struct mliVec a);
struct mliVec mliVec_multiply(const struct mliVec v, const double a);
double mliVec_dot(const struct mliVec a, const struct mliVec b);
struct mliVec mliVec_cross(const struct mliVec a, const struct mliVec b);
struct mliVec mliVec_substract(const struct mliVec a, const struct mliVec b);
struct mliVec mliVec_add(const struct mliVec a, const struct mliVec b);
struct mliVec mliVec_init(const double x, const double y, const double z);
int mliVec_sign3_bitmask(const struct mliVec a, const double epsilon);
struct mliVec mliVec_mean(const struct mliVec *vecs, const uint64_t num_vecs);
void mliVec_set(struct mliVec *a, const uint64_t dim, const double val);
double mliVec_get(const struct mliVec *a, const uint64_t dim);
#endif

/* mliVec_AABB */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIVEC_AABB_H_
#define MLIVEC_AABB_H_


int mliVec_overlap_aabb(
        const struct mliVec a,
        const struct mliVec aabb_lower,
        const struct mliVec aabb_upper);
#endif

/* mliView */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIVIEW_H_
#define MLIVIEW_H_


struct mliView {
        struct mliVec position;
        struct mliVec rotation;
        double field_of_view;
};

struct mliView mliView_look_up_when_possible(
        const struct mliView camin,
        const double rate);
struct mliView mliView_decrease_fov(
        const struct mliView camin,
        const double rate);
struct mliView mliView_increase_fov(
        const struct mliView camin,
        const double rate);
struct mliView mliView_look_down_when_possible(
        const struct mliView camin,
        const double rate);
struct mliView mliView_look_right(
        const struct mliView camin,
        const double rate);
struct mliView mliView_move_up(const struct mliView camin, const double rate);
struct mliView mliView_move_right(
        const struct mliView camin,
        const double rate);
struct mliView mliView_move_forward(
        const struct mliView camin,
        const double rate);
struct mliVec mliView_direction_up(const struct mliView cam);
struct mliVec mliView_direction_right(const struct mliView cam);
struct mliVec mliView_optical_axis(const struct mliView cam);
struct mliHomTraComp mliView_to_HomTraComp(const struct mliView view);
#endif

/* mli_barycentric */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_BARYCENTRIC_H_
#define MLI_BARYCENTRIC_H_


struct mliBarycentrigWeights {
        double a;
        double b;
        double c;
};

struct mliBarycentrigWeights mli_barycentric_weights(
        const struct mliVec a,
        const struct mliVec b,
        const struct mliVec c,
        const struct mliVec t);
#endif

/* mli_cstr */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_CSTR_H_
#define MLI_CSTR_H_


int mli_cstr_ends_with(const char *str, const char *sufix);
int mli_cstr_starts_with(const char *str, const char *prefix);
int mli_cstr_has_prefix_suffix(
        const char *str,
        const char *prefix,
        const char *sufix);

int mli_cstr_split(
        const char *str,
        const char delimiter,
        char *token,
        const uint64_t token_length);
int mli_cstr_is_CRLF(const char *s);
int mli_cstr_is_CR(const char *s);
int mli_cstr_assert_only_NUL_LF_TAB_controls(const char *str);
int mli_cstr_assert_only_NUL_LF_TAB_controls_dbg(
        const char *str,
        const int dbg);

uint64_t mli_cstr_count_chars_up_to(
        const char *str,
        const char c,
        const uint64_t num_chars_to_scan);

int mli_cstr_lines_fprint(
        FILE *f,
        const char *str,
        const uint64_t line,
        const uint64_t line_radius);
void mli_cstr_path_strip_this_dir(char *dst, const char *src);

void mli_cstr_path_basename_without_extension(const char *filename, char *key);
void mli_cstr_strip_spaces(const char *in, char *out);

int mli_cstr_match_templeate(
        const char *s,
        const char *t,
        const char digit_wildcard);

#endif

/* mli_cstr_numbers */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_CSTR_NUMBERS_H_
#define MLI_CSTR_NUMBERS_H_


int mli_cstr_nto_int64(
        int64_t *out,
        const char *s,
        const uint64_t base,
        const uint64_t length);
int mli_cstr_to_int64(int64_t *out, const char *s, const uint64_t base);

int mli_cstr_nto_uint64(
        uint64_t *out,
        const char *s,
        const uint64_t base,
        const uint64_t length);
int mli_cstr_to_uint64(uint64_t *out, const char *s, const uint64_t base);

int mli_cstr_nto_double(double *out, const char *s, const uint64_t length);
int mli_cstr_to_double(double *out, const char *s);

int mli_cstr_print_uint64(
        uint64_t u,
        char *s,
        const uint64_t max_num_chars,
        const uint64_t base,
        const uint64_t min_num_digits);

#endif

/* mli_frame_to_scenery */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FRAME_TO_SCENERY_H_
#define MLI_FRAME_TO_SCENERY_H_

struct mliFrame;
struct mliGeometry;
struct mliGeometryToMaterialMap;

int mliFrame_set_robjects_and_material_map(
        const struct mliFrame *frame,
        struct mliGeometry *geometry,
        struct mliGeometryToMaterialMap *geomap);

int mliFrame_set_robjects_and_material_map_walk(
        const struct mliFrame *frame,
        struct mliGeometry *geometry,
        struct mliGeometryToMaterialMap *geomap,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers);

int mliFrame_estimate_num_robjects_and_total_num_boundary_layers(
        const struct mliFrame *frame,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers);

int mliFrame_estimate_num_robjects_and_total_num_boundary_layers_walk(
        const struct mliFrame *frame,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers);

#endif

/* mli_from_outside_to_inside */
/* -------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FROM_OUTSIDE_TO_INSIDE_H_
#define MLI_FROM_OUTSIDE_TO_INSIDE_H_


int mli_ray_runs_from_outside_to_inside(
        const struct mliVec ray_direction_local,
        const struct mliVec surface_normal_local);
#endif

/* mli_json_jsmn */
/* ------------- */

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
 */

#ifndef JSMN_H_
#define JSMN_H_


/**
 * JSON type identifier. Basic types are:
 *      o Object
 *      o Array
 *      o String
 *      o Other primitive: number, boolean (true/false) or null
 */
enum jsmntype_t {
        JSMN_UNDEFINED = 0,
        JSMN_OBJECT = 1,
        JSMN_ARRAY = 2,
        JSMN_STRING = 3,
        JSMN_PRIMITIVE = 4
};

enum jsmnerr {
        /* Not enough tokens were provided */
        JSMN_ERROR_NOMEM = -1,
        /* Invalid character inside JSON string */
        JSMN_ERROR_INVAL = -2,
        /* The string is not a full JSON packet, more bytes expected */
        JSMN_ERROR_PART = -3
};

/**
 * JSON token description.
 * type         type (object, array, string etc.)
 * start        start position in JSON data string
 * end          end position in JSON data string
 */
struct jsmntok_t {
        enum jsmntype_t type;
        int start;
        int end;
        int size;
};

/**
 * JSON parser. Contains an array of token blocks available. Also stores
 * the string being parsed now and current position in that string.
 */
struct jsmn_parser {
        unsigned int pos;     /* offset in the JSON string */
        unsigned int toknext; /* next token to allocate */
        int toksuper; /* superior token node, e.g. parent object or array */
};

/**
 * Create JSON parser over an array of tokens
 */
int jsmn_parse(
        struct jsmn_parser *parser,
        const char *js,
        const size_t len,
        struct jsmntok_t *tokens,
        const unsigned int num_tokens);
void jsmn_init(struct jsmn_parser *parser);

#endif

/* mli_math */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MATH_H_
#define MLI_MATH_H_


#define MLI_PI 3.14159265358979323846
#define MLI_2PI 6.28318530717958623199
#define MLI_2_OVER_SQRT3 1.1547005383792517
#define MLI_SQRT3_OVER_2 0.8660254037844386
#define MLI_EPSILON 1e-9
#define MLI_NAN 0. / 0.
#define MLI_IS_NAN(a) ((a) != (a))
#define MLI_MIN2(a, b) (((a) < (b)) ? (a) : (b))
#define MLI_MAX2(a, b) (((a) > (b)) ? (a) : (b))
#define MLI_ROUND(num) (num - floor(num) > 0.5) ? ceil(num) : floor(num)
#define MLI_NEAR_INT(x) ((x) > 0 ? (int64_t)((x) + 0.5) : (int64_t)((x)-0.5))

#define MLI_MIN3(a, b, c)                                                      \
        ((((a) < (b)) && ((a) < (c))) ? (a) : (((b) < (c)) ? (b) : (c)))

#define MLI_MAX3(a, b, c)                                                      \
        ((((a) > (b)) && ((a) > (c))) ? (a) : (((b) > (c)) ? (b) : (c)))

#define MLI_ARRAY_SET(arr, val, num)                                           \
        do {                                                                   \
                uint64_t i;                                                    \
                for (i = 0; i < num; i++) {                                    \
                        arr[i] = val;                                          \
                }                                                              \
        } while (0)

#define MLI_UPPER_COMPARE(points, num_points, point_arg, return_idx)           \
        do {                                                                   \
                uint64_t first, last, middle;                                  \
                first = 0u;                                                    \
                last = num_points - 1u;                                        \
                middle = (last - first) / 2;                                   \
                if (num_points == 0) {                                         \
                        return_idx = 0;                                        \
                } else {                                                       \
                        if (point_arg >= points[num_points - 1u]) {            \
                                return_idx = num_points;                       \
                        } else {                                               \
                                while (first < last) {                         \
                                        if (points[middle] > point_arg) {      \
                                                last = middle;                 \
                                        } else {                               \
                                                first = middle + 1u;           \
                                        }                                      \
                                        middle = first + (last - first) / 2;   \
                                }                                              \
                                return_idx = last;                             \
                        }                                                      \
                }                                                              \
        } while (0)

#define MLI_NCPY(src, dst, num)                                                \
        do {                                                                   \
                uint64_t i;                                                    \
                for (i = 0; i < num; i++) {                                    \
                        dst[i] = src[i];                                       \
                }                                                              \
        } while (0)

double mli_std(
        const double vals[],
        const uint64_t size,
        const double vals_mean);
double mli_mean(const double vals[], const uint64_t size);
void mli_linspace(
        const double start,
        const double stop,
        double *points,
        const uint64_t num_points);
void mli_histogram(
        const double *bin_edges,
        const uint64_t num_bin_edges,
        uint64_t *underflow_bin,
        uint64_t *bins,
        uint64_t *overflow_bin,
        const double point);
uint64_t mli_upper_compare_double(
        const double *points,
        const uint64_t num_points,
        const double point_arg);
double mli_square(const double a);
double mli_hypot(const double a, const double b);
double mli_deg2rad(const double angle_in_deg);
double mli_rad2deg(const double angle_in_rad);
double mli_bin_center_in_linear_space(
        const double start,
        const double stop,
        const uint64_t num_bins,
        const uint64_t bin);
double mli_linear_interpolate_1d(
        const double weight,
        const double start,
        const double end);
double mli_linear_interpolate_2d(
        const double xarg,
        const double x0,
        const double y0,
        const double x1,
        const double y1);
double mli_relative_ratio(const double a, const double b);
#endif

/* mli_quadratic_equation */
/* ---------------------- */

/* Copyright 2019 Sebastian A. Mueller */
#ifndef MLI_QUADRATIC_EQUATION_H_
#define MLI_QUADRATIC_EQUATION_H_

int mli_quadratic_equation(
        const double p,
        const double q,
        double *minus_solution,
        double *plus_solution);

#endif

/* mli_random_MT19937 */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RANDOM_MT19937_H_
#define MLI_RANDOM_MT19937_H_


struct mliMT19937 {
        uint32_t N;
        uint32_t M;
        int R;
        int A;
        int F;
        int U;
        int S;
        int B;
        int T;
        int C;
        int L;
        int MASK_LOWER;
        int MASK_UPPER;
        uint32_t mt[624];
        uint16_t index;
};

uint32_t mliMT19937_generate_uint32(struct mliMT19937 *mt);
void mliMT19937_twist(struct mliMT19937 *mt);
struct mliMT19937 mliMT19937_init(const uint32_t seed);
void mliMT19937_reinit(struct mliMT19937 *mt, const uint32_t seed);
void mliMT19937_set_constants(struct mliMT19937 *mt);
#endif

/* mli_random_pcg_variants_32bit_subset */
/* ------------------------------------ */

/**
 *  2021 March 23, Sebastian Achim Mueller
 *
 *  Based on 'pcg_variants.h' written by Melissa O'Neill.
 *
 *  I only kept the version with the 64bit sequence state to generate
 *  32bit numbers.
 *  I dropped 'advance', and 'boundedrand'.
 *  I only kept the seeding and the generation.
 *  I split the original header-only into a source.c and a header.h.
 */

/*
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

#ifndef MLI_RANDOM_PCG_VARIANTS_32BIT_SUBSET_H_INCLUDED
#define MLI_RANDOM_PCG_VARIANTS_32BIT_SUBSET_H_INCLUDED


struct pcg_state_setseq_64 {
        uint64_t state;
        uint64_t inc;
};

void pcg_setseq_64_srandom_r(
        struct pcg_state_setseq_64 *rng,
        uint64_t initstate,
        uint64_t initseq);

uint32_t pcg_setseq_64_xsh_rr_32_random_r(struct pcg_state_setseq_64 *rng);

#endif

/* mli_version */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VERSION_H_
#define MLI_VERSION_H_


#define MLI_VERSION_MAYOR 1
#define MLI_VERSION_MINOR 8
#define MLI_VERSION_PATCH 1

void mli_logo_fprint(FILE *f);
void mli_authors_and_affiliations_fprint(FILE *f);
#endif

/* mli_viewer_Config */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VIEWER_CONFIG_H_
#define MLI_VIEWER_CONFIG_H_


struct mlivrConfig {
        uint32_t random_seed;
        uint64_t preview_num_cols;
        uint64_t preview_num_rows;
        uint64_t export_num_cols;
        uint64_t export_num_rows;
        double step_length;
        struct mliView view;

        double aperture_camera_f_stop_ratio;
        double aperture_camera_image_sensor_width;
};

struct mlivrConfig mlivrConfig_default(void);

#endif

/* mli_viewer_Cursor */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VIEWER_CURSOR_H_
#define MLI_VIEWER_CURSOR_H_


struct mlivrCursor {
        int64_t active;
        uint64_t col;
        uint64_t row;
        uint64_t num_cols;
        uint64_t num_rows;
};

void mlivrCursor_move_up(struct mlivrCursor *cursor);
void mlivrCursor_move_down(struct mlivrCursor *cursor);
void mlivrCursor_move_right(struct mlivrCursor *cursor);
void mlivrCursor_move_left(struct mlivrCursor *cursor);

#endif

/* mli_viewer_toggle_stdin */
/* ----------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VIEWER_TOGGLE_STDIN_H_
#define MLI_VIEWER_TOGGLE_STDIN_H_

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#define HAVE_TERMIOS_H 1
#endif

#ifdef HAVE_TERMIOS_H
struct termios mlivr_non_canonical_stdin(void);
void mlivr_restore_stdin(struct termios *old_terminal);
#endif

#endif

/* mliAABB */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIAABB_H_
#define MLIAABB_H_


struct mliAABB {
        /*
         * Rectangular (A)xis-(A)ligned-(B)ounding-(B)ox
         * oriented w.r.t. the unit-vectors.
         *
         *                     O----------------------O
         *                    /.                     /|
         *                   / .                    / |
         *                  /  .                   /  |
         *                 /   .                  /   |
         *                O----------------------O upper
         *                |    .                 |    |
         *      Z         |    .                 |    |
         *      |       lo|wer O- - - - - - - - -| - -O
         *      |         |   .                  |   /
         *      |         |  .                   |  /
         *      /-----Y   | .                    | /
         *     /          |.                     |/
         *    X           O----------------------O
         *
         *
         */
        struct mliVec lower;
        struct mliVec upper;
};

struct mliAABB mliAABB_set(
        const struct mliVec lower,
        const struct mliVec upper);
struct mliVec mliAABB_center(const struct mliAABB a);
struct mliAABB mliAABB_outermost(
        const struct mliAABB a,
        const struct mliAABB b);
int mliAABB_valid(const struct mliAABB a);
int mliAABB_equal(const struct mliAABB a, const struct mliAABB b);
int mliAABB_overlapp_aabb(const struct mliAABB a, const struct mliAABB b);
#endif

/* mliAtmosphere */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIATMOSPHERE_H_
#define MLIATMOSPHERE_H_


struct mliAtmosphere {
        double sunLatitude;
        double sunHourAngle;

        struct mliVec sunDirection;
        double sunDistance;
        double sunRadius;

        double earthRadius;
        double atmosphereRadius;

        double Height_Rayleigh;
        double Height_Mie;

        struct mliColor beta_Rayleigh;
        struct mliColor beta_Mie;

        uint64_t numSamples;
        uint64_t numSamplesLight;

        double power;
        double altitude;
};

struct mliAtmosphere mliAtmosphere_init(void);
void mliAtmosphere_set_sun_direction(
        struct mliAtmosphere *atmosphere,
        const double sunLatitude,
        const double sunHourAngle);
void mliAtmosphere_increase_latitude(
        struct mliAtmosphere *atmosphere,
        const double increment);
void mliAtmosphere_decrease_latitude(
        struct mliAtmosphere *atmosphere,
        const double increment);
void mliAtmosphere_increase_hours(
        struct mliAtmosphere *atmosphere,
        const double increment);
void mliAtmosphere_decrease_hours(
        struct mliAtmosphere *atmosphere,
        const double increment);
void mliAtmosphere_increase_altitude(
        struct mliAtmosphere *atmosphere,
        const double factor);
void mliAtmosphere_decrease_altitude(
        struct mliAtmosphere *atmosphere,
        const double factor);

struct mliColor mliAtmosphere_query(
        const struct mliAtmosphere *atmosphere,
        const struct mliVec orig,
        const struct mliVec dir);

struct mliColor mliAtmosphere_hit_earth_body(
        const struct mliAtmosphere *atmosphere,
        const struct mliVec orig,
        const struct mliVec dir);

struct mliColor mliAtmosphere_hit_outer_atmosphere(
        const struct mliAtmosphere *atmosphere,
        const struct mliVec orig,
        const struct mliVec dir,
        double tmin,
        double tmax);

struct mliColor mliAtmosphere_compute_depth(
        const struct mliAtmosphere *atmosphere,
        const struct mliVec orig,
        const struct mliVec dir,
        double tmin,
        double tmax);

#endif

/* mliCube */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLICUBE_H_
#define MLICUBE_H_


#define MLI_IS_BIT(var, pos) ((var) & (1 << (pos)))

struct mliCube {
        /*
         * Cubic Oriented-Bounding-Box
         * oriented w.r.t. the unit-vectors.
         */
        struct mliVec lower;
        double edge_length;
};

int mliCube_equal(const struct mliCube a, const struct mliCube b);
struct mliCube mliCube_octree_child_code(
        const struct mliCube cube,
        const uint8_t a);
struct mliCube mliCube_octree_child(
        const struct mliCube cube,
        const uint32_t sx,
        const uint32_t sy,
        const uint32_t sz);
struct mliCube mliCube_outermost_cube(const struct mliAABB a);
struct mliVec mliCube_center(const struct mliCube a);
struct mliAABB mliCube_to_aabb(const struct mliCube a);
struct mliVec mliCube_upper(const struct mliCube a);
#endif

/* mliDynFace */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIDYNARRAY_FACE_H_
#define MLIDYNARRAY_FACE_H_

MLIDYNARRAY_DEFINITON(mli, Face, struct mliFace)
#endif

/* mliDynMap */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIMAP_H_
#define MLIMAP_H_


struct mliDynMapItem {
        char key[MLI_NAME_CAPACITY];
        uint64_t value;
};
MLIDYNARRAY_DEFINITON(mli, Map, struct mliDynMapItem)
int mliDynMap_has(const struct mliDynMap *map, const char *key);
int mliDynMap_insert(struct mliDynMap *map, const char *key, uint64_t value);
int mliDynMap_find(const struct mliDynMap *map, const char *key, uint64_t *idx);
int mliDynMap_get(
        const struct mliDynMap *map,
        const char *key,
        uint64_t *value);

#endif

/* mliDynVec */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIDYNARRAY_VEC_H_
#define MLIDYNARRAY_VEC_H_

MLIDYNARRAY_DEFINITON(mli, Vec, struct mliVec)
#endif

/* mliFresnel */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIFRESNEL_H_
#define MLIFRESNEL_H_


struct mliFresnel {
        struct mliVec incident;
        struct mliVec normal;
        double n_from;
        double n_to;

        double _cosI;
        double _n_from_over_n_to;
        double _sinT2;
        double _cosT;
};

struct mliVec mliFresnel_refraction_direction(const struct mliFresnel fresnel);
struct mliVec mliFresnel_reflection_direction(const struct mliFresnel fresnel);
double mliFresnel_reflection_propability(const struct mliFresnel fresnel);
struct mliFresnel mliFresnel_init(
        const struct mliVec incident,
        const struct mliVec normal,
        const double n_from,
        const double n_to);
#endif

/* mliImage */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIIMAGE_H_
#define MLIIMAGE_H_


struct mliImage {
        uint32_t num_cols;
        uint32_t num_rows;
        struct mliColor *raw;
};

void mliImage_assign_pixel_colors_to_sum_and_exposure_image(
        const struct mliPixels *pixels,
        const struct mliImage *colors,
        struct mliImage *sum_image,
        struct mliImage *exposure_image);
int mliPixels_malloc_from_image_above_threshold(
        struct mliPixels *pixels,
        const struct mliImage *image,
        const float threshold);
void mliPixels_above_threshold(
        const struct mliImage *to_do_image,
        const float threshold,
        struct mliPixels *pixels);
void mliImage_from_sum_and_exposure(
        const struct mliImage *sum,
        const struct mliImage *exposure,
        struct mliImage *out);
void mliImage_luminance_threshold_dilatation(
        const struct mliImage *image,
        const float threshold,
        struct mliImage *out);
void mliImage_sobel(const struct mliImage *image, struct mliImage *out);
int mliImage_scale_down_twice(
        const struct mliImage *source,
        struct mliImage *destination);
struct mliColor mliImage_at(
        const struct mliImage *img,
        const uint32_t col,
        const uint32_t row);
void mliImage_set(
        const struct mliImage *img,
        const uint32_t col,
        const uint32_t row,
        const struct mliColor color);
void mliImage_set_all_pixel(
        const struct mliImage *img,
        const struct mliColor color);
void mliPixels_set_all_from_image(
        struct mliPixels *pixels,
        const struct mliImage *image);
uint32_t mliImage_idx(
        const struct mliImage *img,
        const uint32_t col,
        const uint32_t row);
int mliImage_malloc(
        struct mliImage *img,
        const uint32_t num_cols,
        const uint32_t num_rows);
void mliImage_free(struct mliImage *img);

void mliImage_copy(const struct mliImage *source, struct mliImage *destination);
void mliImage_fabs_difference(
        const struct mliImage *a,
        const struct mliImage *b,
        struct mliImage *out);
struct mliImage mliImage_init(void);

void mliImage_histogram(
        struct mliImage *img,
        const double *col_bin_edges,
        const double *row_bin_edges,
        const double x,
        const double y,
        const struct mliColor weight);
struct mliColor mliImage_max(const struct mliImage *img);
void mliImage_multiply(struct mliImage *img, const struct mliColor color);
void mliImage_divide_pixelwise(
        const struct mliImage *numerator,
        const struct mliImage *denominator,
        struct mliImage *out);
#endif

/* mliImage_ppm */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIIMAGE_PPM_H_
#define MLIIMAGE_PPM_H_


int mliImage_fwrite(const struct mliImage *img, FILE *f);
int mliImage_malloc_fread(struct mliImage *img, FILE *f);

int mliImage_malloc_from_path(struct mliImage *img, const char *path);
int mliImage_write_to_path(const struct mliImage *img, const char *path);
#endif

/* mliImage_print */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIIMAGE_PRINT_H_
#define MLIIMAGE_PRINT_H_


#define MLI_ASCII_MONOCHROME 100
#define MLI_ANSI_ESCAPE_COLOR 101

void mliImage_print(const struct mliImage *img, const uint64_t print_mode);
void mliImage_print_chars(
        const struct mliImage *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols,
        const uint64_t print_mode);

/* Colored ANSI escape sequences */

void mliImage_print_ansi_escape_chars(
        const struct mliImage *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols);

/* Monochrome ASCII chars */

void mliImage_print_ascii_chars(
        const struct mliImage *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols);
#endif

/* mliIntersection */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIINTERSECTION_H_
#define MLIINTERSECTION_H_


struct mliIntersection {
        struct mliGeometryId geometry_id;
        struct mliVec position_local;
        double distance_of_ray;
};

struct mliIntersection mliIntersection_init(void);
#endif

/* mliIntersectionSurfaceNormal */
/* ---------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIINTERSECTIONSURFACENORMAL_H_
#define MLIINTERSECTIONSURFACENORMAL_H_


struct mliIntersectionSurfaceNormal {
        struct mliGeometryId geometry_id;
        struct mliVec position;
        struct mliVec surface_normal;
        struct mliVec position_local;
        struct mliVec surface_normal_local;
        double distance_of_ray;
        int64_t from_outside_to_inside;
};

struct mliIntersectionSurfaceNormal mliIntersectionSurfaceNormal_init(void);

#endif

/* mliIo */
/* ----- */

/* Copyright 2018-2023 Sebastian Achim Mueller */
#ifndef mliIo_H_
#define mliIo_H_

struct mliIo {
        /* memory */
        unsigned char *cstr;

        /* Capacity of the allocated memory */
        uint64_t capacity;

        /* Size of the payload in the allocated memory */
        uint64_t size;

        /* Position of the cursor */
        uint64_t pos;
};

struct mliIo mliIo_init(void);
void mliIo_free(struct mliIo *byt);
int mliIo_malloc_capacity(struct mliIo *byt, const uint64_t capacity);
int mliIo_malloc(struct mliIo *byt);
int64_t mliIo_malloc_cstr(struct mliIo *byt, const char *s);
int mliIo_malloc_from_path(struct mliIo *byt, const char *path);
int mliIo_read_to_path(struct mliIo *byt, const char *path);
int mliIo_putc(struct mliIo *byt, const unsigned char c);
int mliIo_putchar(struct mliIo *byt, const char c);
int mliIo_getc(struct mliIo *byt);
int64_t mliIo_write(
        struct mliIo *byt,
        const void *ptr,
        const uint64_t size,
        const uint64_t count);
int64_t mliIo_read(
        struct mliIo *byt,
        const void *ptr,
        const uint64_t size,
        const uint64_t count);

uint64_t mliIo_ftell(struct mliIo *byt);
void mliIo_rewind(struct mliIo *byt);
int64_t mliIo_printf(struct mliIo *byt, const char *format, ...);

int mli_readline(
        struct mliIo *stream,
        struct mliStr *line,
        const char delimiter);

int mli_path_strip_this_dir(const struct mliStr *src, struct mliStr *dst);
int mli_path_basename(const struct mliStr *src, struct mliStr *dst);
int mli_path_splitext(
        const struct mliStr *src,
        struct mliStr *dst,
        struct mliStr *ext);

int mliStr_convert_line_break_CRLF_CR_to_LF(
        struct mliStr *dst,
        const struct mliStr *src);

int mli_line_viewer_write(
        struct mliIo *f,
        const struct mliStr *text,
        const uint64_t line_number,
        const uint64_t line_radius);

#endif

/* mliMat */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIMAT_H_
#define MLIMAT_H_


struct mliMat {
        double r00;
        double r01;
        double r02;
        double r10;
        double r11;
        double r12;
        double r20;
        double r21;
        double r22;
};

void mliMat_set(struct mliMat *a, uint64_t col, uint64_t row, const double v);
double mliMat_get(const struct mliMat *a, uint64_t col, uint64_t row);
struct mliMat mliMat_unity(void);
int mliMat_equal_margin(
        const struct mliMat a,
        const struct mliMat b,
        const double margin);
struct mliMat mliMat_init_axis_angle(
        const struct mliVec axis,
        const double angle);
struct mliMat mliMat_init_tait_bryan(
        const double rx,
        const double ry,
        const double rz);
struct mliMat mliMat_init_columns(
        const struct mliVec c0,
        const struct mliVec c1,
        const struct mliVec c2);
struct mliMat mliMat_covariance(
        const struct mliVec *vecs,
        const uint64_t num_vecs,
        const struct mliVec vecs_mean);
struct mliMat mliMat_transpose(const struct mliMat m);
struct mliMat mliMat_multiply(const struct mliMat x, const struct mliMat y);
struct mliMat mliMat_minor(const struct mliMat x, const int d);
struct mliMat mliMat_vector_outer_product(const struct mliVec v);
void mliMat_qr_decompose(
        const struct mliMat m,
        struct mliMat *q,
        struct mliMat *r);
int mliMat_has_shurform(const struct mliMat m, const double margin);
void mliMat_find_eigenvalues(
        const struct mliMat a,
        double *e0,
        double *e1,
        double *e2,
        const double margin,
        const uint64_t max_num_iterations);
int mliMat_find_eigenvector_for_eigenvalue(
        struct mliMat a,
        const double eigen_value,
        struct mliVec *eigen_vector,
        const double tolerance);
int mliMat_lup_decompose(struct mliMat *A, int *pivots, const double tolerance);
void mliMat_lup_solve(
        const struct mliMat *A,
        const int *P,
        const struct mliVec *b,
        struct mliVec *x);
#endif

/* mliObject */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIOBJECT_H_
#define MLIOBJECT_H_


struct mliObject {
        uint32_t num_vertices;
        struct mliVec *vertices;

        uint32_t num_vertex_normals;
        struct mliVec *vertex_normals;

        uint32_t num_faces;
        struct mliFace *faces_vertices;
        struct mliFace *faces_vertex_normals;
        uint16_t *faces_materials;

        uint32_t num_materials;
        struct mliName *material_names;
};

int mliObject_malloc(
        struct mliObject *obj,
        const uint64_t num_vertices,
        const uint64_t num_vertex_normals,
        const uint64_t num_faces,
        const uint64_t num_materials);
void mliObject_free(struct mliObject *obj);
struct mliObject mliObject_init(void);
int mliObject_equal(const struct mliObject *a, const struct mliObject *b);
uint32_t mliObject_resolve_material_idx(
        const struct mliObject *obj,
        const uint32_t face_idx);
#endif

/* mliObject_serialize */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIOBJECT_SERIALIZE_H_
#define MLIOBJECT_SERIALIZE_H_


int mliObject_fwrite(const struct mliObject *obj, FILE *f);
int mliObject_malloc_fread(struct mliObject *obj, FILE *f);
#endif

/* mliObject_valid */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIOBJECT_VALID_H_
#define MLIOBJECT_VALID_H_


int mliObject_is_valid(const struct mliObject *obj);
int mliObject_has_valid_vertices(const struct mliObject *obj);
int mliObject_has_valid_faces(const struct mliObject *obj);
int mliObject_has_valid_normals(
        const struct mliObject *obj,
        const double epsilon);
int mliObject_has_valid_materials(const struct mliObject *obj);
int mliObject_num_unused(
        const struct mliObject *obj,
        uint32_t *num_unused_vertices,
        uint32_t *num_unused_vertex_normals,
        uint32_t *num_unused_materials);
#endif

/* mliObject_wavefront */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIOBJECT_WAVEFRONT_H_
#define MLIOBJECT_WAVEFRONT_H_


int mliObject_malloc_from_wavefront(struct mliObject *obj, const char *str);
int mliObject_fprint_to_wavefront(struct mliIo *f, const struct mliObject *obj);
int mliObject_parse_face_line(
        const char *line,
        struct mliFace *faces_vertices,
        struct mliFace *faces_texture_points,
        struct mliFace *faces_vertex_normals,
        int *line_mode);
int mliObject_parse_three_float_line(const char *line, struct mliVec *v);
#endif

/* mliPixelWalk */
/* ------------ */

/* Copyright 2020-2021 Sebastian Achim Mueller */
#ifndef MLIPIXELWALK_H_
#define MLIPIXELWALK_H_


struct mliPixelWalk {
        /*
         * PixelWalk walks over the pixels of a rectangular image in a
         * cache-aware-way with respect to raytracing.
         * The goal is to bundle rays that will go to similar directions.
         * Instead of running fast along one axis of the image, and slow along
         * the other, PixelWalk spreads the walk among both axis by walking
         * small quadratic chunks of pixels.
         */
        uint32_t chunk_size;
        uint32_t num_chunks_row;
        uint32_t num_chunks_col;
        uint32_t chunk_row;
        uint32_t sub_row;
        uint32_t chunk_col;
        uint32_t sub_col;
        uint32_t num_rows;
        uint32_t num_cols;
        uint32_t i;
};

struct mliPixelWalk mliPixelWalk_set(
        const uint32_t num_cols,
        const uint32_t num_rows,
        const uint32_t chunk_size);
struct mliPixel mliPixelWalk_get(const struct mliPixelWalk *walk);
void mliPixelWalk_walk(struct mliPixelWalk *walk);
#endif

/* mliQuaternion */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIQUATERNION_H_
#define MLIQUATERNION_H_


struct mliQuaternion {
        double w;
        double x;
        double y;
        double z;
};

void mliQuaternion_print(const struct mliQuaternion q);
struct mliQuaternion mliQuaternion_set_tait_bryan(
        const double rx,
        const double ry,
        const double rz);
struct mliMat mliQuaternion_to_matrix(const struct mliQuaternion quat);
struct mliQuaternion mliQuaternion_set_rotaxis_and_angle(
        const struct mliVec rot_axis,
        const double angle);
double mliQuaternion_norm(const struct mliQuaternion q);
double mliQuaternion_product_complex_conjugate(const struct mliQuaternion p);
struct mliQuaternion mliQuaternion_product(
        const struct mliQuaternion p,
        const struct mliQuaternion q);
struct mliQuaternion mliQuaternion_complex_conjugate(
        const struct mliQuaternion q);
int mliQuaternion_equal_margin(
        const struct mliQuaternion a,
        const struct mliQuaternion b,
        const double margin);
int mliQuaternion_equal(
        const struct mliQuaternion a,
        const struct mliQuaternion b);
struct mliQuaternion mliQuaternion_set(
        const double w,
        const double x,
        const double y,
        const double z);
struct mliQuaternion mliQuaternion_set_unit_xyz(
        const double x,
        const double y,
        const double z);
#endif

/* mliRay */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIRAY_H_
#define MLIRAY_H_


struct mliRay {
        struct mliVec support;
        struct mliVec direction;
};

struct mliVec mliRay_at(const struct mliRay *ray, const double t);
struct mliRay mliRay_set(
        const struct mliVec support,
        const struct mliVec direction);
int mliRay_sphere_intersection(
        const struct mliVec support,
        const struct mliVec direction,
        const double radius,
        double *minus_solution,
        double *plus_solution);
#endif

/* mliRay_AABB */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIRAY_AABB_H_
#define MLIRAY_AABB_H_


int mliRay_has_overlap_aabb(
        const struct mliRay ray,
        const struct mliAABB aabb,
        double *ray_parameter);
#endif

/* mliTracer */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLITRACER_H_
#define MLITRACER_H_


struct mliScenery;
struct mliPrng;

struct mliTracerConfig {
        uint64_t num_trails_global_light_source;

        int have_atmosphere;
        struct mliAtmosphere atmosphere;

        struct mliColor background_color;
};

struct mliTracerConfig mliTracerConfig_init(void);

struct mliColor mli_trace(
        const struct mliScenery *scenery,
        const struct mliRay ray,
        const struct mliTracerConfig *config,
        struct mliPrng *prng);

struct mliColor mli_trace_with_atmosphere(
        const struct mliScenery *scenery,
        const struct mliRay ray,
        const struct mliTracerConfig *config,
        struct mliPrng *prng);

struct mliColor mli_trace_without_atmosphere(
        const struct mliScenery *scenery,
        const struct mliRay ray,
        const struct mliTracerConfig *config,
        struct mliPrng *prng);

double mli_trace_sun_obstruction(
        const struct mliScenery *scenery,
        const struct mliVec position,
        const struct mliTracerConfig *config,
        struct mliPrng *prng);

double mli_trace_sun_visibility(
        const struct mliScenery *scenery,
        const struct mliVec position,
        const struct mliTracerConfig *config,
        struct mliPrng *prng);

#endif

/* mliTriangle_AABB */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLITRIANGLE_AABB_H_
#define MLITRIANGLE_AABB_H_


struct mliTriangle {
        struct mliVec v1;
        struct mliVec v2;
        struct mliVec v3;
};

struct mliAABB mliTriangle_aabb(
        const struct mliVec a,
        const struct mliVec b,
        const struct mliVec c);
int mliTriangle_has_overlap_aabb(
        const struct mliVec a,
        const struct mliVec b,
        const struct mliVec c,
        const struct mliAABB aabb);
struct mliTriangle mliTriangle_set_in_norm_aabb(
        const struct mliVec a,
        const struct mliVec b,
        const struct mliVec c,
        const struct mliAABB aabb);
int64_t mliTriangle_intersects_norm_aabb(struct mliTriangle t);
int64_t mliTriangle_intersects_point(struct mliTriangle t, struct mliVec p);
int64_t mli_triangle_aabb_check_line(
        struct mliVec p1,
        struct mliVec p2,
        int64_t outcode_diff);
int64_t mli_triangle_aabb_check_point(
        struct mliVec p1,
        struct mliVec p2,
        double alpha,
        int64_t mask);
int64_t mli_triangle_aabb_bevel_3d(struct mliVec p);
int64_t mli_triangle_aabb_bevel_2d(struct mliVec p);
int64_t mli_triangle_aabb_face_plane(struct mliVec p);
#endif

/* mliUserScenery */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIUSERSCENERY_H_
#define MLIUSERSCENERY_H_


struct mliNameMap {
        struct mliDynMap objects;
        struct mliDynMap media;
        struct mliDynMap surfaces;
        struct mliDynMap boundary_layers;
};
struct mliNameMap mliNameMap_init(void);
int mliNameMap_malloc(struct mliNameMap *namemap);
void mliNameMap_free(struct mliNameMap *namemap);

struct mliMaterials;
struct mliArchive;
int mliMaterials_malloc_form_archive(
        struct mliMaterials *materials,
        struct mliNameMap *names,
        const struct mliArchive *archive);
struct mliGeometry;
int mli_set_geometry_objects_and_names_from_archive(
        struct mliGeometry *geometry,
        struct mliDynMap *object_names,
        const struct mliArchive *archive);

struct mliObject;
struct mliFrame;
int mli_check_malloc_root_frame_from_Archive(
        struct mliFrame *root,
        const struct mliArchive *archive,
        const struct mliDynMap *object_names,
        const struct mliObject *objects,
        const struct mliDynMap *boundary_layer_names);
#endif

/* mli_json */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_JSON_H_
#define MLI_JSON_H_


struct mliJson {
        struct mliStr raw;
        uint64_t num_tokens;
        struct jsmntok_t *tokens;
};

int mliJson_debug_to_path(const struct mliJson *json, const char *path);
int mliJson_debug_fprint(FILE *f, const struct mliJson *json);
int mliJson_debug_token_fprint(
        FILE *f,
        const struct mliJson *json,
        const uint64_t token);
uint64_t mliJson_token_by_index(
        const struct mliJson *json,
        const uint64_t start_token_idx,
        const uint64_t child_idx);
int mliJson_token_by_key(
        const struct mliJson *json,
        const uint64_t token,
        const char *key,
        uint64_t *key_token);
int mliJson_token_by_key_eprint(
        const struct mliJson *json,
        const uint64_t token,
        const char *key,
        uint64_t *key_token);
int mliJson_double_by_token(
        const struct mliJson *json,
        const uint64_t token,
        double *val);
int mliJson_double_by_key(
        const struct mliJson *json,
        const uint64_t token,
        double *val,
        const char *key);
int mliJson_int64_by_token(
        const struct mliJson *json,
        const uint64_t token,
        int64_t *return_int64);
int mliJson_uint64_by_token(
        const struct mliJson *json,
        const uint64_t token,
        uint64_t *return_uint64);
int mliJson_cstr_by_token(
        const struct mliJson *json,
        const uint64_t token,
        char *return_string,
        const uint64_t return_string_size);
int mliJson_int64_by_key(
        const struct mliJson *json,
        const uint64_t token,
        int64_t *val,
        const char *key);
int mliJson_uint64_by_key(
        const struct mliJson *json,
        const uint64_t token,
        uint64_t *val,
        const char *key);
int mliJson_cstrcmp(
        const struct mliJson *json,
        const uint64_t token,
        const char *str);
int mliJson_malloc_from_path(struct mliJson *json, const char *path);
int mliJson_malloc_from_cstr(struct mliJson *json, const char *cstr);
int mliJson_malloc_tokens__(struct mliJson *json);
int mliJson_parse_tokens__(struct mliJson *json);
void mliJson_free(struct mliJson *json);
struct mliJson mliJson_init(void);
#endif

/* mli_random_PCG32 */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RANDOM_PCG32_H_
#define MLI_RANDOM_PCG32_H_


/*      Wrapping the pcg implementation by Melissa O'Neill in
 *      pcg_variants_32bit_subset.h
 */

struct mliPCG32 {
        struct pcg_state_setseq_64 state_setseq_64;
};

struct mliPCG32 mliPCG32_init(const uint32_t seed);
uint32_t mliPCG32_generate_uint32(struct mliPCG32 *pcg32);
void mliPCG32_reinit(struct mliPCG32 *pcg32, const uint32_t seed);

#endif

/* mli_random_generator */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RANDOM_GENERATOR_H_
#define MLI_RANDOM_GENERATOR_H_


/**
 *      mliPrng is a transparent container to use different
 *      pseudo-random-number-generators (PRNGs) within merlict.
 *      It defines a minimal interface:
 *
 *      1) (Re)initializing with a seed.
 *      2) Generating the next random number uint32.
 *
 *      Merlict ships with two Prngs:
 *
 *      1) Mersenne Twister 19937
 *      2) PCG32
 *
 *      If you want to use your own, wrapp it here using mliPrng. See below.
 */

union mliPrngStorage {
        struct mliMT19937 mt19937;
        struct mliPCG32 pcg32;
        /* Add your own prng here */
};

struct mliPrng {
        union mliPrngStorage _storage;
        uint32_t (*generate_uint32)(void *);
        void (*reinit)(void *, const uint32_t);
};

uint32_t mliPrng_generate_uint32(struct mliPrng *prng);
void mliPrng_reinit(struct mliPrng *prng, const uint32_t seed);

/**
 *      Mersenne Twister 19937
 *      ----------------------
 */
struct mliPrng mliPrng_init_MT19937(const uint32_t seed);
uint32_t mliPrng_MT19937_generate_uint32(void *mt);
void mliPrng_MT19937_reinit(void *mt, const uint32_t seed);

/**
 *      PCG32
 *      -----
 */
struct mliPrng mliPrng_init_PCG32(const uint32_t seed);
uint32_t mliPrng_PCG32_generate_uint32(void *pcg);
void mliPrng_PCG32_reinit(void *pcg, const uint32_t seed);

/**
 *      Add your own prng here
 *      ----------------------
 */

#endif

/* mli_triangle_intersection */
/* ------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLITRIANGLE_INTERSECTION_H_
#define MLITRIANGLE_INTERSECTION_H_


struct mliVec mliTriangle_interpolate_surface_normal(
        const struct mliVec vertex_normal_a,
        const struct mliVec vertex_normal_b,
        const struct mliVec vertex_normal_c,
        const struct mliBarycentrigWeights weights);

int mliRay_intersects_triangle(
        const struct mliRay ray,
        const struct mliVec vertex_a,
        const struct mliVec vertex_b,
        const struct mliVec vertex_c,
        double *intersection_ray_parameter);

struct mliVec mliTriangle_surface_normal(
        const struct mliVec vertex_a,
        const struct mliVec vertex_b,
        const struct mliVec vertex_c,
        const struct mliVec vertex_normal_a,
        const struct mliVec vertex_normal_b,
        const struct mliVec vertex_normal_c,
        const struct mliVec intersection_position);
#endif

/* mliArchive */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef ARCHIVE_H_
#define ARCHIVE_H_


MLIDYNARRAY_DEFINITON(mli, TextFiles, struct mliStr)

struct mliArchive {
        struct mliDynTextFiles textfiles;
        struct mliDynMap filenames;
};

struct mliArchive mliArchive_init(void);

void mliArchive_free(struct mliArchive *arc);
int mliArchive_malloc(struct mliArchive *arc);
int mliArchive_malloc_fread(struct mliArchive *arc, FILE *f);
int mliArchive_malloc_from_path(struct mliArchive *arc, const char *path);

int mliArchive_push_back(
        struct mliArchive *arc,
        const struct mliStr *filename,
        const struct mliStr *payload);

int mliArchive_has(const struct mliArchive *arc, const char *filename);
int mliArchive_get(
        const struct mliArchive *arc,
        const char *filename,
        struct mliStr **str);
int mliArchive_get_malloc_json(
        const struct mliArchive *arc,
        const char *filename,
        struct mliJson *json);
uint64_t mliArchive_num(const struct mliArchive *arc);
void mliArchive_info_fprint(FILE *f, const struct mliArchive *arc);
void mliArchive_mask_filename_prefix_sufix(
        const struct mliArchive *arc,
        uint64_t *mask,
        const char *prefix,
        const char *suffix);
uint64_t mliArchive_num_filename_prefix_sufix(
        const struct mliArchive *arc,
        const char *prefix,
        const char *sufix);

#endif

/* mliAtmosphere_json */
/* ------------------ */

/* Copyright 2018-2021 Sebastian Achim Mueller */
#ifndef MLIATMOSPHERE_JSON_H_
#define MLIATMOSPHERE_JSON_H_


int mliAtmosphere_from_json_token(
        struct mliAtmosphere *atm,
        const struct mliJson *json,
        const uint64_t tkn);

#endif

/* mliColor_json */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLICOLOR_JSON_H_
#define MLICOLOR_JSON_H_


int mliColor_from_json_token(
        struct mliColor *c,
        const struct mliJson *json,
        const uint64_t token);
#endif

/* mliDynMap_json */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIDYNMAP_JSON_H_
#define MLIDYNMAP_JSON_H_


int mliDynMap_get_value_for_string_from_json(
        const struct mliDynMap *map,
        const struct mliJson *json,
        const uint64_t token_name,
        uint32_t *out_value);
int mliDynMap_insert_key_from_json(
        struct mliDynMap *map,
        const struct mliJson *json,
        const uint64_t token_name,
        const uint64_t value);

#endif

/* mliFunc_json */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIFUNC_JSON_H_
#define MLIFUNC_JSON_H_


int mliFunc_malloc_from_json_token(
        struct mliFunc *func,
        const struct mliJson *json,
        const uint64_t token);
#endif

/* mliHomTra */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIHOMTRA_H_
#define MLIHOMTRA_H_


struct mliHomTraComp {
        struct mliVec translation;
        struct mliQuaternion rotation;
};

struct mliHomTra {
        struct mliVec translation;
        struct mliMat rotation;
};

void mliHomTra_print(const struct mliHomTra h);
struct mliHomTraComp mliHomTraComp_set(
        const struct mliVec translation,
        const struct mliQuaternion rotation);
struct mliHomTraComp mliHomTraComp_sequence(
        const struct mliHomTraComp a,
        const struct mliHomTraComp b);
int mliHomTraComp_equal(
        const struct mliHomTraComp a,
        const struct mliHomTraComp b);
struct mliVec mliHomTra_dir_inverse(
        const struct mliHomTra *t,
        const struct mliVec in);
struct mliVec mliHomTra_dir(const struct mliHomTra *t, const struct mliVec in);
struct mliVec mliHomTra_pos_inverse(
        const struct mliHomTra *t,
        const struct mliVec in);
struct mliVec mliHomTra_pos(const struct mliHomTra *t, const struct mliVec in);
struct mliRay mliHomTra_ray_inverse(
        const struct mliHomTra *t,
        const struct mliRay in);
struct mliRay mliHomTra_ray(const struct mliHomTra *t, const struct mliRay in);
struct mliRay mli_transform_ray_inverse(
        const struct mliMat *rotation,
        const struct mliVec translation,
        const struct mliRay in);
struct mliRay mli_transform_ray(
        const struct mliMat *rotation,
        const struct mliVec translation,
        const struct mliRay in);
struct mliVec mli_transform_position_inverse(
        const struct mliMat *rotation,
        const struct mliVec translation,
        const struct mliVec in);
struct mliVec mli_transform_position(
        const struct mliMat *rotation,
        const struct mliVec translation,
        const struct mliVec in);
struct mliVec mli_transform_orientation_inverse(
        const struct mliMat *rotation,
        const struct mliVec in);
struct mliVec mli_transform_orientation(
        const struct mliMat *rotation,
        const struct mliVec in);
struct mliHomTra mliHomTra_from_compact(const struct mliHomTraComp trafo);
#endif

/* mliMedium_json */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIMEDIUM_JSON_H_
#define MLIMEDIUM_JSON_H_


int mliMedium_malloc_from_json_str(struct mliMedium *med, const char *json_str);
int mliMedium_malloc_from_json_token(
        struct mliMedium *med,
        const struct mliJson *json,
        const uint64_t token);
#endif

/* mliObject_AABB */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIOBJECT_AABB_H_
#define MLIOBJECT_AABB_H_


int mliObject_has_overlap_aabb(
        const struct mliObject *obj,
        const struct mliHomTra local2root,
        const struct mliAABB aabb);

struct mliAABB mliObject_aabb(
        const struct mliObject *obj,
        const struct mliHomTra local2root);

int mliObject_face_in_local_frame_has_overlap_aabb(
        const struct mliObject *obj,
        const uint64_t face_idx,
        const struct mliAABB aabb);

int mliObject_face_in_local_frame_has_overlap_aabb_void(
        const void *obj,
        const uint32_t face_idx,
        const struct mliAABB aabb);

struct mliAABB mliObject_aabb_in_local_frame(const struct mliObject *obj);

#endif

/* mliPhoton */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIPHOTON_H_
#define MLIPHOTON_H_


struct mliPhoton {
        struct mliRay ray;
        double wavelength;
        int64_t id;
};

#endif

/* mliPinHoleCamera */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PIN_HOLE_CAMERA_H_
#define MLI_PIN_HOLE_CAMERA_H_


struct mliPinHoleCamera {
        struct mliVec optical_axis;
        struct mliVec col_axis;
        struct mliVec row_axis;
        struct mliVec principal_point;
        double distance_to_principal_point;
        double row_over_column_pixel_ratio;
};

struct mliPinHoleCamera mliPinHoleCamera_init(
        const double field_of_view,
        const struct mliImage *image,
        const double row_over_column_pixel_ratio);

void mliPinHoleCamera_render_image(
        struct mliPinHoleCamera camera,
        const struct mliHomTraComp camera2root_comp,
        const struct mliScenery *scenery,
        struct mliImage *image,
        const struct mliTracerConfig *tracer_config,
        struct mliPrng *prng);

void mliPinHoleCamera_render_image_with_view(
        const struct mliView view,
        const struct mliScenery *scenery,
        struct mliImage *image,
        const double row_over_column_pixel_ratio,
        const struct mliTracerConfig *tracer_config,
        struct mliPrng *prng);

struct mliRay mliPinHoleCamera_ray_at_row_col(
        const struct mliPinHoleCamera *camera,
        const struct mliImage *image,
        const uint32_t row,
        const uint32_t col);

#endif

/* mliQuaternion_json */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIQUATERNION_JSON_H_
#define MLIQUATERNION_JSON_H_


int mliQuaternion_tait_bryan_from_json(
        struct mliQuaternion *quat,
        const struct mliJson *json,
        const uint64_t token);
int mliQuaternion_axis_angle_from_json(
        struct mliQuaternion *quat,
        const struct mliJson *json,
        const uint64_t token);
int mliQuaternion_quaternion_from_json(
        struct mliQuaternion *quat,
        const struct mliJson *json,
        const uint64_t token);
int mliQuaternion_from_json(
        struct mliQuaternion *quat,
        const struct mliJson *json,
        const uint64_t token);

#endif

/* mliSurface */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLISURFACE_H_
#define MLISURFACE_H_


#define MLI_MATERIAL_PHONG 100u
#define MLI_MATERIAL_TRANSPARENT 102u

struct mliSurface {
        uint32_t material;

        struct mliFunc specular_reflection;
        struct mliFunc diffuse_reflection;

        /*
         *  The color is only relevant for fast rendering of images.
         *  Color will not effect the propagation of photons.
         */
        struct mliColor color;
};

int mliSurface_malloc(
        struct mliSurface *surface,
        const uint32_t num_points_specular_reflection,
        const uint32_t num_points_diffuse_reflection);
void mliSurface_free(struct mliSurface *surface);
struct mliSurface mliSurface_init(void);
int mliSurface_equal(const struct mliSurface *a, const struct mliSurface *b);

int mliSurface_fwrite(const struct mliSurface *srf, FILE *f);
int mliSurface_malloc_fread(struct mliSurface *srf, FILE *f);

int mli_material_type_to_string(const uint32_t type, char *s);
int mli_material_type_from_string(const char *s, uint32_t *id);

int mliSurface_malloc_from_json_str(
        struct mliSurface *surface,
        const char *json_str);
int mliSurface_malloc_from_json_token(
        struct mliSurface *surface,
        const struct mliJson *json,
        const uint64_t token);
#endif

/* mliSurface_json */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLISURFACE_JSON_H_
#define MLISURFACE_JSON_H_


int mliSurface_malloc_from_json_str(
        struct mliSurface *surface,
        const char *json_str);
int mliSurface_malloc_from_json_token(
        struct mliSurface *surface,
        const struct mliJson *json,
        const uint64_t token);
int mli_material_type_from_json_token(
        const struct mliJson *json,
        const uint64_t token,
        uint32_t *material);
#endif

/* mliTracerConfig_json */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLITRACERCONFIG_JSON_H_
#define MLITRACERCONFIG_JSON_H_


int mliTracerConfig_from_json_token(
        struct mliTracerConfig *tc,
        const struct mliJson *json,
        const uint64_t tkn);
#endif

/* mliUserScenery_json */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIUSERSCENERY_JSON_H_
#define MLIUSERSCENERY_JSON_H_


int mliSide_from_json(
        struct mliSide *side,
        const struct mliDynMap *surface_names,
        const struct mliDynMap *medium_names,
        const struct mliJson *json,
        const uint64_t side_token);
int mliBoundaryLayer_from_json(
        struct mliBoundaryLayer *boundary_layer,
        const struct mliDynMap *surface_names,
        const struct mliDynMap *medium_names,
        const struct mliJson *json,
        const uint64_t token_surface);
int mliMaterials_assign_boundary_layers_from_json(
        struct mliMaterials *materials,
        struct mliDynMap *boundary_layer_names,
        const struct mliDynMap *surface_names,
        const struct mliDynMap *medium_names,
        const struct mliJson *json);

#endif

/* mliVec_json */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIVEC_JSON_H_
#define MLIVEC_JSON_H_


int mliVec_from_json_token(
        struct mliVec *v,
        const struct mliJson *json,
        const uint64_t token);
#endif

/* mli_random */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RANDOM_H_
#define MLI_RANDOM_H_


struct mliRandomUniformRange {
        double start;
        double range;
};

struct mliRandomZenithRange {
        double z_min;
        double z_range;
};

struct mliVec mli_random_position_on_disc(
        const double radius,
        struct mliPrng *prng);
struct mliVec mli_random_draw_direction_in_zenith_azimuth_range(
        const struct mliRandomZenithRange zenith,
        const struct mliRandomUniformRange azimuth,
        struct mliPrng *prng);
double mli_random_draw_zenith(
        const struct mliRandomZenithRange range,
        struct mliPrng *prng);
struct mliRandomZenithRange mliRandomZenithRange_set(
        const double min_zenith_distance,
        const double max_zenith_distance);
double mli_random_draw_uniform(
        const struct mliRandomUniformRange uniform_range,
        struct mliPrng *prng);
struct mliRandomUniformRange mliRandomUniformRange_set(
        double start,
        double stop);
double mli_random_normal_Irwin_Hall_approximation(
        struct mliPrng *prng,
        const double mean,
        const double std);
double mli_random_expovariate(struct mliPrng *prng, const double rate);
double mli_random_uniform(struct mliPrng *prng);
struct mliVec mli_random_position_inside_unit_sphere(struct mliPrng *prng);
#endif

/* mliApertureCamera */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIAPTCAM_H_
#define MLIAPTCAM_H_


/*
principal-rays of the thin-lens
===============================
                                        | z
                                        |
        (A)                           --+-- object-distance
          \                             |
         | \\                           |
            \\                          |
         |   \\                         |
              \ \                       |
         |     \ \                      |
                \ \                     |
         |       \  \                   |
                 \   \                  |
         |        \   \                 |
                   \    \               |
         |          \    \              |
                     \    \             |
         |            \     \           |
                       \     \          |
         |             \      \         |
                        \       \       |
         |               \       \      |
                          \       \     |
         |                 \        \   |
                            \        \  |
         |                   \        \ |
                             \         \|
         |                    \       --+--  focal-length
                               \        |\
         |                      \       | \
                                 \      |  \
         |                        \     |   \
                                  \     |     \
         |                         \    |      \
                                    \   |       \
         |                           \  |         \
                                      \ |          \
         |                             \|           \      aperture-plane
   -|----O------------------------------O------------O----------------------|-
          \                             |\                                  |
             \                          | \          |                aperture-
                \                       |  \                           radius
                   \                    |   \        |
                      \                 |    \
                         \              |     \      |
                            \           |     \
                               \        |      \     |
                                  \     |       \
                                     \  |        \   |
                        focal-length  --+--       \
                                        | \       \  |
                                        |    \     \
                                        |       \   \|
    image-sensor-plane                  |          \\
                ------------------------+-----------(P)----------  x/y
                                        |\_ image-sensor-distance
                                        |

1)      Find point P on image-sensor-plane for (row, column).
        With P.z = -image-sensor-distance.
        Add random-scatter in pixel-bin.

2)      Find point A on the object-plane.
        With A.z = +object-distance

3)      Draw random point W on aperture-plane within aperture-radius.

4)      Trace ray(P - W) and assign to pixel (row, column).

*/

struct mliVec mliApertureCamera_pixel_center_on_image_sensor_plane(
        const double image_sensor_width_x,
        const double image_sensor_width_y,
        const double image_sensor_distance,
        const uint64_t num_pixel_x,
        const uint64_t num_pixel_y,
        const uint64_t pixel_x,
        const uint64_t pixel_y);

struct mliVec mliApertureCamera_pixel_support_on_image_sensor_plane(
        const double image_sensor_width_x,
        const double image_sensor_width_y,
        const double image_sensor_distance,
        const uint64_t num_pixel_x,
        const uint64_t num_pixel_y,
        const uint64_t pixel_x,
        const uint64_t pixel_y,
        struct mliPrng *prng);

struct mliVec mliApertureCamera_get_object_point(
        const double focal_length,
        const struct mliVec pixel_support);

double mli_thin_lens_get_object_given_focal_and_image(
        const double focal_length,
        const double image_distance);

double mli_thin_lens_get_image_given_focal_and_object(
        const double focal_length,
        const double object_distance);

double mliApertureCamera_focal_length_given_field_of_view_and_sensor_width(
        const double field_of_view,
        const double image_sensor_width);

struct mliVec mliApertureCamera_ray_support_on_aperture(
        const double aperture_radius,
        struct mliPrng *prng);

struct mliRay mliApertureCamera_get_ray_for_pixel(
        const double focal_length,
        const double aperture_radius,
        const double image_sensor_distance,
        const double image_sensor_width_x,
        const double image_sensor_width_y,
        const uint64_t num_pixel_x,
        const uint64_t num_pixel_y,
        const uint64_t pixel_x,
        const uint64_t pixel_y,
        struct mliPrng *prng);

struct mliApertureCamera {
        double focal_length;
        double aperture_radius;
        double image_sensor_distance;
        double image_sensor_width_x;
        double image_sensor_width_y;
};

struct mliApertureCamera mliApertureCamera_init(void);

int mliApertureCamera_render_image(
        const struct mliApertureCamera camera,
        const struct mliHomTraComp camera2root_comp,
        const struct mliScenery *scenery,
        struct mliImage *image,
        const struct mliTracerConfig *tracer_config,
        struct mliPrng *prng);

void mliApertureCamera_aquire_pixels(
        const struct mliApertureCamera camera,
        const struct mliImage *image,
        const struct mliHomTraComp camera2root_comp,
        const struct mliScenery *scenery,
        const struct mliPixels *pixels_to_do,
        struct mliImage *colors,
        const struct mliTracerConfig *tracer_config,
        struct mliPrng *prng);

void mliApertureCamera_assign_pixel_colors_to_sum_and_exposure_image(
        const struct mliPixels *pixels,
        const struct mliImage *colors,
        struct mliImage *sum_image,
        struct mliImage *exposure_image);

#endif

/* mliDynPhoton */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIDYNPHOTON_H_
#define MLIDYNPHOTON_H_

MLIDYNARRAY_DEFINITON(mli, Photon, struct mliPhoton)
#endif

/* mliFrame */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIFRAME_H_
#define MLIFRAME_H_


#define MLI_FRAME 1000u
#define MLI_OBJECT 1001u

MLIDYNARRAY_DEFINITON(mli, FramePtr, struct mliFrame *)

struct mliFrame {
        uint32_t type;
        uint32_t id;
        struct mliHomTraComp frame2mother;
        struct mliHomTraComp frame2root;
        struct mliFrame *mother;

        struct mliDynFramePtr children;

        uint32_t object;
        struct mliDynUint32 boundary_layers;
};

void mliFrame_set_frame2root(struct mliFrame *f);
void mliFrame_print(struct mliFrame *f);
void mliFrame_print_walk(const struct mliFrame *f, const uint64_t indention);
int mli_string_to_type(const char *s, uint64_t *type);
int mli_type_to_string(const uint64_t type, char *s);
struct mliFrame *mliFrame_add(struct mliFrame *mother, const uint64_t type);
int mliFrame_set_mother_and_child(
        struct mliFrame *mother,
        struct mliFrame *child);
int mliFrame_malloc(struct mliFrame *f, const uint64_t type);
void mliFrame_free(struct mliFrame *f);
struct mliFrame mliFrame_init(void);
#endif

/* mliFrame_json */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIFRAME_JSON_H_
#define MLIFRAME_JSON_H_


int mliFrame_from_json(
        struct mliFrame *mother,
        const struct mliJson *json,
        const uint64_t token_children,
        const struct mliDynMap *object_names,
        const struct mliObject *objects,
        const struct mliDynMap *boundary_layer_names);
int mliFrame_id_from_json_token(
        uint32_t *id,
        const struct mliJson *json,
        const uint64_t token);
int mliFrame_pos_rot_from_json_token(
        struct mliHomTraComp *frame2mother,
        const struct mliJson *json,
        const uint64_t token);
int mliFrame_type_from_json_token(
        uint64_t *type,
        const struct mliJson *json,
        const uint64_t token);
int mliFrame_boundary_layers_form_json_token(
        struct mliDynUint32 *boundary_layers,
        const uint32_t object_idx,
        const struct mliObject *objects,
        const struct mliDynMap *boundary_layer_names,
        const struct mliJson *json,
        const uint64_t token);
int mliFrame_object_reference_form_json_token(
        uint32_t *object_reference,
        const struct mliJson *json,
        const uint64_t token,
        const struct mliDynMap *object_names);
#endif

/* mliMaterials */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIMATERIALS_H_
#define MLIMATERIALS_H_


struct mliMaterialsCapacity {
        uint64_t num_surfaces;
        uint64_t num_media;
        uint64_t num_boundary_layers;
};

struct mliMaterialsCapacity mliMaterialsCapacity_init(void);

struct mliMaterials {
        uint64_t num_surfaces;
        struct mliSurface *surfaces;
        struct mliName *surface_names;

        uint64_t num_media;
        struct mliMedium *media;
        struct mliName *medium_names;

        uint64_t num_boundary_layers;
        struct mliBoundaryLayer *boundary_layers;
        struct mliName *boundary_layer_names;

        uint64_t default_medium;
};

int mliMaterials_malloc(
        struct mliMaterials *res,
        const struct mliMaterialsCapacity rescap);
void mliMaterials_free(struct mliMaterials *res);
struct mliMaterials mliMaterials_init(void);
void mliMaterials_info_fprint(FILE *f, const struct mliMaterials *res);
#endif

/* mliMaterials_equal */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIMATERIALS_EQUAL_H_
#define MLIMATERIALS_EQUAL_H_


int mliMaterials_equal(
        const struct mliMaterials *a,
        const struct mliMaterials *b);
int mliMaterials_surfaces_equal(
        const struct mliMaterials *a,
        const struct mliMaterials *b);
int mliMaterials_media_equal(
        const struct mliMaterials *a,
        const struct mliMaterials *b);
int mliMaterials_boundary_layers_equal(
        const struct mliMaterials *a,
        const struct mliMaterials *b);
#endif

/* mliMaterials_serialize */
/* ---------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIMATERIALS_SERIALIZE_H_
#define MLIMATERIALS_SERIALIZE_H_


int mliMaterials_fwrite(const struct mliMaterials *res, FILE *f);
int mliMaterials_malloc_fread(struct mliMaterials *res, FILE *f);
#endif

/* mliMaterials_valid */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIMATERIALS_VALID_H_
#define MLIMATERIALS_VALID_H_


int mliMaterials_valid(const struct mliMaterials *materials);
int mliMaterials_valid_surfaces(const struct mliMaterials *materials);
int mliMaterials_valid_media(const struct mliMaterials *materials);
int mliMaterials_valid_boundary_layers(const struct mliMaterials *materials);
#endif

/* mliPhotonInteraction */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIPHOTONINTERACTION_H_
#define MLIPHOTONINTERACTION_H_


#define MLI_VACUUM_SPPED_OF_LIGHT 299792458.0
#define MLI_PHOTON_CREATION 101u
#define MLI_PHOTON_ABSORBTION 102u
#define MLI_PHOTON_ABSORBTION_MEDIUM 103u
#define MLI_PHOTON_FRESNEL_REFLECTION 104u
#define MLI_PHOTON_REFRACTION 105u
#define MLI_PHOTON_SPECULAR_REFLECTION 106u
#define MLI_PHOTON_DIFFUSE_REFLECTION 107u

struct mliPhotonInteraction {
        int32_t on_geometry_surface;
        struct mliGeometryId geometry_id;

        struct mliVec position;
        struct mliVec position_local;
        double distance_of_ray;

        uint64_t medium_coming_from;
        uint64_t medium_going_to;

        int32_t from_outside_to_inside;
        int32_t type;
};

int mli_time_of_flight(
        const struct mliMaterials *materials,
        const struct mliPhotonInteraction *phisec,
        const double wavelength,
        double *time_of_flight);
int mli_photoninteraction_type_to_string(const int32_t type, char *s);
#endif

/* mliRenderConfig */
/* --------------- */

/* Copyright 2018-2021 Sebastian Achim Mueller */
#ifndef MLIRENDERCONFIG_H_
#define MLIRENDERCONFIG_H_


struct mliRenderConfig {
        struct mliApertureCamera camera;
        struct mliHomTraComp camera_to_root;
        struct mliTracerConfig tracer;
        uint64_t num_pixel_x;
        uint64_t num_pixel_y;
        uint64_t random_seed;
};

struct mliRenderConfig mliRenderConfig_init(void);

int mliRenderConfig_from_json_token(
        struct mliRenderConfig *cc,
        const struct mliJson *json,
        const uint64_t token);

#endif

/* mli_lambertian_cosine_law */
/* ------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_LAMBERTIAN_COSINE_LAW_H_
#define MLI_LAMBERTIAN_COSINE_LAW_H_


struct mliVec mli_draw_lambertian_direction_wrt_surface_normal(
        struct mliPrng *prng,
        const struct mliVec surface_normal);
struct mliVec mli_draw_lambertian_direction_wrt_z(struct mliPrng *prng);
#endif

/* mli_photon_sources */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PHOTON_SOURCES_H_
#define MLI_PHOTON_SOURCES_H_


int mli_photon_source_point_like_opening_cone_towards_z(
        struct mliDynPhoton *out_photons,
        const double wavelength,
        const double opening_angle,
        const uint64_t num_photons,
        struct mliPrng *prng);
int mli_photon_source_parallel_towards_z_from_xy_disc(
        struct mliDynPhoton *out_photons,
        const double wavelength,
        const double radius,
        const uint64_t num_photons,
        struct mliPrng *prng);
#endif

/* mliGeometry */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRY_H_
#define MLIGEOMETRY_H_


struct mliGeometry {
        uint32_t num_objects;
        struct mliObject *objects;
        struct mliName *object_names;

        uint32_t num_robjects;
        uint32_t *robjects;
        uint32_t *robject_ids;
        struct mliHomTraComp *robject2root;
};

int mliGeometry_malloc(
        struct mliGeometry *geometry,
        const uint32_t num_objects,
        const uint32_t num_robjects);
int mliGeometry_malloc_references(
        struct mliGeometry *geometry,
        const uint32_t num_robjects);
int mliGeometry_malloc_objects(
        struct mliGeometry *geometry,
        const uint32_t num_objects);

void mliGeometry_free(struct mliGeometry *geometry);
void mliGeometry_free_objects(struct mliGeometry *geometry);
void mliGeometry_free_references(struct mliGeometry *geometry);

struct mliGeometry mliGeometry_init(void);
void mliGeometry_init_objects(struct mliGeometry *geometry);
void mliGeometry_init_references(struct mliGeometry *geometry);

void mliGeometry_info_fprint(FILE *f, const struct mliGeometry *geometry);
struct mliBoundaryLayer mliGeometry_object_surfaces(
        const struct mliGeometry *geometry,
        const uint32_t object_idx);
int mliGeometry_warn_objects(const struct mliGeometry *geometry);
#endif

/* mliGeometryToMaterialMap */
/* ------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRYTOMATERIALMAP_H_
#define MLIGEOMETRYTOMATERIALMAP_H_


struct mliGeometryToMaterialMap {
        uint32_t num_robjects;
        uint32_t total_num_boundary_layers;
        uint32_t *boundary_layers;
        uint32_t *first_boundary_layer_in_robject;
};

struct mliGeometryToMaterialMap mliGeometryToMaterialMap_init(void);
int mliGeometryToMaterialMap_malloc(
        struct mliGeometryToMaterialMap *map,
        const uint32_t num_robjects,
        const uint32_t total_num_boundary_layers);
void mliGeometryToMaterialMap_free(struct mliGeometryToMaterialMap *map);

uint32_t mliGeometryToMaterialMap_resolve_idx(
        const struct mliGeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx);

uint32_t mliGeometryToMaterialMap_get(
        const struct mliGeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx);
void mliGeometryToMaterialMap_set(
        const struct mliGeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx,
        const uint32_t boundary_layer_idx);

uint32_t mliGeometryToMaterialMap_num_boundary_layers_in_robject(
        const struct mliGeometryToMaterialMap *map,
        const uint32_t robject_idx);

void mliGeometryToMaterialMap_info_fprint(
        FILE *f,
        const struct mliGeometryToMaterialMap *map);
#endif

/* mliGeometryToMaterialMap_equal */
/* ------------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRYTOMATERIALMAP_EQUAL_H_
#define MLIGEOMETRYTOMATERIALMAP_EQUAL_H_


int mliGeometryToMaterialMap_equal(
        const struct mliGeometryToMaterialMap *a,
        const struct mliGeometryToMaterialMap *b);
#endif

/* mliGeometryToMaterialMap_serialize */
/* ---------------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRYTOMATERIALMAP_SERIALIZE_H_
#define MLIGEOMETRYTOMATERIALMAP_SERIALIZE_H_


int mliGeometryToMaterialMap_malloc_fread(
        struct mliGeometryToMaterialMap *geomap,
        FILE *f);
int mliGeometryToMaterialMap_fwrite(
        const struct mliGeometryToMaterialMap *geomap,
        FILE *f);
#endif

/* mliGeometryToMaterialMap_valid */
/* ------------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRYTOMATERIALMAP_EQUAL_H_
#define MLIGEOMETRYTOMATERIALMAP_EQUAL_H_


int mliGeometryToMaterialMap_valid(
        const struct mliGeometryToMaterialMap *geomap);

int mliGeometryToMaterialMap_valid_wrt_Geometry(
        const struct mliGeometryToMaterialMap *geomap,
        const struct mliGeometry *geometry);
int mliGeometryToMaterialMap_valid_wrt_Materials(
        const struct mliGeometryToMaterialMap *geomap,
        const struct mliMaterials *materials);
#endif

/* mliGeometry_equal */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRY_EQUAL_H_
#define MLIGEOMETRY_EQUAL_H_


int mliGeometry_equal(const struct mliGeometry *a, const struct mliGeometry *b);
int mliGeometry_objects_equal(
        const struct mliGeometry *a,
        const struct mliGeometry *b);
int mliGeometry_object_references_equal(
        const struct mliGeometry *a,
        const struct mliGeometry *b);
#endif

/* mliGeometry_serialize */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRY_SERIALIZE_H_
#define MLIGEOMETRY_SERIALIZE_H_


int mliGeometry_fwrite(const struct mliGeometry *scenery, FILE *f);
int mliGeometry_malloc_fread(struct mliGeometry *scenery, FILE *f);
#endif

/* mliGeometry_valid */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRY_VALID_H_
#define MLIGEOMETRY_VALID_H_


int mliGeometry_valid(const struct mliGeometry *geometry);
int mliGeometry_valid_objects(const struct mliGeometry *geometry);
int mliGeometry_valid_robjects_HomTras(const struct mliGeometry *geometry);
int mliGeometry_valid_object_references(const struct mliGeometry *geometry);
#endif

/* mliTmpOcTree */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLITMPOCTREE_H_
#define MLITMPOCTREE_H_


#define MLI_TMPNODE_FLAT_INDEX_NONE -1

uint64_t mli_guess_octree_depth_based_on_num_objects(
        const uint64_t num_objects);

/*
 * The dynamic node
 * ================
 */

struct mliTmpNode {
        struct mliTmpNode *children[8];
        uint32_t num_objects;
        uint32_t *objects;

        int32_t flat_index;
        int32_t node_index;
        int32_t leaf_index;
};

int mliTmpNode_malloc(struct mliTmpNode *n, const uint32_t num_objects);
void mliTmpNode_free(struct mliTmpNode *n);
struct mliTmpNode mliTmpNode_init(void);
void mliTmpNode_num_nodes_leafs_objects(
        const struct mliTmpNode *root_node,
        uint64_t *num_nodes,
        uint64_t *num_leafs,
        uint64_t *num_object_links);
void mliTmpNode_num_nodes_leafs_objects_walk(
        const struct mliTmpNode *node,
        uint64_t *num_nodes,
        uint64_t *num_leafs,
        uint64_t *num_object_links);
void mliTmpNode_set_flat_index(struct mliTmpNode *root_node);
void mliTmpNode_set_flat_index_walk(
        struct mliTmpNode *node,
        int32_t *flat_index,
        int32_t *node_index,
        int32_t *leaf_index);
int mliTmpNode_exists_and_has_objects(const struct mliTmpNode *node);
void mliTmpNode_print(
        const struct mliTmpNode *node,
        const uint32_t indent,
        const uint32_t child);
int mliTmpNode_num_children(const struct mliTmpNode *node);
int mliTmpNode_malloc_tree_from_bundle(
        struct mliTmpNode *root_node,
        const void *bundle,
        const uint32_t num_items_in_bundle,
        int (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mliAABB),
        const struct mliCube bundle_cube);
int mliTmpNode_add_children(
        struct mliTmpNode *node,
        const void *bundle,
        int (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mliAABB),
        const struct mliCube cube,
        const uint64_t depth,
        const uint64_t max_depth);
uint32_t mliTmpNode_signs_to_child(
        const uint32_t sx,
        const uint32_t sy,
        const uint32_t sz);

/*
 * The dynamic octree
 * ==================
 */

struct mliTmpOcTree {
        struct mliCube cube;
        struct mliTmpNode root;
};

int mliTmpOcTree_malloc_from_bundle(
        struct mliTmpOcTree *octree,
        const void *bundle,
        const uint32_t num_items_in_bundle,
        int (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mliAABB),
        struct mliAABB bundle_aabb);

void mliTmpOcTree_free(struct mliTmpOcTree *octree);
struct mliTmpOcTree mliTmpOcTree_init(void);
void mliTmpOcTree_print(const struct mliTmpOcTree *octree);

#endif

/* mliDynPhotonInteraction */
/* ----------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIDYNPHOTONINTERACTION_H_
#define MLIDYNPHOTONINTERACTION_H_


MLIDYNARRAY_DEFINITON(mli, PhotonInteraction, struct mliPhotonInteraction)

void mliDynPhotonInteraction_print(
        const struct mliDynPhotonInteraction *history,
        const struct mliGeometry *scenery);
#endif

/* mliOcTree */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIOCTREE_H_
#define MLIOCTREE_H_


#define MLI_OCTREE_TYPE_NONE 0
#define MLI_OCTREE_TYPE_NODE 1
#define MLI_OCTREE_TYPE_LEAF 2

struct mliLeafAddress {
        uint32_t first_object_link;
        uint32_t num_object_links;
};

struct mliLeafArray {
        uint64_t num_leafs;
        struct mliLeafAddress *adresses;
        uint64_t num_object_links;
        uint32_t *object_links;
};

struct mliNode {
        uint32_t children[8];
        uint8_t types[8];
};

struct mliOcTree {
        struct mliCube cube;
        uint64_t num_nodes;
        struct mliNode *nodes;
        struct mliLeafArray leafs;
        uint8_t root_type;
};

void mliOcTree_print(const struct mliOcTree *tree);
void mliOcTree_print_walk(
        const struct mliOcTree *tree,
        const int32_t node_idx,
        const uint8_t node_type,
        const uint32_t indent,
        const uint32_t child);
int mliOcTree_equal_payload(
        const struct mliOcTree *tree,
        const struct mliTmpOcTree *tmp_octree);
int mliOcTree_equal_payload_walk(
        const struct mliOcTree *tree,
        const int32_t node_idx,
        const int32_t node_type,
        const struct mliTmpNode *tmp_node);
uint32_t mliOcTree_leaf_object_link(
        const struct mliOcTree *tree,
        const uint64_t leaf,
        const uint64_t object_link);
uint64_t mliOcTree_leaf_num_objects(
        const struct mliOcTree *tree,
        const uint64_t leaf);
uint64_t mliOcTree_node_num_children(
        const struct mliOcTree *tree,
        const uint64_t node_idx);
void mliOcTree_set(struct mliOcTree *tree, const struct mliTmpOcTree *dyntree);
void mliOcTree_set_walk(
        struct mliOcTree *tree,
        const struct mliTmpNode *dynnode,
        uint64_t *object_link_size);
void mliOcTree_set_leaf(
        struct mliOcTree *tree,
        const struct mliTmpNode *dynnode,
        uint64_t *object_link_size);
void mliOcTree_set_node(
        struct mliOcTree *tree,
        const struct mliTmpNode *dynnode);
int mliOcTree_malloc(
        struct mliOcTree *tree,
        const uint64_t num_nodes,
        const uint64_t num_leafs,
        const uint64_t num_object_links);
void mliOcTree_free(struct mliOcTree *tree);
struct mliOcTree mliOcTree_init(void);
struct mliNode mliNode_init(void);
int mliLeafArray_malloc(
        struct mliLeafArray *leafs,
        const uint64_t num_leafs,
        const uint64_t num_object_links);
void mliLeafArray_free(struct mliLeafArray *leafs);
struct mliLeafArray mliLeafArray_init(void);
struct mliLeafAddress mliLeafAddress_init(void);

int mliOcTree_malloc_from_object_wavefront(
        struct mliOcTree *octree,
        const struct mliObject *object);

struct mliGeometryAndAccelerator;
int mliOcTree_malloc_from_Geometry(
        struct mliOcTree *octree,
        const struct mliGeometryAndAccelerator *accgeo,
        const struct mliAABB outermost_aabb);

#endif

/* mliOcTree_equal */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIOCTREE_EQUAL_H_
#define MLIOCTREE_EQUAL_H_


int mliOcTree_equal(const struct mliOcTree *a, const struct mliOcTree *b);
#endif

/* mliOcTree_serialize */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIOCTREE_SERIALIZE_H_
#define MLIOCTREE_SERIALIZE_H_


int mliOcTree_fwrite(const struct mliOcTree *octree, FILE *f);
int mliOcTree_malloc_fread(struct mliOcTree *octree, FILE *f);

int mliOcTree_write_to_path(const struct mliOcTree *octree, const char *path);
int mliOcTree_malloc_from_path(struct mliOcTree *octree, const char *path);
#endif

/* mliOcTree_valid */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIOCTREE_VALID_H_
#define MLIOCTREE_VALID_H_


int mliOcTree_valid(const struct mliOcTree *octree);
int mliOcTree_valid_wrt_links(
        const struct mliOcTree *octree,
        const uint32_t num_links);
#endif

/* mli_ray_octree_traversal */
/* ------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RAY_OCTREE_TRAVERSAL_H_
#define MLI_RAY_OCTREE_TRAVERSAL_H_



#define MLI_RAY_OCTREE_TRAVERSAL_EPSILON 1.0e-307

void mli_ray_octree_traversal(
        const struct mliOcTree *octree,
        const struct mliRay ray,
        void *work,
        void (*work_on_leaf_node)(
                void *,
                const struct mliOcTree *,
                const uint32_t));

void mli_ray_octree_traversal_sub(
        struct mliVec t0,
        struct mliVec t1,
        const struct mliOcTree *octree,
        const int32_t node_idx,
        const int32_t node_type,
        uint8_t permutation,
        void *work,
        void (*work_on_leaf_node)(
                void *,
                const struct mliOcTree *,
                const uint32_t));

int mli_ray_octree_traversal_next_octree_node(
        const struct mliVec tm,
        int x,
        int y,
        int z);

int mli_ray_octree_traversal_first_octree_node(
        const struct mliVec t0,
        const struct mliVec tm);

#endif

/* mliAccelerator */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIACCELERATOR_H_
#define MLIACCELERATOR_H_


struct mliAccelerator {
        uint32_t num_objects;
        struct mliOcTree *object_octrees;

        uint32_t num_robjects;
        struct mliAABB *robject_aabbs;

        struct mliOcTree scenery_octree;
};

struct mliAccelerator mliAccelerator_init(void);

void mliAccelerator_free(struct mliAccelerator *accel);

int mliAccelerator_malloc(
        struct mliAccelerator *accel,
        const uint32_t num_objects,
        const uint32_t num_robjects);

int mliAccelerator_malloc_from_Geometry(
        struct mliAccelerator *accel,
        const struct mliGeometry *geometry);

int mliAccelerator_set_robject_aabbs(
        struct mliAccelerator *accel,
        const struct mliGeometry *geometry);

int mliAccelerator_set_object_octrees(
        struct mliAccelerator *accel,
        const struct mliGeometry *geometry);

void mliAccelerator_info_fprint(FILE *f, const struct mliAccelerator *accel);

struct mliAABB mliAccelerator_outermost_aabb(
        const struct mliAccelerator *accel);

#endif

/* mliAccelerator_equal */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIACCELERATOR_EQUAL_H_
#define MLIACCELERATOR_EQUAL_H_


int mliAccelerator_equal(
        const struct mliAccelerator *a,
        const struct mliAccelerator *b);
#endif

/* mliAccelerator_serialize */
/* ------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIACCELERATOR_SERIALIZE_H_
#define MLIACCELERATOR_SERIALIZE_H_


int mliAccelerator_malloc_fread(struct mliAccelerator *accel, FILE *f);
int mliAccelerator_fwrite(const struct mliAccelerator *accel, FILE *f);
#endif

/* mliAccelerator_valid */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIACCELERATOR_VALID_H_
#define MLIACCELERATOR_VALID_H_


int mliAccelerator_valid(const struct mliAccelerator *accel);
int mliAccelerator_valid_wrt_Geometry(
        const struct mliAccelerator *accel,
        const struct mliGeometry *geometry);
#endif

/* mliGeometryAndAccelerator */
/* ------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRYANDACCELERATOR_H_
#define MLIGEOMETRYANDACCELERATOR_H_


struct mliGeometryAndAccelerator {
        /*
         * A temporary container to allow access to both geometry and its
         * accelerator using only one pointer.
         */
        const struct mliGeometry *geometry;
        const struct mliAccelerator *accelerator;
};

#endif

/* mliGeometry_AABB */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRY_AABB_H_
#define MLIGEOMETRY_AABB_H_


int mliGeometry_robject_has_overlap_aabb_void(
        const void *accgeo,
        const uint32_t robject_idx,
        const struct mliAABB aabb);

int mliGeometry_robject_has_overlap_aabb(
        const struct mliGeometryAndAccelerator *accgeo,
        const uint32_t robject_idx,
        const struct mliAABB aabb);

#endif

/* mliScenery */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLISCENERY_H_
#define MLISCENERY_H_


struct mliScenery {
        struct mliGeometry geometry;
        struct mliAccelerator accelerator;
        struct mliMaterials materials;
        struct mliGeometryToMaterialMap geomap;
};

struct mliScenery mliScenery_init(void);
void mliScenery_free(struct mliScenery *scenery);
void mliScenery_info_fprint(FILE *f, const struct mliScenery *scenery);
#endif

/* mliScenery_equal */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLISCENERY_EQUAL_H_
#define MLISCENERY_EQUAL_H_


int mliScenery_equal(const struct mliScenery *a, const struct mliScenery *b);
#endif

/* mliScenery_minimal_object */
/* ------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLISCENERY_MINIMAL_OBJECT_H_
#define MLISCENERY_MINIMAL_OBJECT_H_


int mliScenery_malloc_minimal_from_wavefront(
        struct mliScenery *scenery,
        const char *path);
#endif

/* mliScenery_serialize */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLISCENERY_SERIALIZE_H_
#define MLISCENERY_SERIALIZE_H_


int mliScenery_fwrite(const struct mliScenery *scenery, FILE *f);
int mliScenery_malloc_fread(struct mliScenery *scenery, FILE *f);

int mliScenery_malloc_from_path(struct mliScenery *scenery, const char *path);
int mliScenery_write_to_path(
        const struct mliScenery *scenery,
        const char *path);
#endif

/* mliScenery_tar */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLISCENERY_TAR_H_
#define MLISCENERY_TAR_H_


int mliScenery_malloc_fread_tar(struct mliScenery *scenery, FILE *f);
int mliScenery_malloc_from_path_tar(
        struct mliScenery *scenery,
        const char *path);
int mliScenery_malloc_from_Archive(
        struct mliScenery *scenery,
        const struct mliArchive *archive);
#endif

/* mliScenery_valid */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLISCENERY_VALID_H_
#define MLISCENERY_VALID_H_


int mliScenery_valid(const struct mliScenery *scenery);
#endif

/* mli_intersection_and_scenery */
/* ---------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_INTERSECTION_AND_SCENERY_H_
#define MLI_INTERSECTION_AND_SCENERY_H_


const struct mliFunc *mli_get_refractive_index_coming_from(
        const struct mliScenery *scenery,
        const struct mliIntersectionSurfaceNormal *isec);
const struct mliFunc *mli_get_refractive_index_going_to(
        const struct mliScenery *scenery,
        const struct mliIntersectionSurfaceNormal *isec);
struct mliSide mli_get_side_going_to(
        const struct mliScenery *scenery,
        const struct mliIntersectionSurfaceNormal *isec);
struct mliSide mli_get_side_coming_from(
        const struct mliScenery *scenery,
        const struct mliIntersectionSurfaceNormal *isec);
#endif

/* mli_photon_propagation */
/* ---------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PHOTON_PROPAGATION_H_
#define MLI_PHOTON_PROPAGATION_H_



struct mliEnv {
        const struct mliScenery *scenery;
        struct mliDynPhotonInteraction *history;
        struct mliPhoton *photon;
        struct mliPrng *prng;
        uint64_t max_interactions;
};

int mli_propagate_photon(
        const struct mliScenery *scenery,
        struct mliDynPhotonInteraction *history,
        struct mliPhoton *photon,
        struct mliPrng *prng,
        const uint64_t max_interactions);
int mli_propagate_photon_work_on_causal_intersection(struct mliEnv *env);
int mli_propagate_photon_distance_until_absorbtion(
        const struct mliFunc *absorbtion_in_medium_passing_through,
        const double wavelength,
        struct mliPrng *prng,
        double *distance_until_absorbtion);
int mli_propagate_photon_interact_with_object(
        struct mliEnv *env,
        const struct mliIntersectionSurfaceNormal *isec);
int mli_propagate_photon_fresnel_refraction_and_reflection(
        struct mliEnv *env,
        const struct mliIntersectionSurfaceNormal *isec);
int mli_propagate_photon_probability_passing_medium_coming_from(
        const struct mliScenery *scenery,
        const struct mliPhoton *photon,
        const struct mliIntersectionSurfaceNormal *isec,
        double *probability_passing);
int mli_propagate_photon_pass_boundary_layer(
        struct mliEnv *env,
        const struct mliIntersectionSurfaceNormal *isec,
        const struct mliFresnel fresnel);
int mli_propagate_photon_phong(
        struct mliEnv *env,
        const struct mliIntersectionSurfaceNormal *isec);
struct mliPhotonInteraction mliPhotonInteraction_from_Intersection(
        const int64_t type,
        const struct mliScenery *scenery,
        const struct mliIntersectionSurfaceNormal *isec);
int mli_propagate_photon_env(struct mliEnv *env);
#endif

/* mli_ray_scenery_query */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RAY_SCENERY_QUERY_H_
#define MLI_RAY_SCENERY_QUERY_H_


int mli_query_intersection(
        const struct mliScenery *scenery,
        const struct mliRay ray_root,
        struct mliIntersection *isec);

int mli_query_intersection_with_surface_normal(
        const struct mliScenery *scenery,
        const struct mliRay ray_root,
        struct mliIntersectionSurfaceNormal *isecsrf);

int mli_query_object_reference(
        const struct mliObject *object,
        const struct mliOcTree *object_octree,
        const struct mliHomTraComp local2root_comp,
        const struct mliRay ray_root,
        struct mliIntersection *isec);

struct mliQueryInnerWork {
        struct mliIntersection *intersection;
        const struct mliObject *object;
        struct mliRay ray_object;
        int has_intersection;
};

struct mliQueryOuterWork {
        struct mliIntersection *intersection;
        const struct mliGeometry *geometry;
        const struct mliAccelerator *accelerator;
        struct mliRay ray_root;
};

void mli_outer_scenery_traversal(
        void *_outer,
        const struct mliOcTree *scenery_octree,
        const uint32_t scenery_octree_leaf_idx);

void mli_inner_object_traversal(
        void *_inner,
        const struct mliOcTree *object_octree,
        const uint32_t object_octree_leaf_idx);

#endif

/* mli_viewer_viewer */
/* ----------------- */

/* Copyright 2019 Sebastian Achim Mueller                                     */
#ifndef MLI_VIEWER_VIEWER_H_
#define MLI_VIEWER_VIEWER_H_


#define MLIVR_ESCAPE_KEY 27
#define MLIVR_SPACE_KEY 32

void mlivr_clear_screen(void);

void mlivr_print_help(void);

void mlivr_print_info_line(
        const struct mliView view,
        const struct mlivrCursor cursor,
        const struct mliTracerConfig tracer_config);

void mlivr_timestamp_now_19chars(char *buffer);

int mlivr_export_image(
        const struct mliScenery *scenery,
        const struct mlivrConfig config,
        const struct mliView view,
        struct mliPrng *prng,
        const struct mliTracerConfig *tracer_config,
        const double object_distance,
        const char *path);

int mlivr_run_interactive_viewer(
        const struct mliScenery *scenery,
        const struct mlivrConfig config);
int mlivr_run_interactive_viewer_try_non_canonical_stdin(
        const struct mliScenery *scenery,
        const struct mlivrConfig config);
#endif

