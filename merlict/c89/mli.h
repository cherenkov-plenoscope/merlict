#include <ctype.h>
#include <errno.h>
#include <inttypes.h>
#include <math.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <termios.h>

/* Cursor */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VIEWER_CURSOR_H_
#define MLI_VIEWER_CURSOR_H_


struct mli_viewer_Cursor {
        int64_t active;
        uint64_t col;
        uint64_t row;
        uint64_t num_cols;
        uint64_t num_rows;
};

void mli_viewer_Cursor_move_up(struct mli_viewer_Cursor *cursor);
void mli_viewer_Cursor_move_down(struct mli_viewer_Cursor *cursor);
void mli_viewer_Cursor_move_right(struct mli_viewer_Cursor *cursor);
void mli_viewer_Cursor_move_left(struct mli_viewer_Cursor *cursor);

#endif

/* EventIo_Header */
/* -------------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */
#ifndef MLI_CORSIKA_EVENTIOHEADER_H_
#define MLI_CORSIKA_EVENTIOHEADER_H_


enum mli_eventio_level { MLI_EVENTIO_TOP_LEVEL = 1, MLI_EVENTIO_SUB_LEVEL = 0 };

struct mliEventIoHeader {
        int is_sync;
        int32_t type;
        int32_t version;
        int user;
        int extended;
        int only_sub_objects;
        uint64_t length;
        int32_t id;
};
struct mliEventIoHeader mliEventIoHeader_init(void);
int mliEventIoHeader_read(struct mliEventIoHeader *header, FILE *f, int level);
void mliEventIoHeader_fprint(const struct mliEventIoHeader head, FILE *f);
#endif

/* avl_Tree */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_AVL_TREE_H_
#define MLI_AVL_TREE_H_


enum mli_avltree_depth_changes {
        MLI_AVL_DEPTH_GREW_BY_ONE = 1,
        MLI_AVL_DEPTH_DID_NOT_CHANGE = 0,
        MLI_AVL_DEPTH_SHRUNK_BY_ONE = -1
};

struct mli_Avl {
        struct mli_Avl *left;
        struct mli_Avl *right;
        int64_t balance;
};

struct mli_AvlTree {
        struct mli_Avl *root;
        int64_t (*compare)(const void *a, const void *b);
};

int mli_AvlTree_insert(struct mli_AvlTree *t, struct mli_Avl *a);
int mli_AvlTree_remove(struct mli_AvlTree *t, struct mli_Avl *a);
int mli_AvlTree_removeroot(struct mli_AvlTree *t);
struct mli_Avl *mli_AvlTree_find(
        struct mli_AvlTree *t,
        const struct mli_Avl *probe);

struct mli_AvlNode {
        struct mli_Avl avl;
        int64_t key;
        int64_t value;
};

struct mli_AvlNode mli_AvlNode_init(void);
int64_t mli_AvlNode_compare(const void *a, const void *b);
void mli_AvlNode_print(struct mli_Avl *a, int m);

#endif

/* bool */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_BOOL_BOOL_H_
#define MLI_BOOL_BOOL_H_

#define mli_bool int

enum mli_bool_states { MLI_FALSE = 0, MLI_TRUE = 1 };

char mli_bool_to_char(const mli_bool self);

#endif

/* chk */
/* --- */

/* Copyright 2018-2021 Sebastian Achim Mueller */
#ifndef CHK_DEBUG_H_
#define CHK_DEBUG_H_


/*
 *  Based on Zed Shawn's awesome Debug Macros from his book:
 *  Learn C the hard way
 */

enum chk_rc_states { CHK_FAIL = 0, CHK_SUCCESS = 1 };

#define chk_rc int

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
                goto chk_error;                                                \
        }

#define chk_msgf(C, MSGFMT)                                                    \
        if (!(C)) {                                                            \
                chk_eprint_head();                                             \
                chk_eprintf MSGFMT;                                            \
                chk_eprintf("\n");                                             \
                errno = 0;                                                     \
                goto chk_error;                                                \
        }

#define chk_bad(MSG)                                                           \
        {                                                                      \
                chk_eprint_line(MSG);                                          \
                errno = 0;                                                     \
                goto chk_error;                                                \
        }

#define chk_badf(MSGFMT)                                                       \
        {                                                                      \
                chk_eprint_head();                                             \
                chk_eprintf MSGFMT;                                            \
                chk_eprintf("\n");                                             \
                errno = 0;                                                     \
                goto chk_error;                                                \
        }

#define chk(C) chk_msg(C, "Not expected.")

#define chk_mem(C) chk_msg((C), "Out of memory.")

#define chk_malloc(PTR, TYPE, NUM)                                             \
        {                                                                      \
                PTR = (TYPE *)malloc(NUM * sizeof(TYPE));                      \
                chk_mem(PTR);                                                  \
        }

#define chk_to_io(PTR, SIZE_OF_TYPE, NUM, F)                                   \
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

#define chk_dbg                                                                \
        {                                                                      \
                /*fprintf(stderr, "%s, %d\n", __FILE__, __LINE__);*/           \
        }

#define chk_warning(MSG)                                                       \
        {                                                                      \
                chk_eprintf("[WARNING] (%s:%d) ", __FILE__, __LINE__);         \
                chk_eprintf("%s", MSG);                                        \
                chk_eprintf("\n");                                             \
        }

#endif

/* color */
/* ----- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_COLOR_H_
#define MLI_COLOR_H_


struct mli_Color {
        float r;
        float g;
        float b;
};

mli_bool mli_Color_equal(const struct mli_Color a, const struct mli_Color b);
mli_bool mli_Color_equal_margin(
        const struct mli_Color a,
        const struct mli_Color b,
        const double epsilon);
struct mli_Color mli_Color_truncate(
        const struct mli_Color color,
        const float start,
        const float stop);
struct mli_Color mli_Color_mean(
        const struct mli_Color colors[],
        const uint32_t num_colors);

struct mli_Color mli_Color_set(const float r, const float g, const float b);
mli_bool mli_Color_is_in_range(
        const struct mli_Color c,
        const float start,
        const float stop);
float mli_Color_luminance(const struct mli_Color self);

struct mli_Color mli_Color_add(
        const struct mli_Color a,
        const struct mli_Color b);
struct mli_Color mli_Color_multiply(const struct mli_Color c, const double f);
#endif

/* color_spectrum */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_COLOR_SPECTRUM_H_
#define MLI_COLOR_SPECTRUM_H_

struct mli_Func;

enum mli_ColorSpectrum_sizes { MLI_COLORSPECTRUM_SIZE = 30 };

double MLI_COLORSPECTRUM_WAVELENGTH_START(void);

double MLI_COLORSPECTRUM_WAVELENGTH_STOP(void);

struct mli_ColorSpectrumBinEdges {
        float values[MLI_COLORSPECTRUM_SIZE + 1];
};

struct mli_ColorSpectrumBinEdges mli_ColorSpectrumBinEdges_init(void);

struct mli_ColorSpectrum {
        float values[MLI_COLORSPECTRUM_SIZE];
};

struct mli_ColorSpectrum mli_ColorSpectrum_init_zeros(void);

struct mli_ColorSpectrum mli_ColorSpectrum_add(
        const struct mli_ColorSpectrum a,
        const struct mli_ColorSpectrum b);

struct mli_ColorSpectrum mli_ColorSpectrum_exp(
        const struct mli_ColorSpectrum a,
        const double factor);

struct mli_ColorSpectrum mli_ColorSpectrum_multiply(
        const struct mli_ColorSpectrum a,
        const struct mli_ColorSpectrum b);

struct mli_ColorSpectrum mli_ColorSpectrum_multiply_scalar(
        const struct mli_ColorSpectrum a,
        const double factor);

double mli_ColorSpectrum_multiply_and_sum(
        const struct mli_ColorSpectrum *a,
        const struct mli_ColorSpectrum *b);

float mli_ColorSpectrum_sum(const struct mli_ColorSpectrum *self);

void mli_ColorSpectrum_set(struct mli_ColorSpectrum *self, const float value);

void mli_ColorSpectrum_set_radiance_of_black_body_W_per_m2_per_sr(
        struct mli_ColorSpectrum *self,
        const double temperature);

chk_rc mli_ColorSpectrum_from_func(
        struct mli_ColorSpectrum *self,
        const struct mli_ColorSpectrumBinEdges *wavelength_bin_edges,
        const struct mli_Func *func);

#endif

/* corsika_version */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_CORSIKA_VERSION_H_
#define MLI_CORSIKA_VERSION_H_

#define MLI_CORSIKA_VERSION_MAYOR 0
#define MLI_CORSIKA_VERSION_MINOR 4
#define MLI_CORSIKA_VERSION_PATCH 1

#endif

/* cstr */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MTL_CSTR_H_
#define MTL_CSTR_H_


mli_bool mli_cstr_ends_with(const char *str, const char *sufix);
mli_bool mli_cstr_starts_with(const char *str, const char *prefix);

mli_bool mli_cstr_is_CRLF(const char *s);
mli_bool mli_cstr_is_CR(const char *s);
mli_bool mli_cstr_assert_only_NUL_LF_TAB_controls(const char *str);
mli_bool mli_cstr_assert_only_NUL_LF_TAB_controls__dbg(
        const char *str,
        const int dbg);

mli_bool mli_cstr_match_templeate(
        const char *s,
        const char *t,
        const char digit_wildcard);

#endif

/* cstr_numbers */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MTL_CSTR_NUMBERS_H_
#define MTL_CSTR_NUMBERS_H_


chk_rc mli_cstr_nto_int64(
        int64_t *out,
        const char *s,
        const uint64_t base,
        const uint64_t length);
chk_rc mli_cstr_to_int64(int64_t *out, const char *s, const uint64_t base);

chk_rc mli_cstr_nto_uint64(
        uint64_t *out,
        const char *s,
        const uint64_t base,
        const uint64_t length);
chk_rc mli_cstr_to_uint64(uint64_t *out, const char *s, const uint64_t base);

chk_rc mli_cstr_nto_double(double *out, const char *s, const uint64_t length);
chk_rc mli_cstr_to_double(double *out, const char *s);

chk_rc mli_cstr_print_uint64(
        uint64_t u,
        char *s,
        const uint64_t max_num_chars,
        const uint64_t base,
        const uint64_t min_num_digits);

#endif

/* func */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FUNC_H_
#define MLI_FUNC_H_


struct mli_Func {
        uint64_t num_points;
        double *x;
        double *y;
};

mli_bool mli_Func_equal(const struct mli_Func a, const struct mli_Func b);

chk_rc mli_Func_evaluate(
        const struct mli_Func *f,
        const double xarg,
        double *out);
mli_bool mli_Func_in_range(const struct mli_Func *f, const double xarg);
mli_bool mli_Func_x_is_strictly_increasing(const struct mli_Func *f);
chk_rc mli_Func_malloc(struct mli_Func *f, const uint64_t num_points);
void mli_Func_free(struct mli_Func *f);
struct mli_Func mli_Func_init(void);
mli_bool mli_Func_is_valid(const struct mli_Func *func);

mli_bool mli_Func_malloc_constant(
        struct mli_Func *self,
        const double start,
        const double stop,
        const double value);
#endif

/* func_fprint */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FUNC_plot_H_
#define MLI_FUNC_plot_H_


struct mli_Func_fprint_Config {
        double x_start;
        double x_stop;
        int x_num;
        double y_start;
        double y_stop;
        int y_num;
};

chk_rc mli_Func_fprint(
        FILE *f,
        const struct mli_Func *func,
        struct mli_Func_fprint_Config plot);

#endif

/* geometry_id */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRYID_H_
#define MLI_GEOMETRYID_H_


struct mli_GeometryId {
        uint32_t robj;
        uint32_t face;
};

struct mli_GeometryId mli_GeometryId_init(void);
#endif

/* geometry_set_from_frame */
/* ----------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef GEOMETRY_SET_FROM_FRAME_H_
#define GEOMETRY_SET_FROM_FRAME_H_

struct mli_Frame;
struct mli_Geometry;
struct mli_GeometryToMaterialMap;

chk_rc mli_Geometry_set_robjects_and_material_map_from_frame(
        const struct mli_Frame *frame,
        struct mli_Geometry *geometry,
        struct mli_GeometryToMaterialMap *geomap);

chk_rc mli_Geometry__set_robjects_and_material_map_from_frame_walk(
        const struct mli_Frame *frame,
        struct mli_Geometry *geometry,
        struct mli_GeometryToMaterialMap *geomap,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers);

#endif

/* image_Pixel */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_IMAGE_PIXEL_H_
#define MLI_IMAGE_PIXEL_H_


struct mli_image_Pixel {
        uint16_t col;
        uint16_t row;
};

struct mli_image_Pixel mli_image_Pixel_set_col_row(
        const uint16_t col,
        const uint16_t row);

void mli_image_Pixel_fprint(FILE *f, const struct mli_image_Pixel *self);

#endif

/* image_chunk */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_IMAGE_CHUNK_H_
#define MLI_IMAGE_CHUNK_H_


struct mli_image_Chunk {
        uint64_t edge_size;
        struct mli_Color *array;
};
struct mli_image_Chunk mli_image_Chunk_init(void);
void mli_image_Chunk_free(struct mli_image_Chunk *self);
chk_rc mli_image_Chunk_malloc(
        struct mli_image_Chunk *self,
        const uint64_t edge_size);
void mli_image_Chunk_set(
        struct mli_image_Chunk *self,
        const uint64_t col,
        const uint64_t row,
        const struct mli_Color color);
struct mli_Color mli_image_Chunk_get(
        const struct mli_image_Chunk *self,
        const uint64_t col,
        const uint64_t row);
struct mli_Color *mli_image_Chunk_get_ptr(
        const struct mli_image_Chunk *self,
        const uint64_t col,
        const uint64_t row);
uint64_t mli_image_Chunk__idx(
        const struct mli_image_Chunk *self,
        const uint64_t col,
        const uint64_t row);

struct mli_image_ChunkGeometry {
        uint64_t num_cols;
        uint64_t num_rows;
        uint64_t chunk_edge_size;
        uint64_t num_chunks_row;
        uint64_t num_chunks_col;
};

struct mli_image_ChunkGeometry mli_image_ChunkGeometry_set(
        const uint64_t num_cols,
        const uint64_t num_rows,
        const uint64_t chunk_edge_size);

mli_bool mli_image_ChunkGeometry_equal(
        const struct mli_image_ChunkGeometry a,
        const struct mli_image_ChunkGeometry b);

#endif

/* io_memory */
/* --------- */

/* Copyright 2018-2023 Sebastian Achim Mueller */
#ifndef MLI_IOMEMORY_H_
#define MLI_IOMEMORY_H_


struct mli_IoMemory {
        /* memory */
        unsigned char *cstr;

        /* Capacity of the allocated memory */
        uint64_t capacity;

        /* Size of the payload in the allocated memory */
        uint64_t size;

        /* Position of the cursor */
        uint64_t pos;
};

struct mli_IoMemory mli_IoMemory_init(void);
chk_rc mli_IoMemory_close(struct mli_IoMemory *self);
chk_rc mli_IoMemory_open(struct mli_IoMemory *self);
size_t mli_IoMemory_write(
        const void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IoMemory *self);
size_t mli_IoMemory_read(
        void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IoMemory *self);
void mli_IoMemory_rewind(struct mli_IoMemory *self);
int64_t mli_IoMemory_tell(struct mli_IoMemory *self);
int64_t mli_IoMemory_seek(
        struct mli_IoMemory *self,
        const int64_t offset,
        const int64_t origin);
int64_t mli_IoMemory_eof(const struct mli_IoMemory *self);

/* internal */
chk_rc mli_IoMemory__malloc(struct mli_IoMemory *self);
chk_rc mli_IoMemory__malloc_capacity(
        struct mli_IoMemory *self,
        const uint64_t capacity);
chk_rc mli_IoMemory__realloc_capacity(
        struct mli_IoMemory *self,
        const uint64_t new_capacity);
chk_rc mli_IoMemory__shrink_to_fit(struct mli_IoMemory *self);
chk_rc mli_IoMemory__write_unsigned_char(
        struct mli_IoMemory *self,
        const unsigned char *c);
chk_rc mli_IoMemory__read_unsigned_char(
        struct mli_IoMemory *self,
        unsigned char *c);
chk_rc mli_IoMemory__write_cstr(struct mli_IoMemory *self, const char *cstr);
#endif

/* json_jsmn */
/* --------- */

/* Copyright (c) 2010 Serge Zaitsev
 *               2018-2020 Sebastian Achim Mueller
 *
 * MIT License
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

/* magicid */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MAGICID_H_
#define MLI_MAGICID_H_


#define MLI_MAGICID_WORD_CAPACITY 52
#define MLI_MAGICID_SIZE MLI_MAGICID_WORD_CAPACITY + 12

struct mli_MagicId {
        char word[MLI_MAGICID_WORD_CAPACITY];
        uint32_t mayor;
        uint32_t minor;
        uint32_t patch;
};

struct mli_MagicId mli_MagicId_init(void);
chk_rc mli_MagicId_set(struct mli_MagicId *magic, const char *word);
chk_rc mli_MagicId_has_word(const struct mli_MagicId *magic, const char *word);
void mli_MagicId_warn_version(const struct mli_MagicId *magic);
#endif

/* math */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MATH_H_
#define MLI_MATH_H_


#define MLI_MATH_PI 3.14159265358979323846
#define MLI_MATH_2PI 6.28318530717958623199
#define MLI_MATH_EPSILON 1e-9
#define MLI_MATH_NAN 0. / 0.
#define MLI_MATH_IS_NAN(a) ((a) != (a))
#define MLI_MATH_MIN2(a, b) (((a) < (b)) ? (a) : (b))
#define MLI_MATH_MAX2(a, b) (((a) > (b)) ? (a) : (b))
#define MLI_MATH_ROUND(num) (num - floor(num) > 0.5) ? ceil(num) : floor(num)
#define MLI_MATH_NEAR_INT(x)                                                   \
        ((x) > 0 ? (int64_t)((x) + 0.5) : (int64_t)((x)-0.5))

#define MLI_MATH_SIGN(x) ((x) == 0 ? 0 : ((x) > 0 ? 1 : -1))

#define MLI_MATH_MIN3(a, b, c)                                                 \
        ((((a) < (b)) && ((a) < (c))) ? (a) : (((b) < (c)) ? (b) : (c)))

#define MLI_MATH_MAX3(a, b, c)                                                 \
        ((((a) > (b)) && ((a) > (c))) ? (a) : (((b) > (c)) ? (b) : (c)))

#define MLI_MATH_ARRAY_SET(arr, val, num)                                      \
        do {                                                                   \
                uint64_t i;                                                    \
                for (i = 0; i < num; i++) {                                    \
                        arr[i] = val;                                          \
                }                                                              \
        } while (0)

#define MLI_MATH_ARRAY_ARGMIN(arr, num, argmin)                                \
        do {                                                                   \
                uint64_t i;                                                    \
                argmin = 0;                                                    \
                for (i = 1; i < num; i++) {                                    \
                        if (arr[i] < arr[argmin]) {                            \
                                argmin = i;                                    \
                        }                                                      \
                }                                                              \
        } while (0)

#define MLI_MATH_ARRAY_ARGMAX(arr, num, argmax)                                \
        do {                                                                   \
                uint64_t i;                                                    \
                argmax = 0;                                                    \
                for (i = 1; i < num; i++) {                                    \
                        if (arr[i] > arr[argmax]) {                            \
                                argmax = i;                                    \
                        }                                                      \
                }                                                              \
        } while (0)

#define MLI_MATH_UPPER_COMPARE(points, num_points, point_arg, return_idx)      \
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

#define MLI_MATH_NCPY(src, dst, num)                                           \
        do {                                                                   \
                uint64_t i;                                                    \
                for (i = 0; i < num; i++) {                                    \
                        dst[i] = src[i];                                       \
                }                                                              \
        } while (0)

#define MLI_MATH_IS_BIT(var, pos) ((var) & (1 << (pos)))

double mli_math_std(
        const double vals[],
        const uint64_t size,
        const double vals_mean);
double mli_math_mean(const double vals[], const uint64_t size);
void mli_math_linspace(
        const double start,
        const double stop,
        double *points,
        const uint64_t num_points);
void mli_math_histogram(
        const double *bin_edges,
        const uint64_t num_bin_edges,
        uint64_t *underflow_bin,
        uint64_t *bins,
        uint64_t *overflow_bin,
        const double point);
uint64_t mli_math_upper_compare_double(
        const double *points,
        const uint64_t num_points,
        const double point_arg);
double mli_math_square(const double a);
double mli_math_hypot(const double a, const double b);
double mli_math_deg2rad(const double angle_in_deg);
double mli_math_rad2deg(const double angle_in_rad);
double mli_math_bin_center_in_linear_space(
        const double start,
        const double stop,
        const uint64_t num_bins,
        const uint64_t bin);
double mli_math_linear_interpolate_1d(
        const double weight,
        const double start,
        const double end);
double mli_math_linear_interpolate_2d(
        const double xarg,
        const double x0,
        const double y0,
        const double x1,
        const double y1);
double mli_math_relative_ratio(const double a, const double b);

double mli_math_interpret_int64_as_double(int64_t i);
int64_t mli_math_interpret_double_as_int64(double d);
#endif

/* math_quadratic_equation */
/* ----------------------- */

/* Copyright 2019 Sebastian A. Mueller */
#ifndef MLI_MATH_QUADRATIC_EQUATION_H_
#define MLI_MATH_QUADRATIC_EQUATION_H_


chk_rc mli_math_quadratic_equation(
        const double p,
        const double q,
        double *minus_solution,
        double *plus_solution);

#endif

/* object_face */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OBJECT_FACE_H_
#define MLI_OBJECT_FACE_H_


struct mli_object_Face {
        uint32_t a;
        uint32_t b;
        uint32_t c;
};

mli_bool mli_object_Face_equal(
        const struct mli_object_Face a,
        const struct mli_object_Face b);
struct mli_object_Face mli_object_Face_set(
        const uint32_t a,
        const uint32_t b,
        const uint32_t c);
#endif

/* pathtracer_path */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PATHTRACER_PATH_H_
#define MLI_PATHTRACER_PATH_H_


struct mli_pathtracer_Path {
        double weight;
        uint64_t num_interactions;
};

struct mli_pathtracer_Path mli_pathtracer_Path_init(void);

#endif

/* physics */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PHYSICS_H_
#define MLI_PHYSICS_H_

double mli_physics_plancks_spectral_radiance_law_W_per_m2_per_sr_per_m(
        const double wavelength,
        const double temperature);

double MLI_PHYSICS_SPEED_OF_LIGHT_M_PER_S(void);
double MLI_PHYSICS_BOLTZMANN_J_PER_K(void);
double MLI_PHYSICS_PLANCK_KG_M2_PER_S(void);

#endif

/* prng_MT19937 */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PRNG_MT19937_H_
#define MLI_PRNG_MT19937_H_


struct mli_prng_MT19937 {
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

uint32_t mli_prng_MT19937_generate_uint32(struct mli_prng_MT19937 *mt);
void mli_prng_MT19937_twist(struct mli_prng_MT19937 *mt);
struct mli_prng_MT19937 mli_prng_MT19937_init(const uint32_t seed);
void mli_prng_MT19937_reinit(struct mli_prng_MT19937 *mt, const uint32_t seed);
void mli_prng_MT19937_set_constants(struct mli_prng_MT19937 *mt);
#endif

/* prng_pcg_variants_32bit_subset */
/* ------------------------------ */

/* Copyright 2014 Melissa O'Neill <oneill@pcg-random.org>
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

/*  2021 March 23, Sebastian Achim Mueller
 *
 *  Based on 'pcg_variants.h' written by Melissa O'Neill.
 *
 *  I only kept the version with the 64bit sequence state to generate
 *  32bit numbers.
 *  I dropped 'advance', and 'boundedrand'.
 *  I only kept the seeding and the generation.
 *  I split the original header-only into a source.c and a header.h.
 */

#ifndef MLI_PRNG_PCG_VARIANTS_32BIT_SUBSET_H_INCLUDED
#define MLI_PRNG_PCG_VARIANTS_32BIT_SUBSET_H_INCLUDED


struct mli_prng_pcg_state_setseq_64 {
        uint64_t state;
        uint64_t inc;
};

void mli_prng_pcg_setseq_64_srandom_r(
        struct mli_prng_pcg_state_setseq_64 *rng,
        uint64_t initstate,
        uint64_t initseq);

uint32_t mli_prng_pcg_setseq_64_xsh_rr_32_random_r(
        struct mli_prng_pcg_state_setseq_64 *rng);

#endif

/* surface_type */
/* ------------ */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_SURFACE_TYPE_H_
#define MLI_SURFACE_TYPE_H_

struct mli_String;

enum mli_surface_type {
        MLI_SURFACE_TYPE_NONE = 0,
        MLI_SURFACE_TYPE_TRANSPARENT = 1000,
        MLI_SURFACE_TYPE_COOKTORRANCE = 5000
};

chk_rc mli_Surface_type_to_string(const uint64_t type, struct mli_String *s);
chk_rc mli_Surface_type_from_string(const struct mli_String *s, uint64_t *id);

#endif

/* testing */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_TESTING_H_
#define MLI_TESTING_H_


#define CHECK_MARGIN(first, second, margin)                                    \
        do {                                                                   \
                if ((fabs((first) - (second)) > margin)) {                     \
                        printf("  In %s, line %d\n", __FILE__, __LINE__);      \
                        printf("  Expected fabs(first - second) < margin.\n"); \
                        printf("  first:  %E\n", first);                       \
                        printf("  second: %E\n", second);                      \
                        printf("  margin: %E\n", margin);                      \
                        goto test_failure;                                     \
                }                                                              \
        } while (0)

#define CHECK(test)                                                            \
        do {                                                                   \
                if (!(test)) {                                                 \
                        printf("  In %s, line %d\n", __FILE__, __LINE__);      \
                        printf("  Expected true\n");                           \
                        goto test_failure;                                     \
                }                                                              \
        } while (0)

#define CASE(msg)                                                              \
        printf("%s: %d, %s\n", __FILE__, __LINE__, msg);                       \
        fflush(stdout);

#endif

/* thin_lens */
/* --------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_THIN_LENS_H_
#define MLI_THIN_LENS_H_

double mli_thin_lens_get_object_given_focal_and_image(
        const double focal_length,
        const double image_distance);

double mli_thin_lens_get_image_given_focal_and_object(
        const double focal_length,
        const double object_distance);
#endif

/* toggle_stdin */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VIEWER_TOGGLE_STDIN_H_
#define MLI_VIEWER_TOGGLE_STDIN_H_

#if defined(__unix__) || (defined(__APPLE__) && defined(__MACH__))
#define HAVE_TERMIOS_H 1
#endif

#ifdef HAVE_TERMIOS_H
struct termios mli_viewer_non_canonical_stdin(void);
void mli_viewer_restore_stdin(struct termios *old_terminal);
#endif

#endif

/* utils */
/* ----- */

/* Copyright 2020 Sebastian A. Mueller */
#ifndef MLI_CORSIKA_UTILS_H_
#define MLI_CORSIKA_UTILS_H_

float mli_chars_to_float(const char *four_char_word);

double mli_corsika_ux_to_cx(const double ux);
double mli_corsika_vy_to_cy(const double vy);
double mli_corsika_wz_to_cz(const double wz);

double mli_corsika_cx_to_ux(const double cx);
double mli_corsika_cy_to_vy(const double cy);
double mli_corsika_cz_to_wz(const double cz);

double mli_corsika_restore_direction_z_component(
        const double x,
        const double y);

enum mli_corsika_block_sizes {
        MLI_CORSIKA_HEADER_SIZE_BYTES = (sizeof(float) * 273),
        MLI_CORSIKA_BUNCH_SIZE_BYTES = (sizeof(float) * 8)
};

enum mli_corsika_runh_types {
        MLI_CORSIKA_RUNH_RUN_NUMBER = 1,
        MLI_CORSIKA_RUNH_SLOPE_OF_ENERGY_SPECTRUM = 15,
        MLI_CORSIKA_RUNH_ENERGY_RANGE_START = 16,
        MLI_CORSIKA_RUNH_ENERGY_RANGE_STOP = 17,
        MLI_CORSIKA_RUNH_NUM_OBSERVATION_LEVELS = 4
};

enum mli_corsika_evth_types {
        MLI_CORSIKA_EVTH_EVENT_NUMBER = 1,
        MLI_CORSIKA_EVTH_RUN_NUMBER = 43,
        MLI_CORSIKA_EVTH_PARTICLE_ID = 2,
        MLI_CORSIKA_EVTH_ENERGY_GEV = 3,
        MLI_CORSIKA_EVTH_ZENITH_RAD = 10,
        MLI_CORSIKA_EVTH_AZIMUTH_RAD = 11,
        MLI_CORSIKA_EVTH_FIRST_INTERACTION_HEIGHT_CM = 6
};

#endif

/* vec */
/* --- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VEC_H_
#define MLI_VEC_H_


struct mli_Vec {
        double x;
        double y;
        double z;
};

void mli_Vec_print(const struct mli_Vec v);
uint32_t mli_Vec_octant(const struct mli_Vec a);
mli_bool mli_Vec_equal(const struct mli_Vec a, const struct mli_Vec b);
mli_bool mli_Vec_equal_margin(
        const struct mli_Vec a,
        const struct mli_Vec b,
        const double distance_margin);
struct mli_Vec mli_Vec_mirror(
        const struct mli_Vec in,
        const struct mli_Vec normal);
double mli_Vec_norm_between(const struct mli_Vec a, const struct mli_Vec b);
double mli_Vec_angle_between(const struct mli_Vec a, const struct mli_Vec b);
struct mli_Vec mli_Vec_normalized(struct mli_Vec a);
double mli_Vec_norm(const struct mli_Vec a);
struct mli_Vec mli_Vec_multiply(const struct mli_Vec v, const double a);
double mli_Vec_dot(const struct mli_Vec a, const struct mli_Vec b);
struct mli_Vec mli_Vec_cross(const struct mli_Vec a, const struct mli_Vec b);
struct mli_Vec mli_Vec_substract(
        const struct mli_Vec a,
        const struct mli_Vec b);
struct mli_Vec mli_Vec_add(const struct mli_Vec a, const struct mli_Vec b);
struct mli_Vec mli_Vec_set(const double x, const double y, const double z);
int64_t mli_Vec_sign3_bitmask(const struct mli_Vec a, const double epsilon);
struct mli_Vec mli_Vec_mean(
        const struct mli_Vec *vecs,
        const uint64_t num_vecs);
void mli_Vec_set_dim(struct mli_Vec *a, const uint64_t dim, const double val);
double mli_Vec_get_dim(const struct mli_Vec *a, const uint64_t dim);
#endif

/* vec_AABB */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VEC_AABB_H_
#define MLI_VEC_AABB_H_


mli_bool mli_Vec_overlap_aabb(
        const struct mli_Vec a,
        const struct mli_Vec aabb_lower,
        const struct mli_Vec aabb_upper);
#endif

/* vector */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VECTOR_H_
#define MLI_VECTOR_H_


#define MLI_VECTOR_DEFINITON(NAME, PAYLOAD_TYPE)                               \
                                                                               \
        struct NAME {                                                          \
                uint64_t capacity;                                             \
                uint64_t size;                                                 \
                PAYLOAD_TYPE *array;                                           \
        };                                                                     \
                                                                               \
        struct NAME NAME##_init(void);                                         \
                                                                               \
        void NAME##_free(struct NAME *self);                                   \
                                                                               \
        int NAME##_malloc(struct NAME *self, const uint64_t capacity);         \
        int NAME##_realloc(struct NAME *self, const uint64_t capacity);        \
        int NAME##_shrink_to_fit(struct NAME *self);                           \
                                                                               \
        int NAME##_push_back(struct NAME *self, PAYLOAD_TYPE item);            \
                                                                               \
        int NAME##_set(                                                        \
                struct NAME *self, const uint64_t at, PAYLOAD_TYPE item);      \
                                                                               \
        int NAME##_get(                                                        \
                const struct NAME *self,                                       \
                const uint64_t at,                                             \
                PAYLOAD_TYPE *item);                                           \
        int NAME##_copy(struct NAME *dst, const struct NAME *src);             \
        int NAME##_copyn(                                                      \
                struct NAME *dst,                                              \
                const struct NAME *src,                                        \
                const uint64_t start,                                          \
                const uint64_t length);

#define MLI_VECTOR_IMPLEMENTATION_MALLOC(NAME, PAYLOAD_TYPE)                   \
        int NAME##_malloc(struct NAME *self, const uint64_t capacity)          \
        {                                                                      \
                NAME##_free(self);                                             \
                self->capacity = MLI_MATH_MAX2(2, capacity);                   \
                self->size = 0;                                                \
                chk_malloc(self->array, PAYLOAD_TYPE, self->capacity);         \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_realloc(struct NAME *self, const uint64_t capacity)         \
        {                                                                      \
                PAYLOAD_TYPE *new_array = (PAYLOAD_TYPE *)realloc(             \
                        (void *)self->array, capacity * sizeof(PAYLOAD_TYPE)); \
                chk_mem(new_array);                                            \
                self->array = new_array;                                       \
                self->capacity = capacity;                                     \
                if (self->capacity < self->size) {                             \
                        self->size = self->capacity;                           \
                }                                                              \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }

#define MLI_VECTOR_IMPLEMENTATION_MALLOC_ZERO_TERMINATION(NAME, PAYLOAD_TYPE)  \
        int NAME##_malloc(struct NAME *self, const uint64_t capacity)          \
        {                                                                      \
                NAME##_free(self);                                             \
                self->capacity = MLI_MATH_MAX2(2, capacity);                   \
                self->size = 0;                                                \
                chk_malloc(self->array, PAYLOAD_TYPE, self->capacity + 1);     \
                memset(self->array,                                            \
                       '\0',                                                   \
                       (self->capacity + 1) * sizeof(PAYLOAD_TYPE));           \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_realloc(struct NAME *self, const uint64_t capacity)         \
        {                                                                      \
                PAYLOAD_TYPE *new_array = (PAYLOAD_TYPE *)realloc(             \
                        (void *)self->array,                                   \
                        (capacity + 1) * sizeof(PAYLOAD_TYPE));                \
                chk_mem(new_array);                                            \
                self->array = new_array;                                       \
                self->capacity = capacity;                                     \
                                                                               \
                if (self->capacity < self->size) {                             \
                        self->size = self->capacity;                           \
                } else {                                                       \
                        int64_t num_fields_after_size =                        \
                                ((int64_t)self->capacity -                     \
                                 (int64_t)self->size);                         \
                        memset(&self->array[self->size],                       \
                               '\0',                                           \
                               num_fields_after_size * sizeof(PAYLOAD_TYPE));  \
                }                                                              \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }

#define MLI_VECTOR_IMPLEMENTATION_PRIMITIVE_FREE(NAME, PAYLOAD_TYPE)           \
        void NAME##_free(struct NAME *self)                                    \
        {                                                                      \
                free(self->array);                                             \
                (*self) = NAME##_init();                                       \
        }

#define MLI_VECTOR_IMPLEMENTATION_PAYLOAD_FREE(                                \
        NAME, PAYLOAD_TYPE, PAYLOAD_FREE)                                      \
        void NAME##_free(struct NAME *self)                                    \
        {                                                                      \
                size_t i;                                                      \
                for (i = 0; i < self->size; i++) {                             \
                        PAYLOAD_FREE(&self->array[i]);                         \
                }                                                              \
                free(self->array);                                             \
                (*self) = NAME##_init();                                       \
        }

#define MLI_VECTOR_IMPLEMENTATION_BASICS(NAME, PAYLOAD_TYPE)                   \
                                                                               \
        struct NAME NAME##_init(void)                                          \
        {                                                                      \
                struct NAME out;                                               \
                out.capacity = 0u;                                             \
                out.size = 0u;                                                 \
                out.array = NULL;                                              \
                return out;                                                    \
        }                                                                      \
                                                                               \
        int NAME##_shrink_to_fit(struct NAME *self)                            \
        {                                                                      \
                return NAME##_realloc(self, self->size);                       \
        }                                                                      \
                                                                               \
        int NAME##_push_back(struct NAME *self, PAYLOAD_TYPE item)             \
        {                                                                      \
                if (self->size == self->capacity) {                            \
                        chk_msg(NAME##_realloc(self, self->capacity * 2),      \
                                "Failed to grow vector.");                     \
                }                                                              \
                                                                               \
                self->array[self->size] = item;                                \
                self->size += 1;                                               \
                                                                               \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_set(                                                        \
                struct NAME *self, const uint64_t at, PAYLOAD_TYPE item)       \
        {                                                                      \
                chk_msg(at < self->size, "Out of range.");                     \
                self->array[at] = item;                                        \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_get(                                                        \
                const struct NAME *self,                                       \
                const uint64_t at,                                             \
                PAYLOAD_TYPE *item)                                            \
        {                                                                      \
                chk_msg(at < self->size, "Out of range.");                     \
                (*item) = self->array[at];                                     \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_copy(struct NAME *dst, const struct NAME *src)              \
        {                                                                      \
                return NAME##_copyn(dst, src, 0, src->size);                   \
        }                                                                      \
                                                                               \
        int NAME##_copyn(                                                      \
                struct NAME *dst,                                              \
                const struct NAME *src,                                        \
                const uint64_t start,                                          \
                const uint64_t length)                                         \
        {                                                                      \
                chk_msg(src->array != NULL, "Expected src to be allocated");   \
                chk_msg(start + length <= src->size,                           \
                        "Expected start + length <= src->size.")               \
                        NAME##_malloc(dst, length);                            \
                memcpy(dst->array,                                             \
                       &src->array[start],                                     \
                       length * sizeof(PAYLOAD_TYPE));                         \
                dst->size = length;                                            \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }

#define MLI_VECTOR_IMPLEMENTATION(NAME, PAYLOAD_TYPE)                          \
        MLI_VECTOR_IMPLEMENTATION_MALLOC(NAME, PAYLOAD_TYPE)                   \
        MLI_VECTOR_IMPLEMENTATION_BASICS(NAME, PAYLOAD_TYPE)                   \
        MLI_VECTOR_IMPLEMENTATION_PRIMITIVE_FREE(NAME, PAYLOAD_TYPE)

#define MLI_VECTOR_IMPLEMENTATION_ZERO_TERMINATION(NAME, PAYLOAD_TYPE)         \
        MLI_VECTOR_IMPLEMENTATION_MALLOC_ZERO_TERMINATION(NAME, PAYLOAD_TYPE)  \
        MLI_VECTOR_IMPLEMENTATION_BASICS(NAME, PAYLOAD_TYPE)                   \
        MLI_VECTOR_IMPLEMENTATION_PRIMITIVE_FREE(NAME, PAYLOAD_TYPE)

#define MLI_VECTOR_IMPLEMENTATION_FREE(NAME, PAYLOAD_TYPE, PAYLOAD_FREE)       \
        MLI_VECTOR_IMPLEMENTATION_MALLOC(NAME, PAYLOAD_TYPE)                   \
        MLI_VECTOR_IMPLEMENTATION_BASICS(NAME, PAYLOAD_TYPE)                   \
        MLI_VECTOR_IMPLEMENTATION_PAYLOAD_FREE(NAME, PAYLOAD_TYPE, PAYLOAD_FREE)

#endif

/* vector_testing */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VECTOR_TESTING_H_
#define MLI_VECTOR_TESTING_H_


#define MTL_VEC_TESTING_DEFINITON(NAME, PAYLOAD_TYPE)                          \
                                                                               \
        int NAME##_test_init(struct NAME *dh);                                 \
                                                                               \
        int NAME##_test_malloc(struct NAME *dh, const uint64_t capacity);      \
                                                                               \
        int NAME##_test_free(struct NAME *dh);

#define MTL_VEC_TESTING_IMPLEMENTATION(NAME, PAYLOAD_TYPE)                     \
                                                                               \
        int NAME##_test_init(struct NAME *dh)                                  \
        {                                                                      \
                chk(dh->capacity == 0u);                                       \
                chk(dh->size == 0u);                                           \
                chk(dh->array == NULL);                                        \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_test_malloc(struct NAME *dh, const uint64_t capacity)       \
        {                                                                      \
                chk(dh->capacity >= dh->size);                                 \
                if (capacity < 2) {                                            \
                        chk(dh->capacity == 2);                                \
                } else {                                                       \
                        chk(dh->capacity == capacity);                         \
                }                                                              \
                chk(dh->array != NULL);                                        \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_test_free(struct NAME *dh) { return NAME##_test_init(dh); }

#endif

/* version */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VERSION_H_
#define MLI_VERSION_H_


#define MLI_VERSION_MAYOR 2
#define MLI_VERSION_MINOR 2
#define MLI_VERSION_PATCH 6

void mli_version_logo_fprint(FILE *f);
void mli_version_authors_and_affiliations_fprint(FILE *f);
#endif

/* view */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VIEW_H_
#define MLI_VIEW_H_


struct mli_View {
        struct mli_Vec position;
        struct mli_Vec rotation;
        double field_of_view;
};

struct mli_View mli_View_look_up_when_possible(
        const struct mli_View camin,
        const double rate);
struct mli_View mli_View_decrease_fov(
        const struct mli_View camin,
        const double rate);
struct mli_View mli_View_increase_fov(
        const struct mli_View camin,
        const double rate);
struct mli_View mli_View_look_down_when_possible(
        const struct mli_View camin,
        const double rate);
struct mli_View mli_View_look_right(
        const struct mli_View camin,
        const double rate);
struct mli_View mli_View_move_up(
        const struct mli_View camin,
        const double rate);
struct mli_View mli_View_move_right(
        const struct mli_View camin,
        const double rate);
struct mli_View mli_View_move_forward(
        const struct mli_View camin,
        const double rate);
struct mli_Vec mli_View_direction_up(const struct mli_View cam);
struct mli_Vec mli_View_direction_right(const struct mli_View cam);
struct mli_Vec mli_View_optical_axis(const struct mli_View cam);
struct mli_HomTraComp mli_View_to_HomTraComp(const struct mli_View view);
#endif

/* Config */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VIEWER_CONFIG_H_
#define MLI_VIEWER_CONFIG_H_


struct mli_viewer_Config {
        uint32_t random_seed;
        uint64_t preview_num_cols;
        uint64_t preview_num_rows;
        uint64_t export_num_cols;
        uint64_t export_num_rows;
        double step_length;
        struct mli_View view;
        double gain;
        double gamma;

        double aperture_camera_f_stop_ratio;
        double aperture_camera_image_sensor_width;
};

struct mli_viewer_Config mli_viewer_Config_default(void);

#endif

/* aabb */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_AABB_H_
#define MLI_AABB_H_


struct mli_AABB {
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
        struct mli_Vec lower;
        struct mli_Vec upper;
};

struct mli_AABB mli_AABB_set(
        const struct mli_Vec lower,
        const struct mli_Vec upper);
struct mli_Vec mli_AABB_center(const struct mli_AABB a);
struct mli_AABB mli_AABB_outermost(
        const struct mli_AABB a,
        const struct mli_AABB b);
mli_bool mli_AABB_valid(const struct mli_AABB a);
mli_bool mli_AABB_equal(const struct mli_AABB a, const struct mli_AABB b);
mli_bool mli_AABB_is_overlapping(
        const struct mli_AABB a,
        const struct mli_AABB b);
mli_bool mli_AABB_is_point_inside(
        const struct mli_AABB a,
        const struct mli_Vec point);
#endif

/* array */
/* ----- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_ARRAY_H_
#define MLI_ARRAY_H_


#define MLI_ARRAY_DEFINITON(NAME, PAYLOAD_TYPE)                                \
                                                                               \
        struct NAME {                                                          \
                uint64_t size;                                                 \
                PAYLOAD_TYPE *array;                                           \
        };                                                                     \
                                                                               \
        struct NAME NAME##_init(void);                                         \
        void NAME##_free(struct NAME *self);                                   \
        int NAME##_malloc(struct NAME *self, const uint64_t size);             \
        int NAME##_realloc(struct NAME *self, const uint64_t size);            \
        int NAME##_set(                                                        \
                struct NAME *self, const uint64_t at, PAYLOAD_TYPE item);      \
        int NAME##_get(                                                        \
                const struct NAME *self,                                       \
                const uint64_t at,                                             \
                PAYLOAD_TYPE *item);                                           \
        int NAME##_copy(struct NAME *dst, const struct NAME *src);             \
        int NAME##_copyn(                                                      \
                struct NAME *dst,                                              \
                const struct NAME *src,                                        \
                const uint64_t start,                                          \
                const uint64_t length);

#define MLI_ARRAY_IMPLEMENTATION_MALLOC(NAME, PAYLOAD_TYPE)                    \
        int NAME##_malloc(struct NAME *self, const uint64_t size)              \
        {                                                                      \
                NAME##_free(self);                                             \
                self->size = size;                                             \
                chk_malloc(self->array, PAYLOAD_TYPE, self->size);             \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_realloc(struct NAME *self, const uint64_t size)             \
        {                                                                      \
                PAYLOAD_TYPE *new_array = (PAYLOAD_TYPE *)realloc(             \
                        (void *)self->array, size * sizeof(PAYLOAD_TYPE));     \
                chk_mem(new_array);                                            \
                self->array = new_array;                                       \
                self->size = size;                                             \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }

#define MLI_ARRAY_IMPLEMENTATION_MALLOC_ZERO_TERMINATION(NAME, PAYLOAD_TYPE)   \
        int NAME##_malloc(struct NAME *self, const uint64_t size)              \
        {                                                                      \
                NAME##_free(self);                                             \
                self->size = size;                                             \
                chk_malloc(self->array, PAYLOAD_TYPE, self->size + 1);         \
                memset(self->array,                                            \
                       '\0',                                                   \
                       (self->size + 1) * sizeof(PAYLOAD_TYPE));               \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_realloc(struct NAME *self, const uint64_t size)             \
        {                                                                      \
                PAYLOAD_TYPE *new_array = (PAYLOAD_TYPE *)realloc(             \
                        (void *)self->array,                                   \
                        (size + 1) * sizeof(PAYLOAD_TYPE));                    \
                chk_mem(new_array);                                            \
                self->array = new_array;                                       \
                self->size = size;                                             \
                memset(self->array,                                            \
                       '\0',                                                   \
                       (self->size + 1) * sizeof(PAYLOAD_TYPE));               \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }

#define MLI_ARRAY_IMPLEMENTATION_PRIMITIVE_FREE(NAME, PAYLOAD_TYPE)            \
        void NAME##_free(struct NAME *self)                                    \
        {                                                                      \
                free(self->array);                                             \
                (*self) = NAME##_init();                                       \
        }

#define MLI_ARRAY_IMPLEMENTATION_PAYLOAD_FREE(                                 \
        NAME, PAYLOAD_TYPE, PAYLOAD_FREE)                                      \
        void NAME##_free(struct NAME *self)                                    \
        {                                                                      \
                size_t i;                                                      \
                for (i = 0; i < self->size; i++) {                             \
                        PAYLOAD_FREE(&self->array[i]);                         \
                }                                                              \
                free(self->array);                                             \
                (*self) = NAME##_init();                                       \
        }

#define MLI_ARRAY_IMPLEMENTATION_BASICS(NAME, PAYLOAD_TYPE)                    \
                                                                               \
        struct NAME NAME##_init(void)                                          \
        {                                                                      \
                struct NAME out;                                               \
                out.size = 0u;                                                 \
                out.array = NULL;                                              \
                return out;                                                    \
        }                                                                      \
                                                                               \
        int NAME##_set(                                                        \
                struct NAME *self, const uint64_t at, PAYLOAD_TYPE item)       \
        {                                                                      \
                chk_msg(at < self->size, "Out of range.");                     \
                self->array[at] = item;                                        \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_get(                                                        \
                const struct NAME *self,                                       \
                const uint64_t at,                                             \
                PAYLOAD_TYPE *item)                                            \
        {                                                                      \
                chk_msg(at < self->size, "Out of range.");                     \
                (*item) = self->array[at];                                     \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }                                                                      \
                                                                               \
        int NAME##_copy(struct NAME *dst, const struct NAME *src)              \
        {                                                                      \
                return NAME##_copyn(dst, src, 0, src->size);                   \
        }                                                                      \
                                                                               \
        int NAME##_copyn(                                                      \
                struct NAME *dst,                                              \
                const struct NAME *src,                                        \
                const uint64_t start,                                          \
                const uint64_t length)                                         \
        {                                                                      \
                chk_msg(src->array != NULL, "Expected src to be allocated");   \
                chk_msg(start + length <= src->size,                           \
                        "Expected start + length <= src->size.")               \
                        NAME##_malloc(dst, length);                            \
                memcpy(dst->array,                                             \
                       &src->array[start],                                     \
                       length * sizeof(PAYLOAD_TYPE));                         \
                return CHK_SUCCESS;                                            \
        chk_error:                                                             \
                return CHK_FAIL;                                               \
        }

#define MLI_ARRAY_IMPLEMENTATION(NAME, PAYLOAD_TYPE)                           \
        MLI_ARRAY_IMPLEMENTATION_MALLOC(NAME, PAYLOAD_TYPE)                    \
        MLI_ARRAY_IMPLEMENTATION_BASICS(NAME, PAYLOAD_TYPE)                    \
        MLI_ARRAY_IMPLEMENTATION_PRIMITIVE_FREE(NAME, PAYLOAD_TYPE)

#define MLI_ARRAY_IMPLEMENTATION_ZERO_TERMINATION(NAME, PAYLOAD_TYPE)          \
        MLI_ARRAY_IMPLEMENTATION_MALLOC_ZERO_TERMINATION(NAME, PAYLOAD_TYPE)   \
        MLI_ARRAY_IMPLEMENTATION_BASICS(NAME, PAYLOAD_TYPE)                    \
        MLI_ARRAY_IMPLEMENTATION_PRIMITIVE_FREE(NAME, PAYLOAD_TYPE)

#define MLI_ARRAY_IMPLEMENTATION_FREE(NAME, PAYLOAD_TYPE, PAYLOAD_FREE)        \
        MLI_ARRAY_IMPLEMENTATION_MALLOC(NAME, PAYLOAD_TYPE)                    \
        MLI_ARRAY_IMPLEMENTATION_BASICS(NAME, PAYLOAD_TYPE)                    \
        MLI_ARRAY_IMPLEMENTATION_PAYLOAD_FREE(NAME, PAYLOAD_TYPE, PAYLOAD_FREE)

#endif

/* array_dummy_testing */
/* ------------------- */

/* Copyright Sebastian Achim Mueller */
#ifndef MLI_ARRAY_DUMMY_TESTING_H_
#define MLI_ARRAY_DUMMY_TESTING_H_


MLI_ARRAY_DEFINITON(mli_ArrayTestingFloat, float)
MLI_ARRAY_DEFINITON(mli_ArrayTestingChar, char)

#endif

/* avl_Dict */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_AVL_DICT_H_
#define MLI_AVL_DICT_H_


struct mli_AvlDict {
        struct mli_AvlTree tree;
        struct mli_AvlNode *nodes;
        uint64_t capacity;
        uint64_t back;
        uint64_t len;
};

struct mli_AvlDict mli_AvlDict_init(void);
void mli_AvlDict_free(struct mli_AvlDict *dict);
chk_rc mli_AvlDict_malloc(struct mli_AvlDict *dict, const uint64_t capacity);

chk_rc mli_AvlDict_set(
        struct mli_AvlDict *dict,
        const int64_t key,
        const int64_t value);
chk_rc mli_AvlDict_pop(struct mli_AvlDict *dict, const int64_t key);
mli_bool mli_AvlDict_has(struct mli_AvlDict *dict, const int64_t key);
chk_rc mli_AvlDict_get(
        struct mli_AvlDict *dict,
        const int64_t key,
        int64_t *value);
void mli_AvlDict_reset(struct mli_AvlDict *dict);

#endif

/* color_spectrum_array */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_COLOR_SPECTRUM_ARRAY_H_
#define MLI_COLOR_SPECTRUM_ARRAY_H_


MLI_ARRAY_DEFINITON(mli_ColorSpectrumArray, struct mli_ColorSpectrum)
#endif

/* color_vector */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_COLOR_VECTOR_H_
#define MLI_COLOR_VECTOR_H_


MLI_VECTOR_DEFINITON(mli_ColorVector, struct mli_Color)

#endif

/* cube */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_CUBE_H_
#define MLI_CUBE_H_


struct mli_Cube {
        /*
         * Cubic Oriented-Bounding-Box
         * oriented w.r.t. the unit-vectors.
         */
        struct mli_Vec lower;
        double edge_length;
};

mli_bool mli_Cube_equal(const struct mli_Cube a, const struct mli_Cube b);
struct mli_Cube mli_Cube_octree_child_code(
        const struct mli_Cube cube,
        const uint8_t a);
struct mli_Cube mli_Cube_octree_child(
        const struct mli_Cube cube,
        const uint32_t sx,
        const uint32_t sy,
        const uint32_t sz);
struct mli_Cube mli_Cube_outermost_cube(const struct mli_AABB a);
struct mli_Vec mli_Cube_center(const struct mli_Cube a);
struct mli_AABB mli_Cube_to_aabb(const struct mli_Cube a);
struct mli_Vec mli_Cube_upper(const struct mli_Cube a);
#endif

/* double_array */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_DOUBLE_ARRAY_H_
#define MLI_DOUBLE_ARRAY_H_

MLI_ARRAY_DEFINITON(mli_DoubleArray, double)
#endif

/* double_vector */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_DOUBLE_VECTOR_H_
#define MLI_DOUBLE_VECTOR_H_

MLI_VECTOR_DEFINITON(mli_DoubleVector, double)
#endif

/* float_array */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FLOAT_ARRAY_H_
#define MLI_FLOAT_ARRAY_H_

MLI_ARRAY_DEFINITON(mli_FloatArray, float)
#endif

/* float_vector */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLID_FLOAT_VECTOR_H_
#define MLID_FLOAT_VECTOR_H_

MLI_VECTOR_DEFINITON(mli_FloatVector, float)
#endif

/* frame_ptr_vector */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FRAME_PTR_VECTOR_H_
#define MLI_FRAME_PTR_VECTOR_H_

struct mli_Frame;
MLI_VECTOR_DEFINITON(mli_FramePtrVector, struct mli_Frame *)
#endif

/* fresnel */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FRESNEL_H_
#define MLI_FRESNEL_H_


struct mli_Fresnel {
        struct mli_Vec incident;
        struct mli_Vec normal;
        double n_from;
        double n_to;

        double _cosI;
        double _n_from_over_n_to;
        double _sinT2;
        double _cosT;
};

struct mli_Vec mli_Fresnel_refraction_direction(
        const struct mli_Fresnel fresnel);
struct mli_Vec mli_Fresnel_reflection_direction(
        const struct mli_Fresnel fresnel);
double mli_Fresnel_reflection_propability(const struct mli_Fresnel fresnel);
struct mli_Fresnel mli_Fresnel_init(
        const struct mli_Vec incident,
        const struct mli_Vec normal,
        const double n_from,
        const double n_to);
#endif

/* from_outside_to_inside */
/* ---------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RAYTRACING_FROM_OUTSIDE_TO_INSIDE_H_
#define MLI_RAYTRACING_FROM_OUTSIDE_TO_INSIDE_H_


mli_bool mli_raytracing_from_outside_to_inside(
        const struct mli_Vec ray_direction_local,
        const struct mli_Vec surface_normal_local);
#endif

/* image_PixelVector */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_IMAGE_PIXEL_VECTOR_H_
#define MLI_IMAGE_PIXEL_VECTOR_H_

MLI_VECTOR_DEFINITON(mli_image_PixelVector, struct mli_image_Pixel)
#endif

/* image_PixelWalk */
/* --------------- */

/* Copyright 2020-2021 Sebastian Achim Mueller */
#ifndef MLI_IMAGE_PIXELWALK_H_
#define MLI_IMAGE_PIXELWALK_H_


struct mli_image_PixelWalk {
        /*
         * PixelWalk walks over the pixels of a rectangular image in a
         * cache-aware-way with respect to raytracing.
         * The goal is to bundle rays that will go to similar directions.
         * Instead of running fast along one axis of the image, and slow along
         * the other, PixelWalk spreads the walk among both axis by walking
         * small quadratic chunks of pixels.
         */
        uint32_t chunk_row;
        uint32_t sub_row;
        uint32_t chunk_col;
        uint32_t sub_col;
        uint32_t i;
};

struct mli_image_PixelWalk mli_image_PixelWalk_init(void);
struct mli_image_Pixel mli_image_PixelWalk_get_Pixel(
        const struct mli_image_PixelWalk *self,
        const struct mli_image_ChunkGeometry *chunk_geometry);
struct mli_image_Pixel mli_image_PixelWalk_get_Pixel(
        const struct mli_image_PixelWalk *self,
        const struct mli_image_ChunkGeometry *chunk_geometry);
void mli_image_PixelWalk_walk(
        struct mli_image_PixelWalk *self,
        const struct mli_image_ChunkGeometry *chunk_geometry);

struct mli_image_PixelWalk mli_image_PixelWalk_from_pixel(
        const struct mli_image_ChunkGeometry *geometry,
        struct mli_image_Pixel pixel);

void mli_image_PixelWalk_fprint(
        FILE *f,
        const struct mli_image_PixelWalk *self);

#endif

/* intersection */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_INTERSECTION_H_
#define MLI_INTERSECTION_H_


struct mli_Intersection {
        struct mli_GeometryId geometry_id;
        struct mli_Vec position_local;
        double distance_of_ray;
};

struct mli_Intersection mli_Intersection_init(void);
#endif

/* intersection_surface_normal */
/* --------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_INTERSECTION_SURFACE_NORMAL_H_
#define MLI_INTERSECTION_SURFACE_NORMAL_H_


struct mli_IntersectionSurfaceNormal {
        struct mli_GeometryId geometry_id;
        struct mli_Vec position;
        struct mli_Vec surface_normal;
        struct mli_Vec position_local;
        struct mli_Vec surface_normal_local;
        double distance_of_ray;
        int64_t from_outside_to_inside;
};

struct mli_IntersectionSurfaceNormal mli_IntersectionSurfaceNormal_init(void);

#endif

/* mat */
/* --- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MAT_H_
#define MLI_MAT_H_


struct mli_Mat {
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

void mli_Mat_set(struct mli_Mat *a, uint64_t col, uint64_t row, const double v);
double mli_Mat_get(const struct mli_Mat *a, uint64_t col, uint64_t row);
struct mli_Mat mli_Mat_unity(void);
mli_bool mli_Mat_equal_margin(
        const struct mli_Mat a,
        const struct mli_Mat b,
        const double margin);
struct mli_Mat mli_Mat_init_axis_angle(
        const struct mli_Vec axis,
        const double angle);
struct mli_Mat mli_Mat_init_tait_bryan(
        const double rx,
        const double ry,
        const double rz);
struct mli_Mat mli_Mat_init_columns(
        const struct mli_Vec c0,
        const struct mli_Vec c1,
        const struct mli_Vec c2);
struct mli_Mat mli_Mat_covariance(
        const struct mli_Vec *vecs,
        const uint64_t num_vecs,
        const struct mli_Vec vecs_mean);
struct mli_Mat mli_Mat_transpose(const struct mli_Mat m);
struct mli_Mat mli_Mat_multiply(const struct mli_Mat x, const struct mli_Mat y);
struct mli_Mat mli_Mat_minor(const struct mli_Mat x, const int d);
struct mli_Mat mli_Mat_vector_outer_product(const struct mli_Vec v);
void mli_Mat_qr_decompose(
        const struct mli_Mat m,
        struct mli_Mat *q,
        struct mli_Mat *r);
mli_bool mli_Mat_has_shurform(const struct mli_Mat m, const double margin);
void mli_Mat_find_eigenvalues(
        const struct mli_Mat a,
        double *e0,
        double *e1,
        double *e2,
        const double margin,
        const uint64_t max_num_iterations);
chk_rc mli_Mat_find_eigenvector_for_eigenvalue(
        struct mli_Mat a,
        const double eigen_value,
        struct mli_Vec *eigen_vector,
        const double tolerance);
chk_rc mli_Mat_lup_decompose(
        struct mli_Mat *A,
        int *pivots,
        const double tolerance);
void mli_Mat_lup_solve(
        const struct mli_Mat *A,
        const int *P,
        const struct mli_Vec *b,
        struct mli_Vec *x);
struct mli_Vec mli_Mat_dot_product(
        const struct mli_Mat *m,
        const struct mli_Vec v);

#endif

/* object_face_vector */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OBJECT_FACE_VECTOR_H_
#define MLI_OBJECT_FACE_VECTOR_H_

MLI_VECTOR_DEFINITON(mli_object_FaceVector, struct mli_object_Face)
#endif

/* octree_overlaps */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OCTREE_OVERLAPS_H_
#define MLI_OCTREE_OVERLAPS_H_


#define mli_octree_OverlapVector mli_Uint32Vector
#define mli_octree_OverlapVector_init mli_Uint32Vector_init
#define mli_octree_OverlapVector_malloc mli_Uint32Vector_malloc
#define mli_octree_OverlapVector_free mli_Uint32Vector_free
#define mli_octree_OverlapVector_push_back mli_Uint32Vector_push_back

#endif

/* pinhole */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_CAMERA_PINHOLE_H_
#define MLI_CAMERA_PINHOLE_H_


struct mli_PathTracer;
struct mli_Image;
struct mli_Prng;

struct mli_camera_PinHole {
        struct mli_Vec optical_axis;
        struct mli_Vec col_axis;
        struct mli_Vec row_axis;
        struct mli_Vec principal_point;
        double distance_to_principal_point;
        double row_over_column_pixel_ratio;
};

struct mli_camera_PinHole mli_camera_PinHole_set(
        const double field_of_view,
        const struct mli_Image *image,
        const double row_over_column_pixel_ratio);

void mli_camera_PinHole_render_image(
        struct mli_camera_PinHole self,
        const struct mli_HomTraComp camera2root_comp,
        const struct mli_PathTracer *pathtracer,
        struct mli_Image *image,
        struct mli_Prng *prng);

void mli_camera_PinHole_render_image_with_view(
        const struct mli_View view,
        const struct mli_PathTracer *pathtracer,
        struct mli_Image *image,
        const double row_over_column_pixel_ratio,
        struct mli_Prng *prng);

struct mli_Ray mli_camera_PinHole_ray_at_row_col(
        const struct mli_camera_PinHole *self,
        const struct mli_Image *image,
        const uint32_t row,
        const uint32_t col);

#endif

/* prng_PCG32 */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PRNG_PCG32_H_
#define MLI_PRNG_PCG32_H_


/*      Wrapping the pcg implementation by Melissa O'Neill in
 *      pcg_variants_32bit_subset.h
 */

struct mli_prng_PCG32 {
        struct mli_prng_pcg_state_setseq_64 state_setseq_64;
};

struct mli_prng_PCG32 mli_prng_PCG32_init(const uint32_t seed);
uint32_t mli_prng_PCG32_generate_uint32(struct mli_prng_PCG32 *pcg32);
void mli_prng_PCG32_reinit(struct mli_prng_PCG32 *pcg32, const uint32_t seed);

#endif

/* quaternion */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_QUATERNION_H_
#define MLI_QUATERNION_H_


struct mli_Quaternion {
        double w;
        double x;
        double y;
        double z;
};

void mli_Quaternion_print(const struct mli_Quaternion q);
struct mli_Quaternion mli_Quaternion_set_tait_bryan(
        const double rx,
        const double ry,
        const double rz);
struct mli_Mat mli_Quaternion_to_matrix(const struct mli_Quaternion quat);
struct mli_Quaternion mli_Quaternion_set_rotaxis_and_angle(
        const struct mli_Vec rot_axis,
        const double angle);
double mli_Quaternion_norm(const struct mli_Quaternion q);
double mli_Quaternion_product_complex_conjugate(const struct mli_Quaternion p);
struct mli_Quaternion mli_Quaternion_product(
        const struct mli_Quaternion p,
        const struct mli_Quaternion q);
struct mli_Quaternion mli_Quaternion_complex_conjugate(
        const struct mli_Quaternion q);
mli_bool mli_Quaternion_equal_margin(
        const struct mli_Quaternion a,
        const struct mli_Quaternion b,
        const double margin);
mli_bool mli_Quaternion_equal(
        const struct mli_Quaternion a,
        const struct mli_Quaternion b);
struct mli_Quaternion mli_Quaternion_set(
        const double w,
        const double x,
        const double y,
        const double z);
struct mli_Quaternion mli_Quaternion_set_unit_xyz(
        const double x,
        const double y,
        const double z);
#endif

/* ray */
/* --- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RAY_H_
#define MLI_RAY_H_


struct mli_Ray {
        struct mli_Vec support;
        struct mli_Vec direction;
};

struct mli_Vec mli_Ray_at(const struct mli_Ray *ray, const double t);
struct mli_Ray mli_Ray_set(
        const struct mli_Vec support,
        const struct mli_Vec direction);
chk_rc mli_Ray_sphere_intersection(
        const struct mli_Vec support,
        const struct mli_Vec direction,
        const double radius,
        double *minus_solution,
        double *plus_solution);
#endif

/* ray_AABB */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RAY_AABB_H_
#define MLI_RAY_AABB_H_


void mli_Ray_aabb_intersections(
        const struct mli_Ray ray,
        const struct mli_AABB aabb,
        double *t_near,
        double *t_far);
mli_bool mli_Ray_aabb_intersections_is_valid_given_near_and_far(
        const double t_near,
        const double t_far);
mli_bool mli_Ray_has_overlap_aabb(
        const struct mli_Ray ray,
        const struct mli_AABB aabb);

#endif

/* string */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_STRING_H_
#define MLI_STRING_H_


MLI_VECTOR_DEFINITON(mli_String, char)

chk_rc mli_String_from_vargs(
        struct mli_String *self,
        const char *format,
        va_list args);
chk_rc mli_String_from_cstr_fromat(
        struct mli_String *self,
        const char *format,
        ...);
chk_rc mli_String_from_cstr(struct mli_String *self, const char *s);

mli_bool mli_String_equal_cstr(const struct mli_String *self, const char *cstr);

mli_bool mli_String_equal(
        const struct mli_String *self,
        const struct mli_String *other);

mli_bool mli_String_ends_with(
        const struct mli_String *self,
        const struct mli_String *suffix);
mli_bool mli_String_ends_with_cstr(
        const struct mli_String *self,
        const char *cstr);

mli_bool mli_String_starts_with(
        const struct mli_String *self,
        const struct mli_String *prefix);
mli_bool mli_String_starts_with_cstr(
        const struct mli_String *self,
        const char *cstr);

mli_bool mli_String_has_prefix_suffix(
        const struct mli_String *self,
        const struct mli_String *prefix,
        const struct mli_String *suffix);

int64_t mli_String_rfind(const struct mli_String *self, const char c);
int64_t mli_String_find(const struct mli_String *self, const char c);
chk_rc mli_String_strip(const struct mli_String *src, struct mli_String *dst);
uint64_t mli_String_countn(
        const struct mli_String *self,
        const char c,
        const uint64_t num_chars_to_scan);
int64_t mli_String_compare(
        const struct mli_String *s1,
        const struct mli_String *s2);
chk_rc mli_String_convert_line_break_CRLF_CR_to_LF(
        struct mli_String *dst,
        const struct mli_String *src);

int64_t mli_String__discover_size(const struct mli_String *self);
mli_bool mli_String_valid(const struct mli_String *self, const size_t min_size);

chk_rc mli_String__find_idx_with_cstr(
        const struct mli_String *names,
        const uint64_t num_names,
        const char *key,
        uint64_t *idx);

#endif

/* string_numbers */
/* -------------- */

/* Copyright 2018-2023 Sebastian Achim Mueller */
#ifndef MTL_STRING_NUMBERS_H_
#define MTL_STRING_NUMBERS_H_

chk_rc mli_String_nto_double(
        double *out,
        const struct mli_String *str,
        const uint64_t expected_num_chars);
chk_rc mli_String_to_double(double *out, const struct mli_String *str);
chk_rc mli_String_nto_int64(
        int64_t *out,
        const struct mli_String *str,
        const uint64_t base,
        const uint64_t expected_num_chars);
chk_rc mli_String_to_int64(
        int64_t *out,
        const struct mli_String *str,
        const uint64_t base);
chk_rc mli_String_nto_uint64(
        uint64_t *out,
        const struct mli_String *str,
        const uint64_t base,
        const uint64_t expected_num_chars);
chk_rc mli_String_to_uint64(
        uint64_t *out,
        const struct mli_String *str,
        const uint64_t base);

chk_rc mli_String_print_uint64(
        const uint64_t u,
        struct mli_String *str,
        const uint64_t base,
        const uint64_t min_num_digits,
        const char leading_char);

chk_rc mli_String_to_uint32(uint32_t *out, const struct mli_String *str);

#endif

/* string_vector */
/* ------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_STRING_VECTOR_H_
#define MLI_STRING_VECTOR_H_

MLI_VECTOR_DEFINITON(mli_StringVector, struct mli_String)
#endif

/* triangle */
/* -------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_TRIANGLE_H_
#define MLI_TRIANGLE_H_


struct mli_Triangle {
        struct mli_Vec v1;
        struct mli_Vec v2;
        struct mli_Vec v3;
};

#endif

/* triangle_aabb */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_TRIANGLE_AABB_H_
#define MLI_TRIANGLE_AABB_H_


struct mli_AABB mli_Triangle_aabb(
        const struct mli_Vec a,
        const struct mli_Vec b,
        const struct mli_Vec c);
mli_bool mli_Triangle_has_overlap_aabb(
        const struct mli_Vec a,
        const struct mli_Vec b,
        const struct mli_Vec c,
        const struct mli_AABB aabb);
struct mli_Triangle mli_Triangle_set_in_norm_aabb(
        const struct mli_Vec a,
        const struct mli_Vec b,
        const struct mli_Vec c,
        const struct mli_AABB aabb);
int64_t mli_Triangle_intersects_norm_aabb(struct mli_Triangle t);
int64_t mli_Triangle_intersects_point(struct mli_Triangle t, struct mli_Vec p);
int64_t mli_triangle_aabb_check_line(
        struct mli_Vec p1,
        struct mli_Vec p2,
        int64_t outcode_diff);
int64_t mli_triangle_aabb_check_point(
        struct mli_Vec p1,
        struct mli_Vec p2,
        double alpha,
        int64_t mask);
int64_t mli_triangle_aabb_bevel_3d(struct mli_Vec p);
int64_t mli_triangle_aabb_bevel_2d(struct mli_Vec p);
int64_t mli_triangle_aabb_face_plane(struct mli_Vec p);
#endif

/* triangle_barycentric */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_TRIANGLE_BARYCENTRIC_H_
#define MLI_TRIANGLE_BARYCENTRIC_H_


struct mli_triangle_BarycentrigWeights {
        double a;
        double b;
        double c;
};

struct mli_triangle_BarycentrigWeights mli_triangle_barycentric_weights(
        const struct mli_Vec a,
        const struct mli_Vec b,
        const struct mli_Vec c,
        const struct mli_Vec t);
#endif

/* triangle_intersection */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_TRIANGLE_INTERSECTION_H_
#define MLI_TRIANGLE_INTERSECTION_H_


struct mli_Vec mli_Triangle_interpolate_surface_normal(
        const struct mli_Vec vertex_normal_a,
        const struct mli_Vec vertex_normal_b,
        const struct mli_Vec vertex_normal_c,
        const struct mli_triangle_BarycentrigWeights weights);

mli_bool mli_Ray_intersects_triangle(
        const struct mli_Ray ray,
        const struct mli_Vec vertex_a,
        const struct mli_Vec vertex_b,
        const struct mli_Vec vertex_c,
        double *intersection_ray_parameter);

struct mli_Vec mli_Triangle_surface_normal(
        const struct mli_Vec vertex_a,
        const struct mli_Vec vertex_b,
        const struct mli_Vec vertex_c,
        const struct mli_Vec vertex_normal_a,
        const struct mli_Vec vertex_normal_b,
        const struct mli_Vec vertex_normal_c,
        const struct mli_Vec intersection_position);
#endif

/* uint32_vector */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_UINT32_VECTOR_H_
#define MLI_UINT32_VECTOR_H_

MLI_VECTOR_DEFINITON(mli_Uint32Vector, uint32_t)
#endif

/* vec_vector */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VEC_VECTOR_H_
#define MLI_VEC_VECTOR_H_

MLI_VECTOR_DEFINITON(mli_VecVector, struct mli_Vec)
#endif

/* vector_dummy_testing */
/* -------------------- */

/* Copyright Sebastian Achim Mueller */
#ifndef MTL_VEC_DUMMY_TESTING_H_
#define MTL_VEC_DUMMY_TESTING_H_


struct mtlDummy {
        float r;
        float g;
        float b;
};

MLI_VECTOR_DEFINITON(mtlDynDummy, struct mtlDummy)
MLI_VECTOR_DEFINITON(mtlDynDummyPtr, struct mtlDummy *)

MTL_VEC_TESTING_DEFINITON(mtlDynDummy, struct mtlDummy)
MTL_VEC_TESTING_DEFINITON(mtlDynDummyPtr, struct mtlDummy *)

MLI_VECTOR_DEFINITON(mtl_VectorChar, char)

#endif

/* EventIo_TelescopeOffset */
/* ----------------------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */
#ifndef MLI_CORSIKA_EVENTIO_TELESCOPEOFFSET_H_
#define MLI_CORSIKA_EVENTIO_TELESCOPEOFFSET_H_


struct mliEventIoTelescopeOffset {
        float toff;
        float xoff;
        float yoff;
        float weight;
};
struct mliEventIoTelescopeOffset mliEventIoTelescopeOffset_init(void);

MLI_ARRAY_DEFINITON(
        mliDynEventIoTelescopeOffset,
        struct mliEventIoTelescopeOffset)

#endif

/* EventIo_TelescopePosition */
/* ------------------------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */
#ifndef MLI_CORSIKA_EVENTIO_TELESCOPEPOSITION_H_
#define MLI_CORSIKA_EVENTIO_TELESCOPEPOSITION_H_


struct mliEventIoTelescopePosition {
        float x;
        float y;
        float z;
        float r;
};
MLI_ARRAY_DEFINITON(
        mliDynEventIoTelescopePosition,
        struct mliEventIoTelescopePosition)

#endif

/* Histogram2d */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_CORSIKA_HISTOGRAM2D_H_
#define MLI_CORSIKA_HISTOGRAM2D_H_


struct mli_corsika_Histogram2d {
        struct mli_AvlDict dict;
};

struct mli_corsika_Histogram2dBin {
        int32_t x;
        int32_t y;
        double value;
};

MLI_VECTOR_DEFINITON(
        mliDynCorsikaHistogram2dBin,
        struct mli_corsika_Histogram2dBin)

struct mli_corsika_Histogram2d mli_corsika_Histogram2d_init(void);

chk_rc mli_corsika_Histogram2d_malloc(
        struct mli_corsika_Histogram2d *hist,
        const uint64_t capacity);

void mli_corsika_Histogram2d_free(struct mli_corsika_Histogram2d *hist);

chk_rc mli_corsika_Histogram2d_assign(
        struct mli_corsika_Histogram2d *hist,
        const int32_t x,
        const int32_t y,
        const double weight);

chk_rc mli_corsika_Histogram2d_flatten(
        const struct mli_corsika_Histogram2d *hist,
        struct mliDynCorsikaHistogram2dBin *dump);

void mli_corsika_Histogram2d_reset(struct mli_corsika_Histogram2d *hist);
uint64_t mli_corsika_Histogram2d_len(
        const struct mli_corsika_Histogram2d *hist);

#endif

/* args */
/* ---- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_ARGS_H_
#define MLI_ARGS_H_


chk_rc mli_StringVector_from_argc_argv(
        struct mli_StringVector *self,
        int argc,
        char *argv[]);

#endif

/* axis_aligned_grid */
/* ----------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_AAG_H_
#define MLI_AAG_H_


struct mli_Idx3 {
        int64_t x;
        int64_t y;
        int64_t z;
};

struct mli_Idx3 mli_Idx3_set(const int64_t x, const int64_t y, const int64_t z);

struct mli_AxisAlignedGrid {
        struct mli_AABB bounds;
        struct mli_Idx3 num_bins;
        struct mli_Vec bin_width;
};

struct mli_AxisAlignedGrid mli_AxisAlignedGrid_set(
        struct mli_AABB bounds,
        struct mli_Idx3 num_bins);

int mli_AxisAlignedGrid_find_voxel_of_first_interaction(
        const struct mli_AxisAlignedGrid *grid,
        const struct mli_Ray *ray,
        struct mli_Idx3 *bin);

#define MLI_AXISALIGNEDGRID_RAY_DOES_NOT_INTERSECT_GRID 0
#define MLI_AXISALIGNEDGRID_RAY_STARTS_INSIDE_GRID 1
#define MLI_AXISALIGNEDGRID_RAY_STARTS_OUTSIDE_GRID_BUT_INTERSECTS 2

struct mli_AxisAlignedGridTraversal {
        const struct mli_AxisAlignedGrid *grid;
        struct mli_Idx3 voxel;
        struct mli_Vec step;
        struct mli_Vec tMax;
        struct mli_Vec tDelta;
        int valid;
};

struct mli_AxisAlignedGridTraversal mli_AxisAlignedGridTraversal_start(
        const struct mli_AxisAlignedGrid *grid,
        const struct mli_Ray *ray);
int mli_AxisAlignedGridTraversal_next(
        struct mli_AxisAlignedGridTraversal *traversal);

void mli_AxisAlignedGridTraversal_fprint(
        FILE *f,
        struct mli_AxisAlignedGridTraversal *traversal);
void mli_Ray_fprint(FILE *f, struct mli_Ray *ray);
#endif

/* color_cie1931 */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLICOLOR_CIE1931_H_
#define MLICOLOR_CIE1931_H_


/* Color detection efficiency by wavelength.
 * The color observer has three color channels (x,y,z).
 * The functions x,y,z define how these channels respond to a
 * specific wavelength of light.
 */

chk_rc mli_cie1931_spectral_matching_curve_x(struct mli_Func *self);
chk_rc mli_cie1931_spectral_matching_curve_y(struct mli_Func *self);
chk_rc mli_cie1931_spectral_matching_curve_z(struct mli_Func *self);

struct mli_Mat mli_cie1931_spectral_matching_xyz_to_rgb(void);
struct mli_Mat mli_cie1931_spectral_matching_rgb_to_xyz(void);

chk_rc mli_cie1931_spectral_radiance_of_black_body_W_per_m2_per_sr_per_m(
        struct mli_Func *self,
        const double wavelength_start,
        const double wavelength_stop,
        const double temperature,
        const uint64_t num_points);

#endif

/* homtra */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_HOMTRA_H_
#define MLI_HOMTRA_H_


struct mli_HomTraComp {
        struct mli_Vec translation;
        struct mli_Quaternion rotation;
};

struct mli_HomTra {
        struct mli_Vec translation;
        struct mli_Mat rotation;
};

void mli_HomTraComp_print(const struct mli_HomTra h);
struct mli_HomTraComp mli_HomTraComp_set(
        const struct mli_Vec translation,
        const struct mli_Quaternion rotation);
struct mli_HomTraComp mli_HomTraComp_sequence(
        const struct mli_HomTraComp a,
        const struct mli_HomTraComp b);
mli_bool mli_HomTraComp_equal(
        const struct mli_HomTraComp a,
        const struct mli_HomTraComp b);
struct mli_Vec mli_HomTraComp_dir_inverse(
        const struct mli_HomTra *t,
        const struct mli_Vec in);
struct mli_Vec mli_HomTraComp_dir(
        const struct mli_HomTra *t,
        const struct mli_Vec in);
struct mli_Vec mli_HomTraComp_pos_inverse(
        const struct mli_HomTra *t,
        const struct mli_Vec in);
struct mli_Vec mli_HomTraComp_pos(
        const struct mli_HomTra *t,
        const struct mli_Vec in);
struct mli_Ray mli_HomTraComp_ray_inverse(
        const struct mli_HomTra *t,
        const struct mli_Ray in);
struct mli_Ray mli_HomTraComp_ray(
        const struct mli_HomTra *t,
        const struct mli_Ray in);
struct mli_Ray mli_transform_ray_inverse(
        const struct mli_Mat *rotation,
        const struct mli_Vec translation,
        const struct mli_Ray in);
struct mli_Ray mli_transform_ray(
        const struct mli_Mat *rotation,
        const struct mli_Vec translation,
        const struct mli_Ray in);
struct mli_Vec mli_transform_position_inverse(
        const struct mli_Mat *rotation,
        const struct mli_Vec translation,
        const struct mli_Vec in);
struct mli_Vec mli_transform_position(
        const struct mli_Mat *rotation,
        const struct mli_Vec translation,
        const struct mli_Vec in);
struct mli_Vec mli_transform_orientation_inverse(
        const struct mli_Mat *rotation,
        const struct mli_Vec in);
struct mli_Vec mli_transform_orientation(
        const struct mli_Mat *rotation,
        const struct mli_Vec in);
struct mli_HomTra mli_HomTraComp_from_compact(
        const struct mli_HomTraComp trafo);
#endif

/* image */
/* ----- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_IMAGE_H_
#define MLI_IMAGE_H_


struct mli_Image {
        struct mli_image_ChunkGeometry geometry;
        struct mli_image_Chunk **chunks;
};

struct mli_Image mli_Image_init(void);

void mli_Image_free(struct mli_Image *self);

chk_rc mli_Image_malloc(
        struct mli_Image *self,
        const uint32_t num_cols,
        const uint32_t num_rows);

chk_rc mli_Image_malloc_same_size(
        struct mli_Image *self,
        const struct mli_Image *other);

chk_rc mli_Image_copy(const struct mli_Image *src, struct mli_Image *dst);

chk_rc mli_Image__malloc(
        struct mli_Image *self,
        const uint32_t num_cols,
        const uint32_t num_rows);

uint64_t mli_Image_num_pixel(const struct mli_Image *self);
uint64_t mli_Image_num_cols(const struct mli_Image *self);
uint64_t mli_Image_num_rows(const struct mli_Image *self);

void mli_Image__set_by_PixelWalk(
        const struct mli_Image *self,
        const struct mli_image_PixelWalk walk,
        const struct mli_Color color);
struct mli_Color mli_Image__get_by_PixelWalk(
        const struct mli_Image *self,
        const struct mli_image_PixelWalk walk);
struct mli_Color *mli_Image__get_ptr_by_PixelWalk(
        const struct mli_Image *self,
        const struct mli_image_PixelWalk walk);

void mli_Image_set_by_Pixel(
        struct mli_Image *self,
        const struct mli_image_Pixel px,
        const struct mli_Color color);
struct mli_Color mli_Image_get_by_Pixel(
        const struct mli_Image *self,
        const struct mli_image_Pixel px);

void mli_Image_set_by_col_row(
        struct mli_Image *self,
        const uint32_t col,
        const uint32_t row,
        const struct mli_Color color);
struct mli_Color mli_Image_get_by_col_row(
        const struct mli_Image *self,
        const uint32_t col,
        const uint32_t row);

chk_rc mli_image_PixelVector_above_threshold(
        const struct mli_Image *to_do_image,
        const float threshold,
        struct mli_image_PixelVector *pixels);

chk_rc mli_Image_luminance_threshold_dilatation(
        const struct mli_Image *self,
        const float threshold,
        const struct mli_Color marker,
        struct mli_Image *out);

chk_rc mli_Image_sobel(const struct mli_Image *src, struct mli_Image *dst);

chk_rc mli_Image_scale_down_twice(
        const struct mli_Image *source,
        struct mli_Image *destination);

void mli_Image_set_all(
        const struct mli_Image *img,
        const struct mli_Color color);
void mli_image_PixelVector_push_back_all_from_image(
        struct mli_image_PixelVector *pixels,
        const struct mli_Image *image);

chk_rc mli_Image_fabs_difference(
        const struct mli_Image *a,
        const struct mli_Image *b,
        struct mli_Image *out);

void mli_Image_histogram(
        struct mli_Image *img,
        const double *col_bin_edges,
        const double *row_bin_edges,
        const double x,
        const double y,
        const struct mli_Color weight);

struct mli_Color mli_Image_max(const struct mli_Image *img);

void mli_Image_multiply(struct mli_Image *img, const struct mli_Color color);
void mli_Image_power(struct mli_Image *self, const struct mli_Color power);

chk_rc mli_Image_divide_pixelwise(
        const struct mli_Image *numerator,
        const struct mli_Image *denominator,
        struct mli_Image *out);
#endif

/* image_print */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_IMAGE_PRINT_H_
#define MLI_IMAGE_PRINT_H_


enum mli_Image_print_modes {
        MLI_IMAGE_PRINT_ASCII_MONOCHROME = 100,
        MLI_IMAGE_PRINT_ASCII_ESCAPE_COLOR = 101
};

void mli_Image_print(const struct mli_Image *img, const uint64_t print_mode);
void mli_Image_print_chars(
        const struct mli_Image *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols,
        const uint64_t print_mode);

/* Colored ANSI escape sequences */

void mli_Image_print_ansi_escape_chars(
        const struct mli_Image *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols);

/* Monochrome ASCII chars */

void mli_Image_print_ascii_chars(
        const struct mli_Image *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols);
#endif

/* io_file */
/* ------- */

/* Copyright 2018-2023 Sebastian Achim Mueller */
#ifndef MLI_IO_FILE_H_
#define MLI_IO_FILE_H_


struct mli_IoFile {
        FILE *cfile;
};

struct mli_IoFile mli_IoFile_init(void);
chk_rc mli_IoFile_close(struct mli_IoFile *self);
chk_rc mli_IoFile_open(
        struct mli_IoFile *self,
        const struct mli_String *filename,
        const struct mli_String *mode);
chk_rc mli_IoFile_adopt_cfile(struct mli_IoFile *self, FILE *cfile);
size_t mli_IoFile_write(
        const void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IoFile *self);
size_t mli_IoFile_read(
        void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IoFile *self);
void mli_IoFile_rewind(struct mli_IoFile *self);
int64_t mli_IoFile_tell(struct mli_IoFile *self);
int64_t mli_IoFile_seek(
        struct mli_IoFile *self,
        const int64_t offset,
        const int64_t origin);
int64_t mli_IoFile_eof(const struct mli_IoFile *self);
int64_t mli_IoFile_flush(struct mli_IoFile *self);

mli_bool mli_IoFile__cfile_is_stdin_or_stdout_stderr(
        const struct mli_IoFile *self);
#endif

/* map */
/* --- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MAP_H_
#define MLI_MAP_H_


struct mliMapItem {
        struct mli_String key;
        uint64_t value;
};

MLI_VECTOR_DEFINITON(mli_MapItemVector, struct mliMapItem)

struct mli_Map {
        struct mli_MapItemVector items;
};

struct mli_Map mli_Map_init(void);
void mli_Map_free(struct mli_Map *map);
chk_rc mli_Map_malloc(struct mli_Map *map);
uint64_t mli_Map_size(const struct mli_Map *map);
mli_bool mli_Map_has(const struct mli_Map *map, const struct mli_String *key);
chk_rc mli_Map_insert(
        struct mli_Map *map,
        const struct mli_String *key,
        uint64_t value);
chk_rc mli_Map_find(
        const struct mli_Map *map,
        const struct mli_String *key,
        uint64_t *idx);
chk_rc mli_Map_get(
        const struct mli_Map *map,
        const struct mli_String *key,
        uint64_t *value);

#endif

/* materials_names */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MATERIALS_NAMES_H_
#define MLI_MATERIALS_NAMES_H_


struct mli_materials_Names {
        struct mli_Map spectra;
        struct mli_Map media;
        struct mli_Map surfaces;
        struct mli_Map boundary_layers;
};
struct mli_materials_Names mli_materials_Names_init(void);
chk_rc mli_materials_Names_malloc(struct mli_materials_Names *namemap);
void mli_materials_Names_free(struct mli_materials_Names *namemap);

#endif

/* medium */
/* ------ */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_MEDIUM_H_
#define MLI_MEDIUM_H_

struct mli_IO;
struct mli_Materials;
struct mli_Map;

struct mli_Medium {
        struct mli_String name;
        uint64_t refraction_spectrum;
        uint64_t absorption_spectrum;
};

struct mli_Medium mli_Medium_init(void);

void mli_Medium_free(struct mli_Medium *self);

mli_bool mli_Medium_valid_wrt_materials(
        const struct mli_Medium *self,
        const struct mli_Materials *materials);

mli_bool mli_Medium_equal(
        const struct mli_Medium *a,
        const struct mli_Medium *b);

chk_rc mli_Medium_to_io(const struct mli_Medium *self, struct mli_IO *f);
chk_rc mli_Medium_from_io(struct mli_Medium *self, struct mli_IO *f);

chk_rc mli_Medium_from_json_string_and_name(
        struct mli_Medium *self,
        const struct mli_Map *spectra_names,
        const struct mli_String *json_string,
        const struct mli_String *name);

#endif

/* medium_array */
/* ------------ */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_MEDIUM_ARRAY_H_
#define MLI_MEDIUM_ARRAY_H_

MLI_ARRAY_DEFINITON(mli_MediumArray, struct mli_Medium)
#endif

/* object */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OBJECT_H_
#define MLI_OBJECT_H_


struct mli_Object {
        uint32_t num_vertices;
        struct mli_Vec *vertices;

        uint32_t num_vertex_normals;
        struct mli_Vec *vertex_normals;

        uint32_t num_faces;
        struct mli_object_Face *faces_vertices;
        struct mli_object_Face *faces_vertex_normals;
        uint16_t *faces_materials;

        uint32_t num_materials;
        struct mli_String *material_names;
};

chk_rc mli_Object_malloc(
        struct mli_Object *obj,
        const uint64_t num_vertices,
        const uint64_t num_vertex_normals,
        const uint64_t num_faces,
        const uint64_t num_materials);
void mli_Object_free(struct mli_Object *obj);
struct mli_Object mli_Object_init(void);
mli_bool mli_Object_equal(
        const struct mli_Object *a,
        const struct mli_Object *b);
uint32_t mli_Object_resolve_material_idx(
        const struct mli_Object *obj,
        const uint32_t face_idx);
#endif

/* object_AABB */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OBJECT_AABB_H_
#define MLI_OBJECT_AABB_H_


mli_bool mli_Object_has_overlap_aabb(
        const struct mli_Object *obj,
        const struct mli_HomTra local2root,
        const struct mli_AABB aabb);

struct mli_AABB mli_Object_aabb(
        const struct mli_Object *obj,
        const struct mli_HomTra local2root);

mli_bool mli_Object_face_in_local_frame_has_overlap_aabb(
        const struct mli_Object *obj,
        const uint64_t face_idx,
        const struct mli_AABB aabb);

mli_bool mli_Object_face_in_local_frame_has_overlap_aabb_void(
        const void *obj,
        const uint32_t face_idx,
        const struct mli_AABB aabb);

struct mli_AABB mli_Object_aabb_in_local_frame(const struct mli_Object *obj);

#endif

/* object_valid */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OBJECT_VALID_H_
#define MLI_OBJECT_VALID_H_


mli_bool mli_Object_is_valid(const struct mli_Object *obj);
mli_bool mli_Object_has_valid_vertices(const struct mli_Object *obj);
mli_bool mli_Object_has_valid_faces(const struct mli_Object *obj);
mli_bool mli_Object_has_valid_normals(
        const struct mli_Object *obj,
        const double epsilon);
mli_bool mli_Object_has_valid_materials(const struct mli_Object *obj);
chk_rc mli_Object_num_unused(
        const struct mli_Object *obj,
        uint32_t *num_unused_vertices,
        uint32_t *num_unused_vertex_normals,
        uint32_t *num_unused_materials);
#endif

/* path */
/* ---- */

/* Copyright Sebastian Achim Mueller */
#ifndef MLI_PATH_H_
#define MLI_PATH_H_

chk_rc mli_path_strip_this_dir(
        const struct mli_String *src,
        struct mli_String *dst);
chk_rc mli_path_basename(const struct mli_String *src, struct mli_String *dst);
chk_rc mli_path_splitext(
        const struct mli_String *src,
        struct mli_String *dst,
        struct mli_String *ext);

#endif

/* photon */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIP_HOTON_H_
#define MLIP_HOTON_H_


struct mli_Photon {
        struct mli_Ray ray;
        double wavelength;
        int64_t id;
};

#endif

/* photon_vector */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PHOTON_VECTOR_H_
#define MLI_PHOTON_VECTOR_H_

MLI_VECTOR_DEFINITON(mli_PhotonVector, struct mli_Photon)
#endif

/* prng */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PRNG_H_
#define MLI_PRNG_H_


/**
 *      mli_Prng is a transparent container to use different
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
 *      If you want to use your own, wrapp it here using mli_Prng. See below.
 */

union mli_PrngStorage {
        struct mli_prng_MT19937 mt19937;
        struct mli_prng_PCG32 pcg32;
        /* Add your own prng here */
};

struct mli_Prng {
        union mli_PrngStorage _storage;
        uint32_t (*generate_uint32)(void *);
        void (*reinit)(void *, const uint32_t);
};

uint32_t mli_Prng_generate_uint32(struct mli_Prng *prng);
void mli_Prng_reinit(struct mli_Prng *prng, const uint32_t seed);

/**
 *      Mersenne Twister 19937
 *      ----------------------
 */
struct mli_Prng mli_Prng_init_MT19937(const uint32_t seed);
uint32_t mli_Prng_MT19937_generate_uint32(void *mt);
void mli_Prng_MT19937_reinit(void *mt, const uint32_t seed);

/**
 *      PCG32
 *      -----
 */
struct mli_Prng mli_Prng_init_PCG32(const uint32_t seed);
uint32_t mli_Prng_PCG32_generate_uint32(void *pcg);
void mli_Prng_PCG32_reinit(void *pcg, const uint32_t seed);

/**
 *      Add your own prng here
 *      ----------------------
 */

/**
 *      API
 *      ---
 */
struct mli_prng_UniformRange {
        double start;
        double range;
};

struct mli_prng_ZenithRange {
        double z_min;
        double z_range;
};

double mli_prng_draw_zenith(
        const struct mli_prng_ZenithRange range,
        struct mli_Prng *prng);
struct mli_prng_ZenithRange mli_prng_ZenithRange_set(
        const double min_zenith_distance,
        const double max_zenith_distance);
double mli_prng_draw_uniform(
        const struct mli_prng_UniformRange uniform_range,
        struct mli_Prng *prng);
struct mli_prng_UniformRange mli_prng_UniformRange_set(
        double start,
        double stop);
double mli_Prng_normal_Irwin_Hall_approximation(
        struct mli_Prng *prng,
        const double mean,
        const double std);
double mli_Prng_expovariate(struct mli_Prng *prng, const double rate);
double mli_Prng_uniform(struct mli_Prng *prng);

#endif

/* vec_random */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VEC_RANDOM_H_
#define MLI_VEC_RANDOM_H_


struct mli_Vec mli_Vec_random_position_on_disc(
        const double radius,
        struct mli_Prng *prng);
struct mli_Vec mli_Vec_random_draw_direction_in_zenith_azimuth_range(
        const struct mli_prng_ZenithRange zenith,
        const struct mli_prng_UniformRange azimuth,
        struct mli_Prng *prng);
struct mli_Vec mli_Vec_random_position_inside_unit_sphere(
        struct mli_Prng *prng);
struct mli_Vec mli_Vec_random_direction_in_hemisphere(
        struct mli_Prng *prng,
        struct mli_Vec normal);
#endif

/* CorsikaPhotonBunch */
/* ------------------ */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */
#ifndef MLI_CORSIKA_CORSIKAPHOTONBUNCH_H_
#define MLI_CORSIKA_CORSIKAPHOTONBUNCH_H_


struct mli_corsika_PhotonBunch {
        /*
         * x in cm
         * y in cm
         * ux = sin(theta) * cos(phi)
         * vy = sin(theta) * sin(phi)
         * time in nanoseconds since first interaction.
         * zem
         * photons
         * wavelength is in nanometer negative if scattered ?!
         *
         * KIT-CORSIKA coordinate-system
         *                /\ z-axis                                            *
         *                |                                                    *
         *                |\ p                                                 *
         *                | \ a                                                *
         *                |  \ r                                               *
         *                |   \ t                                              *
         *                |    \ i                                             *
         *                |     \ c                                            *
         *                |      \ l                                           *
         *                |       \ e                                          *
         *                |        \                                           *
         *                |  theta  \ m                                        *
         *                |       ___\ o                                       *
         *                |___----    \ m      ___                             *
         *                |            \ e       /| y-axis (west)              *
         *                |             \ n    /                               *
         *                |              \ t /                                 *
         *                |               \/u                                  *
         *                |              / \ m                                 *
         *                |            /    \                                  *
         *                |          /       \                                 *
         *                |        /__________\                                *
         *                |      /      ___---/                                *
         *                |    /   __---    /                                  *
         *                |  /__--- phi \ /                                    *
         * _______________|/--__________/______\ x-axis (north)                *
         *               /|                    /                               *
         *             /  |                                                    *
         *           /    |                                                    *
         *         /                                                           *
         *                                                                     *
         *                                                                     *
         * Extensive Air Shower Simulation with CORSIKA, Figure 1, page 114
         * (Version 7.6400 from December 27, 2017)
         *
         * Direction-cosines:
         *      ux = sin(theta) * cos(phi)
         *      vy = sin(theta) * sin(phi)
         * The zenith-angle theta opens relative to the negative z-axis.
         * It is the momentum of the Cherenkov-photon, which is pointing
         * down towards the observation-plane.
         */
        float x_cm;
        float y_cm;
        float ux;
        float vy;
        float time_ns;
        float z_emission_cm;
        float weight_photons;
        float wavelength_nm;
};

MLI_VECTOR_DEFINITON(mliDynCorsikaPhotonBunch, struct mli_corsika_PhotonBunch)

void mli_corsika_PhotonBunch_set_from_raw(
        struct mli_corsika_PhotonBunch *bunch,
        const float *raw);
void mli_corsika_PhotonBunch_to_raw(
        const struct mli_corsika_PhotonBunch *bunch,
        float *raw);

struct mli_Photon mli_corsika_PhotonBunch_to_merlict_photon(
        const struct mli_corsika_PhotonBunch bunch,
        const double production_distance_offset,
        const int64_t id);

struct mli_Vec mli_corsika_photon_direction_of_motion(
        const struct mli_corsika_PhotonBunch bunch);

struct mli_Vec mli_corsika_photon_support_on_observation_level(
        const struct mli_corsika_PhotonBunch bunch);

double mli_corsika_photon_wavelength(
        const struct mli_corsika_PhotonBunch bunch);

double mli_corsika_photon_emission_height(
        const struct mli_corsika_PhotonBunch bunch);

double mli_corsika_photon_relative_arrival_time_on_observation_level(
        const struct mli_corsika_PhotonBunch bunch);

#endif

/* EventIo_Run */
/* ----------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */
#ifndef MLI_CORSIKA_EVENTIO_RUN_H_
#define MLI_CORSIKA_EVENTIO_RUN_H_


struct mliEventIoRun {
        FILE *f;
        struct mliEventIoHeader next_block;
        float corsika_run_header[273];
        struct mli_String corsika_input_card;
        struct mliDynEventIoTelescopePosition telescope_positions;
};
struct mliEventIoRun mliEventIoRun_init(void);
void mliEventIoRun_finalize(struct mliEventIoRun *run);
int mliEventIoRun_begin(struct mliEventIoRun *run, FILE *stream);
int mliEventIoRun_has_still_events_left(struct mliEventIoRun *run);
int mliEventIoRun_next_block(struct mliEventIoRun *run, const int level);
int mliEventIoRun_read_273float32_block(FILE *f, float *block);

#endif

/* EventIo_Telescope */
/* ----------------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */
#ifndef MLI_CORSIKA_EVENTIO_TELESCOPE_H_
#define MLI_CORSIKA_EVENTIO_TELESCOPE_H_


struct mliEventIoTelescope {
        int16_t array_id;
        int16_t telescope_id;
        struct mliEventIoTelescopeOffset offset;
        struct mliDynCorsikaPhotonBunch photon_bunches;
};
struct mliEventIoTelescope mliEventIoTelescope_init(void);
void mliEventIoTelescope_free(struct mliEventIoTelescope *telescope);

MLI_ARRAY_DEFINITON(mliDynEventIoTelescope, struct mliEventIoTelescope)

#endif

/* aperture */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_CAMERA_APERTURE_H_
#define MLI_CAMERA_APERTURE_H_


struct mli_PathTracer;
struct mli_Prng;

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

struct mli_Vec mli_camera_Aperture_pixel_center_on_image_sensor_plane(
        const double image_sensor_width_x,
        const double image_sensor_width_y,
        const double image_sensor_distance,
        const uint64_t num_pixel_x,
        const uint64_t num_pixel_y,
        const uint64_t pixel_x,
        const uint64_t pixel_y);

struct mli_Vec mli_camera_Aperture_pixel_support_on_image_sensor_plane(
        const double image_sensor_width_x,
        const double image_sensor_width_y,
        const double image_sensor_distance,
        const uint64_t num_pixel_x,
        const uint64_t num_pixel_y,
        const uint64_t pixel_x,
        const uint64_t pixel_y,
        struct mli_Prng *prng);

struct mli_Vec mli_camera_Aperture_get_object_point(
        const double focal_length,
        const struct mli_Vec pixel_support);

double mli_camera_Aperture_focal_length_given_field_of_view_and_sensor_width(
        const double field_of_view,
        const double image_sensor_width);

struct mli_Vec mli_camera_Aperture_ray_support_on_aperture(
        const double aperture_radius,
        struct mli_Prng *prng);

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
        struct mli_Prng *prng);

struct mli_camera_Aperture {
        double focal_length;
        double aperture_radius;
        double image_sensor_distance;
        double image_sensor_width_x;
        double image_sensor_width_y;
};

struct mli_camera_Aperture mli_camera_Aperture_init(void);

chk_rc mli_camera_Aperture_render_image(
        const struct mli_camera_Aperture self,
        const struct mli_HomTraComp camera2root_comp,
        const struct mli_PathTracer *pathtracer,
        struct mli_Image *image,
        struct mli_Prng *prng);

void mli_camera_Aperture_aquire_pixels(
        const struct mli_camera_Aperture self,
        const struct mli_Image *image,
        const struct mli_HomTraComp camera2root_comp,
        const struct mli_PathTracer *pathtracer,
        const struct mli_image_PixelVector *pixels_to_do,
        struct mli_ColorVector *colors_to_do,
        struct mli_Prng *prng);

void mli_camera_Aperture_assign_pixel_colors_to_sum_and_exposure_image(
        const struct mli_image_PixelVector *pixels,
        const struct mli_ColorVector *colors,
        struct mli_Image *sum_image,
        struct mli_Image *exposure_image);

#endif

/* frame */
/* ----- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FRAME_H_
#define MLI_FRAME_H_


enum mli_frame_types {
        MLI_FRAME_TYPE_FRAME = 1000,
        MLI_FRAME_TYPE_OBJECT = 1001
};

struct mli_Frame {
        uint32_t type;
        uint32_t id;
        struct mli_HomTraComp frame2mother;
        struct mli_HomTraComp frame2root;
        struct mli_Frame *mother;

        struct mli_FramePtrVector children;

        uint32_t object;
        struct mli_Uint32Vector boundary_layers;
};

void mli_Frame_set_frame2root(struct mli_Frame *f);
void mli_Frame_print(struct mli_Frame *f);
void mli_Frame_print_walk(const struct mli_Frame *f, const uint64_t indention);
chk_rc mli_frame_string_to_type(const char *s, uint64_t *type);
chk_rc mli_frame_type_to_string(const uint64_t type, char *s);
struct mli_Frame *mli_Frame_add(struct mli_Frame *mother, const uint64_t type);
chk_rc mli_Frame_set_mother_and_child(
        struct mli_Frame *mother,
        struct mli_Frame *child);
chk_rc mli_Frame_malloc(struct mli_Frame *f, const uint64_t type);
void mli_Frame_free(struct mli_Frame *f);
struct mli_Frame mli_Frame_init(void);
chk_rc mli_Frame_estimate_num_robjects_and_total_num_boundary_layers(
        const struct mli_Frame *frame,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers);
chk_rc mli_Frame_estimate_num_robjects_and_total_num_boundary_layers_walk(
        const struct mli_Frame *frame,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers);
#endif

/* geometry_from_archive */
/* --------------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRY_FROM_ARCHIVE_H_
#define MLI_GEOMETRY_FROM_ARCHIVE_H_


struct mli_Archive;
struct mli_Geometry;
chk_rc mli_Geometry_from_archive(
        struct mli_Geometry *geometry,
        struct mli_Map *object_names,
        const struct mli_Archive *archive);
#endif

/* io */
/* -- */

/* Copyright 2018-2023 Sebastian Achim Mueller */
#ifndef MLI_STREAM_H_
#define MLI_STREAM_H_


enum mli_io_type_code {
        MLI_IO_TYPE_VOID = 0,
        MLI_IO_TYPE_FILE = 10,
        MLI_IO_TYPE_MEMORY = 20
};

union mli_io_Type {
        struct mli_IoFile file;
        struct mli_IoMemory memory;
};

struct mli_IO {
        enum mli_io_type_code type;
        union mli_io_Type data;
};

struct mli_IO mli_IO_init(void);
chk_rc mli_IO_close(struct mli_IO *self);
chk_rc mli_IO_open_memory(struct mli_IO *self);
chk_rc mli_IO_open_file(
        struct mli_IO *self,
        const struct mli_String *filename,
        const struct mli_String *mode);
chk_rc mli_IO_adopt_file(struct mli_IO *self, FILE *cfile);
chk_rc mli_IO_open_file_cstr(
        struct mli_IO *self,
        const char *filename,
        const char *mode);
size_t mli_IO_write(
        const void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IO *self);
size_t mli_IO_read(
        void *ptr,
        const size_t size,
        const size_t count,
        struct mli_IO *self);
void mli_IO_rewind(struct mli_IO *self);
int64_t mli_IO_tell(struct mli_IO *self);
int64_t mli_IO_seek(
        struct mli_IO *self,
        const int64_t offset,
        const int64_t origin);
int64_t mli_IO_eof(const struct mli_IO *self);
int64_t mli_IO_flush(struct mli_IO *self);

#define chk_IO_write(PTR, SIZE_OF_TYPE, NUM, IO)                               \
        {                                                                      \
                const uint64_t num_written =                                   \
                        mli_IO_write(PTR, SIZE_OF_TYPE, NUM, IO);              \
                chk_msg(num_written == NUM, "Can not write to mli_IO.");       \
        }

#define chk_IO_read(PTR, SIZE_OF_TYPE, NUM, IO)                                \
        {                                                                      \
                const uint64_t num_read =                                      \
                        mli_IO_read(PTR, SIZE_OF_TYPE, NUM, IO);               \
                chk_msg(num_read == NUM, "Can not read from mli_IO.");         \
        }

#endif

/* io_text */
/* ------- */

/* Copyright 2018-2023 Sebastian Achim Mueller */
#ifndef MLI_IO_TEXT_H_
#define MLI_IO_TEXT_H_


int64_t mli_IO_text_getc(struct mli_IO *self);
chk_rc mli_IO_text_putc(struct mli_IO *self, const char c);

chk_rc mli_IO_text_write_cstr(struct mli_IO *self, const char *cstr);
chk_rc mli_IO_text_write_cstr_format(
        struct mli_IO *self,
        const char *format,
        ...);

chk_rc mli_IO_text_read_string(struct mli_IO *self, struct mli_String *str);
chk_rc mli_IO_text_read_line(
        struct mli_IO *self,
        struct mli_String *line,
        const char delimiter);

chk_rc mli_IO_text_write_String(
        struct mli_IO *self,
        const struct mli_String *str);

chk_rc mli_IO_text_write_multi_line_debug_view(
        struct mli_IO *self,
        const struct mli_String *text,
        const uint64_t line_number,
        const uint64_t line_radius);

#endif

/* json */
/* ---- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_JSON_H_
#define MLI_JSON_H_


struct mli_Json {
        struct mli_String raw;
        uint64_t num_tokens;
        struct jsmntok_t *tokens;
};

chk_rc mli_Json_from_string(
        struct mli_Json *self,
        const struct mli_String *str);
chk_rc mli_Json_from_io(struct mli_Json *self, struct mli_IO *io);
void mli_Json_free(struct mli_Json *self);
struct mli_Json mli_Json_init(void);

chk_rc mli_Json_debug_to_io(const struct mli_Json *self, struct mli_IO *io);
chk_rc mli_Json_debug_token_to_io(
        const struct mli_Json *self,
        const uint64_t token,
        struct mli_IO *io);
chk_rc mli_Json_debug_token_fprint(
        FILE *f,
        const struct mli_Json *self,
        const uint64_t token);

uint64_t mli_Json__token_by_index_unsafe(
        const struct mli_Json *json,
        const uint64_t start_token_idx,
        const uint64_t child_idx);
chk_rc mli_Json_token_by_key(
        const struct mli_Json *json,
        const uint64_t token,
        const char *key,
        uint64_t *key_token);
chk_rc mli_Json_token_by_idx(
        const struct mli_Json *json,
        const uint64_t token,
        const uint64_t idx,
        uint64_t *idx_token);

chk_rc mli_Json_token_by_key_eprint(
        const struct mli_Json *json,
        const uint64_t token,
        const char *key,
        uint64_t *key_token);
chk_rc mli_Json_double_by_token(
        const struct mli_Json *json,
        const uint64_t token,
        double *val);
chk_rc mli_Json_double_by_key(
        const struct mli_Json *json,
        const uint64_t token,
        double *val,
        const char *key);
chk_rc mli_Json_int64_by_token(
        const struct mli_Json *json,
        const uint64_t token,
        int64_t *return_int64);
chk_rc mli_Json_uint64_by_token(
        const struct mli_Json *json,
        const uint64_t token,
        uint64_t *return_uint64);
chk_rc mli_Json_cstr_by_token(
        const struct mli_Json *json,
        const uint64_t token,
        char *return_string,
        const uint64_t return_string_size);
chk_rc mli_Json_string_by_token(
        const struct mli_Json *json,
        const uint64_t token,
        struct mli_String *return_string);
chk_rc mli_Json_int64_by_key(
        const struct mli_Json *json,
        const uint64_t token,
        int64_t *val,
        const char *key);
chk_rc mli_Json_uint64_by_key(
        const struct mli_Json *json,
        const uint64_t token,
        uint64_t *val,
        const char *key);
mli_bool mli_Json_cstrcmp(
        const struct mli_Json *json,
        const uint64_t token,
        const char *str);

chk_rc mli_Json__malloc_tokens(struct mli_Json *json);
chk_rc mli_Json__parse_tokens(struct mli_Json *json);
#endif

/* json_walk */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_JSON_WALK_H_
#define MLI_JSON_WALK_H_


struct mli_JsonWalk {
        const struct mli_Json *json;
        uint64_t token;
};

struct mli_JsonWalk mli_JsonWalk_init(void);
struct mli_JsonWalk mli_JsonWalk_set(const struct mli_Json *json);
struct mli_JsonWalk mli_JsonWalk_copy(const struct mli_JsonWalk *self);

chk_rc mli_JsonWalk_to_key(struct mli_JsonWalk *self, const char *key);
chk_rc mli_JsonWalk_to_idx(struct mli_JsonWalk *self, const uint64_t idx);
void mli_JsonWalk_to_root(struct mli_JsonWalk *self);

chk_rc mli_JsonWalk_get_array_size(
        const struct mli_JsonWalk *self,
        uint64_t *size);
chk_rc mli_JsonWalk_get_string(
        const struct mli_JsonWalk *self,
        struct mli_String *val);
chk_rc mli_JsonWalk_get_int64(const struct mli_JsonWalk *self, int64_t *val);
chk_rc mli_JsonWalk_get_double(const struct mli_JsonWalk *self, double *val);

enum jsmntype_t mli_JsonWalk__type(const struct mli_JsonWalk *self);

#endif

/* lambertian */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_LAMBERTIAN_H_
#define MLI_LAMBERTIAN_H_


/*
 * Lambertian cosine law
 */

struct mli_Vec mli_lambertian_cosine_law_draw_direction_wrt_surface_normal(
        struct mli_Prng *prng,
        const struct mli_Vec surface_normal);
struct mli_Vec mli_lambertian_cosine_law_draw_direction_wrt_z(
        struct mli_Prng *prng);
#endif

/* map_json */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MAP_JSON_H_
#define MLI_MAP_JSON_H_


chk_rc mli_Map_get_value_for_string_from_json(
        const struct mli_Map *map,
        const struct mli_Json *json,
        const uint64_t token_name,
        uint32_t *out_value);
chk_rc mli_Map_insert_key_from_json(
        struct mli_Map *map,
        const struct mli_Json *json,
        const uint64_t token_name,
        const uint64_t value);

#endif

/* materials_from_archive */
/* ---------------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_MATERIALS_FORM_ARCHIVE_H_
#define MLI_MATERIALS_FORM_ARCHIVE_H_


struct mli_Materials;
struct mli_Archive;
struct mli_BoundaryLayer;

chk_rc mli_Materials__key_from_filename(
        struct mli_String *key,
        const struct mli_String *filename);

chk_rc mli_Materials_from_Archive__set_spectra(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive);

chk_rc mli_Materials_from_Archive__set_media(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive);

chk_rc mli_Materials_from_Archive__set_surfaces(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive);

chk_rc mli_BoundaryLayer_from_json_string_and_name(
        struct mli_BoundaryLayer *self,
        const struct mli_Map *surface_names,
        const struct mli_Map *media_names,
        const struct mli_String *json_string,
        const struct mli_String *name);

chk_rc mli_Materials_from_Archive__set_boundary_layers(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive);

chk_rc mli_Materials_from_Archive__set_default_medium(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive);

chk_rc mli_Materials_from_Archive(
        struct mli_Materials *materials,
        struct mli_materials_Names *names,
        const struct mli_Archive *archive);
#endif

/* object_serialize */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OBJECT_SERIALIZE_H_
#define MLI_OBJECT_SERIALIZE_H_


chk_rc mli_Object_to_io(const struct mli_Object *obj, struct mli_IO *f);
chk_rc mli_Object_from_io(struct mli_Object *obj, struct mli_IO *f);
#endif

/* object_wavefront */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OBJECT_WAVEFRONT_H_
#define MLI_OBJECT_WAVEFRONT_H_


chk_rc mli_Object_malloc_from_wavefront(
        struct mli_Object *obj,
        struct mli_IO *io);
chk_rc mli_Object_fprint_to_wavefront(
        struct mli_IO *f,
        const struct mli_Object *obj);
chk_rc mli_Object_parse_face_line(
        const struct mli_String *line,
        struct mli_object_Face *faces_vertices,
        struct mli_object_Face *faces_texture_points,
        struct mli_object_Face *faces_vertex_normals,
        int *line_mode);
chk_rc mli_Object_parse_three_float_line(
        const struct mli_String *line,
        struct mli_Vec *v);
#endif

/* photon_source */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PHOTON_SOURCE_H_
#define MLI_PHOTON_SOURCE_H_


chk_rc mli_photon_source_point_like_opening_cone_towards_z(
        struct mli_PhotonVector *out_photons,
        const double wavelength,
        const double opening_angle,
        const uint64_t num_photons,
        struct mli_Prng *prng);
chk_rc mli_photon_source_parallel_towards_z_from_xy_disc(
        struct mli_PhotonVector *out_photons,
        const double wavelength,
        const double radius,
        const uint64_t num_photons,
        struct mli_Prng *prng);
#endif

/* quaternion_json */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_QUATERNION_JSON_H_
#define MLI_QUATERNION_JSON_H_


chk_rc mli_Quaternion_tait_bryan_from_json(
        struct mli_Quaternion *quat,
        const struct mli_Json *json,
        const uint64_t token);
chk_rc mli_Quaternion_axis_angle_from_json(
        struct mli_Quaternion *quat,
        const struct mli_Json *json,
        const uint64_t token);
chk_rc mli_Quaternion_quaternion_from_json(
        struct mli_Quaternion *quat,
        const struct mli_Json *json,
        const uint64_t token);
chk_rc mli_Quaternion_from_json(
        struct mli_Quaternion *quat,
        const struct mli_Json *json,
        const uint64_t token);

#endif

/* ray_voxel_overlap */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_CORSIKA_RAY_VOXEL_OVERLAP_H_
#define MLI_CORSIKA_RAY_VOXEL_OVERLAP_H_


void mli_corsika_overlap_of_ray_with_voxels(
        const struct mli_corsika_PhotonBunch *bunch,
        const struct mli_DoubleVector *x_bin_edges,
        const struct mli_DoubleVector *y_bin_edges,
        const struct mli_DoubleVector *z_bin_edges,
        struct mli_Uint32Vector *x_idxs,
        struct mli_Uint32Vector *y_idxs,
        struct mli_Uint32Vector *z_idxs,
        struct mli_DoubleVector *overlaps);

#endif

/* string_serialize */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_STRING_SERIALIZE_H_
#define MLI_STRING_SERIALIZE_H_


chk_rc mli_String_to_io(const struct mli_String *self, struct mli_IO *f);
chk_rc mli_String_from_io(struct mli_String *self, struct mli_IO *f);

#endif

/* surface_cooktorrance */
/* -------------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_SURFACE_COOKTORRANCE_H_
#define MLI_SURFACE_COOKTORRANCE_H_


struct mli_Map;
struct mli_String;
struct mli_Materials;

struct mli_Surface_CookTorrance {
        uint64_t reflection_spectrum;
        double diffuse_weight;
        double specular_weight;
        double roughness;
};

mli_bool mli_Surface_CookTorrance_equal(
        const struct mli_Surface_CookTorrance *a,
        const struct mli_Surface_CookTorrance *b);

chk_rc mli_Surface_CookTorrance_to_io(
        const struct mli_Surface_CookTorrance *self,
        struct mli_IO *f);
chk_rc mli_Surface_CookTorrance_from_io(
        struct mli_Surface_CookTorrance *self,
        struct mli_IO *f);

chk_rc mli_Surface_CookTorrance_from_json_string(
        struct mli_Surface_CookTorrance *self,
        const struct mli_Map *spectra_names,
        const struct mli_String *json_string);

chk_rc mli_Surface_CookTorrance_valid_wrt_materials(
        const struct mli_Surface_CookTorrance *self,
        const struct mli_Materials *materials);

#endif

/* surface_transparent */
/* ------------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_SURFACE_TRANSPARENT_H_
#define MLI_SURFACE_TRANSPARENT_H_


struct mli_Map;
struct mli_String;
struct mli_Materials;

struct mli_Surface_Transparent {
        uint64_t nothing;
};

mli_bool mli_Surface_Transparent_equal(
        const struct mli_Surface_Transparent *a,
        const struct mli_Surface_Transparent *b);

chk_rc mli_Surface_Transparent_to_io(
        const struct mli_Surface_Transparent *self,
        struct mli_IO *f);
chk_rc mli_Surface_Transparent_from_io(
        struct mli_Surface_Transparent *self,
        struct mli_IO *f);

chk_rc mli_Surface_Transparent_from_json_string(
        struct mli_Surface_Transparent *self,
        const struct mli_Map *spectra_names,
        const struct mli_String *json_string);

chk_rc mli_Surface_Transparent_valid_wrt_materials(
        const struct mli_Surface_Transparent *self,
        const struct mli_Materials *materials);

#endif

/* tar */
/* --- */

/* Copyright (c) 2017 rxi
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

#ifndef MLI_TAR_H_
#define MLI_TAR_H_


#define MLI_TAR_VERSION_MAYOR 1
#define MLI_TAR_VERSION_MINOR 0
#define MLI_TAR_VERSION_PATCH 0

enum mli_tar_file_types {
        MLI_TAR_NORMAL_FILE = '0',
        MLI_TAR_HARD_LINK = '1',
        MLI_TAR_SYMBOLIC_LINK = '2',
        MLI_TAR_CHARACTER_SPECIAL = '3',
        MLI_TAR_BLOCK_SPECIAL = '4',
        MLI_TAR_DIRECTORY = '5',
        MLI_TAR_FIFO = '6'
};

enum mli_tar_name_lengths { MLI_TAR_NAME_LENGTH = 100 };

#define MLI_TAR_OCTAL 8u
#define MLI_TAR_MAX_FILESIZE_OCTAL 8589934592lu /* 8^11 */

/* basics */
/* ====== */
uint64_t mli_Tar_round_up(uint64_t n, uint64_t incr);
chk_rc mli_Tar_field_to_uint(
        uint64_t *out,
        const char *field,
        const uint64_t fieldsize);
chk_rc mli_Tar_uint_to_field(
        const uint64_t value,
        char *field,
        const uint64_t fieldsize);
chk_rc mli_Tar_uint64_to_field12_2001star_base256(uint64_t val, char *field);
chk_rc mli_Tar_field12_to_uint64_2001star_base256(
        const char *field,
        uint64_t *val);

/* header and raw header */
/* ===================== */
struct mli_TarRawHeader {
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

struct mli_TarHeader {
        uint64_t mode;
        uint64_t owner;
        uint64_t size;
        uint64_t mtime;
        uint64_t type;
        char name[MLI_TAR_NAME_LENGTH];
        char linkname[MLI_TAR_NAME_LENGTH];
};

uint64_t mli_TarRawHeader_checksum(const struct mli_TarRawHeader *rh);
mli_bool mli_TarRawHeader_is_null(const struct mli_TarRawHeader *rh);
mli_bool mli_Tar_is_known_file_type(const int file_type);
chk_rc mli_TarRawHeader_from_header(
        struct mli_TarRawHeader *rh,
        const struct mli_TarHeader *h);

struct mli_TarHeader mli_TarHeader_init(void);
chk_rc mli_TarHeader_set_directory(struct mli_TarHeader *h, const char *name);
chk_rc mli_TarHeader_set_normal_file(
        struct mli_TarHeader *h,
        const char *name,
        const uint64_t size);
chk_rc mli_TarHeader_from_raw(
        struct mli_TarHeader *h,
        const struct mli_TarRawHeader *rh);

/* tar */
/* === */
struct mli_Tar {
        struct mli_IO *stream;
        uint64_t pos;
        uint64_t remaining_data;
};

struct mli_Tar mli_Tar_init(void);

chk_rc mli_Tar_read_begin(struct mli_Tar *tar, struct mli_IO *stream);
chk_rc mli_Tar_read_header(struct mli_Tar *tar, struct mli_TarHeader *h);
chk_rc mli_Tar_read_data(struct mli_Tar *tar, void *ptr, uint64_t size);
chk_rc mli_Tar_read_finalize(struct mli_Tar *tar);

chk_rc mli_Tar_write_begin(struct mli_Tar *tar, struct mli_IO *stream);
chk_rc mli_Tar_write_header(struct mli_Tar *tar, const struct mli_TarHeader *h);
chk_rc mli_Tar_write_data(struct mli_Tar *tar, const void *data, uint64_t size);
chk_rc mli_Tar_write_finalize(struct mli_Tar *tar);

#endif

/* tar_io */
/* ------ */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLITARIO_H_
#define MLITARIO_H_


chk_rc mli_Tar_read_data_to_IO(
        struct mli_Tar *tar,
        struct mli_IO *buff,
        const uint64_t size);
chk_rc mli_Tar_write_data_from_IO(
        struct mli_Tar *tar,
        struct mli_IO *buff,
        const uint64_t size);
#endif

/* vec_json */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_VEC_JSON_H_
#define MLI_VEC_JSON_H_


chk_rc mli_Vec_from_json_token(
        struct mli_Vec *v,
        const struct mli_Json *json,
        const uint64_t token);
#endif

/* EventIo_Event */
/* ------------- */

/* Copyright 2016 Sebastian A. Mueller, Dominik Neise */
#ifndef MLI_CORSIKA_EVENTIO_EVENT_H_
#define MLI_CORSIKA_EVENTIO_EVENT_H_


struct mliEventIoEvent {
        float corsika_event_header[273];
        float corsika_event_end[273];
        struct mliDynEventIoTelescope telescopes;
};
struct mliEventIoEvent mliEventIoEvent_init(void);
void mliEventIoEvent_free(struct mliEventIoEvent *event);
int mliEventIoEvent_malloc(
        struct mliEventIoEvent *event,
        uint64_t num_photon_bunches);
int mliEventIoEvent_malloc_from_run(
        struct mliEventIoEvent *event,
        struct mliEventIoRun *run);

#endif

/* EventTape */
/* --------- */

/* Copyright 2020 Sebastian A. Mueller */
#ifndef MLI_CORSIKA_EVENTTAPE_H_
#define MLI_CORSIKA_EVENTTAPE_H_


#define MLI_CORSIKA_EVENTTAPE_VERSION_MAYOR 2
#define MLI_CORSIKA_EVENTTAPE_VERSION_MINOR 1
#define MLI_CORSIKA_EVENTTAPE_VERSION_PATCH 2

struct mliEventTapeWriter {
        struct mli_Tar tar;
        int flush_tar_stream_after_each_file;
        int run_number;
        int event_number;
        int cherenkov_bunch_block_number;
        struct mli_FloatVector buffer;
};
struct mliEventTapeWriter mliEventTapeWriter_init(void);

int mliEventTapeWriter_begin(
        struct mliEventTapeWriter *tio,
        struct mli_IO *stream,
        const uint64_t num_bunches_buffer);
int mliEventTapeWriter_finalize(struct mliEventTapeWriter *tio);

int mliEventTapeWriter_write_runh(
        struct mliEventTapeWriter *tio,
        const float *runh);
int mliEventTapeWriter_write_evth(
        struct mliEventTapeWriter *tio,
        const float *evth);
int mliEventTapeWriter_write_cherenkov_bunch(
        struct mliEventTapeWriter *tio,
        const float *bunch);
int mliEventTapeWriter_flush_cherenkov_bunch_block(
        struct mliEventTapeWriter *tio);

struct mliEventTapeReader {
        uint64_t run_number;

        /* Current event-number */
        uint64_t event_number;

        /* Current cherenkov-block-number inside the current event */
        uint64_t cherenkov_bunch_block_number;

        /* Current bunch-number inside the current cherenkov-block */
        uint64_t block_at;
        uint64_t block_size;
        int has_still_bunches_in_event;

        /* Underlying tape-archive */
        struct mli_Tar tar;

        /* Next file's tar-header in the underlying tape-archive */
        int has_tarh;
        struct mli_TarHeader tarh;
};
struct mliEventTapeReader mliEventTapeReader_init(void);

int mliEventTapeReader_begin(
        struct mliEventTapeReader *tio,
        struct mli_IO *stream);
int mliEventTapeReader_finalize(struct mliEventTapeReader *tio);

int mliEventTapeReader_read_runh(struct mliEventTapeReader *tio, float *runh);
int mliEventTapeReader_read_evth(struct mliEventTapeReader *tio, float *evth);
int mliEventTapeReader_read_cherenkov_bunch(
        struct mliEventTapeReader *tio,
        float *bunch);
int mliEventTapeReader_tarh_is_valid_cherenkov_block(
        const struct mliEventTapeReader *tio);
int mliEventTapeReader_tarh_might_be_valid_cherenkov_block(
        const struct mliEventTapeReader *tio);
int mliEventTapeReader_read_readme_until_runh(struct mliEventTapeReader *tio);

#endif

/* EventTape_testing */
/* ----------------- */

/* Copyright 2020 Sebastian A. Mueller */
#ifndef MLI_CORSIKA_EVENTTAPE_TESTING_H_
#define MLI_CORSIKA_EVENTTAPE_TESTING_H_


void mliEventTape_testing_set_random_corsika_header(
        float *head,
        struct mli_Prng *prng);
void mliEventTape_testing_set_random_RUNH(
        float *runh,
        const float run_number,
        struct mli_Prng *prng);
void mliEventTape_testing_set_random_EVTH(
        float *evth,
        const float event_number,
        const float run_number,
        struct mli_Prng *prng);
void mliEventTape_testing_set_random_bunch(float *bunch, struct mli_Prng *prng);
int mliEventTape_testing_bunches_are_equal(float *b1, float *b2);
int mliEventTape_testing_corsika_headers_are_equal(
        const float *h1,
        const float *h2);
int mliEventTape_testing_write_and_read(
        const char *path,
        const uint64_t num_events,
        const uint64_t buffer_size,
        const float *event_numbers,
        const uint64_t *num_bunches,
        const uint32_t random_seed);

#endif

/* archive */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_ARCHIVE_H_
#define MLI_ARCHIVE_H_


struct mli_Archive {
        struct mli_StringVector textfiles;
        struct mli_Map filenames;
};

struct mli_Archive mli_Archive_init(void);

void mli_Archive_free(struct mli_Archive *self);
chk_rc mli_Archive_malloc(struct mli_Archive *self);
chk_rc mli_Archive_from_io(struct mli_Archive *self, struct mli_IO *f);
chk_rc mli_Archive__from_path_cstr(struct mli_Archive *self, const char *path);
chk_rc mli_Archive_push_back(
        struct mli_Archive *self,
        const struct mli_String *filename,
        const struct mli_String *payload);
mli_bool mli_Archive_has(
        const struct mli_Archive *self,
        const struct mli_String *filename);
chk_rc mli_Archive_get(
        const struct mli_Archive *self,
        const struct mli_String *filename,
        struct mli_String **str);
uint64_t mli_Archive_size(const struct mli_Archive *self);
uint64_t mli_Archive_num_filename_prefix_sufix(
        const struct mli_Archive *self,
        const char *prefix,
        const char *sufix);

#endif

/* boundarylayer */
/* ------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_BOUNDARYLAYER_H_
#define MLI_BOUNDARYLAYER_H_


struct mli_BoundaryLayer_Side {
        uint64_t surface;
        uint64_t medium;
};

struct mli_BoundaryLayer {
        struct mli_BoundaryLayer_Side inner;
        struct mli_BoundaryLayer_Side outer;
        struct mli_String name;
};

void mli_BoundaryLayer_free(struct mli_BoundaryLayer *self);
struct mli_BoundaryLayer mli_BoundaryLayer_init(void);

mli_bool mli_BoundaryLayer_equal(
        const struct mli_BoundaryLayer *a,
        const struct mli_BoundaryLayer *b);

chk_rc mli_BoundaryLayer_to_io(
        const struct mli_BoundaryLayer *self,
        struct mli_IO *f);
chk_rc mli_BoundaryLayer_from_io(
        struct mli_BoundaryLayer *self,
        struct mli_IO *f);

#endif

/* boundarylayer_array */
/* ------------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_BOUNDARYLAYER_ARRAY_H_
#define MLI_BOUNDARYLAYER_ARRAY_H_

MLI_ARRAY_DEFINITON(mli_BoundaryLayerArray, struct mli_BoundaryLayer)
#endif

/* frame_from_archive */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FRAME_FROM_ARCHIVE_H_
#define MLI_FRAME_FROM_ARCHIVE_H_


struct mli_Object;
struct mli_Frame;
chk_rc mli_Frame_from_Archive(
        struct mli_Frame *root,
        const struct mli_Archive *archive,
        const struct mli_Map *object_names,
        const struct mli_Object *objects,
        const struct mli_Map *boundary_layer_names);
#endif

/* frame_json */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FRAME_JSON_H_
#define MLI_FRAME_JSON_H_


chk_rc mli_Frame_from_json(
        struct mli_Frame *mother,
        const struct mli_Json *json,
        const uint64_t token_children,
        const struct mli_Map *object_names,
        const struct mli_Object *objects,
        const struct mli_Map *boundary_layer_names);
chk_rc mli_Frame_id_from_json_token(
        uint32_t *id,
        const struct mli_Json *json,
        const uint64_t token);
chk_rc mli_Frame_pos_rot_from_json_token(
        struct mli_HomTraComp *frame2mother,
        const struct mli_Json *json,
        const uint64_t token);
chk_rc mli_Frame_type_from_json_token(
        uint64_t *type,
        const struct mli_Json *json,
        const uint64_t token);
chk_rc mli_Frame_boundary_layers_form_json_token(
        struct mli_Uint32Vector *boundary_layers,
        const uint32_t object_idx,
        const struct mli_Object *objects,
        const struct mli_Map *boundary_layer_names,
        const struct mli_Json *json,
        const uint64_t token);
chk_rc mli_Frame_object_reference_form_json_token(
        uint32_t *object_reference,
        const struct mli_Json *json,
        const uint64_t token,
        const struct mli_Map *object_names);
#endif

/* func_csv */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FUNC_CSV_H_
#define MLI_FUNC_CSV_H_


chk_rc mli_Func_from_csv(
        struct mli_Func *func,
        struct mli_String *xname,
        struct mli_String *yname,
        struct mli_IO *io);

#endif

/* func_info */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FUNC_INFO_H_
#define MLI_FUNC_INFO_H_


struct mli_FuncInfo {
        struct mli_String x;
        struct mli_String y;
};

struct mli_FuncInfo mli_FuncInfo_init(void);
void mli_FuncInfo_free(struct mli_FuncInfo *self);
chk_rc mli_FuncInfo_malloc(struct mli_FuncInfo *self);

mli_bool mli_FuncInfo_equal(
        const struct mli_FuncInfo *a,
        const struct mli_FuncInfo *b);
chk_rc mli_FuncInfo_to_io(const struct mli_FuncInfo *self, struct mli_IO *f);
chk_rc mli_FuncInfo_from_io(struct mli_FuncInfo *self, struct mli_IO *f);
#endif

/* func_serialize */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_FUNC_SERIALIZE_H_
#define MLI_FUNC_SERIALIZE_H_


chk_rc mli_Func_from_io(struct mli_Func *func, struct mli_IO *f);
chk_rc mli_Func_to_io(const struct mli_Func *func, struct mli_IO *f);
#endif

/* image_ppm */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_IMAGE_PPM_H_
#define MLI_IMAGE_PPM_H_


enum mli_image_ppm_color_depth {
        MLI_IMAGE_PPM_COLOR_DEPTH_8BIT = 8,
        MLI_IMAGE_PPM_COLOR_DEPTH_16BIT = 16
};

chk_rc mli_Image_to_io(
        const struct mli_Image *img,
        struct mli_IO *f,
        const uint64_t color_depth_num_bit);
chk_rc mli_Image_from_io(struct mli_Image *img, struct mli_IO *f);

chk_rc mli_image_ppm__read_color_8bit(
        struct mli_Color *color,
        struct mli_IO *f);
chk_rc mli_image_ppm__write_color_8bit(
        const struct mli_Color color,
        struct mli_IO *f);

chk_rc mli_image_ppm__read_color_16bit(
        struct mli_Color *color,
        struct mli_IO *f);
chk_rc mli_image_ppm__write_color_16bit(
        const struct mli_Color color,
        struct mli_IO *f);

#endif

/* spectrum */
/* -------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_SPECTRUM_H_
#define MLI_SPECTRUM_H_


struct mli_Spectrum {
        struct mli_Func spectrum;
        struct mli_FuncInfo info;
        struct mli_String name;
};
void mli_Spectrum_free(struct mli_Spectrum *self);
struct mli_Spectrum mli_Spectrum_init(void);
mli_bool mli_Spectrum_equal(
        const struct mli_Spectrum *a,
        const struct mli_Spectrum *b);

chk_rc mli_Spectrum_to_io(const struct mli_Spectrum *self, struct mli_IO *f);
chk_rc mli_Spectrum_from_io(struct mli_Spectrum *self, struct mli_IO *f);

chk_rc mli_Spectrum_print_to_io(
        const struct mli_Spectrum *self,
        struct mli_IO *f);

#endif

/* spectrum_array */
/* -------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_SPECTRUM_ARRAY_H_
#define MLI_SPECTRUM_ARRAY_H_

MLI_ARRAY_DEFINITON(mli_SpectrumArray, struct mli_Spectrum)
#endif

/* surface */
/* ------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_SURFACE_H_
#define MLI_SURFACE_H_

struct mli_IO;
struct mli_Map;
struct mli_Materials;

union mli_SurfaceData {
        struct mli_Surface_Transparent transparent;
        struct mli_Surface_CookTorrance cooktorrance;
};

struct mli_Surface {
        struct mli_String name;
        uint64_t type;
        union mli_SurfaceData data;
};
struct mli_Surface mli_Surface_init(void);

void mli_Surface_free(struct mli_Surface *self);

mli_bool mli_Surface_equal(
        const struct mli_Surface *a,
        const struct mli_Surface *b);

chk_rc mli_Surface_to_io(const struct mli_Surface *self, struct mli_IO *f);
chk_rc mli_Surface_from_io(struct mli_Surface *self, struct mli_IO *f);

chk_rc mli_Surface_from_json_string_and_name(
        struct mli_Surface *self,
        const struct mli_Map *spectra_names,
        const struct mli_String *json_string,
        const struct mli_String *name);

chk_rc mli_Surface_valid_wrt_materials(
        const struct mli_Surface *self,
        const struct mli_Materials *materials);

#endif

/* surface_array */
/* ------------- */

/* Copyright 2018-2024 Sebastian Achim Mueller */
#ifndef MLI_SURFACE_ARRAY_H_
#define MLI_SURFACE_ARRAY_H_

MLI_ARRAY_DEFINITON(mli_SurfaceArray, struct mli_Surface)
#endif

/* materials */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MATERIALS_H_
#define MLI_MATERIALS_H_


struct mli_MaterialsCapacity {
        uint64_t num_spectra;
        uint64_t num_surfaces;
        uint64_t num_media;
        uint64_t num_boundary_layers;
};

struct mli_MaterialsCapacity mli_MaterialsCapacity_init(void);

struct mli_Materials {
        struct mli_SpectrumArray spectra;
        struct mli_SurfaceArray surfaces;
        struct mli_MediumArray media;
        struct mli_BoundaryLayerArray boundary_layers;

        uint64_t default_medium;
};

chk_rc mli_Materials_malloc(
        struct mli_Materials *self,
        const struct mli_MaterialsCapacity rescap);
void mli_Materials_free(struct mli_Materials *self);
struct mli_Materials mli_Materials_init(void);
chk_rc mli_Materials_info_fprint(FILE *f, const struct mli_Materials *self);

mli_bool mli_Materials__has_surface_name_cstr(
        const struct mli_Materials *self,
        const char *name);
mli_bool mli_Materials__has_medium_name_cstr(
        const struct mli_Materials *self,
        const char *name);

#endif

/* materials_equal */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MATERIALS_EQUAL_H_
#define MLI_MATERIALS_EQUAL_H_


mli_bool mli_Materials_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b);
mli_bool mli_Materials_spectra_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b);
mli_bool mli_Materials_surfaces_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b);
mli_bool mli_Materials_media_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b);
mli_bool mli_Materials_boundary_layers_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b);
mli_bool mli_Materials_default_medium_equal(
        const struct mli_Materials *a,
        const struct mli_Materials *b);
#endif

/* materials_serialize */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MATERIALS_SERIALIZE_H_
#define MLI_MATERIALS_SERIALIZE_H_


chk_rc mli_Materials_to_io(const struct mli_Materials *self, struct mli_IO *f);
chk_rc mli_Materials_from_io(struct mli_Materials *self, struct mli_IO *f);
#endif

/* materials_valid */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_MATERIALS_VALID_H_
#define MLI_MATERIALS_VALID_H_


mli_bool mli_Materials_valid(const struct mli_Materials *self);
mli_bool mli_Materials_valid_default_medium(const struct mli_Materials *self);
mli_bool mli_Materials_valid_spectra(const struct mli_Materials *self);
mli_bool mli_Materials_valid_surfaces(const struct mli_Materials *self);
mli_bool mli_Materials_valid_media(const struct mli_Materials *self);
mli_bool mli_Materials_valid_boundary_layers(const struct mli_Materials *self);
#endif

/* photon_interaction */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PHOTON_INTERACTION_H_
#define MLI_PHOTON_INTERACTION_H_


enum mli_photon_interaction_types {
        MLI_PHOTON_CREATION = 101u,
        MLI_PHOTON_ABSORPTION = 102u,
        MLI_PHOTON_ABSORPTION_MEDIUM = 103u,
        MLI_PHOTON_FRESNEL_REFLECTION = 104u,
        MLI_PHOTON_REFRACTION = 105u,
        MLI_PHOTON_SPECULAR_REFLECTION = 106u,
        MLI_PHOTON_DIFFUSE_REFLECTION = 107u
};

struct mli_PhotonInteraction {
        int32_t on_geometry_surface;
        struct mli_GeometryId geometry_id;

        struct mli_Vec position;
        struct mli_Vec position_local;
        double distance_of_ray;

        uint64_t medium_coming_from;
        uint64_t medium_going_to;

        int32_t from_outside_to_inside;
        int32_t type;
};

chk_rc mli_photon_time_of_flight(
        const struct mli_Materials *materials,
        const struct mli_PhotonInteraction *phisec,
        const double wavelength,
        double *time_of_flight);
chk_rc mli_photon_interaction_type_to_string(const int32_t type, char *s);
#endif

/* color_materials */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_COLOR_MATERIALS_H_
#define MLI_COLOR_MATERIALS_H_


struct mli_ColorMaterials {
        struct mli_ColorSpectrumBinEdges wavelength_bin_edges;

        struct mli_ColorSpectrum observer_matching_curve_x;
        struct mli_ColorSpectrum observer_matching_curve_y;
        struct mli_ColorSpectrum observer_matching_curve_z;

        struct mli_Mat observer_matching_curve_xyz_to_rgb;

        struct mli_ColorSpectrumArray spectra;
};

struct mli_Vec mli_ColorMaterials_ColorSpectrum_to_xyz(
        const struct mli_ColorMaterials *self,
        const struct mli_ColorSpectrum *spectrum);

struct mli_ColorMaterials mli_ColorMaterials_init(void);
chk_rc mli_ColorMaterials_malloc(
        struct mli_ColorMaterials *self,
        const uint64_t num_spectra);
chk_rc mli_ColorMaterials_set_observer_cie1931(struct mli_ColorMaterials *self);

chk_rc mli_ColorMaterials_malloc_from_Materials(
        struct mli_ColorMaterials *self,
        const struct mli_Materials *materials);

void mli_ColorMaterials_free(struct mli_ColorMaterials *self);
#endif

/* geometry */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRY_H_
#define MLI_GEOMETRY_H_


struct mli_Geometry {
        uint32_t num_objects;
        struct mli_Object *objects;
        struct mli_String *object_names;

        uint32_t num_robjects;
        uint32_t *robjects;
        uint32_t *robject_ids;
        struct mli_HomTraComp *robject2root;
};

chk_rc mli_Geometry_malloc(
        struct mli_Geometry *self,
        const uint32_t num_objects,
        const uint32_t num_robjects);
chk_rc mli_Geometry_malloc_references(
        struct mli_Geometry *self,
        const uint32_t num_robjects);
chk_rc mli_Geometry_malloc_objects(
        struct mli_Geometry *self,
        const uint32_t num_objects);

void mli_Geometry_free(struct mli_Geometry *self);
void mli_Geometry_free_objects(struct mli_Geometry *self);
void mli_Geometry_free_references(struct mli_Geometry *self);

struct mli_Geometry mli_Geometry_init(void);
void mli_Geometry_init_objects(struct mli_Geometry *self);
void mli_Geometry_init_references(struct mli_Geometry *self);

void mli_Geometry_info_fprint(FILE *f, const struct mli_Geometry *self);
struct mli_BoundaryLayer mli_Geometry_object_surfaces(
        const struct mli_Geometry *self,
        const uint32_t object_idx);
chk_rc mli_Geometry_warn_objects(const struct mli_Geometry *self);
#endif

/* geometry_equal */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLIGEOMETRY_EQUAL_H_
#define MLIGEOMETRY_EQUAL_H_


mli_bool mli_Geometry_equal(
        const struct mli_Geometry *a,
        const struct mli_Geometry *b);
mli_bool mli_Geometry_objects_equal(
        const struct mli_Geometry *a,
        const struct mli_Geometry *b);
mli_bool mli_Geometry_object_references_equal(
        const struct mli_Geometry *a,
        const struct mli_Geometry *b);
#endif

/* geometry_serialize */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRY_SERIALIZE_H_
#define MLI_GEOMETRY_SERIALIZE_H_


chk_rc mli_Geometry_to_io(const struct mli_Geometry *scenery, struct mli_IO *f);
chk_rc mli_Geometry_from_io(struct mli_Geometry *scenery, struct mli_IO *f);
#endif

/* geometry_valid */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRY_VALID_H_
#define MLI_GEOMETRY_VALID_H_


mli_bool mli_Geometry_valid(const struct mli_Geometry *geometry);
mli_bool mli_Geometry_valid_objects(const struct mli_Geometry *geometry);
mli_bool mli_Geometry_valid_robjects_HomTras(
        const struct mli_Geometry *geometry);
mli_bool mli_Geometry_valid_object_references(
        const struct mli_Geometry *geometry);
#endif

/* geometrytomaterialmap */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRYTOMATERIALMAP_H_
#define MLI_GEOMETRYTOMATERIALMAP_H_


struct mli_GeometryToMaterialMap {
        uint32_t num_robjects;
        uint32_t total_num_boundary_layers;
        uint32_t *boundary_layers;
        uint32_t *first_boundary_layer_in_robject;
};

struct mli_GeometryToMaterialMap mli_GeometryToMaterialMap_init(void);
chk_rc mli_GeometryToMaterialMap_malloc(
        struct mli_GeometryToMaterialMap *map,
        const uint32_t num_robjects,
        const uint32_t total_num_boundary_layers);
void mli_GeometryToMaterialMap_free(struct mli_GeometryToMaterialMap *map);

uint32_t mli_GeometryToMaterialMap_resolve_idx(
        const struct mli_GeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx);

uint32_t mli_GeometryToMaterialMap_get(
        const struct mli_GeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx);
void mli_GeometryToMaterialMap_set(
        const struct mli_GeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx,
        const uint32_t boundary_layer_idx);

uint32_t mli_GeometryToMaterialMap_num_boundary_layers_in_robject(
        const struct mli_GeometryToMaterialMap *map,
        const uint32_t robject_idx);

void mli_GeometryToMaterialMap_info_fprint(
        FILE *f,
        const struct mli_GeometryToMaterialMap *map);
#endif

/* geometrytomaterialmap_equal */
/* --------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRYTOMATERIALMAP_EQUAL_H_
#define MLI_GEOMETRYTOMATERIALMAP_EQUAL_H_


mli_bool mli_GeometryToMaterialMap_equal(
        const struct mli_GeometryToMaterialMap *a,
        const struct mli_GeometryToMaterialMap *b);
#endif

/* geometrytomaterialmap_serialize */
/* ------------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRYTOMATERIALMAP_SERIALIZE_H_
#define MLI_GEOMETRYTOMATERIALMAP_SERIALIZE_H_


chk_rc mli_GeometryToMaterialMap_from_io(
        struct mli_GeometryToMaterialMap *geomap,
        struct mli_IO *f);
chk_rc mli_GeometryToMaterialMap_to_io(
        const struct mli_GeometryToMaterialMap *geomap,
        struct mli_IO *f);
#endif

/* geometrytomaterialmap_valid */
/* --------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRYTOMATERIALMAP_EQUAL_H_
#define MLI_GEOMETRYTOMATERIALMAP_EQUAL_H_


mli_bool mli_GeometryToMaterialMap_valid(
        const struct mli_GeometryToMaterialMap *geomap);

mli_bool mli_GeometryToMaterialMap_valid_wrt_Geometry(
        const struct mli_GeometryToMaterialMap *geomap,
        const struct mli_Geometry *geometry);
mli_bool mli_GeometryToMaterialMap_valid_wrt_Materials(
        const struct mli_GeometryToMaterialMap *geomap,
        const struct mli_Materials *materials);
#endif

/* octree_tmp */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OCTREE_TMP_H_
#define MLI_OCTREE_TMP_H_


#define MLI_OCTREE_TMPNODE_FLAT_INDEX_NONE -1

uint64_t mli_octree_guess_depth_based_on_num_objects(
        const uint64_t num_objects);

/*
 * The dynamic node
 * ================
 */

struct mli_octree_TmpNode {
        struct mli_octree_TmpNode *children[8];
        uint32_t num_objects;
        uint32_t *objects;

        int32_t flat_index;
        int32_t node_index;
        int32_t leaf_index;
};

chk_rc mli_octree_TmpNode_malloc(
        struct mli_octree_TmpNode *n,
        const uint32_t num_objects);
void mli_octree_TmpNode_free(struct mli_octree_TmpNode *n);
struct mli_octree_TmpNode mli_octree_TmpNode_init(void);
void mli_octree_TmpNode_num_nodes_leafs_objects(
        const struct mli_octree_TmpNode *root_node,
        uint64_t *num_nodes,
        uint64_t *num_leafs,
        uint64_t *num_object_links);
void mli_octree_TmpNode_num_nodes_leafs_objects_walk(
        const struct mli_octree_TmpNode *node,
        uint64_t *num_nodes,
        uint64_t *num_leafs,
        uint64_t *num_object_links);
void mli_octree_TmpNode_set_flat_index(struct mli_octree_TmpNode *root_node);
void mli_octree_TmpNode_set_flat_index_walk(
        struct mli_octree_TmpNode *node,
        int32_t *flat_index,
        int32_t *node_index,
        int32_t *leaf_index);
mli_bool mli_octree_TmpNode_exists_and_has_objects(
        const struct mli_octree_TmpNode *node);
void mli_octree_TmpNode_print(
        const struct mli_octree_TmpNode *node,
        const uint32_t indent,
        const uint32_t child);
uint32_t mli_octree_TmpNode_num_children(const struct mli_octree_TmpNode *node);
chk_rc mli_octree_TmpNode_malloc_tree_from_bundle(
        struct mli_octree_TmpNode *root_node,
        const void *bundle,
        const uint32_t num_items_in_bundle,
        mli_bool (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mli_AABB),
        const struct mli_Cube bundle_cube);
chk_rc mli_octree_TmpNode_add_children(
        struct mli_octree_TmpNode *node,
        const void *bundle,
        mli_bool (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mli_AABB),
        const struct mli_Cube cube,
        const uint64_t depth,
        const uint64_t max_depth);
uint32_t mli_octree_TmpNode_signs_to_child(
        const uint32_t sx,
        const uint32_t sy,
        const uint32_t sz);

/*
 * The dynamic octree
 * ==================
 */

struct mli_octree_TmpOcTree {
        struct mli_Cube cube;
        struct mli_octree_TmpNode root;
};

chk_rc mli_octree_TmpOcTree_malloc_from_bundle(
        struct mli_octree_TmpOcTree *octree,
        const void *bundle,
        const uint32_t num_items_in_bundle,
        mli_bool (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mli_AABB),
        struct mli_AABB bundle_aabb);

void mli_octree_TmpOcTree_free(struct mli_octree_TmpOcTree *octree);
struct mli_octree_TmpOcTree mli_octree_TmpOcTree_init(void);
void mli_octree_TmpOcTree_print(const struct mli_octree_TmpOcTree *octree);

#endif

/* atmosphere */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_ATMOSPHERE_H_
#define MLI_ATMOSPHERE_H_


struct mli_Atmosphere {
        double sunLatitude;
        double sunHourAngle;

        struct mli_Vec sunDirection;
        double sunDistance;
        double sunRadius;

        double earthRadius;
        double atmosphereRadius;

        double Height_Rayleigh;
        double Height_Mie;

        struct mli_ColorSpectrum beta_Rayleigh_spectrum;
        struct mli_ColorSpectrum beta_Mie_spectrum;
        struct mli_ColorSpectrum sun_spectrum;

        uint64_t numSamples;
        uint64_t numSamplesLight;

        double power;
        double altitude;
};

struct mli_Atmosphere mli_Atmosphere_init(void);
void mli_Atmosphere_set_sun_direction(
        struct mli_Atmosphere *self,
        const double sunLatitude,
        const double sunHourAngle);
void mli_Atmosphere_increase_latitude(
        struct mli_Atmosphere *self,
        const double increment);
void mli_Atmosphere_decrease_latitude(
        struct mli_Atmosphere *self,
        const double increment);
void mli_Atmosphere_increase_hours(
        struct mli_Atmosphere *self,
        const double increment);
void mli_Atmosphere_decrease_hours(
        struct mli_Atmosphere *self,
        const double increment);
void mli_Atmosphere_increase_altitude(
        struct mli_Atmosphere *self,
        const double factor);
void mli_Atmosphere_decrease_altitude(
        struct mli_Atmosphere *self,
        const double factor);

void mli_ColorSpectrum_set_beta_rayleigh(struct mli_ColorSpectrum *self);

struct mli_ColorSpectrum mli_Atmosphere_query(
        const struct mli_Atmosphere *self,
        const struct mli_Vec orig,
        const struct mli_Vec dir);

struct mli_ColorSpectrum mli_Atmosphere_hit_earth_body(
        const struct mli_Atmosphere *self,
        const struct mli_Vec orig,
        const struct mli_Vec dir);

struct mli_ColorSpectrum mli_Atmosphere_hit_outer_atmosphere(
        const struct mli_Atmosphere *self,
        const struct mli_Vec orig,
        const struct mli_Vec dir,
        double tmin,
        double tmax);

struct mli_ColorSpectrum mli_Atmosphere_compute_depth(
        const struct mli_Atmosphere *self,
        const struct mli_Vec orig,
        const struct mli_Vec dir,
        double tmin,
        double tmax);

#endif

/* atmosphere_json */
/* --------------- */

/* Copyright 2018-2021 Sebastian Achim Mueller */
#ifndef MLI_ATMOSPHERE_JSON_H_
#define MLI_ATMOSPHERE_JSON_H_


chk_rc mli_Atmosphere_from_json_token(
        struct mli_Atmosphere *atm,
        const struct mli_Json *json,
        const uint64_t tkn);

#endif

/* octree */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OCTREE_H_
#define MLI_OCTREE_H_


enum mli_octree_type {
        MLI_OCTREE_TYPE_NONE = 0,
        MLI_OCTREE_TYPE_NODE = 1,
        MLI_OCTREE_TYPE_LEAF = 2
};

struct mli_octree_LeafAddress {
        uint32_t first_object_link;
        uint32_t num_object_links;
};

struct mli_octree_LeafArray {
        uint64_t num_leafs;
        struct mli_octree_LeafAddress *adresses;
        uint64_t num_object_links;
        uint32_t *object_links;
};

struct mli_octree_Node {
        uint32_t children[8];
        uint8_t types[8];
};

struct mli_OcTree {
        struct mli_Cube cube;
        uint64_t num_nodes;
        struct mli_octree_Node *nodes;
        struct mli_octree_LeafArray leafs;
        uint8_t root_type;
};

void mli_OcTree_print(const struct mli_OcTree *tree);
void mli_OcTree_print_walk(
        const struct mli_OcTree *tree,
        const int32_t node_idx,
        const uint8_t node_type,
        const uint32_t indent,
        const uint32_t child);
mli_bool mli_OcTree_equal_payload(
        const struct mli_OcTree *tree,
        const struct mli_octree_TmpOcTree *tmp_octree);
mli_bool mli_OcTree_equal_payload_walk(
        const struct mli_OcTree *tree,
        const int32_t node_idx,
        const int32_t node_type,
        const struct mli_octree_TmpNode *tmp_node);
uint32_t mli_OcTree_leaf_object_link(
        const struct mli_OcTree *tree,
        const uint64_t leaf,
        const uint64_t object_link);
uint64_t mli_OcTree_leaf_num_objects(
        const struct mli_OcTree *tree,
        const uint64_t leaf);
uint64_t mli_OcTree_node_num_children(
        const struct mli_OcTree *tree,
        const uint64_t node_idx);
void mli_OcTree_set(
        struct mli_OcTree *tree,
        const struct mli_octree_TmpOcTree *dyntree);
void mli_OcTree_set_walk(
        struct mli_OcTree *tree,
        const struct mli_octree_TmpNode *dynnode,
        uint64_t *object_link_size);
void mli_OcTree_set_leaf(
        struct mli_OcTree *tree,
        const struct mli_octree_TmpNode *dynnode,
        uint64_t *object_link_size);
void mli_OcTree_set_node(
        struct mli_OcTree *tree,
        const struct mli_octree_TmpNode *dynnode);
chk_rc mli_OcTree_malloc(
        struct mli_OcTree *tree,
        const uint64_t num_nodes,
        const uint64_t num_leafs,
        const uint64_t num_object_links);
void mli_OcTree_free(struct mli_OcTree *tree);
struct mli_OcTree mli_OcTree_init(void);
struct mli_octree_Node mli_octree_Node_init(void);
chk_rc mli_octree_LeafArray_malloc(
        struct mli_octree_LeafArray *leafs,
        const uint64_t num_leafs,
        const uint64_t num_object_links);
void mli_octree_LeafArray_free(struct mli_octree_LeafArray *leafs);
struct mli_octree_LeafArray mli_octree_LeafArray_init(void);
struct mli_octree_LeafAddress mli_octree_LeafAddress_init(void);

chk_rc mli_OcTree_malloc_from_object_wavefront(
        struct mli_OcTree *octree,
        const struct mli_Object *object);

struct mli_GeometryAndAccelerator;
chk_rc mli_OcTree_malloc_from_Geometry(
        struct mli_OcTree *octree,
        const struct mli_GeometryAndAccelerator *accgeo,
        const struct mli_AABB outermost_aabb);

#endif

/* octree_equal */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OCTREE_EQUAL_H_
#define MLI_OCTREE_EQUAL_H_


mli_bool mli_OcTree_equal(
        const struct mli_OcTree *a,
        const struct mli_OcTree *b);
#endif

/* octree_serialize */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OCTREE_SERIALIZE_H_
#define MLI_OCTREE_SERIALIZE_H_


chk_rc mli_OcTree_to_io(const struct mli_OcTree *octree, struct mli_IO *f);
chk_rc mli_OcTree_from_io(struct mli_OcTree *octree, struct mli_IO *f);

#endif

/* octree_valid */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_OCTREE_VALID_H_
#define MLI_OCTREE_VALID_H_


mli_bool mli_OcTree_valid(const struct mli_OcTree *octree);
mli_bool mli_OcTree_valid_wrt_links(
        const struct mli_OcTree *octree,
        const uint32_t num_links);
#endif

/* pathtracer */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PATHTRACER_H_
#define MLI_PATHTRACER_H_


struct mli_pathtracer_Config;
struct mli_Scenery;
struct mli_Prng;
struct mli_IntersectionSurfaceNormal;
struct mli_IntersectionLayer;

struct mli_PathTracer {
        const struct mli_Scenery *scenery;
        const struct mli_ColorMaterials *scenery_color_materials;
        const struct mli_pathtracer_Config *config;
};

struct mli_PathTracer mli_pathtracer_init(void);

double mli_pathtracer_estimate_sun_obstruction_weight(
        const struct mli_PathTracer *tracer,
        const struct mli_Vec position,
        struct mli_Prng *prng);

double mli_pathtracer_estimate_sun_visibility_weight(
        const struct mli_PathTracer *tracer,
        const struct mli_Vec position,
        struct mli_Prng *prng);

struct mli_Color mli_pathtracer_trace_ray(
        const struct mli_PathTracer *tracer,
        const struct mli_Ray ray,
        struct mli_Prng *prng);

struct mli_ColorSpectrum mli_pathtracer_trace_ambient_background(
        const struct mli_PathTracer *tracer,
        const struct mli_Ray ray);

struct mli_ColorSpectrum mli_pathtracer_trace_ambient_background_atmosphere(
        const struct mli_PathTracer *tracer,
        const struct mli_Ray ray);

struct mli_ColorSpectrum mli_pathtracer_trace_ambient_background_whitebox(
        const struct mli_PathTracer *tracer);

struct mli_ColorSpectrum mli_pathtracer_trace_ambient_sun(
        const struct mli_PathTracer *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        struct mli_Prng *prng);

struct mli_ColorSpectrum mli_pathtracer_trace_ambient_sun_atmosphere(
        const struct mli_PathTracer *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        struct mli_Prng *prng);

struct mli_ColorSpectrum mli_pathtracer_trace_ambient_sun_whitebox(
        const struct mli_PathTracer *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        struct mli_Prng *prng);

struct mli_ColorSpectrum mli_pathtracer_trace_path_to_next_intersection(
        const struct mli_PathTracer *tracer,
        const struct mli_Ray ray,
        struct mli_pathtracer_Path path,
        struct mli_Prng *prng);

struct mli_ColorSpectrum mli_pathtracer_trace_next_intersection(
        const struct mli_PathTracer *tracer,
        const struct mli_Ray ray,
        const struct mli_IntersectionSurfaceNormal *intersection,
        struct mli_pathtracer_Path path,
        struct mli_Prng *prng);

struct mli_ColorSpectrum mli_pathtracer_trace_intersection_cooktorrance(
        const struct mli_PathTracer *tracer,
        const struct mli_Ray ray,
        const struct mli_IntersectionSurfaceNormal *intersection,
        const struct mli_IntersectionLayer *intersection_layer,
        struct mli_pathtracer_Path path,
        struct mli_Prng *prng);

struct mli_ColorSpectrum mli_pathtracer_trace_intersection_transparent(
        const struct mli_PathTracer *tracer,
        const struct mli_Ray ray,
        const struct mli_IntersectionSurfaceNormal *intersection,
        const struct mli_IntersectionLayer *intersection_layer,
        struct mli_pathtracer_Path path,
        struct mli_Prng *prng);

#endif

/* pathtracer_config */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PATHTRACER_CONFIG_H_
#define MLI_PATHTRACER_CONFIG_H_


struct mli_Scenery;
struct mli_Prng;

struct mli_pathtracer_Config {
        uint64_t num_trails_global_light_source;

        int have_atmosphere;
        struct mli_Atmosphere atmosphere;

        struct mli_ColorSpectrum ambient_radiance_W_per_m2_per_sr;
};

struct mli_pathtracer_Config mli_pathtracer_Config_init(void);

#endif

/* pathtracer_config_json */
/* ---------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PATHTRACER_CONFIG_JSON_H_
#define MLI_PATHTRACER_CONFIG_JSON_H_


chk_rc mli_pathtracer_Config_from_json_token(
        struct mli_pathtracer_Config *tc,
        const struct mli_Json *json,
        const uint64_t tkn);
#endif

/* ray_octree_traversal */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RAYTRACING_RAY_OCTREE_TRAVERSAL_H_
#define MLI_RAYTRACING_RAY_OCTREE_TRAVERSAL_H_



#define mli_octree_node int
#define MLI_RAYTRACING_RAY_OCTREE_TRAVERSAL_EPSILON 1.0e-307

void mli_raytracing_ray_octree_traversal(
        const struct mli_OcTree *octree,
        const struct mli_Ray ray,
        void *work,
        void (*work_on_leaf_node)(
                void *,
                const struct mli_OcTree *,
                const uint32_t));

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
                const uint32_t));

mli_octree_node mli_raytracing_ray_octree_traversal_next_octree_node(
        const struct mli_Vec tm,
        mli_octree_node x,
        mli_octree_node y,
        mli_octree_node z);

mli_octree_node mli_raytracing_ray_octree_traversal_first_octree_node(
        const struct mli_Vec t0,
        const struct mli_Vec tm);

#endif

/* accelerator */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_ACCELERATOR_H_
#define MLI_ACCELERATOR_H_


struct mli_Accelerator {
        uint32_t num_objects;
        struct mli_OcTree *object_octrees;

        uint32_t num_robjects;
        struct mli_AABB *robject_aabbs;

        struct mli_OcTree scenery_octree;
};

struct mli_Accelerator mli_Accelerator_init(void);

void mli_Accelerator_free(struct mli_Accelerator *self);

chk_rc mli_Accelerator_malloc(
        struct mli_Accelerator *self,
        const uint32_t num_objects,
        const uint32_t num_robjects);

chk_rc mli_Accelerator_malloc_from_Geometry(
        struct mli_Accelerator *self,
        const struct mli_Geometry *geometry);

chk_rc mli_Accelerator_set_robject_aabbs(
        struct mli_Accelerator *self,
        const struct mli_Geometry *geometry);

chk_rc mli_Accelerator_set_object_octrees(
        struct mli_Accelerator *self,
        const struct mli_Geometry *geometry);

void mli_Accelerator_info_fprint(FILE *f, const struct mli_Accelerator *self);

struct mli_AABB mli_Accelerator_outermost_aabb(
        const struct mli_Accelerator *self);

#endif

/* accelerator_equal */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_ACCELERATOR_EQUAL_H_
#define MLI_ACCELERATOR_EQUAL_H_


mli_bool mli_Accelerator_equal(
        const struct mli_Accelerator *a,
        const struct mli_Accelerator *b);
#endif

/* accelerator_serialize */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_ACCELERATOR_SERIALIZE_H_
#define MLI_ACCELERATOR_SERIALIZE_H_


chk_rc mli_Accelerator_from_io(struct mli_Accelerator *self, struct mli_IO *f);
chk_rc mli_Accelerator_to_io(
        const struct mli_Accelerator *self,
        struct mli_IO *f);
#endif

/* accelerator_valid */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_ACCELERATOR_VALID_H_
#define MLI_ACCELERATOR_VALID_H_


mli_bool mli_Accelerator_valid(const struct mli_Accelerator *self);
mli_bool mli_Accelerator_valid_wrt_Geometry(
        const struct mli_Accelerator *self,
        const struct mli_Geometry *geometry);
#endif

/* geometry_aabb */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRY_AABB_H_
#define MLI_GEOMETRY_AABB_H_


mli_bool mli_Geometry_robject_has_overlap_aabb_void(
        const void *accgeo,
        const uint32_t robject_idx,
        const struct mli_AABB aabb);

mli_bool mli_Geometry_robject_has_overlap_aabb(
        const struct mli_GeometryAndAccelerator *accgeo,
        const uint32_t robject_idx,
        const struct mli_AABB aabb);

#endif

/* geometry_and_accelerator */
/* ------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_GEOMETRY_AND_ACCELERATOR_H_
#define MLI_GEOMETRY_AND_ACCELERATOR_H_


struct mli_GeometryAndAccelerator {
        /*
         * A temporary container to allow access to both geometry and its
         * accelerator using only one pointer.
         */
        const struct mli_Geometry *geometry;
        const struct mli_Accelerator *accelerator;
};

#endif

/* scenery */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_SCENERY_H_
#define MLI_SCENERY_H_


struct mli_Scenery {
        struct mli_Geometry geometry;
        struct mli_Accelerator accelerator;
        struct mli_Materials materials;
        struct mli_GeometryToMaterialMap geomap;
};

struct mli_Scenery mli_Scenery_init(void);
void mli_Scenery_free(struct mli_Scenery *self);
void mli_Scenery_info_fprint(FILE *f, const struct mli_Scenery *self);
uint32_t mli_Scenery_resolve_boundary_layer_idx(
        const struct mli_Scenery *scenery,
        const struct mli_GeometryId geometry_id);
#endif

/* scenery_equal */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_SCENERY_EQUAL_H_
#define MLI_SCENERY_EQUAL_H_


mli_bool mli_Scenery_equal(
        const struct mli_Scenery *a,
        const struct mli_Scenery *b);
#endif

/* scenery_minimal_object */
/* ---------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_SCENERY_MINIMAL_OBJECT_H_
#define MLI_SCENERY_MINIMAL_OBJECT_H_


chk_rc mli_Scenery_malloc_minimal_from_wavefront(
        struct mli_Scenery *self,
        const char *path);
#endif

/* scenery_serialize */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_SCENERY_SERIALIZE_H_
#define MLI_SCENERY_SERIALIZE_H_


chk_rc mli_Scenery_to_io(const struct mli_Scenery *self, struct mli_IO *f);
chk_rc mli_Scenery_from_io(struct mli_Scenery *self, struct mli_IO *f);

chk_rc mli_Scenery_malloc_from_path(struct mli_Scenery *self, const char *path);
chk_rc mli_Scenery_write_to_path(
        const struct mli_Scenery *self,
        const char *path);
#endif

/* scenery_tar */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_SCENERY_TAR_H_
#define MLI_SCENERY_TAR_H_


chk_rc mli_Scenery_from_io_tar(struct mli_Scenery *self, struct mli_IO *f);
chk_rc mli_Scenery__from_path_cstr(struct mli_Scenery *self, const char *path);
chk_rc mli_Scenery_malloc_from_Archive(
        struct mli_Scenery *self,
        const struct mli_Archive *archive);
#endif

/* scenery_valid */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_SCENERY_VALID_H_
#define MLI_SCENERY_VALID_H_


mli_bool mli_Scenery_valid(const struct mli_Scenery *self);
#endif

/* viewer */
/* ------ */

/* Copyright 2019 Sebastian Achim Mueller                                     */
#ifndef MLI_VIEWER_VIEWER_H_
#define MLI_VIEWER_VIEWER_H_


#define mli_key_code int

enum mli_viewer_key_codes {
        MLI_VIEWER_ESCAPE_KEY = 27,
        MLI_VIEWER_SPACE_KEY = 32
};

mli_key_code mli_viewer_get_key(void);

void mli_viewer_clear_screen(void);

void mli_viewer_print_help(void);

void mli_viewer_print_info_line(
        const struct mli_View view,
        const struct mli_viewer_Cursor cursor,
        const struct mli_pathtracer_Config tracer_config,
        const double gamma,
        const double gain);

void mli_viewer_timestamp_now_19chars(char *buffer);

chk_rc mli_viewer_export_image(
        const struct mli_PathTracer *tracer,
        const struct mli_viewer_Config config,
        const struct mli_View view,
        struct mli_Prng *prng,
        const double object_distance,
        const double gamma,
        const double gain,
        const char *path);

chk_rc mli_viewer_run_interactive_viewer(
        const struct mli_Scenery *scenery,
        const struct mli_viewer_Config config);
chk_rc mli_viewer_run_interactive_viewer_try_non_canonical_stdin(
        const struct mli_Scenery *scenery,
        const struct mli_viewer_Config config);
#endif

/* intersection_and_scenery */
/* ------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RAYTRACING_INTERSECTION_AND_SCENERY_H_
#define MLI_RAYTRACING_INTERSECTION_AND_SCENERY_H_


struct mli_IntersectionLayerSide {
        const struct mli_Surface *surface;
        uint64_t surface_idx;
        const struct mli_Medium *medium;
        uint64_t medium_idx;
};

struct mli_IntersectionLayerSide mli_IntersectionLayerSide_init(void);

struct mli_IntersectionLayer {
        struct mli_IntersectionLayerSide side_coming_from;
        struct mli_IntersectionLayerSide side_going_to;
};

struct mli_IntersectionLayer mli_IntersectionLayer_init(void);

struct mli_IntersectionLayer mli_raytracing_get_intersection_layer(
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec);

struct mli_BoundaryLayer_Side mli_raytracing_get_side_going_to(
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec);
struct mli_BoundaryLayer_Side mli_raytracing_get_side_coming_from(
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec);

const struct mli_Func *mli_raytracing_get_refractive_index_coming_from(
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec);
const struct mli_Func *mli_raytracing_get_refractive_index_going_to(
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec);

#endif

/* pathtracer_atmosphere */
/* --------------------- */

/* Copyright 2018-2023 Sebastian Achim Mueller */
#ifndef MLI_PATHTRACER_ATMOSPHERE_H_
#define MLI_PATHTRACER_ATMOSPHERE_H_


struct mli_Prng;
struct mli_PathTracer;
struct mli_pathtracer_Config;
struct mli_IntersectionSurfaceNormal;

struct mli_ColorSpectrum mli_raytracing_color_tone_of_sun(
        const struct mli_pathtracer_Config *config,
        const struct mli_Vec support);
struct mli_ColorSpectrum mli_raytracing_color_tone_of_diffuse_sky(
        const struct mli_PathTracer *tracer,
        const struct mli_IntersectionSurfaceNormal *intersection,
        struct mli_Prng *prng);

#endif

/* photon_interaction_vector */
/* ------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PHOTON_INTERACTION_VECTOR_H_
#define MLI_PHOTON_INTERACTION_VECTOR_H_


MLI_VECTOR_DEFINITON(mli_PhotonInteractionVector, struct mli_PhotonInteraction)

void mli_PhotonInteractionVector_print(
        const struct mli_PhotonInteractionVector *history,
        const struct mli_Scenery *scenery);

chk_rc mli_PhotonInteractionVector_time_of_flight(
        const struct mli_PhotonInteractionVector *history,
        const struct mli_Scenery *scenery,
        const double wavelength,
        double *total_time_of_flight);
#endif

/* photon_propagation */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_PHOTON_PROPAGATION_H_
#define MLI_PHOTON_PROPAGATION_H_



struct mli_PhotonPropagation {
        const struct mli_Scenery *scenery;
        struct mli_PhotonInteractionVector *history;
        struct mli_Photon *photon;
        struct mli_Prng *prng;
        uint64_t max_interactions;
};

chk_rc mli_propagate_photon(
        const struct mli_Scenery *scenery,
        struct mli_PhotonInteractionVector *history,
        struct mli_Photon *photon,
        struct mli_Prng *prng,
        const uint64_t max_interactions);
chk_rc mli_propagate_photon_work_on_causal_intersection(
        struct mli_PhotonPropagation *env);
chk_rc mli_propagate_photon_distance_until_absorption(
        const struct mli_Func *absorption_in_medium_passing_through,
        const double wavelength,
        struct mli_Prng *prng,
        double *distance_until_absorption);
chk_rc mli_propagate_photon_interact_with_object(
        struct mli_PhotonPropagation *env,
        const struct mli_IntersectionSurfaceNormal *isec);
chk_rc mli_propagate_photon_fresnel_refraction_and_reflection(
        struct mli_PhotonPropagation *env,
        const struct mli_IntersectionSurfaceNormal *isec);
chk_rc mli_propagate_photon_probability_passing_medium_coming_from(
        const struct mli_Scenery *scenery,
        const struct mli_Photon *photon,
        const struct mli_IntersectionSurfaceNormal *isec,
        double *probability_passing);
chk_rc mli_propagate_photon_pass_boundary_layer(
        struct mli_PhotonPropagation *env,
        const struct mli_IntersectionSurfaceNormal *isec,
        const struct mli_Fresnel fresnel);
struct mli_PhotonInteraction mliPhotonInteraction_from_Intersection(
        const int64_t type,
        const struct mli_Scenery *scenery,
        const struct mli_IntersectionSurfaceNormal *isec);
chk_rc mli_propagate_photon_env(struct mli_PhotonPropagation *env);
#endif

/* ray_scenery_query */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
#ifndef MLI_RAYTRACING_RAY_SCENERY_QUERY_H_
#define MLI_RAYTRACING_RAY_SCENERY_QUERY_H_


mli_bool mli_raytracing_query_intersection(
        const struct mli_Scenery *scenery,
        const struct mli_Ray ray_root,
        struct mli_Intersection *isec);

mli_bool mli_raytracing_query_intersection_with_surface_normal(
        const struct mli_Scenery *scenery,
        const struct mli_Ray ray_root,
        struct mli_IntersectionSurfaceNormal *isecsrf);

mli_bool mli_raytracing_query_object_reference(
        const struct mli_Object *object,
        const struct mli_OcTree *object_octree,
        const struct mli_HomTraComp local2root_comp,
        const struct mli_Ray ray_root,
        struct mli_Intersection *isec);

struct mli_raytracing_QueryInnerWork {
        struct mli_Intersection *intersection;
        const struct mli_Object *object;
        struct mli_Ray ray_object;
        mli_bool has_intersection;
};

struct mli_raytracing_QueryOuterWork {
        struct mli_Intersection *intersection;
        const struct mli_Geometry *geometry;
        const struct mli_Accelerator *accelerator;
        struct mli_Ray ray_root;
};

void mli_raytracing_outer_scenery_traversal(
        void *_outer,
        const struct mli_OcTree *scenery_octree,
        const uint32_t scenery_octree_leaf_idx);

void mli_raytracing_inner_object_traversal(
        void *_inner,
        const struct mli_OcTree *object_octree,
        const uint32_t object_octree_leaf_idx);

#endif

