#include "merlict_c89.h"

#include <assert.h>
#include <ctype.h>
#include <float.h>
#include <limits.h>
#include <stdlib.h>
#include <time.h>
#include <unistd.h>

/* chk_debug */
/* --------- */

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

/* mliAABB */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliAABB mliAABB_set(const struct mliVec lower, const struct mliVec upper)
{
        struct mliAABB a;
        a.lower = lower;
        a.upper = upper;
        return a;
}

struct mliAABB mliAABB_outermost(const struct mliAABB a, const struct mliAABB b)
{
        struct mliAABB c;
        c.lower.x = MLI_MIN2(a.lower.x, b.lower.x);
        c.lower.y = MLI_MIN2(a.lower.y, b.lower.y);
        c.lower.z = MLI_MIN2(a.lower.z, b.lower.z);
        c.upper.x = MLI_MAX2(a.upper.x, b.upper.x);
        c.upper.y = MLI_MAX2(a.upper.y, b.upper.y);
        c.upper.z = MLI_MAX2(a.upper.z, b.upper.z);
        return c;
}

struct mliVec mliAABB_center(const struct mliAABB a)
{
        struct mliVec sum = mliVec_add(a.upper, a.lower);
        return mliVec_multiply(sum, .5);
}

int mliAABB_valid(const struct mliAABB a)
{
        chk_msg(!MLI_IS_NAN(a.lower.x), "aabb.lower.x is 'nan'.");
        chk_msg(!MLI_IS_NAN(a.lower.y), "aabb.lower.y is 'nan'.");
        chk_msg(!MLI_IS_NAN(a.lower.z), "aabb.lower.z is 'nan'.");

        chk_msg(!MLI_IS_NAN(a.upper.x), "aabb.upper.x is 'nan'.");
        chk_msg(!MLI_IS_NAN(a.upper.y), "aabb.upper.y is 'nan'.");
        chk_msg(!MLI_IS_NAN(a.upper.z), "aabb.upper.z is 'nan'.");

        chk_msg(a.lower.x <= a.upper.x, "Expected lower.x <= upper.x");
        chk_msg(a.lower.y <= a.upper.y, "Expected lower.y <= upper.y");
        chk_msg(a.lower.z <= a.upper.z, "Expected lower.z <= upper.z");
        return 1;
error:
        return 0;
}

int mliAABB_equal(const struct mliAABB a, const struct mliAABB b)
{
        chk_msg(mliVec_equal(a.lower, b.lower),
                "Expected 'lower'-corner to be equal.");
        chk_msg(mliVec_equal(a.upper, b.upper),
                "Expected 'upper'-corner to be equal.");
        return 1;
error:
        return 0;
}

int mliAABB_is_overlapping(const struct mliAABB a, const struct mliAABB b)
{
        const int over_x = (a.upper.x >= b.lower.x) && (b.upper.x >= a.lower.x);
        const int over_y = (a.upper.y >= b.lower.y) && (b.upper.y >= a.lower.y);
        const int over_z = (a.upper.z >= b.lower.z) && (b.upper.z >= a.lower.z);
        return (over_x && over_y) && over_z;
}

/* mliAccelerator */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliAccelerator mliAccelerator_init(void)
{
        struct mliAccelerator accel;

        accel.num_objects = 0u;
        accel.object_octrees = NULL;

        accel.num_robjects = 0u;
        accel.robject_aabbs = NULL;

        accel.scenery_octree = mliOcTree_init();

        return accel;
}

void mliAccelerator_free(struct mliAccelerator *accel)
{
        uint32_t obj;

        mliOcTree_free(&accel->scenery_octree);

        for (obj = 0; obj < accel->num_objects; obj++) {
                mliOcTree_free(&accel->object_octrees[obj]);
        }
        free(accel->object_octrees);

        free(accel->robject_aabbs);
        (*accel) = mliAccelerator_init();
}

int mliAccelerator_malloc(
        struct mliAccelerator *accel,
        const uint32_t num_objects,
        const uint32_t num_robjects)
{
        uint32_t obj, rob;
        mliAccelerator_free(accel);

        accel->num_objects = num_objects;
        chk_malloc(accel->object_octrees, struct mliOcTree, accel->num_objects);
        for (obj = 0; obj < accel->num_objects; obj++) {
                accel->object_octrees[obj] = mliOcTree_init();
        }

        accel->num_robjects = num_robjects;
        chk_malloc(accel->robject_aabbs, struct mliAABB, accel->num_robjects);
        for (rob = 0; rob < accel->num_robjects; rob++) {
                accel->robject_aabbs[rob] = mliAABB_set(
                        mliVec_init(0.0, 0.0, 0.0), mliVec_init(0.0, 0.0, 0.0));
        }

        return 1;
error:
        return 0;
}

int mliAccelerator_set_robject_aabbs(
        struct mliAccelerator *accel,
        const struct mliGeometry *geometry)
{
        uint32_t rob;
        chk_msg(accel->num_robjects == geometry->num_robjects,
                "Expected num_robjects to be equal, but its not.");

        for (rob = 0; rob < accel->num_robjects; rob++) {
                const uint32_t robject = geometry->robjects[rob];
                accel->robject_aabbs[rob] = mliObject_aabb(
                        &geometry->objects[robject],
                        mliHomTra_from_compact(geometry->robject2root[rob]));
        }
        return 1;
error:
        return 0;
}

int mliAccelerator_set_object_octrees(
        struct mliAccelerator *accel,
        const struct mliGeometry *geometry)
{
        uint32_t obj;
        chk_msg(accel->num_objects == geometry->num_objects,
                "Expected num_objects to be equal, but its not.");

        for (obj = 0; obj < accel->num_objects; obj++) {
                chk_msg(mliOcTree_malloc_from_object_wavefront(
                                &accel->object_octrees[obj],
                                &geometry->objects[obj]),
                        "Failed to setup mliOctree for object-wavefront.");
        }

        return 1;
error:
        return 0;
}

int mliAccelerator_malloc_from_Geometry(
        struct mliAccelerator *accel,
        const struct mliGeometry *geometry)
{
        struct mliAABB outermost_aabb;
        struct mliGeometryAndAccelerator accgeo;
        accgeo.accelerator = accel;
        accgeo.geometry = geometry;

        chk_msg(mliAccelerator_malloc(
                        accel, geometry->num_objects, geometry->num_robjects),
                "Failed to malloc mliAccelerator from mliGeometry's "
                "num_robjects");

        chk_msg(mliAccelerator_set_robject_aabbs(accel, geometry),
                "Failed to set AABBs of robjects.");

        chk_msg(mliAccelerator_set_object_octrees(accel, geometry),
                "Failed to setup object octrees.");

        outermost_aabb = mliAccelerator_outermost_aabb(accel);

        chk_msg(mliOcTree_malloc_from_Geometry(
                        &accel->scenery_octree, &accgeo, outermost_aabb),
                "Failed to set up octree across all robjects in geometry.");

        return 1;
error:
        return 0;
}

void mliAccelerator_info_fprint(FILE *f, const struct mliAccelerator *accel)
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

        for (rob = 0; rob < accel->num_robjects; rob++) {
                fprintf(f, "    ");
                fprintf(f, "%5d ", rob);
                fprintf(f, "%9.1f ", accel->robject_aabbs[rob].lower.x);
                fprintf(f, "%9.1f ", accel->robject_aabbs[rob].lower.y);
                fprintf(f, "%9.1f ", accel->robject_aabbs[rob].lower.z);
                fprintf(f, "%9.1f ", accel->robject_aabbs[rob].upper.x);
                fprintf(f, "%9.1f ", accel->robject_aabbs[rob].upper.y);
                fprintf(f, "%9.1f ", accel->robject_aabbs[rob].upper.z);
                fprintf(f, "\n");
        }
}

struct mliAABB mliAccelerator_outermost_aabb(const struct mliAccelerator *accel)
{
        uint32_t rob;
        struct mliAABB aabb;
        if (accel->num_robjects == 0) {
                aabb.lower = mliVec_init(MLI_NAN, MLI_NAN, MLI_NAN);
                aabb.upper = mliVec_init(MLI_NAN, MLI_NAN, MLI_NAN);
                return aabb;
        }
        aabb.lower = accel->robject_aabbs[0].lower;
        aabb.upper = accel->robject_aabbs[0].upper;
        for (rob = 0; rob < accel->num_robjects; rob++) {
                aabb = mliAABB_outermost(aabb, accel->robject_aabbs[rob]);
        }
        return aabb;
}

/* mliAccelerator_equal */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliAccelerator_equal(
        const struct mliAccelerator *a,
        const struct mliAccelerator *b)
{
        uint32_t i = 0u;
        chk_msg(a->num_objects == b->num_objects,
                "Expected num_objects to be equal.");
        for (i = 0; i < a->num_objects; i++) {
                chk_msg(mliOcTree_equal(
                                &a->object_octrees[i], &b->object_octrees[i]),
                        "Expected object_octrees[i] to be equal.");
        }

        chk_msg(a->num_robjects == b->num_robjects,
                "Expected num_robjects to be equal.");
        for (i = 0; i < a->num_robjects; i++) {
                chk_msg(mliAABB_equal(a->robject_aabbs[i], b->robject_aabbs[i]),
                        "Expected robject_aabbs[i] to be equal.");
        }

        chk_msg(mliOcTree_equal(&a->scenery_octree, &b->scenery_octree),
                "Expected scenery_octree to be equal.");

        return 1;
error:
        return 0;
}

/* mliAccelerator_serialize */
/* ------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliAccelerator_fwrite(const struct mliAccelerator *accel, FILE *f)
{
        uint64_t i = 0;

        /* magic identifier */
        struct mliMagicId magic = mliMagicId_init();
        chk(mliMagicId_set(&magic, "mliAccelerator"));
        chk_fwrite(&magic, sizeof(struct mliMagicId), 1u, f);

        /* capacity */
        chk_fwrite(&accel->num_objects, sizeof(uint32_t), 1, f);
        chk_fwrite(&accel->num_robjects, sizeof(uint32_t), 1, f);

        for (i = 0; i < accel->num_objects; i++) {
                mliOcTree_fwrite(&accel->object_octrees[i], f);
        }
        chk_fwrite(
                accel->robject_aabbs,
                sizeof(struct mliAABB),
                accel->num_robjects,
                f);
        mliOcTree_fwrite(&accel->scenery_octree, f);

        return 1;
error:
        return 0;
}

int mliAccelerator_malloc_fread(struct mliAccelerator *accel, FILE *f)
{
        uint64_t i = 0u;
        struct mliMagicId magic;

        uint32_t num_robjects = 0u;
        uint32_t num_objects = 0u;

        /* magic identifier */
        chk_fread(&magic, sizeof(struct mliMagicId), 1u, f);
        chk(mliMagicId_has_word(&magic, "mliAccelerator"));
        mliMagicId_warn_version(&magic);

        /* capacity */
        chk_fread(&num_objects, sizeof(uint32_t), 1u, f);
        chk_fread(&num_robjects, sizeof(uint32_t), 1u, f);

        /* malloc */
        chk_mem(mliAccelerator_malloc(accel, num_objects, num_robjects));

        for (i = 0; i < accel->num_objects; i++) {
                chk_mem(mliOcTree_malloc_fread(&accel->object_octrees[i], f));
        }

        chk_fread(
                accel->robject_aabbs,
                sizeof(struct mliAABB),
                accel->num_robjects,
                f);

        chk_mem(mliOcTree_malloc_fread(&accel->scenery_octree, f));

        return 1;
error:
        return 0;
}

/* mliAccelerator_valid */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliAccelerator_valid(const struct mliAccelerator *accel)
{
        uint32_t i = 0u;
        for (i = 0u; i < accel->num_objects; i++) {
                chk_msg(mliOcTree_valid(&accel->object_octrees[i]),
                        "Expected object_octrees[i] to be valid.");
        }
        for (i = 0u; i < accel->num_robjects; i++) {
                chk_msg(mliAABB_valid(accel->robject_aabbs[i]),
                        "Expected robject_aabbs[i] to be valid.");
        }
        chk_msg(mliOcTree_valid(&accel->scenery_octree),
                "Expected scenery_octree to be valid.");
        return 1;
error:
        return 0;
}

int mliAccelerator_valid_wrt_Geometry(
        const struct mliAccelerator *accel,
        const struct mliGeometry *geometry)
{
        uint32_t i = 0u;
        chk_msg(accel->num_objects == geometry->num_objects,
                "Expected num_objects to be equal.");
        for (i = 0u; i < accel->num_objects; i++) {
                chk_msg(mliOcTree_valid_wrt_links(
                                &accel->object_octrees[i],
                                geometry->objects[i].num_faces),
                        "Expected object_octrees[i] to be valid w.r.t. "
                        "the object's num_faces.");
        }

        chk_msg(accel->num_robjects == geometry->num_robjects,
                "Expected num_robjects to be equal.");
        chk_msg(mliOcTree_valid_wrt_links(
                        &accel->scenery_octree, geometry->num_robjects),
                "Expected scenery_octree to be valid w.r.t. to "
                "geometry's num_robjects.");
        return 1;
error:
        return 0;
}

/* mliApertureCamera */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliApertureCamera mliApertureCamera_init(void)
{
        const double mm = 1e-3;
        struct mliApertureCamera apcam;
        apcam.focal_length = 50.0 * mm;
        apcam.aperture_radius = apcam.focal_length / 2.0;
        apcam.image_sensor_distance = apcam.focal_length;
        apcam.image_sensor_width_x = 36.0 * mm;
        apcam.image_sensor_width_y = 24.0 * mm;
        return apcam;
}

struct mliVec mliApertureCamera_pixel_center_on_image_sensor_plane(
        const double image_sensor_width_x,
        const double image_sensor_width_y,
        const double image_sensor_distance,
        const uint64_t num_pixel_x,
        const uint64_t num_pixel_y,
        const uint64_t pixel_x,
        const uint64_t pixel_y)
{
        struct mliVec pixel_center;
        pixel_center.x = -1.0 * mli_bin_center_in_linear_space(
                                        -0.5 * image_sensor_width_x,
                                        +0.5 * image_sensor_width_x,
                                        num_pixel_x,
                                        pixel_x);
        pixel_center.y = -1.0 * mli_bin_center_in_linear_space(
                                        -0.5 * image_sensor_width_y,
                                        +0.5 * image_sensor_width_y,
                                        num_pixel_y,
                                        pixel_y);
        pixel_center.z = -image_sensor_distance;
        return pixel_center;
}

struct mliVec mliApertureCamera_pixel_support_on_image_sensor_plane(
        const double image_sensor_width_x,
        const double image_sensor_width_y,
        const double image_sensor_distance,
        const uint64_t num_pixel_x,
        const uint64_t num_pixel_y,
        const uint64_t pixel_x,
        const uint64_t pixel_y,
        struct mliPrng *prng)
{
        double pixel_bin_width_x = image_sensor_width_x / (double)num_pixel_x;
        double pixel_bin_width_y = image_sensor_width_y / (double)num_pixel_y;
        struct mliVec support =
                mliApertureCamera_pixel_center_on_image_sensor_plane(
                        image_sensor_width_x,
                        image_sensor_width_y,
                        image_sensor_distance,
                        num_pixel_x,
                        num_pixel_y,
                        pixel_x,
                        pixel_y);
        double rnd_x =
                (mli_random_uniform(prng) * pixel_bin_width_x -
                 0.5 * pixel_bin_width_x);
        double rnd_y =
                (mli_random_uniform(prng) * pixel_bin_width_y -
                 0.5 * pixel_bin_width_y);
        support.x = support.x + rnd_x;
        support.y = support.y + rnd_y;
        return support;
}

struct mliVec mliApertureCamera_get_object_point(
        const double focal_length,
        const struct mliVec pixel_support)
{
        const double object_distance =
                mli_thin_lens_get_object_given_focal_and_image(
                        focal_length, -1.0 * pixel_support.z);
        const double scaleing = object_distance / pixel_support.z;
        return mliVec_init(
                scaleing * pixel_support.x,
                scaleing * pixel_support.y,
                object_distance);
}

struct mliVec mliApertureCamera_ray_support_on_aperture(
        const double aperture_radius,
        struct mliPrng *prng)
{
        /* use a perfect disc as aperture */
        return mli_random_position_on_disc(aperture_radius, prng);
}

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

double mliApertureCamera_focal_length_given_field_of_view_and_sensor_width(
        const double field_of_view,
        const double image_sensor_width)
{
        const double image_sensor_radius = 0.5 * image_sensor_width;
        const double fov_opening_angle = 0.5 * field_of_view;
        return image_sensor_radius / tan(fov_opening_angle);
}

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
        struct mliPrng *prng)
{
        struct mliVec direction;
        struct mliVec aperture_support =
                mliApertureCamera_ray_support_on_aperture(
                        aperture_radius, prng);

        struct mliVec image_sensor_support =
                (mliApertureCamera_pixel_support_on_image_sensor_plane(
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
                direction = mliVec_multiply(image_sensor_support, -1.0);
        } else {
                struct mliVec object_point = mliApertureCamera_get_object_point(
                        focal_length, image_sensor_support);

                direction = mliVec_substract(object_point, aperture_support);
        }

        return mliRay_set(aperture_support, direction);
}

void mliApertureCamera_aquire_pixels(
        const struct mliApertureCamera camera,
        const struct mliImage *image,
        const struct mliHomTraComp camera2root_comp,
        const struct mliScenery *scenery,
        const struct mliPixels *pixels_to_do,
        struct mliImage *colors,
        const struct mliTracerConfig *tracer_config,
        struct mliPrng *prng)
{
        uint64_t i;
        struct mliHomTra camera2root = mliHomTra_from_compact(camera2root_comp);
        for (i = 0; i < pixels_to_do->num_pixels_to_do; i++) {
                struct mliRay ray_wrt_camera =
                        mliApertureCamera_get_ray_for_pixel(
                                camera.focal_length,
                                camera.aperture_radius,
                                camera.image_sensor_distance,
                                camera.image_sensor_width_x,
                                camera.image_sensor_width_y,
                                image->num_cols,
                                image->num_rows,
                                pixels_to_do->pixels[i].col,
                                pixels_to_do->pixels[i].row,
                                prng);

                struct mliRay ray_wrt_root =
                        mliHomTra_ray(&camera2root, ray_wrt_camera);

                struct mliColor set_color =
                        mli_trace(scenery, ray_wrt_root, tracer_config, prng);

                mliImage_set(colors, i, 0u, set_color);
        }

        return;
}

void mliApertureCamera_assign_pixel_colors_to_sum_and_exposure_image(
        const struct mliPixels *pixels,
        const struct mliImage *colors,
        struct mliImage *sum_image,
        struct mliImage *exposure_image)
{
        uint64_t pix;
        for (pix = 0; pix < pixels->num_pixels_to_do; pix++) {
                const uint64_t idx = mliImage_idx(
                        sum_image,
                        pixels->pixels[pix].col,
                        pixels->pixels[pix].row);
                sum_image->raw[idx].r += colors->raw[pix].r;
                sum_image->raw[idx].g += colors->raw[pix].g;
                sum_image->raw[idx].b += colors->raw[pix].b;

                exposure_image->raw[idx].r += 1.;
                exposure_image->raw[idx].g += 1.;
                exposure_image->raw[idx].b += 1.;
        }
}

int mliApertureCamera_render_image(
        const struct mliApertureCamera camera,
        const struct mliHomTraComp camera2root_comp,
        const struct mliScenery *scenery,
        struct mliImage *image,
        const struct mliTracerConfig *tracer_config,
        struct mliPrng *prng)
{
        float noise_threshold = 0.05 * 255.0;
        uint64_t MAX_ITERATIONS = 128;
        uint64_t iteration = 0;

        struct mliColor zero_color = mliColor_set(0.0, 0.0, 0.0);
        struct mliImage sum_image = mliImage_init();
        struct mliImage exposure_image = mliImage_init();
        struct mliImage to_do_image = mliImage_init();
        struct mliImage sobel_image = mliImage_init();
        struct mliImage previous_sobel_image = mliImage_init();
        struct mliImage diff_image = mliImage_init();
        struct mliImage colors = mliImage_init();
        struct mliPixels pixels_to_do = mliPixels_init();

        chk_msg(mliImage_malloc(&sum_image, image->num_cols, image->num_rows),
                "Failed to malloc sum_image.");
        chk_msg(mliImage_malloc(
                        &exposure_image, image->num_cols, image->num_rows),
                "Failed to malloc exposure_image.");
        chk_msg(mliImage_malloc(&to_do_image, image->num_cols, image->num_rows),
                "Failed to malloc to_do_image.");
        chk_msg(mliImage_malloc(&sobel_image, image->num_cols, image->num_rows),
                "Failed to malloc sobel_image.");
        chk_msg(mliImage_malloc(
                        &previous_sobel_image,
                        image->num_cols,
                        image->num_rows),
                "Failed to malloc previous_sobel_image.");
        chk_msg(mliImage_malloc(&diff_image, image->num_cols, image->num_rows),
                "Failed to malloc diff_image.");
        chk_msg(mliImage_malloc(&colors, image->num_cols * image->num_rows, 1),
                "Failed to malloc colors.");

        mliImage_set_all_pixel(image, zero_color);
        mliImage_set_all_pixel(&sum_image, zero_color);
        mliImage_set_all_pixel(&exposure_image, zero_color);
        mliImage_set_all_pixel(&to_do_image, zero_color);
        mliImage_set_all_pixel(&sobel_image, zero_color);
        mliImage_set_all_pixel(&previous_sobel_image, zero_color);
        mliImage_set_all_pixel(&colors, zero_color);

        chk_msg(mliPixels_malloc(
                        &pixels_to_do, image->num_cols * image->num_rows),
                "Failed to malloc pixels_to_do.");

        /*
        initial image
        =============
        */
        mliPixels_set_all_from_image(&pixels_to_do, image);
        pixels_to_do.num_pixels_to_do = image->num_cols * image->num_rows;

        mliApertureCamera_aquire_pixels(
                camera,
                image,
                camera2root_comp,
                scenery,
                &pixels_to_do,
                &colors,
                tracer_config,
                prng);

        mliApertureCamera_assign_pixel_colors_to_sum_and_exposure_image(
                &pixels_to_do, &colors, &sum_image, &exposure_image);

        mliImage_divide_pixelwise(&sum_image, &exposure_image, image);

        mliImage_sobel(image, &sobel_image);

        mliImage_luminance_threshold_dilatation(
                &sobel_image, 128.0, &to_do_image);

        /*printf("\n");*/
        while (1) {
                if (iteration >= MAX_ITERATIONS)
                        break;

                mliPixels_above_threshold(&to_do_image, 0.5, &pixels_to_do);

                if (pixels_to_do.num_pixels_to_do <
                    image->num_rows * image->num_cols / 100.0)
                        break;

                /*printf("loop %3u / %3u, %d,%03d pixel left\n",
                       (uint32_t)iteration + 1,
                       (uint32_t)MAX_ITERATIONS,
                       pixels_to_do.num_pixels_to_do / 1000,
                       pixels_to_do.num_pixels_to_do % 1000);*/
                mliApertureCamera_aquire_pixels(
                        camera,
                        image,
                        camera2root_comp,
                        scenery,
                        &pixels_to_do,
                        &colors,
                        tracer_config,
                        prng);

                mliApertureCamera_assign_pixel_colors_to_sum_and_exposure_image(
                        &pixels_to_do, &colors, &sum_image, &exposure_image);

                mliImage_divide_pixelwise(&sum_image, &exposure_image, image);

                mliImage_copy(&sobel_image, &previous_sobel_image);

                mliImage_sobel(image, &sobel_image);
                mliImage_fabs_difference(
                        &previous_sobel_image, &sobel_image, &diff_image);

                mliImage_set_all_pixel(&to_do_image, zero_color);

                mliImage_luminance_threshold_dilatation(
                        &diff_image, noise_threshold, &to_do_image);

                iteration += 1;
        }

        mliImage_free(&sum_image);
        mliImage_free(&exposure_image);
        mliImage_free(&to_do_image);
        mliImage_free(&sobel_image);
        mliImage_free(&previous_sobel_image);
        mliImage_free(&diff_image);
        mliImage_free(&colors);
        mliPixels_free(&pixels_to_do);
        return 1;
error:
        mliImage_free(&sum_image);
        mliImage_free(&exposure_image);
        mliImage_free(&to_do_image);
        mliImage_free(&sobel_image);
        mliImage_free(&previous_sobel_image);
        mliImage_free(&diff_image);
        mliImage_free(&colors);
        mliPixels_free(&pixels_to_do);
        return 0;
}

/* mliArchive */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

MLIDYNARRAY_IMPLEMENTATION(mli, TextFiles, struct mliStr)

struct mliArchive mliArchive_init(void)
{
        struct mliArchive arc;
        arc.textfiles = mliDynTextFiles_init();
        arc.filenames = mliDynMap_init();
        return arc;
}

void mliArchive_free(struct mliArchive *arc)
{
        uint64_t i;
        for (i = 0; i < arc->textfiles.size; i++) {
                mliStr_free(&arc->textfiles.array[i]);
        }
        mliDynTextFiles_free(&arc->textfiles);
        mliDynMap_free(&arc->filenames);
        (*arc) = mliArchive_init();
}

int mliArchive_malloc(struct mliArchive *arc)
{
        mliArchive_free(arc);
        chk(mliDynTextFiles_malloc(&arc->textfiles, 0u));
        chk(mliDynMap_malloc(&arc->filenames, 0u));
        return 1;
error:
        return 0;
}

int mliArchive_push_back(
        struct mliArchive *arc,
        const struct mliStr *filename,
        const struct mliStr *payload)
{
        uint64_t next;
        chk_msg(filename->length < MLI_TAR_NAME_LENGTH,
                "Expected shorter filename.");
        next = arc->filenames.size;

        /* filename */
        /* ======== */
        chk_msg(mliDynMap_insert(&arc->filenames, filename->cstr, next),
                "Can not insert key.");

        /* payload */
        /* ======= */
        chk_msg(mliDynTextFiles_push_back(&arc->textfiles, mliStr_init()),
                "Can not push back mliStr.");
        chk_msg(mliStr_malloc_copy(&arc->textfiles.array[next], payload),
                "Can not copy payload.");
        return 1;
error:
        return 0;
}

int mliArchive_malloc_fread(struct mliArchive *arc, FILE *f)
{
        struct mliTar tar = mliTar_init();
        struct mliTarHeader tarh = mliTarHeader_init();
        struct mliStr payload = mliStr_init();
        struct mliStr filename = mliStr_init();

        char tarh_name[MLI_TAR_NAME_LENGTH] = {'\0'};

        chk_msg(mliArchive_malloc(arc), "Can not malloc archive.");
        chk_msg(mliTar_read_begin(&tar, f), "Can't begin tar.");

        while (mliTar_read_header(&tar, &tarh)) {

                chk(mliStr_malloc_cstr(&filename, tarh.name));
                chk(mliStr_strip(&filename, &filename));
                chk(mli_path_strip_this_dir(&filename, &filename));

                chk_msg(mliStr_malloc(&payload, tarh.size),
                        "Can not allocate payload.");
                chk_msg(mliTar_read_data(&tar, (void *)payload.cstr, tarh.size),
                        "Failed to read payload from tar into payload.");
                chk_msg(mliStr_convert_line_break_CRLF_CR_to_LF(
                                &payload, &payload),
                        "Failed to replace CRLF and CR linebreaks.");
                chk_msg(mli_cstr_assert_only_NUL_LF_TAB_controls(payload.cstr),
                        "Did not expect control codes other than "
                        "('\\n', '\\t', '\\0') in textfiles.");
                chk_msg(mliArchive_push_back(arc, &filename, &payload),
                        "Can not push back file into archive.");
        }

        chk_msg(mliTar_read_finalize(&tar), "Can't finalize reading tar.");
        mliStr_free(&payload);
        mliStr_free(&filename);
        return 1;
error:
        fprintf(stderr, "tar->filename: '%s'.\n", tarh_name);
        mliStr_free(&payload);
        mliStr_free(&filename);
        mliArchive_free(arc);
        return 0;
}

int mliArchive_malloc_from_path(struct mliArchive *arc, const char *path)
{
        FILE *f = fopen(path, "rb");
        chk_msgf(f != NULL, ("Can't open path '%s'.", path))
                chk_msg(mliArchive_malloc_fread(arc, f),
                        "Can't fread Archive from file.");
        fclose(f);
        return 1;
error:
        return 0;
}

int mliArchive_has(const struct mliArchive *arc, const char *filename)
{
        return mliDynMap_has(&arc->filenames, filename);
}

int mliArchive_get(
        const struct mliArchive *arc,
        const char *filename,
        struct mliStr **str)
{
        uint64_t idx;
        chk(mliDynMap_find(&arc->filenames, filename, &idx));
        (*str) = &arc->textfiles.array[idx];
        return 1;
error:
        return 0;
}

int mliArchive_get_malloc_json(
        const struct mliArchive *arc,
        const char *filename,
        struct mliJson *json)
{
        struct mliStr *text = NULL;

        chk_msg(mliArchive_get(arc, filename, &text),
                "Can not find requested file in archive.");

        chk_msg(mliJson_malloc_from_cstr(json, (char *)text->cstr),
                "Can not parse requested json.");

        return 1;
error:
        return 0;
}

uint64_t mliArchive_num(const struct mliArchive *arc)
{
        return arc->filenames.size;
}

void mliArchive_info_fprint(FILE *f, const struct mliArchive *arc)
{
        uint64_t i;
        for (i = 0; i < arc->textfiles.size; i++) {
                struct mliDynMapItem *map_item = &arc->filenames.array[i];
                fprintf(f,
                        "%u: %s, %u\n",
                        (uint32_t)i,
                        map_item->key,
                        (uint32_t)arc->textfiles.array[i].length);
        }
}

void mliArchive_mask_filename_prefix_sufix(
        const struct mliArchive *arc,
        uint64_t *mask,
        const char *prefix,
        const char *sufix)
{
        uint64_t i = 0u;
        uint64_t match = 0u;
        for (i = 0; i < arc->textfiles.size; i++) {
                struct mliDynMapItem *map_item = &arc->filenames.array[i];

                match = mli_cstr_has_prefix_suffix(
                        map_item->key, prefix, sufix);

                if (match) {
                        mask[i] = 1;
                } else {
                        mask[i] = 0;
                }
        }
}

uint64_t mliArchive_num_filename_prefix_sufix(
        const struct mliArchive *arc,
        const char *prefix,
        const char *sufix)
{
        uint64_t i = 0;
        uint64_t match;
        uint64_t num_matches = 0;
        for (i = 0; i < arc->textfiles.size; i++) {
                struct mliDynMapItem *map_item = &arc->filenames.array[i];

                match = mli_cstr_has_prefix_suffix(
                        map_item->key, prefix, sufix);

                if (match) {
                        num_matches++;
                }
        }
        return num_matches;
}

/* mliAtmosphere */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
/* based on:
 * 2009-2016 Scratchapixel. Distributed under the terms of the
 * GNU General Public License as published by the Free Software Foundation,
 * either version 3 of the License, or (at your option) any later version.
 */


struct mliAtmosphere mliAtmosphere_init(void)
{
        struct mliAtmosphere atm;
        atm.sunLatitude = 0.0;
        atm.sunHourAngle = 12.0;
        mliAtmosphere_set_sun_direction(
                &atm, atm.sunLatitude, atm.sunHourAngle);

        atm.sunDistance = 1.5e11;
        atm.sunRadius = 7e8;

        atm.altitude = 2300.0;
        atm.earthRadius = 6360e3;
        atm.atmosphereRadius = atm.earthRadius + 60e3;

        /* The height for the density to drop by 1 over e */
        atm.Height_Rayleigh = 7994.0;
        atm.Height_Mie = 1200.0;

        atm.beta_Rayleigh = mliColor_set(3.8e-6, 13.5e-6, 33.1e-6);
        atm.beta_Mie = mliColor_multiply(mliColor_set(1.0, 1.0, 1.0), 41e-6);

        atm.numSamples = 16;
        atm.numSamplesLight = 8;

        atm.power = 3000.0;

        return atm;
}

void mliAtmosphere_set_sun_direction(
        struct mliAtmosphere *atmosphere,
        const double sunLatitude,
        const double sunHourAngle)
{
        atmosphere->sunHourAngle = sunHourAngle;
        atmosphere->sunLatitude = sunLatitude;

        {
                const double hours_rad =
                        MLI_PI + 2.0 * MLI_PI * atmosphere->sunHourAngle / 24.0;

                const struct mliHomTraComp tc_latitude = mliHomTraComp_set(
                        mliVec_init(0.0, 0.0, 0.0),
                        mliQuaternion_set_tait_bryan(
                                atmosphere->sunLatitude, 0.0, 0.0));
                const struct mliHomTraComp tc_hour = mliHomTraComp_set(
                        mliVec_init(0.0, 0.0, 0.0),
                        mliQuaternion_set_tait_bryan(0.0, hours_rad, 0.0));
                const struct mliHomTraComp tc =
                        mliHomTraComp_sequence(tc_latitude, tc_hour);
                const struct mliHomTra t = mliHomTra_from_compact(tc);
                const struct mliVec zenith = mliVec_init(0.0, 0.0, 1.0);
                atmosphere->sunDirection = mliHomTra_dir(&t, zenith);
        }
}

struct mliColor mliAtmosphere_compute_depth(
        const struct mliAtmosphere *atmosphere,
        const struct mliVec orig,
        const struct mliVec dir,
        double tmin,
        double tmax)
{
        uint64_t i, j;
        const double segmentLength = (tmax - tmin) / atmosphere->numSamples;
        double tCurrent = tmin;
        struct mliColor sumR = mliColor_set(0.0, 0.0, 0.0);
        struct mliColor sumM = mliColor_set(0.0, 0.0, 0.0);
        /* mie and rayleigh contribution */

        double opticalDepthR = 0.0;
        double opticalDepthM = 0.0;
        /*
         * mu in the paper which is the cosine of the angle between the sun
         * direction and the ray direction
         */
        const double mu = mliVec_dot(dir, atmosphere->sunDirection);
        const double phaseR = 3.f / (16.f * MLI_PI) * (1.0 + mu * mu);
        const double g = 0.76f;
        const double phaseM =
                3.f / (8.f * MLI_PI) *
                (((1.f - g * g) * (1.f + mu * mu)) /
                 ((2.f + g * g) * pow(1.f + g * g - 2.f * g * mu, 1.5f)));

        for (i = 0; i < atmosphere->numSamples; ++i) {
                const struct mliVec samplePosition = mliVec_add(
                        orig,
                        mliVec_multiply(
                                dir, (tCurrent + segmentLength * 0.5f)));
                const double height =
                        (mliVec_norm(samplePosition) - atmosphere->earthRadius);
                /* compute optical depth for light */
                const double hr = segmentLength *
                                  exp(-height / atmosphere->Height_Rayleigh);
                const double hm =
                        segmentLength * exp(-height / atmosphere->Height_Mie);
                double t0Light = 0;
                double t1Light = 0;
                double segmentLengthLight = 0;
                double tCurrentLight = 0;
                double opticalDepthLightR = 0;
                double opticalDepthLightM = 0;

                opticalDepthR += hr;
                opticalDepthM += hm;
                /* light optical depth */
                mliRay_sphere_intersection(
                        samplePosition,
                        atmosphere->sunDirection,
                        atmosphere->atmosphereRadius,
                        &t0Light,
                        &t1Light);
                segmentLengthLight = t1Light / atmosphere->numSamplesLight;

                for (j = 0; j < atmosphere->numSamplesLight; ++j) {
                        const struct mliVec samplePositionLight = mliVec_add(
                                samplePosition,
                                mliVec_multiply(
                                        atmosphere->sunDirection,
                                        tCurrentLight +
                                                0.5 * segmentLengthLight));
                        const double heightLight =
                                mliVec_norm(samplePositionLight) -
                                atmosphere->earthRadius;

                        if (heightLight < 0)
                                break;

                        opticalDepthLightR +=
                                segmentLengthLight *
                                exp(-heightLight / atmosphere->Height_Rayleigh);
                        opticalDepthLightM +=
                                segmentLengthLight *
                                exp(-heightLight / atmosphere->Height_Mie);
                        tCurrentLight += segmentLengthLight;
                }

                if (j == atmosphere->numSamplesLight) {
                        const struct mliColor tau = mliColor_add(
                                mliColor_multiply(
                                        atmosphere->beta_Rayleigh,
                                        opticalDepthR + opticalDepthLightR),
                                mliColor_multiply(
                                        atmosphere->beta_Mie,
                                        1.1 * (opticalDepthM +
                                               opticalDepthLightM)));
                        const struct mliColor attenuation = mliColor_set(
                                exp(-tau.r), exp(-tau.g), exp(-tau.b));
                        sumR = mliColor_add(
                                sumR, mliColor_multiply(attenuation, hr));
                        sumM = mliColor_add(
                                sumM, mliColor_multiply(attenuation, hm));
                }
                tCurrent += segmentLength;
        }

        return mliColor_multiply(
                mliColor_add(
                        mliColor_multiply(
                                mliColor_multiply_elementwise(
                                        sumR, atmosphere->beta_Rayleigh),
                                phaseR),
                        mliColor_multiply(
                                mliColor_multiply_elementwise(
                                        sumM, atmosphere->beta_Mie),
                                phaseM)),
                atmosphere->power);
}

struct mliColor mliAtmosphere_hit_outer_atmosphere(
        const struct mliAtmosphere *atmosphere,
        const struct mliVec orig,
        const struct mliVec dir,
        double tmin,
        double tmax)
{
        double t_minus = -1.0;
        double t_plus = -1.0;
        int has_intersection = mliRay_sphere_intersection(
                orig, dir, atmosphere->atmosphereRadius, &t_minus, &t_plus);

        if (!has_intersection || t_plus < 0.0) {
                return mliColor_set(0.0, 0.0, 0.0);
        }
        if (t_minus > tmin && t_minus > 0) {
                tmin = t_minus;
        }
        if (t_plus < tmax) {
                tmax = t_plus;
        }

        return mliAtmosphere_compute_depth(atmosphere, orig, dir, tmin, tmax);
}

struct mliColor mliAtmosphere_hit_earth_body(
        const struct mliAtmosphere *atmosphere,
        const struct mliVec orig,
        const struct mliVec dir)
{
        double t_minus = DBL_MAX;
        double t_plus = DBL_MAX;
        double t_max = DBL_MAX;
        int intersects_earth_body = mliRay_sphere_intersection(
                orig, dir, atmosphere->earthRadius, &t_minus, &t_plus);

        if (intersects_earth_body && t_minus > 0) {
                t_max = MLI_MAX2(0.0, t_minus);
        }

        return mliAtmosphere_hit_outer_atmosphere(
                atmosphere, orig, dir, 0.0, t_max);
}

struct mliColor mliAtmosphere_query(
        const struct mliAtmosphere *atmosphere,
        const struct mliVec orig,
        const struct mliVec dir)
{
        struct mliVec orig_up = orig;
        orig_up.z += (atmosphere->earthRadius + atmosphere->altitude);

        return mliAtmosphere_hit_earth_body(atmosphere, orig_up, dir);
}

void mliAtmosphere_increase_latitude(
        struct mliAtmosphere *atmosphere,
        const double increment)
{
        if (atmosphere->sunLatitude + increment <= 0.5 * MLI_PI) {
                atmosphere->sunLatitude += increment;
                mliAtmosphere_set_sun_direction(
                        atmosphere,
                        atmosphere->sunLatitude,
                        atmosphere->sunHourAngle);
        }
}

void mliAtmosphere_decrease_latitude(
        struct mliAtmosphere *atmosphere,
        const double increment)
{
        if (atmosphere->sunLatitude - increment >= -0.5 * MLI_PI) {
                atmosphere->sunLatitude -= increment;
                mliAtmosphere_set_sun_direction(
                        atmosphere,
                        atmosphere->sunLatitude,
                        atmosphere->sunHourAngle);
        }
}

void mliAtmosphere_increase_hours(
        struct mliAtmosphere *atmosphere,
        const double increment)
{
        if (atmosphere->sunHourAngle + increment < 24.0) {
                atmosphere->sunHourAngle += increment;
                mliAtmosphere_set_sun_direction(
                        atmosphere,
                        atmosphere->sunLatitude,
                        atmosphere->sunHourAngle);
        }
}

void mliAtmosphere_decrease_hours(
        struct mliAtmosphere *atmosphere,
        const double increment)
{
        if (atmosphere->sunHourAngle - increment >= 0.0) {
                atmosphere->sunHourAngle -= increment;
                mliAtmosphere_set_sun_direction(
                        atmosphere,
                        atmosphere->sunLatitude,
                        atmosphere->sunHourAngle);
        }
}

void mliAtmosphere_increase_altitude(
        struct mliAtmosphere *atmosphere,
        const double factor)
{
        assert(factor > 1.0);
        atmosphere->altitude *= factor;
}

void mliAtmosphere_decrease_altitude(
        struct mliAtmosphere *atmosphere,
        const double factor)
{
        assert(factor < 1.0);
        assert(factor > 0.0);
        atmosphere->altitude *= factor;
}

/* mliAtmosphere_json */
/* ------------------ */

/* Copyright 2018-2021 Sebastian Achim Mueller */

int mliAtmosphere_from_json_token(
        struct mliAtmosphere *atm,
        const struct mliJson *json,
        const uint64_t tkn)
{
        uint64_t beta_rayleigh_tkn;
        uint64_t beta_mie_tkn;

        (*atm) = mliAtmosphere_init();

        chk(mliJson_double_by_key(json, tkn, &atm->sunLatitude, "sunLatitude"));
        chk(mliJson_double_by_key(
                json, tkn, &atm->sunHourAngle, "sunHourAngle"));
        mliAtmosphere_set_sun_direction(
                atm, atm->sunLatitude, atm->sunHourAngle);

        chk(mliJson_double_by_key(json, tkn, &atm->sunDistance, "sunDistance"));
        chk_msg(atm->sunDistance > 0, "Expected atmosphere->sunDistance > 0.");
        chk(mliJson_double_by_key(json, tkn, &atm->sunRadius, "sunRadius"));
        chk_msg(atm->sunRadius > 0, "Expected atmosphere->sunRadius > 0.");

        chk(mliJson_double_by_key(json, tkn, &atm->earthRadius, "earthRadius"));
        chk_msg(atm->earthRadius > 0, "Expected atmosphere->earthRadius > 0.");
        chk(mliJson_double_by_key(
                json, tkn, &atm->atmosphereRadius, "atmosphereRadius"));
        chk_msg(atm->atmosphereRadius > atm->earthRadius,
                "Expected atmosphere->atmosphereRadius > atm->earthRadius.");

        chk(mliJson_double_by_key(
                json, tkn, &atm->Height_Rayleigh, "Height_Rayleigh"));
        chk(mliJson_double_by_key(json, tkn, &atm->Height_Mie, "Height_Mie"));

        chk(mliJson_uint64_by_key(json, tkn, &atm->numSamples, "numSamples"));
        chk_msg(atm->numSamples > 0, "Expected atmosphere->numSamples > 0.");
        chk(mliJson_uint64_by_key(
                json, tkn, &atm->numSamplesLight, "numSamplesLight"));
        chk_msg(atm->numSamplesLight > 0,
                "Expected atmosphere->numSamplesLight > 0.");

        chk(mliJson_token_by_key(
                json, tkn, "beta_Rayleigh", &beta_rayleigh_tkn));
        chk(mliColor_from_json_token(
                &atm->beta_Rayleigh, json, beta_rayleigh_tkn + 1));
        chk_msg(atm->beta_Rayleigh.r > 0.0,
                "Expected atmosphere->beta_Rayleigh.r > 0.");
        chk_msg(atm->beta_Rayleigh.g > 0.0,
                "Expected atmosphere->beta_Rayleigh.g > 0.");
        chk_msg(atm->beta_Rayleigh.b > 0.0,
                "Expected atmosphere->beta_Rayleigh.b > 0.");

        chk(mliJson_token_by_key(json, tkn, "beta_Mie", &beta_mie_tkn));
        chk(mliColor_from_json_token(&atm->beta_Mie, json, beta_mie_tkn + 1));
        chk_msg(atm->beta_Mie.r > 0.0, "Expected atmosphere->beta_Mie.r > 0.");
        chk_msg(atm->beta_Mie.g > 0.0, "Expected atmosphere->beta_Mie.g > 0.");
        chk_msg(atm->beta_Mie.b > 0.0, "Expected atmosphere->beta_Mie.b > 0.");

        chk(mliJson_double_by_key(json, tkn, &atm->power, "power"));
        chk_msg(atm->power > 0, "Expected atmosphere->power > 0.");
        chk(mliJson_double_by_key(json, tkn, &atm->altitude, "altitude"));
        chk_msg(atm->altitude > 0, "Expected atmosphere->altitude > 0.");

        return 1;
error:
        return 0;
}

/* mliBoundaryLayer */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliBoundaryLayer_equal(
        const struct mliBoundaryLayer a,
        const struct mliBoundaryLayer b)
{
        if (a.inner.surface != b.inner.surface)
                return 0;
        if (a.outer.surface != b.outer.surface)
                return 0;
        if (a.inner.medium != b.inner.medium)
                return 0;
        if (a.outer.medium != b.outer.medium)
                return 0;
        return 1;
}

void mliBoundaryLayer_print(const struct mliBoundaryLayer a)
{
        fprintf(stderr,
                "inner %d srf / %d med\nouter %d srf / %d med\n",
                a.inner.surface,
                a.inner.medium,
                a.outer.surface,
                a.outer.medium);
}

/* mliColor */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliColor mliColor_set(const float r, const float g, const float b)
{
        struct mliColor rgb;
        rgb.r = r;
        rgb.g = g;
        rgb.b = b;
        return rgb;
}

struct mliColor mliColor_mix(
        const struct mliColor a,
        const struct mliColor b,
        const float refl)
{
        struct mliColor out;
        out.r = (1.f - refl) * a.r + refl * b.r;
        out.g = (1.f - refl) * a.g + refl * b.g;
        out.b = (1.f - refl) * a.b + refl * b.b;
        return out;
}

struct mliColor mliColor_mean(
        const struct mliColor colors[],
        const uint32_t num_colors)
{
        struct mliColor out = {0., 0., 0.};
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

struct mliColor mliColor_truncate_to_uint8(const struct mliColor color)
{
        struct mliColor out;
        out.r = color.r;
        out.g = color.g;
        out.b = color.b;
        if (out.r > 255.)
                out.r = 255.;
        if (out.r < 0.)
                out.r = 0.;
        if (out.g > 255.)
                out.g = 255.;
        if (out.g < 0.)
                out.g = 0.;
        if (out.b > 255.)
                out.b = 255.;
        if (out.b < 0.)
                out.b = 0.;
        return out;
}

int mliColor_equal(const struct mliColor a, const struct mliColor b)
{
        if (a.r != b.r)
                return 0;
        if (a.g != b.g)
                return 0;
        if (a.b != b.b)
                return 0;
        return 1;
}

int mliColor_is_valid_8bit_range(const struct mliColor c)
{
        if (MLI_IS_NAN(c.r))
                return 0;
        if (MLI_IS_NAN(c.g))
                return 0;
        if (MLI_IS_NAN(c.b))
                return 0;

        if (c.r < 0.0 || c.r >= 256.0)
                return 0;
        if (c.g < 0.0 || c.g >= 256.0)
                return 0;
        if (c.b < 0.0 || c.b >= 256.0)
                return 0;
        return 1;
}

struct mliColor mliColor_add(const struct mliColor u, const struct mliColor v)
{
        return mliColor_set(u.r + v.r, u.g + v.g, u.b + v.b);
}

struct mliColor mliColor_multiply(const struct mliColor c, const double f)
{
        return mliColor_set(c.r * f, c.g * f, c.b * f);
}

struct mliColor mliColor_multiply_elementwise(
        const struct mliColor u,
        const struct mliColor v)
{
        return mliColor_set(u.r * v.r, u.g * v.g, u.b * v.b);
}

/* mliColor_json */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliColor_from_json_token(
        struct mliColor *c,
        const struct mliJson *json,
        const uint64_t token)
{
        struct mliVec v;
        chk_msg(mliVec_from_json_token(&v, json, token),
                "Can not parse json-float-triple to color.");
        c->r = v.x;
        c->g = v.y;
        c->b = v.z;
        return 1;
error:
        return 0;
}

/* mliCube */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliVec mliCube_upper(const struct mliCube a)
{
        return mliVec_add(
                a.lower,
                mliVec_init(a.edge_length, a.edge_length, a.edge_length));
}

struct mliAABB mliCube_to_aabb(const struct mliCube a)
{
        struct mliAABB out;
        out.lower = a.lower;
        out.upper = mliCube_upper(a);
        return out;
}

struct mliVec mliCube_center(const struct mliCube a)
{
        return mliVec_init(
                a.lower.x + a.edge_length * .5,
                a.lower.y + a.edge_length * .5,
                a.lower.z + a.edge_length * .5);
}

struct mliCube mliCube_outermost_cube(const struct mliAABB a)
{
        struct mliCube cube;
        struct mliVec center;
        struct mliVec half_diagonal;
        struct mliVec diff;
        double max_half_length;
        center = mliAABB_center(a);
        diff = mliVec_substract(a.upper, a.lower);
        max_half_length = .5 * MLI_MAX3(diff.x, diff.y, diff.z);
        half_diagonal.x = max_half_length;
        half_diagonal.y = max_half_length;
        half_diagonal.z = max_half_length;
        cube.lower = mliVec_substract(center, half_diagonal);
        cube.edge_length = max_half_length * 2.;
        return cube;
}

struct mliCube mliCube_octree_child(
        const struct mliCube cube,
        const uint32_t sx,
        const uint32_t sy,
        const uint32_t sz)
{
        struct mliCube child;
        struct mliVec length;
        struct mliVec center = mliCube_center(cube);
        length = mliVec_substract(center, cube.lower);
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

struct mliCube mliCube_octree_child_code(
        const struct mliCube cube,
        const uint8_t a)
{
        struct mliCube child;
        struct mliVec length;
        struct mliVec center = mliCube_center(cube);
        length = mliVec_substract(center, cube.lower);
        child.lower = cube.lower;
        child.edge_length = .5 * cube.edge_length;
        if (MLI_IS_BIT(a, 2)) {
                child.lower.x += length.x;
        }
        if (MLI_IS_BIT(a, 1)) {
                child.lower.y += length.y;
        }
        if (MLI_IS_BIT(a, 0)) {
                child.lower.z += length.z;
        }
        return child;
}

int mliCube_equal(const struct mliCube a, const struct mliCube b)
{
        if (a.edge_length != b.edge_length)
                return 0;
        if (!mliVec_equal(a.lower, b.lower))
                return 0;
        return 1;
}

/* mliDynArray */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/* mliDynColor */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLIDYNARRAY_IMPLEMENTATION(mli, Color, struct mliColor)
MLIDYNARRAY_IMPLEMENTATION(mli, ColorPtr, struct mliColor *)

/* mliDynDouble */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLIDYNARRAY_IMPLEMENTATION(mli, Double, double)

/* mliDynFace */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLIDYNARRAY_IMPLEMENTATION(mli, Face, struct mliFace)

/* mliDynFloat */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLIDYNARRAY_IMPLEMENTATION(mli, Float, float)

/* mliDynMap */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

MLIDYNARRAY_IMPLEMENTATION(mli, Map, struct mliDynMapItem)

int mliDynMap_find(const struct mliDynMap *map, const char *key, uint64_t *idx)
{
        uint64_t i;
        for (i = 0; i < map->size; i++) {
                if (strcmp(map->array[i].key, key) == 0) {
                        (*idx) = i;
                        return 1;
                }
        }
        return 0;
}

int mliDynMap_has(const struct mliDynMap *map, const char *key)
{
        uint64_t idx;
        return mliDynMap_find(map, key, &idx);
}

int mliDynMap_insert(struct mliDynMap *map, const char *key, uint64_t value)
{
        struct mliDynMapItem item;
        if (map->size == 0u) {
                chk_msg(mliDynMap_malloc(map, 0u),
                        "Failed to initially malloc dyn-map.");
        }
        chk_msg(strlen(key) < MLI_NAME_CAPACITY, "Key too long.");
        chk_msg(0 == mliDynMap_has(map, key), "Key already in use.");
        memset(item.key, '\0', MLI_NAME_CAPACITY);
        strcpy(item.key, key);
        item.value = value;
        chk_msg(mliDynMap_push_back(map, item), "Failed to mmaloc item.");
        return 1;
error:
        return 0;
}

int mliDynMap_get(const struct mliDynMap *map, const char *key, uint64_t *value)
{
        uint64_t idx;
        chk_msg(mliDynMap_find(map, key, &idx), "Key does not exist.");
        (*value) = map->array[idx].value;
        return 1;
error:
        return 0;
}

/* mliDynMap_json */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliDynMap_insert_key_from_json(
        struct mliDynMap *map,
        const struct mliJson *json,
        const uint64_t token,
        const uint64_t value)
{
        char buff[MLI_NAME_CAPACITY] = {'\0'};
        const uint64_t name_strlen =
                (json->tokens[token].end - json->tokens[token].start);
        chk_msg(name_strlen < sizeof(buff), "Key is too long");
        chk_msg(json->tokens[token].type == JSMN_STRING,
                "Expected key to be of type string.");
        chk_msg(mliJson_cstr_by_token(json, token, buff, name_strlen + 1),
                "Failed to extract string from json.");
        chk_msg(mliDynMap_insert(map, buff, value),
                "Failed to insert name and value into map.");
        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, token);
        return 0;
}

int mliDynMap_get_value_for_string_from_json(
        const struct mliDynMap *map,
        const struct mliJson *json,
        const uint64_t token,
        uint32_t *out_value)
{
        char buff[MLI_NAME_CAPACITY] = {'\0'};
        uint64_t value;
        uint64_t name_strlen =
                (json->tokens[token].end - json->tokens[token].start);
        chk_msg(name_strlen < sizeof(buff), "Key is too long");
        chk_msg(json->tokens[token].type == JSMN_STRING,
                "Expected token to be of type string to be given to mliMap.");
        chk_msg(mliJson_cstr_by_token(json, token, buff, name_strlen + 1),
                "Failed to extract string from json.");
        chk_msg(mliDynMap_get(map, buff, &value),
                "Failed to get value for json-string-key from map.");
        (*out_value) = (uint32_t)value;

        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, token);
        return 0;
}

/* mliDynPhoton */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

MLIDYNARRAY_IMPLEMENTATION(mli, Photon, struct mliPhoton)

/* mliDynPhotonInteraction */
/* ----------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

MLIDYNARRAY_IMPLEMENTATION(mli, PhotonInteraction, struct mliPhotonInteraction)

int mliDynPhotonInteraction_time_of_flight(
        const struct mliDynPhotonInteraction *history,
        const struct mliScenery *scenery,
        const double wavelength,
        double *total_time_of_flight)
{
        uint64_t i;
        (*total_time_of_flight) = 0.0;
        for (i = 0; i < history->size; i++) {
                double time_of_flight = 0.0;
                chk_msg(mli_time_of_flight(
                                &scenery->materials,
                                &history->array[i],
                                wavelength,
                                &time_of_flight),
                        "Can't estimate time-of-flight.");
                (*total_time_of_flight) += time_of_flight;
        }
        return 1;
error:
        return 0;
}

void mliDynPhotonInteraction_print(
        const struct mliDynPhotonInteraction *history,
        const struct mliGeometry *scenery)
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
                struct mliPhotonInteraction phisec = history->array[i];
                printf(" % 3d  ", (int32_t)i);

                if (phisec.on_geometry_surface == 1) {
                        printf("(% 5d;% 5d,% 5d,% 5d)  ",
                               scenery->robject_ids[phisec.geometry_id.robj],
                               phisec.geometry_id.robj,
                               scenery->robjects[phisec.geometry_id.robj],
                               phisec.geometry_id.face);
                } else {
                        printf("            n/a            ");
                }

                printf("[% -.1e,% -.1e,% -.1e]  ",
                       phisec.position.x,
                       phisec.position.y,
                       phisec.position.z);

                mli_photoninteraction_type_to_string(phisec.type, type_string);

                printf("%-12s ", type_string);

                printf("{% 2ld,% 2ld}  ",
                       phisec.medium_coming_from,
                       phisec.medium_going_to);

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
                        printf("n/a  ");
                }
                printf("\n");
        }
}

/* mliDynUint32 */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLIDYNARRAY_IMPLEMENTATION(mli, Uint32, uint32_t)

/* mliDynVec */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */
MLIDYNARRAY_IMPLEMENTATION(mli, Vec, struct mliVec)

/* mliFace */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliFace mliFace_set(const uint32_t a, const uint32_t b, const uint32_t c)
{
        struct mliFace face;
        face.a = a;
        face.b = b;
        face.c = c;
        return face;
}

int mliFace_equal(const struct mliFace a, const struct mliFace b)
{
        if (a.a != b.a)
                return 0;
        if (a.b != b.b)
                return 0;
        if (a.c != b.c)
                return 0;
        return 1;
}

/* mliFrame */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

MLIDYNARRAY_IMPLEMENTATION(mli, FramePtr, struct mliFrame *)

struct mliFrame;

struct mliFrame mliFrame_init(void)
{
        struct mliFrame f;
        f.type = MLI_FRAME;
        f.id = 0u;
        f.frame2mother.translation = mliVec_init(0., 0., 0.);
        f.frame2mother.rotation = mliQuaternion_set_tait_bryan(0., 0., 0.);
        f.frame2root = f.frame2mother;
        f.mother = NULL;
        f.children = mliDynFramePtr_init();
        f.object = 0u;

        f.boundary_layers = mliDynUint32_init();
        return f;
}

void mliFrame_free(struct mliFrame *f)
{
        uint64_t c;
        if (f->type == MLI_FRAME) {
                for (c = 0; c < f->children.size; c++) {
                        struct mliFrame *child = f->children.array[c];
                        mliFrame_free(child);
                }
                mliDynFramePtr_free(&f->children);
        }
        if (f->type == MLI_OBJECT) {
                mliDynUint32_free(&f->boundary_layers);
        }
        (*f) = mliFrame_init();
}

int mliFrame_malloc(struct mliFrame *f, const uint64_t type)
{
        mliFrame_free(f);
        f->type = type;
        if (type == MLI_FRAME) {
                chk_msg(mliDynFramePtr_malloc(&f->children, 0u),
                        "Can not allocate children of frame.");
        }
        if (type == MLI_OBJECT) {
                chk_msg(mliDynUint32_malloc(&f->boundary_layers, 0u),
                        "Failed to malloc frame's boundary_layers.");
        }
        return 1;
error:
        return 0;
}

int mliFrame_set_mother_and_child(
        struct mliFrame *mother,
        struct mliFrame *child)
{
        chk_msg(mother->type == MLI_FRAME,
                "Expected mother to be of type FRAME");
        chk_msg(mliDynFramePtr_push_back(&mother->children, child),
                "Can not push back child-frame.");

        child->mother = (struct mliFrame *)mother;
        return 1;
error:
        return 0;
}

struct mliFrame *mliFrame_add(struct mliFrame *mother, const uint64_t type)
{
        struct mliFrame *child = NULL;
        chk_malloc(child, struct mliFrame, 1u);
        chk_msg(mliFrame_malloc(child, type), "Can not allocate child-frame.");
        chk_msg(mliFrame_set_mother_and_child(mother, child),
                "Can not allocate child-pointer.");
        return child;
error:
        return NULL;
}

int mli_type_to_string(const uint64_t type, char *s)
{
        switch (type) {
        case MLI_FRAME:
                sprintf(s, "frame");
                break;
        case MLI_OBJECT:
                sprintf(s, "object");
                break;
        default:
                chk_bad("Type is unknown.");
                break;
        }
        return 1;
error:
        return 0;
}

int mli_string_to_type(const char *s, uint64_t *type)
{
        if (strcmp(s, "frame") == 0) {
                *type = MLI_FRAME;
        } else if (strcmp(s, "object") == 0) {
                *type = MLI_OBJECT;
        } else {
                chk_bad("Type is unknown.");
        }
        return 1;
error:
        return 0;
}

void mliFrame_print_walk(const struct mliFrame *f, const uint64_t indention)
{
        uint64_t c;
        char type_string[1024];
        mli_type_to_string(f->type, type_string);
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
        if (f->type == MLI_OBJECT) {
                uint32_t ii;
                printf("%*s", (int)indention, "");
                printf("|-boundary_layers [");
                for (ii = 0; ii < f->boundary_layers.size; ii++) {
                        printf("%u,", f->boundary_layers.array[ii]);
                }
                printf("]\n");

                printf("%*s", (int)indention, "");
                printf("|-obj %u\n", f->object);
        }
        for (c = 0; c < f->children.size; c++) {
                const struct mliFrame *child = f->children.array[c];
                mliFrame_print_walk(child, indention + 4);
        }
}

void mliFrame_print(struct mliFrame *f) { mliFrame_print_walk(f, 0u); }

void mliFrame_set_frame2root(struct mliFrame *f)
{
        if (f->mother == NULL) {
                f->frame2root = f->frame2mother;
        } else {
                f->frame2root = mliHomTraComp_sequence(
                        f->frame2mother, f->mother->frame2root);
        }
        if (f->type == MLI_FRAME) {
                uint64_t c;
                for (c = 0; c < f->children.size; c++) {
                        struct mliFrame *child = f->children.array[c];
                        mliFrame_set_frame2root(child);
                }
        }
}

/* mliFrame_json */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliFrame_type_from_json_token(
        uint64_t *type,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t _t;
        int has_obj = mliJson_token_by_key(json, token, "obj", &_t);
        int has_children = mliJson_token_by_key(json, token, "children", &_t);

        if (has_obj && has_children) {
                chk_bad("Frame must not have both keys 'obj', and 'children'.");
        } else if (!has_obj && !has_children) {
                chk_bad("Frame must have either of keys 'obj', or 'children'.");
        } else if (has_obj && !has_children) {
                (*type) = MLI_OBJECT;
        } else if (!has_obj && has_children) {
                (*type) = MLI_FRAME;
        } else {
                chk_bad("Not expected to happen");
        }

        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, token);
        return 0;
}

int mliFrame_id_from_json_token(
        uint32_t *id,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t token_id;
        int64_t _id;
        chk_msg(mliJson_token_by_key(json, token, "id", &token_id),
                "Expected Frame to have key 'id'.");
        chk_msg(mliJson_int64_by_token(json, token_id + 1, &_id),
                "Failed to parse Frame's id.");
        chk_msg(_id >= 0, "Expected Frame's id >= 0.");
        (*id) = _id;

        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, token);
        return 0;
}

int mliFrame_pos_rot_from_json_token(
        struct mliHomTraComp *frame2mother,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t token_pos, token_rot;
        /* pos */
        chk_msg(mliJson_token_by_key(json, token, "pos", &token_pos),
                "Expected Frame to have key 'pos'.");
        chk_msg(mliVec_from_json_token(
                        &frame2mother->translation, json, token_pos + 1),
                "Failed to parse Frame's 'pos' from json.");

        /* rot */
        chk_msg(mliJson_token_by_key(json, token, "rot", &token_rot),
                "Expected Frame to have key 'rot'.");
        chk_msg(mliQuaternion_from_json(
                        &frame2mother->rotation, json, token_rot + 1),
                "Failed to parse Frame's 'rot' from json.");
        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, token);
        return 0;
}

int mliFrame_boundary_layers_form_json_token(
        struct mliDynUint32 *boundary_layers,
        const uint32_t object_idx,
        const struct mliObject *objects,
        const struct mliDynMap *boundary_layer_names,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t token_mtl_key, token_mtl;
        uint64_t material_idx;
        chk_msg(mliJson_token_by_key(json, token, "mtl", &token_mtl_key),
                "Expected 'mtl' in Frame.");
        token_mtl = token_mtl_key + 1;
        chk_msg(json->tokens[token_mtl].type == JSMN_OBJECT,
                "Expected 'mtl' to be a json-object {}.");

        for (material_idx = 0u;
             material_idx < objects[object_idx].num_materials;
             material_idx++) {
                const char *material_key_in_object =
                        objects[object_idx].material_names[material_idx].cstr;

                uint64_t token_material_key = 0u;
                uint32_t boundary_layer_idx = 0u;
                chk_msg(mliJson_token_by_key(
                                json,
                                token_mtl,
                                material_key_in_object,
                                &token_material_key),
                        "Expected object's material-key to be in "
                        "object-reference's mtls in tree.json.");

                chk_msg(mliDynMap_get_value_for_string_from_json(
                                boundary_layer_names,
                                json,
                                token_material_key + 1,
                                &boundary_layer_idx),
                        "Expected boundary-layer to exist in materials.");

                chk_msg(mliDynUint32_push_back(
                                boundary_layers, boundary_layer_idx),
                        "Failed to push-back boundary_layer_idx into "
                        "frame's boundary_layers.");
        }

        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, token);
        return 0;
}

int mliFrame_object_reference_form_json_token(
        uint32_t *object_reference,
        const struct mliJson *json,
        const uint64_t token,
        const struct mliDynMap *object_names)
{
        uint64_t token_obj_key;
        chk_msg(mliJson_token_by_key(json, token, "obj", &token_obj_key),
                "Expected object to have key 'obj'.");
        chk_msg(mliDynMap_get_value_for_string_from_json(
                        object_names,
                        json,
                        token_obj_key + 1,
                        object_reference),
                "Failed to get object-reference 'obj' from map");
        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, token);
        return 0;
}

int mliFrame_from_json(
        struct mliFrame *mother,
        const struct mliJson *json,
        const uint64_t token_children,
        const struct mliDynMap *object_names,
        const struct mliObject *objects,
        const struct mliDynMap *boundary_layer_names)
{
        uint64_t num_children;
        uint64_t c;
        chk_msg(json->tokens[token_children].type == JSMN_ARRAY,
                "Expected Frame's children to be a json-array '[]'.");
        num_children = json->tokens[token_children].size;
        for (c = 0; c < num_children; c++) {
                uint64_t token_child =
                        mliJson_token_by_index(json, token_children, c);
                struct mliFrame *child = NULL;
                uint64_t type;
                uint64_t token_grandchildren;

                chk_msg(mliFrame_type_from_json_token(&type, json, token_child),
                        "Failed to read type of Frame.");

                child = mliFrame_add(mother, type);
                chk_msg(child, "Failed to add child to frame.");

                chk_msg(mliFrame_pos_rot_from_json_token(
                                &child->frame2mother, json, token_child),
                        "Failed to set pos, and rot of Frame from json.");

                chk_msg(mliFrame_id_from_json_token(
                                &child->id, json, token_child),
                        "Failed to set id of Frame from json.");

                switch (type) {
                case MLI_FRAME:
                        chk_msg(mliJson_token_by_key(
                                        json,
                                        token_child,
                                        "children",
                                        &token_grandchildren),
                                "Expected child of type Frame to have "
                                "key 'children'.");
                        chk_msg(mliFrame_from_json(
                                        child,
                                        json,
                                        token_grandchildren + 1,
                                        object_names,
                                        objects,
                                        boundary_layer_names),
                                "Failed to populate grandchildren "
                                "Frames from json.");
                        break;
                case MLI_OBJECT:
                        chk_msg(mliFrame_object_reference_form_json_token(
                                        &child->object,
                                        json,
                                        token_child,
                                        object_names),
                                "Failed to parse object-reference "
                                "from json.");
                        chk_msg(mliFrame_boundary_layers_form_json_token(
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
error:
        return 0;
}

/* mliFresnel */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliFresnel mliFresnel_init(
        const struct mliVec incident,
        const struct mliVec normal,
        const double n_from,
        const double n_to)
{
        struct mliFresnel fresnel;
        fresnel.incident = incident;
        fresnel.normal = normal;
        fresnel.n_from = n_from;
        fresnel.n_to = n_to;

        fresnel._cosI = -1.0 * mliVec_dot(normal, incident);
        fresnel._n_from_over_n_to = n_from / n_to;
        fresnel._sinT2 =
                (fresnel._n_from_over_n_to * fresnel._n_from_over_n_to) *
                (1.0 - (fresnel._cosI * fresnel._cosI));
        fresnel._cosT = sqrt(1.0 - fresnel._sinT2);
        return fresnel;
}

double mliFresnel_reflection_propability(const struct mliFresnel fresnel)
{
        if (fresnel._sinT2 > 1.0) {
                /* total internal reflection */
                return 1.0;
        } else {
                const struct mliFresnel f = fresnel;
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

struct mliVec mliFresnel_reflection_direction(const struct mliFresnel fresnel)
{
        return mliVec_add(
                fresnel.incident,
                mliVec_multiply(fresnel.normal, fresnel._cosI * 2.0));
}

struct mliVec mliFresnel_refraction_direction(const struct mliFresnel fresnel)
{
        return mliVec_add(
                mliVec_multiply(fresnel.incident, fresnel._n_from_over_n_to),
                mliVec_multiply(
                        fresnel.normal,
                        fresnel._n_from_over_n_to * fresnel._cosI -
                                fresnel._cosT));
}

/* mliFunc */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliFunc mliFunc_init(void)
{
        struct mliFunc f;
        f.num_points = 0u;
        f.x = NULL;
        f.y = NULL;
        return f;
}

void mliFunc_free(struct mliFunc *f)
{
        free(f->x);
        free(f->y);
        (*f) = mliFunc_init();
}

int mliFunc_malloc(struct mliFunc *f, const uint32_t num_points)
{
        mliFunc_free(f);
        f->num_points = num_points;
        chk_malloc(f->x, double, f->num_points);
        chk_malloc(f->y, double, f->num_points);
        return 1;
error:
        mliFunc_free(f);
        return 0;
}

int mliFunc_x_is_strictly_increasing(const struct mliFunc *f)
{
        uint32_t i;
        for (i = 1; i < f->num_points; i++) {
                if (f->x[i] <= f->x[i - 1]) {
                        return 0;
                }
        }
        return 1;
}

int mliFunc_evaluate(const struct mliFunc *f, const double xarg, double *out)
{
        double y1, y0, x1, x0;
        uint32_t idx = mli_upper_compare_double(f->x, f->num_points, xarg);
        if (idx == 0) {
                chk_bad("mliFunc argument below lower bound.");
        } else if (idx == f->num_points) {
                chk_bad("mliFunc argument above upper bound.");
        } else {
                y1 = f->y[idx];
                y0 = f->y[idx - 1u];
                x1 = f->x[idx];
                x0 = f->x[idx - 1u];
                (*out) = mli_linear_interpolate_2d(xarg, x0, y0, x1, y1);
        }
        return 1;
error:
        return 0;
}

int mliFunc_fold_numeric(
        const struct mliFunc *a,
        const struct mliFunc *b,
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
                double ra = MLI_NAN;
                double rb = MLI_NAN;
                double x = xmin + (double)i * step_size;
                chk(mliFunc_evaluate(a, x, &ra));
                chk(mliFunc_evaluate(b, x, &rb));
                (*fold) += (ra * rb) * step_size;
        }
        return 1;
error:
        return 0;
}

int mliFunc_equal(const struct mliFunc a, const struct mliFunc b)
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

int mliFunc_is_valid(const struct mliFunc *func)
{
        uint64_t i;
        chk_msg(func->num_points >= 2,
                "Expected function to have at least two points. "
                "Evaluation is not possible when there is no valid range "
                "between two points.");

        for (i = 0; i < func->num_points; i++) {
                chk_msg(!MLI_IS_NAN(func->x[i]),
                        "Expected x-argument to be a real number, "
                        "but it is 'nan'.");
                chk_msg(!MLI_IS_NAN(func->y[i]),
                        "Expected y-value to be a real number, "
                        "but it is 'nan'.");
        }

        chk_msg(mliFunc_x_is_strictly_increasing(func),
                "Expected x-arguments to be strictly increasing, "
                "but they do not.");

        return 1;
error:
        return 0;
}

/* mliFunc_json */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliFunc_malloc_from_json_token(
        struct mliFunc *f,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t p;
        uint64_t point_token;
        chk_msg(json->tokens[token].type == JSMN_ARRAY,
                "Expected function-token to be a json-array.");
        chk_msg(mliFunc_malloc(f, json->tokens[token].size),
                "Can not allocate function.");
        point_token = token + 1;
        for (p = 0; p < f->num_points; p++) {
                chk_msg(json->tokens[point_token].type == JSMN_ARRAY,
                        "Expected function-x-y-pair to be a json-array.");
                chk_msg(json->tokens[point_token].size == 2,
                        "Expected function-x-y-pair to contain exactly 2 "
                        "tokens (x and y).");
                chk_msg(mliJson_double_by_token(
                                json, point_token + 1, &f->x[p]),
                        "Can not parse function-x-value.");
                chk_msg(mliJson_double_by_token(
                                json, point_token + 2, &f->y[p]),
                        "Can not parse function-y-value.");
                point_token += 3;
        }
        return 1;
error:
        mliFunc_free(f);
        return 0;
}

/* mliFunc_serialize */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliFunc_fwrite(const struct mliFunc *func, FILE *f)
{
        struct mliMagicId magic;
        chk(mliMagicId_set(&magic, "mliFunc"));
        chk_fwrite(&magic, sizeof(struct mliMagicId), 1u, f);

        chk_fwrite(&func->num_points, sizeof(uint32_t), 1u, f);
        chk_fwrite(func->x, sizeof(double), func->num_points, f);
        chk_fwrite(func->y, sizeof(double), func->num_points, f);
        return 1;
error:
        return 0;
}

int mliFunc_malloc_fread(struct mliFunc *func, FILE *f)
{
        uint32_t num_points;
        struct mliMagicId magic;
        chk_fread(&magic, sizeof(struct mliMagicId), 1u, f);
        chk(mliMagicId_has_word(&magic, "mliFunc"));
        mliMagicId_warn_version(&magic);

        chk_fread(&num_points, sizeof(uint32_t), 1u, f);
        chk(mliFunc_malloc(func, num_points));
        chk_fread(func->x, sizeof(double), func->num_points, f);
        chk_fread(func->y, sizeof(double), func->num_points, f);
        chk_msg(mliFunc_is_valid(func), "Expected function to be valid.");
        return 1;
error:
        mliFunc_free(func);
        return 0;
}

/* mliGeometry */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mliGeometry_init_objects(struct mliGeometry *geometry)
{
        geometry->num_objects = 0u;
        geometry->objects = NULL;
        geometry->object_names = NULL;
}

void mliGeometry_init_references(struct mliGeometry *geometry)
{
        geometry->num_robjects = 0u;
        geometry->robjects = NULL;
        geometry->robject_ids = NULL;
        geometry->robject2root = NULL;
}

struct mliGeometry mliGeometry_init(void)
{
        struct mliGeometry geometry;
        mliGeometry_init_objects(&geometry);
        mliGeometry_init_references(&geometry);
        return geometry;
}

void mliGeometry_free_objects(struct mliGeometry *geometry)
{
        uint32_t i;
        for (i = 0; i < geometry->num_objects; i++) {
                mliObject_free(&(geometry->objects[i]));
        }
        free(geometry->objects);
        free(geometry->object_names);
        mliGeometry_init_objects(geometry);
}

void mliGeometry_free_references(struct mliGeometry *geometry)
{
        free(geometry->robjects);
        free(geometry->robject_ids);
        free(geometry->robject2root);
        mliGeometry_init_references(geometry);
}

void mliGeometry_free(struct mliGeometry *geometry)
{
        mliGeometry_free_objects(geometry);
        mliGeometry_free_references(geometry);
        (*geometry) = mliGeometry_init();
}

int mliGeometry_malloc_objects(
        struct mliGeometry *geometry,
        const uint32_t num_objects)
{
        uint32_t i;
        geometry->num_objects = num_objects;
        chk_malloc(geometry->objects, struct mliObject, geometry->num_objects);
        chk_malloc(
                geometry->object_names, struct mliName, geometry->num_objects);
        for (i = 0; i < geometry->num_objects; i++) {
                geometry->objects[i] = mliObject_init();
                geometry->object_names[i] = mliName_init();
        }
        return 1;
error:
        mliGeometry_free_objects(geometry);
        return 0;
}

int mliGeometry_malloc_references(
        struct mliGeometry *geometry,
        const uint32_t num_robjects)
{
        geometry->num_robjects = num_robjects;
        chk_malloc(geometry->robjects, uint32_t, geometry->num_robjects);
        chk_malloc(geometry->robject_ids, uint32_t, geometry->num_robjects);
        chk_malloc(
                geometry->robject2root,
                struct mliHomTraComp,
                geometry->num_robjects);
        return 1;
error:
        mliGeometry_free_references(geometry);
        return 0;
}

int mliGeometry_malloc(
        struct mliGeometry *geometry,
        const uint32_t num_objects,
        const uint32_t num_robjects)
{
        mliGeometry_free(geometry);
        chk(mliGeometry_malloc_objects(geometry, num_objects));
        chk(mliGeometry_malloc_references(geometry, num_robjects));
        return 1;
error:
        mliGeometry_free(geometry);
        return 0;
}

void mliGeometry_info_fprint(FILE *f, const struct mliGeometry *geometry)
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

        for (i = 0; i < geometry->num_objects; i++) {
                fprintf(f, "    ");
                fprintf(f, "%5d ", i);
                fprintf(f, "%24s ", geometry->object_names[i].cstr);
                fprintf(f, "%8d ", geometry->objects[i].num_vertices);
                fprintf(f, "%8d ", geometry->objects[i].num_vertex_normals);
                fprintf(f, "%8d ", geometry->objects[i].num_faces);
                fprintf(f, "%8d ", geometry->objects[i].num_materials);
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

        for (rob = 0; rob < geometry->num_robjects; rob++) {
                fprintf(f, "    ");
                fprintf(f, "%6d ", rob);
                fprintf(f, "%6d ", geometry->robjects[rob]);
                fprintf(f, "%6d ", geometry->robject_ids[rob]);
                fprintf(f,
                        "% 6.1f,",
                        geometry->robject2root[rob].translation.x);
                fprintf(f,
                        "% 6.1f,",
                        geometry->robject2root[rob].translation.y);
                fprintf(f,
                        "% 6.1f ",
                        geometry->robject2root[rob].translation.z);
                fprintf(f, " ");
                fprintf(f, "% 6.1f,", geometry->robject2root[rob].rotation.x);
                fprintf(f, "% 6.1f,", geometry->robject2root[rob].rotation.y);
                fprintf(f, "% 6.1f;", geometry->robject2root[rob].rotation.z);
                fprintf(f, "% 6.1f ", geometry->robject2root[rob].rotation.w);
                fprintf(f, "\n");
        }
}

int mliGeometry_warn_objects(const struct mliGeometry *geometry)
{
        uint32_t o;
        for (o = 0; o < geometry->num_objects; o++) {
                uint32_t v, vn, mtl;
                chk_msg(mliObject_num_unused(
                                &(geometry->objects)[o], &v, &vn, &mtl),
                        "Can't estimate num unused v/vn/mtl in object.");
                if (v > 0 || vn > 0 || mtl > 0) {
                        fprintf(stderr,
                                "[WARNING] Object '%s' at [%u] ",
                                geometry->object_names[o].cstr,
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
error:
        return 0;
}

/* mliGeometryAndAccelerator */
/* ------------------------- */



/* mliGeometryId */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliGeometryId mliGeometryId_init(void)
{
        struct mliGeometryId id;
        id.robj = 0u;
        id.face = 0u;
        return id;
}

/* mliGeometryToMaterialMap */
/* ------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliGeometryToMaterialMap mliGeometryToMaterialMap_init(void)
{
        struct mliGeometryToMaterialMap map;
        map.num_robjects = 0u;
        map.total_num_boundary_layers = 0u;
        map.boundary_layers = NULL;
        map.first_boundary_layer_in_robject = NULL;
        return map;
}

void mliGeometryToMaterialMap_free(struct mliGeometryToMaterialMap *map)
{
        free(map->first_boundary_layer_in_robject);
        free(map->boundary_layers);
        (*map) = mliGeometryToMaterialMap_init();
}

int mliGeometryToMaterialMap_malloc(
        struct mliGeometryToMaterialMap *map,
        const uint32_t num_robjects,
        const uint32_t total_num_boundary_layers)
{
        mliGeometryToMaterialMap_free(map);
        map->num_robjects = num_robjects;
        map->total_num_boundary_layers = total_num_boundary_layers;

        chk_malloc(
                map->boundary_layers, uint32_t, map->total_num_boundary_layers);
        chk_malloc(
                map->first_boundary_layer_in_robject,
                uint32_t,
                map->num_robjects);
        MLI_ARRAY_SET(
                map->first_boundary_layer_in_robject, 0, map->num_robjects);
        MLI_ARRAY_SET(map->boundary_layers, 0, map->total_num_boundary_layers);
        return 1;
error:
        return 0;
}

uint32_t mliGeometryToMaterialMap_resolve_idx(
        const struct mliGeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx)
{
        uint32_t idx = map->first_boundary_layer_in_robject[robject_idx];
        idx += material_idx;
        return idx;
}

uint32_t mliGeometryToMaterialMap_get(
        const struct mliGeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx)
{
        return map->boundary_layers[mliGeometryToMaterialMap_resolve_idx(
                map, robject_idx, material_idx)];
}

void mliGeometryToMaterialMap_set(
        const struct mliGeometryToMaterialMap *map,
        const uint32_t robject_idx,
        const uint32_t material_idx,
        const uint32_t boundary_layer_idx)
{
        map->boundary_layers[mliGeometryToMaterialMap_resolve_idx(
                map, robject_idx, material_idx)] = boundary_layer_idx;
}

uint32_t mliGeometryToMaterialMap_num_boundary_layers_in_robject(
        const struct mliGeometryToMaterialMap *map,
        const uint32_t robject_idx)
{
        const uint32_t start =
                mliGeometryToMaterialMap_resolve_idx(map, robject_idx, 0u);
        uint32_t end = start;
        if (robject_idx + 1 < map->num_robjects) {
                end = mliGeometryToMaterialMap_resolve_idx(
                        map, robject_idx + 1, 0u);
        } else {
                end = map->total_num_boundary_layers;
        }
        return end - start;
}

void mliGeometryToMaterialMap_info_fprint(
        FILE *f,
        const struct mliGeometryToMaterialMap *map)
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
                bdl_start = mliGeometryToMaterialMap_resolve_idx(map, robj, 0u);
                num_bdls =
                        mliGeometryToMaterialMap_num_boundary_layers_in_robject(
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

/* mliGeometryToMaterialMap_equal */
/* ------------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliGeometryToMaterialMap_equal(
        const struct mliGeometryToMaterialMap *a,
        const struct mliGeometryToMaterialMap *b)
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
error:
        return 0;
}

/* mliGeometryToMaterialMap_serialize */
/* ---------------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliGeometryToMaterialMap_malloc_fread(
        struct mliGeometryToMaterialMap *geomap,
        FILE *f)
{
        uint32_t num_robjects = 0u;
        uint32_t total_num_boundary_layers = 0u;
        struct mliMagicId magic;

        /* magic identifier */
        chk_fread(&magic, sizeof(struct mliMagicId), 1u, f);
        chk(mliMagicId_has_word(&magic, "mliGeometryToMaterialMap"));
        mliMagicId_warn_version(&magic);

        /* payload */
        chk_fread(&num_robjects, sizeof(uint32_t), 1u, f);
        chk_fread(&total_num_boundary_layers, sizeof(uint32_t), 1u, f);
        chk_msg(mliGeometryToMaterialMap_malloc(
                        geomap, num_robjects, total_num_boundary_layers),
                "Failed to malloc mliGeometryToMaterialMap.");
        chk_fread(
                geomap->boundary_layers,
                sizeof(uint32_t),
                geomap->total_num_boundary_layers,
                f);
        chk_fread(
                geomap->first_boundary_layer_in_robject,
                sizeof(uint32_t),
                geomap->num_robjects,
                f);
        return 1;
error:
        return 0;
}

int mliGeometryToMaterialMap_fwrite(
        const struct mliGeometryToMaterialMap *geomap,
        FILE *f)
{
        struct mliMagicId magic;

        /* magic identifier */
        chk(mliMagicId_set(&magic, "mliGeometryToMaterialMap"));
        chk_fwrite(&magic, sizeof(struct mliMagicId), 1u, f);

        /* payload */
        chk_fwrite(&geomap->num_robjects, sizeof(uint32_t), 1, f);
        chk_fwrite(&geomap->total_num_boundary_layers, sizeof(uint32_t), 1, f);
        chk_fwrite(
                geomap->boundary_layers,
                sizeof(uint32_t),
                geomap->total_num_boundary_layers,
                f);
        chk_fwrite(
                geomap->first_boundary_layer_in_robject,
                sizeof(uint32_t),
                geomap->num_robjects,
                f);
        return 1;
error:
        return 0;
}

/* mliGeometryToMaterialMap_valid */
/* ------------------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliGeometryToMaterialMap_valid(
        const struct mliGeometryToMaterialMap *geomap)
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
error:
        return 0;
}

int mliGeometryToMaterialMap_valid_wrt_Geometry(
        const struct mliGeometryToMaterialMap *geomap,
        const struct mliGeometry *geometry)
{
        uint32_t robj = 0u;
        uint32_t total_num_boundary_layers = 0u;
        chk_msg(geomap->num_robjects == geometry->num_robjects,
                "Expected num_robjects to be equal in Geometry and GeoMtlMap.");

        for (robj = 0u; robj < geomap->num_robjects; robj++) {
                const uint32_t obj = geometry->robjects[robj];
                const uint32_t obj_num_materials =
                        geometry->objects[obj].num_materials;
                chk_msg(mliGeometryToMaterialMap_num_boundary_layers_in_robject(
                                geomap, robj) == obj_num_materials,
                        "Expected robject to have same num boundary-layers.");
                total_num_boundary_layers += obj_num_materials;
        }
        chk_msg(total_num_boundary_layers == geomap->total_num_boundary_layers,
                "Expected total_num_boundary_layers to match the Geometry.");

        return 1;
error:
        return 0;
}

int mliGeometryToMaterialMap_valid_wrt_Materials(
        const struct mliGeometryToMaterialMap *geomap,
        const struct mliMaterials *materials)
{
        uint32_t i = 0u;
        for (i = 0u; i < geomap->total_num_boundary_layers; i++) {
                chk_msg(geomap->boundary_layers[i] <
                                materials->num_boundary_layers,
                        "Expected geomap's boundary_layers[i] to refer to "
                        "a valid boundary_layer in Materials.");
        }
        return 1;
error:
        return 0;
}

/* mliGeometry_AABB */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliGeometry_robject_has_overlap_aabb(
        const struct mliGeometryAndAccelerator *accgeo,
        const uint32_t robject_idx,
        const struct mliAABB aabb)
{
        const struct mliAABB robject_aabb =
                accgeo->accelerator->robject_aabbs[robject_idx];
        if (mliAABB_is_overlapping(aabb, robject_aabb)) {
                /* Only test the object's faces when its own aabb has an
                 * overlap with the aabb to test for.
                 */
                const uint32_t obj_idx =
                        accgeo->geometry->robjects[robject_idx];
                const struct mliObject *obj_ptr =
                        &accgeo->geometry->objects[obj_idx];
                const struct mliHomTra robj2root = mliHomTra_from_compact(
                        accgeo->geometry->robject2root[robject_idx]);
                return mliObject_has_overlap_aabb(obj_ptr, robj2root, aabb);
        } else {
                return 0;
        }
}

int mliGeometry_robject_has_overlap_aabb_void(
        const void *accgeo,
        const uint32_t robject_idx,
        const struct mliAABB aabb)
{
        return mliGeometry_robject_has_overlap_aabb(
                (const struct mliGeometryAndAccelerator *)accgeo,
                robject_idx,
                aabb);
}

/* mliGeometry_equal */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliGeometry_objects_equal(
        const struct mliGeometry *a,
        const struct mliGeometry *b)
{
        uint32_t i = 0u;
        chk_msg(a->num_objects == b->num_objects,
                "Expected num_objects to be equal.");

        for (i = 0; i < a->num_objects; i++) {
                chk_msg(mliObject_equal(&a->objects[i], &b->objects[i]),
                        "Expected object to be equal.");
                chk_msg(mliName_equal(&a->object_names[i], &b->object_names[i]),
                        "Expected object_name to be equal.");
        }
        return 1;
error:
        fprintf(stderr, "In geometry.object[%u]\n", i);
        return 0;
}

int mliGeometry_object_references_equal(
        const struct mliGeometry *a,
        const struct mliGeometry *b)
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

                chk_msg(mliHomTraComp_equal(
                                a->robject2root[rob], b->robject2root[rob]),
                        "Expected homogenous transformation of "
                        "object-references to be equal");
        }
        return 1;
error:
        fprintf(stderr, "In geometry.object_reference[%lu]\n", rob);
        return 0;
}

int mliGeometry_equal(const struct mliGeometry *a, const struct mliGeometry *b)
{
        chk_msg(mliGeometry_objects_equal(a, b),
                "Expected objects to be equal.");
        chk_msg(mliGeometry_object_references_equal(a, b),
                "Expected object-references to be equal.");
        return 1;
error:
        return 0;
}

/* mliGeometry_serialize */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliGeometry_malloc_fread(struct mliGeometry *geometry, FILE *f)
{
        uint32_t i;
        uint32_t num_objects = 0u;
        uint32_t num_robjects = 0u;
        struct mliMagicId magic;

        /* magic identifier */
        chk_fread(&magic, sizeof(struct mliMagicId), 1u, f);
        chk(mliMagicId_has_word(&magic, "mliGeometry"));
        mliMagicId_warn_version(&magic);

        /* payload */
        chk_fread(&num_objects, sizeof(uint32_t), 1u, f);
        chk_fread(&num_robjects, sizeof(uint32_t), 1u, f);

        chk_msg(mliGeometry_malloc(geometry, num_objects, num_robjects),
                "Failed to malloc robjects in mliGeometry.");

        for (i = 0; i < geometry->num_objects; i++) {
                chk_msg(mliObject_malloc_fread(&geometry->objects[i], f),
                        "Failed to read object into geometry.");
        }
        chk_fread(
                geometry->object_names,
                sizeof(struct mliName),
                geometry->num_objects,
                f);

        chk_fread(
                geometry->robjects,
                sizeof(uint32_t),
                geometry->num_robjects,
                f);
        chk_fread(
                geometry->robject_ids,
                sizeof(uint32_t),
                geometry->num_robjects,
                f);
        chk_fread(
                geometry->robject2root,
                sizeof(struct mliHomTraComp),
                geometry->num_robjects,
                f);

        return 1;
error:
        return 0;
}

int mliGeometry_fwrite(const struct mliGeometry *geometry, FILE *f)
{
        uint32_t i;
        struct mliMagicId magic;

        /* magic identifier */
        chk(mliMagicId_set(&magic, "mliGeometry"));
        chk_fwrite(&magic, sizeof(struct mliMagicId), 1u, f);

        /* payload */
        chk_fwrite(&geometry->num_objects, sizeof(uint32_t), 1, f);
        chk_fwrite(&geometry->num_robjects, sizeof(uint32_t), 1, f);

        for (i = 0; i < geometry->num_objects; i++) {
                chk_msg(mliObject_fwrite(&geometry->objects[i], f),
                        "Failed to write objects.");
        }
        chk_fwrite(
                geometry->object_names,
                sizeof(struct mliName),
                geometry->num_objects,
                f);

        chk_fwrite(
                geometry->robjects,
                sizeof(uint32_t),
                geometry->num_robjects,
                f);
        chk_fwrite(
                geometry->robject_ids,
                sizeof(uint32_t),
                geometry->num_robjects,
                f);
        chk_fwrite(
                geometry->robject2root,
                sizeof(struct mliHomTraComp),
                geometry->num_robjects,
                f);

        return 1;
error:
        return 0;
}

/* mliGeometry_valid */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliGeometry_valid_objects(const struct mliGeometry *geometry)
{
        uint32_t i;
        for (i = 0; i < geometry->num_objects; i++) {
                chk_msg(mliObject_is_valid(&geometry->objects[i]),
                        "Expected object to be valid.");
                chk_msg(mliName_valid(&geometry->object_names[i]),
                        "Expected object_name to be valid.");
        }
        return 1;
error:
        fprintf(stderr, "In geometry.objects[%u]\n", i);
        return 0;
}

int mliGeometry_valid_robjects_HomTras(const struct mliGeometry *geometry)
{
        uint32_t i;
        for (i = 0; i < geometry->num_robjects; i++) {
                const struct mliVec t = geometry->robject2root[i].translation;
                const struct mliQuaternion q =
                        geometry->robject2root[i].rotation;
                chk_msg(!MLI_IS_NAN(t.x), "translation.x is 'nan'.");
                chk_msg(!MLI_IS_NAN(t.y), "translation.y is 'nan'.");
                chk_msg(!MLI_IS_NAN(t.z), "translation.z is 'nan'.");

                chk_msg(!MLI_IS_NAN(q.w), "quaternion.w is 'nan'.");
                chk_msg(!MLI_IS_NAN(q.x), "quaternion.x is 'nan'.");
                chk_msg(!MLI_IS_NAN(q.y), "quaternion.y is 'nan'.");
                chk_msg(!MLI_IS_NAN(q.z), "quaternion.z is 'nan'.");
        }
        return 1;
error:
        fprintf(stderr, "In geometry.robject2root[%u]\n", i);
        return 0;
}

int mliGeometry_valid_object_references(const struct mliGeometry *geometry)
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
error:
        return 0;
}

int mliGeometry_valid(const struct mliGeometry *geometry)
{
        chk_msg(mliGeometry_valid_objects(geometry),
                "Expected objects to be valid.");
        chk_msg(mliGeometry_valid_robjects_HomTras(geometry),
                "Expected robject transformations to be free of 'nan'.");
        chk_msg(mliGeometry_valid_object_references(geometry),
                "Expected object-references to be valid.");
        return 1;
error:
        return 0;
}

/* mliHomTra */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliHomTraComp mliHomTraComp_set(
        const struct mliVec translation,
        const struct mliQuaternion rotation)
{
        struct mliHomTraComp comp;
        comp.translation = translation;
        comp.rotation = rotation;
        return comp;
}

struct mliHomTra mliHomTra_from_compact(const struct mliHomTraComp trafo)
{
        struct mliHomTra t;
        t.translation = trafo.translation;
        t.rotation = mliQuaternion_to_matrix(trafo.rotation);
        return t;
}

struct mliVec mli_transform_orientation(
        const struct mliMat *rotation,
        const struct mliVec ori)
{
        struct mliVec out;
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

struct mliVec mli_transform_orientation_inverse(
        const struct mliMat *rotation,
        const struct mliVec ori)
{
        struct mliVec out;
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

struct mliVec mli_transform_position(
        const struct mliMat *rotation,
        const struct mliVec translation,
        const struct mliVec pos)
{
        struct mliVec out;
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

struct mliVec mli_transform_position_inverse(
        const struct mliMat *rotation,
        const struct mliVec translation,
        const struct mliVec pos)
{
        struct mliVec out;
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

struct mliRay mli_transform_ray(
        const struct mliMat *rotation,
        const struct mliVec translation,
        const struct mliRay in)
{
        struct mliRay out;
        out.support = mli_transform_position(rotation, translation, in.support);
        out.direction = mli_transform_orientation(rotation, in.direction);
        return out;
}

struct mliRay mli_transform_ray_inverse(
        const struct mliMat *rotation,
        const struct mliVec translation,
        const struct mliRay in)
{
        struct mliRay out;
        out.support = mli_transform_position_inverse(
                rotation, translation, in.support);
        out.direction =
                mli_transform_orientation_inverse(rotation, in.direction);
        return out;
}

struct mliRay mliHomTra_ray(const struct mliHomTra *t, const struct mliRay in)
{
        return mli_transform_ray(&t->rotation, t->translation, in);
}

struct mliRay mliHomTra_ray_inverse(
        const struct mliHomTra *t,
        const struct mliRay in)
{
        return mli_transform_ray_inverse(&t->rotation, t->translation, in);
}

struct mliVec mliHomTra_pos(const struct mliHomTra *t, const struct mliVec in)
{
        return mli_transform_position(&t->rotation, t->translation, in);
}

struct mliVec mliHomTra_pos_inverse(
        const struct mliHomTra *t,
        const struct mliVec in)
{
        return mli_transform_position_inverse(&t->rotation, t->translation, in);
}

struct mliVec mliHomTra_dir(const struct mliHomTra *t, const struct mliVec in)
{
        return mli_transform_orientation(&t->rotation, in);
}

struct mliVec mliHomTra_dir_inverse(
        const struct mliHomTra *t,
        const struct mliVec in)
{
        return mli_transform_orientation_inverse(&t->rotation, in);
}

int mliHomTraComp_equal(
        const struct mliHomTraComp a,
        const struct mliHomTraComp b)
{
        if (!mliVec_equal(a.translation, b.translation))
                return 0;
        if (!mliQuaternion_equal(a.rotation, b.rotation))
                return 0;
        return 1;
}

struct mliHomTraComp mliHomTraComp_sequence(
        const struct mliHomTraComp a,
        const struct mliHomTraComp b)
{
        struct mliHomTra b_;
        struct mliHomTraComp s;
        b_ = mliHomTra_from_compact(b);
        s.translation = mliHomTra_pos(&b_, a.translation);
        s.rotation = mliQuaternion_product(b.rotation, a.rotation);
        return s;
}

/* mliImage */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliImage mliImage_init(void)
{
        struct mliImage img;
        img.num_cols = 0u;
        img.num_rows = 0u;
        img.raw = NULL;
        return img;
}

void mliImage_free(struct mliImage *img)
{
        free(img->raw);
        (*img) = mliImage_init();
}

int mliImage_malloc(
        struct mliImage *img,
        const uint32_t num_cols,
        const uint32_t num_rows)
{
        mliImage_free(img);
        img->num_cols = num_cols;
        img->num_rows = num_rows;
        chk_malloc(img->raw, struct mliColor, img->num_cols * img->num_rows);
        return 1;
error:
        mliImage_free(img);
        return 0;
}

uint32_t mliImage_idx(
        const struct mliImage *img,
        const uint32_t col,
        const uint32_t row)
{
        return col * img->num_rows + row;
}

void mliImage_set(
        const struct mliImage *img,
        const uint32_t col,
        const uint32_t row,
        const struct mliColor color)
{
        img->raw[mliImage_idx(img, col, row)] = color;
}

void mliImage_set_all_pixel(
        const struct mliImage *img,
        const struct mliColor color)
{
        uint64_t row, col;
        for (row = 0; row < img->num_rows; row++) {
                for (col = 0; col < img->num_cols; col++) {
                        img->raw[mliImage_idx(img, col, row)] = color;
                }
        }
}

struct mliColor mliImage_at(
        const struct mliImage *img,
        const uint32_t col,
        const uint32_t row)
{
        struct mliColor out;
        out = img->raw[mliImage_idx(img, col, row)];
        return out;
}

int mliImage_scale_down_twice(
        const struct mliImage *source,
        struct mliImage *destination)
{
        uint64_t row, col, sr, sc;
        chk_msg(destination->num_cols * 2u == source->num_cols,
                "Expected destination.num_cols*2u == source.num_cols");
        chk_msg(destination->num_rows * 2u == source->num_rows,
                "Expected destination.num_rows*2u == source.num_rows");
        for (row = 0; row < destination->num_rows; row++) {
                for (col = 0; col < destination->num_cols; col++) {
                        struct mliColor mix[4];
                        sr = row * 2u;
                        sc = col * 2u;
                        mix[0] = mliImage_at(source, sc + 0, sr + 0);
                        mix[1] = mliImage_at(source, sc + 0, sr + 1);
                        mix[2] = mliImage_at(source, sc + 1, sr + 0);
                        mix[3] = mliImage_at(source, sc + 1, sr + 1);
                        mliImage_set(
                                destination, col, row, mliColor_mean(mix, 4));
                }
        }
        return 1;
error:
        return 0;
}

void mliImage_sobel(const struct mliImage *image, struct mliImage *out)
{
        uint64_t idx_cm1_rp1;
        uint64_t idx_cm1_rp0;
        uint64_t idx_cm1_rm1;

        uint64_t idx_cp1_rp1;
        uint64_t idx_cp1_rp0;
        uint64_t idx_cp1_rm1;

        uint64_t idx_cp0_rp1;
        uint64_t idx_cp0_rm1;

        uint64_t col, row;
        for (col = 1; col < image->num_cols - 1; col++) {
                for (row = 1; row < image->num_rows - 1; row++) {
                        double xr = 0;
                        double xg = 0;
                        double xb = 0;
                        double yr = 0;
                        double yg = 0;
                        double yb = 0;
                        uint64_t idx;

                        idx_cm1_rp1 = mliImage_idx(image, col - 1, row + 1);
                        idx_cm1_rp0 = mliImage_idx(image, col - 1, row);
                        idx_cm1_rm1 = mliImage_idx(image, col - 1, row - 1);

                        idx_cp1_rp1 = mliImage_idx(image, col + 1, row + 1);
                        idx_cp1_rp0 = mliImage_idx(image, col + 1, row);
                        idx_cp1_rm1 = mliImage_idx(image, col + 1, row - 1);

                        idx_cp0_rp1 = mliImage_idx(image, col, row + 1);
                        idx_cp0_rm1 = mliImage_idx(image, col, row - 1);

                        idx = mliImage_idx(out, col, row);

                        xr += -1. * image->raw[idx_cm1_rp1].r;
                        xg += -1. * image->raw[idx_cm1_rp1].g;
                        xb += -1. * image->raw[idx_cm1_rp1].b;

                        xr += -2. * image->raw[idx_cm1_rp0].r;
                        xg += -2. * image->raw[idx_cm1_rp0].g;
                        xb += -2. * image->raw[idx_cm1_rp0].b;

                        xr += -1. * image->raw[idx_cm1_rm1].r;
                        xg += -1. * image->raw[idx_cm1_rm1].g;
                        xb += -1. * image->raw[idx_cm1_rm1].b;

                        xr += +1. * image->raw[idx_cp1_rp1].r;
                        xg += +1. * image->raw[idx_cp1_rp1].g;
                        xb += +1. * image->raw[idx_cp1_rp1].b;

                        xr += +2. * image->raw[idx_cp1_rp0].r;
                        xg += +2. * image->raw[idx_cp1_rp0].g;
                        xb += +2. * image->raw[idx_cp1_rp0].b;

                        xr += +1. * image->raw[idx_cp1_rm1].r;
                        xg += +1. * image->raw[idx_cp1_rm1].g;
                        xb += +1. * image->raw[idx_cp1_rm1].b;

                        yr += -1. * image->raw[idx_cm1_rp1].r;
                        yg += -1. * image->raw[idx_cm1_rp1].g;
                        yb += -1. * image->raw[idx_cm1_rp1].b;

                        yr += -2. * image->raw[idx_cp0_rp1].r;
                        yg += -2. * image->raw[idx_cp0_rp1].g;
                        yb += -2. * image->raw[idx_cp0_rp1].b;

                        yr += -1. * image->raw[idx_cp1_rp1].r;
                        yg += -1. * image->raw[idx_cp1_rp1].g;
                        yb += -1. * image->raw[idx_cp1_rp1].b;

                        yr += +1. * image->raw[idx_cp1_rm1].r;
                        yg += +1. * image->raw[idx_cp1_rm1].g;
                        yb += +1. * image->raw[idx_cp1_rm1].b;

                        yr += +2. * image->raw[idx_cp0_rm1].r;
                        yg += +2. * image->raw[idx_cp0_rm1].g;
                        yb += +2. * image->raw[idx_cp0_rm1].b;

                        yr += +1. * image->raw[idx_cm1_rm1].r;
                        yg += +1. * image->raw[idx_cm1_rm1].g;
                        yb += +1. * image->raw[idx_cm1_rm1].b;

                        out->raw[idx].r = mli_hypot(xr, yr);
                        out->raw[idx].g = mli_hypot(xg, yg);
                        out->raw[idx].b = mli_hypot(xb, yb);
                }
        }
}

void mliImage_luminance_threshold_dilatation(
        const struct mliImage *image,
        const float threshold,
        struct mliImage *out)
{
        const int32_t rows = image->num_rows;
        const int32_t cols = image->num_cols;
        int32_t col, row;
        const struct mliColor color_max = mliColor_set(255., 255., 255.);
        for (row = 0; row < rows; row++) {
                for (col = 0; col < cols; col++) {
                        const struct mliColor color_at =
                                mliImage_at(image, col, row);
                        const float luminance =
                                (color_at.r + color_at.g + color_at.b);
                        if (luminance > threshold) {
                                int32_t orow, ocol;
                                for (orow = -1; orow < 2; orow++) {
                                        for (ocol = -1; ocol < 2; ocol++) {
                                                if (row + orow >= 0 &&
                                                    col + ocol >= 0 &&
                                                    row + orow < rows &&
                                                    col + ocol < cols) {
                                                        mliImage_set(
                                                                out,
                                                                col + ocol,
                                                                row + orow,
                                                                color_max);
                                                }
                                        }
                                }
                        }
                }
        }
}

void mliImage_from_sum_and_exposure(
        const struct mliImage *sum,
        const struct mliImage *exposure,
        struct mliImage *out)
{
        uint64_t pix;
        for (pix = 0u; pix < out->num_rows * out->num_cols; pix++) {
                out->raw[pix].r = sum->raw[pix].r / exposure->raw[pix].r;
                out->raw[pix].g = sum->raw[pix].g / exposure->raw[pix].g;
                out->raw[pix].b = sum->raw[pix].b / exposure->raw[pix].b;
        }
}

void mliPixels_set_all_from_image(
        struct mliPixels *pixels,
        const struct mliImage *image)
{
        struct mliPixelWalk walk =
                mliPixelWalk_set(image->num_cols, image->num_rows, 16u);
        uint64_t i;
        const uint64_t num_pixel = image->num_rows * image->num_cols;

        assert(image->num_cols * image->num_rows == pixels->num_pixels);
        for (i = 0; i < num_pixel; i++) {
                pixels->pixels[i] = mliPixelWalk_get(&walk);
                mliPixelWalk_walk(&walk);
        }
}

void mliPixels_above_threshold(
        const struct mliImage *image,
        const float threshold,
        struct mliPixels *pixels)
{
        struct mliPixelWalk walk =
                mliPixelWalk_set(image->num_cols, image->num_rows, 16u);
        uint64_t i;
        const uint64_t num_pixel = image->num_rows * image->num_cols;

        pixels->num_pixels_to_do = 0;
        for (i = 0; i < num_pixel; i++) {
                struct mliPixel px = mliPixelWalk_get(&walk);
                double lum = 0.0;
                struct mliColor c = mliImage_at(image, px.col, px.row);
                lum = c.r + c.g + c.b;
                if (lum > threshold) {
                        pixels->pixels[pixels->num_pixels_to_do] = px;
                        pixels->num_pixels_to_do += 1;
                }
                mliPixelWalk_walk(&walk);
        }
}

void mliImage_assign_pixel_colors_to_sum_and_exposure_image(
        const struct mliPixels *pixels,
        const struct mliImage *colors,
        struct mliImage *sum_image,
        struct mliImage *exposure_image)
{
        uint64_t pix;
        for (pix = 0u; pix < pixels->num_pixels; pix++) {
                const uint64_t idx = mliImage_idx(
                        sum_image,
                        pixels->pixels[pix].col,
                        pixels->pixels[pix].row);
                sum_image->raw[idx].r += colors->raw[idx].r;
                sum_image->raw[idx].g += colors->raw[idx].g;
                sum_image->raw[idx].b += colors->raw[idx].b;
                exposure_image->raw[idx].r += 1.;
                exposure_image->raw[idx].g += 1.;
                exposure_image->raw[idx].b += 1.;
        }
}

void mliImage_copy(const struct mliImage *source, struct mliImage *destination)
{
        uint64_t pix;
        assert(source->num_rows == destination->num_rows);
        assert(source->num_cols == destination->num_cols);
        for (pix = 0u; pix < source->num_rows * source->num_cols; pix++) {
                destination->raw[pix] = source->raw[pix];
        }
}

void mliImage_fabs_difference(
        const struct mliImage *a,
        const struct mliImage *b,
        struct mliImage *out)
{
        uint64_t pix;
        assert(a->num_cols == b->num_cols);
        assert(a->num_rows == b->num_rows);
        assert(a->num_cols == out->num_cols);
        assert(a->num_rows == out->num_rows);
        for (pix = 0; pix < out->num_cols * out->num_rows; pix++) {
                out->raw[pix].r = fabs(a->raw[pix].r - b->raw[pix].r);
                out->raw[pix].g = fabs(a->raw[pix].g - b->raw[pix].g);
                out->raw[pix].b = fabs(a->raw[pix].b - b->raw[pix].b);
        }
}

void mliImage_histogram(
        struct mliImage *img,
        const double *col_bin_edges,
        const double *row_bin_edges,
        const double col_val,
        const double row_val,
        const struct mliColor weight)
{
        uint64_t col_upper_idx, row_upper_idx;
        int valid_col, valid_row;

        col_upper_idx = mli_upper_compare_double(
                col_bin_edges, img->num_cols + 1, col_val);

        row_upper_idx = mli_upper_compare_double(
                row_bin_edges, img->num_rows + 1, row_val);

        valid_col = col_upper_idx > 0 && col_upper_idx < img->num_cols + 1;
        valid_row = row_upper_idx > 0 && row_upper_idx < img->num_rows + 1;

        if (valid_col && valid_row) {
                const uint32_t col_idx = col_upper_idx - 1;
                const uint32_t row_idx = row_upper_idx - 1;
                const uint64_t pix = mliImage_idx(img, col_idx, row_idx);
                img->raw[pix] = mliColor_add(img->raw[pix], weight);
        }
}

struct mliColor mliImage_max(const struct mliImage *img)
{
        struct mliColor max = mliColor_set(-FLT_MAX, -FLT_MAX, -FLT_MAX);
        uint64_t pix;
        for (pix = 0u; pix < img->num_rows * img->num_cols; pix++) {
                struct mliColor c = img->raw[pix];
                if (c.r > max.r) {
                        max.r = c.r;
                }
                if (c.g > max.g) {
                        max.g = c.g;
                }
                if (c.b > max.b) {
                        max.b = c.b;
                }
        }
        return max;
}

void mliImage_multiply(struct mliImage *img, const struct mliColor color)
{
        uint64_t pix;
        for (pix = 0u; pix < img->num_rows * img->num_cols; pix++) {
                img->raw[pix].r *= color.r;
                img->raw[pix].g *= color.g;
                img->raw[pix].b *= color.b;
        }
}

void mliImage_divide_pixelwise(
        const struct mliImage *numerator,
        const struct mliImage *denominator,
        struct mliImage *out)
{
        uint64_t p;
        for (p = 0; p < out->num_rows * out->num_cols; p++) {
                out->raw[p].r = numerator->raw[p].r / denominator->raw[p].r;
                out->raw[p].g = numerator->raw[p].g / denominator->raw[p].g;
                out->raw[p].b = numerator->raw[p].b / denominator->raw[p].b;
        }
}

/* mliImage_ppm */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliImage_malloc_fread(struct mliImage *img, FILE *f)
{
        char line[1024];
        uint32_t num_commen_lines = 0;
        uint32_t num_cols;
        uint32_t num_rows;
        uint32_t col;
        uint32_t row;
        chk_msg(fgets(line, 1024, f), "Can't read header-line.") chk_msg(
                strcmp(line, "P6\n") == 0, "Expected starts with 'P6'.");
        while (1) {
                chk_msg(num_commen_lines < 1024, "Expected < 1024 lines.");
                chk_msg(fgets(line, 1024, f), "Can't read header-line.");
                if (line[0] == '#') {
                        num_commen_lines += 1u;
                } else {
                        break;
                }
        }
        num_cols = atoi(line);
        chk_msg(fgets(line, 1024, f), "Can't read header-line.");
        num_rows = atoi(line);
        chk_msg(fgets(line, 1024, f), "Can't read header-line.");
        chk_msg(strcmp(line, "255\n") == 0, "Expected 8bit range '255'.");
        chk_mem(mliImage_malloc(img, num_cols, num_rows));
        for (row = 0; row < img->num_rows; row++) {
                for (col = 0; col < img->num_cols; col++) {
                        uint8_t r, g, b;
                        struct mliColor color;
                        chk_fread(&r, sizeof(uint8_t), 1u, f);
                        chk_fread(&g, sizeof(uint8_t), 1u, f);
                        chk_fread(&b, sizeof(uint8_t), 1u, f);
                        color.r = (float)r;
                        color.g = (float)g;
                        color.b = (float)b;
                        mliImage_set(img, col, row, color);
                }
        }
        return 1;
error:
        mliImage_free(img);
        return 0;
}

int mliImage_malloc_from_path(struct mliImage *img, const char *path)
{
        FILE *f;
        f = fopen(path, "rb");
        chk_msgf(f, ("Can't open path '%s'.", path));
        chk_msg(mliImage_malloc_fread(img, f), "Can't read ppm from file.");
        chk_msg(!feof(f), "Unexpected end-of-file.");
        chk_msg(!ferror(f), "File error.");
        fclose(f);
        return 1;
error:
        mliImage_free(img);
        fclose(f);
        return 0;
}

int mliImage_fwrite(const struct mliImage *img, FILE *f)
{

        uint32_t col;
        uint32_t row;
        chk(fprintf(f, "P6\n"));
        chk(fprintf(f, "# merlict_c89\n"));
        chk(
                fprintf(f,
                        "# MLI_VERSION %d.%d.%d\n",
                        MLI_VERSION_MAYOR,
                        MLI_VERSION_MINOR,
                        MLI_VERSION_PATCH));
        chk(fprintf(f, "%d\n", img->num_cols));
        chk(fprintf(f, "%d\n", img->num_rows));
        chk(fprintf(f, "255\n"));
        for (row = 0; row < img->num_rows; row++) {
                for (col = 0; col < img->num_cols; col++) {
                        struct mliColor color = mliImage_at(img, col, row);
                        struct mliColor out = mliColor_truncate_to_uint8(color);
                        uint8_t r = (uint8_t)out.r;
                        uint8_t g = (uint8_t)out.g;
                        uint8_t b = (uint8_t)out.b;
                        chk_fwrite(&r, sizeof(uint8_t), 1u, f);
                        chk_fwrite(&g, sizeof(uint8_t), 1u, f);
                        chk_fwrite(&b, sizeof(uint8_t), 1u, f);
                }
        }
        return 1;
error:
        return 0;
}

int mliImage_write_to_path(const struct mliImage *img, const char *path)
{
        FILE *f;
        f = fopen(path, "wb");
        chk_msgf(f, ("Can't open path '%s'.", path));
        chk_msg(mliImage_fwrite(img, f), "Can't write ppm to file.");
        chk_msg(!feof(f), "Unexpected end-of-file.");
        chk_msg(!ferror(f), "File error.");
        fclose(f);
        return 1;
error:
        fclose(f);
        return 0;
}

/* mliImage_print */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mliImage_print(const struct mliImage *img, const uint64_t print_mode)
{
        const uint64_t num_symbols = 0;
        mliImage_print_chars(img, NULL, NULL, NULL, num_symbols, print_mode);
}

void mliImage_print_chars(
        const struct mliImage *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols,
        const uint64_t print_mode)
{
        if (print_mode == MLI_ANSI_ESCAPE_COLOR) {
                mliImage_print_ansi_escape_chars(
                        img, symbols, rows, cols, num_symbols);
        } else {
                mliImage_print_ascii_chars(
                        img, symbols, rows, cols, num_symbols);
                if (print_mode != MLI_ASCII_MONOCHROME) {
                        fprintf(stderr,
                                "Do not know print_mode %u\n",
                                (uint32_t)print_mode);
                }
        }
        return;
}

void mliImage_print_ansi_escape_chars(
        const struct mliImage *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols)
{
        uint32_t col, row, sym;
        char symbol;
        for (row = 0; row < img->num_rows; row++) {
                for (col = 0; col < img->num_cols; col++) {
                        struct mliColor color = mliImage_at(img, col, row);
                        struct mliColor out = mliColor_truncate_to_uint8(color);
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

void mliImage_print_ascii_chars(
        const struct mliImage *img,
        const char *symbols,
        const uint64_t *rows,
        const uint64_t *cols,
        const uint64_t num_symbols)
{
        uint32_t col, row, sym;
        char symbol;
        char chars_with_ascending_fill[16] = {' ',
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
        for (row = 0; row < img->num_rows; row++) {
                for (col = 0; col < img->num_cols; col++) {
                        struct mliColor color = mliImage_at(img, col, row);
                        struct mliColor out = mliColor_truncate_to_uint8(color);
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

/* mliIntersection */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliIntersection mliIntersection_init(void)
{
        struct mliIntersection psec;
        psec.geometry_id = mliGeometryId_init();
        psec.position_local = mliVec_init(0.0, 0.0, 0.0);
        psec.distance_of_ray = DBL_MAX;
        return psec;
}

/* mliIntersectionSurfaceNormal */
/* ---------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliIntersectionSurfaceNormal mliIntersectionSurfaceNormal_init(void)
{
        struct mliIntersectionSurfaceNormal isec;
        isec.geometry_id = mliGeometryId_init();
        isec.position = mliVec_init(0.0, 0.0, 0.0);
        isec.surface_normal = mliVec_init(0.0, 0.0, 1.0);
        isec.position_local = mliVec_init(0.0, 0.0, 0.0);
        isec.surface_normal_local = mliVec_init(0.0, 0.0, 1.0);
        isec.distance_of_ray = DBL_MAX;
        isec.from_outside_to_inside = 1;
        return isec;
}

/* mliIo */
/* ----- */

/* Copyright 2018-2023 Sebastian Achim Mueller */

struct mliIo mliIo_init(void)
{
        struct mliIo byt;
        byt.cstr = NULL;
        byt.capacity = 0u;
        byt.size = 0u;
        byt.pos = 0u;
        return byt;
}

void mliIo_free(struct mliIo *byt)
{
        free(byt->cstr);
        (*byt) = mliIo_init();
}

int mliIo_malloc_capacity(struct mliIo *byt, const uint64_t capacity)
{
        mliIo_free(byt);
        byt->capacity = MLI_MAX2(2u, capacity);
        byt->size = 0u;
        chk_malloc(byt->cstr, unsigned char, byt->capacity);
        memset(byt->cstr, '\0', byt->capacity);
        return 1;
error:
        return 0;
}

int mliIo_malloc(struct mliIo *byt)
{
        chk(mliIo_malloc_capacity(byt, 0u));
        return 1;
error:
        return 0;
}

int mliIo_putc(struct mliIo *byt, const unsigned char c)
{
        /* 'puts' a single byte (char) to the BytesIo-buffer */
        const uint64_t new_size = byt->size + 1u;

        if (byt->cstr == NULL) {
                chk(mliIo_malloc(byt));
        }

        if (new_size >= byt->capacity) {
                const uint64_t min_new_capacity =
                        MLI_MAX2(new_size, 2 * byt->capacity);

                byt->capacity = min_new_capacity;
                byt->cstr = (unsigned char *)realloc(
                        (void *)byt->cstr, byt->capacity * sizeof(char));
                memset(&byt->cstr[byt->size], '\0', byt->capacity - byt->size);
                chk_mem(byt->cstr);
        }
        byt->cstr[byt->size] = c;
        byt->size = new_size;
        byt->pos += 1;

        return 1;
error:
        return 0;
}

int mliIo_putchar(struct mliIo *byt, const char c)
{
        return mliIo_putc(byt, (unsigned char)c);
}

int mliIo_getc(struct mliIo *byt)
{
        int c;
        if (byt->pos >= byt->size) {
                return EOF;
        }
        c = byt->cstr[byt->pos];
        byt->pos += 1;
        return c;
}

int64_t mliIo_write(
        struct mliIo *byt,
        const void *ptr,
        const uint64_t size,
        const uint64_t count)
{
        const uint64_t block_size = size * count;
        unsigned char *block = (unsigned char *)ptr;

        uint64_t i;
        for (i = 0u; i < block_size; i++) {
                chk_msg(mliIo_putc(byt, block[i]), "Failed to put byte");
        }
        return (i + 1u) / size;
error:
        return -1;
}

int64_t mliIo_read(
        struct mliIo *byt,
        const void *ptr,
        const uint64_t size,
        const uint64_t count)
{
        const uint64_t block_size = size * count;
        unsigned char *block = (unsigned char *)ptr;

        uint64_t i;
        for (i = 0u; i < block_size; i++) {
                int c = mliIo_getc(byt);

                if (c == EOF) {
                        return EOF;
                } else {
                        block[i] = (unsigned char)c;
                }
        }
        return (i + 1u) / size;
}

uint64_t mliIo_ftell(struct mliIo *byt) { return byt->pos; }

void mliIo_rewind(struct mliIo *byt) { byt->pos = 0u; }

int64_t mliIo_printf(struct mliIo *byt, const char *format, ...)
{
        uint64_t i;
        char tmp[256];
        va_list args;
        memset(tmp, '\0', sizeof(tmp));

        va_start(args, format);

        /* DANGER, THIS WILL 100% BUFFEROVERFLOW  */
        /* PROBLEM IS THERE IS NO vsnprintf (note the 'n') in c89 */
        /* SO SEFAULT-PARADE IT IS UNTIL I WRITE MY OWN snprintf */
        vsprintf(tmp, format, args);

        for (i = 0; i < 256; i++) {
                if (tmp[i] == '\0') {
                        break;
                }
                chk(mliIo_putchar(byt, tmp[i]));
        }
        va_end(args);
        return (int64_t)i;
error:
        va_end(args);
        return -1;
}

int mliStr_convert_line_break_CRLF_CR_to_LF(
        struct mliStr *dst,
        const struct mliStr *src)
{
        uint64_t i = 0;
        struct mliIo sdst = mliIo_init();
        chk(mliIo_malloc(&sdst));

        while (i < src->length) {
                if (mli_cstr_is_CRLF((char *)&src->cstr[i])) {
                        chk(mliIo_putchar(&sdst, '\n'));
                        i += 2;
                } else if (mli_cstr_is_CR((char *)&src->cstr[i])) {
                        chk(mliIo_putchar(&sdst, '\n'));
                        i += 1;
                } else {
                        chk(mliIo_putchar(&sdst, src->cstr[i]));
                        i += 1;
                }
        }

        chk(mliStr_malloc(dst, sdst.size));
        strncpy(dst->cstr, (char *)sdst.cstr, sdst.size);

        mliIo_free(&sdst);
        return 1;
error:
        mliIo_free(&sdst);
        mliStr_free(dst);
        return 0;
}

int mliIo_malloc_from_path(struct mliIo *byt, const char *path)
{
        int c = EOF;
        FILE *f = fopen(path, "rt");
        chk_msg(f, "Failed to open file.");
        chk_msg(mliIo_malloc(byt), "Can not malloc string.");
        c = getc(f);
        while (c != EOF) {
                chk_msg(mliIo_putchar(byt, c), "Failed to push back char.");
                c = getc(f);
        }
        fclose(f);
        return 1;
error:
        if (f != NULL)
                fclose(f);
        mliIo_free(byt);
        return 0;
}

int mliIo_read_to_path(struct mliIo *byt, const char *path)
{
        int c = 1;
        FILE *f = fopen(path, "w");
        chk_msg(f != NULL, "Failed to open path.");
        while (c != EOF) {
                c = mliIo_getc(byt);
                chk(fputc(c, f));
        }
        fclose(f);
        return 1;
error:
        fclose(f);
        return 0;
}

int64_t mliIo_malloc_cstr(struct mliIo *byt, const char *s)
{
        const uint64_t slen = strlen(s);
        uint64_t i;
        for (i = 0; i < slen; i++) {
                chk_msg(mliIo_putchar(byt, s[i]), "Failed to push back char");
        }
        return i;
error:
        return 0;
}

int mli_readline(
        struct mliIo *stream,
        struct mliStr *line,
        const char delimiter)
{
        struct mliIo buf = mliIo_init();
        chk(mliIo_malloc(&buf));

        while (stream->pos < stream->size) {
                const int c = mliIo_getc(stream);
                if (c == '\0') {
                        break;
                } else if (c == delimiter) {
                        break;
                } else {
                        chk(mliIo_putchar(&buf, c));
                }
        }

        mliIo_rewind(&buf);
        chk(mliStr_malloc(line, buf.size));
        strcpy(line->cstr, (char *)buf.cstr);

        mliIo_free(&buf);
        return 1;
error:
        mliIo_free(&buf);
        return 0;
}

int mli_path_strip_this_dir(const struct mliStr *src, struct mliStr *dst)
{
        uint64_t i = 0;
        struct mliStr cpysrc = mliStr_init();
        chk_msg(src->cstr, "Expected src-string to be allocated.");
        chk_msg(mliStr_malloc_copy(&cpysrc, src), "Can not copy input.");
        mliStr_free(dst);

        if (cpysrc.cstr == NULL) {
                return 1;
        }

        while (i + 1 < cpysrc.length) {
                if (cpysrc.cstr[i] == '.' && cpysrc.cstr[i + 1] == '/') {
                        i += 2;
                } else {
                        break;
                }
        }
        chk(mliStr_malloc(dst, cpysrc.length - i));
        strcpy(dst->cstr, &cpysrc.cstr[i]);
        mliStr_free(&cpysrc);
        return 1;
error:
        mliStr_free(&cpysrc);
        mliStr_free(dst);
        return 0;
}

int mli_path_basename(const struct mliStr *src, struct mliStr *dst)
{
        int64_t pos_last_del = -1;
        mliStr_free(dst);
        chk_msg(src->cstr != NULL, "Expected src-path to be allocated");

        pos_last_del = mliStr_rfind(src, '/');

        if (pos_last_del < 0) {
                chk(mliStr_mallocf(dst, src->cstr));
        } else {
                chk(mliStr_mallocf(dst, &src->cstr[pos_last_del + 1]));
        }
        return 1;
error:
        mliStr_free(dst);
        return 0;
}

int mli_path_splitext(
        const struct mliStr *src,
        struct mliStr *dst,
        struct mliStr *ext)
{
        int64_t p = -1;
        int64_t d = -1;
        struct mliStr tmp = mliStr_init();
        chk_msg(src->cstr != NULL, "Expected src-path to be allocated");
        chk(mliStr_malloc_copy(&tmp, src));

        mliStr_free(dst);
        mliStr_free(ext);

        p = mliStr_rfind(&tmp, '.');
        d = mliStr_rfind(&tmp, '/');

        if (p <= 0 || d > p ||
            ((d + 1 == p) && (p + 1 < (int64_t)tmp.length))) {
                chk(mliStr_mallocf(dst, tmp.cstr));
                chk(mliStr_mallocf(ext, ""));
        } else {
                chk(mliStr_malloc(dst, p));
                strncpy(dst->cstr, tmp.cstr, p);
                chk(mliStr_malloc(ext, tmp.length - p));
                strncpy(ext->cstr, &tmp.cstr[p + 1], tmp.length - p);
        }

        mliStr_free(&tmp);
        return 1;
error:
        mliStr_free(&tmp);
        mliStr_free(dst);
        mliStr_free(ext);
        return 0;
}

int mli_line_viewer_write_line_match(
        struct mliIo *f,
        const int64_t line_number,
        const int64_t line_number_of_interest)
{
        chk(mliIo_printf(f, "% 6d", (int32_t)line_number));
        if (line_number == line_number_of_interest) {
                chk(mliIo_printf(f, "->|  "));
        } else {
                chk(mliIo_printf(f, "  |  "));
        }
        return 1;
error:
        return 0;
}

int mli_line_viewer_write(
        struct mliIo *f,
        const struct mliStr *text,
        const uint64_t line_number,
        const uint64_t line_radius)
{
        int64_t _line_number = (int64_t)line_number;
        int64_t _line_radius = (int64_t)line_radius;
        int64_t line_start = MLI_MAX2(_line_number - _line_radius, 1);
        int64_t line_stop = line_number + line_radius;
        int64_t line = 1;
        uint64_t i = 0;

        chk_msg(line_radius > 1, "Expected line_radius > 1.");

        chk(mliIo_printf(f, "  line     text\n"));
        chk(mliIo_printf(f, "        |\n"));

        while (i < text->length && text->cstr[i]) {
                int prefix = (line + 1 >= line_start) && (line < line_stop);
                int valid = (line >= line_start) && (line <= line_stop);
                if (text->cstr[i] == '\n') {
                        line++;
                }
                if (prefix && i == 0) {
                        chk(mli_line_viewer_write_line_match(
                                f, line, _line_number));
                }
                if (valid) {
                        chk(mliIo_putchar(f, text->cstr[i]));
                }
                if (prefix && text->cstr[i] == '\n') {
                        chk(mli_line_viewer_write_line_match(
                                f, line, _line_number));
                }
                i++;
        }
        chk(mliIo_putchar(f, '\n'));

        return 1;
error:
        return 0;
}

/* mliMagicId */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliMagicId mliMagicId_init(void)
{
        struct mliMagicId magic;
        memset(magic.word, '\0', sizeof(magic.word));
        magic.mayor = MLI_VERSION_MAYOR;
        magic.minor = MLI_VERSION_MINOR;
        magic.patch = MLI_VERSION_PATCH;
        return magic;
}

int mliMagicId_set(struct mliMagicId *magic, const char *word)
{
        uint64_t i, len;
        (*magic) = mliMagicId_init();
        chk_msg(strlen(word) < sizeof(magic->word),
                "Expected magic word to be shorter.");

        len = MLI_MIN2(strlen(word), sizeof(magic->word));

        for (i = 0; i < len; i++) {
                magic->word[i] = word[i];
        }
        while (i < sizeof(magic->word)) {
                magic->word[i] = '\0';
                i += 1;
        }
        return 1;
error:
        return 0;
}

int mliMagicId_has_word(const struct mliMagicId *magic, const char *word)
{
        uint64_t i, len;
        chk_msg(strlen(word) < sizeof(magic->word),
                "Expected magic word to be shorter.");

        len = MLI_MIN2(strlen(word), sizeof(magic->word));

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
error:
        return 0;
}

void mliMagicId_warn_version(const struct mliMagicId *magic)
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

/* mliMat */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mliMat_set(struct mliMat *a, uint64_t col, uint64_t row, const double v)
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

double mliMat_get(const struct mliMat *a, uint64_t col, uint64_t row)
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

struct mliMat mliMat_unity(void)
{
        struct mliMat u;
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

int mliMat_equal_margin(
        const struct mliMat a,
        const struct mliMat b,
        const double margin)
{
        int i, j;
        for (i = 0; i < 3; i++) {
                for (j = 0; j < 3; j++) {
                        const double diff = fabs(
                                mliMat_get(&a, i, j) - mliMat_get(&b, i, j));
                        if (diff > margin) {
                                return 0;
                        }
                }
        }
        return 1;
}

struct mliMat mliMat_init_tait_bryan(
        const double rx,
        const double ry,
        const double rz)
{
        struct mliMat rot;
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

struct mliMat mliMat_init_axis_angle(
        const struct mliVec axis,
        const double angle)
{
        struct mliMat rot;
        const double norm = mliVec_norm(axis);
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

struct mliMat mliMat_init_columns(
        const struct mliVec c0,
        const struct mliVec c1,
        const struct mliVec c2)
{
        struct mliMat m;
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

struct mliMat mliMat_covariance(
        const struct mliVec *vecs,
        const uint64_t num_vecs,
        const struct mliVec vecs_mean)
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
        struct mliVec v;
        struct mliMat cov;

        for (i = 0; i < num_vecs; i++) {
                v = mliVec_substract(vecs[i], vecs_mean);
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

struct mliMat mliMat_transpose(const struct mliMat m)
{
        struct mliMat t;
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

struct mliMat mliMat_multiply(const struct mliMat x, const struct mliMat y)
{
        struct mliMat p;
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

struct mliMat mliMat_minor(const struct mliMat x, const int d)
{
        struct mliMat m = x;
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

struct mliMat mliMat_vector_outer_product(const struct mliVec v)
{
        struct mliMat x;
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

void mliMat_qr_decompose(
        const struct mliMat m,
        struct mliMat *q,
        struct mliMat *r)
{
        /* housholder */
        /* ========== */

        struct mliMat q0, q1;
        struct mliMat z = m;

        struct mliVec e;
        struct mliVec x;
        double xnorm;

        /* column 0 */
        /* -------- */
        z = mliMat_minor(z, 0);
        x = mliVec_init(z.r00, z.r10, z.r20);
        xnorm = mliVec_norm(x);
        if (m.r00 > 0) {
                xnorm = -xnorm;
        };

        e = mliVec_init(xnorm, 0.0, 0.0);
        e = mliVec_add(x, e);
        e = mliVec_normalized(e);

        q0 = mliMat_vector_outer_product(e);
        z = mliMat_multiply(q0, z);

        /* column 1 */
        /* -------- */
        z = mliMat_minor(z, 1);
        x = mliVec_init(z.r01, z.r11, z.r21);
        xnorm = mliVec_norm(x);
        if (m.r11 > 0) {
                xnorm = -xnorm;
        };

        e = mliVec_init(0.0, xnorm, 0.0);
        e = mliVec_add(x, e);
        e = mliVec_normalized(e);

        q1 = mliMat_vector_outer_product(e);
        z = mliMat_multiply(q1, z);

        /* finalize */
        /* ======== */
        (*q) = q0;
        (*r) = mliMat_multiply(q0, m);

        (*q) = mliMat_multiply(q1, *q);

        (*r) = mliMat_multiply(*q, m);
        (*q) = mliMat_transpose(*q);
}

int mliMat_has_shurform(const struct mliMat m, const double margin)
{
        if (fabs(m.r10) > margin)
                return 0;
        if (fabs(m.r20) > margin)
                return 0;
        if (fabs(m.r21) > margin)
                return 0;
        return 1;
}

void mliMat_find_eigenvalues(
        const struct mliMat a,
        double *e0,
        double *e1,
        double *e2,
        const double margin,
        const uint64_t max_num_iterations)
{
        struct mliMat ak = a;
        struct mliMat qq = mliMat_unity();
        uint64_t k = 0;

        while (k < max_num_iterations) {
                struct mliMat q, r;
                k += 1;
                mliMat_qr_decompose(ak, &q, &r);
                ak = mliMat_multiply(r, q);
                qq = mliMat_multiply(qq, q);

                if (mliMat_has_shurform(ak, margin)) {
                        break;
                }
        }
        (*e0) = ak.r00;
        (*e1) = ak.r11;
        (*e2) = ak.r22;
}

int mliMat_find_eigenvector_for_eigenvalue(
        struct mliMat A,
        const double eigen_value,
        struct mliVec *eigen_vector,
        const double tolerance)
{
        int pivots[4] = {0, 0, 0, 0};
        struct mliVec right_hand_side = {1.0, 1.0, 1.0};

        chk_msg(tolerance > 0.0, "Expected tolerance > 0.0");

        A.r00 -= eigen_value;
        A.r11 -= eigen_value;
        A.r22 -= eigen_value;

        chk_msg(mliMat_lup_decompose(&A, pivots, tolerance),
                "Can not decompose LU-matices");
        mliMat_lup_solve(&A, pivots, &right_hand_side, eigen_vector);
        (*eigen_vector) = mliVec_normalized(*eigen_vector);

        return 1;
error:
        return 0;
}

int mliMat_lup_decompose(struct mliMat *A, int *P, const double tolerance)
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
                        absA = fabs(mliMat_get(A, k, i));
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
                        tmp[0] = mliMat_get(A, i, 0);
                        tmp[1] = mliMat_get(A, i, 1);
                        tmp[2] = mliMat_get(A, i, 2);

                        mliMat_set(A, i, 0, mliMat_get(A, imax, 0));
                        mliMat_set(A, i, 1, mliMat_get(A, imax, 1));
                        mliMat_set(A, i, 2, mliMat_get(A, imax, 2));

                        mliMat_set(A, imax, 0, tmp[0]);
                        mliMat_set(A, imax, 1, tmp[1]);
                        mliMat_set(A, imax, 2, tmp[2]);

                        /* counting pivots starting from 3 (for determinant) */
                        P[3]++;
                }

                for (j = i + 1; j < 3; j++) {
                        double tmp = mliMat_get(A, j, i) / mliMat_get(A, i, i);
                        mliMat_set(A, j, i, tmp);

                        for (k = i + 1; k < 3; k++) {
                                double tmp =
                                        (mliMat_get(A, j, k) -
                                         mliMat_get(A, j, i) *
                                                 mliMat_get(A, i, k));
                                mliMat_set(A, j, k, tmp);
                        }
                }
        }

        return 1; /* decomposition done */
}

void mliMat_lup_solve(
        const struct mliMat *A,
        const int *P,
        const struct mliVec *b,
        struct mliVec *x)
{
        /*
        Solve linear equation A*x = b for x.

        Parameters
        ----------
        A and P
                Filled in mliMat_lup_decompose.
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
                mliVec_set(x, i, mliVec_get(b, P[i]));
                for (k = 0; k < i; k++) {
                        const double tmp =
                                (mliVec_get(x, i) -
                                 mliMat_get(A, i, k) * mliVec_get(x, k));
                        mliVec_set(x, i, tmp);
                }
        }

        i = N - 1;
        for (idown = 0; idown < N; idown++) {
                double tmp;
                for (k = i + 1; k < N; k++) {
                        const double tmp =
                                (mliVec_get(x, i) -
                                 mliMat_get(A, i, k) * mliVec_get(x, k));
                        mliVec_set(x, i, tmp);
                }
                tmp = mliVec_get(x, i) / mliMat_get(A, i, i);
                mliVec_set(x, i, tmp);
                i--;
        }
}

/* mliMaterials */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliMaterialsCapacity mliMaterialsCapacity_init(void)
{
        struct mliMaterialsCapacity cap;
        cap.num_surfaces = 0;
        cap.num_media = 0;
        cap.num_boundary_layers = 0;
        return cap;
}

struct mliMaterials mliMaterials_init(void)
{
        struct mliMaterials res;
        res.default_medium = 0u;

        res.num_media = 0u;
        res.media = NULL;
        res.medium_names = NULL;

        res.num_surfaces = 0u;
        res.surfaces = NULL;
        res.surface_names = NULL;

        res.num_boundary_layers = 0u;
        res.boundary_layers = NULL;
        res.boundary_layer_names = NULL;
        return res;
}

void mliMaterials_free(struct mliMaterials *res)
{
        uint64_t i;
        for (i = 0; i < res->num_media; i++) {
                mliMedium_free(&(res->media[i]));
        }
        free(res->media);
        free(res->medium_names);

        for (i = 0; i < res->num_surfaces; i++) {
                mliSurface_free(&(res->surfaces[i]));
        }
        free(res->surfaces);
        free(res->surface_names);

        free(res->boundary_layers);
        free(res->boundary_layer_names);

        (*res) = mliMaterials_init();
}

int mliMaterials_malloc(
        struct mliMaterials *res,
        const struct mliMaterialsCapacity rescap)
{
        uint64_t i;
        mliMaterials_free(res);
        res->num_surfaces = rescap.num_surfaces;
        res->num_media = rescap.num_media;
        res->num_boundary_layers = rescap.num_boundary_layers;

        chk_malloc(res->media, struct mliMedium, res->num_media);
        chk_malloc(res->medium_names, struct mliName, res->num_media);
        for (i = 0; i < res->num_media; i++) {
                res->media[i] = mliMedium_init();
                res->medium_names[i] = mliName_init();
        }

        chk_malloc(res->surfaces, struct mliSurface, res->num_surfaces);
        chk_malloc(res->surface_names, struct mliName, res->num_surfaces);
        for (i = 0; i < res->num_surfaces; i++) {
                res->surfaces[i] = mliSurface_init();
                res->surface_names[i] = mliName_init();
        }

        chk_malloc(
                res->boundary_layers,
                struct mliBoundaryLayer,
                res->num_boundary_layers);
        chk_malloc(
                res->boundary_layer_names,
                struct mliName,
                res->num_boundary_layers);
        for (i = 0; i < res->num_boundary_layers; i++) {
                res->boundary_layer_names[i] = mliName_init();
        }

        return 1;
error:
        mliMaterials_free(res);
        return 0;
}

void mliMaterials_info_fprint(FILE *f, const struct mliMaterials *res)
{
        uint32_t i = 0;
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
        fprintf(f, "%12s ", "absorbtion");
        fprintf(f, "%12s ", "refraction");
        fprintf(f, "%12s ", "default");
        fprintf(f, "\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        for (i = 0; i < res->num_media; i++) {
                fprintf(f, "    ");
                fprintf(f, "% 3d ", i);
                fprintf(f, "%24s ", res->medium_names[i].cstr);
                fprintf(f, "%12d ", res->media[i].absorbtion.num_points);
                fprintf(f, "%12d ", res->media[i].refraction.num_points);
                if (i == res->default_medium) {
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
        fprintf(f, "    ");
        fprintf(f, "%48s ", "spec.");
        fprintf(f, "%6s ", "diff.");
        fprintf(f, "%10s ", "color");
        fprintf(f, "\n");
        fprintf(f, "    ");
        fprintf(f, "%3s ", "#");
        fprintf(f, "%24s ", "name");
        fprintf(f, "%12s ", "material");
        fprintf(f, "%6s ", "refl.");
        fprintf(f, "%6s ", "refl.");
        fprintf(f, "%10s ", "r,g,b");
        fprintf(f, "\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        for (i = 0; i < res->num_surfaces; i++) {
                fprintf(f, "    ");
                fprintf(f, "% 3d ", i);
                fprintf(f, "%24s ", res->surface_names[i].cstr);
                if (res->surfaces[i].material == MLI_MATERIAL_TRANSPARENT) {
                        fprintf(f, "%12s ", "transparent");
                } else if (res->surfaces[i].material == MLI_MATERIAL_PHONG) {
                        fprintf(f, "%12s ", "Phong");
                } else {
                        fprintf(f, "%12s ", "UNKNOWN");
                }
                fprintf(f,
                        "%6d ",
                        res->surfaces[i].specular_reflection.num_points);
                fprintf(f,
                        "%6d ",
                        res->surfaces[i].diffuse_reflection.num_points);
                fprintf(f,
                        "   %3d,%3d,%3d ",
                        (int)res->surfaces[i].color.r,
                        (int)res->surfaces[i].color.g,
                        (int)res->surfaces[i].color.b);
                fprintf(f, "\n");
        }
        fprintf(f, "\n");

        fprintf(f, "    boundary_layers\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");
        fprintf(f, "    ");
        fprintf(f, "%36s ", "inner");
        fprintf(f, "%17s ", "outer");
        fprintf(f, "\n");
        fprintf(f, "    ");
        fprintf(f, "%3s ", "#");
        fprintf(f, "%24s ", "name");
        fprintf(f, "%8s ", "medium");
        fprintf(f, "%8s ", "surface");
        fprintf(f, "%8s ", "medium");
        fprintf(f, "%8s ", "surface");
        fprintf(f, "\n");
        fprintf(f, "    ");
        for (i = 0; i < 70; i++) {
                fprintf(f, "-");
        }
        fprintf(f, "\n");

        for (i = 0; i < res->num_boundary_layers; i++) {
                fprintf(f, "    ");
                fprintf(f, "% 3d ", i);
                fprintf(f, "%24s ", res->boundary_layer_names[i].cstr);

                fprintf(f, "%8d ", res->boundary_layers[i].inner.medium);
                fprintf(f, "%8d ", res->boundary_layers[i].inner.surface);

                fprintf(f, "%8d ", res->boundary_layers[i].outer.medium);
                fprintf(f, "%8d ", res->boundary_layers[i].outer.surface);

                fprintf(f, "\n");
        }
}

/* mliMaterials_equal */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliMaterials_media_equal(
        const struct mliMaterials *a,
        const struct mliMaterials *b)
{
        uint32_t i = 0u;
        chk_msg(a->num_media == b->num_media, "Different number of media.");
        for (i = 0; i < a->num_media; i++) {
                chk_msg(mliMedium_equal(&a->media[i], &b->media[i]),
                        "Different Medium.");
                chk_msg(mliName_equal(&a->medium_names[i], &b->medium_names[i]),
                        "Different medium-name.");
        }
        return 1;
error:
        fprintf(stderr, "In materials.media[%u].\n", i);
        return 0;
}

int mliMaterials_surfaces_equal(
        const struct mliMaterials *a,
        const struct mliMaterials *b)
{
        uint32_t i = 0u;
        chk_msg(a->num_surfaces == b->num_surfaces,
                "Different number of surfaces.");
        for (i = 0; i < a->num_surfaces; i++) {
                chk_msg(mliSurface_equal(&a->surfaces[i], &b->surfaces[i]),
                        "Different Surface.");
                chk_msg(mliName_equal(
                                &a->surface_names[i], &b->surface_names[i]),
                        "Different surface-name.");
        }
        return 1;
error:
        fprintf(stderr, "In materials.surfaces[%u].\n", i);
        return 0;
}

int mliMaterials_boundary_layers_equal(
        const struct mliMaterials *a,
        const struct mliMaterials *b)
{
        uint32_t i = 0u;
        chk_msg(a->num_boundary_layers == b->num_boundary_layers,
                "Different number of boundary_layers.");
        for (i = 0; i < a->num_boundary_layers; i++) {
                chk_msg(mliBoundaryLayer_equal(
                                a->boundary_layers[i], b->boundary_layers[i]),
                        "Different boundary_layer.");
                chk_msg(mliName_equal(
                                &a->boundary_layer_names[i],
                                &b->boundary_layer_names[i]),
                        "Different boundary_layer-name.");
        }
        return 1;
error:
        fprintf(stderr, "In materials.boundary_layers[%u].\n", i);
        return 0;
}

int mliMaterials_equal(
        const struct mliMaterials *a,
        const struct mliMaterials *b)
{
        chk_msg(a->default_medium == b->default_medium,
                "Different default_medium.");
        chk_msg(mliMaterials_media_equal(a, b), "Different media.");
        chk_msg(mliMaterials_surfaces_equal(a, b), "Different surfaces.");
        chk_msg(mliMaterials_boundary_layers_equal(a, b),
                "Different boundary_layers.");
        return 1;
error:
        return 0;
}

/* mliMaterials_serialize */
/* ---------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliMaterials_fwrite(const struct mliMaterials *res, FILE *f)
{
        uint64_t i;
        struct mliMagicId magic = mliMagicId_init();

        /* magic identifier */
        chk(mliMagicId_set(&magic, "mliMaterials"));
        chk_fwrite(&magic, sizeof(struct mliMagicId), 1u, f);

        chk_fwrite(&res->num_media, sizeof(uint64_t), 1u, f);
        chk_fwrite(&res->num_surfaces, sizeof(uint64_t), 1u, f);
        chk_fwrite(&res->num_boundary_layers, sizeof(uint64_t), 1u, f);

        for (i = 0; i < res->num_media; i++) {
                chk_fwrite(
                        &res->medium_names[i].cstr,
                        sizeof(char),
                        MLI_NAME_CAPACITY,
                        f);
                chk(mliMedium_fwrite(&res->media[i], f));
        }
        for (i = 0; i < res->num_surfaces; i++) {
                chk_fwrite(
                        &res->surface_names[i].cstr,
                        sizeof(char),
                        MLI_NAME_CAPACITY,
                        f);
                chk(mliSurface_fwrite(&res->surfaces[i], f));
        }

        chk_fwrite(
                res->boundary_layers,
                sizeof(struct mliBoundaryLayer),
                res->num_boundary_layers,
                f);
        chk_fwrite(
                res->boundary_layer_names,
                sizeof(struct mliName),
                res->num_boundary_layers,
                f);

        chk_fwrite(&res->default_medium, sizeof(uint64_t), 1, f);

        return 1;
error:
        return 0;
}

int mliMaterials_malloc_fread(struct mliMaterials *res, FILE *f)
{
        uint64_t i;
        struct mliMagicId magic;
        struct mliMaterialsCapacity cap;

        /* magic identifier */
        chk_fread(&magic, sizeof(struct mliMagicId), 1u, f);
        chk(mliMagicId_has_word(&magic, "mliMaterials"));
        mliMagicId_warn_version(&magic);

        chk_fread(&cap.num_media, sizeof(uint64_t), 1u, f);
        chk_fread(&cap.num_surfaces, sizeof(uint64_t), 1u, f);
        chk_fread(&cap.num_boundary_layers, sizeof(uint64_t), 1u, f);

        chk(mliMaterials_malloc(res, cap));

        /* payload */
        for (i = 0; i < res->num_media; i++) {
                chk_fread(
                        &res->medium_names[i].cstr,
                        sizeof(char),
                        MLI_NAME_CAPACITY,
                        f);
                chk_msg(mliMedium_malloc_fread(&res->media[i], f),
                        "Failed to fread Medium.");
        }
        for (i = 0; i < res->num_surfaces; i++) {
                chk_fread(
                        &res->surface_names[i].cstr,
                        sizeof(char),
                        MLI_NAME_CAPACITY,
                        f);
                chk_msg(mliSurface_malloc_fread(&res->surfaces[i], f),
                        "Failed to fread Surface.");
        }

        chk_fread(
                res->boundary_layers,
                sizeof(struct mliBoundaryLayer),
                res->num_boundary_layers,
                f);
        chk_fread(
                res->boundary_layer_names,
                sizeof(struct mliName),
                res->num_boundary_layers,
                f);

        chk_fread(&res->default_medium, sizeof(uint64_t), 1, f);

        return 1;
error:
        return 0;
}

/* mliMaterials_valid */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliMaterials_valid_media(const struct mliMaterials *materials)
{
        uint32_t i = 0u;
        for (i = 0; i < materials->num_media; i++) {
                chk_msg(mliFunc_is_valid(&materials->media[i].refraction),
                        "Expected refraction of medium to be valid.");
                chk_msg(mliFunc_is_valid(&materials->media[i].absorbtion),
                        "Expected refraction of medium to be valid.");
        }
        return 1;
error:
        fprintf(stderr, "In materials.media[%u]\n", i);
        return 0;
}

int mliMaterials_valid_surfaces(const struct mliMaterials *materials)
{
        uint32_t i = 0u;
        for (i = 0; i < materials->num_surfaces; i++) {
                chk_msg((materials->surfaces[i].material ==
                         MLI_MATERIAL_PHONG) ||
                                (materials->surfaces[i].material ==
                                 MLI_MATERIAL_TRANSPARENT),
                        "Material-type is unknown.");
                chk_msg(mliFunc_is_valid(
                                &materials->surfaces[i].specular_reflection),
                        "Expected specular_reflection of surface to be valid.");
                chk_msg(mliFunc_is_valid(
                                &materials->surfaces[i].diffuse_reflection),
                        "Expected diffuse_reflection of surface to be valid.");
                chk_msg(mliColor_is_valid_8bit_range(
                                materials->surfaces[i].color),
                        "Expected 0.0 <= color < 256.0.");
        }
        return 1;
error:
        fprintf(stderr, "In materials.surface[%u]\n", i);
        return 0;
}

int mliMaterials_valid_boundary_layers(const struct mliMaterials *materials)
{
        uint32_t i = 0u;
        for (i = 0; i < materials->num_boundary_layers; i++) {
                chk_msg(materials->boundary_layers[i].inner.surface <
                                materials->num_surfaces,
                        "inner.surface is invalid.");
                chk_msg(materials->boundary_layers[i].outer.surface <
                                materials->num_surfaces,
                        "outer.surface is invalid.");
                chk_msg(materials->boundary_layers[i].inner.medium <
                                materials->num_media,
                        "inner.medium is invalid.");
                chk_msg(materials->boundary_layers[i].outer.medium <
                                materials->num_media,
                        "outer.medium is invalid.");

                chk_msg(mliName_valid(&materials->boundary_layer_names[i]),
                        "Name is invalid.");
        }
        return 1;
error:
        fprintf(stderr, "In materials.boundary_layers[%u]\n", i);
        return 0;
}

int mliMaterials_valid(const struct mliMaterials *materials)
{
        chk_msg(materials->default_medium <= materials->num_media,
                "Expected default-medium to reference a valid medium.");
        chk_msg(mliMaterials_valid_media(materials),
                "Expected media to be valid.");
        chk_msg(mliMaterials_valid_surfaces(materials),
                "Expected surfaces to be valid.");
        chk_msg(mliMaterials_valid_boundary_layers(materials),
                "Expected boundary_layers to be valid.");
        return 1;
error:
        return 0;
}

/* mliMedium */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliMedium mliMedium_init(void)
{
        struct mliMedium medium;
        medium.refraction = mliFunc_init();
        medium.absorbtion = mliFunc_init();
        return medium;
}

void mliMedium_free(struct mliMedium *medium)
{
        mliFunc_free(&medium->refraction);
        mliFunc_free(&medium->absorbtion);
        (*medium) = mliMedium_init();
}

int mliMedium_equal(const struct mliMedium *a, const struct mliMedium *b)
{
        if (!mliFunc_equal(a->refraction, b->refraction))
                return 0;
        if (!mliFunc_equal(a->absorbtion, b->absorbtion))
                return 0;
        return 1;
}

/* mliMedium_json */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliMedium_malloc_from_json_str(struct mliMedium *med, const char *json_str)
{
        struct mliJson json = mliJson_init();
        chk_msg(mliJson_malloc_from_cstr(&json, json_str),
                "Failed to read json_str to malloc medium.");
        chk_msg(mliMedium_malloc_from_json_token(med, &json, 0),
                "Failed to malloc medium from json.");
        mliJson_free(&json);
        return 1;
error:
        return 0;
}

int mliMedium_malloc_from_json_token(
        struct mliMedium *med,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t refraction_token;
        uint64_t absorbtion_token;
        chk_msg(mliJson_token_by_key(
                        json, token, "refraction", &refraction_token),
                "Expected medium to have key 'refraction', but it does not.");
        chk_msg(mliFunc_malloc_from_json_token(
                        &med->refraction, json, refraction_token + 1),
                "Failed to read medium's refraction from json.");
        chk_msg(mliJson_token_by_key(
                        json, token, "absorbtion", &absorbtion_token),
                "Expected medium to have key 'absorbtion', but it does not.");
        chk_msg(mliFunc_malloc_from_json_token(
                        &med->absorbtion, json, absorbtion_token + 1),
                "Failed to read medium's absorbtion from json.");
        return 1;
error:
        return 0;
}

/* mliMedium_serialize */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliMedium_fwrite(const struct mliMedium *med, FILE *f)
{
        struct mliMagicId magic = mliMagicId_init();
        /* magic identifier */
        chk(mliMagicId_set(&magic, "mliMedium"));
        chk_fwrite(&magic, sizeof(struct mliMagicId), 1u, f);

        chk(mliFunc_fwrite(&med->refraction, f));
        chk(mliFunc_fwrite(&med->absorbtion, f));

        return 1;
error:
        return 0;
}

int mliMedium_malloc_fread(struct mliMedium *med, FILE *f)
{
        struct mliMagicId magic;
        /* magic identifier */
        chk_fread(&magic, sizeof(struct mliMagicId), 1u, f);
        chk(mliMagicId_has_word(&magic, "mliMedium"));
        mliMagicId_warn_version(&magic);

        chk(mliFunc_malloc_fread(&med->refraction, f));
        chk(mliFunc_malloc_fread(&med->absorbtion, f));

        return 1;
error:
        return 0;
}

/* mliName */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliName mliName_init(void)
{
        struct mliName name;
        memset(name.cstr, '\0', sizeof(name.cstr));
        return name;
}

int mliName_valid(const struct mliName *name)
{
        if (name->cstr[sizeof(name->cstr) - 1] == '\0') {
                return 1;
        } else {
                return 0;
        }
}

int mliName_equal(const struct mliName *a, const struct mliName *b)
{
        uint32_t i;
        for (i = 0; i < sizeof(a->cstr); i++) {
                if (a->cstr[i] != b->cstr[i]) {
                        return 0;
                }
        }
        return 1;
}

int mliName_find_idx(
        const struct mliName *names,
        const uint64_t num_names,
        const char *key,
        uint64_t *idx)
{
        uint64_t i;
        (*idx) = 0u;
        for (i = 0; i < num_names; i++) {
                if (strcmp(names[i].cstr, key) == 0) {
                        (*idx) = i;
                        return 1;
                }
        }
        return 0;
}

/* mliObject */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliObject mliObject_init(void)
{
        struct mliObject obj;
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

void mliObject_free(struct mliObject *obj)
{
        free(obj->vertices);
        free(obj->vertex_normals);
        free(obj->faces_vertices);
        free(obj->faces_vertex_normals);
        free(obj->faces_materials);
        free(obj->material_names);
        *obj = mliObject_init();
}

int mliObject_malloc(
        struct mliObject *obj,
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

        mliObject_free(obj);
        obj->num_vertices = num_vertices;
        obj->num_vertex_normals = num_vertex_normals;
        obj->num_faces = num_faces;
        obj->num_materials = num_materials;
        chk_malloc(obj->vertices, struct mliVec, obj->num_vertices);
        chk_malloc(obj->vertex_normals, struct mliVec, obj->num_vertex_normals);
        chk_malloc(obj->faces_vertices, struct mliFace, obj->num_faces);
        chk_malloc(obj->faces_vertex_normals, struct mliFace, obj->num_faces);
        chk_malloc(obj->faces_materials, uint16_t, obj->num_faces);
        chk_malloc(obj->material_names, struct mliName, obj->num_materials);
        for (i = 0; i < obj->num_materials; i++) {
                obj->material_names[i] = mliName_init();
        }
        return 1;
error:
        return 0;
}

int mliObject_equal(const struct mliObject *a, const struct mliObject *b)
{
        uint64_t i;
        chk(a->num_vertices == b->num_vertices);
        chk(a->num_vertex_normals == b->num_vertex_normals);
        chk(a->num_faces == b->num_faces);
        chk(a->num_materials == b->num_materials);

        for (i = 0; i < a->num_vertices; i++) {
                chk(mliVec_equal(a->vertices[i], b->vertices[i]));
        }
        for (i = 0; i < a->num_vertex_normals; i++) {
                chk(mliVec_equal(a->vertex_normals[i], b->vertex_normals[i]));
        }
        for (i = 0; i < a->num_faces; i++) {
                chk(mliFace_equal(a->faces_vertices[i], b->faces_vertices[i]));
                chk(mliFace_equal(
                        a->faces_vertex_normals[i],
                        b->faces_vertex_normals[i]));
                chk(a->faces_materials[i] == b->faces_materials[i]);
        }
        for (i = 0; i < a->num_materials; i++) {
                chk(mliName_equal(
                        &a->material_names[i], &b->material_names[i]));
        }
        return 1;
error:
        return 0;
}

uint32_t mliObject_resolve_material_idx(
        const struct mliObject *obj,
        const uint32_t face_idx)
{
        assert(face_idx < obj->num_faces);
        return obj->faces_materials[face_idx];
}

/* mliObject_AABB */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliObject_face_in_local_frame_has_overlap_aabb(
        const struct mliObject *obj,
        const uint64_t face_idx,
        const struct mliAABB aabb)
{
        struct mliFace face;
        if (face_idx >= obj->num_faces) {
                return 0;
        }
        face = obj->faces_vertices[face_idx];
        if (mliTriangle_has_overlap_aabb(
                    obj->vertices[face.a],
                    obj->vertices[face.b],
                    obj->vertices[face.c],
                    aabb)) {
                return 1;
        }
        return 0;
}

int mliObject_has_overlap_aabb(
        const struct mliObject *obj,
        const struct mliHomTra local2root,
        const struct mliAABB aabb)
{
        uint64_t face_idx;
        for (face_idx = 0; face_idx < obj->num_faces; face_idx++) {
                const struct mliFace face = obj->faces_vertices[face_idx];

                const struct mliVec a_local = obj->vertices[face.a];
                const struct mliVec b_local = obj->vertices[face.b];
                const struct mliVec c_local = obj->vertices[face.c];

                const struct mliVec a_root =
                        mliHomTra_pos(&local2root, a_local);
                const struct mliVec b_root =
                        mliHomTra_pos(&local2root, b_local);
                const struct mliVec c_root =
                        mliHomTra_pos(&local2root, c_local);

                if (mliTriangle_has_overlap_aabb(
                            a_root, b_root, c_root, aabb)) {
                        return 1;
                }
        }
        return 0;
}

int mliObject_face_in_local_frame_has_overlap_aabb_void(
        const void *obj,
        const uint32_t face_idx,
        const struct mliAABB aabb)
{
        return mliObject_face_in_local_frame_has_overlap_aabb(
                (const struct mliObject *)obj, (uint64_t)face_idx, aabb);
}

struct mliAABB mliObject_aabb(
        const struct mliObject *obj,
        const struct mliHomTra local2root)
{
        struct mliAABB aabb;
        struct mliVec seed_vertex_local;
        struct mliVec seed_vertex_root;
        struct mliFace seed_face;
        uint64_t face_idx;

        if (obj->num_faces == 0) {
                aabb.lower = mliVec_init(MLI_NAN, MLI_NAN, MLI_NAN);
                aabb.upper = mliVec_init(MLI_NAN, MLI_NAN, MLI_NAN);
                return aabb;
        }

        seed_face = obj->faces_vertices[0];
        seed_vertex_local = obj->vertices[seed_face.a];
        seed_vertex_root = mliHomTra_pos(&local2root, seed_vertex_local);

        aabb.lower = seed_vertex_root;
        aabb.upper = seed_vertex_root;
        for (face_idx = 0; face_idx < obj->num_faces; face_idx++) {
                const struct mliFace face = obj->faces_vertices[face_idx];
                const struct mliVec a_local = obj->vertices[face.a];
                const struct mliVec b_local = obj->vertices[face.b];
                const struct mliVec c_local = obj->vertices[face.c];
                const struct mliVec a_root =
                        mliHomTra_pos(&local2root, a_local);
                const struct mliVec b_root =
                        mliHomTra_pos(&local2root, b_local);
                const struct mliVec c_root =
                        mliHomTra_pos(&local2root, c_local);

                aabb.lower.x = MLI_MIN2(aabb.lower.x, a_root.x);
                aabb.lower.y = MLI_MIN2(aabb.lower.y, a_root.y);
                aabb.lower.z = MLI_MIN2(aabb.lower.z, a_root.z);

                aabb.upper.x = MLI_MAX2(aabb.upper.x, a_root.x);
                aabb.upper.y = MLI_MAX2(aabb.upper.y, a_root.y);
                aabb.upper.z = MLI_MAX2(aabb.upper.z, a_root.z);

                aabb.lower.x = MLI_MIN2(aabb.lower.x, b_root.x);
                aabb.lower.y = MLI_MIN2(aabb.lower.y, b_root.y);
                aabb.lower.z = MLI_MIN2(aabb.lower.z, b_root.z);

                aabb.upper.x = MLI_MAX2(aabb.upper.x, b_root.x);
                aabb.upper.y = MLI_MAX2(aabb.upper.y, b_root.y);
                aabb.upper.z = MLI_MAX2(aabb.upper.z, b_root.z);

                aabb.lower.x = MLI_MIN2(aabb.lower.x, c_root.x);
                aabb.lower.y = MLI_MIN2(aabb.lower.y, c_root.y);
                aabb.lower.z = MLI_MIN2(aabb.lower.z, c_root.z);

                aabb.upper.x = MLI_MAX2(aabb.upper.x, c_root.x);
                aabb.upper.y = MLI_MAX2(aabb.upper.y, c_root.y);
                aabb.upper.z = MLI_MAX2(aabb.upper.z, c_root.z);
        }
        return aabb;
}

struct mliAABB mliObject_aabb_in_local_frame(const struct mliObject *obj)
{
        struct mliHomTra unity;
        unity.translation = mliVec_init(0.0, 0.0, 0.0);
        unity.rotation = mliMat_unity();
        return mliObject_aabb(obj, unity);
}

/* mliObject_serialize */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliObject_fwrite(const struct mliObject *obj, FILE *f)
{
        struct mliMagicId magic;
        chk(mliMagicId_set(&magic, "mliObject"));
        chk_fwrite(&magic, sizeof(struct mliMagicId), 1u, f);

        chk_fwrite(&obj->num_vertices, sizeof(uint32_t), 1u, f);
        chk_fwrite(&obj->num_vertex_normals, sizeof(uint32_t), 1u, f);
        chk_fwrite(&obj->num_faces, sizeof(uint32_t), 1u, f);
        chk_fwrite(&obj->num_materials, sizeof(uint32_t), 1u, f);

        chk_fwrite(obj->vertices, sizeof(struct mliVec), obj->num_vertices, f);

        chk_fwrite(
                obj->vertex_normals,
                sizeof(struct mliVec),
                obj->num_vertex_normals,
                f);

        chk_fwrite(
                obj->faces_vertices, sizeof(struct mliFace), obj->num_faces, f);
        chk_fwrite(
                obj->faces_vertex_normals,
                sizeof(struct mliFace),
                obj->num_faces,
                f);
        chk_fwrite(obj->faces_materials, sizeof(uint16_t), obj->num_faces, f);
        chk_fwrite(
                obj->material_names,
                sizeof(struct mliName),
                obj->num_materials,
                f);

        return 1;
error:
        return 0;
}

int mliObject_malloc_fread(struct mliObject *obj, FILE *f)
{
        uint32_t num_vertices;
        uint32_t num_vertex_normals;
        uint32_t num_faces;
        uint32_t num_materials;
        struct mliMagicId magic;
        chk_fread(&magic, sizeof(struct mliMagicId), 1u, f);
        chk(mliMagicId_has_word(&magic, "mliObject"));
        mliMagicId_warn_version(&magic);

        chk_fread(&num_vertices, sizeof(uint32_t), 1u, f);
        chk_fread(&num_vertex_normals, sizeof(uint32_t), 1u, f);
        chk_fread(&num_faces, sizeof(uint32_t), 1u, f);
        chk_fread(&num_materials, sizeof(uint32_t), 1u, f);

        chk(mliObject_malloc(
                obj,
                num_vertices,
                num_vertex_normals,
                num_faces,
                num_materials));

        chk_fread(obj->vertices, sizeof(struct mliVec), obj->num_vertices, f);
        chk_fread(
                obj->vertex_normals,
                sizeof(struct mliVec),
                obj->num_vertex_normals,
                f);

        chk_fread(
                obj->faces_vertices, sizeof(struct mliFace), obj->num_faces, f);
        chk_fread(
                obj->faces_vertex_normals,
                sizeof(struct mliFace),
                obj->num_faces,
                f);
        chk_fread(obj->faces_materials, sizeof(uint16_t), obj->num_faces, f);
        chk_fread(
                obj->material_names,
                sizeof(struct mliName),
                obj->num_materials,
                f);

        chk_msg(mliObject_has_valid_faces(obj),
                "A face refers to a not existing vertex/vertex_normal.");
        return 1;
error:
        mliObject_free(obj);
        return 0;
}

/* mliObject_valid */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliObject_is_valid(const struct mliObject *obj)
{
        chk_msg(mliObject_has_valid_vertices(obj), "Bad vertex.");
        chk_msg(mliObject_has_valid_faces(obj), "Expected faces to be valid.");
        chk_msg(mliObject_has_valid_normals(obj, MLI_EPSILON),
                "Bad vertex-normal.");
        chk_msg(mliObject_has_valid_materials(obj), "Bad material.");
        return 1;
error:
        return 0;
}

int mliObject_has_valid_faces(const struct mliObject *obj)
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
error:
        fprintf(stderr, "In obj.faces[%u]\n", i);
        return 0;
}

int mliObject_has_valid_vertices(const struct mliObject *obj)
{
        uint32_t i = 0;
        for (i = 0; i < obj->num_vertices; i++) {
                chk_msg(!MLI_IS_NAN(obj->vertices[i].x), "X is 'nan'.");
                chk_msg(!MLI_IS_NAN(obj->vertices[i].y), "Y is 'nan'.");
                chk_msg(!MLI_IS_NAN(obj->vertices[i].z), "Z is 'nan'.");
        }
        return 1;
error:
        fprintf(stderr, "In obj.vertices[%u]\n", i);
        return 0;
}

int mliObject_has_valid_normals(
        const struct mliObject *obj,
        const double epsilon)
{
        uint32_t i = 0;
        for (i = 0; i < obj->num_vertex_normals; i++) {
                double norm;
                chk_msg(!MLI_IS_NAN(obj->vertex_normals[i].x), "X is 'nan'.");
                chk_msg(!MLI_IS_NAN(obj->vertex_normals[i].y), "Y is 'nan'.");
                chk_msg(!MLI_IS_NAN(obj->vertex_normals[i].z), "Z is 'nan'.");

                norm = mliVec_norm(obj->vertex_normals[i]);
                chk_msg(fabs(norm - 1.0) <= epsilon,
                        "Expected vertex_normals to be normalized.");
        }
        return 1;
error:
        fprintf(stderr, "In obj.vertex_normals[%u]\n", i);
        return 0;
}

int mliObject_has_valid_materials(const struct mliObject *obj)
{
        uint32_t i = 0;
        for (i = 0; i < obj->num_materials; i++) {
                chk_msg(mliName_valid(&obj->material_names[i]),
                        "Expected material_name to be '\\0' terminated.");
                chk_msg(strlen(obj->material_names[i].cstr) > 0,
                        "Expected strlen(material_name) > 0.");
        }
        return 1;
error:
        fprintf(stderr, "In obj.materials[%u]\n", i);
        return 0;
}

int mliObject_num_unused(
        const struct mliObject *obj,
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
        MLI_ARRAY_SET(v, 0, obj->num_vertices);

        chk_malloc(vn, uint64_t, obj->num_vertex_normals);
        MLI_ARRAY_SET(vn, 0, obj->num_vertex_normals);

        chk_malloc(mtl, uint64_t, obj->num_materials);
        MLI_ARRAY_SET(mtl, 0, obj->num_materials);

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
error:
        return 0;
}

/* mliObject_wavefront */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

#define MLI_WAVEFRONT_FACE_LINE_V 7
#define MLI_WAVEFRONT_FACE_LINE_V_VN 37
#define MLI_WAVEFRONT_FACE_LINE_V_VT_VN 25

#define MLI_WAVEFRONT_LINE_BUFF_LENGTH 64

int mliObject_is_face_line_toggle(const int state)
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

int mliStr_to_uint32(uint32_t *out, const struct mliStr *str)
{
        uint64_t u = 0;
        chk(mliStr_to_uint64(&u, str, 10));
        (*out) = (uint32_t)u;
        return 1;
error:
        return 0;
}

int mliObject_parse_face_line(
        const char *line,
        struct mliFace *faces_vertices,
        struct mliFace *faces_texture_points,
        struct mliFace *faces_vertex_normals,
        int *line_mode)
{
        struct mliFace *v = faces_vertices;
        struct mliFace *vt = faces_texture_points;
        struct mliFace *vn = faces_vertex_normals;
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

        struct mliStr wuff = mliStr_init();
        chk(mliStr_malloc(&wuff, MAX_NUM_CHARS));
        wuff.length = 0;

        while (state != final_state) {
                chk_msg(i <= MAX_NUM_CHARS, "Expected less chars in line.");
                c = line[i];

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

                if (mliObject_is_face_line_toggle(state)) {
                        wuff.cstr[wuff.length] = c;
                        wuff.length++;
                } else if (mliObject_is_face_line_toggle(old_state)) {
                        uint64_t r;
                        for (r = wuff.length; r < MAX_NUM_CHARS; r++) {
                                wuff.cstr[r] = '\0';
                        }

                        switch (old_state) {
                        /* MLI_WAVEFRONT_FACE_LINE_V                          */
                        case 3:
                                chk(mliStr_to_uint32(&v->a, &wuff));
                                break;
                        case 5:
                                chk(mliStr_to_uint32(&v->b, &wuff));
                                break;
                        case 7:
                                chk(mliStr_to_uint32(&v->c, &wuff));
                                break;

                        /* MLI_WAVEFRONT_FACE_LINE_V_VN                       */
                        /*    3                        v->a                   */
                        case 27:
                                chk(mliStr_to_uint32(&vn->a, &wuff));
                                break;
                        case 29:
                                chk(mliStr_to_uint32(&v->b, &wuff));
                                break;
                        case 32:
                                chk(mliStr_to_uint32(&vn->b, &wuff));
                                break;
                        case 34:
                                chk(mliStr_to_uint32(&v->c, &wuff));
                                break;
                        case 37:
                                chk(mliStr_to_uint32(&vn->c, &wuff));
                                break;

                        /* MLI_WAVEFRONT_FACE_LINE_V_VT_VN                    */
                        /*    3                        v->a                   */
                        case 11:
                                chk(mliStr_to_uint32(&vt->a, &wuff));
                                break;
                        case 13:
                                chk(mliStr_to_uint32(&vn->a, &wuff));
                                break;
                        case 15:
                                chk(mliStr_to_uint32(&v->b, &wuff));
                                break;
                        case 17:
                                chk(mliStr_to_uint32(&vt->b, &wuff));
                                break;
                        case 19:
                                chk(mliStr_to_uint32(&vn->b, &wuff));
                                break;
                        case 21:
                                chk(mliStr_to_uint32(&v->c, &wuff));
                                break;
                        case 23:
                                chk(mliStr_to_uint32(&vt->c, &wuff));
                                break;
                        case 25:
                                chk(mliStr_to_uint32(&vn->c, &wuff));
                                break;

                        default:
                                break;
                        }
                        wuff.length = 0;
                }

                old_state = state;
                i++;
        }
        mliStr_free(&wuff);
        return 1;
error:
        mliStr_free(&wuff);
        return 0;
}

int mliObject_is_vert_line_toggle(const int state)
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

int mliObject_parse_three_float_line(const char *line, struct mliVec *v)
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

        struct mliStr wuff = mliStr_init();
        chk(mliStr_malloc(&wuff, MAX_NUM_CHARS));
        wuff.length = 0;

        while (state != final_state) {
                chk_msg(i <= MAX_NUM_CHARS, "Expected less chars in line.");
                c = line[i];

                if (state == error_state) {
                        fprintf(stderr,
                                "[ERROR] Can not parse line '%s'\n",
                                line);
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

                if (mliObject_is_vert_line_toggle(state)) {
                        wuff.cstr[wuff.length] = c;
                        wuff.length++;
                } else if (mliObject_is_vert_line_toggle(old_state)) {
                        uint64_t r;
                        for (r = wuff.length; r < MAX_NUM_CHARS; r++) {
                                wuff.cstr[r] = '\0';
                        }

                        switch (old_state) {
                        case 2:
                                chk(mliStr_to_double(&v->x, &wuff));
                                break;
                        case 4:
                                chk(mliStr_to_double(&v->y, &wuff));
                                break;
                        case 6:
                                chk(mliStr_to_double(&v->z, &wuff));
                                break;
                        default:
                                break;
                        }
                        wuff.length = 0;
                }

                old_state = state;
                i++;
        }
        return 1;
error:
        return 0;
}

int mliObject_parse_face_vertices_and_normals(
        const char *line,
        struct mliFace *fv,
        struct mliFace *fvn)
{
        int line_mode = -1;
        struct mliFace tmp_fvt;

        chk_msg(mliObject_parse_face_line(line, fv, &tmp_fvt, fvn, &line_mode),
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
error:
        return 0;
}

int mliObject_malloc_from_wavefront(struct mliObject *obj, const char *str)
{
        uint64_t i = 0u;
        uint64_t p = 0u;
        uint64_t line_number = 0u;
        char line[2 * MLI_NAME_CAPACITY];
        uint64_t line_length = 0u;
        uint64_t mtl = 0u;
        const uint64_t debug_line_radius = 5u;

        /* init dyn */
        struct mliDynVec v = mliDynVec_init();
        struct mliDynVec vn = mliDynVec_init();

        struct mliDynFace fv = mliDynFace_init();
        struct mliDynFace fvn = mliDynFace_init();
        struct mliDynUint32 fm = mliDynUint32_init();

        struct mliDynMap material_names = mliDynMap_init();

        memset(line, '\0', sizeof(line));

        /* malloc dyn */
        chk(mliDynVec_malloc(&v, 0u));
        chk(mliDynVec_malloc(&vn, 0u));

        chk(mliDynFace_malloc(&fv, 0u));
        chk(mliDynFace_malloc(&fvn, 0u));
        chk(mliDynUint32_malloc(&fm, 0u));

        chk(mliDynMap_malloc(&material_names, 0u));

        /* parse wavefront into dyn */
        while (1) {
                line_number += 1;
                chk_msg(line_number < 1000 * 1000 * 1000,
                        "Expected less than 1e9 lines in wavefront-file. "
                        "Something went wrong.");

                line_length = mli_cstr_split(&str[p], '\n', line, sizeof(line));

                chk_msg(line_length < sizeof(line), "Line is too long.");

                if (line_length > 0) {
                        if (mli_cstr_starts_with(line, "vn ")) {
                                struct mliVec tmp_vn;
                                chk_msg(mliObject_parse_three_float_line(
                                                &line[2], &tmp_vn),
                                        "Can not parse vertex-normal-line.");
                                chk_msg(mliVec_dot(tmp_vn, tmp_vn) > 0.0,
                                        "vn can not be normalized.") tmp_vn =
                                        mliVec_normalized(tmp_vn);

                                chk(mliDynVec_push_back(&vn, tmp_vn));
                        } else if (mli_cstr_starts_with(line, "v ")) {
                                struct mliVec tmp_v;
                                chk_msg(mliObject_parse_three_float_line(
                                                &line[1], &tmp_v),
                                        "Can not parse vertex-line.");
                                chk(mliDynVec_push_back(&v, tmp_v));
                        } else if (mli_cstr_starts_with(line, "f ")) {
                                struct mliFace tmp_fv;
                                struct mliFace tmp_fvn;
                                chk_msg(material_names.size > 0,
                                        "Expected 'usemtl' before first "
                                        "face 'f'.");
                                chk_msg(mliObject_parse_face_vertices_and_normals(
                                                line, &tmp_fv, &tmp_fvn),
                                        "Failed to parse face-line.");
                                chk(mliDynFace_push_back(&fv, tmp_fv));
                                chk(mliDynFace_push_back(&fvn, tmp_fvn));
                                chk(mliDynUint32_push_back(&fm, mtl));
                        } else if (mli_cstr_starts_with(line, "usemtl ")) {
                                const char *mtl_key = &line[7];
                                if (!mliDynMap_has(&material_names, mtl_key)) {
                                        chk(mliDynMap_insert(
                                                &material_names, mtl_key, 0));
                                }
                                chk(mliDynMap_find(
                                        &material_names, mtl_key, &mtl));
                        }
                } /* line_length > 0 */

                if (str[p + line_length] == '\0') {
                        break;
                }
                p += line_length + 1;
        }

        /* copy dyn into static mliObject */
        chk_msg(fv.size == fvn.size,
                "Expected num. vertex-indices == num. vertex-normal-indices.");
        chk_msg(mliObject_malloc(
                        obj, v.size, vn.size, fv.size, material_names.size),
                "Failed to malloc mliObject from file.");

        MLI_NCPY(v.array, obj->vertices, v.size);
        MLI_NCPY(vn.array, obj->vertex_normals, vn.size);
        MLI_NCPY(fv.array, obj->faces_vertices, fv.size);
        MLI_NCPY(fvn.array, obj->faces_vertex_normals, fvn.size);
        MLI_NCPY(fm.array, obj->faces_materials, fm.size);
        for (i = 0; i < material_names.size; i++) {
                memcpy(obj->material_names[i].cstr,
                       material_names.array[i].key,
                       MLI_NAME_CAPACITY);
        }

        chk_msg(mliObject_is_valid(obj), "Expected object to be valid.");

        /* free dyn */
        mliDynVec_free(&v);
        mliDynVec_free(&vn);

        mliDynFace_free(&fv);
        mliDynFace_free(&fvn);
        mliDynUint32_free(&fm);

        mliDynMap_free(&material_names);

        return 1;
error:
        mliObject_free(obj);

        /* free dyn */
        mliDynVec_free(&v);
        mliDynVec_free(&vn);

        mliDynFace_free(&fv);
        mliDynFace_free(&fvn);
        mliDynUint32_free(&fm);

        mliDynMap_free(&material_names);

        mli_cstr_lines_fprint(stderr, str, line_number, debug_line_radius);

        return 0;
}

int mliObject_fprint_to_wavefront(struct mliIo *f, const struct mliObject *obj)
{
        uint32_t i, mtl, face;
        chk(mliIo_printf(f, "# vertices\n"));
        for (i = 0; i < obj->num_vertices; i++) {
                chk(mliIo_printf(
                        f,
                        "v %.6f %.6f %.6f\n",
                        obj->vertices[i].x,
                        obj->vertices[i].y,
                        obj->vertices[i].z));
        }

        chk(mliIo_printf(f, "# vertex normals\n"));
        for (i = 0; i < obj->num_vertex_normals; i++) {
                chk(mliIo_printf(
                        f,
                        "vn %.6f %.6f %.6f\n",
                        obj->vertex_normals[i].x,
                        obj->vertex_normals[i].y,
                        obj->vertex_normals[i].z));
        }

        chk(mliIo_printf(f, "# faces\n"));
        for (face = 0; face < obj->num_faces; face++) {
                if ((face == 0) || (mtl != obj->faces_materials[face])) {
                        mtl = obj->faces_materials[face];
                        chk(mliIo_printf(
                                f,
                                "usemtl %s\n",
                                obj->material_names[mtl].cstr));
                }

                chk(mliIo_printf(
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
error:
        return 0;
}

/* mliOcTree */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/*
 * mliOcTree
 * =========
 *
 * A cache-aware mliOcTree to accelerate the traversal.
 * mliOcTree is created from the dynamic, mliTmpOcTree.
 *
 * Structure of mliTmpOcTree
 * -------------------------
 *
 *                               mliTmpNode(0,-)
 *                                    |-objects[a,b,c,d,e,f,g,h]
 *                     _______________|________________
 *                    |                                |
 *                mliTmpNode(1,-)                  mliTmpNode(2,-)
 *                    |objects[a,b,c,d]                |objects[e,f,g,h]
 *             _______|_______                 ________|________
 *            |               |               |                 |
 *        mliTmpNode(3,0) mliTmpNode(4,1) mliTmpNode(5,2)   mliTmpNode(6,3)
 *            |objects[a,b,c] |objects[c,d]  |objects[e,f]     |objects[f,g,h]
 *
 *      mliTmpNode link to each other via pointer.
 *      Each mliTmpNode is allocated separately.
 *
 *      advantages:
 *              easy to grow and fill
 *      disadvantage:
 *              memory is scattered, difficult to serialize
 *
 * Structure of mliOcTree
 * ----------------------
 *
 *                object_links w.r.t. mliGeometry ->
 *          leafs +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
 *            |---|  a  |  b  |  c  |  c  |  d  |  e  |  f  |  f  |  g  |  h  |
 *            |   +-----+-----+-----+-----+-----+-----+-----+-----+-----+-----+
 * mliOcTree  |   0     1     2     3     4     5     6     7     8     9     10
 *  |         |   ^                 ^           ^           ^
 *  |         |
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
 *            |  1  |  0  |  2  |                  w.r.t. mliOcTree.nodes
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

struct mliLeafAddress mliLeafAddress_init(void)
{
        struct mliLeafAddress address;
        address.first_object_link = 0u;
        address.num_object_links = 0u;
        return address;
}

struct mliLeafArray mliLeafArray_init(void)
{
        struct mliLeafArray leafs;
        leafs.num_leafs = 0;
        leafs.adresses = NULL;
        leafs.num_object_links = 0u;
        leafs.object_links = NULL;
        return leafs;
}

void mliLeafArray_free(struct mliLeafArray *leafs)
{
        free(leafs->object_links);
        free(leafs->adresses);
        *leafs = mliLeafArray_init();
}

int mliLeafArray_malloc(
        struct mliLeafArray *leafs,
        const uint64_t num_leafs,
        const uint64_t num_object_links)
{
        uint64_t i;
        mliLeafArray_free(leafs);
        leafs->num_leafs = num_leafs;
        chk_malloc(leafs->adresses, struct mliLeafAddress, leafs->num_leafs);
        for (i = 0; i < leafs->num_leafs; i++) {
                leafs->adresses[i] = mliLeafAddress_init();
        }
        leafs->num_object_links = num_object_links;
        chk_malloc(leafs->object_links, uint32_t, leafs->num_object_links);
        for (i = 0; i < leafs->num_object_links; i++) {
                leafs->object_links[i] = 0u;
        }
        return 1;
error:
        return 0;
}

struct mliNode mliNode_init(void)
{
        uint64_t c = 0;
        struct mliNode node;
        for (c = 0; c < 8u; c++) {
                node.children[c] = 0u;
                node.types[c] = MLI_OCTREE_TYPE_NONE;
        }
        return node;
}

struct mliOcTree mliOcTree_init(void)
{
        struct mliOcTree tree;
        tree.cube.lower = mliVec_init(0., 0., 0.);
        tree.cube.edge_length = 0.;
        tree.num_nodes = 0u;
        tree.nodes = NULL;
        tree.leafs = mliLeafArray_init();
        tree.root_type = MLI_OCTREE_TYPE_NONE;
        return tree;
}

void mliOcTree_free(struct mliOcTree *tree)
{
        free(tree->nodes);
        mliLeafArray_free(&tree->leafs);
        *tree = mliOcTree_init();
}

int mliOcTree_malloc(
        struct mliOcTree *tree,
        const uint64_t num_nodes,
        const uint64_t num_leafs,
        const uint64_t num_object_links)
{
        uint64_t i;
        mliOcTree_free(tree);
        tree->num_nodes = num_nodes;
        chk_malloc(tree->nodes, struct mliNode, tree->num_nodes);
        for (i = 0; i < tree->num_nodes; i++)
                tree->nodes[i] = mliNode_init();
        chk(mliLeafArray_malloc(&tree->leafs, num_leafs, num_object_links));
        return 1;
error:
        return 0;
}

void mliOcTree_set_node(
        struct mliOcTree *tree,
        const struct mliTmpNode *dynnode)
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

void mliOcTree_set_leaf(
        struct mliOcTree *tree,
        const struct mliTmpNode *dynnode,
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

void mliOcTree_set_walk(
        struct mliOcTree *tree,
        const struct mliTmpNode *dynnode,
        uint64_t *object_link_size)
{
        uint64_t c;
        if (dynnode->node_index >= 0) {
                mliOcTree_set_node(tree, dynnode);
        } else if (dynnode->leaf_index >= 0) {
                mliOcTree_set_leaf(tree, dynnode, object_link_size);
        }

        for (c = 0; c < 8u; c++) {
                if (dynnode->children[c] != NULL) {
                        mliOcTree_set_walk(
                                tree, dynnode->children[c], object_link_size);
                }
        }
}

void mliOcTree_set(struct mliOcTree *tree, const struct mliTmpOcTree *dyntree)
{
        uint64_t object_link_size = 0u;
        tree->cube = dyntree->cube;
        if (mliTmpNode_num_children(&dyntree->root) > 0) {
                tree->root_type = MLI_OCTREE_TYPE_NODE;
        } else {
                tree->root_type = MLI_OCTREE_TYPE_LEAF;
        }

        mliOcTree_set_walk(tree, &dyntree->root, &object_link_size);
}

uint64_t mliOcTree_node_num_children(
        const struct mliOcTree *tree,
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

uint64_t mliOcTree_leaf_num_objects(
        const struct mliOcTree *tree,
        const uint64_t leaf)
{
        return tree->leafs.adresses[leaf].num_object_links;
}

uint32_t mliOcTree_leaf_object_link(
        const struct mliOcTree *tree,
        const uint64_t leaf,
        const uint64_t object_link)
{
        uint64_t i = tree->leafs.adresses[leaf].first_object_link + object_link;
        return tree->leafs.object_links[i];
}

int mliOcTree_equal_payload_walk(
        const struct mliOcTree *tree,
        const int32_t node_idx,
        const int32_t node_type,
        const struct mliTmpNode *tmp_node)
{
        if (node_type == MLI_OCTREE_TYPE_LEAF) {
                /* leaf */
                uint64_t leaf_idx;
                uint64_t obj;
                chk_msg(mliTmpNode_num_children(tmp_node) == 0,
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
                                chk_msg(mliOcTree_equal_payload_walk(
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
error:
        return 0;
}

int mliOcTree_equal_payload(
        const struct mliOcTree *tree,
        const struct mliTmpOcTree *tmp_octree)
{
        int32_t root_node_idx = 0;
        int32_t root_node_type = MLI_OCTREE_TYPE_NODE;
        chk_msg(mliCube_equal(tree->cube, tmp_octree->cube),
                "Cubes are not equal");
        chk_msg(mliOcTree_equal_payload_walk(
                        tree, root_node_idx, root_node_type, &tmp_octree->root),
                "Tree is not equal");
        return 1;
error:
        return 0;
}

void mliOcTree_print_walk(
        const struct mliOcTree *tree,
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
                                mliOcTree_print_walk(
                                        tree,
                                        child_node_idx,
                                        child_node_type,
                                        indent + 2,
                                        c);
                        }
                }
        }
}

void mliOcTree_print(const struct mliOcTree *tree)
{
        printf("__ mliOctree __\n");
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
                mliOcTree_print_walk(tree, root_idx, tree->root_type, 0u, 0u);
        }
}

int mliOcTree_malloc_from_object_wavefront(
        struct mliOcTree *octree,
        const struct mliObject *object)
{
        uint64_t num_nodes = 0;
        uint64_t num_leafs = 0;
        uint64_t num_object_links = 0;
        struct mliTmpOcTree tmp_octree = mliTmpOcTree_init();
        chk_msg(mliTmpOcTree_malloc_from_bundle(
                        &tmp_octree,
                        (const void *)object,
                        object->num_faces,
                        mliObject_face_in_local_frame_has_overlap_aabb_void,
                        mliObject_aabb_in_local_frame(object)),
                "Failed to create dynamic, and temporary TmpOcTree "
                "from mliObject");
        mliTmpNode_set_flat_index(&tmp_octree.root);
        mliTmpNode_num_nodes_leafs_objects(
                &tmp_octree.root, &num_nodes, &num_leafs, &num_object_links);

        chk_msg(mliOcTree_malloc(
                        octree, num_nodes, num_leafs, num_object_links),
                "Failed to allocate cache-aware octree from dynamic octree.");
        mliOcTree_set(octree, &tmp_octree);
        mliTmpOcTree_free(&tmp_octree);

        return 1;
error:
        return 0;
}

int mliOcTree_malloc_from_Geometry(
        struct mliOcTree *octree,
        const struct mliGeometryAndAccelerator *accgeo,
        const struct mliAABB outermost_aabb)
{
        uint64_t num_nodes = 0;
        uint64_t num_leafs = 0;
        uint64_t num_object_links = 0;
        struct mliTmpOcTree tmp_octree = mliTmpOcTree_init();
        chk_msg(mliTmpOcTree_malloc_from_bundle(
                        &tmp_octree,
                        (const void *)accgeo,
                        accgeo->geometry->num_robjects,
                        mliGeometry_robject_has_overlap_aabb_void,
                        outermost_aabb),
                "Failed to create dynamic, and temporary TmpOcTree "
                "from scenery(Geometry, Accelerator)");
        mliTmpNode_set_flat_index(&tmp_octree.root);
        mliTmpNode_num_nodes_leafs_objects(
                &tmp_octree.root, &num_nodes, &num_leafs, &num_object_links);

        chk_msg(mliOcTree_malloc(
                        octree, num_nodes, num_leafs, num_object_links),
                "Failed to allocate cache-aware octree from dynamic octree.");
        mliOcTree_set(octree, &tmp_octree);
        mliTmpOcTree_free(&tmp_octree);

        return 1;
error:
        return 0;
}

/* mliOcTree_equal */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliOcTree_equal(const struct mliOcTree *a, const struct mliOcTree *b)
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
error:
        return 0;
}

/* mliOcTree_serialize */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliOcTree_fwrite(const struct mliOcTree *octree, FILE *f)
{
        struct mliMagicId magic;

        /* magic identifier */
        chk(mliMagicId_set(&magic, "mliOcTree"));
        chk_fwrite(&magic, sizeof(struct mliMagicId), 1u, f);

        /* capacity */
        chk_fwrite(&octree->num_nodes, sizeof(uint64_t), 1u, f);
        chk_fwrite(&octree->leafs.num_leafs, sizeof(uint64_t), 1u, f);
        chk_fwrite(&octree->leafs.num_object_links, sizeof(uint64_t), 1u, f);

        /* nodes */
        chk_fwrite(octree->nodes, sizeof(struct mliNode), octree->num_nodes, f);

        /* leaf addresses */
        chk_fwrite(
                octree->leafs.adresses,
                sizeof(struct mliLeafAddress),
                octree->leafs.num_leafs,
                f);

        /* leaf links */
        chk_fwrite(
                octree->leafs.object_links,
                sizeof(uint32_t),
                octree->leafs.num_object_links,
                f);

        /* mliCube */
        chk_fwrite(&octree->cube.lower.x, sizeof(double), 1u, f);
        chk_fwrite(&octree->cube.lower.y, sizeof(double), 1u, f);
        chk_fwrite(&octree->cube.lower.z, sizeof(double), 1u, f);
        chk_fwrite(&octree->cube.edge_length, sizeof(double), 1u, f);

        /* root type */
        chk_fwrite(&octree->root_type, sizeof(uint8_t), 1u, f);

        return 1;
error:
        return 0;
}

int mliOcTree_malloc_fread(struct mliOcTree *octree, FILE *f)
{
        uint64_t num_nodes;
        uint64_t num_leafs;
        uint64_t num_object_links;
        struct mliMagicId magic;

        /* magic identifier */
        chk_fread(&magic, sizeof(struct mliMagicId), 1u, f);
        chk(mliMagicId_has_word(&magic, "mliOcTree"));
        mliMagicId_warn_version(&magic);

        /* capacity */
        chk_fread(&num_nodes, sizeof(uint64_t), 1u, f);
        chk_fread(&num_leafs, sizeof(uint64_t), 1u, f);
        chk_fread(&num_object_links, sizeof(uint64_t), 1u, f);

        chk_msg(mliOcTree_malloc(
                        octree, num_nodes, num_leafs, num_object_links),
                "Can not malloc octree from file.");

        chk_fread(octree->nodes, sizeof(struct mliNode), octree->num_nodes, f);
        chk_fread(
                octree->leafs.adresses,
                sizeof(struct mliLeafAddress),
                octree->leafs.num_leafs,
                f);
        chk_fread(
                octree->leafs.object_links,
                sizeof(uint32_t),
                octree->leafs.num_object_links,
                f);

        /* mliCube */
        chk_fread(&octree->cube.lower.x, sizeof(double), 1u, f);
        chk_fread(&octree->cube.lower.y, sizeof(double), 1u, f);
        chk_fread(&octree->cube.lower.z, sizeof(double), 1u, f);
        chk_fread(&octree->cube.edge_length, sizeof(double), 1u, f);

        /* root type */
        chk_fread(&octree->root_type, sizeof(uint8_t), 1u, f);

        return 1;
error:
        return 0;
}

int mliOcTree_write_to_path(const struct mliOcTree *octree, const char *path)
{
        FILE *f;
        f = fopen(path, "w");
        chk_msg(f != NULL, "Can not open octree-file for writing.");

        chk_msg(mliOcTree_fwrite(octree, f), "Can not write octree to file.");

        fclose(f);
        return 1;
error:
        if (f != NULL) {
                fclose(f);
        }
        return 0;
}

int mliOcTree_malloc_from_path(struct mliOcTree *octree, const char *path)
{
        FILE *f;
        f = fopen(path, "r");
        chk_msg(f != NULL, "Can not open octree-file for reading.");

        chk_msg(mliOcTree_malloc_fread(octree, f),
                "Can not read octree to file.");

        fclose(f);
        return 1;
error:
        if (f != NULL) {
                fclose(f);
        }
        return 0;
}

/* mliOcTree_valid */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliOcTree_valid(const struct mliOcTree *octree)
{
        uint32_t n, c;
        chk_msg(!MLI_IS_NAN(octree->cube.lower.x), "cube.lower.x is 'nan'.");
        chk_msg(!MLI_IS_NAN(octree->cube.lower.y), "cube.lower.y is 'nan'.");
        chk_msg(!MLI_IS_NAN(octree->cube.lower.z), "cube.lower.z is 'nan'.");
        chk_msg(!MLI_IS_NAN(octree->cube.edge_length),
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
error:
        return 0;
}

int mliOcTree_valid_wrt_links(
        const struct mliOcTree *octree,
        const uint32_t num_links)
{
        uint32_t n;
        for (n = 0u; n < octree->leafs.num_object_links; n++) {
                chk_msg(octree->leafs.object_links[n] < num_links,
                        "Expected object_links[n] <  num_links.");
        }
        return 1;
error:
        return 0;
}

/* mliOctOverlaps */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/* mliPhoton */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/* mliPhotonInteraction */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_photoninteraction_type_to_string(const int32_t type, char *s)
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
error:
        return 0;
}

int mli_time_of_flight(
        const struct mliMaterials *materials,
        const struct mliPhotonInteraction *phisec,
        const double wavelength,
        double *time_of_flight)
{
        double refractive_index;
        chk_msg(mliFunc_evaluate(
                        &materials->media[phisec->medium_coming_from]
                                 .refraction,
                        wavelength,
                        &refractive_index),
                "Failed to eval. refraction for wavelength.");

        (*time_of_flight) = (refractive_index * phisec->distance_of_ray) /
                            MLI_VACUUM_SPPED_OF_LIGHT;
        return 1;
error:
        return 0;
}

/* mliPinHoleCamera */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliPinHoleCamera mliPinHoleCamera_init(
        const double field_of_view,
        const struct mliImage *image,
        const double row_over_column_pixel_ratio)
{
        struct mliPinHoleCamera sensor;

        assert(field_of_view > 0.);
        assert(field_of_view < mli_deg2rad(180.));

        sensor.optical_axis = mliVec_init(0.0, 0.0, 1.0);
        sensor.col_axis = mliVec_init(1.0, 0.0, 0.0);
        sensor.row_axis = mliVec_init(0.0, row_over_column_pixel_ratio, 0.0);

        sensor.distance_to_principal_point =
                ((0.5 * image->num_cols) / tan(0.5 * field_of_view));
        sensor.principal_point = mliVec_multiply(
                sensor.optical_axis, sensor.distance_to_principal_point);
        return sensor;
}

struct mliRay mliPinHoleCamera_ray_at_row_col(
        const struct mliPinHoleCamera *camera,
        const struct mliImage *image,
        const uint32_t row,
        const uint32_t col)
{
        struct mliVec pin_hole_position = mliVec_init(0.0, 0.0, 0.0);
        int row_idx_on_sensor = row - image->num_rows / 2;
        int col_idx_on_sensor = col - image->num_cols / 2;
        struct mliVec s_row =
                mliVec_multiply(camera->row_axis, row_idx_on_sensor);
        struct mliVec s_col =
                mliVec_multiply(camera->col_axis, col_idx_on_sensor);
        struct mliVec sensor_intersection = camera->principal_point;
        sensor_intersection = mliVec_add(sensor_intersection, s_row);
        sensor_intersection = mliVec_add(sensor_intersection, s_col);
        return mliRay_set(pin_hole_position, sensor_intersection);
}

void mliPinHoleCamera_render_image(
        struct mliPinHoleCamera camera,
        const struct mliHomTraComp camera2root_comp,
        const struct mliScenery *scenery,
        struct mliImage *image,
        const struct mliTracerConfig *tracer_config,
        struct mliPrng *prng)
{
        struct mliPixelWalk walk =
                mliPixelWalk_set(image->num_cols, image->num_rows, 16u);
        struct mliHomTra camera2root = mliHomTra_from_compact(camera2root_comp);
        uint32_t i;
        const uint32_t num_pixel = image->num_rows * image->num_cols;

        for (i = 0; i < num_pixel; i++) {
                struct mliPixel px = mliPixelWalk_get(&walk);
                struct mliRay ray_wrt_camera = mliPinHoleCamera_ray_at_row_col(
                        &camera, image, px.row, px.col);

                struct mliRay ray_wrt_root =
                        mliHomTra_ray(&camera2root, ray_wrt_camera);

                struct mliColor color =
                        mli_trace(scenery, ray_wrt_root, tracer_config, prng);
                mliImage_set(image, px.col, px.row, color);
                mliPixelWalk_walk(&walk);
        }
}

void mliPinHoleCamera_render_image_with_view(
        const struct mliView view,
        const struct mliScenery *scenery,
        struct mliImage *image,
        const double row_over_column_pixel_ratio,
        const struct mliTracerConfig *tracer_config,
        struct mliPrng *prng)
{
        struct mliPinHoleCamera camera = mliPinHoleCamera_init(
                view.field_of_view, image, row_over_column_pixel_ratio);
        struct mliHomTraComp camera2root_comp = mliView_to_HomTraComp(view);
        mliPinHoleCamera_render_image(
                camera, camera2root_comp, scenery, image, tracer_config, prng);
}

/* mliPixelWalk */
/* ------------ */

/* Copyright 2020-2021 Sebastian Achim Mueller */

struct mliPixelWalk mliPixelWalk_set(
        const uint32_t num_cols,
        const uint32_t num_rows,
        const uint32_t chunk_size)
{
        struct mliPixelWalk walk;
        uint32_t full_row, full_col, rest_row, rest_col;

        walk.num_cols = num_cols;
        walk.num_rows = num_rows;
        walk.chunk_size = chunk_size;

        walk.i = 0;

        full_row = num_rows / chunk_size;
        full_col = num_cols / chunk_size;
        rest_row = num_rows % chunk_size;
        rest_col = num_cols % chunk_size;

        walk.num_chunks_row = rest_row ? full_row + 1 : full_row;
        walk.num_chunks_col = rest_col ? full_col + 1 : full_col;

        walk.chunk_row = 0u;
        walk.sub_row = 0u;

        walk.chunk_col = 0u;
        walk.sub_col = 0u;

        return walk;
}

struct mliPixel mliPixelWalk_get(const struct mliPixelWalk *walk)
{
        struct mliPixel px;
        px.row = walk->chunk_row * walk->chunk_size + walk->sub_row;
        px.col = walk->chunk_col * walk->chunk_size + walk->sub_col;
        return px;
}

void mliPixelWalk_walk(struct mliPixelWalk *walk)
{
        struct mliPixel px;

        walk->sub_row += 1u;
        px = mliPixelWalk_get(walk);
        if (walk->sub_row < walk->chunk_size && px.row < walk->num_rows) {
                return;
        }

        walk->sub_row = 0u;
        walk->sub_col += 1u;
        px = mliPixelWalk_get(walk);
        if (walk->sub_col < walk->chunk_size && px.col < walk->num_cols) {
                return;
        }

        walk->sub_col = 0u;
        walk->chunk_row += 1u;
        if (walk->chunk_row < walk->num_chunks_row) {
                return;
        }

        walk->chunk_row = 0u;
        walk->chunk_col += 1u;
        if (walk->chunk_col < walk->num_chunks_col) {
                return;
        }

        walk->chunk_col = 0;
        return;
}

/* mliPixels */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliPixels mliPixels_init(void)
{
        struct mliPixels pix;
        pix.num_pixels_to_do = 0u;
        pix.num_pixels = 0u;
        pix.pixels = NULL;
        return pix;
}

void mliPixels_free(struct mliPixels *pix)
{
        free(pix->pixels);
        *pix = mliPixels_init();
}

int mliPixels_malloc(struct mliPixels *pix, const uint32_t num_pixels)
{
        mliPixels_free(pix);
        pix->num_pixels = num_pixels;
        chk_malloc(pix->pixels, struct mliPixel, pix->num_pixels);
        return 1;
error:
        mliPixels_free(pix);
        return 0;
}

/* mliQuaternion */
/* ------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliQuaternion mliQuaternion_set(
        const double w,
        const double x,
        const double y,
        const double z)
{
        struct mliQuaternion o;
        o.w = w;
        o.x = x;
        o.y = y;
        o.z = z;
        return o;
}

struct mliQuaternion mliQuaternion_set_unit_xyz(
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
        return mliQuaternion_set(w, x, y, z);
}

int mliQuaternion_equal(
        const struct mliQuaternion a,
        const struct mliQuaternion b)
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

int mliQuaternion_equal_margin(
        const struct mliQuaternion a,
        const struct mliQuaternion b,
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

struct mliQuaternion mliQuaternion_complex_conjugate(
        const struct mliQuaternion q)
{
        struct mliQuaternion c;
        c.w = q.w;
        c.x = -q.x;
        c.y = -q.y;
        c.z = -q.z;
        return c;
}

struct mliQuaternion mliQuaternion_product(
        const struct mliQuaternion p,
        const struct mliQuaternion q)
{
        struct mliQuaternion pq;
        const struct mliVec P = mliVec_init(p.x, p.y, p.z);
        const struct mliVec Q = mliVec_init(q.x, q.y, q.z);
        const struct mliVec P_cross_Q = mliVec_cross(P, Q);
        pq.w = p.w * q.w - mliVec_dot(P, Q);
        pq.x = p.w * Q.x + q.w * P.x + P_cross_Q.x;
        pq.y = p.w * Q.y + q.w * P.y + P_cross_Q.y;
        pq.z = p.w * Q.z + q.w * P.z + P_cross_Q.z;
        return pq;
}

double mliQuaternion_product_complex_conjugate(const struct mliQuaternion p)
{
        return p.w * p.w + p.x * p.x + p.y * p.y + p.z * p.z;
}

double mliQuaternion_norm(const struct mliQuaternion q)
{
        return sqrt(mliQuaternion_product_complex_conjugate(q));
}

struct mliQuaternion mliQuaternion_set_rotaxis_and_angle(
        const struct mliVec rot_axis,
        const double angle)
{
        const struct mliVec normed_rot_axis = mliVec_normalized(rot_axis);
        struct mliQuaternion quat;
        const double angle_half = .5 * angle;
        const double sin_angle_half = sin(angle_half);
        quat.w = cos(angle_half);
        quat.x = normed_rot_axis.x * sin_angle_half;
        quat.y = normed_rot_axis.y * sin_angle_half;
        quat.z = normed_rot_axis.z * sin_angle_half;
        return quat;
}

struct mliMat mliQuaternion_to_matrix(const struct mliQuaternion quat)
{
        struct mliMat o;
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

struct mliQuaternion mliQuaternion_set_tait_bryan(
        const double rx,
        const double ry,
        const double rz)
{
        const struct mliQuaternion qz =
                mliQuaternion_set_rotaxis_and_angle(mliVec_init(0, 0, 1), -rz);
        const struct mliQuaternion qy =
                mliQuaternion_set_rotaxis_and_angle(mliVec_init(0, 1, 0), -ry);
        const struct mliQuaternion qx =
                mliQuaternion_set_rotaxis_and_angle(mliVec_init(1, 0, 0), -rx);
        const struct mliQuaternion qz_qy = mliQuaternion_product(qz, qy);
        return mliQuaternion_product(qz_qy, qx);
}

/* mliQuaternion_json */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliQuaternion_tait_bryan_from_json(
        struct mliQuaternion *quat,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t token_xyz;
        struct mliVec xyz_deg;
        chk_msg(mliJson_token_by_key(json, token, "xyz_deg", &token_xyz),
                "Expected tait_bryan to have key 'xyz_deg'.");
        chk_msg(mliVec_from_json_token(&xyz_deg, json, token_xyz + 1),
                "Failed to parse tait_bryan's 'xyz_deg' from json.");
        *quat = mliQuaternion_set_tait_bryan(
                mli_deg2rad(xyz_deg.x),
                mli_deg2rad(xyz_deg.y),
                mli_deg2rad(xyz_deg.z));
        return 1;
error:
        return 0;
}

int mliQuaternion_axis_angle_from_json(
        struct mliQuaternion *quat,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t token_axis, token_angle;
        double angle_deg;
        struct mliVec axis;
        chk_msg(mliJson_token_by_key(json, token, "axis", &token_axis),
                "Expected axis_angle to have key 'axis'.");
        chk_msg(mliVec_from_json_token(&axis, json, token_axis + 1),
                "Failed to parse axis_angle's 'axis' from json.");
        chk_msg(mliJson_token_by_key(json, token, "angle_deg", &token_angle),
                "Expected axis_angle to have key 'angle_deg'.");
        chk_msg(mliJson_double_by_token(json, token_angle + 1, &angle_deg),
                "Failed to parse axis_angle's 'angle_deg' from json.");
        *quat = mliQuaternion_set_rotaxis_and_angle(
                axis, mli_deg2rad(angle_deg));
        return 1;
error:
        return 0;
}

int mliQuaternion_quaternion_from_json(
        struct mliQuaternion *quat,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t token_xyz;
        struct mliVec q;
        chk_msg(mliJson_token_by_key(json, token, "xyz", &token_xyz),
                "Expected quaternion to have key 'xyz'.");
        chk_msg(mliVec_from_json_token(&q, json, token_xyz + 1),
                "Failed to parse quaternion's 'xyz' from json.");
        *quat = mliQuaternion_set_unit_xyz(q.x, q.y, q.z);
        chk_msg(fabs(mliQuaternion_norm(*quat) - 1.) < 1e-6,
                "Expected norm(quaternion) < 1e-6. Expected unit-quaternion.");
        return 1;
error:
        return 0;
}

int mliQuaternion_from_json(
        struct mliQuaternion *quat,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t token_repr = 0u;
        uint64_t token_repr_value = 0u;
        chk_msg(mliJson_token_by_key(json, token, "repr", &token_repr),
                "Expected 'rot' to have key 'repr'.");
        token_repr_value = token_repr + 1;

        if (mliJson_cstrcmp(json, token_repr_value, "tait_bryan")) {
                chk_msg(mliQuaternion_tait_bryan_from_json(quat, json, token),
                        "Failed to parse tait_bryan rotation.");
        } else if (mliJson_cstrcmp(json, token_repr_value, "axis_angle")) {
                chk_msg(mliQuaternion_axis_angle_from_json(quat, json, token),
                        "Failed to parse axis_angle rotation.");
        } else if (mliJson_cstrcmp(json, token_repr_value, "quaternion")) {
                chk_msg(mliQuaternion_quaternion_from_json(quat, json, token),
                        "Failed to parse quaternion rotation.");
        } else {
                chk_bad("Unknown representation 'repr' in rotation.");
        }
        return 1;
error:
        return 0;
}

/* mliRay */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliRay mliRay_set(
        const struct mliVec support,
        const struct mliVec direction)
{
        struct mliRay ray;
        ray.support = support;
        ray.direction = mliVec_multiply(direction, 1. / mliVec_norm(direction));
        return ray;
}

struct mliVec mliRay_at(const struct mliRay *ray, const double t)
{
        struct mliVec out;
        out.x = ray->support.x + t * ray->direction.x;
        out.y = ray->support.y + t * ray->direction.y;
        out.z = ray->support.z + t * ray->direction.z;
        return out;
}

int mliRay_sphere_intersection(
        const struct mliVec support,
        const struct mliVec direction,
        const double radius,
        double *minus_solution,
        double *plus_solution)
{
        const double sup_times_dir = mliVec_dot(support, direction);
        const double dir_times_dir = mliVec_dot(direction, direction);
        const double sup_times_sup = mliVec_dot(support, support);
        const double radius_square = radius * radius;

        const double p = 2.0 * (sup_times_dir / dir_times_dir);
        const double q = sup_times_sup / dir_times_dir - radius_square;

        return mli_quadratic_equation(p, q, minus_solution, plus_solution);
}

/* mliRay_AABB */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliRay_has_overlap_aabb(
        const struct mliRay ray,
        const struct mliAABB aabb,
        double *ray_parameter)
{
        const double frac_x = 1. / ray.direction.x;
        const double frac_y = 1. / ray.direction.y;
        const double frac_z = 1. / ray.direction.z;

        const double t1 = (aabb.lower.x - ray.support.x) * frac_x;
        const double t2 = (aabb.upper.x - ray.support.x) * frac_x;
        const double t3 = (aabb.lower.y - ray.support.y) * frac_y;
        const double t4 = (aabb.upper.y - ray.support.y) * frac_y;
        const double t5 = (aabb.lower.z - ray.support.z) * frac_z;
        const double t6 = (aabb.upper.z - ray.support.z) * frac_z;

        const double tmin = MLI_MAX2(
                MLI_MAX2(MLI_MIN2(t1, t2), MLI_MIN2(t3, t4)), MLI_MIN2(t5, t6));
        const double tmax = MLI_MIN2(
                MLI_MIN2(MLI_MAX2(t1, t2), MLI_MAX2(t3, t4)), MLI_MAX2(t5, t6));
        if (tmax < 0) {
                (*ray_parameter) = tmax;
                return 0;
        }

        /* if tmin > tmax, ray doesn't intersect AABB */
        if (tmin > tmax) {
                (*ray_parameter) = tmax;
                return 0;
        }

        (*ray_parameter) = tmin;
        return 1;
}

/* mliRenderConfig */
/* --------------- */

/* Copyright 2018-2021 Sebastian Achim Mueller */

struct mliRenderConfig mliRenderConfig_init(void)
{
        struct mliRenderConfig c;
        c.camera = mliApertureCamera_init();
        c.camera_to_root.translation = mliVec_init(0.0, 0.0, 0.0);
        c.camera_to_root.rotation = mliQuaternion_set_tait_bryan(0.0, 0.0, 0.0);
        c.tracer = mliTracerConfig_init();
        c.num_pixel_x = 64;
        c.num_pixel_y = 48;
        c.random_seed = 0;
        return c;
}

int mliRenderConfig_camera_from_json(
        struct mliRenderConfig *cc,
        const struct mliJson *json,
        const uint64_t tkn)
{
        chk(mliFrame_pos_rot_from_json_token(&cc->camera_to_root, json, tkn));
        chk(mliJson_double_by_key(
                json,
                tkn,
                &cc->camera.image_sensor_width_x,
                "image_sensor_width_x"));
        chk(mliJson_double_by_key(
                json,
                tkn,
                &cc->camera.image_sensor_width_y,
                "image_sensor_width_y"));
        chk(mliJson_double_by_key(
                json, tkn, &cc->camera.focal_length, "focal_length"));
        chk(mliJson_double_by_key(
                json, tkn, &cc->camera.aperture_radius, "aperture_radius"));
        chk(mliJson_double_by_key(
                json,
                tkn,
                &cc->camera.image_sensor_distance,
                "image_sensor_distance"));

        chk_msg(cc->camera.image_sensor_width_x > 0.0,
                "Expected camera.image_sensor_width_x > 0.");
        chk_msg(cc->camera.image_sensor_width_y > 0.0,
                "Expected camera.image_sensor_width_y > 0.");
        chk_msg(cc->camera.focal_length > 0.0,
                "Expected camera.focal_length > 0.");
        chk_msg(cc->camera.aperture_radius > 0.0,
                "Expected camera.aperture_radius > 0.");
        chk_msg(cc->camera.image_sensor_distance > 0.0,
                "Expected camera.image_sensor_distance > 0.");
        chk_msg(cc->camera.image_sensor_distance >= cc->camera.focal_length,
                "Expected "
                "camera.image_sensor_distance >= camera.focal_length.");
        return 1;
error:
        return 0;
}

int mliRenderConfig_image_from_json(
        struct mliRenderConfig *cc,
        const struct mliJson *json,
        const uint64_t tkn)
{
        chk(mliJson_uint64_by_key(json, tkn, &cc->num_pixel_x, "num_pixel_x"));
        chk(mliJson_uint64_by_key(json, tkn, &cc->num_pixel_y, "num_pixel_y"));
        chk_msg(cc->num_pixel_x > 0, "Expected image.num_pixel_x > 0.");
        chk_msg(cc->num_pixel_y > 0, "Expected image.num_pixel_y > 0.");
        return 1;
error:
        return 0;
}

int mliRenderConfig_from_json(
        struct mliRenderConfig *cc,
        const struct mliJson *json,
        const uint64_t tkn)
{
        uint64_t camtkn;
        uint64_t imgtkn;
        uint64_t tratkn;
        chk(mliJson_uint64_by_key(json, tkn, &cc->random_seed, "random_seed"));
        chk(mliJson_token_by_key_eprint(json, tkn, "camera", &camtkn));
        chk(mliJson_token_by_key_eprint(json, tkn, "image", &imgtkn));
        chk(mliJson_token_by_key_eprint(json, tkn, "tracer", &tratkn));

        chk(mliRenderConfig_camera_from_json(cc, json, camtkn + 1));
        chk(mliRenderConfig_image_from_json(cc, json, imgtkn + 1));
        chk(mliTracerConfig_from_json_token(&cc->tracer, json, tratkn + 1));

        return 1;
error:
        return 0;
}

/* mliScenery */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliScenery mliScenery_init(void)
{
        struct mliScenery scenery;
        scenery.geometry = mliGeometry_init();
        scenery.accelerator = mliAccelerator_init();
        scenery.materials = mliMaterials_init();
        scenery.geomap = mliGeometryToMaterialMap_init();
        return scenery;
}

void mliScenery_free(struct mliScenery *scenery)
{
        mliGeometry_free(&scenery->geometry);
        mliAccelerator_free(&scenery->accelerator);
        mliMaterials_free(&scenery->materials);
        mliGeometryToMaterialMap_free(&scenery->geomap);
}

void mliScenery_info_fprint(FILE *f, const struct mliScenery *scenery)
{
        mliGeometry_info_fprint(f, &scenery->geometry);
        fprintf(f, "\n");
        mliAccelerator_info_fprint(f, &scenery->accelerator);
        fprintf(f, "\n");
        mliMaterials_info_fprint(f, &scenery->materials);
        fprintf(f, "\n");
        mliGeometryToMaterialMap_info_fprint(f, &scenery->geomap);
        fprintf(f, "\n");
}

/* mliScenery_equal */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliScenery_equal(const struct mliScenery *a, const struct mliScenery *b)
{
        chk_msg(mliGeometry_equal(&a->geometry, &b->geometry),
                "Expected geometry to be valid.");
        chk_msg(mliMaterials_equal(&a->materials, &b->materials),
                "Expected materials to be valid.");
        chk_msg(mliAccelerator_equal(&a->accelerator, &b->accelerator),
                "Expected accelerator to be valid.");
        chk_msg(mliGeometryToMaterialMap_equal(&a->geomap, &b->geomap),
                "Expected geomap to be valid.");
        return 1;
error:
        return 0;
}

/* mliScenery_minimal_object */
/* ------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliColor mli_random_color(struct mliPrng *prng)
{
        struct mliRandomUniformRange uniform_8bit_range;
        uniform_8bit_range.start = 0.0;
        uniform_8bit_range.range = 255.0;
        return mliColor_set(
                mli_random_draw_uniform(uniform_8bit_range, prng),
                mli_random_draw_uniform(uniform_8bit_range, prng),
                mli_random_draw_uniform(uniform_8bit_range, prng));
}

int mliScenery_malloc_minimal_from_wavefront(
        struct mliScenery *scenery,
        const char *path)
{
        uint32_t i, total_num_boundary_layers;
        struct mliPrng prng = mliPrng_init_MT19937(1u);
        struct mliIo str = mliIo_init();
        struct mliMaterialsCapacity mtlcap = mliMaterialsCapacity_init();

        mliScenery_free(scenery);

        chk_msg(mliGeometry_malloc(&scenery->geometry, 1u, 1u),
                "Failed to malloc geometry.");

        /* set object */
        chk_msg(mliIo_malloc_from_path(&str, path), "Failed to read file.");
        chk_msg(mli_cstr_assert_only_NUL_LF_TAB_controls((char *)str.cstr),
                "Expected object-wavefront file to be free of "
                "control characters, except [NUL, TAB, LF].");
        chk_msg(mliObject_malloc_from_wavefront(
                        &scenery->geometry.objects[0], (char *)str.cstr),
                "Failed to malloc wavefront-object from string.");
        mliIo_free(&str);
        sprintf(scenery->geometry.object_names[0].cstr, "default-object");

        /* set reference */
        scenery->geometry.robjects[0] = 0u;
        scenery->geometry.robject_ids[0] = 0u;
        scenery->geometry.robject2root[0] = mliHomTraComp_set(
                mliVec_init(0.0, 0.0, 0.0),
                mliQuaternion_set_tait_bryan(0.0, 0.0, 0.0));

        /* materials */
        total_num_boundary_layers = scenery->geometry.objects[0].num_materials;

        mtlcap.num_media = 1u;
        mtlcap.num_boundary_layers = total_num_boundary_layers;
        mtlcap.num_surfaces = total_num_boundary_layers;

        chk_msg(mliMaterials_malloc(&scenery->materials, mtlcap),
                "Failed to malloc materials.");

        sprintf(scenery->materials.medium_names[0].cstr, "vacuum");

        chk(mliFunc_malloc(&scenery->materials.media[0].refraction, 2));
        scenery->materials.media[0].refraction.x[0] = 200e-9;
        scenery->materials.media[0].refraction.x[1] = 1200e-9;
        scenery->materials.media[0].refraction.y[0] = 1.0;
        scenery->materials.media[0].refraction.y[1] = 1.0;

        chk(mliFunc_malloc(&scenery->materials.media[0].absorbtion, 2));
        scenery->materials.media[0].absorbtion.x[0] = 200e-9;
        scenery->materials.media[0].absorbtion.x[1] = 1200e-9;
        scenery->materials.media[0].absorbtion.y[0] = 0.0;
        scenery->materials.media[0].absorbtion.y[1] = 0.0;

        for (i = 0u; i < total_num_boundary_layers; i++) {
                scenery->materials.surfaces[i].material = MLI_MATERIAL_PHONG;
                scenery->materials.surfaces[i].color = mli_random_color(&prng);

                chk(mliFunc_malloc(
                        &scenery->materials.surfaces[i].specular_reflection,
                        2));
                scenery->materials.surfaces[i].specular_reflection.x[0] =
                        200e-9;
                scenery->materials.surfaces[i].specular_reflection.x[1] =
                        1200e-9;
                scenery->materials.surfaces[i].specular_reflection.y[0] = 0.0;
                scenery->materials.surfaces[i].specular_reflection.y[1] = 0.0;

                chk(mliFunc_malloc(
                        &scenery->materials.surfaces[i].diffuse_reflection, 2));
                scenery->materials.surfaces[i].diffuse_reflection.x[0] = 200e-9;
                scenery->materials.surfaces[i].diffuse_reflection.x[1] =
                        1200e-9;
                scenery->materials.surfaces[i].diffuse_reflection.y[0] = 1.0;
                scenery->materials.surfaces[i].diffuse_reflection.y[1] = 1.0;

                sprintf(scenery->materials.surface_names[i].cstr,
                        "surface_%06u",
                        i);

                scenery->materials.boundary_layers[i].inner.medium = 0u;
                scenery->materials.boundary_layers[i].outer.medium = 0u;
                scenery->materials.boundary_layers[i].inner.surface = i;
                scenery->materials.boundary_layers[i].outer.surface = i;
                sprintf(scenery->materials.boundary_layer_names[i].cstr,
                        "boundary_layer_%06u",
                        i);
        }

        chk_msg(mliGeometryToMaterialMap_malloc(
                        &scenery->geomap,
                        scenery->geometry.num_robjects,
                        total_num_boundary_layers),
                "Failed to malloc geometry to materials map.");

        /* set map */
        for (i = 0u; i < total_num_boundary_layers; i++) {
                mliGeometryToMaterialMap_set(&scenery->geomap, 0u, i, i);
        }

        chk_msg(mliAccelerator_malloc_from_Geometry(
                        &scenery->accelerator, &scenery->geometry),
                "Failed to malloc accelerator from geometry.");

        chk_msg(mliScenery_valid(scenery), "Expected scenery to be valid.");
        return 1;
error:
        return 0;
}

/* mliScenery_serialize */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliScenery_fwrite(const struct mliScenery *scenery, FILE *f)
{
        struct mliMagicId magic;
        chk(mliMagicId_set(&magic, "mliScenery"));
        chk_fwrite(&magic, sizeof(struct mliMagicId), 1u, f);

        chk(mliGeometry_fwrite(&scenery->geometry, f));
        chk(mliAccelerator_fwrite(&scenery->accelerator, f));
        chk(mliMaterials_fwrite(&scenery->materials, f));
        chk(mliGeometryToMaterialMap_fwrite(&scenery->geomap, f));
        return 1;
error:
        return 0;
}

int mliScenery_malloc_fread(struct mliScenery *scenery, FILE *f)
{
        struct mliMagicId magic;

        mliScenery_free(scenery);

        chk_fread(&magic, sizeof(struct mliMagicId), 1u, f);
        chk(mliMagicId_has_word(&magic, "mliScenery"));
        mliMagicId_warn_version(&magic);

        chk(mliGeometry_malloc_fread(&scenery->geometry, f));
        chk(mliAccelerator_malloc_fread(&scenery->accelerator, f));
        chk(mliMaterials_malloc_fread(&scenery->materials, f));
        chk(mliGeometryToMaterialMap_malloc_fread(&scenery->geomap, f));
        return 1;
error:
        mliScenery_free(scenery);
        return 0;
}

int mliScenery_malloc_from_path(struct mliScenery *scenery, const char *path)
{
        FILE *f;
        f = fopen(path, "r");
        chk_msg(f != NULL, "Can not open file for reading.");
        chk_msg(mliScenery_malloc_fread(scenery, f), "Can not read file.");
        fclose(f);
        return 1;
error:
        if (f != NULL) {
                fclose(f);
        }
        mliScenery_free(scenery);
        return 0;
}

int mliScenery_write_to_path(const struct mliScenery *scenery, const char *path)
{
        FILE *f;
        f = fopen(path, "w");
        chk_msg(f != NULL, "Can not open file for writing.");
        chk_msg(mliScenery_fwrite(scenery, f), "Failed to write to file.");
        fclose(f);
        return 1;
error:
        if (f != NULL) {
                fclose(f);
        }
        return 0;
}

/* mliScenery_tar */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliScenery_malloc_fread_tar(struct mliScenery *scenery, FILE *f)
{
        struct mliArchive archive = mliArchive_init();
        chk_msg(mliArchive_malloc_fread(&archive, f),
                "Can't read archive from file.");
        chk_msg(mliScenery_malloc_from_Archive(scenery, &archive),
                "Can't malloc Scenery from Archive.");
        mliArchive_free(&archive);
        return 1;
error:
        return 0;
}

int mliScenery_malloc_from_path_tar(
        struct mliScenery *scenery,
        const char *path)
{
        FILE *f = fopen(path, "rb");
        chk_msgf(f != NULL, ("Can't open path '%s'.", path))
                chk_msg(mliScenery_malloc_fread_tar(scenery, f),
                        "Can't fread Scenery from file.");
        fclose(f);
        return 1;
error:
        return 0;
}

int mliScenery_malloc_from_Archive(
        struct mliScenery *scenery,
        const struct mliArchive *archive)
{
        uint64_t num_robjects = 0u;
        uint64_t num_objects = 0u;
        uint64_t total_num_boundary_layers = 0u;

        struct mliNameMap material_names = mliNameMap_init();
        struct mliDynMap object_names = mliDynMap_init();
        struct mliFrame root = mliFrame_init();

        chk_msg(mliMaterials_malloc_form_archive(
                        &scenery->materials, &material_names, archive),
                "Failed to malloc materials.");

        chk(mliDynMap_malloc(&object_names, 0u));

        num_objects = mliArchive_num_filename_prefix_sufix(
                archive, "geometry/objects/", ".obj");

        chk_msg(mliGeometry_malloc_objects(&scenery->geometry, num_objects),
                "Failed to malloc geometry.objects.");

        chk_msg(mli_set_geometry_objects_and_names_from_archive(
                        &scenery->geometry, &object_names, archive),
                "Failed to malloc geometry.objects.");

        chk_msg(mli_check_malloc_root_frame_from_Archive(
                        &root,
                        archive,
                        &object_names,
                        scenery->geometry.objects,
                        &material_names.boundary_layers),
                "Failed to malloc and populate tree of frames.");

        chk_msg(mliFrame_estimate_num_robjects_and_total_num_boundary_layers(
                        &root, &num_robjects, &total_num_boundary_layers),
                "Can not estimate num_robjects from tree of frames.");

        chk_msg(mliGeometryToMaterialMap_malloc(
                        &scenery->geomap,
                        num_robjects,
                        total_num_boundary_layers),
                "Failed to malloc geometry to materials map.");

        chk_msg(mliGeometry_malloc_references(&scenery->geometry, num_robjects),
                "Failed to malloc geometry.references.");

        chk_msg(mliFrame_set_robjects_and_material_map(
                        &root, &scenery->geometry, &scenery->geomap),
                "Can not set robjects.");

        mliNameMap_free(&material_names);
        mliDynMap_free(&object_names);
        mliFrame_free(&root);

        chk_msg(mliAccelerator_malloc_from_Geometry(
                        &scenery->accelerator, &scenery->geometry),
                "Failed to malloc accelerator from geometry.");

        chk_msg(mliScenery_valid(scenery), "Expected scenery to be valid.");
        chk_msg(mliGeometry_warn_objects(&scenery->geometry),
                "Failed to warn about objects.");

        return 1;
error:
        return 0;
}

/* mliScenery_valid */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliScenery_valid(const struct mliScenery *scenery)
{
        /* check in itself */
        chk_msg(mliMaterials_valid(&scenery->materials),
                "Expected materials to be valid.");
        chk_msg(mliGeometry_valid(&scenery->geometry),
                "Expected geometry to be valid.");
        chk_msg(mliAccelerator_valid(&scenery->accelerator),
                "Expected accelerator to be valid");
        chk_msg(mliGeometryToMaterialMap_valid(&scenery->geomap),
                "Expected geometry-to-materials-map to be valid.");

        /* check interplay */
        chk_msg(mliAccelerator_valid_wrt_Geometry(
                        &scenery->accelerator, &scenery->geometry),
                "Expected accelerator to be valid w.r.t. geometry.");
        chk_msg(mliGeometryToMaterialMap_valid_wrt_Geometry(
                        &scenery->geomap, &scenery->geometry),
                "Expected geomap to be valid w.r.t. geometry.");
        chk_msg(mliGeometryToMaterialMap_valid_wrt_Materials(
                        &scenery->geomap, &scenery->materials),
                "Expected geomap to be valid w.r.t. materials.");
        return 1;
error:
        return 0;
}

/* mliStr */
/* ------ */

/* Copyright 2018-2023 Sebastian Achim Mueller */

struct mliStr mliStr_init(void)
{
        struct mliStr str;
        str.length = 0u;
        str.cstr = NULL;
        return str;
}

void mliStr_free(struct mliStr *str)
{
        free(str->cstr);
        (*str) = mliStr_init();
}

int mliStr_malloc(struct mliStr *str, const uint64_t length)
{
        mliStr_free(str);
        str->length = length;
        chk_malloc(str->cstr, char, str->length + 1);
        memset(str->cstr, '\0', str->length + 1);
        return 1;
error:
        return 0;
}

int mliStr_malloc_copy(struct mliStr *str, const struct mliStr *src)
{
        return mliStr_malloc_copyn(str, src, 0, src->length);
}

int mliStr_malloc_copyn(
        struct mliStr *str,
        const struct mliStr *src,
        const uint64_t start,
        const uint64_t length)
{
        chk_msg(src->cstr != NULL, "Expected src to be allocated");
        chk_msg(start + length <= src->length,
                "Expected start + length < src->length.")
                mliStr_malloc(str, length);
        strncpy(str->cstr, &src->cstr[start], length);
        return 1;
error:
        return 0;
}

int mliStr_mallocf(struct mliStr *str, const char *format, ...)
{
        struct mliStr tmp = mliStr_init();
        va_list args;

        chk(mliStr_malloc(&tmp, 10 * strlen(format)));
        va_start(args, format);
        vsprintf(tmp.cstr, format, args);
        chk(mliStr_malloc(str, strlen(tmp.cstr)));
        strncpy(str->cstr, tmp.cstr, str->length);
        va_end(args);
        mliStr_free(&tmp);
        return 1;
error:
        va_end(args);
        mliStr_free(&tmp);
        return 0;
}

int mliStr_malloc_cstr(struct mliStr *str, const char *s)
{
        chk(mliStr_malloc(str, strlen(s)));
        strncpy(str->cstr, s, str->length);
        return 1;
error:
        return 0;
}

int mliStr_ends_with(const struct mliStr *str, const struct mliStr *suffix)
{
        if (!str->cstr || !suffix->cstr) {
                return 0;
        }
        if (suffix->length > str->length) {
                return 0;
        }
        return strncmp(str->cstr + str->length - suffix->length,
                       suffix->cstr,
                       suffix->length) == 0;
}

int mliStr_starts_with(const struct mliStr *str, const struct mliStr *prefix)
{
        if (!str->cstr || !prefix->cstr) {
                return 0;
        }
        if (prefix->length > str->length) {
                return 0;
        }
        return strncmp(str->cstr, prefix->cstr, prefix->length) == 0;
}

int mliStr_has_prefix_suffix(
        const struct mliStr *str,
        const struct mliStr *prefix,
        const struct mliStr *suffix)
{
        uint64_t has_pre = 1;
        uint64_t has_suf = 1;
        if (prefix->cstr != NULL) {
                has_pre = mliStr_starts_with(str, prefix);
        }

        if (suffix->cstr != NULL) {
                has_suf = mliStr_ends_with(str, suffix);
        }

        if (has_pre == 1 && has_suf == 1) {
                return 1;
        } else {
                return 0;
        }
}

int64_t mliStr_rfind(const struct mliStr *str, const char c)
{
        int64_t i;
        for (i = str->length - 1; i >= 0; i--) {
                if (str->cstr[i] == '\0') {
                        return -1;
                } else if (str->cstr[i] == c) {
                        return i;
                }
        }
        return -1;
}

int64_t mliStr_find(const struct mliStr *str, const char c)
{
        int64_t i;
        for (i = 0; i < (int64_t)str->length; i++) {
                if (str->cstr[i] == '\0') {
                        return -1;
                } else if (str->cstr[i] == c) {
                        return i;
                }
        }
        return -1;
}

int mliStr_match_templeate(
        const struct mliStr *s,
        const struct mliStr *t,
        const char digit_wildcard)
{
        uint64_t i;
        if (s->length != t->length) {
                return 0;
        }
        for (i = 0; i < s->length; i++) {
                if (t->cstr[i] == digit_wildcard) {
                        if (!isdigit(s->cstr[i])) {
                                return 0;
                        }
                } else {
                        if (s->cstr[i] != t->cstr[i]) {
                                return 0;
                        }
                }
        }
        return 1;
}

int mliStr_strip(const struct mliStr *src, struct mliStr *dst)
{
        int64_t start = 0;
        int64_t stop = 0;
        int64_t len = -1;
        struct mliStr cpysrc = mliStr_init();

        chk_msg(src->cstr, "Expected src-string to be allocated.");
        chk_msg(mliStr_malloc_copy(&cpysrc, src), "Can not copy input.");
        mliStr_free(dst);

        while (start < (int64_t)cpysrc.length && isspace(cpysrc.cstr[start])) {
                start += 1;
        }

        stop = cpysrc.length - 1;
        while (stop >= 0 && isspace(cpysrc.cstr[stop])) {
                stop -= 1;
        }

        len = stop - start;

        if (len < 0) {
                chk(mliStr_mallocf(dst, ""));
        } else {
                chk(mliStr_malloc(dst, len + 1));
                strncpy(dst->cstr, &cpysrc.cstr[start], len + 1);
        }
        mliStr_free(&cpysrc);
        return 1;
error:
        mliStr_free(&cpysrc);
        mliStr_free(dst);
        return 0;
}

uint64_t mliStr_countn(
        const struct mliStr *str,
        const char match,
        const uint64_t num_chars_to_scan)
{
        uint64_t i = 0;
        uint64_t count = 0u;
        while (i < str->length && i < num_chars_to_scan) {
                if (str->cstr[i] == match) {
                        count++;
                }
                i++;
        }
        return count;
}

/* mliStr_numbers */
/* -------------- */

/* Copyright 2018-2023 Sebastian Achim Mueller */

int mliStr_nto_double(
        double *out,
        const struct mliStr *str,
        const uint64_t expected_num_chars)
{
        char *end;
        uint64_t actual_num_chars = 0u;
        double l;
        chk_msg(str->cstr != NULL, "Expected str to be allocated.");
        chk_msg(!(str->cstr[0] == '\0' || isspace(str->cstr[0])),
                "Can not convert string to double, bad string.");
        errno = 0;
        l = strtod(str->cstr, &end);
        chk_msg(errno != ERANGE,
                "Can not convert string to double, over-, under-flow.");
        chk_msg(end != NULL, "Can not convert string to double.");

        actual_num_chars = end - str->cstr;
        chk_msg(actual_num_chars == expected_num_chars,
                "double has not the expected number of chars.");
        *out = l;
        return 1;
error:
        return 0;
}

int mliStr_to_double(double *out, const struct mliStr *str)
{
        chk_msg(mliStr_nto_double(out, str, str->length),
                "Can not convert mliStr to double.");
        return 1;
error:
        return 0;
}

int mliStr_nto_int64(
        int64_t *out,
        const struct mliStr *str,
        const uint64_t base,
        const uint64_t expected_num_chars)
{
        char *end;
        uint64_t actual_num_chars = 0u;
        int64_t l;
        chk_msg(str->cstr != NULL, "Expected str to be allocated.");
        chk_msg(!(str->cstr[0] == '\0' || isspace(str->cstr[0])),
                "Can not convert string to int64, bad string.");
        errno = 0;
        l = strtol(str->cstr, &end, base);
        chk_msg(errno != ERANGE,
                "Can not convert string to int64, over-, under-flow.");
        chk_msg(end != NULL, "Can not convert string to int64, bad string.");
        actual_num_chars = end - str->cstr;
        chk_msg(actual_num_chars == expected_num_chars,
                "Integer has not the expected number of chars.");
        *out = l;
        return 1;
error:
        return 0;
}

int mliStr_to_int64(int64_t *out, const struct mliStr *str, const uint64_t base)
{
        chk_msg(mliStr_nto_int64(out, str, base, str->length),
                "Can not convert string to int64.");
        return 1;
error:
        return 0;
}

int mliStr_nto_uint64(
        uint64_t *out,
        const struct mliStr *str,
        const uint64_t base,
        const uint64_t expected_num_chars)
{
        int64_t tmp;
        chk(mliStr_nto_int64(&tmp, str, base, expected_num_chars));
        chk_msg(tmp >= 0, "Expected a positive integer.");
        (*out) = tmp;
        return 1;
error:
        return 0;
}

int mliStr_to_uint64(
        uint64_t *out,
        const struct mliStr *str,
        const uint64_t base)
{
        int64_t tmp;
        chk(mliStr_to_int64(&tmp, str, base));
        chk_msg(tmp >= 0, "Expected a positive integer.");
        (*out) = tmp;
        return 1;
error:
        return 0;
}

int mliStr_reverse_print_uint64(
        const uint64_t u,
        struct mliStr *str,
        const uint64_t base)
{
        char literals[] = {'0',
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

        chk_msg(base <= 16, "Expected base <= 16");
        chk_msg(base > 1, "Expected base > 1");
        mliStr_free(str);

        do {
                remainder = quotient % base;
                quotient = quotient / base;
                remainder32 = (uint32_t)remainder;
                tmp[digs] = literals[remainder32];
                digs++;
                chk_msg(digs < 127, "Exceeded max_num_chars.");
        } while (quotient > 0u);

        chk(mliStr_malloc(str, digs));
        strncpy(str->cstr, tmp, digs);
        return 1;
error:
        mliStr_free(str);
        return 0;
}

int mliStr_print_uint64(
        const uint64_t u,
        struct mliStr *str,
        const uint64_t base,
        const uint64_t min_num_digits,
        const char leading_char)
{
        struct mliStr tmp = mliStr_init();
        int64_t pos = 0;
        int64_t i = 0;
        int64_t length = 0;
        int64_t num_leading = 0;
        int64_t MAX_NUM_CHARS = 128;

        chk_msg(base <= 16, "Expected base <= 16");
        chk_msg(base > 1, "Expected base > 1");

        chk(mliStr_reverse_print_uint64(u, &tmp, base));

        num_leading = min_num_digits - tmp.length;
        if (num_leading < 0) {
                num_leading = 0;
        }
        length = num_leading + tmp.length;
        chk(mliStr_malloc(str, length));
        chk_msg(length < MAX_NUM_CHARS, "Exceeded max_num_chars.");

        pos = 0;
        for (i = 0; i < num_leading; i++) {
                str->cstr[pos] = leading_char;
                pos++;
        }

        for (i = 0; i < (int64_t)tmp.length; i++) {
                str->cstr[pos] = tmp.cstr[(int64_t)tmp.length - i - 1];
                pos++;
        }

        mliStr_free(&tmp);
        return 1;
error:
        mliStr_free(&tmp);
        return 0;
}

/* mliSurface */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliSurface mliSurface_init(void)
{
        struct mliSurface surface;
        surface.material = MLI_MATERIAL_PHONG;
        surface.specular_reflection = mliFunc_init();
        surface.diffuse_reflection = mliFunc_init();
        surface.color = mliColor_set(0.0, 0.0, 0.0);
        return surface;
}

void mliSurface_free(struct mliSurface *surface)
{
        mliFunc_free(&surface->specular_reflection);
        mliFunc_free(&surface->diffuse_reflection);
        (*surface) = mliSurface_init();
}

int mliSurface_malloc(
        struct mliSurface *surface,
        const uint32_t num_points_specular_reflection,
        const uint32_t num_points_diffuse_reflection)
{
        mliSurface_free(surface);
        chk(mliFunc_malloc(
                &surface->specular_reflection, num_points_specular_reflection));
        chk(mliFunc_malloc(
                &surface->diffuse_reflection, num_points_diffuse_reflection));
        return 1;
error:
        return 0;
}

int mliSurface_equal(const struct mliSurface *a, const struct mliSurface *b)
{
        if (a->material != b->material)
                return 0;
        if (!mliFunc_equal(a->specular_reflection, b->specular_reflection))
                return 0;
        if (!mliFunc_equal(a->diffuse_reflection, b->diffuse_reflection))
                return 0;
        if (!mliColor_equal(a->color, b->color))
                return 0;
        return 1;
}

int mli_material_type_to_string(const uint32_t type, char *s)
{
        switch (type) {
        case MLI_MATERIAL_PHONG:
                sprintf(s, "Phong");
                break;
        case MLI_MATERIAL_TRANSPARENT:
                sprintf(s, "transparent");
                break;
        default:
                chk_bad("material-type-id is unknown.");
        }
        return 1;
error:
        return 0;
}

int mli_material_type_from_string(const char *s, uint32_t *id)
{
        if (0 == strcmp(s, "Phong")) {
                (*id) = MLI_MATERIAL_PHONG;
                return 1;
        } else if (0 == strcmp(s, "transparent")) {
                (*id) = MLI_MATERIAL_TRANSPARENT;
                return 1;
        } else {
                chk_bad("material-type-string is unknown.");
        }
        return 1;
error:
        return 0;
}

int mliSurface_fwrite(const struct mliSurface *srf, FILE *f)
{
        struct mliMagicId magic = mliMagicId_init();
        /* magic identifier */
        chk(mliMagicId_set(&magic, "mliSurface"));
        chk_fwrite(&magic, sizeof(struct mliMagicId), 1u, f);

        chk_fwrite(&srf->material, sizeof(uint32_t), 1u, f);
        chk(mliFunc_fwrite(&srf->specular_reflection, f));
        chk(mliFunc_fwrite(&srf->diffuse_reflection, f));
        chk_fwrite(&srf->color, sizeof(struct mliColor), 1u, f);

        return 1;
error:
        return 0;
}

int mliSurface_malloc_fread(struct mliSurface *srf, FILE *f)
{
        struct mliMagicId magic;
        /* magic identifier */
        chk_fread(&magic, sizeof(struct mliMagicId), 1u, f);
        chk_msg(mliMagicId_has_word(&magic, "mliSurface"),
                "Expected MagicID 'mliSurface'.");
        mliMagicId_warn_version(&magic);

        chk_fread(&srf->material, sizeof(uint32_t), 1u, f);
        chk(mliFunc_malloc_fread(&srf->specular_reflection, f));
        chk(mliFunc_malloc_fread(&srf->diffuse_reflection, f));
        chk_fread(&srf->color, sizeof(struct mliColor), 1u, f);

        return 1;
error:
        return 0;
}

/* mliSurface_json */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliSurface_malloc_from_json_str(
        struct mliSurface *surface,
        const char *json_str)
{
        struct mliJson json = mliJson_init();
        chk_msg(mliJson_malloc_from_cstr(&json, json_str),
                "Failed to read json_str to malloc surface.");
        chk_msg(mliSurface_malloc_from_json_token(surface, &json, 0),
                "Failed to malloc surface from json.");
        mliJson_free(&json);
        return 1;
error:
        return 0;
}

int mli_material_type_from_json_token(
        const struct mliJson *json,
        const uint64_t token,
        uint32_t *material)
{
        char buff[MLI_NAME_CAPACITY] = {'\0'};
        const uint64_t name_strlen =
                (json->tokens[token + 1].end - json->tokens[token + 1].start);
        chk_msg(name_strlen < sizeof(buff), "Value of 'name' is too long");
        chk_msg(json->tokens[token + 1].type == JSMN_STRING,
                "Expected 'name' to be of type string.");
        chk_msg(mliJson_cstr_by_token(json, token + 1, buff, name_strlen + 1),
                "Failed to extract string from json.");
        chk_msg(mli_material_type_from_string(buff, material),
                "Failed to parse material type from json-string.");
        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, token + 1);
        return 0;
}

int mliSurface_malloc_from_json_token(
        struct mliSurface *surface,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t token_specref;
        uint64_t token_diffref;
        uint64_t token_color;
        uint64_t token_material;

        /* material */
        /* -------- */
        chk_msg(mliJson_token_by_key(json, token, "material", &token_material),
                "Expected json-surface-item to contain key 'material'.");
        chk_msg(json->tokens[token_material].type == JSMN_STRING,
                "Expected surface's material to be of type string.");
        chk_msg(mli_material_type_from_json_token(
                        json, token_material, &surface->material),
                "Failed to get material-idx from map for string from json");

        /* specular_reflection */
        /* ------------------- */
        chk_msg(mliJson_token_by_key(
                        json, token, "specular_reflection", &token_specref),
                "Expected surface to have key 'specular_reflection'.");
        chk_msg(mliFunc_malloc_from_json_token(
                        &surface->specular_reflection, json, token_specref + 1),
                "Failed to read surface's specular_reflection from json.");

        /* diffuse_reflection */
        /* ------------------ */
        chk_msg(mliJson_token_by_key(
                        json, token, "diffuse_reflection", &token_diffref),
                "Expected surface to have key 'diffuse_reflection'.");
        chk_msg(mliFunc_malloc_from_json_token(
                        &surface->diffuse_reflection, json, token_diffref + 1),
                "Failed to read surface's diffuse_reflection from json.");

        /* color */
        /* ----- */
        chk_msg(mliJson_token_by_key(json, token, "color", &token_color),
                "Expected surface to have key 'color'.");
        chk_msg(mliColor_from_json_token(
                        &surface->color, json, token_color + 1),
                "Failed to read surface's color from json.");

        return 1;
error:
        return 0;
}

/* mliTar */
/* ------ */

/**
 * Copyright (c) 2017 rxi
 * Copyright (c) 2019 Sebastian A. Mueller
 *                    Max-Planck-Institute for nuclear-physics, Heidelberg
 */


/*                             basics                                         */
/* ========================================================================== */

uint64_t mliTar_round_up(uint64_t n, uint64_t incr)
{
        return n + (incr - n % incr) % incr;
}

int mliTar_field_to_uint(
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
error:
        return 0;
}

int mliTar_uint_to_field(
        const uint64_t val,
        char *field,
        const uint64_t fieldsize)
{
        chk(mli_cstr_print_uint64(
                val, field, fieldsize, MLI_TAR_OCTAL, fieldsize - 1));
        return 1;
error:
        return 0;
}

int mliTar_uint64_to_field12_2001star_base256(uint64_t val, char *field)
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
error:
        return 0;
}

int mliTar_field12_to_uint64_2001star_base256(const char *field, uint64_t *val)
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
error:
        return 0;
}

/*                               raw header                                   */
/* ========================================================================== */

uint64_t mliTarRawHeader_checksum(const struct mliTarRawHeader *rh)
{
        uint64_t i;
        unsigned char *p = (unsigned char *)rh;
        uint64_t res = 256;
        for (i = 0; i < offsetof(struct mliTarRawHeader, checksum); i++) {
                res += p[i];
        }
        for (i = offsetof(struct mliTarRawHeader, type); i < sizeof(*rh); i++) {
                res += p[i];
        }
        return res;
}

int mliTarRawHeader_is_null(const struct mliTarRawHeader *rh)
{
        uint64_t i = 0u;
        unsigned char *p = (unsigned char *)rh;
        for (i = 0; i < sizeof(struct mliTarRawHeader); i++) {
                if (p[i] != '\0') {
                        return 0;
                }
        }
        return 1;
}

int mliTarRawHeader_from_header(
        struct mliTarRawHeader *rh,
        const struct mliTarHeader *h)
{
        uint64_t chksum;

        /* Load header into raw header */
        memset(rh, 0, sizeof(*rh));
        chk_msg(mliTar_uint_to_field(h->mode, rh->mode, sizeof(rh->mode)),
                "bad mode");
        chk_msg(mliTar_uint_to_field(h->owner, rh->owner, sizeof(rh->owner)),
                "bad owner");
        if (h->size >= MLI_TAR_MAX_FILESIZE_OCTAL) {
                chk_msg(mliTar_uint64_to_field12_2001star_base256(
                                h->size, rh->size),
                        "bad size, mode: base-256");
        } else {
                chk_msg(mliTar_uint_to_field(
                                h->size, rh->size, sizeof(rh->size)),
                        "bad size, mode: base-octal");
        }
        chk_msg(mliTar_uint_to_field(h->mtime, rh->mtime, sizeof(rh->mtime)),
                "bad mtime");
        rh->type = h->type ? h->type : MLI_TAR_NORMAL_FILE;
        memcpy(rh->name, h->name, sizeof(rh->name));
        memcpy(rh->linkname, h->linkname, sizeof(rh->linkname));

        /* Calculate and write checksum */
        chksum = mliTarRawHeader_checksum(rh);
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
error:
        return 0;
}

/*                                  header                                    */
/* ========================================================================== */

struct mliTarHeader mliTarHeader_init(void)
{
        struct mliTarHeader h;
        h.mode = 0;
        h.owner = 0;
        h.size = 0;
        h.mtime = 0;
        h.type = 0;
        memset(h.name, '\0', sizeof(h.name));
        memset(h.linkname, '\0', sizeof(h.linkname));
        return h;
}

int mliTarHeader_set_directory(struct mliTarHeader *h, const char *name)
{
        (*h) = mliTarHeader_init();
        chk_msg(strlen(name) < sizeof(h->name), "Dirname is too long.");
        memcpy(h->name, name, strlen(name));
        h->type = MLI_TAR_DIRECTORY;
        h->mode = 0775;
        return 1;
error:
        return 0;
}

int mliTarHeader_set_normal_file(
        struct mliTarHeader *h,
        const char *name,
        const uint64_t size)
{
        (*h) = mliTarHeader_init();
        chk_msg(strlen(name) < sizeof(h->name), "Filename is too long.");
        memcpy(h->name, name, strlen(name));
        h->size = size;
        h->type = MLI_TAR_NORMAL_FILE;
        h->mode = 0664;
        return 1;
error:
        return 0;
}

int mliTarHeader_from_raw(
        struct mliTarHeader *h,
        const struct mliTarRawHeader *rh)
{
        uint64_t chksum_actual, chksum_expected;
        chksum_actual = mliTarRawHeader_checksum(rh);

        /* Build and compare checksum */
        chk_msg(mliTar_field_to_uint(
                        &chksum_expected, rh->checksum, sizeof(rh->checksum)),
                "bad checksum string.");
        chk_msg(chksum_actual == chksum_expected, "bad checksum.");

        /* Load raw header into header */
        chk_msg(mliTar_field_to_uint(&h->mode, rh->mode, sizeof(rh->mode)),
                "bad mode");
        chk_msg(mliTar_field_to_uint(&h->owner, rh->owner, sizeof(rh->owner)),
                "bad owner");
        if (rh->size[0] == -128) {
                chk_msg(mliTar_field12_to_uint64_2001star_base256(
                                rh->size, &h->size),
                        "bad size, mode: base-256");
        } else {
                chk_msg(mliTar_field_to_uint(
                                &h->size, rh->size, sizeof(rh->size)),
                        "bad size, mode: base-octal");
        }
        chk_msg(mliTar_field_to_uint(&h->mtime, rh->mtime, sizeof(rh->mtime)),
                "bad mtime");
        h->type = rh->type;
        memcpy(h->name, rh->name, sizeof(h->name));
        memcpy(h->linkname, rh->linkname, sizeof(h->linkname));

        return 1;
error:
        return 0;
}

/* tar */
/* === */

struct mliTar mliTar_init(void)
{
        struct mliTar out;
        out.stream = NULL;
        out.pos = 0u;
        out.remaining_data = 0u;
        return out;
}

/*                                 read                                       */
/* ========================================================================== */

int mliTar_read_begin(struct mliTar *tar, FILE *stream)
{
        chk_msg(tar->stream == NULL,
                "Can't begin reading tar. "
                "tar is either still open or not initialized.");
        (*tar) = mliTar_init();
        tar->stream = stream;
        chk_msg(tar->stream, "Can't begin reading tar. Tar->stream is NULL.");
        return 1;
error:
        return 0;
}

int mliTar_tread(struct mliTar *tar, void *data, const uint64_t size)
{
        int64_t res = fread(data, 1, size, tar->stream);
        chk_msg(res >= 0, "Failed reading from tar.");
        chk_msg((uint64_t)res == size, "Failed reading from tar.");
        tar->pos += size;
        return 1;
error:
        return 0;
}

int mliTar_read_header(struct mliTar *tar, struct mliTarHeader *h)
{
        struct mliTarRawHeader rh;

        chk_msg(mliTar_tread(tar, &rh, sizeof(rh)),
                "Failed to read raw header");

        if (mliTarRawHeader_is_null(&rh)) {
                (*h) = mliTarHeader_init();
                return 0;
        }

        chk_msg(mliTarHeader_from_raw(h, &rh), "Failed to parse raw header.");
        tar->remaining_data = h->size;
        return 1;
error:
        return 0;
}

int mliTar_read_data(struct mliTar *tar, void *ptr, uint64_t size)
{
        chk_msg(tar->remaining_data >= size,
                "Expect size to be read >= remaining_data");
        chk_msg(mliTar_tread(tar, ptr, size), "Failed to read payload-data.");
        tar->remaining_data -= size;

        if (tar->remaining_data == 0) {
                uint64_t i;
                const uint64_t next_record = mliTar_round_up(tar->pos, 512);
                const uint64_t padding_size = next_record - tar->pos;
                char padding;

                for (i = 0; i < padding_size; i++) {
                        chk_msg(mliTar_tread(tar, &padding, 1),
                                "Failed to read padding-block "
                                "to reach next record.");
                }
        }

        return 1;
error:
        return 0;
}

int mliTar_read_finalize(struct mliTar *tar)
{
        struct mliTarHeader h = mliTarHeader_init();
        chk_msg(mliTar_read_header(tar, &h) == 0,
                "Failed to read the 2nd final block of zeros.");
        chk(h.mode == 0);
        chk(h.owner == 0);
        chk(h.size == 0);
        chk(h.mtime == 0);
        chk(h.type == 0);
        chk(h.name[0] == '\0');
        chk(h.linkname[0] == '\0');
        return 1;
error:
        return 0;
}

/*                                  write                                     */
/* ========================================================================== */

int mliTar_write_begin(struct mliTar *tar, FILE *stream)
{
        chk_msg(tar->stream == NULL,
                "Can't begin writing tar. "
                "tar is either still open or not initialized.");
        (*tar) = mliTar_init();
        tar->stream = stream;
        chk_msg(tar->stream, "Can't begin writing tar. Tar->stream is NULL.");
        return 1;
error:
        return 0;
}

int mliTar_twrite(struct mliTar *tar, const void *data, const uint64_t size)
{
        int64_t res = fwrite(data, 1, size, tar->stream);
        chk_msg(res >= 0, "Failed writing to tar.");
        chk_msg((uint64_t)res == size, "Failed writing to tar.");
        tar->pos += size;
        return 1;
error:
        return 0;
}

int mliTar_write_header(struct mliTar *tar, const struct mliTarHeader *h)
{
        struct mliTarRawHeader rh;
        chk_msg(mliTarRawHeader_from_header(&rh, h),
                "Failed to make raw-header");
        tar->remaining_data = h->size;
        chk_msg(mliTar_twrite(tar, &rh, sizeof(rh)), "Failed to write header.");
        return 1;
error:
        return 0;
}

int mliTar_write_null_bytes(struct mliTar *tar, uint64_t n)
{
        uint64_t i;
        char nul = '\0';
        for (i = 0; i < n; i++) {
                chk_msg(mliTar_twrite(tar, &nul, 1), "Failed to write nulls");
        }
        return 1;
error:
        return 0;
}

int mliTar_write_data(struct mliTar *tar, const void *data, uint64_t size)
{
        chk_msg(tar->remaining_data >= size,
                "Expect tar->remaining_data >= size to be written.");
        chk_msg(mliTar_twrite(tar, data, size),
                "Failed to write payload-data.");
        tar->remaining_data -= size;

        if (tar->remaining_data == 0) {
                const uint64_t next_record = mliTar_round_up(tar->pos, 512);
                const uint64_t padding_size = next_record - tar->pos;
                chk_msg(mliTar_write_null_bytes(tar, padding_size),
                        "Failed to write padding zeros.");
        }
        return 1;
error:
        return 0;
}

int mliTar_write_finalize(struct mliTar *tar)
{
        chk_msg(mliTar_write_null_bytes(
                        tar, sizeof(struct mliTarRawHeader) * 2),
                "Failed to write two final null records.");
        return 1;
error:
        return 0;
}

/* mliTmpOcTree */
/* ------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

uint64_t mli_guess_octree_depth_based_on_num_objects(const uint64_t num_objects)
{
        return 3u + (uint64_t)ceil(log((double)num_objects) / log(8.0));
}

/*
 * The dynamic node
 * ================
 */

struct mliTmpNode mliTmpNode_init(void)
{
        struct mliTmpNode n;
        uint64_t c;
        for (c = 0; c < 8u; c++) {
                n.children[c] = NULL;
        }
        n.num_objects = 0u;
        n.objects = NULL;
        n.flat_index = MLI_TMPNODE_FLAT_INDEX_NONE;
        n.node_index = MLI_TMPNODE_FLAT_INDEX_NONE;
        n.leaf_index = MLI_TMPNODE_FLAT_INDEX_NONE;
        return n;
}

void mliTmpNode_free(struct mliTmpNode *n)
{
        uint32_t c;
        for (c = 0; c < 8u; c++)
                if (n->children[c] != NULL)
                        mliTmpNode_free(n->children[c]);
        free(n->objects);
}

int mliTmpNode_malloc(struct mliTmpNode *n, const uint32_t num_objects)
{
        mliTmpNode_free(n);
        n->num_objects = num_objects;
        chk_malloc(n->objects, uint32_t, n->num_objects);
        return 1;
error:
        return 0;
}

uint32_t mliTmpNode_signs_to_child(
        const uint32_t sx,
        const uint32_t sy,
        const uint32_t sz)
{
        return 4 * sx + 2 * sy + 1 * sz;
}

int mliTmpNode_add_children(
        struct mliTmpNode *node,
        const void *bundle,
        int (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mliAABB),
        const struct mliCube cube,
        const uint64_t depth,
        const uint64_t max_depth)
{
        uint32_t c;
        uint32_t sx, sy, sz, obj;
        struct mliCube child_cubes[8];
        struct mliOctOverlap overlap[8];

        if (node->num_objects <= 32u) {
                return 1;
        }

        if (depth == max_depth) {
                return 1;
        }

        for (c = 0u; c < 8u; c++) {
                overlap[c] = mliOctOverlap_init();
                chk(mliOctOverlap_malloc(&overlap[c], node->num_objects));
        }

        /* sense possible children */
        for (sx = 0u; sx < 2u; sx++) {
                for (sy = 0u; sy < 2u; sy++) {
                        for (sz = 0u; sz < 2u; sz++) {
                                const uint32_t child =
                                        mliTmpNode_signs_to_child(sx, sy, sz);
                                child_cubes[child] =
                                        mliCube_octree_child(cube, sx, sy, sz);
                                for (obj = 0u; obj < node->num_objects; obj++) {
                                        const uint32_t object_idx =
                                                node->objects[obj];
                                        if (item_in_bundle_has_overlap_aabb(
                                                    bundle,
                                                    object_idx,
                                                    mliCube_to_aabb(
                                                            child_cubes
                                                                    [child]))) {
                                                chk(mliOctOverlap_push_back(
                                                        &overlap[child],
                                                        object_idx));
                                        }
                                }
                        }
                }
        }

        for (c = 0; c < 8u; c++) {
                chk_malloc(node->children[c], struct mliTmpNode, 1u);
                (*node->children[c]) = mliTmpNode_init();
                chk(mliTmpNode_malloc(node->children[c], overlap[c].size));
                MLI_NCPY(
                        overlap[c].array,
                        node->children[c]->objects,
                        overlap[c].size);
        }

        for (c = 0u; c < 8u; c++) {
                mliOctOverlap_free(&overlap[c]);
        }

        for (c = 0; c < 8u; c++) {
                mliTmpNode_add_children(
                        node->children[c],
                        bundle,
                        item_in_bundle_has_overlap_aabb,
                        child_cubes[c],
                        depth + 1u,
                        max_depth);
        }

        return 1;
error:
        return 0;
}

int mliTmpNode_malloc_tree_from_bundle(
        struct mliTmpNode *root_node,
        const void *bundle,
        const uint32_t num_items_in_bundle,
        int (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mliAABB),
        const struct mliCube bundle_cube)
{
        uint32_t idx, start_depth, max_depth;
        start_depth = 0u;
        max_depth = mli_guess_octree_depth_based_on_num_objects(
                num_items_in_bundle);

        chk_msg(mliTmpNode_malloc(root_node, num_items_in_bundle),
                "Failed to allocate root-node in dynamic octree.");

        for (idx = 0; idx < root_node->num_objects; idx++) {
                root_node->objects[idx] = idx;
        }
        mliTmpNode_add_children(
                root_node,
                bundle,
                item_in_bundle_has_overlap_aabb,
                bundle_cube,
                start_depth,
                max_depth);
        return 1;
error:
        return 0;
}

int mliTmpNode_num_children(const struct mliTmpNode *node)
{
        uint32_t c, num = 0;
        for (c = 0u; c < 8u; c++)
                if (node->children[c] != NULL)
                        num++;
        return num;
}

void mliTmpNode_print(
        const struct mliTmpNode *node,
        const uint32_t indent,
        const uint32_t child)
{
        uint32_t i;
        uint32_t c;
        uint32_t num_children = mliTmpNode_num_children(node);
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
                        mliTmpNode_print(node->children[c], indent + 2, c);
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

int mliTmpNode_exists_and_has_objects(const struct mliTmpNode *node)
{
        if (node != NULL) {
                if (node->num_objects > 0u) {
                        return 1;
                }
        }
        return 0;
}

void mliTmpNode_set_flat_index_walk(
        struct mliTmpNode *node,
        int32_t *flat_index,
        int32_t *node_index,
        int32_t *leaf_index)
{
        uint64_t c;
        for (c = 0u; c < 8u; c++) {
                if (mliTmpNode_exists_and_has_objects(node->children[c])) {
                        (*flat_index)++;
                        node->children[c]->flat_index = *flat_index;

                        if (mliTmpNode_num_children(node->children[c]) == 0) {
                                node->children[c]->leaf_index = *leaf_index;
                                (*leaf_index)++;
                        } else {
                                (*node_index)++;
                                node->children[c]->node_index = *node_index;
                        }
                }
        }
        for (c = 0u; c < 8u; c++) {
                if (mliTmpNode_exists_and_has_objects(node->children[c])) {
                        mliTmpNode_set_flat_index_walk(
                                node->children[c],
                                flat_index,
                                node_index,
                                leaf_index);
                }
        }
}

void mliTmpNode_set_flat_index(struct mliTmpNode *root_node)
{
        int32_t flat_index = 0;
        int32_t node_index = 0;
        int32_t leaf_index = 0;
        root_node->flat_index = flat_index;

        if (mliTmpNode_num_children(root_node) == 0) {
                root_node->leaf_index = leaf_index;
        } else {
                root_node->node_index = node_index;
        }

        mliTmpNode_set_flat_index_walk(
                root_node, &flat_index, &node_index, &leaf_index);
}

/*
 * Find the number of valid nodes in dynamic tree
 */

void mliTmpNode_num_nodes_leafs_objects_walk(
        const struct mliTmpNode *node,
        uint64_t *num_nodes,
        uint64_t *num_leafs,
        uint64_t *num_object_links)
{
        uint64_t c;
        if (node->node_index != MLI_TMPNODE_FLAT_INDEX_NONE) {
                (*num_nodes)++;
        }
        if (node->leaf_index != MLI_TMPNODE_FLAT_INDEX_NONE) {
                (*num_leafs)++;
                (*num_object_links) += node->num_objects;
        }
        for (c = 0; c < 8u; c++) {
                if (node->children[c] != NULL) {
                        mliTmpNode_num_nodes_leafs_objects_walk(
                                node->children[c],
                                num_nodes,
                                num_leafs,
                                num_object_links);
                }
        }
}

void mliTmpNode_num_nodes_leafs_objects(
        const struct mliTmpNode *root_node,
        uint64_t *num_nodes,
        uint64_t *num_leafs,
        uint64_t *num_object_links)
{
        *num_nodes = 0;
        *num_leafs = 0;
        *num_object_links = 0;
        mliTmpNode_num_nodes_leafs_objects_walk(
                root_node, num_nodes, num_leafs, num_object_links);
}

/*
 * The dynamic octree
 * ==================
 */

struct mliTmpOcTree mliTmpOcTree_init(void)
{
        struct mliTmpOcTree octree;
        octree.cube.lower = mliVec_init(0., 0., 0.);
        octree.cube.edge_length = 0.;
        octree.root = mliTmpNode_init();
        return octree;
}

void mliTmpOcTree_free(struct mliTmpOcTree *octree)
{
        mliTmpNode_free(&octree->root);
}

int mliTmpOcTree_malloc_from_bundle(
        struct mliTmpOcTree *octree,
        const void *bundle,
        const uint32_t num_items_in_bundle,
        int (*item_in_bundle_has_overlap_aabb)(
                const void *,
                const uint32_t,
                const struct mliAABB),
        struct mliAABB bundle_aabb)
{
        mliTmpOcTree_free(octree);
        octree->cube = mliCube_outermost_cube(bundle_aabb);
        chk_msg(mliTmpNode_malloc_tree_from_bundle(
                        &octree->root,
                        bundle,
                        num_items_in_bundle,
                        item_in_bundle_has_overlap_aabb,
                        octree->cube),
                "Failed to allocate dynamic octree from bundle.");
        return 1;
error:
        return 0;
}

void mliTmpOcTree_print(const struct mliTmpOcTree *octree)
{
        mliTmpNode_print(&octree->root, 0u, 0u);
}

/* mliTracer */
/* --------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

double mli_trace_sun_visibility(
        const struct mliScenery *scenery,
        const struct mliVec position,
        const struct mliTracerConfig *config,
        struct mliPrng *prng)
{
        return (1.0 -
                mli_trace_sun_obstruction(scenery, position, config, prng));
}

double mli_trace_sun_obstruction(
        const struct mliScenery *scenery,
        const struct mliVec position,
        const struct mliTracerConfig *config,
        struct mliPrng *prng)
{
        uint64_t i;
        double num_obstructions = 0.0;

        for (i = 0; i < config->num_trails_global_light_source; i++) {
                struct mliVec pos_in_source = mliVec_add(
                        mliVec_multiply(
                                config->atmosphere.sunDirection,
                                config->atmosphere.sunDistance),
                        mliVec_multiply(
                                mli_random_position_inside_unit_sphere(prng),
                                config->atmosphere.sunRadius));

                struct mliRay line_of_sight_to_source = mliRay_set(
                        position, mliVec_substract(pos_in_source, position));

                struct mliIntersection isec;

                const int has_intersection = mli_query_intersection(
                        scenery, line_of_sight_to_source, &isec);

                if (has_intersection) {
                        num_obstructions += 1.0;
                }
        }

        return num_obstructions / config->num_trails_global_light_source;
}

struct mliColor mli_trace_to_intersection(
        const struct mliTracerConfig *config,
        const struct mliIntersectionSurfaceNormal *intersection,
        const struct mliScenery *scenery,
        struct mliPrng *prng)
{
        struct mliColor color;
        struct mliSide side;
        double theta;
        double lambert_factor;

        const double sun_visibility = mli_trace_sun_visibility(
                scenery, intersection->position, config, prng);

        side = mli_get_side_coming_from(scenery, intersection);
        color = scenery->materials.surfaces[side.surface].color;

        theta = mliVec_angle_between(
                config->atmosphere.sunDirection, intersection->surface_normal);
        lambert_factor = fabs(cos(theta));

        color.r = color.r * 0.5 * (1.0 + sun_visibility * lambert_factor);
        color.g = color.g * 0.5 * (1.0 + sun_visibility * lambert_factor);
        color.b = color.b * 0.5 * (1.0 + sun_visibility * lambert_factor);

        return color;
}

struct mliColor mli_trace_without_atmosphere(
        const struct mliScenery *scenery,
        const struct mliRay ray,
        const struct mliTracerConfig *config,
        struct mliPrng *prng)
{
        struct mliIntersectionSurfaceNormal intersection =
                mliIntersectionSurfaceNormal_init();

        if (mli_query_intersection_with_surface_normal(
                    scenery, ray, &intersection)) {
                return mli_trace_to_intersection(
                        config, &intersection, scenery, prng);
        } else {
                return config->background_color;
        }
}

struct mliColor mli_trace(
        const struct mliScenery *scenery,
        const struct mliRay ray,
        const struct mliTracerConfig *config,
        struct mliPrng *prng)
{
        if (config->have_atmosphere) {
                return mli_trace_with_atmosphere(scenery, ray, config, prng);
        } else {
                return mli_trace_without_atmosphere(scenery, ray, config, prng);
        }
}

struct mliTracerConfig mliTracerConfig_init(void)
{
        struct mliTracerConfig config;
        config.background_color = mliColor_set(128.0, 128.0, 128.0);
        config.num_trails_global_light_source = 3;
        config.have_atmosphere = 0;
        config.atmosphere = mliAtmosphere_init();
        return config;
}

/* mliTracerConfig_json */
/* -------------------- */

/* Copyright 2018-2021 Sebastian Achim Mueller */

int mliTracerConfig_from_json_token(
        struct mliTracerConfig *tc,
        const struct mliJson *json,
        const uint64_t tkn)
{
        uint64_t bgctkn;
        uint64_t atmtkn;
        uint64_t have_atmosphere;
        chk(mliJson_uint64_by_key(
                json,
                tkn,
                &tc->num_trails_global_light_source,
                "num_trails_global_light_source"));
        chk_msg(tc->num_trails_global_light_source > 0,
                "Expected num_trails_global_light_source > 0.");

        chk(mliJson_uint64_by_key(
                json, tkn, &have_atmosphere, "have_atmosphere"));
        tc->have_atmosphere = (int)have_atmosphere;

        chk(mliJson_token_by_key(json, tkn, "background_color", &bgctkn));
        chk(mliColor_from_json_token(&tc->background_color, json, bgctkn + 1));

        chk(mliJson_token_by_key(json, tkn, "atmosphere", &atmtkn));
        chk(mliAtmosphere_from_json_token(&tc->atmosphere, json, atmtkn + 1));

        return 1;
error:
        return 0;
}

/* mliTracer_atmosphere */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliVec mli_random_direction_in_hemisphere(
        struct mliPrng *prng,
        struct mliVec normal)
{
        struct mliVec rnd_dir;
        do {
                rnd_dir = mli_random_position_inside_unit_sphere(prng);
        } while (mliVec_dot(normal, rnd_dir) <= 0.0);
        return mliVec_normalized(rnd_dir);
}

struct mliColor mli_trace_color_tone_of_sun(
        const struct mliTracerConfig *config,
        const struct mliVec support)
{
        struct mliColor sun_color = mliColor_set(1.0, 1.0, 1.0);
        double width_atmosphere = config->atmosphere.atmosphereRadius -
                                  config->atmosphere.earthRadius;

        if (config->atmosphere.altitude < width_atmosphere) {
                struct mliColor color_close_to_sun = mliAtmosphere_query(
                        &config->atmosphere,
                        support,
                        config->atmosphere.sunDirection);

                double f = config->atmosphere.sunDirection.z;
                double max = MLI_MAX3(
                        color_close_to_sun.r,
                        color_close_to_sun.g,
                        color_close_to_sun.b);

                color_close_to_sun.r /= max;
                color_close_to_sun.g /= max;
                color_close_to_sun.b /= max;

                return mliColor_add(
                        mliColor_multiply(sun_color, f),
                        mliColor_multiply(color_close_to_sun, (1.0 - f)));
        } else {
                return sun_color;
        }
}

struct mliColor mli_trace_color_tone_of_diffuse_sky(
        const struct mliTracerConfig *config,
        const struct mliIntersectionSurfaceNormal *intersection,
        const struct mliScenery *scenery,
        struct mliPrng *prng)
{
        int i;
        struct mliColor sky = mliColor_set(0.0, 0.0, 0.0);
        struct mliRay obstruction_ray;
        struct mliVec facing_surface_normal;
        struct mliIntersection isec;
        int has_direct_view_to_sky = 0;
        int num_samples = 5;

        facing_surface_normal =
                intersection->from_outside_to_inside
                        ? intersection->surface_normal
                        : mliVec_multiply(intersection->surface_normal, -1.0);

        for (i = 0; i < num_samples; i++) {
                struct mliVec rnd_dir = mli_random_direction_in_hemisphere(
                        prng, facing_surface_normal);

                obstruction_ray.support = intersection->position;
                obstruction_ray.direction = rnd_dir;

                has_direct_view_to_sky = !mli_query_intersection(
                        scenery, obstruction_ray, &isec);

                if (has_direct_view_to_sky) {
                        struct mliColor sample = mliAtmosphere_query(
                                &config->atmosphere,
                                intersection->position,
                                rnd_dir);

                        double theta = mliVec_angle_between(
                                rnd_dir, facing_surface_normal);
                        double lambert_factor = fabs(cos(theta));

                        sky = mliColor_add(
                                sky, mliColor_multiply(sample, lambert_factor));
                }
        }

        return mliColor_multiply(sky, 1.0 / (255.0 * num_samples));
}

struct mliColor mli_trace_to_intersection_atmosphere(
        const struct mliTracerConfig *config,
        const struct mliIntersectionSurfaceNormal *intersection,
        const struct mliScenery *scenery,
        struct mliPrng *prng)
{
        struct mliColor color;
        struct mliColor tone;
        struct mliSide side;
        double theta;
        double lambert_factor;

        const double sun_visibility = mli_trace_sun_visibility(
                scenery, intersection->position, config, prng);

        if (sun_visibility > 0.0) {
                tone = mli_trace_color_tone_of_sun(
                        config, intersection->position);
                tone = mliColor_multiply(tone, sun_visibility);
        } else {
                tone = mli_trace_color_tone_of_diffuse_sky(
                        config, intersection, scenery, prng);
        }

        side = mli_get_side_coming_from(scenery, intersection);
        color = scenery->materials.surfaces[side.surface].color;

        theta = mliVec_angle_between(
                config->atmosphere.sunDirection, intersection->surface_normal);
        lambert_factor = fabs(cos(theta));

        color = mliColor_multiply(color, lambert_factor);

        return mliColor_multiply_elementwise(color, tone);
}

struct mliColor mli_trace_with_atmosphere(
        const struct mliScenery *scenery,
        const struct mliRay ray,
        const struct mliTracerConfig *config,
        struct mliPrng *prng)
{
        struct mliIntersectionSurfaceNormal intersection =
                mliIntersectionSurfaceNormal_init();
        struct mliColor out;

        if (mli_query_intersection_with_surface_normal(
                    scenery, ray, &intersection)) {
                out = mli_trace_to_intersection_atmosphere(
                        config, &intersection, scenery, prng);
        } else {
                out = mliAtmosphere_query(
                        &config->atmosphere, ray.support, ray.direction);
        }
        return out;
}

/* mliTriangle_AABB */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

#define MLI_INSIDE 0
#define MLI_OUTSIDE 1

/* Voorhies, Douglas,
 * Triangle-Cube Intersection,
 * Graphics Gems III, p. 236-239, code: p. 521-526 */

/* Which of the six face-plane(s) is point P outside of? */

int64_t mli_triangle_aabb_face_plane(struct mliVec p)
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

int64_t mli_triangle_aabb_bevel_2d(struct mliVec p)
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

int64_t mli_triangle_aabb_bevel_3d(struct mliVec p)
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
        struct mliVec p1,
        struct mliVec p2,
        double alpha,
        int64_t mask)
{
        struct mliVec plane_point;
        plane_point.x = mli_linear_interpolate_1d(alpha, p1.x, p2.x);
        plane_point.y = mli_linear_interpolate_1d(alpha, p1.y, p2.y);
        plane_point.z = mli_linear_interpolate_1d(alpha, p1.z, p2.z);
        return (mli_triangle_aabb_face_plane(plane_point) & mask);
}

/* Compute intersection of P1 --> P2 line segment with face planes */
/* Then test intersection point to see if it is on cube face       */
/* Consider only face planes in "outcode_diff"                     */
/* Note: Zero bits in "outcode_diff" means face line is outside of */

int64_t mli_triangle_aabb_check_line(
        struct mliVec p1,
        struct mliVec p2,
        int64_t outcode_diff)
{
        if ((0x01 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (0.5f - p1.x) / (p2.x - p1.x), 0x3e) ==
                    MLI_INSIDE)
                        return (MLI_INSIDE);
        if ((0x02 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (-0.5f - p1.x) / (p2.x - p1.x), 0x3d) ==
                    MLI_INSIDE)
                        return (MLI_INSIDE);
        if ((0x04 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (0.5f - p1.y) / (p2.y - p1.y), 0x3b) ==
                    MLI_INSIDE)
                        return (MLI_INSIDE);
        if ((0x08 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (-0.5f - p1.y) / (p2.y - p1.y), 0x37) ==
                    MLI_INSIDE)
                        return (MLI_INSIDE);
        if ((0x10 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (0.5f - p1.z) / (p2.z - p1.z), 0x2f) ==
                    MLI_INSIDE)
                        return (MLI_INSIDE);
        if ((0x20 & outcode_diff) != 0)
                if (mli_triangle_aabb_check_point(
                            p1, p2, (-0.5f - p1.z) / (p2.z - p1.z), 0x1f) ==
                    MLI_INSIDE)
                        return (MLI_INSIDE);
        return (MLI_OUTSIDE);
}

/* Test if 3D point is inside 3D triangle */

int64_t mliTriangle_intersects_point(struct mliTriangle t, struct mliVec p)
{
        int64_t sign12_bitmask, sign23_bitmask, sign31_bitmask;
        struct mliVec vect12, vect23, vect31, vect1h, vect2h, vect3h;
        struct mliVec cross12_1p, cross23_2p, cross31_3p;

        /* First, a quick bounding-box test:                               */
        /* If P is outside triangle bbox, there cannot be an intersection. */

        if (p.x > MLI_MAX3(t.v1.x, t.v2.x, t.v3.x))
                return (MLI_OUTSIDE);
        if (p.y > MLI_MAX3(t.v1.y, t.v2.y, t.v3.y))
                return (MLI_OUTSIDE);
        if (p.z > MLI_MAX3(t.v1.z, t.v2.z, t.v3.z))
                return (MLI_OUTSIDE);
        if (p.x < MLI_MIN3(t.v1.x, t.v2.x, t.v3.x))
                return (MLI_OUTSIDE);
        if (p.y < MLI_MIN3(t.v1.y, t.v2.y, t.v3.y))
                return (MLI_OUTSIDE);
        if (p.z < MLI_MIN3(t.v1.z, t.v2.z, t.v3.z))
                return (MLI_OUTSIDE);

        /* For each triangle side, make a vector out of it by subtracting
         * vertexes; */
        /* make another vector from one vertex to point P. */
        /* The crossproduct of these two vectors is orthogonal to both and the
         */
        /* signs of its X,Y,Z components indicate whether P was to the inside or
         */
        /* to the outside of this triangle side. */

        vect12 = mliVec_substract(t.v1, t.v2);
        vect1h = mliVec_substract(t.v1, p);
        cross12_1p = mliVec_cross(vect12, vect1h);
        /*sign12_bitmask = MLI_SIGN3(cross12_1p);*/
        sign12_bitmask = mliVec_sign3_bitmask(cross12_1p, MLI_EPSILON);
        /* Extract X,Y,Z signs as 0..7 or 0...63 integer */

        vect23 = mliVec_substract(t.v2, t.v3);
        vect2h = mliVec_substract(t.v2, p);
        cross23_2p = mliVec_cross(vect23, vect2h);
        sign23_bitmask = mliVec_sign3_bitmask(cross23_2p, MLI_EPSILON);

        vect31 = mliVec_substract(t.v3, t.v1);
        vect3h = mliVec_substract(t.v3, p);
        cross31_3p = mliVec_cross(vect31, vect3h);
        sign31_bitmask = mliVec_sign3_bitmask(cross31_3p, MLI_EPSILON);

        /* If all three crossproduct vectors agree in their component signs,  */
        /* then the point must be inside all three.                           */
        /* P cannot be MLI_OUTSIDE all three sides simultaneously. */

        /* this is the old test; with the revised MLI_SIGN3() macro, the test
         * needs to be revised. */
        return ((sign12_bitmask & sign23_bitmask & sign31_bitmask) == 0)
                       ? MLI_OUTSIDE
                       : MLI_INSIDE;
}

/**********************************************/
/* This is the main algorithm procedure.      */
/* Triangle t is compared with a unit cube,   */
/* centered on the origin.                    */
/* It returns MLI_INSIDE (0) or MLI_OUTSIDE(1) if t   */
/* intersects or does not intersect the cube. */
/**********************************************/

int64_t mliTriangle_intersects_norm_aabb(struct mliTriangle t)
{
        int64_t v1_test, v2_test, v3_test;
        double d, denom;
        struct mliVec vect12, vect13, norm;
        struct mliVec hitpp, hitpn, hitnp, hitnn;

        /* First compare all three vertexes with all six face-planes */
        /* If any vertex is inside the cube, return immediately!     */

        if ((v1_test = mli_triangle_aabb_face_plane(t.v1)) == MLI_INSIDE)
                return (MLI_INSIDE);
        if ((v2_test = mli_triangle_aabb_face_plane(t.v2)) == MLI_INSIDE)
                return (MLI_INSIDE);
        if ((v3_test = mli_triangle_aabb_face_plane(t.v3)) == MLI_INSIDE)
                return (MLI_INSIDE);

        /* If all three vertexes were outside of one or more face-planes, */
        /* return immediately with a trivial rejection!                   */

        if ((v1_test & v2_test & v3_test) != 0)
                return (MLI_OUTSIDE);

        /* Now do the same trivial rejection test for the 12 edge planes */

        v1_test |= mli_triangle_aabb_bevel_2d(t.v1) << 8;
        v2_test |= mli_triangle_aabb_bevel_2d(t.v2) << 8;
        v3_test |= mli_triangle_aabb_bevel_2d(t.v3) << 8;
        if ((v1_test & v2_test & v3_test) != 0)
                return (MLI_OUTSIDE);

        /* Now do the same trivial rejection test for the 8 corner planes */

        v1_test |= mli_triangle_aabb_bevel_3d(t.v1) << 24;
        v2_test |= mli_triangle_aabb_bevel_3d(t.v2) << 24;
        v3_test |= mli_triangle_aabb_bevel_3d(t.v3) << 24;
        if ((v1_test & v2_test & v3_test) != 0)
                return (MLI_OUTSIDE);

        /* If vertex 1 and 2, as a pair, cannot be trivially rejected */
        /* by the above tests, then see if the v1-->v2 triangle edge  */
        /* intersects the cube.  Do the same for v1-->v3 and v2-->v3. */
        /* Pass to the intersection algorithm the "OR" of the outcode */
        /* bits, so that only those cube faces which are spanned by   */
        /* each triangle edge need be tested.                         */

        if ((v1_test & v2_test) == 0)
                if (mli_triangle_aabb_check_line(
                            t.v1, t.v2, v1_test | v2_test) == MLI_INSIDE)
                        return (MLI_INSIDE);
        if ((v1_test & v3_test) == 0)
                if (mli_triangle_aabb_check_line(
                            t.v1, t.v3, v1_test | v3_test) == MLI_INSIDE)
                        return (MLI_INSIDE);
        if ((v2_test & v3_test) == 0)
                if (mli_triangle_aabb_check_line(
                            t.v2, t.v3, v2_test | v3_test) == MLI_INSIDE)
                        return (MLI_INSIDE);

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

        vect12 = mliVec_substract(t.v1, t.v2);
        vect13 = mliVec_substract(t.v1, t.v3);
        norm = mliVec_cross(vect12, vect13);

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
        if (fabs(denom = (norm.x + norm.y + norm.z)) > MLI_EPSILON)
        /* skip parallel diagonals to the plane; division by 0 can occur */
        {
                hitpp.x = hitpp.y = hitpp.z = d / denom;
                if (fabs(hitpp.x) <= 0.5)
                        if (mliTriangle_intersects_point(t, hitpp) ==
                            MLI_INSIDE)
                                return (MLI_INSIDE);
        }
        if (fabs(denom = (norm.x + norm.y - norm.z)) > MLI_EPSILON) {
                hitpn.z = -(hitpn.x = hitpn.y = d / denom);
                if (fabs(hitpn.x) <= 0.5)
                        if (mliTriangle_intersects_point(t, hitpn) ==
                            MLI_INSIDE)
                                return (MLI_INSIDE);
        }
        if (fabs(denom = (norm.x - norm.y + norm.z)) > MLI_EPSILON) {
                hitnp.y = -(hitnp.x = hitnp.z = d / denom);
                if (fabs(hitnp.x) <= 0.5)
                        if (mliTriangle_intersects_point(t, hitnp) ==
                            MLI_INSIDE)
                                return (MLI_INSIDE);
        }
        if (fabs(denom = (norm.x - norm.y - norm.z)) > MLI_EPSILON) {
                hitnn.y = hitnn.z = -(hitnn.x = d / denom);
                if (fabs(hitnn.x) <= 0.5)
                        if (mliTriangle_intersects_point(t, hitnn) ==
                            MLI_INSIDE)
                                return (MLI_INSIDE);
        }

        /* No edge touched the cube; no cube diagonal touched the triangle. */
        /* We're done...there was no intersection.                          */

        return (MLI_OUTSIDE);
}

struct mliTriangle mliTriangle_set_in_norm_aabb(
        const struct mliVec a,
        const struct mliVec b,
        const struct mliVec c,
        const struct mliAABB aabb)
{
        struct mliTriangle tri;
        struct mliVec aabb_center = mliAABB_center(aabb);
        const double inv_scale_x = 1.0 / (aabb.upper.x - aabb.lower.x);
        const double inv_scale_y = 1.0 / (aabb.upper.y - aabb.lower.y);
        const double inv_scale_z = 1.0 / (aabb.upper.z - aabb.lower.z);
        /* translate */
        tri.v1 = mliVec_substract(a, aabb_center);
        tri.v2 = mliVec_substract(b, aabb_center);
        tri.v3 = mliVec_substract(c, aabb_center);
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

int mliTriangle_has_overlap_aabb(
        const struct mliVec a,
        const struct mliVec b,
        const struct mliVec c,
        const struct mliAABB aabb)
{
        struct mliTriangle tri = mliTriangle_set_in_norm_aabb(a, b, c, aabb);
        if (mliTriangle_intersects_norm_aabb(tri) == MLI_INSIDE)
                return 1;
        else
                return 0;
}

struct mliAABB mliTriangle_aabb(
        const struct mliVec a,
        const struct mliVec b,
        const struct mliVec c)
{
        struct mliAABB aabb;
        aabb.lower.x = MLI_MIN3(a.x, b.x, c.x);
        aabb.lower.y = MLI_MIN3(a.y, b.y, c.y);
        aabb.lower.z = MLI_MIN3(a.z, b.z, c.z);
        aabb.upper.x = MLI_MAX3(a.x, b.x, c.x);
        aabb.upper.y = MLI_MAX3(a.y, b.y, c.y);
        aabb.upper.z = MLI_MAX3(a.z, b.z, c.z);
        return aabb;
}

/* mliUserScenery */
/* -------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliNameMap mliNameMap_init(void)
{
        struct mliNameMap nm;
        nm.media = mliDynMap_init();
        nm.surfaces = mliDynMap_init();
        nm.boundary_layers = mliDynMap_init();
        return nm;
}

int mliNameMap_malloc(struct mliNameMap *namemap)
{
        mliNameMap_free(namemap);
        chk_mem(mliDynMap_malloc(&namemap->media, 0u));
        chk_mem(mliDynMap_malloc(&namemap->surfaces, 0u));
        chk_mem(mliDynMap_malloc(&namemap->boundary_layers, 0u));
        return 1;
error:
        return 0;
}

void mliNameMap_free(struct mliNameMap *namemap)
{
        mliDynMap_free(&namemap->media);
        mliDynMap_free(&namemap->surfaces);
        mliDynMap_free(&namemap->boundary_layers);
}

int mli_set_geometry_objects_and_names_from_archive(
        struct mliGeometry *geometry,
        struct mliDynMap *object_names,
        const struct mliArchive *archive)
{
        uint64_t arc_idx = 0u;
        uint64_t obj_idx = 0u;
        char key[MLI_NAME_CAPACITY];

        /* objects */
        obj_idx = 0u;
        for (arc_idx = 0u; arc_idx < mliArchive_num(archive); arc_idx++) {
                if (mli_cstr_has_prefix_suffix(
                            archive->filenames.array[arc_idx].key,
                            "geometry/objects/",
                            ".obj")) {
                        chk_msg(obj_idx < geometry->num_objects,
                                "Expected less objects in archive.");

                        memset(key, '\0', sizeof(key));
                        mli_cstr_path_basename_without_extension(
                                archive->filenames.array[arc_idx].key, key);
                        mli_cstr_path_basename_without_extension(key, key);
                        chk_msg(mliDynMap_insert(object_names, key, obj_idx),
                                "Failed to insert object-filename into map.");
                        chk_msg(mliObject_malloc_from_wavefront(
                                        &geometry->objects[obj_idx],
                                        (char *)archive->textfiles
                                                .array[arc_idx]
                                                .cstr),
                                "Failed to parse wave-front-object.");
                        memcpy(geometry->object_names[obj_idx].cstr,
                               key,
                               MLI_NAME_CAPACITY);

                        obj_idx += 1u;
                }
        }

        return 1;
error:
        return 0;
}

int mliMaterials_malloc_form_archive(
        struct mliMaterials *materials,
        struct mliNameMap *names,
        const struct mliArchive *archive)
{
        uint64_t i = 0u;
        uint64_t med_idx = 0;
        uint64_t srf_idx = 0;
        uint64_t arc_idx = 0;
        char key[MLI_NAME_CAPACITY];

        struct mliStr *default_medium_text = NULL;
        struct mliJson boundary_layers_json = mliJson_init();
        struct mliMaterialsCapacity cap = mliMaterialsCapacity_init();

        /* free */

        mliMaterials_free(materials);
        mliNameMap_free(names);

        chk(mliNameMap_malloc(names));

        /* estimate capacity */
        /* boundary_layers */

        chk_msg(mliArchive_get_malloc_json(
                        archive,
                        "materials/boundary_layers.json",
                        &boundary_layers_json),
                "Failed to parse 'materials/boundary_layers.json'.");

        chk_msg(boundary_layers_json.tokens[0].type == JSMN_OBJECT,
                "Expected key 'boundary_layers' to be a json-object.");
        cap.num_boundary_layers = boundary_layers_json.tokens[0].size;

        cap.num_media = mliArchive_num_filename_prefix_sufix(
                archive, "materials/media/", ".json");

        cap.num_surfaces = mliArchive_num_filename_prefix_sufix(
                archive, "materials/surfaces/", ".json");

        chk_msg(mliMaterials_malloc(materials, cap),
                "Can not malloc materials.");

        /* set fields */
        /* ---------- */

        /* media */
        med_idx = 0u;
        for (arc_idx = 0u; arc_idx < mliArchive_num(archive); arc_idx++) {
                if (mli_cstr_has_prefix_suffix(
                            archive->filenames.array[arc_idx].key,
                            "materials/media/",
                            ".json")) {
                        chk_msg(mliMedium_malloc_from_json_str(
                                        &materials->media[med_idx],
                                        (char *)archive->textfiles
                                                .array[arc_idx]
                                                .cstr),
                                "Failed to parse media json from "
                                "file.");

                        memset(key, '\0', sizeof(key));
                        mli_cstr_path_basename_without_extension(
                                archive->filenames.array[arc_idx].key, key);
                        mli_cstr_path_basename_without_extension(key, key);

                        chk_msg(mliDynMap_insert(&names->media, key, med_idx),
                                "Failed to insert media-name into map.");

                        memcpy(materials->medium_names[med_idx].cstr,
                               names->media.array[med_idx].key,
                               MLI_NAME_CAPACITY);

                        med_idx += 1u;
                }
        }

        /* surfaces */
        srf_idx = 0u;
        for (arc_idx = 0u; arc_idx < mliArchive_num(archive); arc_idx++) {
                if (mli_cstr_has_prefix_suffix(
                            archive->filenames.array[arc_idx].key,
                            "materials/surfaces/",
                            ".json")) {
                        chk_msg(mliSurface_malloc_from_json_str(
                                        &materials->surfaces[srf_idx],
                                        (char *)archive->textfiles
                                                .array[arc_idx]
                                                .cstr),
                                "Failed to parse surface json from "
                                "file.");

                        memset(key, '\0', sizeof(key));
                        mli_cstr_path_basename_without_extension(
                                archive->filenames.array[arc_idx].key, key);
                        mli_cstr_path_basename_without_extension(key, key);

                        chk_msg(mliDynMap_insert(
                                        &names->surfaces, key, srf_idx),
                                "Failed to insert surface-name into map.");

                        memcpy(materials->surface_names[srf_idx].cstr,
                               names->surfaces.array[srf_idx].key,
                               MLI_NAME_CAPACITY);

                        srf_idx += 1u;
                }
        }

        /* boundary_layers */
        chk_msg(mliMaterials_assign_boundary_layers_from_json(
                        materials,
                        &names->boundary_layers,
                        &names->surfaces,
                        &names->media,
                        &boundary_layers_json),
                "Failed to copy boundary_layers from materials.json.");
        for (i = 0; i < materials->num_boundary_layers; i++) {
                memcpy(materials->boundary_layer_names[i].cstr,
                       names->boundary_layers.array[i].key,
                       MLI_NAME_CAPACITY);
        }

        mliJson_free(&boundary_layers_json);

        /* default medium */

        chk_msg(mliArchive_get(
                        archive,
                        "materials/default_medium.txt",
                        &default_medium_text),
                "Can not find 'materials/default_medium.txt' in scenery.");

        memset(key, '\0', sizeof(key));
        mli_cstr_strip_spaces((char *)default_medium_text->cstr, key);

        chk_msg(mliDynMap_get(&names->media, key, &materials->default_medium),
                "Failed to assign the 'default_medium'.");

        return 1;
error:
        mliJson_free(&boundary_layers_json);
        mliMaterials_free(materials);
        mliNameMap_free(names);
        return 0;
}

int mli_check_malloc_root_frame_from_Archive(
        struct mliFrame *root,
        const struct mliArchive *archive,
        const struct mliDynMap *object_names,
        const struct mliObject *objects,
        const struct mliDynMap *boundary_layer_names)
{
        uint64_t token = 0u;
        struct mliJson tree_json = mliJson_init();
        chk_msg(mliArchive_get_malloc_json(
                        archive, "geometry/relations.json", &tree_json),
                "Failed to parse 'geometry/relations.json'.");
        chk_msg(mliJson_token_by_key(&tree_json, 0, "children", &token),
                "Expected 'tree.json' to have key 'children'.");
        chk_msg(mliFrame_malloc(root, MLI_FRAME), "Can not malloc root-frame.");
        chk_msg(mliFrame_from_json(
                        root,
                        &tree_json,
                        token + 1,
                        object_names,
                        objects,
                        boundary_layer_names),
                "Failed to populate tree of Frames from 'tree.json'.");
        mliJson_free(&tree_json);

        /* init transformations */
        mliFrame_set_frame2root(root);

        return 1;
error:
        mliJson_free(&tree_json);
        return 0;
}

/* mliUserScenery_json */
/* ------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliMaterials_assign_boundary_layers_from_json(
        struct mliMaterials *materials,
        struct mliDynMap *boundary_layer_names,
        const struct mliDynMap *surface_names,
        const struct mliDynMap *medium_names,
        const struct mliJson *json)
{
        uint64_t token = 0;
        uint64_t s;

        chk_msg(json->tokens[token].type == JSMN_OBJECT,
                "Expected boundary_layers.json to be a json-object.");

        chk_msg(materials->num_boundary_layers ==
                        (uint32_t)json->tokens[token].size,
                "Expected num_boundary_layers to match "
                "json-object.size.");

        for (s = 0; s < materials->num_boundary_layers; s++) {
                uint64_t token_s_name = mliJson_token_by_index(json, token, s);
                uint64_t token_s = token_s_name + 1;

                chk_msg(json->tokens[token_s_name].type == JSMN_STRING,
                        "Expected boundary_layer to be a String.");

                chk_msg(mliDynMap_insert_key_from_json(
                                boundary_layer_names, json, token_s_name, s),
                        "Failed to insert boundary_layer's name into map.");

                chk_msg(json->tokens[token_s].type == JSMN_OBJECT,
                        "Expected boundary_layer to be of type object {}.");

                chk_msg(mliBoundaryLayer_from_json(
                                &materials->boundary_layers[s],
                                surface_names,
                                medium_names,
                                json,
                                token_s),
                        "Failed to copy boundary_layer from json.");
        }

        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, token);
        return 0;
}

int mliSide_from_json(
        struct mliSide *side,
        const struct mliDynMap *surface_names,
        const struct mliDynMap *medium_names,
        const struct mliJson *json,
        const uint64_t side_token)
{
        uint64_t token_medium, token_surface;

        chk_msg(mliJson_token_by_key(
                        json, side_token + 1, "medium", &token_medium),
                "Expected key 'medium' in side.");
        chk_msg(mliDynMap_get_value_for_string_from_json(
                        medium_names, json, token_medium + 1, &side->medium),
                "Failed to get medium-idx from map");

        chk_msg(mliJson_token_by_key(
                        json, side_token + 1, "surface", &token_surface),
                "Expected key 'surface' in side.");
        chk_msg(mliDynMap_get_value_for_string_from_json(
                        surface_names, json, token_surface + 1, &side->surface),
                "Failed to get surface-idx from map");

        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, side_token + 1);
        return 0;
}

int mliBoundaryLayer_from_json(
        struct mliBoundaryLayer *boundary_layer,
        const struct mliDynMap *surface_names,
        const struct mliDynMap *medium_names,
        const struct mliJson *json,
        const uint64_t token_surface)
{
        uint64_t token_inner_side, token_outer_side;
        chk_msg(mliJson_token_by_key(
                        json, token_surface, "inner", &token_inner_side),
                "Expected key 'inner' in surface.");
        chk_msg(mliJson_token_by_key(
                        json, token_surface, "outer", &token_outer_side),
                "Expected key 'outer' in surface.");

        chk_msg(mliSide_from_json(
                        &boundary_layer->inner,
                        surface_names,
                        medium_names,
                        json,
                        token_inner_side),
                "Failed to parse inner side.");
        chk_msg(mliSide_from_json(
                        &boundary_layer->outer,
                        surface_names,
                        medium_names,
                        json,
                        token_outer_side),
                "Failed to parse outer side.");
        return 1;
error:
        mliJson_debug_token_fprint(stderr, json, token_surface);
        return 0;
}

/* mliVec */
/* ------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliVec mliVec_init(const double x, const double y, const double z)
{
        struct mliVec out;
        out.x = x;
        out.y = y;
        out.z = z;
        return out;
}

struct mliVec mliVec_add(const struct mliVec a, const struct mliVec b)
{
        struct mliVec out;
        out.x = a.x + b.x;
        out.y = a.y + b.y;
        out.z = a.z + b.z;
        return out;
}

struct mliVec mliVec_substract(const struct mliVec a, const struct mliVec b)
{
        struct mliVec out;
        out.x = a.x - b.x;
        out.y = a.y - b.y;
        out.z = a.z - b.z;
        return out;
}

struct mliVec mliVec_cross(const struct mliVec a, const struct mliVec b)
{
        struct mliVec out;
        out.x = (a.y * b.z - a.z * b.y);
        out.y = (a.z * b.x - a.x * b.z);
        out.z = (a.x * b.y - a.y * b.x);
        return out;
}

double mliVec_dot(const struct mliVec a, const struct mliVec b)
{
        return a.x * b.x + a.y * b.y + a.z * b.z;
}

struct mliVec mliVec_multiply(const struct mliVec v, const double a)
{
        struct mliVec out;
        out.x = v.x * a;
        out.y = v.y * a;
        out.z = v.z * a;
        return out;
}

double mliVec_norm(const struct mliVec a) { return sqrt(mliVec_dot(a, a)); }

struct mliVec mliVec_normalized(struct mliVec a)
{
        return mliVec_multiply(a, 1. / mliVec_norm(a));
}

double mliVec_angle_between(const struct mliVec a, const struct mliVec b)
{
        struct mliVec a_normalized = mliVec_multiply(a, 1. / mliVec_norm(a));
        struct mliVec b_normalized = mliVec_multiply(b, 1. / mliVec_norm(b));
        return acos(mliVec_dot(a_normalized, b_normalized));
}

double mliVec_norm_between(const struct mliVec a, const struct mliVec b)
{
        return mliVec_norm(mliVec_substract(a, b));
}

struct mliVec mliVec_mirror(const struct mliVec in, const struct mliVec normal)
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
        struct mliVec out;
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

int mliVec_equal_margin(
        const struct mliVec a,
        const struct mliVec b,
        const double distance_margin)
{
        struct mliVec diff;
        double distance_squared;
        diff = mliVec_substract(a, b);
        distance_squared = mliVec_dot(diff, diff);
        return distance_squared <= distance_margin * distance_margin;
}

int mliVec_equal(const struct mliVec a, const struct mliVec b)
{
        if (fabs(a.x - b.x) > DBL_EPSILON)
                return 0;
        if (fabs(a.y - b.y) > DBL_EPSILON)
                return 0;
        if (fabs(a.z - b.z) > DBL_EPSILON)
                return 0;
        return 1;
}

uint32_t mliVec_octant(const struct mliVec a)
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

int mliVec_sign3_bitmask(const struct mliVec a, const double epsilon)
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

struct mliVec mliVec_mean(const struct mliVec *vecs, const uint64_t num_vecs)
{
        uint64_t i;
        struct mliVec mean = mliVec_init(0.0, 0.0, 0.0);
        for (i = 0; i < num_vecs; i++) {
                mean = mliVec_add(mean, vecs[i]);
        }
        return mliVec_multiply(mean, (1.0 / num_vecs));
}

void mliVec_set(struct mliVec *a, const uint64_t dim, const double v)
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

double mliVec_get(const struct mliVec *a, const uint64_t dim)
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

/* mliVec_AABB */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliVec_overlap_aabb(
        const struct mliVec a,
        const struct mliVec aabb_lower,
        const struct mliVec aabb_upper)
{
        if (a.x >= aabb_lower.x && a.x <= aabb_upper.x && a.y >= aabb_lower.y &&
            a.y <= aabb_upper.y && a.z >= aabb_lower.z && a.z <= aabb_upper.z) {
                return 1;
        } else {
                return 0;
        }
}

/* mliVec_json */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliVec_from_json_token(
        struct mliVec *v,
        const struct mliJson *json,
        const uint64_t token)
{
        chk_msg(json->tokens[token].type == JSMN_ARRAY,
                "Expected vec-token to be a json-array.");
        chk_msg(json->tokens[token].size == 3,
                "Expected vec-token to contain exactly 3 tokens.");
        chk_msg(mliJson_double_by_token(json, token + 1, &v->x),
                "Can not parse mliVec-x-value.");
        chk_msg(mliJson_double_by_token(json, token + 2, &v->y),
                "Can not parse mliVec y-value.");
        chk_msg(mliJson_double_by_token(json, token + 3, &v->z),
                "Can not parse mliVec z-value.");
        return 1;
error:
        return 0;
}

/* mliView */
/* ------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliVec mliView_optical_axis(const struct mliView cam)
{
        struct mliMat rotation = mliMat_init_tait_bryan(
                cam.rotation.x, cam.rotation.y, cam.rotation.z);
        return mli_transform_orientation(&rotation, mliVec_init(0., 0., 1.));
}

struct mliVec mliView_direction_right(const struct mliView cam)
{
        struct mliMat rotation = mliMat_init_tait_bryan(
                cam.rotation.x, cam.rotation.y, cam.rotation.z);
        return mli_transform_orientation(&rotation, mliVec_init(1., 0., 0.));
}

struct mliVec mliView_direction_up(const struct mliView cam)
{
        struct mliMat rotation = mliMat_init_tait_bryan(
                cam.rotation.x, cam.rotation.y, cam.rotation.z);
        return mli_transform_orientation(&rotation, mliVec_init(0., 1., 0.));
}

struct mliView mliView_move_forward(
        const struct mliView camin,
        const double rate)
{
        struct mliView camout = camin;
        struct mliVec optical_axis = mliView_optical_axis(camin);
        camout.position = mliVec_add(
                camout.position, mliVec_multiply(optical_axis, rate));
        return camout;
}

struct mliView mliView_move_right(const struct mliView camin, const double rate)
{
        struct mliView camout = camin;
        struct mliVec direction_right = mliView_direction_right(camout);
        camout.position = mliVec_add(
                camout.position, mliVec_multiply(direction_right, rate));
        return camout;
}

struct mliView mliView_move_up(const struct mliView camin, const double rate)
{
        struct mliView camout = camin;
        camout.position.z += rate;
        return camout;
}

struct mliView mliView_look_right(const struct mliView camin, const double rate)
{
        struct mliView camout = camin;
        const double diff = camin.field_of_view * rate;
        camout.rotation.z = fmod(camout.rotation.z - diff, (2. * MLI_PI));
        return camout;
}

struct mliView mliView_look_down_when_possible(
        const struct mliView camin,
        const double rate)
{
        struct mliView camout = camin;
        const double diff = camin.field_of_view * rate;
        const double next_rotation_x = camout.rotation.x + diff;
        const int fals_forward_over = next_rotation_x > MLI_PI;
        if (fals_forward_over) {
                camout.rotation.x = MLI_PI;
        } else {
                camout.rotation.x = next_rotation_x;
        }
        return camout;
}

struct mliView mliView_increase_fov(
        const struct mliView camin,
        const double rate)
{
        struct mliView camout = camin;
        if (camout.field_of_view * rate > mli_deg2rad(170)) {
                camout.field_of_view = mli_deg2rad(170);
        } else {
                camout.field_of_view *= rate;
        }
        return camout;
}

struct mliView mliView_decrease_fov(
        const struct mliView camin,
        const double rate)
{
        struct mliView camout = camin;
        if (camout.field_of_view / rate < mli_deg2rad(.1)) {
                camout.field_of_view = mli_deg2rad(.1);
        } else {
                camout.field_of_view /= rate;
        }
        return camout;
}

struct mliView mliView_look_up_when_possible(
        const struct mliView camin,
        const double rate)
{
        struct mliView camout = camin;
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

struct mliHomTraComp mliView_to_HomTraComp(const struct mliView view)
{
        struct mliHomTraComp view2root_comp;
        view2root_comp.translation = view.position;
        view2root_comp.rotation = mliQuaternion_set_tait_bryan(
                view.rotation.x, view.rotation.y, view.rotation.z);
        return view2root_comp;
}

/* mli_barycentric */
/* --------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliBarycentrigWeights mli_barycentric_weights(
        const struct mliVec a,
        const struct mliVec b,
        const struct mliVec c,
        const struct mliVec t)
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
        struct mliBarycentrigWeights weights;

        const struct mliVec ab = mliVec_substract(b, a);
        const struct mliVec ac = mliVec_substract(c, a);

        const double ab_ab = mliVec_dot(ab, ab);
        const double ab_ac = mliVec_dot(ab, ac);
        const double ac_ac = mliVec_dot(ac, ac);

        const double abc_area_pow2 = ab_ab * ac_ac - ab_ac * ab_ac;

        const struct mliVec at = mliVec_substract(t, a);
        const double at_ab = mliVec_dot(at, ab);
        const double at_ac = mliVec_dot(at, ac);

        const double atc_area_pow2 = ab_ab * at_ac - ab_ac * at_ab;
        const double abt_area_pow2 = ac_ac * at_ab - ab_ac * at_ac;

        weights.c = atc_area_pow2 / abc_area_pow2;
        weights.b = abt_area_pow2 / abc_area_pow2;
        weights.a = 1.0 - weights.b - weights.c;

        return weights;
}

/* mli_cstr */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

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

int mli_cstr_has_prefix_suffix(
        const char *str,
        const char *prefix,
        const char *sufix)
{
        uint64_t has_pre = 1;
        uint64_t has_suf = 1;
        if (prefix != NULL) {
                has_pre = mli_cstr_starts_with(str, prefix);
        }

        if (sufix != NULL) {
                has_suf = mli_cstr_ends_with(str, sufix);
        }

        if (has_pre == 1 && has_suf == 1) {
                return 1;
        } else {
                return 0;
        }
}

int mli_cstr_split(
        const char *str,
        const char delimiter,
        char *token,
        const uint64_t token_length)
{
        uint64_t i = 0;
        memset(token, '\0', token_length);
        for (i = 0; i < token_length; i++) {
                if (str[i] == '\0') {
                        break;
                } else if (str[i] == delimiter) {
                        break;
                } else {
                        token[i] = str[i];
                }
        }
        return i;
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

uint64_t mli_cstr_count_chars_up_to(
        const char *str,
        const char c,
        const uint64_t num_chars_to_scan)
{
        uint64_t i = 0;
        uint64_t count = 0u;
        while (str[i] != '\0' && i < num_chars_to_scan) {
                if (str[i] == c) {
                        count++;
                }
                i++;
        }
        return count;
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
error:
        return 0;
}

int mli_cstr_lines_fprint(
        FILE *f,
        const char *text,
        const uint64_t line_number,
        const uint64_t line_radius)
{
        int64_t _line_number = (int64_t)line_number;
        int64_t _line_radius = (int64_t)line_radius;
        int64_t line_start = MLI_MAX2(_line_number - _line_radius, 1);
        int64_t line_stop = line_number + line_radius;
        int64_t line = 1;
        int64_t i = 0;

        chk_msg(line_radius > 1, "Expected line_radius > 1.");

        chk(fprintf(f, "  line     text\n"));
        chk(fprintf(f, "        |\n"));

        while (text[i]) {
                int prefix = (line + 1 >= line_start) && (line < line_stop);
                int valid = (line >= line_start) && (line <= line_stop);
                if (text[i] == '\n') {
                        line++;
                }
                if (prefix && i == 0) {
                        chk(mli_fprint_line_match(f, line, _line_number));
                }
                if (valid) {
                        chk(putc(text[i], f));
                }
                if (prefix && text[i] == '\n') {
                        chk(mli_fprint_line_match(f, line, _line_number));
                }
                i++;
        }
        chk(putc('\n', f));

        return 1;
error:
        return 0;
}

void mli_cstr_path_strip_this_dir(char *dst, const char *src)
{
        const char *_src = &src[0];
        memset(dst, '\0', strlen(src));
        while (mli_cstr_starts_with(_src, "./") && _src[0] != '\0') {
                _src += 2;
        }
        strcpy(dst, _src);
}

void mli_cstr_path_basename_without_extension(const char *filename, char *key)
{
        uint64_t i = 0u;
        uint64_t o = 0u;

        while (1) {
                if (filename[i] == '\0') {
                        goto finalize;
                }
                if (filename[i] == '/') {
                        i += 1;
                        break;
                }
                i += 1;
        }

        while (1) {
                if (filename[i] == '\0') {
                        goto finalize;
                }
                if (filename[i] == '.') {
                        i += 1;
                        break;
                }
                key[o] = filename[i];
                i += 1;
                o += 1;
        }

finalize:
        key[o] = '\0';
}

void mli_cstr_strip_spaces(const char *in, char *out)
{
        uint64_t i = 0u;
        uint64_t o = 0u;
        while (in[i] && isspace(in[i])) {
                i += 1;
        }
        while (in[i] && !isspace(in[i])) {
                out[o] = in[i];
                i += 1;
                o += 1;
        }
        out[o] = '\0';
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

/* mli_cstr_numbers */
/* ---------------- */

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
error:
        return 0;
}

int mli_cstr_to_int64(int64_t *out, const char *s, const uint64_t base)
{
        chk_msg(mli_cstr_nto_int64(out, s, base, strlen(s)),
                "Can not convert string to int64.");
        return 1;
error:
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
error:
        return 0;
}

int mli_cstr_to_uint64(uint64_t *out, const char *s, const uint64_t base)
{
        int64_t tmp;
        chk(mli_cstr_to_int64(&tmp, s, base));
        chk_msg(tmp >= 0, "Expected a positive integer.");
        (*out) = tmp;
        return 1;
error:
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
error:
        return 0;
}

int mli_cstr_to_double(double *out, const char *s)
{
        chk_msg(mli_cstr_nto_double(out, s, strlen(s)),
                "Can not convert string to float64.");
        return 1;
error:
        return 0;
}

int mli_cstr_print_uint64(
        uint64_t u,
        char *s,
        const uint64_t max_num_chars,
        const uint64_t base,
        const uint64_t min_num_digits)
{
        char literals[] = {'0',
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
error:
        return 0;
}

/* mli_frame_to_scenery */
/* -------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliFrame_estimate_num_robjects_and_total_num_boundary_layers_walk(
        const struct mliFrame *frame,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers)
{
        uint64_t c;
        switch (frame->type) {
        case MLI_FRAME:
                for (c = 0; c < frame->children.size; c++) {
                        chk(mliFrame_estimate_num_robjects_and_total_num_boundary_layers_walk(
                                frame->children.array[c],
                                num_robjects,
                                total_num_boundary_layers));
                }
                break;
        case MLI_OBJECT:
                (*num_robjects) += 1;
                (*total_num_boundary_layers) += frame->boundary_layers.size;
                break;
        default:
                chk_bad("Expected either type 'frame' or 'object'.");
                break;
        }
        return 1;
error:
        return 0;
}

int mliFrame_estimate_num_robjects_and_total_num_boundary_layers(
        const struct mliFrame *frame,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers)
{
        (*num_robjects) = 0u;
        (*total_num_boundary_layers) = 0u;
        chk_msg(mliFrame_estimate_num_robjects_and_total_num_boundary_layers_walk(
                        frame, num_robjects, total_num_boundary_layers),
                "Failed to walk tree of frames to estimate "
                "num_robjects and total_num_boundary_layers.");
        return 1;
error:
        return 0;
}

int mliFrame_set_robjects_and_material_map_walk(
        const struct mliFrame *frame,
        struct mliGeometry *geometry,
        struct mliGeometryToMaterialMap *geomap,
        uint64_t *num_robjects,
        uint64_t *total_num_boundary_layers)
{
        uint64_t c;
        uint64_t material_idx;
        uint64_t robject_idx;
        switch (frame->type) {
        case MLI_FRAME:
                for (c = 0; c < frame->children.size; c++) {
                        chk(mliFrame_set_robjects_and_material_map_walk(
                                frame->children.array[c],
                                geometry,
                                geomap,
                                num_robjects,
                                total_num_boundary_layers));
                }
                break;
        case MLI_OBJECT:
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
                        mliGeometryToMaterialMap_set(
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
error:
        return 0;
}

int mliFrame_set_robjects_and_material_map(
        const struct mliFrame *frame,
        struct mliGeometry *geometry,
        struct mliGeometryToMaterialMap *geomap)
{
        uint64_t num_robjects = 0u;
        uint64_t total_num_boundary_layers = 0u;
        chk_msg(mliFrame_set_robjects_and_material_map_walk(
                        frame,
                        geometry,
                        geomap,
                        &num_robjects,
                        &total_num_boundary_layers),
                "Failed to walk tree of frames to set "
                "robjects and material map.");
        return 1;
error:
        return 0;
}

/* mli_from_outside_to_inside */
/* -------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_ray_runs_from_outside_to_inside(
        const struct mliVec ray_direction_local,
        const struct mliVec surface_normal_local)
{
        const double proj =
                mliVec_dot(surface_normal_local, ray_direction_local);
        if (proj < 0.)
                return 1;
        else
                return 0;
}

/* mli_intersection_and_scenery */
/* ---------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

uint32_t mliScenery_resolve_boundary_layer_idx(
        const struct mliScenery *scenery,
        const struct mliGeometryId geometry_id)
{
        const uint32_t robject_idx = geometry_id.robj;
        const uint32_t object_idx = scenery->geometry.robjects[robject_idx];
        const uint32_t face_idx = geometry_id.face;
        const uint32_t obj_mtl_idx = mliObject_resolve_material_idx(
                &scenery->geometry.objects[object_idx], face_idx);
        const uint32_t boundary_layer_idx = mliGeometryToMaterialMap_get(
                &scenery->geomap, robject_idx, obj_mtl_idx);
        return boundary_layer_idx;
}

struct mliSide mli_get_side_coming_from(
        const struct mliScenery *scenery,
        const struct mliIntersectionSurfaceNormal *isec)
{
        struct mliBoundaryLayer layer =
                scenery->materials
                        .boundary_layers[mliScenery_resolve_boundary_layer_idx(
                                scenery, isec->geometry_id)];
        if (isec->from_outside_to_inside)
                return layer.outer;
        else
                return layer.inner;
}

struct mliSide mli_get_side_going_to(
        const struct mliScenery *scenery,
        const struct mliIntersectionSurfaceNormal *isec)
{
        struct mliBoundaryLayer layer =
                scenery->materials
                        .boundary_layers[mliScenery_resolve_boundary_layer_idx(
                                scenery, isec->geometry_id)];
        if (isec->from_outside_to_inside)
                return layer.inner;
        else
                return layer.outer;
}

const struct mliFunc *mli_get_refractive_index_going_to(
        const struct mliScenery *scenery,
        const struct mliIntersectionSurfaceNormal *isec)
{
        const struct mliSide going_to = mli_get_side_going_to(scenery, isec);
        return &scenery->materials.media[going_to.medium].refraction;
}

const struct mliFunc *mli_get_refractive_index_coming_from(
        const struct mliScenery *scenery,
        const struct mliIntersectionSurfaceNormal *isec)
{
        const struct mliSide coming_from =
                mli_get_side_coming_from(scenery, isec);
        return &scenery->materials.media[coming_from.medium].refraction;
}

/* mli_json */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliJson mliJson_init(void)
{
        struct mliJson j;
        j.raw = mliStr_init();
        j.num_tokens = 0u;
        j.tokens = NULL;
        return j;
}

void mliJson_free(struct mliJson *json)
{
        mliStr_free(&json->raw);
        free(json->tokens);
        (*json) = mliJson_init();
}

int mliJson_malloc_tokens__(struct mliJson *json)
{
        struct jsmntok_t default_token = {JSMN_UNDEFINED, 0, 0, 0};
        chk_msg(&json->raw.cstr != NULL, "Expected raw cstr to be malloced.");
        json->num_tokens = json->raw.length / 2;
        chk_malloc(json->tokens, struct jsmntok_t, json->num_tokens);
        MLI_ARRAY_SET(json->tokens, default_token, json->num_tokens);
        return 1;
error:
        return 0;
}

int mliJson_parse_tokens__(struct mliJson *json)
{
        int64_t num_tokens_parsed;
        struct jsmn_parser parser;

        chk_msg(&json->tokens != NULL, "Expected tokens to be malloced.");
        jsmn_init(&parser);
        num_tokens_parsed = jsmn_parse(
                &parser,
                json->raw.cstr,
                json->raw.length,
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
error:
        return 0;
}

int mliJson_malloc_from_cstr(struct mliJson *json, const char *cstr)
{
        mliJson_free(json);
        chk_msg(mliStr_malloc_cstr(&json->raw, cstr), "Can't copy cstr.");
        chk_msg(mliJson_malloc_tokens__(json), "Can't malloc Json's tokens.");
        chk_msg(mliJson_parse_tokens__(json), "Can't parse Json into tokens.");
        return 1;
error:
        mliJson_free(json);
        return 0;
}

int mliJson_malloc_from_path(struct mliJson *json, const char *path)
{
        struct mliIo ff = mliIo_init();
        mliJson_free(json);
        chk_msg(mliIo_malloc_from_path(&ff, path),
                "Failed to read file into Json's Str.");
        chk_msg(mliStr_malloc_cstr(&json->raw, (char *)ff.cstr),
                "Failed to copy cstr.");
        mliIo_free(&ff);
        chk_msg(mliJson_malloc_tokens__(json), "Can't malloc Json's tokens.");
        chk_msg(mliJson_parse_tokens__(json), "Can't parse Json into tokens.");
        return 1;
error:
        mliIo_free(&ff);
        mliJson_free(json);
        return 0;
}

int mliJson_cstr_by_token(
        const struct mliJson *json,
        const uint64_t token,
        char *return_string,
        const uint64_t return_string_size)
{
        const struct jsmntok_t t = json->tokens[token];
        const uint64_t actual_length = t.end - t.start;
        chk_msg(actual_length < return_string_size,
                "Expected return_string_size to be sufficiently large for "
                "json-string, but it is not.");
        memcpy(return_string, json->raw.cstr + t.start, actual_length);
        return_string[actual_length] = '\0';
        return 1;
error:
        return 0;
}

int mliJson_int64_by_token(
        const struct mliJson *json,
        const uint64_t token,
        int64_t *return_int64)
{
        const struct jsmntok_t t = json->tokens[token];
        const uint64_t token_length = t.end - t.start;
        chk_msg(t.type == JSMN_PRIMITIVE,
                "Json int64 expected json-token-to be JSMN_PRIMITIVE.");
        chk_msg(mli_cstr_nto_int64(
                        return_int64,
                        (char *)&json->raw.cstr[t.start],
                        10,
                        token_length),
                "Can't parse int64.");
        return 1;
error:
        return 0;
}

int mliJson_uint64_by_token(
        const struct mliJson *json,
        const uint64_t token,
        uint64_t *val)
{
        int64_t tmp;
        chk(mliJson_int64_by_token(json, token, &tmp));
        chk_msg(tmp >= 0, "Expected value to be unsigned.");
        (*val) = (uint64_t)tmp;
        return 1;
error:
        return 0;
}

int mliJson_int64_by_key(
        const struct mliJson *json,
        const uint64_t token,
        int64_t *val,
        const char *key)
{
        uint64_t token_n;
        chk(mliJson_token_by_key_eprint(json, token, key, &token_n));
        chk_msgf(
                mliJson_int64_by_token(json, token_n + 1, val),
                ("Can't parse value of '%s' into int64.", key));
        return 1;
error:
        return 0;
}

int mliJson_uint64_by_key(
        const struct mliJson *json,
        const uint64_t token,
        uint64_t *val,
        const char *key)
{
        int64_t tmp;
        chk(mliJson_int64_by_key(json, token, &tmp, key));
        chk_msg(tmp >= 0, "Expected value to be unsigned.");
        (*val) = (uint64_t)tmp;
        return 1;
error:
        return 0;
}

int mliJson_double_by_token(
        const struct mliJson *json,
        const uint64_t token,
        double *val)
{
        const struct jsmntok_t t = json->tokens[token];
        const uint64_t token_length = t.end - t.start;
        chk_msg(t.type == JSMN_PRIMITIVE,
                "Json float64 expected json-token-to be JSMN_PRIMITIVE.");
        chk_msg(mli_cstr_nto_double(
                        val, (char *)&json->raw.cstr[t.start], token_length),
                "Can't parse double.");
        return 1;
error:
        return 0;
}

int mliJson_double_by_key(
        const struct mliJson *json,
        const uint64_t token,
        double *val,
        const char *key)
{
        uint64_t token_n;
        chk(mliJson_token_by_key_eprint(json, token, key, &token_n));
        chk_msgf(
                mliJson_double_by_token(json, token_n + 1, val),
                ("Can't parse value of '%s' into double.", key));

        return 1;
error:
        return 0;
}

int mliJson_cstrcmp(
        const struct mliJson *json,
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
                const char token_char = (char)json->raw.cstr[t.start + i];
                const char str_char = str[i];
                if (token_char != str_char) {
                        return 0;
                }
        }
        return 1;
}

int mliJson_token_by_key(
        const struct mliJson *json,
        const uint64_t token,
        const char *key,
        uint64_t *key_token)
{
        int64_t found = 0;
        int64_t child = 0;
        int64_t subchild_balance = 0;
        int64_t idx = token + 1;

        while (child < json->tokens[token].size) {
                if (mliJson_cstrcmp(json, idx, key)) {
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

int mliJson_token_by_key_eprint(
        const struct mliJson *json,
        const uint64_t token,
        const char *key,
        uint64_t *key_token)
{
        chk_msgf(
                mliJson_token_by_key(json, token, key, key_token),
                ("Expected key '%s' in json.", key));
        return 1;
error:
        return 0;
}

uint64_t mliJson_token_by_index(
        const struct mliJson *json,
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

int mliJson_debug_token_fprint(
        FILE *f,
        const struct mliJson *json,
        const uint64_t token)
{
        uint64_t i = 0u;
        struct jsmntok_t t = json->tokens[token];
        uint32_t token_size = t.end - t.start;
        uint64_t line_number = 1u + mliStr_countn(&json->raw, '\n', t.start);
        chk(fprintf(f, "line: %u, ", (uint32_t)line_number));
        chk(fprintf(f, "token: %u, ", (uint32_t)token));
        chk(fprintf(f, "type: %d, ", t.type));
        chk(fprintf(f, "children: %d, ", t.size));
        chk(fprintf(f, "chars: (%d -> %d, %d)\n", t.start, t.end, token_size));
        for (i = 0; i < token_size; i++) {
                chk(fputc((char)json->raw.cstr[t.start + i], f));
        }
        chk(fprintf(f, "\n"));
        return 1;
error:
        return 0;
}

int mliJson_debug_fprint(FILE *f, const struct mliJson *json)
{
        uint64_t i;
        for (i = 0; i < json->num_tokens; i++) {
                chk_msg(mliJson_debug_token_fprint(f, json, i),
                        "Failed to write json-token debug-info to file.");
        }
        return 1;
error:
        return 0;
}

int mliJson_debug_to_path(const struct mliJson *json, const char *path)
{
        FILE *f;
        f = fopen(path, "wt");
        chk_msg(f != NULL, "Failed to open file for Json debug output.");
        chk_msg(mliJson_debug_fprint(f, json), "Failed to fprint debug.");
        fclose(f);
        return 1;
error:
        fclose(f);
        return 0;
}

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

/* mli_lambertian_cosine_law */
/* ------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliVec mli_draw_lambertian_direction_wrt_z(struct mliPrng *prng)
{
        double azimuth;
        double sin_theta, cos_theta;
        azimuth = MLI_2PI * mli_random_uniform(prng);
        sin_theta = mli_random_uniform(prng);
        cos_theta = sqrt(1.0 - sin_theta * sin_theta);
        return mliVec_init(
                sin_theta * cos(azimuth), sin_theta * sin(azimuth), cos_theta);
}

struct mliVec mli_draw_lambertian_direction_wrt_surface_normal(
        struct mliPrng *prng,
        const struct mliVec surface_normal)
{
        const struct mliVec z = mliVec_init(0, 0, 1);
        const struct mliVec lambertian_wrt_z =
                mli_draw_lambertian_direction_wrt_z(prng);
        const double rho = mliVec_angle_between(z, surface_normal);
        if (rho > 0.0) {
                const struct mliMat rot = mliMat_init_axis_angle(
                        mliVec_cross(z, surface_normal), -1.0 * rho);
                return mli_transform_orientation(&rot, lambertian_wrt_z);
        } else {
                return lambertian_wrt_z;
        }
}

/* mli_math */
/* -------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

double mli_rad2deg(const double angle_in_rad)
{
        return 180. * angle_in_rad / MLI_PI;
}

double mli_deg2rad(const double angle_in_deg)
{
        return angle_in_deg * (1. / 180.) * MLI_PI;
}

double mli_hypot(const double a, const double b) { return sqrt(a * a + b * b); }

double mli_square(const double a) { return a * a; }

/*
 *  parameters
 *  ----------
 *      points          Sorted array in ascending order.
 *      num_points      Number of points.
 *      point_arg       The point to find the upper-bound for.
 */
uint64_t mli_upper_compare_double(
        const double *points,
        const uint64_t num_points,
        const double point_arg)
{
        uint64_t upper_index = 0;
        MLI_UPPER_COMPARE(points, num_points, point_arg, upper_index);
        return upper_index;
}

void mli_histogram(
        const double *bin_edges,
        const uint64_t num_bin_edges,
        uint64_t *underflow_bin,
        uint64_t *bins,
        uint64_t *overflow_bin,
        const double point)
{
        uint64_t idx_upper =
                mli_upper_compare_double(bin_edges, num_bin_edges, point);
        if (idx_upper == 0) {
                (*underflow_bin) += 1u;
        } else if (idx_upper == num_bin_edges) {
                (*overflow_bin) += 1u;
        } else {
                bins[idx_upper - 1] += 1u;
        }
}

void mli_linspace(
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

double mli_mean(const double vals[], const uint64_t size)
{
        uint64_t i;
        double sum = 0;
        for (i = 0; i < size; i++) {
                sum = sum + vals[i];
        }
        return sum / (double)size;
}

double mli_std(const double vals[], const uint64_t size, const double vals_mean)
{
        uint64_t i;
        double s = 0.;
        for (i = 0; i < size; i++) {
                s = s + (vals[i] - vals_mean) * (vals[i] - vals_mean);
        }
        return sqrt(s / (double)size);
}

double mli_bin_center_in_linear_space(
        const double start,
        const double stop,
        const uint64_t num_bins,
        const uint64_t bin)
{
        const double width = stop - start;
        const double bin_width = width / (double)num_bins;
        return start + bin * bin_width + 0.5 * bin_width;
}

double mli_linear_interpolate_1d(
        const double weight,
        const double start,
        const double end)
{
        return start + weight * (end - start);
}

double mli_linear_interpolate_2d(
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

double mli_relative_ratio(const double a, const double b)
{
        return fabs(a - b) / (0.5 * (a + b));
}

/* mli_photon_propagation */
/* ---------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliPhotonInteraction mliPhotonInteraction_from_Intersection(
        const int64_t type,
        const struct mliScenery *scenery,
        const struct mliIntersectionSurfaceNormal *isec)
{
        struct mliPhotonInteraction phia;

        struct mliSide side_coming_from, side_going_to;

        side_coming_from = mli_get_side_coming_from(scenery, isec);
        side_going_to = mli_get_side_going_to(scenery, isec);

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

int mli_propagate_photon_phong(
        struct mliEnv *env,
        const struct mliIntersectionSurfaceNormal *isec)
{
        double specular;
        double diffuse;
        double rnd;
        struct mliSide side_coming_from =
                mli_get_side_coming_from(env->scenery, isec);

        chk_msg(mliFunc_evaluate(
                        &env->scenery->materials
                                 .surfaces[side_coming_from.surface]
                                 .diffuse_reflection,
                        env->photon->wavelength,
                        &diffuse),
                "Failed to eval. diffuse reflection for wavelength.");
        chk_msg(mliFunc_evaluate(
                        &env->scenery->materials
                                 .surfaces[side_coming_from.surface]
                                 .specular_reflection,
                        env->photon->wavelength,
                        &specular),
                "Failed to eval. specular reflection for wavelength.");
        rnd = mli_random_uniform(env->prng);
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
                chk(mliDynPhotonInteraction_push_back(
                        env->history,
                        mliPhotonInteraction_from_Intersection(
                                MLI_PHOTON_DIFFUSE_REFLECTION,
                                env->scenery,
                                isec)));
                env->photon->ray = mliRay_set(
                        isec->position,
                        mli_draw_lambertian_direction_wrt_surface_normal(
                                env->prng, isec->surface_normal));
                chk_msg(mli_propagate_photon_env(env),
                        "Failed to continue after diffuse reflection phong.");
        } else if (rnd < (specular + diffuse)) {
                chk(mliDynPhotonInteraction_push_back(
                        env->history,
                        mliPhotonInteraction_from_Intersection(
                                MLI_PHOTON_SPECULAR_REFLECTION,
                                env->scenery,
                                isec)));
                env->photon->ray = mliRay_set(
                        isec->position,
                        mliVec_mirror(
                                env->photon->ray.direction,
                                isec->surface_normal));
                chk_msg(mli_propagate_photon_env(env),
                        "Failed to continue after specular reflection phong.");
        } else {
                chk(mliDynPhotonInteraction_push_back(
                        env->history,
                        mliPhotonInteraction_from_Intersection(
                                MLI_PHOTON_ABSORBTION, env->scenery, isec)));
        }
        return 1;
error:
        return 0;
}

int mli_propagate_photon_pass_boundary_layer(
        struct mliEnv *env,
        const struct mliIntersectionSurfaceNormal *isec,
        const struct mliFresnel fresnel)
{
        chk(mliDynPhotonInteraction_push_back(
                env->history,
                mliPhotonInteraction_from_Intersection(
                        MLI_PHOTON_REFRACTION, env->scenery, isec)));
        env->photon->ray = mliRay_set(
                isec->position, mliFresnel_refraction_direction(fresnel));
        chk_msg(mli_propagate_photon_env(env),
                "Failed to continue after passing boundary layer");
        return 1;
error:
        return 0;
}

int mli_propagate_photon_probability_passing_medium_coming_from(
        const struct mliScenery *scenery,
        const struct mliPhoton *photon,
        const struct mliIntersectionSurfaceNormal *isec,
        double *probability_passing)
{
        double one_over_e_way;
        const struct mliSide side_coming_from =
                mli_get_side_coming_from(scenery, isec);
        chk_msg(mliFunc_evaluate(
                        &scenery->materials.media[side_coming_from.medium]
                                 .absorbtion,
                        photon->wavelength,
                        &one_over_e_way),
                "Photon's wavelength is out of range to "
                "evaluate absorbtion in medium coming from");
        (*probability_passing) = exp(-isec->distance_of_ray / one_over_e_way);
        return 1;
error:
        return 0;
}

int mli_propagate_photon_fresnel_refraction_and_reflection(
        struct mliEnv *env,
        const struct mliIntersectionSurfaceNormal *isec)
{
        struct mliFresnel fresnel;
        double n_going_to;
        double n_coming_from;
        double reflection_propability;
        struct mliVec facing_surface_normal;
        chk_msg(mliFunc_evaluate(
                        mli_get_refractive_index_going_to(env->scenery, isec),
                        env->photon->wavelength,
                        &n_going_to),
                "Failed to eval. refraction going to for wavelength.");
        chk_msg(mliFunc_evaluate(
                        mli_get_refractive_index_coming_from(
                                env->scenery, isec),
                        env->photon->wavelength,
                        &n_coming_from),
                "Failed to eval. refraction coming from for wavelength.");
        facing_surface_normal =
                isec->from_outside_to_inside
                        ? isec->surface_normal
                        : mliVec_multiply(isec->surface_normal, -1.0);
        fresnel = mliFresnel_init(
                env->photon->ray.direction,
                facing_surface_normal,
                n_coming_from,
                n_going_to);
        reflection_propability = mliFresnel_reflection_propability(fresnel);
        if (reflection_propability > mli_random_uniform(env->prng)) {
                chk(mliDynPhotonInteraction_push_back(
                        env->history,
                        mliPhotonInteraction_from_Intersection(
                                MLI_PHOTON_FRESNEL_REFLECTION,
                                env->scenery,
                                isec)));
                env->photon->ray = mliRay_set(
                        isec->position,
                        mliFresnel_reflection_direction(fresnel));
                chk_msg(mli_propagate_photon_env(env),
                        "Failed to continue after reflection");
        } else {
                chk_msg(mli_propagate_photon_pass_boundary_layer(
                                env, isec, fresnel),
                        "Failed to pass boundary");
        }
        return 1;
error:
        return 0;
}

int mli_propagate_photon_interact_with_object(
        struct mliEnv *env,
        const struct mliIntersectionSurfaceNormal *isec)
{
        struct mliSurface surface_coming_from;
        const struct mliSide side_coming_from =
                mli_get_side_coming_from(env->scenery, isec);
        surface_coming_from =
                env->scenery->materials.surfaces[side_coming_from.surface];
        switch (surface_coming_from.material) {
        case MLI_MATERIAL_TRANSPARENT:
                chk_msg(mli_propagate_photon_fresnel_refraction_and_reflection(
                                env, isec),
                        "Failed Fresnel.");
                break;
        case MLI_MATERIAL_PHONG:
                chk_msg(mli_propagate_photon_phong(env, isec),
                        "Failed Phong-material.");
                break;
        default:
                chk_bad("Unkown material of surface.");
                break;
        }
        return 1;
error:
        return 0;
}

int mli_propagate_photon_distance_until_absorbtion(
        const struct mliFunc *absorbtion_in_medium_passing_through,
        const double wavelength,
        struct mliPrng *prng,
        double *distance_until_absorbtion)
{
        double one_over_e_way;
        chk_msg(mliFunc_evaluate(
                        absorbtion_in_medium_passing_through,
                        wavelength,
                        &one_over_e_way),
                "Failed to eval. absorbtion for wavelength.");
        (*distance_until_absorbtion) =
                mli_random_expovariate(prng, 1. / one_over_e_way);
        return 1;
error:
        return 0;
}

int mli_propagate_photon_work_on_causal_intersection(struct mliEnv *env)
{
        int ray_does_intersect_surface = 0;
        double distance_until_absorbtion = 0.0;
        struct mliIntersectionSurfaceNormal next_intersection;
        struct mliFunc *absorbtion_in_medium_passing_through;
        struct mliPhotonInteraction phia;

        ray_does_intersect_surface = mli_query_intersection_with_surface_normal(
                env->scenery, env->photon->ray, &next_intersection);

        if (ray_does_intersect_surface) {
                int photon_is_absorbed_before_reaching_surface;
                struct mliSide side_coming_from;

                side_coming_from = mli_get_side_coming_from(
                        env->scenery, &next_intersection);
                absorbtion_in_medium_passing_through =
                        &env->scenery->materials.media[side_coming_from.medium]
                                 .absorbtion;
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
                        phia.geometry_id = mliGeometryId_init();
                        phia.from_outside_to_inside = 1;

                        phia.medium_coming_from = side_coming_from.medium;
                        phia.medium_going_to = side_coming_from.medium;

                        chk(mliDynPhotonInteraction_push_back(
                                env->history, phia));
                }

                if (photon_is_absorbed_before_reaching_surface) {
                        /* absorbtion in medium */
                        phia.type = MLI_PHOTON_ABSORBTION_MEDIUM;
                        phia.position = mliRay_at(
                                &env->photon->ray, distance_until_absorbtion);
                        ;
                        phia.position_local = phia.position;
                        phia.distance_of_ray = distance_until_absorbtion;
                        phia.on_geometry_surface = 0;
                        phia.geometry_id = mliGeometryId_init();
                        phia.from_outside_to_inside = 1;

                        phia.medium_coming_from = side_coming_from.medium;
                        phia.medium_going_to = side_coming_from.medium;

                        chk(mliDynPhotonInteraction_push_back(
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

                absorbtion_in_medium_passing_through =
                        &env->scenery->materials.media[default_medium]
                                 .absorbtion;
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
                        phia.geometry_id = mliGeometryId_init();
                        phia.from_outside_to_inside = 1;

                        phia.medium_coming_from = default_medium;
                        phia.medium_going_to = default_medium;

                        chk(mliDynPhotonInteraction_push_back(
                                env->history, phia));
                }

                /* absorbtion in medium */
                phia.type = MLI_PHOTON_ABSORBTION_MEDIUM;
                phia.position =
                        mliRay_at(&env->photon->ray, distance_until_absorbtion);
                phia.position_local = phia.position;
                phia.distance_of_ray = distance_until_absorbtion;
                phia.on_geometry_surface = 0;
                phia.geometry_id = mliGeometryId_init();
                phia.from_outside_to_inside = 1;

                phia.medium_coming_from = default_medium;
                phia.medium_going_to = default_medium;

                chk(mliDynPhotonInteraction_push_back(env->history, phia));
        }

        return 1;
error:
        return 0;
}

int mli_propagate_photon_env(struct mliEnv *env)
{
        if (env->max_interactions > env->history->size) {
                chk_msg(mli_propagate_photon_work_on_causal_intersection(env),
                        "Failed to work on intersection.");
        }
        return 1;
error:
        return 0;
}

int mli_propagate_photon(
        const struct mliScenery *scenery,
        struct mliDynPhotonInteraction *history,
        struct mliPhoton *photon,
        struct mliPrng *prng,
        const uint64_t max_interactions)
{
        struct mliEnv env;
        env.scenery = scenery;
        env.history = history;
        env.photon = photon;
        env.prng = prng;
        env.max_interactions = max_interactions;
        chk(mli_propagate_photon_env(&env));
        return 1;
error:
        return 0;
}

/* mli_photon_sources */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mli_photon_source_parallel_towards_z_from_xy_disc(
        struct mliDynPhoton *out_photons,
        const double wavelength,
        const double radius,
        const uint64_t num_photons,
        struct mliPrng *prng)
{
        uint64_t i;
        const struct mliVec direction = mliVec_init(0., 0., 1.);
        for (i = 0; i < num_photons; i++) {
                struct mliPhoton ph;
                ph.ray.support = mli_random_position_on_disc(radius, prng);
                ph.ray.direction = direction;
                ph.wavelength = wavelength;
                ph.id = i;
                chk(mliDynPhoton_push_back(out_photons, ph));
        }
        return 1;
error:
        return 0;
}

int mli_photon_source_point_like_opening_cone_towards_z(
        struct mliDynPhoton *out_photons,
        const double wavelength,
        const double opening_angle,
        const uint64_t num_photons,
        struct mliPrng *prng)
{
        uint64_t i;
        struct mliRandomUniformRange azimuth =
                mliRandomUniformRange_set(0.0, 2.0 * MLI_PI);
        struct mliRandomZenithRange zenith =
                mliRandomZenithRange_set(0.0, opening_angle);
        for (i = 0; i < num_photons; i++) {
                struct mliVec direction =
                        mli_random_draw_direction_in_zenith_azimuth_range(
                                zenith, azimuth, prng);
                struct mliPhoton ph;
                ph.ray.support = mliVec_init(0., 0., 0.);
                ph.ray.direction = direction;
                ph.wavelength = wavelength;
                ph.id = i;
                chk(mliDynPhoton_push_back(out_photons, ph));
        }
        return 1;
error:
        return 0;
}

/* mli_quadratic_equation */
/* ---------------------- */

/* Copyright 2019 Sebastian A. Mueller */

int mli_quadratic_equation(
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

/* mli_random */
/* ---------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

double mli_random_uniform(struct mliPrng *prng)
{
        uint32_t rn_int = mliPrng_generate_uint32(prng);
        const double rn = (double)rn_int;
        const double max_uint32 = (double)UINT32_MAX;
        return rn / max_uint32;
}

double mli_random_expovariate(struct mliPrng *prng, const double rate)
{
        /*      Sampling from a poisson distribution */
        return -log(mli_random_uniform(prng)) / rate;
}

double mli_random_normal_Irwin_Hall_approximation(
        struct mliPrng *prng,
        const double mean,
        const double std)
{
        uint64_t i;
        double sum_of_12 = 0.;
        double std1 = 0.;
        for (i = 0; i < 12; i++) {
                sum_of_12 += mli_random_uniform(prng);
        }
        std1 = sum_of_12 - 6.;
        return mean + std1 * std;
}

/*
        uniform linear range
        ====================
        Draw a uniform distribution within a limited range.
*/
struct mliRandomUniformRange mliRandomUniformRange_set(
        double start,
        double stop)
{
        struct mliRandomUniformRange p;
        p.start = start;
        p.range = stop - start;
        assert(p.range >= 0.0);
        return p;
}

double mli_random_draw_uniform(
        const struct mliRandomUniformRange uniform_range,
        struct mliPrng *prng)
{
        return uniform_range.range * mli_random_uniform(prng) +
               uniform_range.start;
}

/*
        zenith range
        ============
        Draw zenith-distances for an even distribution of points on a sphere.
*/
struct mliRandomZenithRange mliRandomZenithRange_set(
        const double min_zenith_distance,
        const double max_zenith_distance)
{
        struct mliRandomZenithRange zp;
        zp.z_min = (cos(min_zenith_distance) + 1.0) / 2.0;
        zp.z_range = (cos(max_zenith_distance) + 1.0) / 2.0 - zp.z_min;
        return zp;
}

double mli_random_draw_zenith(
        const struct mliRandomZenithRange range,
        struct mliPrng *prng)
{
        const double z =
                (range.z_range * mli_random_uniform(prng)) + range.z_min;
        return acos(2.0 * z - 1.0);
}

/*
        direction
        ==========
*/
struct mliVec mli_random_draw_direction_in_zenith_azimuth_range(
        const struct mliRandomZenithRange zenith,
        const struct mliRandomUniformRange azimuth,
        struct mliPrng *prng)
{
        const double az = mli_random_draw_uniform(azimuth, prng);
        const double zd = mli_random_draw_zenith(zenith, prng);
        const double sin_zd = sin(zd);
        return mliVec_init(sin_zd * cos(az), sin_zd * sin(az), cos(zd));
}

struct mliVec mli_random_position_on_disc(
        const double radius,
        struct mliPrng *prng)
{
        const double r = sqrt(mli_random_uniform(prng)) * radius;
        const double azimuth = mli_random_uniform(prng) * MLI_2PI;
        return mliVec_init(r * cos(azimuth), r * sin(azimuth), 0.0);
}

struct mliVec mli_random_position_inside_unit_sphere(struct mliPrng *prng)
{
        /* rejection sampling */
        struct mliVec pos;
        do {
                pos.x = -1.0 + 2.0 * mli_random_uniform(prng);
                pos.y = -1.0 + 2.0 * mli_random_uniform(prng);
                pos.z = -1.0 + 2.0 * mli_random_uniform(prng);
        } while ((pos.x * pos.x + pos.y * pos.y + pos.z * pos.z) > 1.0);
        return pos;
}

/* mli_random_MT19937 */
/* ------------------ */

/* Copyright 2018-2020 Sebastian Achim Mueller */

/*
 *      Adopted from https://en.wikipedia.org/wiki/Mersenne_Twister
 */
void mliMT19937_set_constants(struct mliMT19937 *mt)
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

void mliMT19937_reinit(struct mliMT19937 *mt, const uint32_t seed)
{
        uint32_t i;
        mliMT19937_set_constants(mt);
        mt->mt[0] = seed;
        for (i = 1; i < mt->N; i++) {
                mt->mt[i] =
                        (mt->F * (mt->mt[i - 1] ^ (mt->mt[i - 1] >> 30)) + i);
        }
        mt->index = mt->N;
}

struct mliMT19937 mliMT19937_init(const uint32_t seed)
{
        struct mliMT19937 mt;
        mliMT19937_reinit(&mt, seed);
        return mt;
}

void mliMT19937_twist(struct mliMT19937 *mt)
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

uint32_t mliMT19937_generate_uint32(struct mliMT19937 *mt)
{
        uint32_t y;
        int i = mt->index;
        if (mt->index >= mt->N) {
                mliMT19937_twist(mt);
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

/* mli_random_PCG32 */
/* ---------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mliPCG32 mliPCG32_init(const uint32_t seed)
{
        struct mliPCG32 pcg32;
        pcg_setseq_64_srandom_r(&pcg32.state_setseq_64, seed, 0u);
        return pcg32;
}

uint32_t mliPCG32_generate_uint32(struct mliPCG32 *pcg32)
{
        return pcg_setseq_64_xsh_rr_32_random_r(&pcg32->state_setseq_64);
}

void mliPCG32_reinit(struct mliPCG32 *pcg32, const uint32_t seed)
{
        pcg_setseq_64_srandom_r(&pcg32->state_setseq_64, seed, 0u);
}

/* mli_random_generator */
/* -------------------- */


uint32_t mliPrng_generate_uint32(struct mliPrng *prng)
{
        return prng->generate_uint32(&prng->_storage);
}

void mliPrng_reinit(struct mliPrng *prng, const uint32_t seed)
{
        prng->reinit(&prng->_storage, seed);
}

/**
 *      Mersenne Twister 19937
 *      ----------------------
 */

struct mliPrng mliPrng_init_MT19937(const uint32_t seed)
{
        struct mliPrng prng;
        prng._storage.mt19937 = mliMT19937_init(seed);
        prng.generate_uint32 = mliPrng_MT19937_generate_uint32;
        prng.reinit = mliPrng_MT19937_reinit;
        return prng;
}

uint32_t mliPrng_MT19937_generate_uint32(void *mt)
{
        return mliMT19937_generate_uint32((struct mliMT19937 *)mt);
}

void mliPrng_MT19937_reinit(void *mt, const uint32_t seed)
{
        mliMT19937_reinit((struct mliMT19937 *)mt, seed);
}

/**
 *      PCG32
 *      -----
 */

struct mliPrng mliPrng_init_PCG32(const uint32_t seed)
{
        struct mliPrng prng;
        prng._storage.pcg32 = mliPCG32_init(seed);
        prng.generate_uint32 = mliPrng_PCG32_generate_uint32;
        prng.reinit = mliPrng_PCG32_reinit;
        return prng;
}

uint32_t mliPrng_PCG32_generate_uint32(void *pcg)
{
        return mliPCG32_generate_uint32((struct mliPCG32 *)pcg);
}

void mliPrng_PCG32_reinit(void *pcg, const uint32_t seed)
{
        mliPCG32_reinit((struct mliPCG32 *)pcg, seed);
}

/* mli_random_pcg_variants_32bit_subset */
/* ------------------------------------ */

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


#define PCG_DEFAULT_MULTIPLIER_64 6364136223846793005U
#define PCG_DEFAULT_INCREMENT_64 1442695040888963407U

/**     Rotate helper functions.
 */

uint32_t pcg_rotr_32(uint32_t value, unsigned int rot)
{
        return (value >> rot) | (value << ((-rot) & 31));
}

/**     Output functions.  These are the core of the PCG generation scheme.
 *
 *      XSH RR
 */

uint32_t pcg_output_xsh_rr_64_32(uint64_t state)
{
        return pcg_rotr_32(((state >> 18u) ^ state) >> 27u, state >> 59u);
}

/**     Functions to advance the underlying LCG.
 *      These functions are considered semi-private.
 *      There is rarely a good reason to call them directly.
 */

void pcg_setseq_64_step_r(struct pcg_state_setseq_64 *rng)
{
        rng->state = rng->state * PCG_DEFAULT_MULTIPLIER_64 + rng->inc;
}

/**     Seed the RNG state.
 */

void pcg_setseq_64_srandom_r(
        struct pcg_state_setseq_64 *rng,
        uint64_t initstate,
        uint64_t initseq)
{
        rng->state = 0U;
        rng->inc = (initseq << 1u) | 1u;
        pcg_setseq_64_step_r(rng);
        rng->state += initstate;
        pcg_setseq_64_step_r(rng);
}

/**     Now, finally we provide
 *      a random_r function that provides a random number of the appropriate
 *      type (using the full range of the type).
 *
 *      Generation functions for XSH RR
 */

uint32_t pcg_setseq_64_xsh_rr_32_random_r(struct pcg_state_setseq_64 *rng)
{
        uint64_t oldstate = rng->state;
        pcg_setseq_64_step_r(rng);
        return pcg_output_xsh_rr_64_32(oldstate);
}

/* mli_ray_octree_traversal */
/* ------------------------ */

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

int mli_ray_octree_traversal_first_octree_node(
        const struct mliVec t0,
        const struct mliVec tm)
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

int mli_ray_octree_traversal_next_octree_node(
        const struct mliVec tm,
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
                const uint32_t))
{
        int32_t proc_node;
        struct mliVec tm, nt0, nt1;

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

        proc_node = mli_ray_octree_traversal_first_octree_node(t0, tm);

        do {
                switch (proc_node) {
                case 0: {
                        nt0 = mliVec_init(t0.x, t0.y, t0.z);
                        nt1 = mliVec_init(tm.x, tm.y, tm.z);
                        mli_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx].children[permutation],
                                octree->nodes[node_idx].types[permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node = mli_ray_octree_traversal_next_octree_node(
                                nt1, 4, 2, 1);
                        break;
                }
                case 1: {
                        nt0 = mliVec_init(t0.x, t0.y, tm.z);
                        nt1 = mliVec_init(tm.x, tm.y, t1.z);
                        mli_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[1 ^ permutation],
                                octree->nodes[node_idx].types[1 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node = mli_ray_octree_traversal_next_octree_node(
                                nt1, 5, 3, 8);
                        break;
                }
                case 2: {
                        nt0 = mliVec_init(t0.x, tm.y, t0.z);
                        nt1 = mliVec_init(tm.x, t1.y, tm.z);
                        mli_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[2 ^ permutation],
                                octree->nodes[node_idx].types[2 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node = mli_ray_octree_traversal_next_octree_node(
                                nt1, 6, 8, 3);
                        break;
                }
                case 3: {
                        nt0 = mliVec_init(t0.x, tm.y, tm.z);
                        nt1 = mliVec_init(tm.x, t1.y, t1.z);
                        mli_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[3 ^ permutation],
                                octree->nodes[node_idx].types[3 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node = mli_ray_octree_traversal_next_octree_node(
                                nt1, 7, 8, 8);
                        break;
                }
                case 4: {
                        nt0 = mliVec_init(tm.x, t0.y, t0.z);
                        nt1 = mliVec_init(t1.x, tm.y, tm.z);
                        mli_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[4 ^ permutation],
                                octree->nodes[node_idx].types[4 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node = mli_ray_octree_traversal_next_octree_node(
                                nt1, 8, 6, 5);
                        break;
                }
                case 5: {
                        nt0 = mliVec_init(tm.x, t0.y, tm.z);
                        nt1 = mliVec_init(t1.x, tm.y, t1.z);
                        mli_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[5 ^ permutation],
                                octree->nodes[node_idx].types[5 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node = mli_ray_octree_traversal_next_octree_node(
                                nt1, 8, 7, 8);
                        break;
                }
                case 6: {
                        nt0 = mliVec_init(tm.x, tm.y, t0.z);
                        nt1 = mliVec_init(t1.x, t1.y, tm.z);
                        mli_ray_octree_traversal_sub(
                                nt0,
                                nt1,
                                octree,
                                octree->nodes[node_idx]
                                        .children[6 ^ permutation],
                                octree->nodes[node_idx].types[6 ^ permutation],
                                permutation,
                                work,
                                work_on_leaf_node);
                        proc_node = mli_ray_octree_traversal_next_octree_node(
                                nt1, 8, 8, 7);
                        break;
                }
                case 7: {
                        nt0 = mliVec_init(tm.x, tm.y, tm.z);
                        nt1 = mliVec_init(t1.x, t1.y, t1.z);
                        mli_ray_octree_traversal_sub(
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

void mli_ray_octree_traversal(
        const struct mliOcTree *octree,
        const struct mliRay ray,
        void *work,
        void (*work_on_leaf_node)(
                void *,
                const struct mliOcTree *,
                const uint32_t))
{
        struct mliVec t0, t1;
        struct mliVec div;
        struct mliRay ray_wrt_octree;
        struct mliVec cube_upper, cube_size;
        struct mliCube cube;
        int32_t octree_root_node, octree_root_type;
        uint8_t permutation = 0;
        cube = octree->cube;
        cube_upper = mliCube_upper(cube);
        octree_root_node = 0u;
        octree_root_type = octree->root_type;
        cube_size = mliVec_add(cube.lower, cube_upper);

        ray_wrt_octree = ray;

        if (ray_wrt_octree.direction.x < 0) {
                ray_wrt_octree.support.x =
                        -ray_wrt_octree.support.x + cube_size.x;
                ray_wrt_octree.direction.x = -ray_wrt_octree.direction.x;
                permutation |= 4;
        } else if (ray_wrt_octree.direction.x == 0.0) {
                ray_wrt_octree.direction.x = MLI_RAY_OCTREE_TRAVERSAL_EPSILON;
        }

        if (ray_wrt_octree.direction.y < 0) {
                ray_wrt_octree.support.y =
                        -ray_wrt_octree.support.y + cube_size.y;
                ray_wrt_octree.direction.y = -ray_wrt_octree.direction.y;
                permutation |= 2;
        } else if (ray_wrt_octree.direction.y == 0.0) {
                ray_wrt_octree.direction.y = MLI_RAY_OCTREE_TRAVERSAL_EPSILON;
        }

        if (ray_wrt_octree.direction.z < 0) {
                ray_wrt_octree.support.z =
                        -ray_wrt_octree.support.z + cube_size.z;
                ray_wrt_octree.direction.z = -ray_wrt_octree.direction.z;
                permutation |= 1;
        } else if (ray_wrt_octree.direction.z == 0.0) {
                ray_wrt_octree.direction.z = MLI_RAY_OCTREE_TRAVERSAL_EPSILON;
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

        if (MLI_MAX3(t0.x, t0.y, t0.z) < MLI_MIN3(t1.x, t1.y, t1.z)) {
                mli_ray_octree_traversal_sub(
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

/* mli_ray_scenery_query */
/* --------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mli_inner_object_traversal(
        void *_inner,
        const struct mliOcTree *object_octree,
        const uint32_t object_octree_leaf_idx)
{
        /* traverse faces in an object-wavefront */
        struct mliQueryInnerWork *inner = (struct mliQueryInnerWork *)_inner;

        uint32_t f;
        const uint32_t num_faces_in_object_leaf = mliOcTree_leaf_num_objects(
                object_octree, object_octree_leaf_idx);

        struct mliIntersection tmp_isec = mliIntersection_init();

        for (f = 0; f < num_faces_in_object_leaf; f++) {

                uint32_t face_idx = mliOcTree_leaf_object_link(
                        object_octree, object_octree_leaf_idx, f);

                struct mliFace fv = inner->object->faces_vertices[face_idx];

                int32_t hit = mliRay_intersects_triangle(
                        inner->ray_object,
                        inner->object->vertices[fv.a],
                        inner->object->vertices[fv.b],
                        inner->object->vertices[fv.c],
                        &tmp_isec.distance_of_ray);

                tmp_isec.geometry_id.face = face_idx;

                if (hit) {
                        tmp_isec.position_local = mliRay_at(
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

int mli_query_object_reference(
        const struct mliObject *object,
        const struct mliOcTree *object_octree,
        const struct mliHomTraComp robject2root_comp,
        const struct mliRay ray_root,
        struct mliIntersection *isec)
{
        struct mliHomTra robject2root =
                mliHomTra_from_compact(robject2root_comp);

        struct mliQueryInnerWork inner;
        inner.has_intersection = 0;
        inner.intersection = isec;
        inner.ray_object = mliHomTra_ray_inverse(&robject2root, ray_root);
        inner.object = object;

        mli_ray_octree_traversal(
                object_octree,
                inner.ray_object,
                (void *)&inner,
                mli_inner_object_traversal);

        return inner.has_intersection;
}

void mli_outer_scenery_traversal(
        void *_outer,
        const struct mliOcTree *scenery_octree,
        const uint32_t scenery_octree_leaf_idx)
{
        /* traverse object-wavefronts in a scenery */
        struct mliQueryOuterWork *outer = (struct mliQueryOuterWork *)_outer;

        uint32_t ro;
        const uint32_t num_robjects_in_scenery_leaf =
                mliOcTree_leaf_num_objects(
                        scenery_octree, scenery_octree_leaf_idx);

        struct mliIntersection tmp_isec = mliIntersection_init();

        for (ro = 0; ro < num_robjects_in_scenery_leaf; ro++) {

                uint32_t robject_idx = mliOcTree_leaf_object_link(
                        scenery_octree, scenery_octree_leaf_idx, ro);
                uint32_t object_idx = outer->geometry->robjects[robject_idx];

                int32_t hit = mli_query_object_reference(
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

int mli_query_intersection(
        const struct mliScenery *scenery,
        const struct mliRay ray_root,
        struct mliIntersection *isec)
{
        struct mliQueryOuterWork outer;

        (*isec) = mliIntersection_init();

        outer.intersection = isec;
        outer.geometry = &scenery->geometry;
        outer.accelerator = &scenery->accelerator;
        outer.ray_root = ray_root;

        mli_ray_octree_traversal(
                &scenery->accelerator.scenery_octree,
                ray_root,
                (void *)&outer,
                mli_outer_scenery_traversal);

        if (isec->distance_of_ray == DBL_MAX) {
                return 0;
        } else {
                return 1;
        }
}

int mli_query_intersection_with_surface_normal(
        const struct mliScenery *scenery,
        const struct mliRay ray_root,
        struct mliIntersectionSurfaceNormal *isecsrf)
{
        struct mliIntersection isec = mliIntersection_init();

        const int has_intersection =
                mli_query_intersection(scenery, ray_root, &isec);

        if (has_intersection) {
                uint32_t robject_idx = isec.geometry_id.robj;
                uint32_t object_idx =
                        scenery->geometry.robjects[isec.geometry_id.robj];
                uint32_t face_idx = isec.geometry_id.face;

                struct mliHomTra robject2root = mliHomTra_from_compact(
                        scenery->geometry.robject2root[robject_idx]);
                struct mliRay ray_object =
                        mliHomTra_ray_inverse(&robject2root, ray_root);

                struct mliObject *obj = &scenery->geometry.objects[object_idx];

                struct mliFace fv = obj->faces_vertices[face_idx];
                struct mliFace fvn = obj->faces_vertex_normals[face_idx];

                (*isecsrf) = mliIntersectionSurfaceNormal_init();
                isecsrf->distance_of_ray = isec.distance_of_ray;
                isecsrf->geometry_id = isec.geometry_id;
                isecsrf->position = mliRay_at(&ray_root, isec.distance_of_ray);
                isecsrf->position_local = isec.position_local;

                /* find surface-normal */
                isecsrf->surface_normal_local = mliTriangle_surface_normal(
                        obj->vertices[fv.a],
                        obj->vertices[fv.b],
                        obj->vertices[fv.c],
                        obj->vertex_normals[fvn.a],
                        obj->vertex_normals[fvn.b],
                        obj->vertex_normals[fvn.c],
                        isecsrf->position_local);

                isecsrf->surface_normal = mliHomTra_dir(
                        &robject2root, isecsrf->surface_normal_local);

                isecsrf->from_outside_to_inside =
                        mli_ray_runs_from_outside_to_inside(
                                ray_object.direction,
                                isecsrf->surface_normal_local);

                return 1;
        } else {
                return 0;
        }
}

/* mli_triangle_intersection */
/* ------------------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

int mliRay_intersects_triangle(
        const struct mliRay ray,
        const struct mliVec vertex_a,
        const struct mliVec vertex_b,
        const struct mliVec vertex_c,
        double *intersection_ray_parameter)
{
        /* Moeller-Trumbore-intersection-algorithm */
        struct mliVec edge1;
        struct mliVec edge2;
        struct mliVec h, s, q;
        double a, f, u, v, t;
        edge1 = mliVec_substract(vertex_b, vertex_a);
        edge2 = mliVec_substract(vertex_c, vertex_a);
        h = mliVec_cross(ray.direction, edge2);
        a = mliVec_dot(edge1, h);

        if (a > -MLI_EPSILON && a < MLI_EPSILON)
                return 0; /* This ray is parallel to this triangle. */
        f = 1.0 / a;
        s = mliVec_substract(ray.support, vertex_a);
        u = f * mliVec_dot(s, h);
        if (u < 0.0 || u > 1.0)
                return 0;
        q = mliVec_cross(s, edge1);
        v = f * mliVec_dot(ray.direction, q);
        if (v < 0.0 || u + v > 1.0)
                return 0;
        /* At this stage we can compute t to find out where the intersection */
        /* point is on the line. */
        t = f * mliVec_dot(edge2, q);
        if (t > MLI_EPSILON) {
                (*intersection_ray_parameter) = t;
                return 1;
        } else {
                /* This means that there is a line intersection but not a */
                /* ray intersection. */
                return 0;
        }
}

struct mliVec mliTriangle_interpolate_surface_normal(
        const struct mliVec vertex_normal_a,
        const struct mliVec vertex_normal_b,
        const struct mliVec vertex_normal_c,
        const struct mliBarycentrigWeights weights)
{
        return mliVec_init(
                vertex_normal_a.x * weights.a + vertex_normal_b.x * weights.b +
                        vertex_normal_c.x * weights.c,

                vertex_normal_a.y * weights.a + vertex_normal_b.y * weights.b +
                        vertex_normal_c.y * weights.c,

                vertex_normal_a.z * weights.a + vertex_normal_b.z * weights.b +
                        vertex_normal_c.z * weights.c);
}

struct mliVec mliTriangle_surface_normal(
        const struct mliVec vertex_a,
        const struct mliVec vertex_b,
        const struct mliVec vertex_c,
        const struct mliVec vertex_normal_a,
        const struct mliVec vertex_normal_b,
        const struct mliVec vertex_normal_c,
        const struct mliVec intersection_position)
{
        struct mliVec surface_normal;
        struct mliBarycentrigWeights normal_weights = mli_barycentric_weights(
                vertex_a, vertex_b, vertex_c, intersection_position);

        surface_normal = mliTriangle_interpolate_surface_normal(
                vertex_normal_a,
                vertex_normal_b,
                vertex_normal_c,
                normal_weights);

        surface_normal = mliVec_normalized(surface_normal);

        return surface_normal;
}

/* mli_version */
/* ----------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

void mli_logo_fprint(FILE *f)
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

void mli_authors_and_affiliations_fprint(FILE *f)
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

/* mli_viewer_Config */
/* ----------------- */

/* Copyright 2018-2020 Sebastian Achim Mueller */

struct mlivrConfig mlivrConfig_default(void)
{
        struct mlivrConfig cfg;
        cfg.random_seed = 0u;

        cfg.preview_num_cols = 160u;
        cfg.preview_num_rows = 90u / 2;

        cfg.export_num_cols = 1280u;
        cfg.export_num_rows = 720u;

        cfg.step_length = 1.0;

        cfg.view.position.x = 0.;
        cfg.view.position.y = 0.;
        cfg.view.position.z = 0.;

        cfg.view.rotation.x = mli_deg2rad(90.);
        cfg.view.rotation.y = 0.;
        cfg.view.rotation.z = 0.;

        cfg.view.field_of_view = mli_deg2rad(80.);

        cfg.aperture_camera_f_stop_ratio = 0.95;
        cfg.aperture_camera_image_sensor_width = 64e-3;
        return cfg;
}

/* mli_viewer_Cursor */
/* ----------------- */

/* Copyright 2019 Sebastian Achim Mueller                                     */

void mlivrCursor_move_up(struct mlivrCursor *cursor)
{
        if (cursor->row != 0)
                cursor->row -= 1;
}

void mlivrCursor_move_down(struct mlivrCursor *cursor)
{
        if (cursor->row + 1u < cursor->num_rows)
                cursor->row += 1;
}

void mlivrCursor_move_right(struct mlivrCursor *cursor)
{
        if (cursor->col != 0)
                cursor->col -= 1;
}

void mlivrCursor_move_left(struct mlivrCursor *cursor)
{
        if (cursor->col + 1u < cursor->num_cols)
                cursor->col += 1;
}

/* mli_viewer_toggle_stdin */
/* ----------------------- */


#ifdef __unix__


struct termios mlivr_non_canonical_stdin(void)
{
        struct termios old_terminal;
        struct termios new_terminal;
        tcgetattr(STDIN_FILENO, &old_terminal);
        new_terminal = old_terminal;
        new_terminal.c_lflag &= ~(ICANON);
        tcsetattr(STDIN_FILENO, TCSANOW, &new_terminal);
        return old_terminal;
}

void mlivr_restore_stdin(struct termios *old_terminal)
{
        tcsetattr(STDIN_FILENO, TCSANOW, old_terminal);
}

#endif /* __unix__ */

/* mli_viewer_viewer */
/* ----------------- */

/* Copyright 2019 Sebastian Achim Mueller                                     */

void mlivr_clear_screen(void)
{
        uint64_t n = 20;
        while (n) {
                putchar('\n');
                n--;
        }
}

void mlivr_print_help(void)
{
        mlivr_clear_screen();
        mli_logo_fprint(stdout);
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
        printf("    increace          [  m  ]     super sampling    [  b  ]\n");
        printf("    decreace          [  n  ]     color/monochrome  [  g  ]\n");
        printf("\n");
        printf("  Atmosphere                    Sun\n");
        printf("    on/off            [  0  ]     later daytime     [  9  ]\n");
        printf("    - altitude        [  4  ]     earlier daytime   [  8  ]\n");
        printf("    + altitude        [  5  ]     + latitude        [  7  ]\n");
        printf("                                  - latitude        [  6  ]\n");
        printf("\n");
        mli_authors_and_affiliations_fprint(stdout);
}

void mlivr_print_info_line(
        const struct mliView view,
        const struct mlivrCursor cursor,
        const struct mliTracerConfig tracer_config)
{
        printf("Help 'h', "
               "Cam: "
               "pos[% -.2e, % -.2e, % -.2e]m, "
               "rot[% -.1f, % -.1f, % -.1f]deg, "
               "fov %.2fdeg, ",
               view.position.x,
               view.position.y,
               view.position.z,
               mli_rad2deg(view.rotation.x),
               mli_rad2deg(view.rotation.y),
               mli_rad2deg(view.rotation.z),
               mli_rad2deg(view.field_of_view));
        printf("Sun: lat % 3.0fdeg, %02d:%02dh, alt % 3.1fkm",
               mli_rad2deg(tracer_config.atmosphere.sunLatitude),
               (int)(tracer_config.atmosphere.sunHourAngle),
               (int)(tracer_config.atmosphere.sunHourAngle * 60) % 60,
               tracer_config.atmosphere.altitude * 1e-3);
        if (cursor.active) {
                printf(", Cursor[%3ld, %3ld]pix", cursor.col, cursor.row);
        }
        printf(".\n");
}

void mlivr_timestamp_now_19chars(char *buffer)
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

int mlivr_get_key(void)
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

int mlivr_export_image(
        const struct mliScenery *scenery,
        const struct mlivrConfig config,
        const struct mliView view,
        struct mliPrng *prng,
        const struct mliTracerConfig *tracer_config,
        const double object_distance,
        const char *path)
{
        struct mliImage full = mliImage_init();
        struct mliHomTraComp camera2root_comp;
        struct mliApertureCamera apcam = mliApertureCamera_init();

        const double image_ratio =
                ((double)config.export_num_cols /
                 (double)config.export_num_rows);
        chk_mem(mliImage_malloc(
                &full, config.export_num_cols, config.export_num_rows));
        camera2root_comp = mliView_to_HomTraComp(view);
        apcam.focal_length =
                mliApertureCamera_focal_length_given_field_of_view_and_sensor_width(
                        view.field_of_view,
                        config.aperture_camera_image_sensor_width);
        apcam.aperture_radius = 0.5 * (apcam.focal_length /
                                       config.aperture_camera_f_stop_ratio);
        apcam.image_sensor_distance =
                mli_thin_lens_get_image_given_focal_and_object(
                        apcam.focal_length, object_distance);
        apcam.image_sensor_width_x = config.aperture_camera_image_sensor_width;
        apcam.image_sensor_width_y = apcam.image_sensor_width_x / image_ratio;
        mliApertureCamera_render_image(
                apcam, camera2root_comp, scenery, &full, tracer_config, prng);
        chk_msg(mliImage_write_to_path(&full, path), "Failed to write ppm.");
        mliImage_free(&full);
        return 1;
error:
        return 0;
}

int mlivr_run_interactive_viewer_try_non_canonical_stdin(
        const struct mliScenery *scenery,
        const struct mlivrConfig config)
{
#ifdef HAVE_TERMIOS_H
        struct termios old_terminal = mlivr_non_canonical_stdin();
#endif
        int rc = mlivr_run_interactive_viewer(scenery, config);

#ifdef HAVE_TERMIOS_H
        mlivr_restore_stdin(&old_terminal);
#endif
        return rc;
}

int mlivr_run_interactive_viewer(
        const struct mliScenery *scenery,
        const struct mlivrConfig config)
{
        struct mliPrng prng = mliPrng_init_MT19937(config.random_seed);
        struct mliTracerConfig tracer_config = mliTracerConfig_init();
        char path[1024];
        int key;
        int super_resolution = 0;
        struct mlivrCursor cursor;
        uint64_t num_screenshots = 0;
        uint64_t print_mode = MLI_ASCII_MONOCHROME;
        char timestamp[20];
        struct mliView view = config.view;
        struct mliImage img = mliImage_init();
        struct mliImage img2 = mliImage_init();
        const double row_over_column_pixel_ratio = 2.0;
        int update_image = 1;
        int print_help = 0;
        int print_scenery_info = 0;
        int has_probing_intersection = 0;
        struct mliIntersectionSurfaceNormal probing_intersection;

        mlivr_timestamp_now_19chars(timestamp);
        chk_mem(mliImage_malloc(
                &img, config.preview_num_cols, config.preview_num_rows));
        chk_mem(mliImage_malloc(
                &img2,
                config.preview_num_cols * 2u,
                config.preview_num_rows * 2u));

        cursor.active = 0;
        cursor.col = config.preview_num_cols / 2;
        cursor.row = config.preview_num_rows / 2;
        cursor.num_cols = config.preview_num_cols;
        cursor.num_rows = config.preview_num_rows;
        goto show_image;

        while ((key = mlivr_get_key()) != MLIVR_ESCAPE_KEY) {
                update_image = 1;
                print_help = 0;
                print_scenery_info = 0;
                if (cursor.active) {
                        update_image = 0;
                        super_resolution = 0;
                        switch (key) {
                        case 'i':
                                mlivrCursor_move_up(&cursor);
                                break;
                        case 'k':
                                mlivrCursor_move_down(&cursor);
                                break;
                        case 'l':
                                mlivrCursor_move_left(&cursor);
                                break;
                        case 'j':
                                mlivrCursor_move_right(&cursor);
                                break;
                        case 'c':
                                cursor.active = !cursor.active;
                                break;
                        case 'h':
                                print_help = 1;
                                break;
                        case MLIVR_SPACE_KEY:
                                sprintf(path,
                                        "%s_%06lu.ppm",
                                        timestamp,
                                        num_screenshots);
                                num_screenshots++;
                                chk(mlivr_export_image(
                                        scenery,
                                        config,
                                        view,
                                        &prng,
                                        &tracer_config,
                                        probing_intersection.distance_of_ray,
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
                                view = mliView_move_forward(
                                        view, config.step_length);
                                break;
                        case 's':
                                view = mliView_move_forward(
                                        view, -config.step_length);
                                break;
                        case 'a':
                                view = mliView_move_right(
                                        view, -config.step_length);
                                break;
                        case 'd':
                                view = mliView_move_right(
                                        view, config.step_length);
                                break;
                        case 'q':
                                view = mliView_move_up(
                                        view, config.step_length);
                                break;
                        case 'e':
                                view = mliView_move_up(
                                        view, -config.step_length);
                                break;
                        case 'i':
                                view = mliView_look_up_when_possible(view, .05);
                                break;
                        case 'k':
                                view = mliView_look_down_when_possible(
                                        view, .05);
                                break;
                        case 'l':
                                view = mliView_look_right(view, -.05);
                                break;
                        case 'j':
                                view = mliView_look_right(view, .05);
                                break;
                        case 'n':
                                view = mliView_decrease_fov(view, 1.05);
                                break;
                        case 'm':
                                view = mliView_increase_fov(view, 1.05);
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
                        case MLIVR_SPACE_KEY:
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
                                mliAtmosphere_decrease_altitude(
                                        &tracer_config.atmosphere, 0.9);
                                break;
                        case '5':
                                mliAtmosphere_increase_altitude(
                                        &tracer_config.atmosphere, 1.1);
                                break;
                        case '6':
                                mliAtmosphere_decrease_latitude(
                                        &tracer_config.atmosphere,
                                        mli_deg2rad(2.0));
                                break;
                        case '7':
                                mliAtmosphere_increase_latitude(
                                        &tracer_config.atmosphere,
                                        mli_deg2rad(2.0));
                                break;
                        case '8':
                                mliAtmosphere_decrease_hours(
                                        &tracer_config.atmosphere, 0.1);
                                break;
                        case '9':
                                mliAtmosphere_increase_hours(
                                        &tracer_config.atmosphere, 0.1);
                                break;
                        case '0':
                                tracer_config.have_atmosphere =
                                        !tracer_config.have_atmosphere;
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
                                mliPinHoleCamera_render_image_with_view(
                                        view,
                                        scenery,
                                        &img2,
                                        row_over_column_pixel_ratio,
                                        &tracer_config,
                                        &prng);
                                mliImage_scale_down_twice(&img2, &img);
                        } else {
                                mliPinHoleCamera_render_image_with_view(
                                        view,
                                        scenery,
                                        &img,
                                        row_over_column_pixel_ratio,
                                        &tracer_config,
                                        &prng);
                        }
                }
                mlivr_clear_screen();
                if (cursor.active) {
                        char symbols[1];
                        uint64_t rows[1];
                        uint64_t cols[1];
                        const uint64_t num_symbols = 1u;
                        symbols[0] = 'X';
                        rows[0] = cursor.row;
                        cols[0] = cursor.col;
                        mliImage_print_chars(
                                &img,
                                symbols,
                                rows,
                                cols,
                                num_symbols,
                                print_mode);
                        {
                                struct mliPinHoleCamera pin_hole_camera =
                                        mliPinHoleCamera_init(
                                                view.field_of_view,
                                                &img,
                                                row_over_column_pixel_ratio);

                                struct mliHomTraComp camera2root_comp =
                                        mliView_to_HomTraComp(view);
                                struct mliHomTra camera2root =
                                        mliHomTra_from_compact(
                                                camera2root_comp);
                                struct mliRay probing_ray_wrt_camera;
                                struct mliRay probing_ray_wrt_root;

                                probing_ray_wrt_camera =
                                        mliPinHoleCamera_ray_at_row_col(
                                                &pin_hole_camera,
                                                &img,
                                                cursor.row,
                                                cursor.col);
                                probing_ray_wrt_root = mliHomTra_ray(
                                        &camera2root, probing_ray_wrt_camera);

                                probing_intersection =
                                        mliIntersectionSurfaceNormal_init();

                                has_probing_intersection =
                                        mli_query_intersection_with_surface_normal(
                                                scenery,
                                                probing_ray_wrt_root,
                                                &probing_intersection);
                        }
                } else {
                        mliImage_print(&img, print_mode);
                }
                mlivr_print_info_line(view, cursor, tracer_config);
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
                        mlivr_print_help();
                }
                if (print_scenery_info) {
                        mlivr_clear_screen();
                        mliScenery_info_fprint(stdout, scenery);
                }
        }

        mliImage_free(&img);
        mliImage_free(&img2);
        return 1;
error:
        mliImage_free(&img);
        mliImage_free(&img2);
        return 0;
}

