/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_drot_avx2.c
 *  @brief Plane rotations in AVX2.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* Application of 2x2 Plane Rotation on two vectors */
int fla_drot_avx2(integer *n, doublereal *dx, integer *incx, doublereal *dy, integer *incy, doublereal *c__, doublereal *s)
{
    integer i__1;

    integer i__;
    integer ix, iy;

    __m256d vd4_c, vd4_s;
    __m256d vd4_idx0, vd4_idx1, vd4_odx0, vd4_odx1;
    __m256d vd4_idy0, vd4_idy1, vd4_ody0, vd4_ody1;

    __m128d vd2_c, vd2_s;
    __m128d vd2_cs, vd2_sc;
    __m128d vd2_idx0, vd2_idx1;
    __m128d vd2_idy0, vd2_idy1;
    __m128d vd2_odx0, vd2_ody0;

    vd2_c = _mm_loaddup_pd((double const *) c__);
    vd2_s = _mm_loaddup_pd((double const *) s);

    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if (*n <= 0)
    {
        return 0;
    }
    i__1 = *n;
    i__ = 1;
    if (*incx == 1 && *incy == 1)
    {
        goto L20;
    }
/*       code for unequal increments or equal increments not equal */
/*         to 1 */
    ix = 1;
    iy = 1;
    if (*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if (*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }

    vd2_cs = _mm_unpackhi_pd(vd2_s, vd2_c);
    vd2_sc = _mm_unpackhi_pd(vd2_c, vd2_s);
    if (i__1 >= 0x02) /* 2 iterations at once using SIMD */
    {
        for ( ; i__ <= (i__1 - 1); i__ += 2)
        {
            /* load input vectors */
            vd2_idx0 = _mm_loaddup_pd((double const *) &dx[ix]);
            vd2_idy0 = _mm_loaddup_pd((double const *) &dy[iy]);
            vd2_idx1 = _mm_loaddup_pd((double const *) &dx[ix + *incx]);
            vd2_idy1 = _mm_loaddup_pd((double const *) &dy[iy + *incy]);

            /* apply the plane rotation matrix  */
            vd2_idx0 = _mm_mul_pd(vd2_sc, vd2_idx0);
            vd2_idx0 = _mm_fmsubadd_pd(vd2_cs, vd2_idy0, vd2_idx0);
            vd2_idx1 = _mm_mul_pd(vd2_sc, vd2_idx1);
            vd2_idx1 = _mm_fmsubadd_pd(vd2_cs, vd2_idy1, vd2_idx1);

            /* store the outputs */
            _mm_storel_pd((double *) &dx[ix], vd2_idx0);
            _mm_storeh_pd((double *) &dy[iy], vd2_idx0);
            _mm_storel_pd((double *) &dx[ix + *incx], vd2_idx1);
            _mm_storeh_pd((double *) &dy[iy + *incy], vd2_idx1);

            ix += 2 * *incx;
            iy += 2 * *incy;
        }
    }
    if (i__1 & 0x01) /* last iteration */
    {
        /* load input vectors */
        vd2_idx0 = _mm_loaddup_pd((double const *) &dx[ix]);
        vd2_idy0 = _mm_loaddup_pd((double const *) &dy[iy]);

        /* apply the plane rotation matrix  */
        vd2_odx0 = _mm_mul_pd(vd2_sc, vd2_idx0);
        vd2_ody0 = _mm_fmsubadd_pd(vd2_cs, vd2_idy0, vd2_odx0);

        /* store the outputs */
        _mm_storel_pd((double *) &dx[ix], vd2_ody0);
        _mm_storeh_pd((double *) &dy[iy], vd2_ody0);
    }
    return 0;

/*       code for both increments equal to 1 */

L20:
    vd4_c = _mm256_broadcastsd_pd(vd2_c);
    vd4_s = _mm256_broadcastsd_pd(vd2_s);
    if (i__1 >= 0x08)
    {
        for ( ; i__ <= (i__1 - 7); i__ += 8)
        {
            /* load input vectors */
            vd4_idx0 = _mm256_loadu_pd((double const *) &dx[i__]);
            vd4_idx1 = _mm256_loadu_pd((double const *) &dx[i__ + 4]);
            vd4_idy0 = _mm256_loadu_pd((double const *) &dy[i__]);
            vd4_idy1 = _mm256_loadu_pd((double const *) &dy[i__ + 4]);

            /* apply the plane rotation matrix  */
            vd4_odx0 = _mm256_mul_pd(vd4_c, vd4_idx0);
            vd4_odx1 = _mm256_mul_pd(vd4_c, vd4_idx1);
            vd4_ody0 = _mm256_mul_pd(vd4_s, vd4_idx0);
            vd4_ody1 = _mm256_mul_pd(vd4_s, vd4_idx1);

            vd4_odx0 = _mm256_fmadd_pd(vd4_s, vd4_idy0, vd4_odx0);
            vd4_odx1 = _mm256_fmadd_pd(vd4_s, vd4_idy1, vd4_odx1);
            vd4_ody0 = _mm256_fmsub_pd(vd4_c, vd4_idy0, vd4_ody0);
            vd4_ody1 = _mm256_fmsub_pd(vd4_c, vd4_idy1, vd4_ody1);

            /* store the outputs */
            _mm256_storeu_pd((double *) &dx[i__], vd4_odx0);
            _mm256_storeu_pd((double *) &dx[i__ + 4], vd4_odx1);
            _mm256_storeu_pd((double *) &dy[i__], vd4_ody0);
            _mm256_storeu_pd((double *) &dy[i__ + 4], vd4_ody1);
        }
    }
    if (i__1 & 0x04)
    {
        /* load input vectors */
        vd4_idx0 = _mm256_loadu_pd((double const *) &dx[i__]);
        vd4_idy0 = _mm256_loadu_pd((double const *) &dy[i__]);

        /* apply the plane rotation matrix  */
        vd4_odx0 = _mm256_mul_pd(vd4_c, vd4_idx0);
        vd4_ody0 = _mm256_mul_pd(vd4_s, vd4_idx0);

        vd4_odx0 = _mm256_fmadd_pd(vd4_s, vd4_idy0, vd4_odx0);
        vd4_ody0 = _mm256_fmsub_pd(vd4_c, vd4_idy0, vd4_ody0);

        /* store the outputs */
        _mm256_storeu_pd((double *) &dx[i__], vd4_odx0);
        _mm256_storeu_pd((double *) &dy[i__], vd4_ody0);

        i__ += 4;
    }
    if (i__1 & 0x02)
    {
        /* load input vectors */
        vd2_idx0 = _mm_loadu_pd((double const *) &dx[i__]);
        vd2_idy0 = _mm_loadu_pd((double const *) &dy[i__]);

        /* apply the plane rotation matrix  */
        vd2_odx0 = _mm_mul_pd(vd2_c, vd2_idx0);
        vd2_ody0 = _mm_mul_pd(vd2_s, vd2_idx0);

        vd2_odx0 = _mm_fmadd_pd(vd2_s, vd2_idy0, vd2_odx0);
        vd2_ody0 = _mm_fmsub_pd(vd2_c, vd2_idy0, vd2_ody0);

        /* store the outputs */
        _mm_storeu_pd((double *) &dx[i__], vd2_odx0);
        _mm_storeu_pd((double *) &dy[i__], vd2_ody0);

        i__ += 2;
    }
    if (i__1 & 0x01)
    {
        /* load input vectors */
        vd2_idx0 = _mm_loaddup_pd((double const *) &dx[i__]);
        vd2_idy0 = _mm_loaddup_pd((double const *) &dy[i__]);

        /* apply the plane rotation matrix  */
        vd2_odx0 = _mm_mul_pd(vd2_c, vd2_idx0);
        vd2_ody0 = _mm_mul_pd(vd2_s, vd2_idx0);

        vd2_odx0 = _mm_fmadd_pd(vd2_s, vd2_idy0, vd2_odx0);
        vd2_ody0 = _mm_fmsub_pd(vd2_c, vd2_idy0, vd2_ody0);

        /* store the outputs */
        _mm_store_sd((double *) &dx[i__], vd2_odx0);
        _mm_store_sd((double *) &dy[i__], vd2_ody0);
    }
    return 0;
}
#endif
