/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_srot_avx2.c
 *  @brief Plane rotations in AVX2.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* Application of 2x2 Plane Rotation on two vectors */
int fla_srot_avx2(aocl_int64_t *n, real *dx, aocl_int64_t *incx, real *dy, aocl_int64_t *incy,
                  real *c__, real *s)
{
    aocl_int64_t i__1;

    aocl_int64_t i__;
    aocl_int64_t ix, iy;

    real stemp;

    __m256 vs8_c, vs8_s;
    __m256 vs8_idx0, vs8_idx1, vs8_odx0, vs8_odx1;
    __m256 vs8_idy0, vs8_idy1, vs8_ody0, vs8_ody1;

    __m128 vs4_c, vs4_s;
    __m128 vs4_idx0, vs4_idy0;
    __m128 vs4_odx0, vs4_ody0;

    /* Parameter adjustments */
    --dy;
    --dx;

    /* Function Body */
    if(*n <= 0)
    {
        return 0;
    }
    i__1 = *n;
    i__ = 1;
    if(*incx == 1 && *incy == 1)
    {
        goto L20;
    }
    /*       code for unequal increments or equal increments not equal */
    /*         to 1 */
    ix = 1;
    iy = 1;
    if(*incx < 0)
    {
        ix = (-(*n) + 1) * *incx + 1;
    }
    if(*incy < 0)
    {
        iy = (-(*n) + 1) * *incy + 1;
    }
    for(; i__ <= i__1; i__++)
    {
        stemp = *c__ * dx[ix] + *s * dy[iy];
        dy[iy] = *c__ * dy[iy] - *s * dx[ix];
        dx[ix] = stemp;
        ix += *incx;
        iy += *incy;
    }
    return 0;

    /*       code for both increments equal to 1 */

L20:
    vs4_c = _mm_load1_ps((real const *)c__);
    vs4_s = _mm_load1_ps((real const *)s);
    vs8_c = _mm256_broadcastss_ps(vs4_c);
    vs8_s = _mm256_broadcastss_ps(vs4_s);
    if(i__1 >= 0x10)
    {
        for(; i__ <= (i__1 - 15); i__ += 16)
        {
            /* load input vectors */
            vs8_idx0 = _mm256_loadu_ps((real const *)&dx[i__]);
            vs8_idx1 = _mm256_loadu_ps((real const *)&dx[i__ + 8]);
            vs8_idy0 = _mm256_loadu_ps((real const *)&dy[i__]);
            vs8_idy1 = _mm256_loadu_ps((real const *)&dy[i__ + 8]);

            /* apply the plane rotation matrix  */
            vs8_odx0 = _mm256_mul_ps(vs8_c, vs8_idx0);
            vs8_odx1 = _mm256_mul_ps(vs8_c, vs8_idx1);
            vs8_ody0 = _mm256_mul_ps(vs8_s, vs8_idx0);
            vs8_ody1 = _mm256_mul_ps(vs8_s, vs8_idx1);

            vs8_odx0 = _mm256_fmadd_ps(vs8_s, vs8_idy0, vs8_odx0);
            vs8_odx1 = _mm256_fmadd_ps(vs8_s, vs8_idy1, vs8_odx1);
            vs8_ody0 = _mm256_fmsub_ps(vs8_c, vs8_idy0, vs8_ody0);
            vs8_ody1 = _mm256_fmsub_ps(vs8_c, vs8_idy1, vs8_ody1);

            /* store the outputs */
            _mm256_storeu_ps((real *)&dx[i__], vs8_odx0);
            _mm256_storeu_ps((real *)&dx[i__ + 8], vs8_odx1);
            _mm256_storeu_ps((real *)&dy[i__], vs8_ody0);
            _mm256_storeu_ps((real *)&dy[i__ + 8], vs8_ody1);
        }
    }
    if(i__1 & 0x08)
    {
        /* load input vectors */
        vs8_idx0 = _mm256_loadu_ps((real const *)&dx[i__]);
        vs8_idy0 = _mm256_loadu_ps((real const *)&dy[i__]);

        /* apply the plane rotation matrix  */
        vs8_odx0 = _mm256_mul_ps(vs8_c, vs8_idx0);
        vs8_ody0 = _mm256_mul_ps(vs8_s, vs8_idx0);

        vs8_odx0 = _mm256_fmadd_ps(vs8_s, vs8_idy0, vs8_odx0);
        vs8_ody0 = _mm256_fmsub_ps(vs8_c, vs8_idy0, vs8_ody0);

        /* store the outputs */
        _mm256_storeu_ps((real *)&dx[i__], vs8_odx0);
        _mm256_storeu_ps((real *)&dy[i__], vs8_ody0);

        i__ += 8;
    }
    if(i__1 & 0x04)
    {
        /* load input vectors */
        vs4_idx0 = _mm_loadu_ps((real const *)&dx[i__]);
        vs4_idy0 = _mm_loadu_ps((real const *)&dy[i__]);

        /* apply the plane rotation matrix  */
        vs4_odx0 = _mm_mul_ps(vs4_c, vs4_idx0);
        vs4_ody0 = _mm_mul_ps(vs4_s, vs4_idx0);

        vs4_odx0 = _mm_fmadd_ps(vs4_s, vs4_idy0, vs4_odx0);
        vs4_ody0 = _mm_fmsub_ps(vs4_c, vs4_idy0, vs4_ody0);

        /* store the outputs */
        _mm_storeu_ps((real *)&dx[i__], vs4_odx0);
        _mm_storeu_ps((real *)&dy[i__], vs4_ody0);

        i__ += 4;
    }
    /* remaining 1 to 3 elements */
    for(; i__ <= i__1; i__++)
    {
        stemp = *c__ * dx[i__] + *s * dy[i__];
        dy[i__] = *c__ * dy[i__] - *s * dx[i__];
        dx[i__] = stemp;
    }
    return 0;
}
#endif
