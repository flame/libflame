/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_dhrot3_avx2.c
 *  @brief 3x3 Householder in AVX2.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* Application of 3x3 Householder reflector on a 3xn matrix */
int fla_dhrot3_avx2(integer *n, doublereal *a, integer *lda, doublereal *v, doublereal *tau)
{
    integer ir;
    doublereal mtau = -*tau;
    __m256d vd4_a0, vd4_a1, vd4_a2;
    __m256d vd4_v1, vd4_v2;
    __m256d vd4_p1, vd4_p2;
    __m256d vd4_tmp, vd4_tau;

    __m128d vd2_a0, vd2_a1, vd2_a2;
    __m128d vd2_v1, vd2_v2;
    __m128d vd2_p1, vd2_p2;
    __m128d vd2_tmp, vd2_tau;

    vd2_v1  = _mm_loaddup_pd((double const *) &v[1]);
    vd2_v2  = _mm_loaddup_pd((double const *) &v[2]);
    vd2_tau = _mm_loaddup_pd((double const *) &mtau);

    ir = 1;
    if (*n >= 4) /* 4 iterations at once through SIMD */
    {
        vd4_v1  = _mm256_broadcastsd_pd(vd2_v1);
        vd4_v2  = _mm256_broadcastsd_pd(vd2_v2);
        vd4_tau = _mm256_broadcastsd_pd(vd2_tau);

        for ( ; ir <= (*n - 3); ir += 4)
        {
            /* load column elements a0, a1, a2 */
            vd4_a0 = _mm256_loadu_pd((double const *) &a[ir + 1 * *lda]);
            vd4_a1 = _mm256_loadu_pd((double const *) &a[ir + 2 * *lda]);
            vd4_a2 = _mm256_loadu_pd((double const *) &a[ir + 3 * *lda]);

            vd4_p1 = _mm256_mul_pd(vd4_a1, vd4_v1); /* a1 * v1 */
            vd4_p2 = _mm256_mul_pd(vd4_a2, vd4_v2); /* a2 * v2 */

            vd4_p1 = _mm256_add_pd(vd4_p1, vd4_a0); /* a0 + a1 * v1 */
            vd4_p1 = _mm256_add_pd(vd4_p1, vd4_p2); /*  " + a2 * v2 */

            vd4_tmp = _mm256_mul_pd(vd4_p1, vd4_tau); /* tmp = " * tau */
            vd4_p1 = _mm256_mul_pd(vd4_tmp, vd4_v1); /* - tmp * v1 */
            vd4_p2 = _mm256_mul_pd(vd4_tmp, vd4_v2); /* - tmp * v2 */

            vd4_a0 = _mm256_add_pd(vd4_tmp, vd4_a0); /* a0 -= tmp */
            vd4_a1 = _mm256_add_pd(vd4_p1, vd4_a1); /* a1 -= tmp * v1 */
            vd4_a2 = _mm256_add_pd(vd4_p2, vd4_a2); /* a2 -= tmp * v2 */

            /* store the column elements */
            _mm256_storeu_pd((double *) &a[ir + 1 * *lda], vd4_a0);
            _mm256_storeu_pd((double *) &a[ir + 2 * *lda], vd4_a1);
            _mm256_storeu_pd((double *) &a[ir + 3 * *lda], vd4_a2);
        }
    }
    if (*n & 0x02) /* 2 remaining iterations */
    {
        /* load column elements a0, a1, a2 */
        vd2_a0 = _mm_loadu_pd((double const *) &a[ir + 1 * *lda]);
        vd2_a1 = _mm_loadu_pd((double const *) &a[ir + 2 * *lda]);
        vd2_a2 = _mm_loadu_pd((double const *) &a[ir + 3 * *lda]);

        vd2_p1 = _mm_mul_pd(vd2_a1, vd2_v1); /* a1 * v1 */
        vd2_p2 = _mm_mul_pd(vd2_a2, vd2_v2); /* a2 * v2 */

        vd2_p1 = _mm_add_pd(vd2_p1, vd2_a0); /* a0 + a1 * v1 */
        vd2_p1 = _mm_add_pd(vd2_p1, vd2_p2); /*  " + a2 * v2 */

        vd2_tmp = _mm_mul_pd(vd2_p1, vd2_tau); /* tmp = " * tau */
        vd2_p1 = _mm_mul_pd(vd2_tmp, vd2_v1); /* - tmp * v1 */
        vd2_p2 = _mm_mul_pd(vd2_tmp, vd2_v2); /* - tmp * v2 */

        vd2_a0 = _mm_add_pd(vd2_tmp, vd2_a0); /* a0 -= tmp */
        vd2_a1 = _mm_add_pd(vd2_p1, vd2_a1); /* a1 -= tmp * v1 */
        vd2_a2 = _mm_add_pd(vd2_p2, vd2_a2); /* a2 -= tmp * v2 */

        /* store the column elements */
        _mm_storeu_pd((double *) &a[ir + 1 * *lda], vd2_a0);
        _mm_storeu_pd((double *) &a[ir + 2 * *lda], vd2_a1);
        _mm_storeu_pd((double *) &a[ir + 3 * *lda], vd2_a2);

        ir += 2;
    }
    if (*n & 0x01) /* last iteration */
    {
        vd2_a0 = _mm_loaddup_pd((double const *) &a[ir + 1 * *lda]);
        vd2_a1 = _mm_loaddup_pd((double const *) &a[ir + 2 * *lda]);
        vd2_a2 = _mm_loaddup_pd((double const *) &a[ir + 3 * *lda]);

        vd2_p1 = _mm_mul_pd(vd2_a1, vd2_v1); /* a1 * v1 */
        vd2_p2 = _mm_mul_pd(vd2_a2, vd2_v2); /* a2 * v2 */

        vd2_p1 = _mm_add_pd(vd2_p1, vd2_a0); /* a0 + a1 * v1 */
        vd2_p1 = _mm_add_pd(vd2_p1, vd2_p2); /*  " + a2 * v2 */

        vd2_tmp = _mm_mul_pd(vd2_p1, vd2_tau); /* tmp = " * tau */
        vd2_p1 = _mm_mul_pd(vd2_tmp, vd2_v1); /* - tmp * v1 */
        vd2_p2 = _mm_mul_pd(vd2_tmp, vd2_v2); /* - tmp * v2 */

        vd2_a0 = _mm_add_pd(vd2_tmp, vd2_a0); /* a0 -= tmp */
        vd2_a1 = _mm_add_pd(vd2_p1, vd2_a1); /* a1 -= tmp * v1 */
        vd2_a2 = _mm_add_pd(vd2_p2, vd2_a2); /* a2 -= tmp * v2 */

        /* store the column elements */
        _mm_storel_pd((double *) &a[ir + 1 * *lda], vd2_a0);
        _mm_storel_pd((double *) &a[ir + 2 * *lda], vd2_a1);
        _mm_storel_pd((double *) &a[ir + 3 * *lda], vd2_a2);
    }

    return 0;
}
#endif
