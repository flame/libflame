/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_zscal_ix1_avx2.c
 *  @brief ZSCAL scales a vector by a scalar constant using AVX2 intrinsics.
 *         The vector elements are assumed to be contiguosly stored in memory.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT
int fla_zscal_ix1_avx2(integer *n, doublecomplex *alpha, doublecomplex *x)
{
    /* Local variables */
    integer i__1, i;
    --x;
    i__1 = *n;
    if (i__1 <= 0)
    {
        return 0;
    }

    __m256d srmm, simm, sirmm, srimm;
    __m128d srm, sim, sirm, srim;
    __m256d xmm0, xrmm0, ximm0, ximm1, xmm1, xrmm1;
    __m256d oxm0, oxm1;
    __m128d xmm, oxm, xrmm, ximm;

    /* load scale factor in 256 bit register */
    srmm = _mm256_broadcast_sd((double const *) &alpha->r);
    simm = _mm256_broadcast_sd((double const *) &alpha->i);
    sirmm = _mm256_shuffle_pd(srmm, simm, 0xA);
    srimm = _mm256_shuffle_pd(simm, srmm, 0x5);

    /* load scale factor in 128 bit register */
    srm = _mm_loaddup_pd((double const *) &alpha->r);
    sim = _mm_loaddup_pd((double const *) &alpha->i);
    sirm = _mm_shuffle_pd(srm, sim, 0x2);
    srim = _mm_shuffle_pd(sim, srm, 0x1);

    /* Code for increments equal to 1 only */
    for (i = 1; i <= (i__1 - 3); i += 4)
    {
        /* load complex inputs */
        xmm0   = _mm256_loadu_pd((double const *) &x[i]);
        xmm1   = _mm256_loadu_pd((double const *) &x[i + 2]);

        /* shuffle the loaded inputs */
        xrmm0 = _mm256_movedup_pd(xmm0);
        ximm0 = _mm256_unpackhi_pd(xmm0, xmm0);
        xrmm1 = _mm256_movedup_pd(xmm1);
        ximm1 = _mm256_unpackhi_pd(xmm1, xmm1);

        /* performs the scaling */
        oxm0 = _mm256_mul_pd(srimm, ximm0);
        oxm0 = _mm256_fmaddsub_pd(sirmm, xrmm0, oxm0);
        oxm1 = _mm256_mul_pd(srimm, ximm1);
        oxm1 = _mm256_fmaddsub_pd(sirmm, xrmm1, oxm1);

        /* store the results */
        _mm256_storeu_pd((double *) &x[i], oxm0);
        _mm256_storeu_pd((double *) &x[i + 2], oxm1);
    }

    /* remainder iterations */
    for ( ; i <= i__1; ++i)
    {
        /* load inputs */
        xmm  = _mm_loadu_pd((double const *) &x[i]);

        /* shuffle inputs */
        xrmm = _mm_movedup_pd(xmm);
        ximm = _mm_unpackhi_pd(xmm, xmm);

        /* performs scaling */
        oxm = _mm_mul_pd(srim, ximm);
        oxm = _mm_fmaddsub_pd(sirm, xrmm, oxm);

        /* store result */
        _mm_storeu_pd((double *) &x[i], oxm);
    }
    return 0;
}
#endif
