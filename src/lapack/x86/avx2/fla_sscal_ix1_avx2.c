/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_sscal_ix1_avx2.c
 *  @brief SSCAL scales a vector by a constant in AVX2.
 *  */

#include "FLAME.h"

#ifdef FLA_ENABLE_AMD_OPT
int fla_sscal_ix1_avx2(integer *n, real *alpha, real *x)
{
    /* Local variables */
    integer i__1, i;
    i__1 = *n;
    if (i__1 <= 0)
    {
        return 0;
    }

    __m256 vd4_x[4], vd4_o[4], vd4_alpha;
    __m128 vd2_x, vd2_o, vd2_alpha;

    /* load scale factor in 256 bit register */
    vd4_alpha = _mm256_broadcast_ss((float const *) alpha);

    /* load scale factor in 128 bit register */
    vd2_alpha = _mm_broadcast_ss((float const *) alpha);

    for ( i = 0; (i + 63) < i__1; i += 64)
    {
        // Load the input values.
        vd4_x[0] = _mm256_loadu_ps( (float const *) &x[i]);
        vd4_x[1] = _mm256_loadu_ps( (float const *) &x[i + 8]);
        vd4_x[2] = _mm256_loadu_ps( (float const *) &x[i + 16]);
        vd4_x[3] = _mm256_loadu_ps( (float const *) &x[i + 24]);

        // perform : x := alpha * x;
        vd4_o[0] = _mm256_mul_ps( vd4_alpha, vd4_x[0] );
        vd4_o[1] = _mm256_mul_ps( vd4_alpha, vd4_x[1] );
        vd4_o[2] = _mm256_mul_ps( vd4_alpha, vd4_x[2] );
        vd4_o[3] = _mm256_mul_ps( vd4_alpha, vd4_x[3] );
        
        // Store the output.
        _mm256_storeu_ps( ((float *) &x[i]), vd4_o[0]);
        _mm256_storeu_ps( ((float *) &x[i + 8]), vd4_o[1]);
        _mm256_storeu_ps( ((float *) &x[i + 16]), vd4_o[2]);
        _mm256_storeu_ps( ((float *) &x[i + 24]), vd4_o[3]);

        // Load the input values.
        vd4_x[0] = _mm256_loadu_ps( (float const *) &x[i + 32]);
        vd4_x[1] = _mm256_loadu_ps( (float const *) &x[i + 40]);
        vd4_x[2] = _mm256_loadu_ps( (float const *) &x[i + 48]);
        vd4_x[3] = _mm256_loadu_ps( (float const *) &x[i + 56]);

        // perform : x := alpha * x;
        vd4_o[0] = _mm256_mul_ps( vd4_alpha, vd4_x[0] );
        vd4_o[1] = _mm256_mul_ps( vd4_alpha, vd4_x[1] );
        vd4_o[2] = _mm256_mul_ps( vd4_alpha, vd4_x[2] );
        vd4_o[3] = _mm256_mul_ps( vd4_alpha, vd4_x[3] );

        // Store the output.
        _mm256_storeu_ps( ((float *) &x[i + 32]), vd4_o[0]);
        _mm256_storeu_ps( ((float *) &x[i + 40]), vd4_o[1]);
        _mm256_storeu_ps( ((float *) &x[i + 48]), vd4_o[2]);
        _mm256_storeu_ps( ((float *) &x[i + 56]), vd4_o[3]);

    }

    for ( ; (i + 31) < i__1; i += 32)
    {
        // Load the input values.
        vd4_x[0] = _mm256_loadu_ps( (float const *) &x[i]);
        vd4_x[1] = _mm256_loadu_ps( (float const *) &x[i + 8]);
        vd4_x[2] = _mm256_loadu_ps( (float const *) &x[i + 16]);
        vd4_x[3] = _mm256_loadu_ps( (float const *) &x[i + 24]);

        // perform : x := alpha * x;
        vd4_o[0] = _mm256_mul_ps( vd4_alpha, vd4_x[0] );
        vd4_o[1] = _mm256_mul_ps( vd4_alpha, vd4_x[1] );
        vd4_o[2] = _mm256_mul_ps( vd4_alpha, vd4_x[2] );
        vd4_o[3] = _mm256_mul_ps( vd4_alpha, vd4_x[3] );

        // Store the output.
        _mm256_storeu_ps( ((float *) &x[i]), vd4_o[0]);
        _mm256_storeu_ps( ((float *) &x[i + 8]), vd4_o[1]);
        _mm256_storeu_ps( ((float *) &x[i + 16]), vd4_o[2]);
        _mm256_storeu_ps( ((float *) &x[i + 24]), vd4_o[3]);

    }

    for ( ; (i + 15) < i__1; i += 16)
    {
        /* load complex inputs */
        vd4_x[0] = _mm256_loadu_ps((float const *) &x[i]);
        vd4_x[1] = _mm256_loadu_ps((float const *) &x[i + 8]);
        
        /* performs the scaling */
        vd4_o[0] = _mm256_mul_ps(vd4_x[0], vd4_alpha);
        vd4_o[1] = _mm256_mul_ps(vd4_x[1], vd4_alpha);
        
        /* store the results */
        _mm256_storeu_ps((float *) &x[i], vd4_o[0]);
        _mm256_storeu_ps((float *) &x[i + 8], vd4_o[1]);
    }

    /* Code for increments equal to 1 only */
    for ( ; (i + 7) < i__1; i += 8)
    {
        /* load complex inputs */
        vd4_x[0] = _mm256_loadu_ps((float const *) &x[i]);
        
        /* performs the scaling */
        vd4_o[0] = _mm256_mul_ps(vd4_x[0], vd4_alpha);
        
        /* store the results */
        _mm256_storeu_ps((float *) &x[i], vd4_o[0]);
    }

    /* remainder iterations */
    for ( ; (i + 0) < i__1; i+= 1)
    {
        /* load complex inputs */
        vd2_x = _mm_load_ss((float const *) &x[i]);
        
        /* performs the scaling */
        vd2_o = _mm_mul_ss(vd2_x, vd2_alpha);
        
        /* store the results */
        _mm_store_ss((float *) &x[i], vd2_o);
    }
    return 0;
}
#endif