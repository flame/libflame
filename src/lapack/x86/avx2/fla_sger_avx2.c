/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_sger_avx2.c
 *  @brief SGER perfroms the rank 1 operation.
 *  */

#include "FLAME.h"

#ifdef FLA_ENABLE_AMD_OPT
int fla_sger_avx2(integer *m, integer *n, real *alpha, real *x, integer *incx, real *y, integer *incy, real *a, integer *lda)
{
    /* Local Variables */
    integer i__1, i__2, j, jy, i, kx, ix, a_offset;
    real temp, *p;
    i__1 = *n;
    i__2 = *m;
    --x;
    --y;
    a_offset = 1 + (*lda) * 1;
    a -= a_offset;
    __m256 vd4_x[4], vd4_o[4], vd4_a[4], vd4_alpha;
    __m128 vd2_x, vd2_o, vd2_a, vd2_alpha;

    if( *m == 0 || *n == 0 || *alpha == 0.f )
    {
        return 0;
    }

    if (*incy > 0)
    {
        jy = 1;
    }
    else
    {
        jy = 1 - (*n - 1) * *incy;
    }

    if (*incx == 1)
    {
        for (j = 1; j <= i__1; ++j)
        {
            if (y[jy] != 0.f)
            {
                temp = *alpha * y[jy];
                p = &temp;
                /* load scale factor in 256 bit register */
                vd4_alpha = _mm256_broadcast_ss((float *const) p);
                /* load scale factor in 128 bit register */
                vd2_alpha = _mm_broadcast_ss((float *const) p);
                for (i = 1; (i + 63) <= i__2; i += 64)
                {
                    // Load the inputs
                    vd4_x[0] = _mm256_loadu_ps( (float const *) &x[i]);
                    vd4_x[1] = _mm256_loadu_ps( (float const *) &x[i + 8]);
                    vd4_x[2] = _mm256_loadu_ps( (float const *) &x[i + 16]);
                    vd4_x[3] = _mm256_loadu_ps( (float const *) &x[i + 24]);

                    vd4_a[0] = _mm256_loadu_ps( (float const *) &a[(i + 0) + j*(*lda)]);
                    vd4_a[1] = _mm256_loadu_ps( (float const *) &a[(i + 8) + j*(*lda)]);
                    vd4_a[2] = _mm256_loadu_ps( (float const *) &a[(i + 16) + j*(*lda)]);
                    vd4_a[3] = _mm256_loadu_ps( (float const *) &a[(i + 24) + j*(*lda)]);

                    // perform : a := alpha * x + a
                    vd4_o[0] = _mm256_fmadd_ps( vd4_alpha, vd4_x[0], vd4_a[0]);
                    vd4_o[1] = _mm256_fmadd_ps( vd4_alpha, vd4_x[1], vd4_a[1]);
                    vd4_o[2] = _mm256_fmadd_ps( vd4_alpha, vd4_x[2], vd4_a[2]);
                    vd4_o[3] = _mm256_fmadd_ps( vd4_alpha, vd4_x[3], vd4_a[3]);

                    // Store the inputs
                    _mm256_storeu_ps( ((float *) &a[(i + 0) + j*(*lda)]), vd4_o[0]);
                    _mm256_storeu_ps( ((float *) &a[(i + 8) + j*(*lda)]), vd4_o[1]);
                    _mm256_storeu_ps( ((float *) &a[(i + 16) + j*(*lda)]), vd4_o[2]);
                    _mm256_storeu_ps( ((float *) &a[(i + 24) + j*(*lda)]), vd4_o[3]);

                    // Load the inputs
                    vd4_x[0] = _mm256_loadu_ps( (float const *) &x[i + 32]);
                    vd4_x[1] = _mm256_loadu_ps( (float const *) &x[i + 40]);
                    vd4_x[2] = _mm256_loadu_ps( (float const *) &x[i + 48]);
                    vd4_x[3] = _mm256_loadu_ps( (float const *) &x[i + 56]);

                    vd4_a[0] = _mm256_loadu_ps( (float const *) &a[(i + 32) + j*(*lda)]);
                    vd4_a[1] = _mm256_loadu_ps( (float const *) &a[(i + 40) + j*(*lda)]);
                    vd4_a[2] = _mm256_loadu_ps( (float const *) &a[(i + 48) + j*(*lda)]);
                    vd4_a[3] = _mm256_loadu_ps( (float const *) &a[(i + 56) + j*(*lda)]);

                    // perform : a := alpha * x + a
                    vd4_o[0] = _mm256_fmadd_ps( vd4_alpha, vd4_x[0], vd4_a[0]);
                    vd4_o[1] = _mm256_fmadd_ps( vd4_alpha, vd4_x[1], vd4_a[1]);
                    vd4_o[2] = _mm256_fmadd_ps( vd4_alpha, vd4_x[2], vd4_a[2]);
                    vd4_o[3] = _mm256_fmadd_ps( vd4_alpha, vd4_x[3], vd4_a[3]);

                    // Store the inputs
                    _mm256_storeu_ps( ((float *) &a[(i + 32) + j*(*lda)]), vd4_o[0]);
                    _mm256_storeu_ps( ((float *) &a[(i + 40) + j*(*lda)]), vd4_o[1]);
                    _mm256_storeu_ps( ((float *) &a[(i + 48) + j*(*lda)]), vd4_o[2]);
                    _mm256_storeu_ps( ((float *) &a[(i + 56) + j*(*lda)]), vd4_o[3]);
                }

                for ( ; (i + 31) <= i__2; i += 32)
                {
                    // Load the inputs
                    vd4_x[0] = _mm256_loadu_ps( (float const *) &x[i]);
                    vd4_x[1] = _mm256_loadu_ps( (float const *) &x[i + 8]);
                    vd4_x[2] = _mm256_loadu_ps( (float const *) &x[i + 16]);
                    vd4_x[3] = _mm256_loadu_ps( (float const *) &x[i + 24]);

                    vd4_a[0] = _mm256_loadu_ps( (float const *) &a[(i + 0) + j*(*lda)]);
                    vd4_a[1] = _mm256_loadu_ps( (float const *) &a[(i + 8) + j*(*lda)]);
                    vd4_a[2] = _mm256_loadu_ps( (float const *) &a[(i + 16) + j*(*lda)]);
                    vd4_a[3] = _mm256_loadu_ps( (float const *) &a[(i + 24) + j*(*lda)]);

                    // perform : a := alpha * x + a
                    vd4_o[0] = _mm256_fmadd_ps( vd4_alpha, vd4_x[0], vd4_a[0]);
                    vd4_o[1] = _mm256_fmadd_ps( vd4_alpha, vd4_x[1], vd4_a[1]);
                    vd4_o[2] = _mm256_fmadd_ps( vd4_alpha, vd4_x[2], vd4_a[2]);
                    vd4_o[3] = _mm256_fmadd_ps( vd4_alpha, vd4_x[3], vd4_a[3]);

                    // Store the inputs
                    _mm256_storeu_ps( ((float *) &a[(i + 0) + j*(*lda)]), vd4_o[0]);
                    _mm256_storeu_ps( ((float *) &a[(i + 8) + j*(*lda)]), vd4_o[1]);
                    _mm256_storeu_ps( ((float *) &a[(i + 16) + j*(*lda)]), vd4_o[2]);
                    _mm256_storeu_ps( ((float *) &a[(i + 24) + j*(*lda)]), vd4_o[3]);
                }
                for ( ; (i + 15) <= i__2; i += 16)
                {
                    // Load the inputs
                    vd4_x[0] = _mm256_loadu_ps( (float const *) &x[i]);
                    vd4_x[1] = _mm256_loadu_ps( (float const *) &x[i + 8]);

                    vd4_a[0] = _mm256_loadu_ps( (float const *) &a[(i + 0) + j*(*lda)]);
                    vd4_a[1] = _mm256_loadu_ps( (float const *) &a[(i + 8) + j*(*lda)]);

                    // perform : a := alpha * x + a
                    vd4_o[0] = _mm256_fmadd_ps( vd4_alpha, vd4_x[0], vd4_a[0]);
                    vd4_o[1] = _mm256_fmadd_ps( vd4_alpha, vd4_x[1], vd4_a[1]);

                    // Store the inputs
                    _mm256_storeu_ps( ((float *) &a[(i + 0) + j*(*lda)]), vd4_o[0]);
                    _mm256_storeu_ps( ((float *) &a[(i + 8) + j*(*lda)]), vd4_o[1]);
                }
                for ( ; (i + 7) <= i__2; i += 8)
                {
                    // Load the inputs
                    vd4_x[0] = _mm256_loadu_ps( (float const *) &x[i]);

                    // perform : a := alpha * x + a
                    vd4_a[0] = _mm256_loadu_ps( (float const *) &a[(i + 0) + j*(*lda)]);

                    vd4_o[0] = _mm256_fmadd_ps( vd4_alpha, vd4_x[0], vd4_a[0]);

                    // Store the inputs
                    _mm256_storeu_ps( ((float *) &a[(i + 0) + j*(*lda)]), vd4_o[0]);
                }
                for ( ; (i + 3) <= i__2; i += 4)
                {
                    // Load the inputs
                    vd2_x = _mm_loadu_ps( (float const *) &x[i]);

                    vd2_a = _mm_loadu_ps( (float const *) &a[(i + 0) + j*(*lda)]);

                    // perform : a := alpha * x + a
                    vd2_o = _mm_fmadd_ps( vd2_alpha, vd2_x, vd2_a);

                    // Store the inputs
                    _mm_storeu_ps( ((float *) &a[(i + 0) + j*(*lda)]), vd2_o);
                }
                for( ; i <= i__2; ++i)
                {
                    vd2_x = _mm_load_ss( (float const*) &x[i]);

                    vd2_a = _mm_load_ss( (float const*) &a[(i + 0) + j*(*lda)]);

                    vd2_o = _mm_fmadd_ss( vd2_alpha, vd2_x, vd2_a);
                    
                    _mm_store_ss( ((float *) &a[(i + 0) + j*(*lda)]), vd2_o);
                }
            }
            jy += *incy;
        }
    }
    else
    {
        if (*incx > 0)
        {
            kx = 1;
        }
        else
        {
            kx = 1 - (*m - 1) * *incx;
        }
        for (j = 1;
                j <= i__1;
                ++j)
        {
            if (y[jy] != 0.f)
            {
                temp = *alpha * y[jy];
                ix = kx;
                for (i = 1;
                        i <= i__2;
                        ++i)
                {
                    a[i + j * *(lda)] += x[ix] * temp;
                    ix += *incx;
                }
            }
            jy += *incy;
        }
    }
    return 0;
}
#endif