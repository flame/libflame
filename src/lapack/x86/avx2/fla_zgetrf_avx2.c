/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

#include "FLAME.h"

#ifdef FLA_ENABLE_AMD_OPT

/*
 * LU with partial pivoting for tiny matrices
 *
 * All the computations are done inline without using
 * corresponding BLAS APIs to reduce function overheads.
 */
int FLA_LU_piv_small_z_avx2( integer *m, integer *n, doublecomplex *a, integer *lda, integer *ipiv, integer *info)
{
    integer mi, ni;
    integer i, j, i_1, i_2, i_3, i_4, i_5, i_6, i_7;
    doublereal max_val, t_val, z_val, x_val, y_val;
    doublecomplex *acur, *apiv, *asrc;
    doublecomplex z__1, y__1 = {1, 0};
    integer p_idx;
    integer min_m_n = fla_min(*m, *n);
    __m256d alpha;
    __m256d bv[2], bv_p[2], temp[8], xv0[2], xv1[2], yv0[2], yv1[2];
    __m256d neg = _mm256_setr_pd(1.0, -1.0, 1.0, -1.0);
#ifndef _WIN32
    double _Complex pinv;
#endif

    for( i = 0; i < min_m_n; i++ )
    {
        mi = *m - i;
        ni = *n - i;

        acur = &a[i + *lda * i];

        // Find the pivot element
        max_val = 0;
        p_idx = i;
        for( i_1 = 0; i_1 < mi; i_1++ )
        {
            t_val = f2c_abs(acur[i_1].r) + f2c_abs(acur[i_1].i);
            if( t_val > max_val )
            {
                max_val = t_val;
                p_idx = i + i_1;
            }
        }

        apiv = a + p_idx;
        asrc = a + i;
        ipiv[i] = p_idx + 1;

        // Swap rows
        if( apiv[*lda * i].r != 0. || apiv[*lda * i].i != 0. )
        {
            if( p_idx != i )
            {
                for( i_1 = 0; i_1 < *n ; i_1++ )
                {
                    i_2 = i_1 * *lda;
                    t_val = apiv[i_2].r;
                    z_val = apiv[i_2].i;
                    apiv[i_2].r = asrc[i_2].r;
                    apiv[i_2].i = asrc[i_2].i;
                    asrc[i_2].r = t_val;
                    asrc[i_2].i = z_val;
                 }
            }

            /*----------------unblocked LU algorithm-------------------------

            A00 | a01    A02
            ----|-----------
            a10 | a11    a12
            A20 | a21    A22

            alpha = 1 / a11
            a21 := a21 * alpha
            A22 := A22 - a21 * a12
            ------------------------------------------------------------------*/

            // Calculate scalefactors (a21) & update trailing matrix
#ifndef _WIN32
            pinv = 1.0 / ((*acur).r + I * (*acur).i);
            z__1.r = creal(pinv);
            z__1.i = cimag(pinv);
#else
            dladiv_(&y__1.r, &y__1.i, &acur->r, &acur->i, &z__1.r, &z__1.i);
#endif

            // Load alpha from memory
            alpha = _mm256_set_pd(z__1.i, z__1.r, z__1.i, z__1.r);

            // Updates 4 rows of trailing matrix per iteration
            for(i_1 = 1; i_1 < mi - 3; i_1+=4 )
            {
                /*-----------Trailing matrix update for LU factorisation-----------

                A00 | a01    A02
                ----|-----------
                a10 | alpha  x
                A20 | b      Y

                b := alpha * b
                Y := Y - b * x


                SIMD algorithm:

                alpha = aR1  aI1  aR1  aI1
                bv    = bR1  bI1  bR2  bI2
                bv_p  = bI1  bR1  bI2  bR2
                xv    = xR1  xI1  xR1  xI1
                yv    = yR1  yI1  yR2  yI2

                step 1 => b := alpha * b
                    bv = alpha * bv
                    bv_p = alpha * (-bv_p)
                    bv = bv - bv_p
                    bv_p = shuffle(bv)

                step 2 => Y := Y - b * x
                    temp = xv * bv
                    temp1 = xv * (-bv_p)
                    xv = temp0 - temp1
                    yv = yv - xv
                ------------------------------------------------------------------*/

                i_2 = i_1 + 2;

                // Load alpha from memory
                bv[0] = _mm256_loadu_pd((double const *) &acur[i_1].r);
                bv[1] = _mm256_loadu_pd((double const *) &acur[i_2].r);
                bv_p[0] = _mm256_permute_pd(bv[0], 0x5);
                bv_p[1] = _mm256_permute_pd(bv[1], 0x5);

                // b := alpha * b
                bv[0] = _mm256_mul_pd(alpha, bv[0]);
                bv_p[0] = _mm256_mul_pd(bv_p[0], neg); 
                bv_p[0] = _mm256_mul_pd(alpha, bv_p[0]);
                bv[1] = _mm256_mul_pd(alpha, bv[1]);
                bv_p[1] = _mm256_mul_pd(bv_p[1], neg); 
                bv_p[1] = _mm256_mul_pd(alpha, bv_p[1]);

                bv[0] = _mm256_hsub_pd(bv[0], bv_p[0]);
                bv[1] = _mm256_hsub_pd(bv[1], bv_p[1]);

                _mm256_storeu_pd ((double *) &acur[i_1].r, bv[0]);
                _mm256_storeu_pd ((double *) &acur[i_2].r, bv[1]);

                bv_p[0] = _mm256_permute_pd(bv[0], 0x5);
                bv_p[1] = _mm256_permute_pd(bv[1], 0x5);
                bv_p[0] = _mm256_mul_pd(bv_p[0], neg);
                bv_p[1] = _mm256_mul_pd(bv_p[1], neg);

                for( j = 1; j < ni - 1; j = j + 2 )
                {
                    i_3 = j * *lda;
                    i_2 = i_1 + i_3;

                    i_4 = (j + 1) * *lda;
                    i_5 = i_1 + i_4;

                    i_6 = i_2 + 2;
                    i_7 = i_5 + 2;

                    // Load x from memory
                    xv0[0] = _mm256_set_pd(acur[i_3].i, acur[i_3].r, acur[i_3].i, acur[i_3].r);
                    xv1[0] = _mm256_set_pd(acur[i_4].i, acur[i_4].r, acur[i_4].i, acur[i_4].r);

                    // Y := Y - b * x
                    temp[0] = _mm256_mul_pd(xv0[0], bv[0]);
                    temp[1] = _mm256_mul_pd(xv0[0], bv_p[0]);
                    temp[2] = _mm256_mul_pd(xv0[0], bv[1]);
                    temp[3] = _mm256_mul_pd(xv0[0], bv_p[1]);
                    temp[4] = _mm256_mul_pd(xv1[0], bv[0]);
                    temp[5] = _mm256_mul_pd(xv1[0], bv_p[0]);
                    temp[6] = _mm256_mul_pd(xv1[0], bv[1]);
                    temp[7] = _mm256_mul_pd(xv1[0], bv_p[1]);

                    xv0[0] = _mm256_hsub_pd(temp[0], temp[1]);
                    xv0[1] = _mm256_hsub_pd(temp[2], temp[3]);
                    xv1[0] = _mm256_hsub_pd(temp[4], temp[5]);
                    xv1[1] = _mm256_hsub_pd(temp[6], temp[7]);

                    yv0[0] = _mm256_loadu_pd((double const *) &acur[i_2].r);
                    yv0[1] = _mm256_loadu_pd((double const *) &acur[i_6].r);
                    yv1[0] = _mm256_loadu_pd((double const *) &acur[i_5].r);
                    yv1[1] = _mm256_loadu_pd((double const *) &acur[i_7].r);

                    yv0[0] = _mm256_sub_pd(yv0[0], xv0[0]);
                    yv0[1] = _mm256_sub_pd(yv0[1], xv0[1]);
                    yv1[0] = _mm256_sub_pd(yv1[0], xv1[0]);
                    yv1[1] = _mm256_sub_pd(yv1[1], xv1[1]);

                    _mm256_storeu_pd ((double *) &acur[i_2].r, yv0[0]);
                    _mm256_storeu_pd ((double *) &acur[i_6].r, yv0[1]);
                    _mm256_storeu_pd ((double *) &acur[i_5].r, yv1[0]);
                    _mm256_storeu_pd ((double *) &acur[i_7].r, yv1[1]);
                }
                if(ni - j > 0)
                {
                    i_3 = j * *lda;
                    i_2 = i_1 + i_3;

                    i_4 = i_2 + 2;

                    // Load x from memory
                    xv0[0] = _mm256_set_pd(acur[i_3].i, acur[i_3].r, acur[i_3].i, acur[i_3].r);

                    // Y := Y - b * x
                    temp[0] = _mm256_mul_pd(xv0[0], bv[0]);
                    temp[1] = _mm256_mul_pd(xv0[0], bv_p[0]);
                    temp[2] = _mm256_mul_pd(xv0[0], bv[1]);
                    temp[3] = _mm256_mul_pd(xv0[0], bv_p[1]);

                    xv0[0] = _mm256_hsub_pd(temp[0], temp[1]);
                    xv0[1] = _mm256_hsub_pd(temp[2], temp[3]);

                    yv0[0] = _mm256_loadu_pd((double const *) &acur[i_2].r);
                    yv0[1] = _mm256_loadu_pd((double const *) &acur[i_4].r);

                    yv0[0] = _mm256_sub_pd(yv0[0], xv0[0]);
                    yv0[1] = _mm256_sub_pd(yv0[1], xv0[1]);

                    _mm256_storeu_pd ((double *) &acur[i_2].r, yv0[0]);
                    _mm256_storeu_pd ((double *) &acur[i_4].r, yv0[1]);
                }
            }

            // Updates 1 row of trailing matrix per iteration
            for( ; i_1 < mi; i_1++ )
            {
                t_val = acur[i_1].r;
                acur[i_1].r = (t_val * z__1.r - acur[i_1].i * z__1.i);
                acur[i_1].i = (t_val * z__1.i + acur[i_1].i * z__1.r);

                t_val = acur[i_1].r;
                z_val = acur[i_1].i;

                for( j = 1; j < ni; j++ )
                {
                    i_3 = j * *lda;
                    i_2 = i_1 + i_3;

                    acur[i_2].r = acur[i_2].r - t_val * acur[i_3].r + z_val * acur[i_3].i;
                    acur[i_2].i = acur[i_2].i - t_val * acur[i_3].i - z_val * acur[i_3].r;
                }
            }
        }
        else
        {
            *info = ( *info == 0 ) ? p_idx + 1 : *info;
        }
    }
    return *info;
}

#endif