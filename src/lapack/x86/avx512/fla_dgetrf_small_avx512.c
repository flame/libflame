/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#ifdef FLA_ENABLE_AMD_OPT

/*
 * LU with partial pivoting for tiny matrices
 *
 * All the computations are done inline without using
 * corresponding BLAS APIs to reduce function overheads.
 */
integer fla_dgetrf_small_avx512( integer *m, integer *n,
                                   doublereal *a, integer *lda,
                                   integer *ipiv,
                                   integer *info)
{
    integer mi, ni;
    integer i, j, i_1, lda_t, b_off, y_off;

    doublereal p_val, max_val, t_val;
    doublereal *acur, *apiv, *asrc;
    integer p_idx;
    __m512d result[8], tempY[8], tempb[8], tempx, p_val4;
    integer min_m_n = fla_min(*m, *n);
    lda_t = *lda;

    for( i = 0; i < min_m_n; i++ )
    {
        mi = *m - i;
        ni = *n - i;

        acur = &a[i + lda_t * i];

        /* Find the pivot element */
        max_val = 0;
        p_idx = i;
        for( i_1 = 0; i_1 < mi; i_1++ )
        {
            t_val = acur[i_1];
            t_val = ( t_val < 0.0 ) ? -t_val : t_val;
            if( t_val > max_val )
            {
                max_val = t_val;
                p_idx = i + i_1;
            }
        }

        apiv = a + p_idx;
        asrc = a + i;
        ipiv[i] = p_idx + 1;

        /* Swap rows and calculate a column of L */
        if( max_val != 0.0 )
        {
            /* Swap entire rows */
            if( p_idx != i)
            {
                for( i_1 = 0; i_1 < *n; i_1++ )
                {
                    t_val = apiv[i_1 * lda_t];
                    apiv[i_1 * *lda] = asrc[i_1 * lda_t];
                    asrc[i_1 * *lda] = t_val;
                }
            }

            /* Calculate scalefactors (L)  & update trailing matrix */
            p_val = *acur;
            p_val = 1 / p_val;
            p_val4 = _mm512_set1_pd(p_val);
            for( i_1 = 1; i_1 < mi-7; i_1+=8 )
            {
                tempx = _mm512_mul_pd(_mm512_loadu_pd(&acur[i_1]), p_val4);
                _mm512_storeu_pd(&acur[i_1], tempx);

                for( j = 1; j < ni-7; j+=8 )
                {
                    b_off = j * lda_t;
                    y_off = i_1 + j * lda_t;
                    tempb[0] = _mm512_set1_pd(acur[b_off]);
                    tempb[1] = _mm512_set1_pd(acur[b_off + 1*lda_t]);
                    tempb[2] = _mm512_set1_pd(acur[b_off + 2*lda_t]);
                    tempb[3] = _mm512_set1_pd(acur[b_off + 3*lda_t]);
                    tempb[4] = _mm512_set1_pd(acur[b_off + 4*lda_t]);
                    tempb[5] = _mm512_set1_pd(acur[b_off + 5*lda_t]);
                    tempb[6] = _mm512_set1_pd(acur[b_off + 6*lda_t]);
                    tempb[7] = _mm512_set1_pd(acur[b_off + 7*lda_t]);
                    tempY[0] = _mm512_loadu_pd(&acur[y_off]);
                    tempY[1] = _mm512_loadu_pd(&acur[y_off + 1*lda_t]);
                    tempY[2] = _mm512_loadu_pd(&acur[y_off + 2*lda_t]);
                    tempY[3] = _mm512_loadu_pd(&acur[y_off + 3*lda_t]);
                    tempY[4] = _mm512_loadu_pd(&acur[y_off + 4*lda_t]);
                    tempY[5] = _mm512_loadu_pd(&acur[y_off + 5*lda_t]);
                    tempY[6] = _mm512_loadu_pd(&acur[y_off + 6*lda_t]);
                    tempY[7] = _mm512_loadu_pd(&acur[y_off + 7*lda_t]);
                    /* Y := Y - b * x */
                    result[0] = _mm512_fnmadd_pd(tempb[0], tempx, tempY[0]);
                    result[1] = _mm512_fnmadd_pd(tempb[1], tempx, tempY[1]);
                    result[2] = _mm512_fnmadd_pd(tempb[2], tempx, tempY[2]);
                    result[3] = _mm512_fnmadd_pd(tempb[3], tempx, tempY[3]);
                    result[4] = _mm512_fnmadd_pd(tempb[4], tempx, tempY[4]);
                    result[5] = _mm512_fnmadd_pd(tempb[5], tempx, tempY[5]);
                    result[6] = _mm512_fnmadd_pd(tempb[6], tempx, tempY[6]);
                    result[7] = _mm512_fnmadd_pd(tempb[7], tempx, tempY[7]);
                    _mm512_storeu_pd(&acur[y_off], result[0]);
                    _mm512_storeu_pd(&acur[y_off + 1 * lda_t], result[1]);
                    _mm512_storeu_pd(&acur[y_off + 2 * lda_t], result[2]);
                    _mm512_storeu_pd(&acur[y_off + 3 * lda_t], result[3]);
                    _mm512_storeu_pd(&acur[y_off + 4 * lda_t], result[4]);
                    _mm512_storeu_pd(&acur[y_off + 5 * lda_t], result[5]);
                    _mm512_storeu_pd(&acur[y_off + 6 * lda_t], result[6]);
                    _mm512_storeu_pd(&acur[y_off + 7 * lda_t], result[7]);
                }
                /* remining inner loop updation*/
                for (; j < ni; j++)
                {
                    b_off = j * lda_t;
                    y_off = i_1 + j * lda_t;
                    tempb[0] = _mm512_set1_pd(acur[b_off]);
                    tempY[0] = _mm512_loadu_pd(&acur[y_off]);
                    result[0] = _mm512_fnmadd_pd(tempb[0], tempx, tempY[0]);
                    _mm512_storeu_pd(&acur[y_off], result[0]);
                }
            }
            /* remining outer loop last 7 iteration updat with avx2 */
            fla_lu_piv_small_d_update_tr_matrix_avx2(i_1, mi, ni, acur, *lda);
        }
        else
        {
            *info = ( *info == 0 ) ? p_idx + 1 : *info;
        }
    }
    return *info;
}
#endif
