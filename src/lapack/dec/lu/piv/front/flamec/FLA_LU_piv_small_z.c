/*
    Copyright (c) 2022 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"
#include "FLA_f2c.h"



/*
 * LU with partial pivoting for tiny matrices
 *
 * All the computations are done inline without using
 * corresponding BLAS APIs to reduce function overheads.
 */
integer FLA_LU_piv_small_z_var0( integer *m, integer *n,
                                   doublecomplex *a, integer *lda,
                                   integer *ipiv,
                                   integer *info)
{
    integer mi, ni;
    integer i, j, i_1, i_2, i_3;
    doublereal max_val, t_val, z_val;
    doublecomplex *acur, *apiv, *asrc;
    integer p_idx;
    integer min_m_n = fla_min(*m, *n);
#ifndef _WIN32
    doublecomplex z__1;
    double _Complex pinv;
#else
    doublereal piv_r, piv_i;
    doublereal pinv;
#endif

    for( i = 0; i < min_m_n; i++ )
    {
        mi = *m - i;
        ni = *n - i;

        acur = &a[i + *lda * i];

        /* Find the pivot element */
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

        /* Swap rows and calculate a column of L */
        if( apiv[*lda * i].r != 0. || apiv[*lda * i].i != 0. )
        {
            /* Swap entire rows */
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

            /* Calculate scalefactors (L) & update trailing matrix */

#ifndef _WIN32
            pinv = 1.0 / ((*acur).r + I * (*acur).i);
            z__1.r = creal(pinv);
            z__1.i = cimag(pinv);
#else
            piv_r = (*acur).r;
            piv_i = (*acur).i;
            pinv = piv_r * piv_r + piv_i * piv_i;
#endif

            for( i_1 = 1; i_1 < mi; i_1++ )
            {
                t_val = acur[i_1].r;
#ifndef _WIN32
                acur[i_1].r = (t_val * z__1.r - acur[i_1].i * z__1.i);
                acur[i_1].i = (t_val * z__1.i + acur[i_1].i * z__1.r);
#else
                acur[i_1].r = (acur[i_1].i * piv_i + t_val * piv_r) / pinv;
                acur[i_1].i = (acur[i_1].i * piv_r - t_val * piv_i) / pinv;
#endif
                t_val = acur[i_1].r;
                z_val = acur[i_1].i;

                for( j = 1; j < ni; j++ )
                {
                    i_2 = i_1 + j * *lda;
                    i_3 = j * *lda;

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
