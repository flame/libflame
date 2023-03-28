/*
    Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"
#include "FLA_f2c.h"

#define FLA_LU_SMALL_BLOCK_SIZE 16

static doublecomplex z__1 = { -1, 0}, c_b1 = {1.,0.};
static integer c__1 = 1;

/*
 * LU with partial pivoting for tiny matrices
 *
 * All the computations are done inline without using
 * corresponding BLAS APIs to reduce function overheads.
 */
integer FLA_LU_piv_small_z_var0( integer *m, integer *n,
                                   dcomplex *a, integer *lda,
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

    *info = 0;

    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -4;
    }

    if (*info != 0)
    {
        return 0;
    }

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
            pinv = 1.0 / ((*acur).r + (I * (*acur).i));
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


/* LU factorization recursive variant*/
integer FLA_LU_piv_z_var0(integer *m, integer *n, doublecomplex *a, integer *lda, integer *ipiv, integer *info)
{
    integer a_dim1, i__1, i__2, i__, n1, n2;
    integer iinfo;

    /* Adjust dimension of the matrix */
    a_dim1 = *lda;

    /* Function Body */
    *info = 0;
    if (*m < 0)
    {
        *info = -1;
    }
    else if (*n < 0)
    {
        *info = -2;
    }
    else if (*lda < fla_max(1,*m))
    {
        *info = -4;
    }

    if (*info != 0)
    {
        return 0;
    }

    /* Quick return if possible */
    if (*m == 0 || *n == 0)
    {
        return 0;
    }

    if (*m <= FLA_LU_SMALL_BLOCK_SIZE && *n <= FLA_LU_SMALL_BLOCK_SIZE)
    {
       fla_zgetrf_small_avx2(m, n, a, lda, ipiv, &iinfo);

       if (*info == 0 && iinfo > 0)
       {
           *info = iinfo;
       }
    }
    else if (*m <= FLA_LU_SMALL_BLOCK_SIZE || *n <= FLA_LU_SMALL_BLOCK_SIZE)
    {
        lapack_zgetf2(m, n, a, lda, ipiv, &iinfo);

        if (*info == 0 && iinfo > 0)
        {
            *info = iinfo;
        }
    }
    else
    {
        /* Use recursive code */

        /* calculate n1 and n2 for recursive call*/
        n1 = fla_min(*m,*n) / 2;
        n2 = *n - n1;

        /*        [ A11 ] */
        /* Factor [ --- ] */
        /*        [ A21 ] */
        FLA_LU_piv_z_var0(m, &n1, a, lda, ipiv, &iinfo);

        /* Update info if necessary*/
        if (*info == 0 && iinfo > 0)
        {
            *info = iinfo;
        }

        /*                       [ A12 ] */
        /* Apply interchanges to [ --- ] */
        /*                       [ A22 ] */
        zlaswp_(&n2, &a[(n1 ) * a_dim1 ], lda, &c__1, &n1, ipiv, & c__1);

        /* Solve A12 */
        ztrsm_("L", "L", "N", "U", &n1, &n2, &c_b1, a, lda, &a[(n1 ) * a_dim1], lda);

        /* Update A22 */
        i__1 = *m - n1;
        zgemm_("N", "N", &i__1, &n2, &n1, &z__1, &a[n1], lda, &a[ (n1) * a_dim1 ], lda, &c_b1, &a[n1 +  (n1) * a_dim1], lda);

        /* Factor A22 */
        i__1 = *m - n1;
        FLA_LU_piv_z_var0(&i__1, &n2, &a[n1 + (n1) * a_dim1], lda, &ipiv[n1], &iinfo);

        /* Adjust INFO and the pivot indices */
        if (*info == 0 && iinfo > 0)
        {
            *info = iinfo + n1;
        }

        /* Update pivot indices*/
        i__1 = fla_min(*m,*n);
        for (i__ = n1; i__ < i__1; ++i__)
        {
            ipiv[i__] += n1;
        }

        /* Apply interchanges to A21 */
        i__1 = n1 + 1;
        i__2 = fla_min(*m,*n);
        zlaswp_(&n1, a, lda, &i__1, &i__2, ipiv, &c__1);
    }

    return 0;
}