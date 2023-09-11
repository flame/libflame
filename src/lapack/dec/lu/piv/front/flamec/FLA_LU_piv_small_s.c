/*
    Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"
#include "fla_lapack_x86_common.h"

/*
 * LU with partial pivoting for tiny matrices
 *
 * All the computations are done inline without using
 * corresponding BLAS APIs to reduce function overheads.
 */
integer FLA_LU_piv_small_s_var0( integer *m, integer *n,
                                   real *a, integer *lda,
                                   integer *ipiv,
                                   integer *info)
{
    integer mi, ni;
    integer i, j, i_1;

    real p_val, max_val, t_val;
    real *acur, *apiv, *asrc;
    integer p_idx;
    integer min_m_n = fla_min(*m, *n);

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
            t_val = acur[i_1];
            t_val = ( t_val < 0 ) ? -t_val : t_val;
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
        if( max_val != 0 )
        {
            /* Swap entire rows */
            if( p_idx != i)
            {
                for( i_1 = 0; i_1 < *n; i_1++ )
                {
                    t_val = apiv[i_1 * *lda];
                    apiv[i_1 * *lda] = asrc[i_1 * *lda];
                    asrc[i_1 * *lda] = t_val;
                }
            }

            /* Calculate scalefactors (L)  & update trailing matrix */
            p_val = *acur;
            for( i_1 = 1; i_1 < mi; i_1++ )
            {
                acur[i_1] = t_val = acur[i_1] / p_val;
                for( j = 1; j < ni; j++ )
                {
                    acur[i_1 + j * *lda] = acur[i_1 + j * *lda] - acur[j * *lda] * t_val;
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

/*
 * LU with partial pivoting for small matrices
 *
 * This is an unblocked variant making use of BLAS APIs
 */
/* Further Optimizations possible:
 *
 *  TODO: AVX optimizations to be done 
 */
integer FLA_LU_piv_small_s_var1( integer *m, integer *n, 
                                  real *a, integer *lda,
                                  integer *ipiv,
                                  integer *info)
{
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real d__1;
    integer c__1 = 1;
    real c_n1 = -1.;

    /* Local variables */
    integer i__, j, jp;
    extern real slamch_(char *);
    extern integer isamax_(integer *, real *, integer *);
    extern /* Subroutine */ int xerbla_(const char *srname, const integer *info, ftnlen srname_len);
    real sfmin;
    
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipiv;

    /* Compute machine safe minimum */
    sfmin = slamch_("S");

    i__1 = fla_min(*m,*n);
    for ( j = 1; j <= i__1; ++j )
    {
        /* Find pivot and test for singularity. */
        i__2 = *m - j + 1;
        jp = j - 1 + isamax_(&i__2, &a[j + j * a_dim1], &c__1);
        ipiv[j] = jp;
        if (a[jp + j * a_dim1] != 0.)
        {
            /*Apply the interchange to columns 1:N. */
            if (jp != j)
            {
                sswap_(n, &a[j + a_dim1], lda, &a[jp + a_dim1], lda);
            }
            /*Compute elements J+1:M of J-th column. */
            if (j < *m)
            {
                d__1 = a[j + j * a_dim1];
                d__1 = (d__1 < 0) ? -d__1 : d__1;
                if (d__1 >= sfmin)
                {
                    i__2 = *m - j;
                    d__1 = 1. / a[j + j * a_dim1];
                    fla_sscal(&i__2, &d__1, &a[j + 1 + j * a_dim1], &c__1);
                }
                else
                {
                    i__2 = *m - j;
                    for (i__ = 1; i__ <= i__2; ++i__)
                    {
                        a[j + i__ + j * a_dim1] /= a[j + j * a_dim1];
                    }
                }
            }
        }
        else if ( *info == 0 )
        {
            *info = j;
        }
        if ( j < fla_min( *m, *n ) )
        {
            /* Update trailing submatrix. */
            i__2 = *m - j;
            i__3 = *n - j;
            fla_sger(&i__2, &i__3, &c_n1, &a[j + 1 + j * a_dim1], &c__1, &a[j + (
                   j + 1) * a_dim1], lda, &a[j + 1 + (j + 1) * a_dim1], lda);
        }
    }
    return *info;
}