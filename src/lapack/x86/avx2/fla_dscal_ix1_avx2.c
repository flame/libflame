/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_dscal_ix1_avx2.c
 *  @brief scales a vector by a constant
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#ifdef FLA_ENABLE_AMD_OPT

int fla_dscal_ix1_avx2(integer *n, doublereal *da, doublereal *dx, integer *incx)
{
    /* Parameter adjustments */
    --dx;
    /* Function Body */
    if (*n <= 0)
    {
        return 0;
    }

    integer i, i__1;
    doublereal d__1;
    i__1 = *n;
    d__1 = *da;

    /* Load scaling factor alpha*/
    __m256d alphav = _mm256_set1_pd(d__1);

    for (i = 1; i <= (i__1 - 3); i += 4)
    {
        /* Load the input values */
        __m256d x0v = _mm256_loadu_pd((double const *) &dx[i]);
        /* perform alpha * x  */
        x0v = _mm256_mul_pd( alphav, x0v );
        /* Store the output */
        _mm256_storeu_pd((double *) &dx[i], x0v);
    }
    /* Remainder iterations */
    if((i__1-i) >= 2)
    {
        for ( ; i <= (i__1-1); i += 2 )
        {
            dx[i] *= d__1;
            dx[i+1] *= d__1;
        }
        for ( ; i <= i__1; ++i )
        {
            dx[i] *= d__1;
        }
    }
    else
    {
        for ( ; i <= i__1; ++i )
        {
            dx[i] *= d__1;
        }
    }
    return 0;
}
#endif
