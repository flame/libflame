/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_dgeqrf_small_avx2.c
 *  @brief QR for small inputs in AVX2.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

static integer c__1 = 1;
/* QR for small sizes */
int fla_dgeqrf_small_avx2(integer *m, integer *n,
                          doublereal *a, integer *lda,
                          doublereal *tau, doublereal *work)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    integer min_m_n;

    /* Adjust pointers */
    a -= (1 + *lda * 1);
    tau--;
    work--;

    min_m_n = fla_min(*m, *n);
    for (i = 1; i <= min_m_n; i++)
    {
        slen = *m - i;
        /* input address */
        const doublereal *iptr = (const doublereal *) &a[i + 1 + i * *lda - 1];
        integer has_outliers = 0;

        if (slen <= 0)
        {
            tau[i] = 0;
        }
        else if (slen < 4)
        {
            FLA_ELEM_REFLECTOR_GENERATE_DSMALL(i, m, n, tau);
            FLA_ELEM_REFLECTOR_APPLY_DSMALL(i, m, n, tau);
        }
        else
        {
            FLA_ELEM_REFLECTOR_GENERATE_DLARGE(i, m, n, tau);
            FLA_ELEM_REFLECTOR_APPLY_DLARGE(i, m, n, tau);
        }
    }
    return 0;
}
#endif
