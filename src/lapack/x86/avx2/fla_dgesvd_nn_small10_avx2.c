/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file fla_dgesvd_nn_small10_avx2.c
 *  @brief DGESVD Small path (path 10)
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

double d_sign(doublereal *, doublereal *);

void fla_dgesvd_nn_small10_avx2(integer *m, integer *n,
                                doublereal *a, integer *lda,
                                doublereal *s,
                                doublereal *work,
                                integer *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    doublereal d__1;
    doublereal *tau, *tauq, *taup;
    doublereal *e;
    doublereal dum[1];

    integer c__0 = 0;
    integer c__1 = 1;

    integer iu, ie, iwork;
    integer itau, itauq, itaup;
    integer i__1, rlen, knt;

    integer ldu_val = 0;
    integer *ldu = &ldu_val;

    /* indices for partitioning work buffer */
    iu = 1;
    itau = iu + *m * *ldu;
    ie = itau + *m;
    itauq = ie + *m;
    itaup = itauq + *m;
    iwork = itaup + *m;

    /* parameter adjustments */
    a -= (1 + *lda);
    --s;
    --work;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* Upper Bidiagonalization */
    FLA_BIDIAGONALIZE_SMALL(*m, *n);

    /* Compute Singular Values */
    lapack_dbdsqr("U", n, &c__0, &c__0, &c__0, &s[1], &e[1],
                  NULL, &c__1,
                  NULL, &c__1,
                  dum, &c__1,
                  &work[iwork], info);

    return;
}
#endif
