/******************************************************************************
 * Copyright (C) 2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_sgesvd_nn_small1T_avx2.c
 *  @brief SGESVD Small path (path 1T)
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#ifdef FLA_ENABLE_AMD_OPT

void fla_sgesvd_nn_small1T_avx2(aocl_int64_t *m, aocl_int64_t *n, real *a, aocl_int64_t *lda,
                                real *s, real *work, aocl_int64_t *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_SSMALL();

    real d__1;
    real *tau, *tauq, *taup;
    real *e;

    aocl_int64_t c__0 = 0;
    aocl_int64_t c__1 = 1;

    aocl_int64_t ie;
    aocl_int64_t itauq, itaup;
    aocl_int64_t rlen, knt;

    /* indices for partitioning work buffer */
    ie = 1;
    itauq = ie + *m;
    itaup = itauq + *m;

    /* parameter adjustments */
    a -= (1 + *lda);
    --s;
    --work;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* Upper Bidiagonalization */
    FLA_BIDIAGONALIZE_SSMALL(*m, *m, a, lda, tauq, taup, s, e);

    /* Compute Singular Values and Vectors */
    lapack_sbdsqr_small("U", m, &c__0, &c__0, &s[1], &e[1], NULL, &c__1, NULL, &c__1, info);

    return;
}
#endif
