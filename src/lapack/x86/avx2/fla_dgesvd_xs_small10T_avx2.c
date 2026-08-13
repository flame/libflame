/******************************************************************************
 * Copyright (C) 2023-2026, Advanced Micro Devices, Inc. All rights reserved.
 *******************************************************************************/

/*! @file fla_dgesvd_small6_avx2_.c
 *  @brief DGESVD Small path (path 6)
 *  without the LQ Factorization.
 *  */

#include "FLAME.h"
#include "fla_lapack_avx2_kernels.h"

#if FLA_ENABLE_AMD_OPT

/* SVD for small fat-matrices
 */
void fla_dgesvd_xs_small10T_avx2(aocl_int64_t *m, aocl_int64_t *n, doublereal *a, aocl_int64_t *lda,
                                 doublereal *s, doublereal *u, aocl_int64_t *ldu, doublereal *vt,
                                 aocl_int64_t *ldvt, doublereal *work, aocl_int64_t *info)
{
    /* Declare and init local variables */
    FLA_GEQRF_INIT_DSMALL();

    aocl_int64_t ie;
    aocl_int64_t itauq, itaup;
    aocl_int64_t rlen, knt;
    aocl_int64_t tm, tn;
    aocl_int64_t c__1 = 1;

    doublereal *tau, *tauq, *taup;
    doublereal *e;
    doublereal *iptr;
    doublereal stau, d__1;

    doublereal *ta, *ts;

    /* indices for partitioning work buffer */
    ie = 1;
    itauq = ie + *m;
    itaup = itauq + *m;

    /* parameter adjustments */
    a -= (1 + *lda);
    u -= (1 + *ldu);
    vt -= (1 + *ldvt);
    --s;
    --work;

    /* work buffer distribution */
    e = &work[ie - 1];
    tauq = &work[itauq - 1];
    taup = &work[itaup - 1];

    /* Lower Bidiagonalization */
    {
        /* Annihilate first row elements to the right of the diagonal */
        rlen = *n - 1;
        slen = *m - 1;
        iptr = a + 1;
        tau = taup;
        FLA_LARF_GEN_DSMALL_ROW(1, m, n, iptr, lda, tau);
        s[1] = beta;
        FLA_LARF_APPLY_SMALL_ROW(1, m, n, iptr, lda, tau);

        /* Upper Bidiagonalize the matrix excluding the first row */
        tm = *m - 1;
        tn = *n;
        ta = a + 1;
        tau = taup + 1;
        ts = s + 1;
        FLA_BIDIAGONALIZE_DSMALL(tm, tn, ta, lda, tauq, tau, e, ts);
    }

    /* Generate Qr (from bidiag) in vt */
    FLA_GESVD_FORM_VT_LOWER_SMALL(m, n, a, lda, taup, vt, ldvt);

    /* Generate Ql (from bidiag) in u from a */
    FLA_GESVD_FORM_U_LOWER_SMALL(m, a, lda, tauq, u, ldu);

    /* Compute Singular Values and Vetcors */
    lapack_dbdsqr_small("L", m, n, m, &s[1], &e[1], &vt[1 + *ldvt], ldvt, &u[1 + *ldu], ldu, info);

    return;
}
#endif
