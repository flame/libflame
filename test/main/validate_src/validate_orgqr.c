/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_orgqr.c
 *  @brief Defines validate function of ORGQR() to use in test suite.
 *  */

#include "test_common.h"

void validate_orgqr(integer m,
    integer n,
    void *A,
    integer lda,
    void* Q,
    void *R,
    void* work,
    integer datatype,
    double* residual,
    integer* info)
{
    integer k;
    *info = 0;

    if( m == n)
        k = m;
    else
        k = n;

    switch( datatype )
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            norm_A = fla_lapack_slange("1", &k, &k, R, &k, work);
            sgemm_("T", "N", &k, &k, &m, &s_n_one, Q, &m, A, &lda, &s_one, R, &k);

            norm = fla_lapack_slange("1", &k, &k, R, &k, work);
            resid1 = norm/(eps * norm_A * (float)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m, k, m);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            norm_A = fla_lapack_dlange("1", &k, &k, R, &k, work);
            dgemm_("T", "N", &k, &k, &m, &d_n_one, Q, &m, A, &lda, &d_one, R, &k);

            norm = fla_lapack_dlange("1", &k, &k, R, &k, work);
            resid1 = norm/(eps * norm_A * (double)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m, k, m);
   
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            norm_A = fla_lapack_clange("1", &k, &k, R, &k, work);
            cgemm_("C", "N", &k, &k, &m, &c_n_one, Q, &m, A, &lda, &c_one, R, &k);

            norm = fla_lapack_clange("1", &k, &k, R, &k, work);
            resid1 = norm/(eps * norm_A * (double)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m, k, m);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            norm_A = fla_lapack_zlange("1", &k, &k, R, &k, work);
            zgemm_("C", "N", &k, &k, &m, &z_n_one, Q, &m, A, &lda, &z_one, R, &k);

            norm = fla_lapack_zlange("1", &k, &k, R, &k, work);
            resid1 = norm/(eps * norm_A * (double)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m, k, m);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }
}
