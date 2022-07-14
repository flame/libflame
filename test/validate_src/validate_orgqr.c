/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_orgqr.c
 *  @brief Defines validate function of ORGQR() to use in test suite.
 *  */

#include "test_common.h"

void validate_orgqr(integer m,
    integer n,
    void *A,
    void* Q,
    void *R,
    void* work,
    integer datatype,
    double* residual)
{
    integer k;

    if( m == n)
        k = m;
    else
        k = n;

    switch( datatype )
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = slamch_("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            norm_A = slange_("1", &k, &k, R, &k, work);
            sgemm_("T", "N", &k, &k, &m, &s_n_one, Q, &m, A, &m, &s_one, R, &k);

            norm = slange_("1", &k, &k, R, &k, work);
            resid1 = norm/(eps * norm_A * (float)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m, k);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = dlamch_("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            norm_A = dlange_("1", &k, &k, R, &k, work);
            dgemm_("T", "N", &k, &k, &m, &d_n_one, Q, &m, A, &m, &d_one, R, &k);

            norm = dlange_("1", &k, &k, R, &k, work);
            resid1 = norm/(eps * norm_A * (double)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m, k);
   
            *residual = (double)max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = slamch_("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            norm_A = clange_("1", &k, &k, R, &k, work);
            cgemm_("C", "N", &k, &k, &m, &c_n_one, Q, &m, A, &m, &c_one, R, &k);

            norm = clange_("1", &k, &k, R, &k, work);
            resid1 = norm/(eps * norm_A * (double)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m, k);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = dlamch_("P");

            /* Test 1
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            norm_A = zlange_("1", &k, &k, R, &k, work);
            zgemm_("C", "N", &k, &k, &m, &z_n_one, Q, &m, A, &m, &z_one, R, &k);

            norm = zlange_("1", &k, &k, R, &k, work);
            resid1 = norm/(eps * norm_A * (double)k);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m, k);

            *residual = (double)max(resid1, resid2);
            break;
        }
    }
}
