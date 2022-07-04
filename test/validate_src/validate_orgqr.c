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
    void *I = NULL;

    /* Create Identity matrix to validate orthogonal property of matrix Q */
    create_matrix(datatype, &I, m, m);

    switch( datatype )
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2;

            eps = slamch_("P");

            /* Test 1
               compute norm(I - Q*Q') / (N * EPS)*/
            slaset_("full", &m, &m, &s_zero, &s_one, I, &m);
            sgemm_("N", "T", &m, &m, &m, &s_n_one, Q, &m, Q, &m, &s_one, I, &m);

            norm = slange_("1", &m, &m, I, &m, work);
            resid1 = norm/(eps * (float)n);

                        /* Test 2
                           compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
                        sgemm_("T", "N", &m, &n, &m, &s_n_one, Q, &m, A, &m, &s_one, R, &m);

                        norm_A = slange_("1", &m, &n, A, &m, work);
                        norm = slange_("1", &m, &n, R, &m, work);

                        resid2 = norm/(eps * norm_A * (float)n);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2;

            eps = dlamch_("P");

            /* Test 1
               compute norm(I - Q*Q') / (N * EPS)*/
            dlaset_("full", &m, &m, &d_zero, &d_one, I, &m);
            dgemm_("N", "T", &m, &m, &m, &d_n_one, Q, &m, Q, &m, &d_one, I, &m);

            norm = dlange_("1", &m, &m, I, &m, work);
            resid1 = norm/(eps * (double)n);

            /* Test 2
                           compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
                        dgemm_("T", "N", &m, &n, &m, &d_n_one, Q, &m, A, &m, &d_one, R, &m);

                        norm_A = dlange_("1", &m, &n, A, &m, work);
                        norm = dlange_("1", &m, &n, R, &m, work);

                        resid2 = norm/(eps * norm_A * (double)n);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2;

            eps = slamch_("P");

            /* Test 1
               compute norm(I - Q*Q') / (N * EPS)*/
            claset_("full", &m, &m, &c_zero, &c_one, I, &m);
            cgemm_("N", "C", &m, &m, &m, &c_n_one, Q, &m, Q, &m, &c_one, I, &m);

            norm = clange_("1", &m, &m, I, &m, work);
            resid1 = norm/(eps * (float)n);

             /* Test 2
                           compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
                        cgemm_("C", "N", &m, &n, &m, &c_n_one, Q, &m, A, &m, &c_one, R, &m);

                        norm_A = clange_("1", &m, &n, A, &m, work);
                        norm = clange_("1", &m, &n, R, &m, work);

                        resid2 = norm/(eps * norm_A * (float)n);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2;

            eps = dlamch_("P");

            /* Test 1
               compute norm(I - Q*Q') / (N * EPS)*/
            zlaset_("full", &m, &m, &z_zero, &z_one, I, &m);
            zgemm_("N", "C", &m, &m, &m, &z_n_one, Q, &m, Q, &m, &z_one, I, &m);

            norm = zlange_("1", &m, &m, I, &m, work);
            resid1 = norm/(eps * (double)n);

                        /* Test 2
                           compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
                        zgemm_("C", "N", &m, &n, &m, &z_n_one, Q, &m, A, &m, &z_one, R, &m);

                        norm_A = zlange_("1", &m, &n, A, &m, work);
                        norm = zlange_("1", &m, &n, R, &m, work);

                        resid2 = norm/(eps * norm_A * (double)n);

            *residual = (double)max(resid1, resid2);

            break;
        }
    }

    /* Free up buffers */
    free_matrix( I );
}