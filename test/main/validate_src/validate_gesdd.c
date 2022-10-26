/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_gesdd.c
 *  @brief Defines validate function of GESDD() to use in test suite.
 *  */

#include "test_common.h"

void validate_gesdd(char *jobz, integer m, integer n, void* A, void* A_test, void* s, void* U, void* V, integer datatype, double *residual)
{
    void *sigma = NULL, *Usigma = NULL;
    void *work = NULL;

    create_matrix(datatype, &sigma, m, n);
    create_matrix(datatype, &Usigma, m, n);
    reset_matrix(datatype, m, n, Usigma, m);

    diagonalize_vector(datatype, s, sigma, m, n, m);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2, resid3;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (N * norm(A) * EPS)*/
            norm_A = slange_("1", &m, &n, A, &m, work);
            sgemm_("N", "N", &m, &n, &m, &s_one, U, &m, sigma, &m, &s_zero, Usigma, &m);
            sgemm_("N", "N", &m, &n, &n, &s_one, Usigma, &m, V, &n, &s_n_one, A, &m);
            norm = slange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, U, m, m);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            resid3 = (float)check_orthogonality(datatype, V, n, n);

            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2, resid3;
            eps = dlamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (N * norm(A) * EPS)*/
            norm_A = dlange_("1", &m, &n, A, &m, work);
            dgemm_("N", "N", &m, &n, &m, &d_one, U, &m, sigma, &m, &d_zero, Usigma, &m);
            dgemm_("N", "N", &m, &n, &n, &d_one, Usigma, &m, V, &n, &d_n_one, A, &m);
            norm = dlange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (double)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, U, m, m);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            resid3 = check_orthogonality(datatype, V, n, n);

            *residual = (double)max(resid1, max(resid2, resid3));
             break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2, resid3;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (N * norm(A) * EPS)*/
            norm_A = clange_("1", &m, &n, A, &m, work);
            cgemm_("N", "N", &m, &n, &m, &c_one, U, &m, sigma, &m, &c_zero, Usigma, &m);
            cgemm_("N", "N", &m, &n, &n, &c_one, Usigma, &m, V, &n, &c_n_one, A, &m);
            norm = clange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, U, m, m);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            resid3 = (float)check_orthogonality(datatype, V, n, n);

            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2, resid3;
            eps = dlamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (N * norm(A) * EPS)*/
            norm_A = zlange_("1", &m, &n, A, &m, work);
            zgemm_("N", "N", &m, &n, &m, &z_one, U, &m, sigma, &m, &z_zero, Usigma, &m);
            zgemm_("N", "N", &m, &n, &n, &z_one, Usigma, &m, V, &n, &z_n_one, A, &m);
            norm = zlange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (double)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, U, m, m);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            resid3 = check_orthogonality(datatype, V, n, n);

            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }
    }
    free_matrix(sigma);
    free_matrix(Usigma);
}
