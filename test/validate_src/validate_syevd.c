/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_syevd.c
 *  @brief Defines validate function of SYEVD() to use in test suite.
 *  */

#include "test_common.h"

void validate_syevd(char* jobz, char* uplo, integer n, void* A, void* A_test, void* w, integer datatype, double* residual)
{
    void *lambda = NULL, *zlambda = NULL;
    void *I = NULL, *Z = NULL;
    void *work = NULL;

    create_matrix(datatype, &lambda, n, n);
    create_matrix(datatype, &zlambda, n, n);
    create_matrix(datatype, &Z, n, n);

    reset_matrix(datatype, n, n, zlambda, n);
    reset_matrix(datatype, n, n, Z, n);

    create_matrix(datatype, &I, n, n);
    reset_matrix(datatype, n, n, I, n);

    copy_matrix(datatype, "full", n, n, A_test, n, Z, n);

    diagonalize_vector(datatype, w, lambda, n, n, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = slange_("1", &n, &n, A, &n, work);
            sgemm_("N", "N", &n, &n, &n, &s_one, Z, &n, lambda, &n, &s_zero, zlambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, zlambda, &n, Z, &n, &s_n_one, A, &n);
            norm = slange_("1", &n, &n, A, &n, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Z*Z') / (N * EPS)*/
            slaset_("full", &n, &n, &s_zero, &s_one, I, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, Z, &n, Z, &n, &s_n_one, I, &n);
            norm = slange_("1", &n, &n, I, &n, work);
            resid2 = norm/(eps * (float)n);
            *residual = (double)max(resid1, resid2);
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = dlamch_("P");

            /* Test 1
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = dlange_("1", &n, &n, A, &n, work);
            dgemm_("N", "N", &n, &n, &n, &d_one, Z, &n, lambda, &n, &d_zero, zlambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, zlambda, &n, Z, &n, &d_n_one, A, &n);
            norm = dlange_("1", &n, &n, A, &n, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Z*Z') / (N * EPS)*/
            dlaset_("full", &n, &n, &d_zero, &d_one, I, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, Z, &n, Z, &n, &d_n_one, I, &n);
            norm = dlange_("1", &n, &n, I, &n, work);
            resid2 = norm/(eps * (float)n);
            *residual = (double)max(resid1, resid2);
            break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = clange_("1", &n, &n, A, &n, work);
            cgemm_("N", "N", &n, &n, &n, &c_one, Z, &n, lambda, &n, &c_zero, zlambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, zlambda, &n, Z, &n, &c_n_one, A, &n);
            norm = clange_("1", &n, &n, A, &n, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Z*Z') / (N * EPS)*/
            claset_("full", &n, &n, &c_zero, &c_one, I, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, Z, &n, Z, &n, &c_n_one, I, &n);
            norm = clange_("1", &n, &n, I, &n, work);
            resid2 = norm/(eps * (float)n);
            *residual = (double)max(resid1, resid2);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = dlamch_("P");

          /* Test 1
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = zlange_("1", &n, &n, A, &n, work);
            zgemm_("N", "N", &n, &n, &n, &z_one, Z, &n, lambda, &n, &z_zero, zlambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, zlambda, &n, Z, &n, &z_n_one, A, &n);
            norm = zlange_("1", &n, &n, A, &n, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Z*Z') / (N * EPS)*/
            zlaset_("full", &n, &n, &z_zero, &z_one, I, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, Z, &n, Z, &n, &z_n_one, I, &n);
            norm = zlange_("1", &n, &n, I, &n, work);
            resid2 = norm/(eps * (float)n);
            *residual = (double)max(resid1, resid2);
            break;
        }
    }
    free_matrix(lambda);
    free_matrix(zlambda);
    free_matrix(I);
    free_matrix(Z);
}