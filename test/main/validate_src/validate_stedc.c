/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_stedc.c
 *  @brief Defines validate function of STEDC() to use in test suite.
 *  */

#include "test_common.h"

/* This function will validate STEDC() output eigenvectors and orthogonal
   matrices only if compz != N, as output will not be generated 
   if compz = N.*/
void validate_stedc(char compz, integer n, void* D_test, void* Z_input, void* Z, integer ldz, integer datatype, double* residual)
{
    void *lambda = NULL, *zlambda = NULL;
    void *I = NULL;
    void *work = NULL;

    if (compz == 'N') {
        *residual = 0.0;
        return;
    }
    create_matrix(datatype, &lambda, n, n);
    create_matrix(datatype, &zlambda, n, n);
    reset_matrix(datatype, n, n, zlambda, n);
    diagonalize_vector(datatype, D_test, lambda, n, n, n);

    create_matrix(datatype, &I, n, n);
    reset_matrix(datatype, n, n, I, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1 = 0, resid2 = 0;

            eps = slamch_("P");
            /* Test 1 - Check for Eigen vectors and Eigen values.
               compute norm(A - (Z * lambda * Z')) / (N * norm(A) * EPS)*/
            norm_A = slange_("1", &n, &n, Z_input, &ldz, work);
            sgemm_("N", "N", &n, &n, &n, &s_one, Z, &ldz, lambda, &n, &s_zero, zlambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, zlambda, &n, Z, &ldz, &s_n_one, Z_input, &ldz);
            norm = slange_("1", &n, &n, Z_input, &ldz, work);
            resid1 = norm/(eps * norm_A * (float)n);
            /* Test 2 - Check for orthogonality of matrix.
               compute norm(I - Z*Z') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Z, n, n, ldz);
            *residual = (double)max(resid1, resid2);
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps, resid1 = 0, resid2 = 0;

            eps = dlamch_("P");
            /* Test 1 - Check for Eigen vectors and Eigen values.
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = dlange_("1", &n, &n, Z_input, &ldz, work);
            dgemm_("N", "N", &n, &n, &n, &d_one, Z, &ldz, lambda, &n, &d_zero, zlambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, zlambda, &n, Z, &ldz, &d_n_one, Z_input, &ldz);
            norm = dlange_("1", &n, &n, Z_input, &ldz, work);
            resid1 = norm/(eps * norm_A * (float)n);
            /* Test 2 - Check for orthogonality of matrix.
               compute norm(I - Z*Z') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Z, n, n, ldz);
            *residual = (double)max(resid1, resid2);
            break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1 = 0, resid2 = 0;

            eps = slamch_("P");
            /* Test 1 - Check for Eigen vectors and Eigen values.
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = clange_("1", &n, &n, Z_input, &ldz, work);
            cgemm_("N", "N", &n, &n, &n, &c_one, Z, &ldz, lambda, &n, &c_zero, zlambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, zlambda, &n, Z, &ldz, &c_n_one, Z_input, &ldz);
            norm = clange_("1", &n, &n, Z_input, &ldz, work);
            resid1 = norm/(eps * norm_A * (float)n);
            /* Test 2 - Check for orthogonality of matrix.
               compute norm(I - Z*Z') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Z, n, n, ldz);
            *residual = (double)max(resid1, resid2);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1 = 0, resid2 = 0;

            eps = dlamch_("P");
            /* Test 1 - Check for Eigen vectors and Eigen values.
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = zlange_("1", &n, &n, Z_input, &ldz, work);
            zgemm_("N", "N", &n, &n, &n, &z_one, Z, &ldz, lambda, &n, &z_zero, zlambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, zlambda, &n, Z, &ldz, &z_n_one, Z_input, &ldz);
            norm = zlange_("1", &n, &n, Z_input, &ldz, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2 - Check for orthogonality of matrix.
               compute norm(I - Z*Z') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Z, n, n, ldz);
            *residual = (double)max(resid1, resid2);
            break;
        }
    }
    free_matrix(lambda);
    free_matrix(zlambda);
    free_matrix(I);
}