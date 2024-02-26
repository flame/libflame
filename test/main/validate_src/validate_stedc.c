/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_stedc.c
 *  @brief Defines validate function of STEDC() to use in test suite.
 *  */

#include "test_common.h"

/* This function will validate STEDC() output eigenvectors and orthogonal
   matrices only if compz != N, as output will not be generated 
   if compz = N.*/
void validate_stedc(char compz, integer n, void* D_test, void* Z_input, void* Z, integer ldz, integer datatype, double* residual, integer* info)
{
   if(n == 0)
        return;
    void *lambda = NULL, *zlambda = NULL;
    void *a_temp = NULL;
    void *work = NULL;
    *info = 0;

    if (compz == 'N')
        return;

    create_matrix(datatype, &lambda, n, n);
    create_matrix(datatype, &zlambda, n, n);
    reset_matrix(datatype, n, n, zlambda, n);
    diagonalize_vector(datatype, D_test, lambda, n, n, n);

    create_matrix(datatype, &a_temp, n, n);
    reset_matrix(datatype, n, n, a_temp, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1 = 0, resid2 = 0;

            eps = fla_lapack_slamch("P");
            /* Test 1 - Check for Eigen vectors and Eigen values.
               compute norm(A - (Z * lambda * Z')) / (N * norm(A) * EPS)*/
            norm_A = fla_lapack_slange("1", &n, &n, Z_input, &ldz, work);
            sgemm_("N", "N", &n, &n, &n, &s_one, Z, &ldz, lambda, &n, &s_zero, zlambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, zlambda, &n, Z, &ldz, &s_n_one, Z_input, &ldz);
            norm = fla_lapack_slange("1", &n, &n, Z_input, &ldz, work);
            resid1 = norm/(eps * norm_A * (float)n);
            /* Test 2 - Check for orthogonality of matrix.
               compute norm(I - Z*Z') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Z, n, n, ldz);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps, resid1 = 0, resid2 = 0;

            eps = fla_lapack_dlamch("P");
            /* Test 1 - Check for Eigen vectors and Eigen values.
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = fla_lapack_dlange("1", &n, &n, Z_input, &ldz, work);
            dgemm_("N", "N", &n, &n, &n, &d_one, Z, &ldz, lambda, &n, &d_zero, zlambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, zlambda, &n, Z, &ldz, &d_n_one, Z_input, &ldz);
            norm = fla_lapack_dlange("1", &n, &n, Z_input, &ldz, work);
            resid1 = norm/(eps * norm_A * (float)n);
            /* Test 2 - Check for orthogonality of matrix.
               compute norm(I - Z*Z') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Z, n, n, ldz);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1 = 0, resid2 = 0;

            eps = fla_lapack_slamch("P");
            /* Test 1 - Check for Eigen vectors and Eigen values.
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = fla_lapack_clange("1", &n, &n, Z_input, &ldz, work);
            cgemm_("N", "N", &n, &n, &n, &c_one, Z, &ldz, lambda, &n, &c_zero, zlambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, zlambda, &n, Z, &ldz, &c_n_one, Z_input, &ldz);
            norm = fla_lapack_clange("1", &n, &n, Z_input, &ldz, work);
            resid1 = norm/(eps * norm_A * (float)n);
            /* Test 2 - Check for orthogonality of matrix.
               compute norm(I - Z*Z') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Z, n, n, ldz);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1 = 0, resid2 = 0;

            eps = fla_lapack_dlamch("P");
            /* Test 1 - Check for Eigen vectors and Eigen values.
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = fla_lapack_zlange("1", &n, &n, Z_input, &ldz, work);
            zgemm_("N", "N", &n, &n, &n, &z_one, Z, &ldz, lambda, &n, &z_zero, zlambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, zlambda, &n, Z, &ldz, &z_n_one, Z_input, &ldz);
            norm = fla_lapack_zlange("1", &n, &n, Z_input, &ldz, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2 - Check for orthogonality of matrix.
               compute norm(I - Z*Z') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Z, n, n, ldz);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }
    free_matrix(lambda);
    free_matrix(zlambda);
    free_matrix(a_temp);
}
