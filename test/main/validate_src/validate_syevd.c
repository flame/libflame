/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_syevd.c
 *  @brief Defines validate function of SYEVD() to use in test suite.
 *  */

#include "test_common.h"

void validate_syevd(char* jobz, integer n, void* A, void* A_test, integer lda, void* w, integer datatype, double* residual, integer* info)
{
    *info = 0;
    if(*jobz != 'N')
    {
        void *lambda = NULL, *zlambda = NULL, *Z = NULL;
        void *work = NULL;

        create_matrix(datatype, &lambda, n, n);
        create_matrix(datatype, &zlambda, n, n);
        create_matrix(datatype, &Z, lda, n);

        reset_matrix(datatype, n, n, zlambda, n);
        reset_matrix(datatype, n, n, Z, lda);

        copy_matrix(datatype, "full", n, n, A_test, lda, Z, lda);

        diagonalize_vector(datatype, w, lambda, n, n, n);

        switch(datatype)
        {
            case FLOAT:
            {
                float norm, norm_A, eps, resid1, resid2;
                eps = fla_lapack_slamch("P");

                /* Test 1
                   compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
                norm_A = fla_lapack_slange("1", &n, &n, A, &lda, work);
                sgemm_("N", "N", &n, &n, &n, &s_one, Z, &lda, lambda, &n, &s_zero, zlambda, &n);
                sgemm_("N", "T", &n, &n, &n, &s_one, zlambda, &n, Z, &lda, &s_n_one, A, &lda);
                norm = fla_lapack_slange("1", &n, &n, A, &lda, work);
                resid1 = norm/(eps * norm_A * (float)n);

                /* Test 2
                   compute norm(I - Z'*Z) / (N * EPS)*/
                resid2 = (float)check_orthogonality(datatype, Z, n, n, lda);

                *residual = (double)fla_max(resid1, resid2);
                break;
            }

            case DOUBLE:
            {
                double norm, norm_A, eps, resid1, resid2;
                eps = fla_lapack_dlamch("P");

                /* Test 1
                   compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
                norm_A = fla_lapack_dlange("1", &n, &n, A, &lda, work);
                dgemm_("N", "N", &n, &n, &n, &d_one, Z, &lda, lambda, &n, &d_zero, zlambda, &n);
                dgemm_("N", "T", &n, &n, &n, &d_one, zlambda, &n, Z, &lda, &d_n_one, A, &lda);
                norm = fla_lapack_dlange("1", &n, &n, A, &lda, work);
                resid1 = norm/(eps * norm_A * (float)n);

                /* Test 2
                   compute norm(I - Z'*Z) / (N * EPS)*/
                resid2 = check_orthogonality(datatype, Z, n, n, lda);

                *residual = (double)fla_max(resid1, resid2);
                break;
            }

            case COMPLEX:
            {
                float norm, norm_A, eps, resid1, resid2;
                eps = fla_lapack_slamch("P");

                /* Test 1
                   compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
                norm_A = fla_lapack_clange("1", &n, &n, A, &lda, work);
                cgemm_("N", "N", &n, &n, &n, &c_one, Z, &lda, lambda, &n, &c_zero, zlambda, &n);
                cgemm_("N", "C", &n, &n, &n, &c_one, zlambda, &n, Z, &lda, &c_n_one, A, &lda);
                norm = fla_lapack_clange("1", &n, &n, A, &lda, work);
                resid1 = norm/(eps * norm_A * (float)n);

                /* Test 2
                   compute norm(I - Z'*Z) / (N * EPS)*/
                resid2 = (float)check_orthogonality(datatype, Z, n, n, lda);

                *residual = (double)fla_max(resid1, resid2);
                break;
            }

            case DOUBLE_COMPLEX:
            {
                double norm, norm_A, eps, resid1, resid2;
                eps = fla_lapack_dlamch("P");

                /* Test 1
                   compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
                norm_A = fla_lapack_zlange("1", &n, &n, A, &lda, work);
                zgemm_("N", "N", &n, &n, &n, &z_one, Z, &lda, lambda, &n, &z_zero, zlambda, &n);
                zgemm_("N", "C", &n, &n, &n, &z_one, zlambda, &n, Z, &lda, &z_n_one, A, &lda);
                norm = fla_lapack_zlange("1", &n, &n, A, &lda, work);
                resid1 = norm/(eps * norm_A * (float)n);

                /* Test 2
                   compute norm(I - Z'*Z) / (N * EPS)*/
                resid2 = check_orthogonality(datatype, Z, n, n, lda);

                *residual = (double)fla_max(resid1, resid2);
                break;
            }
        }
        free_matrix(lambda);
        free_matrix(zlambda);
        free_matrix(Z);
    }
}
