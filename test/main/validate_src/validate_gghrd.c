/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_gghrd.c
 *  @brief Defines validate function of GGHRD() to use in test suite.
 *  */

#include "test_common.h"

void validate_gghrd(char* compq,
    char* compz,
    integer n,
    void* A,
    void* A_test,
    integer lda,
    void* B,
    void* B_test,
    integer ldb,
    void* Q,
    void* Q_test,
    integer ldq,
    void* Z,
    void* Z_test,
    integer ldz,
    integer datatype,
    double* residual,
    integer *info)
{
    if(n == 0)
        return;
    if (*compz == 'N' || *compq == 'N')
        return;

    void *work = NULL, *lambda = NULL, *alambda = NULL;
    *info = 0;

    create_matrix(datatype, &lambda, n, n);
    create_matrix(datatype, &alambda, n, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, norm_T, eps, resid1, resid2, resid3, resid4;
            double res_max, res_max1;
            eps = fla_lapack_slamch("P");

            /* Test 1
                | A - Q H Z**T  | / ( |A| n ulp ) */
            norm_A = fla_lapack_slange("1", &n, &n, A, &lda, work);
            sgemm_("N", "N", &n, &n, &n, &s_one, Q_test, &ldq, A_test, &lda, &s_zero, lambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Z_test, &ldz, &s_zero, alambda, &n);
            sgemm_("T", "N", &n, &n, &n, &s_one, Q, &ldq, alambda, &n, &s_zero, lambda, &n);
            sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, Z, &ldz, &s_n_one, A, &lda);
            norm = fla_lapack_slange("1", &n, &n, A, &lda, work);
            resid1 = norm /( eps * norm_A * (float)n);

            /* Test 2
                | B - Q T Z**T  | / ( |B| n ulp ) */
            norm_T = fla_lapack_slange("1", &n, &n, B, &ldb, work);
            sgemm_("N", "N", &n, &n, &n, &s_one, Q_test, &ldq, B_test, &ldb, &s_zero, lambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Z_test, &ldz, &s_zero, alambda, &n);
            sgemm_("T", "N", &n, &n, &n, &s_one, Q, &ldq, alambda, &n, &s_zero, lambda, &n);
            sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, Z, &ldz, &s_n_one, B, &ldb);
            norm = fla_lapack_slange("1", &n, &n, B, &ldb, work);
            resid2 = norm /( eps * norm_T * (float)n);

            /* Test 3
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid3 = (float)check_orthogonality(datatype, Z, n, n, ldz);

            /* Test 4
                compute norm(I - Q'*Q) / (N * EPS)*/
            resid4 = (float)check_orthogonality(datatype, Q, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, norm_B, eps, resid1, resid2, resid3, resid4;
            double res_max, res_max1;
            eps = fla_lapack_dlamch("P");

            /* Test 1
                | A - Q H Z**T  | / ( |A| n ulp ) */
            norm_A = fla_lapack_dlange("1", &n, &n, A, &lda, work);
            dgemm_("N", "N", &n, &n, &n, &d_one, Q_test, &ldq, A_test, &lda, &d_zero, lambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Z_test, &ldz, &d_zero, alambda, &n);
            dgemm_("T", "N", &n, &n, &n, &d_one, Q, &ldq, alambda, &n, &d_zero, lambda, &n);
            dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, Z, &ldz, &d_n_one, A, &lda);
            norm = fla_lapack_dlange("1", &n, &n, A, &lda, work);
            resid1 = norm /( eps * norm_A * (float)n);

            /* Test 2
                | B - Q T Z**T  | / ( |B| n ulp ) */
            norm_B = fla_lapack_dlange("1", &n, &n, B, &ldb, work);
            dgemm_("N", "N", &n, &n, &n, &d_one, Q_test, &ldq, B_test, &ldb, &d_zero, lambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Z_test, &ldz, &d_zero, alambda, &n);
            dgemm_("T", "N", &n, &n, &n, &d_one, Q, &ldq, alambda, &n, &d_zero, lambda, &n);
            dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, Z, &ldz, &d_n_one, B, &ldb);
            norm = fla_lapack_dlange("1", &n, &n, B, &ldb, work);
            resid2 = norm /( eps * norm_B * (float)n);

            /* Test 3
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid3 = check_orthogonality(datatype, Z, n, n, ldz);

            /* Test 4
                compute norm(I - Q'*Q) / (N * EPS)*/
            resid4 = check_orthogonality(datatype, Q, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            break;
        }        
        case COMPLEX:
        {
            float norm, norm_A, norm_B, eps, resid1, resid2, resid3, resid4;
            double res_max, res_max1;
            eps = fla_lapack_slamch("P");

            /* Test 1
                | A - Q H Z**T  | / ( |A| n ulp ) */
            norm_A = fla_lapack_clange("1", &n, &n, A, &lda, work);
            cgemm_("N", "N", &n, &n, &n, &c_one, Q_test, &ldq, A_test, &lda, &c_zero, lambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Z_test, &ldz, &c_zero, alambda, &n);
            cgemm_("C", "N", &n, &n, &n, &c_one, Q, &ldq, alambda, &n, &c_zero, lambda, &n);
            cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, Z, &ldz, &c_n_one, A, &lda);
            norm = fla_lapack_clange("1", &n, &n, A, &lda, work);
            resid1 = norm /( eps * norm_A * (float)n);
            

            /* Test 2
                | B - Q T Z**T  | / ( |B| n ulp ) */
            norm_B = fla_lapack_clange("1", &n, &n, B, &ldb, work);
            cgemm_("N", "N", &n, &n, &n, &c_one, Q_test, &ldq, B_test, &ldb, &c_zero, lambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Z_test, &ldz, &c_zero, alambda, &n);
            cgemm_("C", "N", &n, &n, &n, &c_one, Q, &ldq, alambda, &n, &c_zero, lambda, &n);
            cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, Z, &ldz, &c_n_one, B, &ldb);
            norm = fla_lapack_clange("1", &n, &n, B, &ldb, work);
            resid2 = norm /( eps * norm_B * (float)n);

            /* Test 3
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid3 = (float)check_orthogonality(datatype, Z, n, n, ldz);

            /* Test 4
                compute norm(I - Q'*Q) / (N * EPS)*/
            resid4 = (float)check_orthogonality(datatype, Q, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            break;   
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, norm_B, eps, resid1, resid2, resid3, resid4;
            double res_max, res_max1;
            eps = fla_lapack_dlamch("P");

            /* Test 1
                | A - Q H Z**T  | / ( |A| n ulp ) */
            norm_A = fla_lapack_zlange("1", &n, &n, A, &lda, work);
            zgemm_("N", "N", &n, &n, &n, &z_one, Q_test, &ldq, A_test, &lda, &z_zero, lambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Z_test, &ldz, &z_zero, alambda, &n);
            zgemm_("C", "N", &n, &n, &n, &z_one, Q, &ldq, alambda, &n, &z_zero, lambda, &n);
            zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, Z, &ldz, &z_n_one, A, &lda);
            norm = fla_lapack_zlange("1", &n, &n, A, &lda, work);
            resid1 = norm /( eps * norm_A * (float)n);

            /* Test 2
                | B - Q T Z**T  | / ( |B| n ulp ) */
            norm_B = fla_lapack_zlange("1", &n, &n, B, &ldb, work);
            zgemm_("N", "N", &n, &n, &n, &z_one, Q_test, &ldq, B_test, &ldb, &z_zero, lambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Z_test, &ldz, &z_zero, alambda, &n);
            zgemm_("C", "N", &n, &n, &n, &z_one, Q, &ldq, alambda, &n, &z_zero, lambda, &n);
            zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, Z, &ldz, &z_n_one, B, &ldb);
            norm = fla_lapack_zlange("1", &n, &n, B, &ldb, work);
            resid2 = norm /( eps * norm_B * (float)n);

            /* Test 3
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid3 = check_orthogonality(datatype, Z, n, n, ldz);

            /* Test 4
                compute norm(I - Q'*Q) / (N * EPS)*/
            resid4 = check_orthogonality(datatype, Q, n, n, ldq);

            res_max = (double)fla_max(resid1, resid2);
            res_max1 = (double)fla_max(resid3, resid4);
            *residual = (double)fla_max(res_max, res_max1);
            break;
        }
    }
    free_matrix(lambda);
    free_matrix(alambda);
}
