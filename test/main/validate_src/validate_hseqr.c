/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_hseqr.c
 *  @brief Defines validate function of HSEQR() to use in test suite.
 *  */

#include "test_common.h"

void validate_hseqr(char* job, char* compz,
    integer n,
    void* H,
    void* H_test,
    integer ldh,
    void* Z,
    void* Z_test,
    integer ldz,
    integer datatype,
    double* residual,
    integer *info)
{
    if(n == 0)
        return;
    if (*job == 'E' || *compz == 'N')
        return;

    void *zlambda = NULL, *work = NULL, *lambda = NULL;
    *info = 0;

    create_matrix(datatype, &zlambda, n, n);
    reset_matrix(datatype, n, n, zlambda, n);
    create_matrix(datatype, &lambda, n, n);
    reset_matrix(datatype, n, n, lambda, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_H, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            norm_H = fla_lapack_slange("1", &n, &n, H, &ldh, work);
            sgemm_("T", "N", &n, &n, &n, &s_one, Z, &ldz, Z_test, &ldz, &s_zero, lambda, &n);
            sgemm_("N", "N", &n, &n, &n, &s_one, lambda, &n, H_test, &ldh, &s_zero, zlambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, zlambda, &n, lambda, &n, &s_n_one, H, &ldh);
            norm = fla_lapack_slange("1", &n, &n, H, &ldh, work);
            resid1 = norm /( eps * norm_H * (float)n);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, lambda, n, n, n);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_H, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            norm_H = fla_lapack_dlange("1", &n, &n, H, &ldh, work);
            dgemm_("T", "N", &n, &n, &n, &d_one, Z, &ldz, Z_test, &ldz, &d_zero, lambda, &n);
            dgemm_("N", "N", &n, &n, &n, &d_one, lambda, &n, H_test, &ldh, &d_zero, zlambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, zlambda, &n, lambda, &n, &d_n_one, H, &ldh);
            norm = fla_lapack_dlange("1", &n, &n, H, &ldh, work);
            resid1 = norm/ ( norm_H * (float)n * eps);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid2 = check_orthogonality(datatype, lambda, n, n, n);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_H, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            norm_H = fla_lapack_clange("1", &n, &n, H, &ldh, work);
            cgemm_("C", "N", &n, &n, &n, &c_one, Z, &ldz, Z_test, &ldz, &c_zero, lambda, &n);
            cgemm_("N", "N", &n, &n, &n, &c_one, lambda, &n, H_test, &ldh, &c_zero, zlambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, zlambda, &n, lambda, &n, &c_n_one, H, &ldh);
            norm = fla_lapack_clange("1", &n, &n, H, &ldh, work);
            resid1 = norm/( norm_H * (float)n * eps);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, lambda, n, n, n);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_H, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            /* Test 1
                lambda = Z' * ZU
                compute norm(H - (lambda * T * lambda')) / (V * norm(H) * eps)*/
            norm_H = fla_lapack_zlange("1", &n, &n, H, &ldh, work);
            zgemm_("C", "N", &n, &n, &n, &z_one, Z, &ldz, Z_test, &ldz, &z_zero, lambda, &n);
            zgemm_("N", "N", &n, &n, &n, &z_one, lambda, &n, H_test, &ldh, &z_zero, zlambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, zlambda, &n, lambda, &n, &z_n_one, H, &ldh);
            norm = fla_lapack_zlange("1", &n, &n, H, &ldh, work);
            resid1 = norm/( norm_H * (float)n * eps);

            /* Test 2
                compute norm(I - Z'*Z) / (N * EPS)*/
            resid2 = check_orthogonality(datatype, lambda, n, n, n);
            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }
    free_matrix(zlambda);
    free_matrix(lambda);
}
