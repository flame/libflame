/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_gehrd.c
 *  @brief Defines validate function of GEHRD() to use in test suite.
 *  */

#include "test_common.h"

void validate_gehrd(integer n,
    integer ilo,
    integer ihi,
    void* A,
    void* A_test,
    integer lda,
    void* tau,
    integer datatype,
    double* residual,
    integer *info)
{
    if(n == 0)
        return;
    void *Q = NULL, *work = NULL, *lambda = NULL;
    integer lwork = -1;
    *info = 0;
    create_matrix(datatype, &Q, lda, n);
    create_matrix(datatype, &lambda, n, n);
    copy_matrix(datatype, "full", n, n, A_test, lda, Q, lda);
    create_vector(datatype, &work, 1);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");

            /* Test 1 
                | A - Q H Q**T | / ( |A| n ulp ) */

            norm_A = fla_lapack_slange("1", &n, &n, A, &lda, work);
            fla_lapack_sorghr(&n, &ilo, &ihi, NULL, &lda, NULL, work, &lwork, info);
            if(*info == 0)
            {
                lwork = get_work_value(datatype, work);
                free_vector(work);
            }
            else
                break;
            create_vector(datatype, &work, lwork);
            fla_lapack_sorghr(&n, &ilo, &ihi, Q, &lda, tau, work, &lwork, info);
            if(*info < 0)
                break;
            extract_upper_hessenberg_matrix(datatype, n, A_test, lda);
            sgemm_("N", "N", &n, &n, &n, &s_one, Q, &lda, A_test, &lda, &s_zero, lambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, lambda, &n, Q, &lda, &s_n_one, A, &lda);
            norm = fla_lapack_slange("1", &n, &n, A, &lda, work);
            resid1 = norm /( eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, n, n, lda);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            /* Test 1 
                | A - Q H Q**T | / ( |A| n ulp ) */

            norm_A = fla_lapack_dlange("1", &n, &n, A, &lda, work);
            fla_lapack_dorghr(&n, &ilo, &ihi, NULL, &lda, NULL, work, &lwork, info);
            if(*info == 0)
            {
                lwork = get_work_value(datatype, work);
                free_vector(work);
            }
            else
                break;
            create_vector(datatype, &work, lwork);
            fla_lapack_dorghr(&n, &ilo, &ihi, Q, &lda, tau, work, &lwork, info);
            if(*info < 0)
                break;
            extract_upper_hessenberg_matrix(datatype, n, A_test, lda);
            dgemm_("N", "N", &n, &n, &n, &d_one, Q, &lda, A_test, &lda, &d_zero, lambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, lambda, &n, Q, &lda, &d_n_one, A, &lda);
            norm = fla_lapack_dlange("1", &n, &n, A, &lda, work);
            resid1 = norm /( eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n, n, lda);

            *residual = (double)fla_max(resid1, resid2);

            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_slamch("P");

            /* Test 1 
                | A - Q H Q**T | / ( |A| n ulp ) */

            norm_A = fla_lapack_clange("1", &n, &n, A, &lda, work);
            fla_lapack_cunghr(&n, &ilo, &ihi, NULL, &lda, NULL, work, &lwork, info);
            if(*info == 0)
            {
                lwork = get_work_value(datatype,work);
                free_vector(work);
            }
            else
                break;
            create_vector(datatype, &work, lwork);
            fla_lapack_cunghr(&n, &ilo, &ihi, Q, &lda, tau, work, &lwork, info);
            if(*info < 0)
                break;
            extract_upper_hessenberg_matrix(datatype, n, A_test, lda);
            cgemm_("N", "N", &n, &n, &n, &c_one, Q, &lda, A_test, &lda, &c_zero, lambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, lambda, &n, Q, &lda, &c_n_one, A, &lda);
            norm = fla_lapack_clange("1", &n, &n, A, &lda, work);
            resid1 = norm /( eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, n, n, lda);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = fla_lapack_dlamch("P");

            /* Test 1 
                | A - Q H Q**T | / ( |A| n ulp ) */

            norm_A = fla_lapack_zlange("1", &n, &n, A, &lda, work);
            fla_lapack_zunghr(&n, &ilo, &ihi, NULL, &lda, NULL, work, &lwork, info);
            if(*info == 0)
            {
                lwork = get_work_value(datatype,work);
                free_vector(work);
            }
            else
                break;
            create_vector(datatype, &work, lwork);
            fla_lapack_zunghr(&n, &ilo, &ihi, Q, &lda, tau, work, &lwork, info);
            if(*info < 0)
                break;
            extract_upper_hessenberg_matrix(datatype, n, A_test, lda);
            zgemm_("N", "N", &n, &n, &n, &z_one, Q, &lda, A_test, &lda, &z_zero, lambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, lambda, &n, Q, &lda, &z_n_one, A, &lda);
            norm = fla_lapack_zlange("1", &n, &n, A, &lda, work);
            resid1 = norm /( eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n, n, lda);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }
    free_vector(work);
    free_matrix(Q);
    free_matrix(lambda);
}
