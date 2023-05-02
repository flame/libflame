/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_gesdd.c
 *  @brief Defines validate function of GESDD() to use in test suite.
 *  */

#include "test_common.h"

void validate_gesdd(char *jobz, integer m, integer n, void* A, void* A_test, integer lda, void* s, void* U, integer ldu, void* V, integer ldvt, integer datatype, double *residual, integer* info)
{
    if(m == 0 || n == 0)
      return;
    void *sigma = NULL, *Usigma = NULL;
    void *work = NULL;
    *info = 0;

    create_matrix(datatype, &sigma, m, n);
    create_matrix(datatype, &Usigma, m, n);
    reset_matrix(datatype, m, n, Usigma, m);

    diagonalize_vector(datatype, s, sigma, m, n, m);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2, resid3, resid4;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (N * norm(A) * EPS)*/
            norm_A = fla_lapack_slange("1", &m, &n, A, &lda, work);
            sgemm_("N", "N", &m, &n, &m, &s_one, U, &ldu, sigma, &m, &s_zero, Usigma, &m);
            sgemm_("N", "N", &m, &n, &n, &s_one, Usigma, &m, V, &ldvt, &s_n_one, A, &lda);
            norm = fla_lapack_slange("1", &m, &n, A, &lda, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, U, m, m, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            resid3 = (float)check_orthogonality(datatype, V, n, n, ldvt);

            /* Test 4
               Test to Check order of Singular values of SVD (positive and non-decreasing) */
            resid4 = (float)svd_check_order( datatype, s, m, n, *residual );

            *residual = (double)fla_max(fla_max(resid1, fla_max(resid2, resid3)), resid4);
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2, resid3, resid4;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (N * norm(A) * EPS)*/
            norm_A = fla_lapack_dlange("1", &m, &n, A, &lda, work);
            dgemm_("N", "N", &m, &n, &m, &d_one, U, &ldu, sigma, &m, &d_zero, Usigma, &m);
            dgemm_("N", "N", &m, &n, &n, &d_one, Usigma, &m, V, &ldvt, &d_n_one, A, &lda);
            norm = fla_lapack_dlange("1", &m, &n, A, &lda, work);
            resid1 = norm/(eps * norm_A * (double)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, U, m, m, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            resid3 = check_orthogonality(datatype, V, n, n, ldvt);
            
            /* Test 4
               Test to Check order of Singular values of SVD (positive and non-decreasing) */
            resid4 = svd_check_order( datatype, s, m, n, *residual );

            *residual = (double)fla_max(fla_max(resid1, fla_max(resid2, resid3)), resid4);
             break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2, resid3, resid4;
            eps = fla_lapack_slamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (N * norm(A) * EPS)*/
            norm_A = fla_lapack_clange("1", &m, &n, A, &lda, work);
            cgemm_("N", "N", &m, &n, &m, &c_one, U, &ldu, sigma, &m, &c_zero, Usigma, &m);
            cgemm_("N", "N", &m, &n, &n, &c_one, Usigma, &m, V, &ldvt, &c_n_one, A, &lda);
            norm = fla_lapack_clange("1", &m, &n, A, &lda, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, U, m, m, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            resid3 = (float)check_orthogonality(datatype, V, n, n, ldvt);

            /* Test 4
               Test to Check order of Singular values of SVD (positive and non-decreasing) */
            resid4 = (float)svd_check_order( datatype, s, m, n, *residual );

            *residual = (double)fla_max(fla_max(resid1, fla_max(resid2, resid3)), resid4);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2, resid3, resid4;
            eps = fla_lapack_dlamch("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (N * norm(A) * EPS)*/
            norm_A = fla_lapack_zlange("1", &m, &n, A, &lda, work);
            zgemm_("N", "N", &m, &n, &m, &z_one, U, &ldu, sigma, &m, &z_zero, Usigma, &m);
            zgemm_("N", "N", &m, &n, &n, &z_one, Usigma, &m, V, &ldvt, &z_n_one, A, &lda);
            norm = fla_lapack_zlange("1", &m, &n, A, &lda, work);
            resid1 = norm/(eps * norm_A * (double)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, U, m, m, ldu);

            /* Test 3
               compute norm(I - V*V') / (N * EPS)*/
            resid3 = check_orthogonality(datatype, V, n, n, ldvt);

           /* Test 4
              Test to Check order of Singular values of SVD  (positive and non-decreasing) */
            resid4 = svd_check_order( datatype, s, m, n, *residual );

            *residual = (double)fla_max(fla_max(resid1, fla_max(resid2, resid3)), resid4);
            break;
        }
    }
    free_matrix(sigma);
    free_matrix(Usigma);
}
