/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_gelqf.c
 *  @brief Defines validate function of GELQF() to use in test suite.
 *  */

#include "test_common.h"

void validate_gelqf(integer m_A,
    integer n_A,
    void *A,
    void *A_test,
    integer lda,
    void *T_test,
    integer datatype,
    double* residual,
    integer* info)
{
    void *Q = NULL, *L = NULL, *work = NULL;
    integer min_A;
    integer lwork = -1;
    *info = 0;

    min_A = fla_min(m_A, n_A);

    /* Create Q and L matrices. */
    create_matrix(datatype, &Q, n_A, n_A);
    create_matrix(datatype, &L, m_A, n_A);
    reset_matrix(datatype, n_A, n_A, Q, n_A);
    reset_matrix(datatype, m_A, n_A, L, m_A);

    /* Extract L matrix and elementary reflectors from the input/output matrix parameter A_test. */
    copy_matrix(datatype, "full", min_A, n_A, A_test, lda, Q, n_A);
    copy_matrix(datatype, "Lower", m_A, min_A, A_test, lda, L, m_A);

    switch( datatype )
    {
        case FLOAT:
        {
            float twork;
            float norm, norm_A, eps, resid1, resid2;

            /* sorglq api generates the Q martrix using the elementary reflectors and scalar
               factor values */
            fla_lapack_sorglq(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, info);
            if (*info < 0)
               break;

            lwork = twork;
            create_vector(datatype,  &work, lwork);

            fla_lapack_sorglq(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, info);
            if(*info < 0)
               break;

            /* Test 1
               compute norm(L - A*Q') / (V * norm(A) * EPS) */
            sgemm_("N", "T", &m_A, &n_A, &n_A, &s_n_one, A, &lda, Q, &n_A, &s_one, L, &m_A);

            norm_A = fla_lapack_slange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_slange("1", &m_A, &n_A, L, &m_A, work);

            eps = fla_lapack_slamch("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS) */
            resid2 = (float)check_orthogonality(datatype, Q, n_A, n_A, n_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double twork;
            double norm, norm_A, eps, resid1, resid2;

            /* dorglq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_dorglq(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, info);
            if(*info < 0)
               break;

            lwork = twork;
            create_vector(datatype,  &work, lwork);

            fla_lapack_dorglq(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, info);
            if(*info < 0)
               break;

            /* Test 1
               compute norm(L - A*Q') / (V * norm(A) * EPS)*/
            dgemm_("N", "T", &m_A, &n_A, &n_A, &d_n_one, A, &lda, Q, &n_A, &d_one, L, &m_A);

            norm_A = fla_lapack_dlange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_dlange("1", &m_A, &n_A, L, &m_A, work);

            eps = fla_lapack_dlamch("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n_A, n_A, n_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            scomplex twork;
            float norm, norm_A, eps, resid1, resid2;

            /* corglq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_cunglq(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, info);
            if(*info < 0)
               break;

            lwork = twork.real;
            create_vector(datatype,  &work, lwork);

            fla_lapack_cunglq(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, info);
            if(*info < 0)
               break;

            /* Test 1
               compute norm(L - A*Q') / (V * norm(A) * EPS)*/
            cgemm_("N", "C", &m_A, &n_A, &n_A, &c_n_one, A, &lda, Q, &n_A, &c_one, L, &m_A);

            norm_A = fla_lapack_clange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_clange("1", &m_A, &n_A, L, &m_A, work);

            eps = fla_lapack_slamch("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, n_A, n_A, n_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex twork;
            double norm, norm_A, eps, resid1, resid2;

            /* zorglq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            fla_lapack_zunglq(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, info);
            if(*info < 0)
               break;

            lwork = twork.real;
            create_vector(datatype, &work, lwork);

            fla_lapack_zunglq(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, info);
            if(*info < 0)
               break;

            /* Test 1
               compute norm(L - Q'*A) / (V * norm(A) * EPS)*/
            zgemm_("N", "C", &m_A, &n_A, &n_A, &z_n_one, A, &lda, Q, &n_A, &z_one, L, &m_A);

            norm_A = fla_lapack_zlange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_zlange("1", &m_A, &n_A, L, &m_A, work);

            eps = fla_lapack_dlamch("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n_A, n_A, n_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }

    // Free up buffers
    free_matrix( L );
    free_matrix( Q );
    free_vector( work );
}

