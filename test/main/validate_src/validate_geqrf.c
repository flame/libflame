/******************************************************************************
* Copyright (C) 2022-2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_geqrf.c
 *  @brief Defines validate function of GEQRF() to use in test suite.
 *  */

#include "test_common.h"

void validate_geqrf(integer m_A,
    integer n_A,
    void *A,
    void *A_test,
    integer lda,
    void *T_test,
    integer datatype,
    double* residual,
    integer* info)
{
    void *Q = NULL, *R = NULL, *work = NULL;
    integer min_A;
    integer lwork = -1;
    *info = 0;

    min_A = fla_min(m_A, n_A);

    // Create Q and R matrices.
    create_matrix(datatype, &Q, m_A, m_A);
    create_matrix(datatype, &R, m_A, n_A);
    reset_matrix(datatype, m_A, m_A, Q, m_A);
    reset_matrix(datatype, m_A, n_A, R, m_A);

    // Extract R matrix and elementary reflectors from the input/output matrix parameter A_test.
    copy_matrix(datatype, "full", m_A, min_A, A_test, lda, Q, m_A);
    copy_matrix(datatype, "Upper", min_A, n_A, A_test, lda, R, m_A);

    switch( datatype )
    {
        case FLOAT:
        {
            float twork;
            float norm, norm_A, eps, resid1, resid2;

            /* sorgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            fla_lapack_sorgqr(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, info);
            if (*info < 0)
               break;

            lwork = twork;
            create_vector(datatype,  &work, lwork);

            fla_lapack_sorgqr(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, info);
            if(*info < 0)
               break;

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            sgemm_("T", "N", &m_A, &n_A, &m_A, &s_n_one, Q, &m_A, A, &lda, &s_one, R, &m_A);

            norm_A = fla_lapack_slange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_slange("1", &m_A, &n_A, R, &m_A, work);

            eps = fla_lapack_slamch("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m_A, m_A, m_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double twork;
            double norm, norm_A, eps, resid1, resid2;

            /* dorgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            fla_lapack_dorgqr(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, info);
            if(*info < 0)
               break;

            lwork = twork;
            create_vector(datatype,  &work, lwork);

            fla_lapack_dorgqr(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, info);
            if(*info < 0)
               break;

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            dgemm_("T", "N", &m_A, &n_A, &m_A, &d_n_one, Q, &m_A, A, &lda, &d_one, R, &m_A);

            norm_A = fla_lapack_dlange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_dlange("1", &m_A, &n_A, R, &m_A, work);

            eps = fla_lapack_dlamch("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m_A, m_A, m_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            scomplex twork;
            float norm, norm_A, eps, resid1, resid2;

            /* corgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            fla_lapack_cungqr(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, info);
            if(*info < 0)
               break;

            lwork = twork.real;
            create_vector(datatype,  &work, lwork);

            fla_lapack_cungqr(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, info);
            if(*info < 0)
               break;

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            cgemm_("C", "N", &m_A, &n_A, &m_A, &c_n_one, Q, &m_A, A, &lda, &c_one, R, &m_A);

            norm_A = fla_lapack_clange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_clange("1", &m_A, &n_A, R, &m_A, work);

            eps = fla_lapack_slamch("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m_A, m_A, m_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex twork;
            double norm, norm_A, eps, resid1, resid2;

            /* zorgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            fla_lapack_zungqr(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, info);
            if(*info < 0)
               break;
            
            lwork = twork.real;
            create_vector(datatype, &work, lwork);

            fla_lapack_zungqr(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, info);
            if(*info < 0)
               break;
            
            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            zgemm_("C", "N", &m_A, &n_A, &m_A, &z_n_one, Q, &m_A, A, &lda, &z_one, R, &m_A);

            norm_A = fla_lapack_zlange("1", &m_A, &n_A, A, &lda, work);
            norm = fla_lapack_zlange("1", &m_A, &n_A, R, &m_A, work);

            eps = fla_lapack_dlamch("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m_A, m_A, m_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }

    // Free up buffers
    free_matrix( R );
    free_matrix( Q );
    free_vector( work );
}