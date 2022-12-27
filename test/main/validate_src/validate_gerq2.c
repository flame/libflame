/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_gerq2.c
 *  @brief Defines validate function of GERQ2() to use in test suite.
 *  */

#include "test_common.h"

void validate_gerq2(integer m_A,
    integer n_A,
    void *A,
    void *A_test,
    integer lda,
    void *T_test,
    integer datatype,
    double* residual)
{
    void *R = NULL, *Q = NULL, *work = NULL;
    integer min_A, diff_A;
    integer lwork = -1, tinfo;

    min_A = fla_min(m_A, n_A);

    if(m_A <= n_A)
        diff_A = n_A - m_A;
    else
        diff_A = m_A - n_A;

    // Create R and Q matrices.
    create_matrix(datatype, &R, m_A, n_A);
    create_matrix(datatype, &Q, n_A, n_A);
    reset_matrix(datatype, m_A, n_A, R, m_A);
    reset_matrix(datatype, n_A, n_A, Q, n_A);

    // Extract R matrix and elementary reflectors from the input/output matrix parameter A_test.
    if(m_A <= n_A)
    {
        copy_matrix(datatype, "Upper", m_A, m_A, get_m_ptr(datatype, A_test, 0, diff_A, lda), lda, get_m_ptr(datatype, R, 0, diff_A, m_A), m_A);
        copy_matrix(datatype, "full", m_A, n_A, A_test, lda, get_m_ptr(datatype, Q, diff_A, 0, n_A), n_A);
    }
    else
    {
        copy_matrix(datatype, "full", diff_A, n_A, A_test, lda, R, m_A);
        copy_matrix(datatype, "Upper", n_A, n_A, get_m_ptr(datatype, A_test, diff_A, 0, lda), lda, get_m_ptr(datatype, R, diff_A, 0, m_A), m_A);
        copy_matrix(datatype, "full", n_A, n_A, get_m_ptr(datatype, A_test, diff_A, 0, lda), lda, Q, n_A);
    }

    switch( datatype )
    {
        case FLOAT:
        {
            float twork;
            float norm, norm_A, eps, resid1, resid2;

            /* sorgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            sorgrq_(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, &tinfo);

            lwork = twork;
            create_vector(datatype,  &work, lwork);

            sorgrq_(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            sgemm_("N", "T", &m_A, &n_A, &n_A, &s_n_one, A, &lda, Q, &n_A, &s_one, R, &m_A);

            norm_A = slange_("1", &m_A, &n_A, A, &lda, work);
            norm = slange_("1", &m_A, &n_A, R, &m_A, work);

            eps = slamch_("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q'*Q) / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, n_A, n_A, n_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double twork;
            double norm, norm_A, eps, resid1, resid2;

            /* dorgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            dorgrq_(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, &tinfo);

            lwork = twork;
            create_vector(datatype,  &work, lwork);

            dorgrq_(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            dgemm_("N", "T", &m_A, &n_A, &n_A, &d_n_one, A, &lda, Q, &n_A, &d_one, R, &m_A);

            norm_A = dlange_("1", &m_A, &n_A, A, &lda, work);
            norm = dlange_("1", &m_A, &n_A, R, &m_A, work);

            eps = dlamch_("P");

            resid1 = norm/(eps * norm_A * (double)n_A);


            /* Test 2
               compute norm(I - Q'*Q) / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n_A, n_A, n_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            scomplex twork;
            float norm, norm_A, eps, resid1, resid2;

            /* dorgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            cungrq_(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, &tinfo);

            lwork = twork.real;
            create_vector(datatype,  &work, lwork);

            cungrq_(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            cgemm_("N", "C", &m_A, &n_A, &n_A, &c_n_one, A, &lda, Q, &n_A, &c_one, R, &m_A);
            norm_A = clange_("1", &m_A, &n_A, A, &lda, work);
            norm = clange_("1", &m_A, &n_A, R, &m_A, work);

            eps = slamch_("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q'*Q) / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, n_A, n_A, n_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex twork;
            double norm, norm_A, eps, resid1, resid2;

            /* dorgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            zungrq_(&n_A, &n_A, &min_A, NULL, &n_A, NULL, &twork, &lwork, &tinfo);

            lwork = twork.real;
            create_vector(datatype,  &work, lwork);

            zungrq_(&n_A, &n_A, &min_A, Q, &n_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            zgemm_("N", "C", &m_A, &n_A, &n_A, &z_n_one, A, &lda, Q, &n_A, &z_one, R, &m_A);
            norm_A = zlange_("1", &m_A, &n_A, A, &lda, work);
            norm = zlange_("1", &m_A, &n_A, R, &m_A, work);

            eps = dlamch_("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q'*Q) / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, n_A, n_A, n_A);

            *residual = (double)fla_max(resid1, resid2);
            break;
        }
    }

    // Free up buffers
    free_matrix( R );
    free_matrix( Q );
    free_vector( work );
}
