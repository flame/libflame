/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file validate_geqp3.c
 *  @brief Defines validate function of GEQP3() to use in test suite.
 *  */

#include "test_common.h"

void validate_geqp3(integer m_A, integer n_A,
    void *A,
    void *A_test,
    integer *jpvt,
    void *T_test,
    integer datatype,
    double* residual)
{
    void *Q = NULL, *R = NULL, *work = NULL;
    integer cs_A, min_A;
    integer lwork = -1, tinfo, TRUE = 1;

    cs_A = m_A;
    min_A = min(m_A, n_A);

    // Create Q and R matrices.
    create_matrix(datatype, &Q, m_A, m_A);
    create_matrix(datatype, &R, m_A, n_A);
    reset_matrix(datatype, m_A, m_A, Q, m_A);
    reset_matrix(datatype, m_A, n_A, R, cs_A);

    // Extract R matrix and elementary reflectors from the input/output matrix parameter A_test.
    if(m_A <= n_A)
    {
        copy_matrix(datatype, "full", m_A, m_A, A_test, m_A, Q, m_A);
        copy_matrix(datatype, "Upper", m_A, n_A, A_test, m_A, R, m_A);
    }
    else
    {
        copy_matrix(datatype, "full", m_A, n_A, get_m_ptr(datatype, A_test, 1, 0, m_A), m_A, get_m_ptr(datatype, Q, 1, 0, m_A), m_A);
        copy_matrix(datatype, "Upper", n_A, n_A, A_test, m_A, R, m_A);
    }

    switch( datatype )
    {
        case FLOAT:
        {
            float twork;
            float norm, norm_A, eps, resid1, resid2;
            eps = slamch_("P");

            /* permute A using the permuted vector jpvt to get (A * P) */
            slapmt_(&TRUE, &m_A, &n_A, A, &m_A, jpvt);
            norm_A = slange_("1", &m_A, &n_A, A, &m_A, work);

            /* sorgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            sorgqr_(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &tinfo);
            lwork = (integer)twork;
            create_vector(datatype, &work, lwork);
            sorgqr_(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm(((Q * R) - (A * P)) / (N * norm(A) * EPS)*/
            sgemm_("N", "N", &m_A, &n_A, &m_A, &s_one, Q, &m_A, R, &m_A, &s_n_one, A, &m_A);
            norm = slange_("1", &m_A, &n_A, A, &m_A, work);

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m_A, n_A);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double twork;
            double norm, norm_A, eps, resid1, resid2;

            eps = dlamch_("P");
            /* permute A using the permuted vector jpvt to get (A * P) */
            dlapmt_(&TRUE, &m_A, &n_A, A, &m_A, jpvt);
            norm_A = dlange_("1", &m_A, &n_A, A, &m_A, work);

            /* dorgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            dorgqr_(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &tinfo);
            lwork = (integer)twork;
            create_vector(datatype,  &work, lwork);

            dorgqr_(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm((Q * R) - (A * P)) / (N * norm(A) * EPS)*/
            dgemm_("N", "N", &m_A, &n_A, &m_A, &d_one, Q, &m_A, R, &m_A, &d_n_one, A, &m_A);

            norm = dlange_("1", &m_A, &n_A, A, &m_A, work);
            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m_A, n_A); 

            *residual = (double)max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            scomplex twork;
            float norm, norm_A, eps, resid1, resid2;

            eps = slamch_("P");
            /* permute A using the permuted vector jpvt to get (A * P) */
            clapmt_(&TRUE, &m_A, &n_A, A, &m_A, jpvt);
            norm_A = clange_("1", &m_A, &n_A, A, &m_A, work);

            /* corgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            cungqr_(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &tinfo);

            lwork = (integer)twork.real;
            create_vector(datatype,  &work, lwork);

            cungqr_(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm((Q * R) - (A * P)) / (N * norm(A) * EPS)*/
            cgemm_("N", "N", &m_A, &n_A, &m_A, &c_one, Q, &m_A, R, &m_A, &c_n_one, A, &m_A);

            norm = clange_("1", &m_A, &n_A, A, &m_A, work);
            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = (float)check_orthogonality(datatype, Q, m_A, n_A); 

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex twork;
            double norm, norm_A, eps, resid1, resid2;

            eps = dlamch_("P");
            /* permute A using the permuted vector jpvt to get (A * P) */
            zlapmt_(&TRUE, &m_A, &n_A, A, &m_A, jpvt);
            norm_A = zlange_("1", &m_A, &n_A, A, &m_A, work);

            /* zorgrq api generates the Q martrix using the elementary reflectors and scalar
               factor values*/
            zungqr_(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &tinfo);

            lwork = (integer)twork.real;
            create_vector(datatype, &work, lwork);

            zungqr_(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm((Q * R) - (A * P)) / (N * norm(A) * EPS)*/
            zgemm_("N", "N", &m_A, &n_A, &m_A, &z_n_one, Q, &m_A, R, &m_A, &z_one, A, &m_A);

            norm = zlange_("1", &m_A, &n_A, A, &m_A, work);
            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            resid2 = check_orthogonality(datatype, Q, m_A, n_A);  

            *residual = (double)max(resid1, resid2);
            break;
        }
    }
    // Free up buffers
    free_matrix( R );
    free_matrix( Q );
    free_vector( work );
}
