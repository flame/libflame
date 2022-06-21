/******************************************************************************
* Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/

/*! @file main.cc
 *  @brief Defines functions to use in linear solver APIs of test suite.
 *  */

#include "test_common.h"

void validate_geqrf(integer m_A,
    integer n_A,
    void *A,
    void *A_test,
    void *T_test,
    integer datatype,
    double* residual)
{
    void *Q = NULL, *R = NULL, *Ibuff = NULL, *work = NULL;
    integer cs_A, min_A;
    integer lwork = -1, tinfo;

    cs_A = m_A;
    min_A = min(m_A, n_A);

    // Create Q and R matrices.
    create_matrix(datatype, &Q, m_A, m_A);
    create_matrix(datatype, &R, m_A, n_A);
    reset_matrix(datatype, m_A, m_A, Q, m_A);
    reset_matrix(datatype, m_A, n_A, R, cs_A);

    // Create Identity matrix to validate orthogonal property of matrix Q
    create_matrix(datatype, &Ibuff, m_A, m_A);

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

            /* sorgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            sorgqr_(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &tinfo);

            lwork = twork;
            create_vector(datatype,  &work, lwork);

            sorgqr_(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            sgemm_("T", "N", &m_A, &n_A, &m_A, &s_n_one, Q, &m_A, A, &m_A, &s_one, R, &m_A);

            norm_A = slange_("1", &m_A, &n_A, A, &m_A, work);
            norm = slange_("1", &m_A, &n_A, R, &m_A, work);

            eps = slamch_("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            slaset_("full", &m_A, &m_A, &s_zero, &s_one, Ibuff, &m_A);
            sgemm_("N", "T", &m_A, &m_A, &m_A, &s_n_one, Q, &m_A, Q, &m_A, &s_one, Ibuff, &m_A);

            norm = slange_("1", &m_A, &m_A, Ibuff, &m_A, work);
            resid2 = norm/(eps * (float)n_A);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double twork;
            double norm, norm_A, eps, resid1, resid2;

            /* dorgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            dorgqr_(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &tinfo);

            lwork = twork;
            create_vector(datatype,  &work, lwork);

            dorgqr_(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            dgemm_("T", "N", &m_A, &n_A, &m_A, &d_n_one, Q, &m_A, A, &m_A, &d_one, R, &m_A);

            norm_A = dlange_("1", &m_A, &n_A, A, &m_A, work);
            norm = dlange_("1", &m_A, &n_A, R, &m_A, work);

            eps = dlamch_("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            dlaset_("full", &m_A, &m_A, &d_zero, &d_one, Ibuff, &m_A);
            dgemm_("N", "T", &m_A, &m_A, &m_A, &d_n_one, Q, &m_A, Q, &m_A, &d_one, Ibuff, &m_A);

            norm = dlange_("1", &m_A, &m_A, Ibuff, &m_A, work);
            resid2 = norm/(eps * (double)n_A);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            scomplex twork;
            float norm, norm_A, eps, resid1, resid2;

            /* corgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            cungqr_(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &tinfo);

            lwork = twork.real;
            create_vector(datatype,  &work, lwork);

            cungqr_(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            cgemm_("C", "N", &m_A, &n_A, &m_A, &c_n_one, Q, &m_A, A, &m_A, &c_one, R, &m_A);

            norm_A = clange_("1", &m_A, &n_A, A, &m_A, work);
            norm = clange_("1", &m_A, &n_A, R, &m_A, work);

            eps = slamch_("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            claset_("full", &m_A, &m_A, &c_zero, &c_one, Ibuff, &m_A);
            cgemm_("N", "C", &m_A, &m_A, &m_A, &c_n_one, Q, &m_A, Q, &m_A, &c_one, Ibuff, &m_A);

            norm = clange_("1", &m_A, &m_A, Ibuff, &m_A, work);
            resid2 = norm/(eps * (float)n_A);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            dcomplex twork;
            double norm, norm_A, eps, resid1, resid2;

            /* zorgrq api generates the Q martrix using the elementary reflectors and scalar 
               factor values*/
            zungqr_(&m_A, &m_A, &min_A, NULL, &m_A, NULL, &twork, &lwork, &tinfo);

            lwork = twork.real;
            create_vector(datatype, &work, lwork);

            zungqr_(&m_A, &m_A, &min_A, Q, &m_A, T_test, work, &lwork, &tinfo);

            /* Test 1
               compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
            zgemm_("C", "N", &m_A, &n_A, &m_A, &z_n_one, Q, &m_A, A, &m_A, &z_one, R, &m_A);

            norm_A = zlange_("1", &m_A, &n_A, A, &m_A, work);
            norm = zlange_("1", &m_A, &n_A, R, &m_A, work);

            eps = dlamch_("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            zlaset_("full", &m_A, &m_A, &z_zero, &z_one, Ibuff, &m_A);
            zgemm_("N", "C", &m_A, &m_A, &m_A, &z_n_one, Q, &m_A, Q, &m_A, &z_one, Ibuff, &m_A);

            norm = zlange_("1", &m_A, &m_A, Ibuff, &m_A, work);
            resid2 = norm/(eps * (double)n_A);

            *residual = (double)max(resid1, resid2);

            break;
        }
    }

    // Free up buffers
    free_matrix( R );
    free_matrix( Q );
    free_matrix( Ibuff );
    free_vector( work );
}

void validate_gerq2(integer m_A,
    integer n_A,
    void *A,
    void *A_test,
    void *T_test,
    integer datatype,
    double* residual)
{
    void *R = NULL, *Q = NULL, *Ibuff = NULL, *work = NULL;
    integer cs_A, min_A, diff_A;
    integer lwork = -1, tinfo;

    cs_A = m_A;
    min_A = min(m_A, n_A);

    if(m_A <= n_A)
        diff_A = n_A - m_A;
    else
        diff_A = m_A - n_A;

    // Create R and Q matrices.
    create_matrix(datatype, &R, m_A, n_A);
    create_matrix(datatype, &Q, n_A, n_A);
    reset_matrix(datatype, m_A, n_A, R, cs_A);
    reset_matrix(datatype, n_A, n_A, Q, n_A);

    // Create Identity matrix to validate orthogonal property of matrix Q
    create_matrix(datatype, &Ibuff, n_A, n_A);

    // Extract R matrix and elementary reflectors from the input/output matrix parameter A_test.
    if(m_A <= n_A)
    {
        copy_matrix(datatype, "Upper", m_A, m_A, get_m_ptr(datatype, A_test, 0, diff_A, m_A), m_A, get_m_ptr(datatype, R, 0, diff_A, m_A), m_A);
        copy_matrix(datatype, "full", m_A, n_A, A_test, m_A, get_m_ptr(datatype, Q, diff_A, 0, n_A), n_A);
    }
    else
    {
        copy_matrix(datatype, "full", diff_A, n_A, A_test, m_A, R, m_A);
        copy_matrix(datatype, "Upper", n_A, n_A, get_m_ptr(datatype, A_test, diff_A, 0, m_A), m_A, get_m_ptr(datatype, R, diff_A, 0, m_A), m_A);
        copy_matrix(datatype, "full", n_A, n_A, get_m_ptr(datatype, A_test, diff_A, 0, m_A), m_A, Q, n_A);
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
            sgemm_("N", "T", &m_A, &n_A, &n_A, &s_n_one, A, &m_A, Q, &n_A, &s_one, R, &m_A);

            norm_A = slange_("1", &m_A, &n_A, A, &m_A, work);
            norm = slange_("1", &m_A, &n_A, R, &m_A, work);

            eps = slamch_("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            slaset_("full", &n_A, &n_A, &s_zero, &s_one, Ibuff, &n_A);
            sgemm_("N", "T", &n_A, &n_A, &n_A, &s_n_one, Q, &n_A, Q, &n_A, &s_one, Ibuff, &n_A);

            norm = slange_("1", &n_A, &n_A, Ibuff, &n_A, work);
            resid2 = norm/(eps * (float)n_A);

            *residual = (double)max(resid1, resid2);
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
            dgemm_("N", "T", &m_A, &n_A, &n_A, &d_n_one, A, &m_A, Q, &n_A, &d_one, R, &m_A);

            norm_A = dlange_("1", &m_A, &n_A, A, &m_A, work);
            norm = dlange_("1", &m_A, &n_A, R, &m_A, work);

            eps = dlamch_("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            dlaset_("full", &n_A, &n_A, &d_zero, &d_one, Ibuff, &n_A);
            dgemm_("N", "T", &n_A, &n_A, &n_A, &d_n_one, Q, &n_A, Q, &n_A, &d_one, Ibuff, &n_A);

            norm = dlange_("1", &n_A, &n_A, Ibuff, &n_A, work);
            resid2 = norm/(eps * (double)n_A);

            *residual = (double)max(resid1, resid2);
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
            cgemm_("N", "C", &m_A, &n_A, &n_A, &c_n_one, A, &m_A, Q, &n_A, &c_one, R, &m_A);
            norm_A = clange_("1", &m_A, &n_A, A, &m_A, work);
            norm = clange_("1", &m_A, &n_A, R, &m_A, work);

            eps = slamch_("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            claset_("full", &n_A, &n_A, &c_zero, &c_one, Ibuff, &n_A);
            cgemm_("N", "C", &n_A, &n_A, &n_A, &c_n_one, Q, &n_A, Q, &n_A, &c_one, Ibuff, &n_A);

            norm = clange_("1", &n_A, &n_A, Ibuff, &n_A, work);
            resid2 = norm/(eps * (float)n_A);

            *residual = (double)max(resid1, resid2);
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
            zgemm_("N", "C", &m_A, &n_A, &n_A, &z_n_one, A, &m_A, Q, &n_A, &z_one, R, &m_A);
            norm_A = zlange_("1", &m_A, &n_A, A, &m_A, work);
            norm = zlange_("1", &m_A, &n_A, R, &m_A, work);

            eps = dlamch_("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            zlaset_("full", &n_A, &n_A, &z_zero, &z_one, Ibuff, &n_A);
            zgemm_("N", "C", &n_A, &n_A, &n_A, &z_n_one, Q, &n_A, Q, &n_A, &z_one, Ibuff, &n_A);

            norm = zlange_("1", &n_A, &n_A, Ibuff, &n_A, work);
            resid2 = norm/(eps * (double)n_A);

            *residual = (double)max(resid1, resid2);
            break;
        }
    }

    // Free up buffers
    free_matrix( R );
    free_matrix( Q );
    free_matrix( Ibuff );
    free_vector( work );
}

void validate_gerqf(integer m_A,
    integer n_A,
    void *A,
    void *A_test,
    void *T_test,
    integer datatype,
    double* residual)
{
    void *R, *Q, *Ibuff, *work;
    integer cs_A, min_A, diff_A;
    integer lwork = -1, tinfo;

    cs_A = m_A;
    min_A = min(m_A, n_A);

    if(m_A <= n_A)
        diff_A = n_A - m_A;
    else
        diff_A = m_A - n_A;

    // Create R and Q matrices.
    create_matrix(datatype, &R, m_A, n_A);
    create_matrix(datatype, &Q, n_A, n_A);
    reset_matrix(datatype, m_A, n_A, R, cs_A);
    reset_matrix(datatype, n_A, n_A, Q, n_A);

    // Create Identity matrix to validate orthogonal property of matrix Q
    create_matrix(datatype, &Ibuff, n_A, n_A);

    // Extract R matrix and elementary reflectors from the input/output matrix parameter A_test.
    if(m_A <= n_A)
    {
        copy_matrix(datatype, "Upper", m_A, m_A, get_m_ptr(datatype, A_test, 0, diff_A, m_A), m_A, get_m_ptr(datatype, R, 0, diff_A, m_A), m_A);
        copy_matrix(datatype, "full", m_A, n_A, A_test, m_A, get_m_ptr(datatype, Q, diff_A, 0, n_A), n_A);
    }
    else
    {
        copy_matrix(datatype, "full", diff_A, n_A, A_test, m_A, R, m_A);
        copy_matrix(datatype, "Upper", n_A, n_A, get_m_ptr(datatype, A_test, diff_A, 0, m_A), m_A, get_m_ptr(datatype, R, diff_A, 0, m_A), m_A);
        copy_matrix(datatype, "full", n_A, n_A, get_m_ptr(datatype, A_test, diff_A, 0, m_A), m_A, Q, n_A);
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
               compute norm(R - Q'*A) / (N * norm(A) * EPS)*/
            sgemm_("N", "T", &m_A, &n_A, &n_A, &s_n_one, A, &m_A, Q, &n_A, &s_one, R, &m_A);

            norm_A = slange_("1", &m_A, &n_A, A, &m_A, work);
            norm = slange_("1", &m_A, &n_A, R, &m_A, work);

            // Get machine precision
            eps = slamch_("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            slaset_("full", &n_A, &n_A, &s_zero, &s_one, Ibuff, &n_A);
            sgemm_("N", "T", &n_A, &n_A, &n_A, &s_n_one, Q, &n_A, Q, &n_A, &s_one, Ibuff, &n_A);

            norm = slange_("1", &n_A, &n_A, Ibuff, &n_A, work);
            resid2 = norm/(eps * (float)n_A);

            *residual = (double)max(resid1, resid2);
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
            dgemm_("N", "T", &m_A, &n_A, &n_A, &d_n_one, A, &m_A, Q, &n_A, &d_one, R, &m_A);

            norm_A = dlange_("1", &m_A, &n_A, A, &m_A, work);
            norm = dlange_("1", &m_A, &n_A, R, &m_A, work);

            eps = dlamch_("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            dlaset_("full", &n_A, &n_A, &d_zero, &d_one, Ibuff, &n_A);
            dgemm_("N", "T", &n_A, &n_A, &n_A, &d_n_one, Q, &n_A, Q, &n_A, &d_one, Ibuff, &n_A);

            norm = dlange_("1", &n_A, &n_A, Ibuff, &n_A, work);
            resid2 = norm/(eps * (double)n_A);

            *residual = (double)max(resid1, resid2);
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
            cgemm_("N", "C", &m_A, &n_A, &n_A, &c_n_one, A, &m_A, Q, &n_A, &c_one, R, &m_A);
            norm_A = clange_("1", &m_A, &n_A, A, &m_A, work);
            norm = clange_("1", &m_A, &n_A, R, &m_A, work);

            eps = slamch_("P");

            resid1 = norm/(eps * norm_A * (float)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            claset_("full", &n_A, &n_A, &c_zero, &c_one, Ibuff, &n_A);
            cgemm_("N", "C", &n_A, &n_A, &n_A, &c_n_one, Q, &n_A, Q, &n_A, &c_one, Ibuff, &n_A);

            norm = clange_("1", &n_A, &n_A, Ibuff, &n_A, work);
            resid2 = norm/(eps * (float)n_A);

            *residual = (double)max(resid1, resid2);
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
            zgemm_("N", "C", &m_A, &n_A, &n_A, &z_n_one, A, &m_A, Q, &n_A, &z_one, R, &m_A);
            norm_A = zlange_("1", &m_A, &n_A, A, &m_A, work);
            norm = zlange_("1", &m_A, &n_A, R, &m_A, work);

            eps = dlamch_("P");

            resid1 = norm/(eps * norm_A * (double)n_A);

            /* Test 2
               compute norm(I - Q*Q') / (N * EPS)*/
            zlaset_("full", &n_A, &n_A, &z_zero, &z_one, Ibuff, &n_A);
            zgemm_("N", "C", &n_A, &n_A, &n_A, &z_n_one, Q, &n_A, Q, &n_A, &z_one, Ibuff, &n_A);

            norm = zlange_("1", &n_A, &n_A, Ibuff, &n_A, work);
            resid2 = norm/(eps * (double)n_A);

            *residual = (double)max(resid1, resid2);
            break;
        }
    }

    // Free up buffers
    free_matrix( R );
    free_matrix( Q );
    free_matrix( Ibuff );
    free_vector( work );
}


void validate_getrf(integer m_A,
    integer n_A,
    void* A,
    void* A_test,/*AFACT*/
    integer* IPIV,
    integer datatype,
    double* residual)
{

    /* System generated locals */
    integer m_n_vector, min_A;
    void *L, *U, *T, *work, *B, *B_test;
    integer nrhs=1, info;

    m_n_vector = m_A * n_A;
    min_A = min(m_A, n_A);
    create_vector(datatype, &B, m_A);
    create_vector(datatype, &B_test, m_A);
    create_matrix(datatype, &L, m_A, m_A);
    create_matrix(datatype, &U, m_A, n_A);
    reset_matrix(datatype, m_A, n_A, U, m_A);
    create_matrix(datatype, &T, m_A, n_A);
    create_vector(datatype, &work, 2 * m_A);

    rand_vector(datatype, B, m_A, 1);
    copy_vector(datatype, m_A, B, 1, B_test, 1);

    /* Lower triangular matrix should be sqare matrix m x m */
    /* For m==i OR  m < n OR m > n -->  A(mxn) = L(mxm) * U(mxn) */
    copy_matrix(datatype, "Lower", m_A, m_A, A_test, m_A, L, m_A);
    copy_matrix(datatype, "Upper", m_A, n_A, A_test, m_A, U, m_A);


    switch (datatype)
    {
    case FLOAT:
    {
        float norm, norm_A, norm_B, eps, resid1, resid2;
        eps = slamch_("Epsilon");
        /* Test 1 */
        if (m_A == n_A)
        {
            norm_B = snrm2_(&m_A, B, &i_one);
            /* Compute X by passing A and B */
            sgetrs_("N", &m_A, &nrhs, A_test, &m_A, IPIV, B_test, &m_A, &info);
            /* Compute AX-B */
            sgemv_("N", &m_A, &m_A, &s_one, A, &m_A, B_test, &i_one, &s_n_one, B, &i_one);
            norm = snrm2_(&m_A, B, &i_one);
            resid1 = norm / (float)m_A / norm_B / eps;
        }
        else
        {
            resid1 = 0.0;
        }
        /* Test 2 */

        /* Unity diagonal elements to Lower triangular matrix */
        slaset_("U", &m_A, &m_A, &s_zero, &s_one, L, &m_A);
        norm_A = slange_("1", &m_A, &n_A, A, &m_A, work);
        /* T = L * U  */
        sgemm_("N", "N", &m_A, &n_A, &m_A, &s_one, L, &m_A, U, &m_A, &s_zero, T, &m_A);
        /*  Row interchanges based on IPIV values */
        slaswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
        /* T - A --> L*U - A */
        saxpy_(&m_n_vector, &s_n_one, A, &i_one, T, &i_one);
        /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
        norm = slange_("1", &m_A, &n_A, T, &m_A, work);

        resid2 = norm / (float)n_A / norm_A / eps;
        *residual = (float)max(resid1, resid2);
    }

    case DOUBLE:
    {
        double norm, norm_A, norm_B, eps, resid1, resid2;

        eps = slamch_("Epsilon");
        /* Test 1 */
        if (m_A == n_A)
        {
            norm_B = dnrm2_(&m_A, B, &i_one);
            /* Compute X by passing A and B */
            dgetrs_("N", &m_A, &nrhs, A_test, &m_A, IPIV, B_test, &m_A, &info);
            /* Compute AX-B */
            dgemv_("N", &m_A, &m_A, &d_one, A, &m_A, B_test, &i_one, &d_n_one, B, &i_one);
            norm = dnrm2_(&m_A, B, &i_one);
            resid1 = norm / (double)m_A / norm_B / eps;
        }
        else
        {
            resid1 = 0.0;
        }
        /* Test 2 */
        /* Unity diagonal elements to Lower triangular matrix */
        dlaset_("U",&m_A, &m_A, &d_zero, &d_one, L, &m_A);
        norm_A = dlange_("1", &m_A, &n_A, A, &m_A, work);
        /* T = L * U  */
        dgemm_("N", "N", &m_A, &n_A, &m_A, &d_one, L, &m_A, U, &m_A, &d_zero, T, &m_A);
        /*  Row interchanges based on IPIV values*/
        dlaswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
        /* T - A --> L*U - A */
        daxpy_(&m_n_vector, &d_n_one, A, &i_one, T, &i_one);
        /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
        norm = dlange_("1", &m_A, &n_A, T, &m_A, work);

        resid2 = norm / (double)n_A / norm_A / eps;
        *residual = (double)max(resid1, resid2);

        break;
    }
    case COMPLEX:
    {
        float norm, norm_A, norm_B, eps, resid1, resid2;

        eps = slamch_("Epsilon");
        /* Test 1 */
        if (m_A == n_A)
        {
            norm_B = snrm2_(&m_A, B, &i_one);
            /* Compute X by passing A and B */
            cgetrs_("N", &m_A, &nrhs, A_test, &m_A, IPIV, B_test, &m_A, &info);
            /* Compute AX-B */
            cgemv_("N", &m_A, &m_A, &c_one, A, &m_A, B_test, &i_one, &c_n_one, B, &i_one);
            norm = snrm2_(&m_A, B, &i_one);
            resid1 = norm / (float)m_A / norm_B / eps;
        }
        else
        {
            resid1 = 0.0;
        }
        /* Test 2 */

        /* Unity diagonal elements to Lower triangular matrix */
        claset_("U", &m_A, &m_A, &c_zero, &c_one, L, &m_A);

        norm_A = clange_("1", &m_A, &n_A, A, &m_A, work);
        /* T = L * U  */
        cgemm_("N", "N", &m_A, &n_A, &m_A, &c_one, L, &m_A, U, &m_A, &c_zero, T, &m_A);
        /*  Row interchanges based on IPIV values*/
        claswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
        /* T - A --> L*U - A */
        caxpy_(&m_n_vector, &c_n_one, A, &i_one, T, &i_one);
        /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
        norm = clange_("1", &m_A, &n_A, T, &m_A, work);

        resid2 = norm / (float)n_A / norm_A / eps;
        *residual = (float)max(resid1, resid2);
    }
    case DOUBLE_COMPLEX:
    {
        double norm, norm_A, norm_B, eps, resid1, resid2;

        eps = slamch_("Epsilon");
        /* Test 1 */
        if (m_A == n_A)
        {
            norm_B = dnrm2_(&m_A, B, &i_one);
            /* Compute X by passing A and B */
            zgetrs_("N", &m_A, &nrhs, A_test, &m_A, IPIV, B_test, &m_A, &info);
            /* Compute AX-B */
            zgemv_("N", &m_A, &m_A, &z_one, A, &m_A, B_test, &i_one, &z_n_one, B, &i_one);
            norm = dnrm2_(&m_A, B, &i_one);
            resid1 = norm / (double)m_A / norm_B / eps;
        }
        else
        {
            resid1 = 0.0;
        }
        /* Test 2 */

        /* Unity diagonal elements to Lower triangular matrix */
        zlaset_("U", &m_A, &m_A, &z_zero, &z_one, L, &m_A);

        norm_A = zlange_("1", &m_A, &n_A, A, &m_A, work);
        /* T = L * U  */
        zgemm_("N", "N", &m_A, &n_A, &m_A, &z_one, L, &m_A, U, &m_A, &z_zero, T, &m_A);
        /*  Row interchanges based on IPIV values*/
        zlaswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
        /* T - A --> L*U - A */
        zaxpy_(&m_n_vector, &z_n_one, A, &i_one, T, &i_one);
        /* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
        norm = zlange_("1", &m_A, &n_A, T, &m_A, work);

        resid2 = norm / (double)n_A / norm_A / eps;
        *residual = (double)max(resid1, resid2);
    }
    }

    // Free up buffers
    free_matrix(L);
    free_matrix(U);
    free_matrix(T);
    free_vector(work);
    free_vector(B);
    free_vector(B_test);

}

void validate_getri(integer m_A,
    integer n_A,
    void* A,
    void* A_inv,
    integer* IPIV,
    integer datatype,
    double* residual)
{

    /* System generated locals */
    void *I, * work;
    
    /* Create Identity matrix */
    create_matrix(datatype, &I, n_A, n_A);
    create_vector(datatype, &work, 2 * m_A);

    switch (datatype)
    {
    case FLOAT:
    {
        float res, norm, norm_A, norm_A_I, eps;

        eps = slamch_("Epsilon");
        norm_A = slange_("1", &m_A, &m_A, A, &m_A, work);
        norm_A_I = slange_("1", &m_A, &m_A, A_inv, &m_A, work);
        res = (1 / norm_A) / norm_A_I;
        /* compute I - A' * A */
        slaset_("full", &m_A, &m_A, &s_zero, &s_one, I, &m_A);
        /* compute I - A' * A */
        sgemm_("N", "N", &m_A, &n_A, &m_A, &s_n_one, A_inv, &m_A, A, &m_A, &s_one, I, &m_A);

        norm = slange_("1", &m_A, &m_A, I, &m_A, work);
        /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
        *residual = (norm * res) / eps / (double)n_A;
    }

    case DOUBLE:
    {
        double res, norm, norm_A, norm_A_I, eps;

        eps = dlamch_("Epsilon");
        norm_A = dlange_("1", &m_A, &m_A, A, &m_A, work);
        norm_A_I = dlange_("1", &m_A, &m_A, A_inv, &m_A, work);
        res = (1 / norm_A) / norm_A_I;
        /* compute I - A' * A */
        dlaset_("full", &m_A, &m_A, &d_zero, &d_one, I, &m_A);
        /* compute I - A' * A */
        dgemm_("N", "N", &m_A, &n_A, &m_A, &d_n_one, A_inv, &m_A, A, &m_A, &d_one, I, &m_A);

        norm = dlange_("1", &m_A, &m_A, I, &m_A, work);
        /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
        *residual = (norm * res) / eps / (double)n_A;

        break;
    }
    case COMPLEX:
    {
        float res, norm, norm_A, norm_A_I, eps;

        eps = slamch_("Epsilon");
        norm_A = clange_("1", &m_A, &m_A, A, &m_A, work);
        norm_A_I = clange_("1", &m_A, &m_A, A_inv, &m_A, work);
        res = (1 / norm_A) / norm_A_I;
        /* compute I - A' * A */
        claset_("full", &m_A, &m_A, &c_zero, &c_one, I, &m_A);
        /* compute I - A' * A */
        cgemm_("N", "N", &m_A, &n_A, &m_A, &c_n_one, A_inv, &m_A, A, &m_A, &c_one, I, &m_A);

        norm = clange_("1", &m_A, &m_A, I, &m_A, work);
        /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
        *residual = (norm * res) / eps / (double)n_A;
    }
    case DOUBLE_COMPLEX:
    {
        double res, norm, norm_A, norm_A_I, eps;

        eps = dlamch_("Epsilon");
        norm_A = zlange_("1", &m_A, &m_A, A, &m_A, work);
        norm_A_I = zlange_("1", &m_A, &m_A, A_inv, &m_A, work);
        res = (1 / norm_A) / norm_A_I;
        /* compute I - A' * A */
        zlaset_("full", &m_A, &m_A, &z_zero, &z_one, I, &m_A);
        /* compute I - A' * A */
        zgemm_("N", "N", &m_A, &n_A, &m_A, &z_n_one, A_inv, &m_A, A, &m_A, &z_one, I, &m_A);

        norm = zlange_("1", &m_A, &m_A, I, &m_A, work);
        /* Compute norm(I - A'*A) / (N * norm(A) * norm(AINV) * EPS)*/
        *residual = (norm * res) / eps / (double)n_A;
    }
    }

    // Free up buffers
    free_vector(work);
    free_vector(I);
}

void validate_potrs(char *uplo, integer m,
    void *A,
    void *A_test,
    integer datatype,
    void *x,
    void *b,
    double* residual)
{
    integer incx = 1, incy = 1;

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = snrm2_(&m, b, &incx);
            eps = slamch_("P");

            /* Compute Ax-b */
            sgemv_("N", &m, &m, &s_one, A, &m, x, &incx, &s_n_one, b, &incy);
            norm = snrm2_(&m, b, &incx);

            resid = norm/(eps * norm_b * (float)m);

            *residual = (double)resid;
            break;
        }
        case DOUBLE:
        {
            double norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = dnrm2_(&m, b, &incx);
            eps = dlamch_("P");

            /* Compute Ax-b */
            dgemv_("N", &m, &m, &d_one, A, &m, x, &incx, &d_n_one, b, &incy);
            norm = dnrm2_(&m, b, &incx);

            resid = norm/(eps * norm_b * (double)m);

            *residual = (double)resid;
            break;
        }
        case COMPLEX:
        {
            float norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = scnrm2_(&m, b, &incx);
            eps = slamch_("P");

            /* Compute Ax-b */
            cgemv_("N", &m, &m, &c_one, A, &m, x, &incx, &c_n_one, b, &incy);
            norm = scnrm2_(&m, b, &incx);

            resid = norm/(eps * norm_b * (float)m);

            *residual = (double)resid;
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_b, norm, eps, resid;

            /* Test 1 */
            norm_b = dznrm2_(&m, b, &incx);
            eps = dlamch_("P");

            /* Compute Ax-b */
            zgemv_("N", &m, &m, &z_one, A, &m, x, &incx, &z_n_one, b, &incy);
            norm = dznrm2_(&m, b, &incx);

            resid = norm/(eps * norm_b * (double)m);

            *residual = (double)resid;
            break;
        }
    }
}

void validate_orgqr(integer m,
    integer n,
    void *A,
    void* Q,
    void *R,
    void* work,
    integer datatype,
    double* residual)
{
    void *I = NULL;

    /* Create Identity matrix to validate orthogonal property of matrix Q */
    create_matrix(datatype, &I, m, m);

    switch( datatype )
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2;

            eps = slamch_("P");

            /* Test 1
               compute norm(I - Q*Q') / (N * EPS)*/
            slaset_("full", &m, &m, &s_zero, &s_one, I, &m);
            sgemm_("N", "T", &m, &m, &m, &s_n_one, Q, &m, Q, &m, &s_one, I, &m);

            norm = slange_("1", &m, &m, I, &m, work);
            resid1 = norm/(eps * (float)n);

                        /* Test 2
                           compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
                        sgemm_("T", "N", &m, &n, &m, &s_n_one, Q, &m, A, &m, &s_one, R, &m);

                        norm_A = slange_("1", &m, &n, A, &m, work);
                        norm = slange_("1", &m, &n, R, &m, work);

                        resid2 = norm/(eps * norm_A * (float)n);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2;

            eps = dlamch_("P");

            /* Test 1
               compute norm(I - Q*Q') / (N * EPS)*/
            dlaset_("full", &m, &m, &d_zero, &d_one, I, &m);
            dgemm_("N", "T", &m, &m, &m, &d_n_one, Q, &m, Q, &m, &d_one, I, &m);

            norm = dlange_("1", &m, &m, I, &m, work);
            resid1 = norm/(eps * (double)n);

            /* Test 2
                           compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
                        dgemm_("T", "N", &m, &n, &m, &d_n_one, Q, &m, A, &m, &d_one, R, &m);

                        norm_A = dlange_("1", &m, &n, A, &m, work);
                        norm = dlange_("1", &m, &n, R, &m, work);

                        resid2 = norm/(eps * norm_A * (double)n);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2;

            eps = slamch_("P");

            /* Test 1
               compute norm(I - Q*Q') / (N * EPS)*/
            claset_("full", &m, &m, &c_zero, &c_one, I, &m);
            cgemm_("N", "C", &m, &m, &m, &c_n_one, Q, &m, Q, &m, &c_one, I, &m);

            norm = clange_("1", &m, &m, I, &m, work);
            resid1 = norm/(eps * (float)n);

             /* Test 2
                           compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
                        cgemm_("C", "N", &m, &n, &m, &c_n_one, Q, &m, A, &m, &c_one, R, &m);

                        norm_A = clange_("1", &m, &n, A, &m, work);
                        norm = clange_("1", &m, &n, R, &m, work);

                        resid2 = norm/(eps * norm_A * (float)n);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2;

            eps = dlamch_("P");

            /* Test 1
               compute norm(I - Q*Q') / (N * EPS)*/
            zlaset_("full", &m, &m, &z_zero, &z_one, I, &m);
            zgemm_("N", "C", &m, &m, &m, &z_n_one, Q, &m, Q, &m, &z_one, I, &m);

            norm = zlange_("1", &m, &m, I, &m, work);
            resid1 = norm/(eps * (double)n);

                        /* Test 2
                           compute norm(R - Q'*A) / (V * norm(A) * EPS)*/
                        zgemm_("C", "N", &m, &n, &m, &z_n_one, Q, &m, A, &m, &z_one, R, &m);

                        norm_A = zlange_("1", &m, &n, A, &m, work);
                        norm = zlange_("1", &m, &n, R, &m, work);

                        resid2 = norm/(eps * norm_A * (double)n);

            *residual = (double)max(resid1, resid2);

            break;
        }
    }

    /* Free up buffers */
    free_matrix( I );
}

void validate_potrf(char *uplo, integer m, void *A, void *A_test, integer datatype, double* residual)
{
    void *b = NULL, *x = NULL;
    void *x_test = NULL, *b_test = NULL;
    void *work = NULL;
    void *A_save = NULL;
    void *buff_A = NULL, *buff_B = NULL;
    integer info;
    integer nrhs = 1, incx = 1, incy = 1;
    char trans_A, trans_B;

    /* Create Matrix buff_A and buff_B for computing LL'*/
    create_matrix(datatype, &buff_A, m, m);
    create_matrix(datatype, &buff_B, m, m);

    /* Reset the matrix to avoid junk entries */
    reset_matrix(datatype, m, m, buff_A, m);
    reset_matrix(datatype, m, m, buff_B, m);

    create_matrix(datatype, &A_save, m, m);

    /* Create vector to compute Ax-b */
    create_vector(datatype, &b, m);
    create_vector(datatype, &x, m);
    create_vector(datatype, &b_test, m);
    create_vector(datatype, &x_test, m);

    /* Generate random vector b */
    rand_vector(datatype, b, m, 1);

    /* Copy lower or upper triangular matrix based on uplo to buff_A, buff_B */
    copy_matrix(datatype, uplo, m, m, A_test, m, buff_A, m);
    copy_matrix(datatype, uplo, m, m, A_test, m, buff_B, m);
    copy_matrix(datatype, "full", m, m, A, m, A_save, m);

    copy_vector(datatype, m, b, 1, b_test, 1);

    /* set transpose flag based on uplo */
    set_transpose(datatype, uplo, &trans_A, &trans_B);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm_b, norm, eps, resid1, resid2;
            float norm_A;

            /* Test 1 */
            norm_A = slange_("1", &m, &m, A, &m, work);

            /* Compute LL'-A */
            sgemm_(&trans_A, &trans_B, &m, &m, &m, &s_one, buff_A, &m, buff_B, &m, &s_n_one, A, &m);

            norm = slange_("1", &m, &m, A, &m, work);
            eps = slamch_("P");

            resid1 = norm/(eps * norm_A * (float)m);

            /* Test 2 */
            copy_matrix(datatype, "full", m, m, A_save, m, A, m);
            norm_b = snrm2_(&m, b, &incx);

            /* Find x to compute Ax-b */
            spotrs_(uplo, &m, &nrhs, A_test, &m, b_test, &m, &info);
            copy_vector(datatype, m, b_test, 1, x,  1);

            /* Compute Ax-b */
            sgemv_("N", &m, &m, &s_one, A, &m, x, &incx, &s_n_one, b, &incy);
            norm = snrm2_(&m, b, &incx);

            resid2 = norm/(eps * norm_b * (float)m);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE:
        {
            double norm_b, norm, eps, resid1, resid2;
            double norm_A;

            /* Test 1 */
            norm_A = dlange_("1", &m, &m, A, &m, work);

            /* compute L*L'-A */
            dgemm_(&trans_A, &trans_B, &m, &m, &m, &d_one, buff_A, &m, buff_B, &m, &d_n_one, A, &m);

            norm = dlange_("1", &m, &m, A, &m, work);
            eps = dlamch_("P");

            resid1 = norm/(eps * norm_A * (double)m);

            /* Test 2 */
            copy_matrix(datatype, "full", m, m, A_save, m, A, m);
            norm_b = dnrm2_(&m, b, &incx);

            /*Compute Ax-b Linear equations and find x */
            dpotrs_(uplo, &m, &nrhs, A_test, &m, b_test, &m, &info);
            copy_vector(datatype, m, b_test, 1, x,  1);

            /* Compute Ax-b */
            dgemv_("N", &m, &m, &d_one, A, &m, x, &incx, &d_n_one, b, &incy);
            norm = dnrm2_(&m, b, &incx);

            resid2 = norm/(eps * norm_b * (double)m);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case COMPLEX:
        {
            float norm_b, norm, eps, resid1, resid2;
            float norm_A;

            /* Test 1 */
            norm_A = clange_("1", &m, &m, A, &m, work);

            /* compute L*L'-A */
            cgemm_(&trans_A, &trans_B, &m, &m, &m, &c_one, buff_A, &m, buff_B, &m, &c_n_one, A, &m);

            norm = clange_("1", &m, &m, A, &m, work);
            eps = slamch_("P");

            resid1 = norm/(eps * norm_A * (float)m);

            /* Test 2 */
            copy_matrix(datatype, "full", m, m, A_save, m, A, m);
            norm_b = scnrm2_(&m, b, &incx);

            /*Find x to compute Ax-b */
            cpotrs_(uplo, &m, &nrhs, A_test, &m, b_test, &m, &info);
            copy_vector(datatype, m, b_test, 1, x,  1);

            /* Compute Ax-b */
            cgemv_("N", &m, &m, &c_one, A, &m, x, &incx, &c_n_one, b, &incy);
            norm = scnrm2_(&m, b, &incx);

            resid2 = norm/(eps * norm_b * (float)m);

            *residual = (double)max(resid1, resid2);
            break;
        }
        case DOUBLE_COMPLEX:
        {
            double norm_b, norm, eps, resid1, resid2;
            double norm_A;

            /* Test 1 */
            norm_A = zlange_("1", &m, &m, A, &m, work);

            /* compute L*L'-A */
            zgemm_(&trans_A, &trans_B, &m, &m, &m, &z_one, buff_A, &m, buff_B, &m, &z_n_one, A, &m);

            norm = zlange_("1", &m, &m, A, &m, work);
            eps = dlamch_("P");

            resid1 = norm/(eps * norm_A * (double)m);

            /* Test 2 */
            copy_matrix(datatype, "full", m, m, A_save, m, A, m);
            norm_b = dznrm2_(&m, b, &incx);

            /* Find x to compute Ax-b */
            zpotrs_(uplo, &m, &nrhs, A_test, &m, b_test, &m, &info);
            copy_vector(datatype, m, b_test, 1, x,  1);

            /* Compute Ax-b */
            zgemv_("N", &m, &m, &z_one, A, &m, x, &incx, &z_n_one, b, &incy);

            norm = dznrm2_(&m, b, &incx);
            eps = dlamch_("P");

            resid2 = norm/(eps * norm_b * (double)m);

            *residual = (double)max(resid1, resid2);
            break;
        }
    }
    /* Free up buffers */
    free_matrix( b );
    free_matrix( x );
    free_matrix( x_test );
    free_matrix( b_test );
    free_matrix( A_save );
    free_matrix( buff_A );
    free_matrix( buff_B );
}


void validate_getrs(char *trans,
    integer m,
    integer n,
    void* A,
    void* B,
    void* X,
    integer datatype,
    double* residual)
{

    switch (datatype)
    {
    case FLOAT:
    {
        float norm_b, norm, eps, resid;

        /* Test 1 */
        norm_b = snrm2_(&m, B, &i_one);
        eps = slamch_("P");

        /* Compute Ax-b */
        sgemv_(trans, &m, &m, &s_one, A, &m, X, &i_one, &s_n_one, B, &i_one);
        norm = snrm2_(&m, B, &i_one);

        resid = norm / (eps * norm_b * (float)m);

        *residual = (double)resid;
        break;
    }
    case DOUBLE:
    {
        double norm_b, norm, eps, resid;

        /* Test 1 */
        norm_b = dnrm2_(&m, B, &i_one);
        eps = dlamch_("P");

        /* Compute Ax-b */
        dgemv_(trans, &m, &m, &d_one, A, &m, X, &i_one, &d_n_one, B, &i_one);
        norm = dnrm2_(&m, B, &i_one);

        resid = norm / (eps * norm_b * (double)m);

        *residual = (double)resid;
        break;
    }
    case COMPLEX:
    {
        float norm_b, norm, eps, resid;

        /* Test 1 */
        norm_b = scnrm2_(&m, B, &i_one);
        eps = slamch_("P");

        /* Compute Ax-b */
        cgemv_(trans, &m, &m, &c_one, A, &m, X, &i_one, &c_n_one, B, &i_one);
        norm = scnrm2_(&m, B, &i_one);

        resid = norm / (eps * norm_b * (float)m);

        *residual = (double)resid;
        break;
    }
    case DOUBLE_COMPLEX:
    {
        double norm_b, norm, eps, resid;

        /* Test 1 */
        norm_b = dznrm2_(&m, B, &i_one);
        eps = dlamch_("P");

        /* Compute Ax-b */
        zgemv_(trans, &m, &m, &z_one, A, &m, X, &i_one, &z_n_one, B, &i_one);
        norm = dznrm2_(&m, B, &i_one);

        resid = norm / (eps * norm_b * (double)m);

        *residual = (double)resid;
        break;
    }
    }
}

void validate_gesdd(char *jobz, integer m, integer n, void* A, void* A_test, void* s, void* U, void* V, integer datatype, double *residual)
{
    void *sigma = NULL, *Usigma = NULL;
    void *Ibuff_U = NULL, *Ibuff_V = NULL;
    void *work = NULL;

    create_matrix(datatype, &sigma, m, n);
    create_matrix(datatype, &Usigma, m, n);
    reset_matrix(datatype, m, n, Usigma, m);

    create_matrix(datatype, &Ibuff_U, m, m);
    create_matrix(datatype, &Ibuff_V, n, n);
    reset_matrix(datatype, m, m, Ibuff_U, m);
    reset_matrix(datatype, n, n, Ibuff_V, n);

    diagonalize_vector(datatype, s, sigma, m, n, m);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2, resid3;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = slange_("1", &m, &n, A, &m, work);
            sgemm_("N", "N", &m, &n, &m, &s_one, U, &m, sigma, &m, &s_zero, Usigma, &m);
            sgemm_("N", "N", &m, &n, &n, &s_one, Usigma, &m, V, &n, &s_n_one, A, &m);
            norm = slange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            slaset_("full", &m, &m, &s_zero, &s_one, Ibuff_U, &m);
            sgemm_("N", "T", &m, &m, &m, &s_one, U, &m, U, &m, &s_n_one, Ibuff_U, &m);
            norm = slange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (float)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            slaset_("full", &n, &n, &s_zero, &s_one, Ibuff_V, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, V, &n, V, &n, &s_n_one, Ibuff_V, &n);
            norm = slange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (float)n);
            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2, resid3;
            eps = dlamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = dlange_("1", &m, &n, A, &m, work);
            dgemm_("N", "N", &m, &n, &m, &d_one, U, &m, sigma, &m, &d_zero, Usigma, &m);
            dgemm_("N", "N", &m, &n, &n, &d_one, Usigma, &m, V, &n, &d_n_one, A, &m);
            norm = dlange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (double)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            dlaset_("full", &m, &m, &d_zero, &d_one, Ibuff_U, &m);
            dgemm_("N", "T", &m, &m, &m, &d_one, U, &m, U, &m, &d_n_one, Ibuff_U, &m);
            norm = dlange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (double)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            dlaset_("full", &n, &n, &d_zero, &d_one, Ibuff_V, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, V, &n, V, &n, &d_n_one, Ibuff_V, &n);
            norm = dlange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (double)n);
            *residual = (double)max(resid1, max(resid2, resid3));
             break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2, resid3;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = clange_("1", &m, &n, A, &m, work);
            cgemm_("N", "N", &m, &n, &m, &c_one, U, &m, sigma, &m, &c_zero, Usigma, &m);
            cgemm_("N", "N", &m, &n, &n, &c_one, Usigma, &m, V, &n, &c_n_one, A, &m);
            norm = clange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            claset_("full", &m, &m, &c_zero, &c_one, Ibuff_U, &m);
            cgemm_("N", "C", &m, &m, &m, &c_one, U, &m, U, &m, &c_n_one, Ibuff_U, &m);
            norm = clange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (float)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            claset_("full", &n, &n, &c_zero, &c_one, Ibuff_V, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, V, &n, V, &n, &c_n_one, Ibuff_V, &n);
            norm = clange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (float)n);
            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2, resid3;
            eps = dlamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = zlange_("1", &m, &n, A, &m, work);
            zgemm_("N", "N", &m, &n, &m, &z_one, U, &m, sigma, &m, &z_zero, Usigma, &m);
            zgemm_("N", "N", &m, &n, &n, &z_one, Usigma, &m, V, &n, &z_n_one, A, &m);
            norm = zlange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (double)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            zlaset_("full", &m, &m, &z_zero, &z_one, Ibuff_U, &m);
            zgemm_("N", "C", &m, &m, &m, &z_one, U, &m, U, &m, &z_n_one, Ibuff_U, &m);
            norm = zlange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (double)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            zlaset_("full", &n, &n, &z_zero, &z_one, Ibuff_V, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, V, &n, V, &n, &z_n_one, Ibuff_V, &n);
            norm = zlange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (float)n);
            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }
    }
    free_matrix(sigma);
    free_matrix(Usigma);
    free_matrix(Ibuff_U);
    free_matrix(Ibuff_V);
    return;
}

void validate_syevd(char* jobz, char* uplo, integer n, void* A, void* A_test, void* w, integer datatype, double* residual)
{
    void *lambda = NULL, *zlambda = NULL;
    void *I = NULL, *Z = NULL;
    void *work = NULL;

    create_matrix(datatype, &lambda, n, n);
    create_matrix(datatype, &zlambda, n, n);
    create_matrix(datatype, &Z, n, n);

    reset_matrix(datatype, n, n, zlambda, n);
    reset_matrix(datatype, n, n, Z, n);

    create_matrix(datatype, &I, n, n);
    reset_matrix(datatype, n, n, I, n);

    copy_matrix(datatype, "full", n, n, A_test, n, Z, n);

    diagonalize_vector(datatype, w, lambda, n, n, n);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = slange_("1", &n, &n, A, &n, work);
            sgemm_("N", "N", &n, &n, &n, &s_one, Z, &n, lambda, &n, &s_zero, zlambda, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, zlambda, &n, Z, &n, &s_n_one, A, &n);
            norm = slange_("1", &n, &n, A, &n, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Z*Z') / (N * EPS)*/
            slaset_("full", &n, &n, &s_zero, &s_one, I, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, Z, &n, Z, &n, &s_n_one, I, &n);
            norm = slange_("1", &n, &n, I, &n, work);
            resid2 = norm/(eps * (float)n);
            *residual = (double)max(resid1, resid2);
            break;
        }

        case DOUBLE:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = dlamch_("P");

            /* Test 1
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = dlange_("1", &n, &n, A, &n, work);
            dgemm_("N", "N", &n, &n, &n, &d_one, Z, &n, lambda, &n, &d_zero, zlambda, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, zlambda, &n, Z, &n, &d_n_one, A, &n);
            norm = dlange_("1", &n, &n, A, &n, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Z*Z') / (N * EPS)*/
            dlaset_("full", &n, &n, &d_zero, &d_one, I, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, Z, &n, Z, &n, &d_n_one, I, &n);
            norm = dlange_("1", &n, &n, I, &n, work);
            resid2 = norm/(eps * (float)n);
            *residual = (double)max(resid1, resid2);
            break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = clange_("1", &n, &n, A, &n, work);
            cgemm_("N", "N", &n, &n, &n, &c_one, Z, &n, lambda, &n, &c_zero, zlambda, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, zlambda, &n, Z, &n, &c_n_one, A, &n);
            norm = clange_("1", &n, &n, A, &n, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Z*Z') / (N * EPS)*/
            claset_("full", &n, &n, &c_zero, &c_one, I, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, Z, &n, Z, &n, &c_n_one, I, &n);
            norm = clange_("1", &n, &n, I, &n, work);
            resid2 = norm/(eps * (float)n);
            *residual = (double)max(resid1, resid2);
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2;
            eps = dlamch_("P");

          /* Test 1
               compute norm(A - (Z * lambda * Z')) / (V * norm(A) * EPS)*/
            norm_A = zlange_("1", &n, &n, A, &n, work);
            zgemm_("N", "N", &n, &n, &n, &z_one, Z, &n, lambda, &n, &z_zero, zlambda, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, zlambda, &n, Z, &n, &z_n_one, A, &n);
            norm = zlange_("1", &n, &n, A, &n, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - Z*Z') / (N * EPS)*/
            zlaset_("full", &n, &n, &z_zero, &z_one, I, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, Z, &n, Z, &n, &z_n_one, I, &n);
            norm = zlange_("1", &n, &n, I, &n, work);
            resid2 = norm/(eps * (float)n);
            *residual = (double)max(resid1, resid2);
            break;
        }
    }
    free_matrix(lambda);
    free_matrix(zlambda);
    free_matrix(I);
    free_matrix(Z);
    return;
}

void validate_gesvd(char *jobu, char *jobvt, integer m, integer n, void* A, void* A_test, void* s, void* U, void* V, integer datatype, double *residual)
{
    void *sigma = NULL, *Usigma = NULL;
    void *Ibuff_U = NULL, *Ibuff_V = NULL;
    void *work = NULL;

    create_matrix(datatype, &sigma, m, n);
    create_matrix(datatype, &Usigma, m, n);
    reset_matrix(datatype, m, n, Usigma, m);

    create_matrix(datatype, &Ibuff_U, m, m);
    create_matrix(datatype, &Ibuff_V, n, n);
    reset_matrix(datatype, m, m, Ibuff_U, m);
    reset_matrix(datatype, n, n, Ibuff_V, n);

    diagonalize_vector(datatype, s, sigma, m, n, m);

    switch(datatype)
    {
        case FLOAT:
        {
            float norm, norm_A, eps, resid1, resid2, resid3;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = slange_("1", &m, &n, A, &m, work);
            sgemm_("N", "N", &m, &n, &m, &s_one, U, &m, sigma, &m, &s_zero, Usigma, &m);
            sgemm_("N", "N", &m, &n, &n, &s_one, Usigma, &m, V, &n, &s_n_one, A, &m);
            norm = slange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            slaset_("full", &m, &m, &s_zero, &s_one, Ibuff_U, &m);
            sgemm_("N", "T", &m, &m, &m, &s_one, U, &m, U, &m, &s_n_one, Ibuff_U, &m);
            norm = slange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (float)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            slaset_("full", &n, &n, &s_zero, &s_one, Ibuff_V, &n);
            sgemm_("N", "T", &n, &n, &n, &s_one, V, &n, V, &n, &s_n_one, Ibuff_V, &n);
            norm = slange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (float)n);
            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }

        case DOUBLE:
       {
            double norm, norm_A, eps, resid1, resid2, resid3;
            eps = dlamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = dlange_("1", &m, &n, A, &m, work);
            dgemm_("N", "N", &m, &n, &m, &d_one, U, &m, sigma, &m, &d_zero, Usigma, &m);
            dgemm_("N", "N", &m, &n, &n, &d_one, Usigma, &m, V, &n, &d_n_one, A, &m);
            norm = dlange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (double)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            dlaset_("full", &m, &m, &d_zero, &d_one, Ibuff_U, &m);
            dgemm_("N", "T", &m, &m, &m, &d_one, U, &m, U, &m, &d_n_one, Ibuff_U, &m);
            norm = dlange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (double)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            dlaset_("full", &n, &n, &d_zero, &d_one, Ibuff_V, &n);
            dgemm_("N", "T", &n, &n, &n, &d_one, V, &n, V, &n, &d_n_one, Ibuff_V, &n);
            norm = dlange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (double)n);
            *residual = (double)max(resid1, max(resid2, resid3));
             break;
        }

        case COMPLEX:
        {
            float norm, norm_A, eps, resid1, resid2, resid3;
            eps = slamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = clange_("1", &m, &n, A, &m, work);
            cgemm_("N", "N", &m, &n, &m, &c_one, U, &m, sigma, &m, &c_zero, Usigma, &m);
            cgemm_("N", "N", &m, &n, &n, &c_one, Usigma, &m, V, &n, &c_n_one, A, &m);
            norm = clange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (float)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            claset_("full", &m, &m, &c_zero, &c_one, Ibuff_U, &m);
            cgemm_("N", "C", &m, &m, &m, &c_one, U, &m, U, &m, &c_n_one, Ibuff_U, &m);
            norm = clange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (float)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            claset_("full", &n, &n, &c_zero, &c_one, Ibuff_V, &n);
            cgemm_("N", "C", &n, &n, &n, &c_one, V, &n, V, &n, &c_n_one, Ibuff_V, &n);
            norm = clange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (float)n);
            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }

        case DOUBLE_COMPLEX:
        {
            double norm, norm_A, eps, resid1, resid2, resid3;
            eps = dlamch_("P");

            /* Test 1
               compute norm(A - (U*sigma*Vt)) / (V * norm(A) * EPS)*/
            norm_A = zlange_("1", &m, &n, A, &m, work);
            zgemm_("N", "N", &m, &n, &m, &z_one, U, &m, sigma, &m, &z_zero, Usigma, &m);
            zgemm_("N", "N", &m, &n, &n, &z_one, Usigma, &m, V, &n, &z_n_one, A, &m);
            norm = zlange_("1", &m, &n, A, &m, work);
            resid1 = norm/(eps * norm_A * (double)n);

            /* Test 2
               compute norm(I - U*U') / (N * EPS)*/
            zlaset_("full", &m, &m, &z_zero, &z_one, Ibuff_U, &m);
            zgemm_("N", "C", &m, &m, &m, &z_one, U, &m, U, &m, &z_n_one, Ibuff_U, &m);
            norm = zlange_("1", &m, &m, Ibuff_U, &m, work);
            resid2 = norm/(eps * (double)m);

            /* Test 3
               compute norm(I - V'*V) / (N * EPS)*/
            zlaset_("full", &n, &n, &z_zero, &z_one, Ibuff_V, &n);
            zgemm_("N", "C", &n, &n, &n, &z_one, V, &n, V, &n, &z_n_one, Ibuff_V, &n);
            norm = zlange_("1", &n, &n, Ibuff_V, &n, work);
            resid3 = norm/(eps * (float)n);
            *residual = (double)max(resid1, max(resid2, resid3));
            break;
        }
    }
    free_matrix(sigma);
    free_matrix(Usigma);
    free_matrix(Ibuff_U);
    free_matrix(Ibuff_V);
    return;
}
