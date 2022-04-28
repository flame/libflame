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
	void *L, *U, *T, *work;
	
	m_n_vector = m_A * n_A;
	min_A = min(m_A, n_A);
	create_matrix(datatype, &L, m_A, m_A);
	create_matrix(datatype, &U, m_A, n_A);
	reset_matrix(datatype, m_A, n_A, U, m_A);
	create_matrix(datatype, &T, m_A, n_A);
	create_vector(datatype, &work, 2 * m_A);
	/* Lower triangular matrix should be sqare matrix m x m */
	/* For m==i OR  m < n OR m > n -->  A(mxn) = L(mxm) * U(mxn) */
	copy_matrix(datatype, "Lower", m_A, m_A, A_test, m_A, L, m_A);
	copy_matrix(datatype, "Upper", m_A, n_A, A_test, m_A, U, m_A);


	switch (datatype)
	{
	case FLOAT:
	{
		float norm, norm_A, eps;

		/* Unity diagonal elements to Lower triangular matrix */
		slaset_("U", &m_A, &m_A, &s_zero, &s_one, L, &m_A);
		norm_A = slange_("1", &m_A, &n_A, A, &m_A, work);
		eps = slamch_("Epsilon");
		/* T = L * U  */
		sgemm_("N", "N", &m_A, &n_A, &m_A, &s_one, L, &m_A, U, &m_A, &s_zero, T, &m_A);
		/*  Row interchanges based on IPIV values */
		slaswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
		/* T - A --> L*U - A */
		saxpy_(&m_n_vector, &s_n_one, A, &i_one, T, &i_one);
		/* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
		norm = slange_("1", &m_A, &n_A, T, &m_A, work);

		*residual = norm / (float)n_A / norm_A / eps;
	}

	case DOUBLE:
	{
		double norm, norm_A, eps;

		/* Unity diagonal elements to Lower triangular matrix */
		dlaset_("U",&m_A, &m_A, &d_zero, &d_one, L, &m_A);
		norm_A = dlange_("1", &m_A, &n_A, A, &m_A, work);
		eps = dlamch_("Epsilon");
		/* T = L * U  */
		dgemm_("N", "N", &m_A, &n_A, &m_A, &d_one, L, &m_A, U, &m_A, &d_zero, T, &m_A);
		/*  Row interchanges based on IPIV values*/
		dlaswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
		/* T - A --> L*U - A */
		daxpy_(&m_n_vector, &d_n_one, A, &i_one, T, &i_one);
		/* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
		norm = dlange_("1", &m_A, &n_A, T, &m_A, work);

		*residual = norm / (double)n_A / norm_A / eps;

		break;
	}
	case COMPLEX:
	{
		float norm, norm_A, eps;

		/* Unity diagonal elements to Lower triangular matrix */
		claset_("U", &m_A, &m_A, &c_zero, &c_one, L, &m_A);

		norm_A = clange_("1", &m_A, &n_A, A, &m_A, work);
		eps = slamch_("Epsilon");
		/* T = L * U  */
		cgemm_("N", "N", &m_A, &n_A, &m_A, &c_one, L, &m_A, U, &m_A, &c_zero, T, &m_A);
		/*  Row interchanges based on IPIV values*/
		claswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
		/* T - A --> L*U - A */
		caxpy_(&m_n_vector, &c_n_one, A, &i_one, T, &i_one);
		/* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
		norm = clange_("1", &m_A, &n_A, T, &m_A, work);

		*residual = norm / (float)n_A / norm_A /eps;
	}
	case DOUBLE_COMPLEX:
	{
		double norm, norm_A, eps;

		/* Unity diagonal elements to Lower triangular matrix */
		zlaset_("U", &m_A, &m_A, &z_zero, &z_one, L, &m_A);

		norm_A = zlange_("1", &m_A, &n_A, A, &m_A, work);
		eps = dlamch_("Epsilon");
		/* T = L * U  */
		zgemm_("N", "N", &m_A, &n_A, &m_A, &z_one, L, &m_A, U, &m_A, &z_zero, T, &m_A);
		/*  Row interchanges based on IPIV values*/
		zlaswp_(&n_A, T, &m_A, &i_one, &min_A, IPIV, &i_n_one);
		/* T - A --> L*U - A */
		zaxpy_(&m_n_vector, &z_n_one, A, &i_one, T, &i_one);
		/* Compute norm( L*U - A ) / ( N * norm(A) * EPS ) */
		norm = zlange_("1", &m_A, &n_A, T, &m_A, work);

		*residual = norm / (double)n_A / norm_A / eps;
	}
	}

	// Free up buffers
	free_matrix(L);
	free_matrix(U);
	free_matrix(T);
	free_vector(work);
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
