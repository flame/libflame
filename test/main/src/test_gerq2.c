/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_libflame.h"
#include "test_common.h"
#include "test_prototype.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1


// Local prototypes.
void fla_test_gerq2_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer  pci,
									integer  n_repeats, double* perf, double* t, double* residual);
void prepare_gerq2_run(integer m_A, integer n_A, void *A, void *T, integer datatype, integer n_repeats, double* time_min_);
inline void invoke_gerq2(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info);
void validate_gerq2(integer m_A, integer n_A, void *A, void *A_test, void *T_test, integer datatype, double* residual);


void fla_test_gerq2(test_params_t *params)
{
	char* op_str = "RQ factorization with unblocked algorithm";
	char* front_str = "GERQ2";
	char* lapack_str = "LAPACK";
	char* pc_str[NUM_PARAM_COMBOS] = { "" };

	fla_test_output_info("--- %s ---\n", op_str);
	fla_test_output_info("\n");
	fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS,
							params, LIN, fla_test_gerq2_experiment);
}


void fla_test_gerq2_experiment(test_params_t *params,
	integer  datatype,
	integer  p_cur,
	integer  q_cur,
	integer  pci,
	integer  n_repeats,
	double* perf,
	double* t,
	double* residual)
{
	integer m, n, cs_A;
	void *A = NULL, *A_test = NULL, *T = NULL;
	double time_min = 1e9;

	// Get input matrix dimensions.
	m = p_cur;
	n = q_cur;
	cs_A = m;

	// Create input matrix parameters
	create_matrix(datatype, &A, m, n);
	create_vector(datatype, &T, min(m,n));

	// Initialize input matrix A with random numbers
	rand_matrix(datatype, A, m, n, cs_A);

	// Make a copy of input matrix A. This is required to validate the API functionality.
	create_matrix(datatype, &A_test, m, n);
	copy_matrix(datatype, "full", m, n, A, cs_A, A_test, cs_A);

	prepare_gerq2_run(m, n, A_test, T, datatype, n_repeats, &time_min);

	// execution time
	*t = time_min;

	// performance computation
	// 2mn^2 - (2/3)n^3 flops
	if(m >= n)
		*perf = (double)((2.0 * m * n * n) - (( 2.0 / 3.0 ) * n * n * n )) / time_min / FLOPS_PER_UNIT_PERF;
	else
		*perf = (double)((2.0 * n * m * m) - (( 2.0 / 3.0 ) * m * m * m )) / time_min / FLOPS_PER_UNIT_PERF;
	if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
		*perf *= 4.0;

	// output validation
	validate_gerq2(m, n, A, A_test, T, datatype, residual);

	// Free up buffers
	free_matrix(A);
	free_matrix(A_test);
	free_vector(T);
}


void prepare_gerq2_run(integer m_A, integer n_A,
	void *A,
	void *T,
	integer datatype,
	integer n_repeats,
	double* time_min_)
{
	integer cs_A, min_A, i;
	void *A_save = NULL, *T_test = NULL, *work = NULL;
	integer info = 0;
	double time_min = 1e9, exe_time;

	cs_A = m_A;
	min_A = min(m_A, n_A);

	/* Make a copy of the input matrix A. Same input values will be passed in
	   each itertaion.*/
	create_matrix(datatype, &A_save, m_A, n_A);
	copy_matrix(datatype, "full", m_A, n_A, A, cs_A, A_save, cs_A);

	for (i = 0; i < n_repeats; ++i)
	{
		/* Restore input matrix A value and allocate memory to output buffers
		   for each iteration*/
		copy_matrix(datatype, "full", m_A, n_A, A_save, cs_A, A, cs_A);
		create_vector(datatype, &T_test, min_A);
		create_vector(datatype, &work, cs_A);

		exe_time = fla_test_clock();

		// call to gerq2 API
		invoke_gerq2(datatype, &m_A, &n_A, A, &cs_A, T_test, work, &info);

		exe_time = fla_test_clock() - exe_time;

		// Get the best execution time
		time_min = min(time_min, exe_time);

		// Make a copy of the output buffers. This is required to validate the API functionality.
		copy_vector(datatype, min_A, T_test, 1, T, 1);

		// Free up the output buffers
		free_vector(work);
		free_vector(T_test);
	}

	*time_min_ = time_min;

	free_matrix(A_save);
}


inline void invoke_gerq2(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info)
{
	switch(datatype)
	{
		case FLOAT:
		{
			sgerq2_(m, n, a, lda, tau, work, info);
			break;
		}
		
		case DOUBLE:
		{
			dgerq2_(m, n, a, lda, tau, work, info);
			break;
		}

		case COMPLEX:
		{
			cgerq2_(m, n, a, lda, tau, work, info);
			break;
		}

		case DOUBLE_COMPLEX:
		{
			zgerq2_(m, n, a, lda, tau, work, info);
			break;
		}
	}
}


void validate_gerq2(integer m_A,
	integer n_A,
	void *A,
	void *A_test,
	void *T_test,
	integer datatype,
	double* residual)
{
	void *R = NULL, *Q = NULL, *I = NULL, *work = NULL;
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
	create_matrix(datatype, &I, n_A, n_A);

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
			slaset_("full", &n_A, &n_A, &s_zero, &s_one, I, &n_A);
			sgemm_("N", "T", &n_A, &n_A, &n_A, &s_n_one, Q, &n_A, Q, &n_A, &s_one, I, &n_A);

			norm = slange_("1", &n_A, &n_A, I, &n_A, work);
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
			dlaset_("full", &n_A, &n_A, &d_zero, &d_one, I, &n_A);
			dgemm_("N", "T", &n_A, &n_A, &n_A, &d_n_one, Q, &n_A, Q, &n_A, &d_one, I, &n_A);

			norm = dlange_("1", &n_A, &n_A, I, &n_A, work);
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
			claset_("full", &n_A, &n_A, &c_zero, &c_one, I, &n_A);
			cgemm_("N", "C", &n_A, &n_A, &n_A, &c_n_one, Q, &n_A, Q, &n_A, &c_one, I, &n_A);

			norm = clange_("1", &n_A, &n_A, I, &n_A, work);
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
			zlaset_("full", &n_A, &n_A, &z_zero, &z_one, I, &n_A);
			zgemm_("N", "C", &n_A, &n_A, &n_A, &z_n_one, Q, &n_A, Q, &n_A, &z_one, I, &n_A);

			norm = zlange_("1", &n_A, &n_A, I, &n_A, work);
			resid2 = norm/(eps * (double)n_A);

			*residual = (double)max(resid1, resid2);
			break;
		}
	}

	// Free up buffers
	free_matrix( R );
	free_matrix( Q );
	free_matrix( I );
	free_vector( work );
}