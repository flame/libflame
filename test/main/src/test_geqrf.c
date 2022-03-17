/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_libflame.h"
#include "test_common.h"
#include "test_prototype.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1


// Local prototypes.
void fla_test_geqrf_experiment(test_params_t *params, integer datatype, integer  p_cur, integer  q_cur, integer  pci, integer  n_repeats,
									double* perf, double* t,double* residual);
void prepare_geqrf_run(integer m_A, integer n_A, void *A, void *T, integer datatype, integer n_repeats, double* time_min_);
inline void invoke_geqrf(integer datatype, integer* m, integer* n, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
void validate_geqrf(integer m_A, integer n_A, void *A, void *A_test, void *T_test, integer datatype, double* residual);



void fla_test_geqrf(test_params_t *params)
{
	char* op_str = "RQ factorization";
	char* front_str = "GEQRF";
	char* lapack_str = "LAPACK";
	char* pc_str[NUM_PARAM_COMBOS] = { ""};

	fla_test_output_info("--- %s ---\n", op_str);
	fla_test_output_info("\n");
	fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS,
							params, LIN, fla_test_geqrf_experiment);
}

void fla_test_geqrf_experiment(test_params_t *params,
	int  datatype,
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

	prepare_geqrf_run(m, n, A_test, T, datatype, n_repeats, &time_min);

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
	validate_geqrf(m, n, A, A_test, T, datatype, residual);

	// Free up the buffers
	free_matrix(A);
	free_matrix(A_test);
	free_vector(T);
}


void prepare_geqrf_run(integer m_A, integer n_A,
	void *A,
	void *T,
	integer datatype,
	integer n_repeats,
	double* time_min_)
{
	integer cs_A, min_A, i;
	void *A_save = NULL, *T_test = NULL, *work = NULL;
	integer lwork = -1, info = 0;
	double time_min = 1e9, exe_time;

	cs_A = m_A;
	min_A = min(m_A, n_A);

	/* Make a copy of the input matrix A. Same input values will be passed in
	   each itertaion.*/
	create_matrix(datatype, &A_save, m_A, n_A);
	copy_matrix(datatype, "full", m_A, n_A, A, cs_A, A_save, cs_A);

	/* Make a workspace query the first time. This will provide us with
	   and ideal workspace size based on internal block size.*/
	create_vector(datatype, &work, 1);

	// call to  geqrf API
	invoke_geqrf(datatype, &m_A, &n_A, NULL, &cs_A, NULL, work, &lwork, &info);

	// Get work size
	lwork = get_work_value( datatype, work );

	/* Output buffers will be freshly allocated for each iterations, free up 
	   the current output buffers.*/ 
	free_vector(work);

	for (i = 0; i < n_repeats; ++i)
	{
		/* Restore input matrix A value and allocate memory to output buffers
		   for each iteration*/
		copy_matrix(datatype, "full", m_A, n_A, A_save, cs_A, A, cs_A);

		// T_test vector will hold the scalar factors of the elementary reflectors.
		create_vector(datatype, &T_test, min_A);

		// Create work buffer
		create_matrix(datatype, &work, lwork, 1);

		exe_time = fla_test_clock();

		// Call to  gerqf API
		invoke_geqrf(datatype, &m_A, &n_A, A, &cs_A, T_test, work, &lwork, &info);

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


inline void invoke_geqrf(integer datatype, integer* m, integer* n, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info)
{
	switch(datatype)
	{
		case FLOAT:
		{
			sgeqrf_(m, n, a, lda, tau, work, lwork, info);
			break;
		}
		
		case DOUBLE:
		{
			dgeqrf_(m, n, a, lda, tau, work, lwork, info);
			break;
		}

		case COMPLEX:
		{
			cgeqrf_(m, n, a, lda, tau, work, lwork, info);
			break;
		}

		case DOUBLE_COMPLEX:
		{
			zgeqrf_(m, n, a, lda, tau, work, lwork, info);
			break;
		}
	}
}


void validate_geqrf(integer m_A,
	integer n_A,
	void *A,
	void *A_test,
	void *T_test,
	int datatype,
	double* residual)
{
	void *Q = NULL, *R = NULL, *I = NULL, *work = NULL;
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
	create_matrix(datatype, &I, m_A, m_A);

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
			slaset_("full", &m_A, &m_A, &s_zero, &s_one, I, &m_A);
			sgemm_("N", "T", &m_A, &m_A, &m_A, &s_n_one, Q, &m_A, Q, &m_A, &s_one, I, &m_A);

			norm = slange_("1", &m_A, &m_A, I, &m_A, work);
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
			dlaset_("full", &m_A, &m_A, &d_zero, &d_one, I, &m_A);
			dgemm_("N", "T", &m_A, &m_A, &m_A, &d_n_one, Q, &m_A, Q, &m_A, &d_one, I, &m_A);

			norm = dlange_("1", &m_A, &m_A, I, &m_A, work);
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
			claset_("full", &m_A, &m_A, &c_zero, &c_one, I, &m_A);
			cgemm_("N", "C", &m_A, &m_A, &m_A, &c_n_one, Q, &m_A, Q, &m_A, &c_one, I, &m_A);

			norm = clange_("1", &m_A, &m_A, I, &m_A, work);
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
			zlaset_("full", &m_A, &m_A, &z_zero, &z_one, I, &m_A);
			zgemm_("N", "C", &m_A, &m_A, &m_A, &z_n_one, Q, &m_A, Q, &m_A, &z_one, I, &m_A);

			norm = zlange_("1", &m_A, &m_A, I, &m_A, work);
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
