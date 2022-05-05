/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_libflame.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1

// Local prototypes.
void fla_test_gerqf_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer pci,
									integer n_repeats, double* perf, double* t, double* residual);
void prepare_gerqf_run(integer m_A, integer n_A, void *A, void *T, integer datatype, integer n_repeats, double* time_min_);
void invoke_gerqf(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);

void fla_test_gerqf(test_params_t *params)
{
	char* op_str = "RQ factorization";
	char* front_str = "GERQF";
	char* lapack_str = "LAPACK";
	char* pc_str[NUM_PARAM_COMBOS] = { "" };

	fla_test_output_info("--- %s ---\n", op_str);
	fla_test_output_info("\n");
	fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS,
							params, LIN, fla_test_gerqf_experiment);

}


void fla_test_gerqf_experiment(test_params_t *params,
	integer  datatype,
	integer  p_cur,
	integer  q_cur,
	integer pci,
	integer n_repeats,
	double* perf,
	double* t,
	double* residual)
{
	integer m, n, cs_A;
	void *A, *A_test, *T;
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

	prepare_gerqf_run(m, n, A_test, T, datatype, n_repeats, &time_min);

	// Execution time
	*t = time_min;

	// performance computation
	// 2mn^2 - (2/3)n^3 flops
	if(m >= n)
		*perf = (double)((2.0 * m * n * n) - (( 2.0 / 3.0 ) * n * n * n )) / time_min / FLOPS_PER_UNIT_PERF;
	else
		*perf = (double)((2.0 * n * m * m) - (( 2.0 / 3.0 ) * m * m * m )) / time_min / FLOPS_PER_UNIT_PERF;
	if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
		*perf *= 4.0;

	// Output validation
	validate_gerqf(m, n, A, A_test, T, datatype, residual);

	// Free up the buffers
	free_matrix(A);
	free_matrix(A_test);
	free_vector(T);
}


void prepare_gerqf_run(integer m_A,
	integer n_A,
	void *A,
	void *T,
	integer datatype,
	integer n_repeats,
	double* time_min_)
{
	integer cs_A, min_A, i;
	void *A_save, *T_test, *work;
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

	// call to  gerqf API
	invoke_gerqf(datatype, &m_A, &n_A, NULL, &cs_A, NULL, work, &lwork, &info);

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

		// create work buffer
		create_matrix(datatype, &work, lwork, 1);

		exe_time = fla_test_clock();

		// Call to  gerqf API
		invoke_gerqf(datatype, &m_A, &n_A, A, &cs_A, T_test, work, &lwork, &info);

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



void invoke_gerqf(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info)
{
	switch(datatype)
	{
		case FLOAT:
		{
			sgerqf_(m, n, a, lda, tau, work, lwork, info);
			break;
		}
		case DOUBLE:
		{
			dgerqf_(m, n, a, lda, tau, work, lwork, info);
			break;
		}
		case COMPLEX:
		{
			cgerqf_(m, n, a, lda, tau, work, lwork, info);
			break;
		}
		case DOUBLE_COMPLEX:
		{
			zgerqf_(m, n, a, lda, tau, work, lwork, info);
			break;
		}
	}
}
