/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_libflame.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1

/* Local prototypes */
void fla_test_getrf_experiment(test_params_t *params, integer  datatype, integer  p_cur, integer  q_cur, integer pci,
									integer n_repeats, double* perf, double* t, double* residual);
void prepare_getrf_run(integer m_A, integer n_A, void *A, integer* ipiv, integer datatype, integer n_repeats, double* time_min_);
inline void invoke_getrf(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *ipiv, integer *info);

void fla_test_getrf(test_params_t *params)
{
	char* op_str = "LU factorization";
	char* front_str = "GETRF";
	char* lapack_str = "LAPACK";
	char* pc_str[NUM_PARAM_COMBOS] = { "" };

	fla_test_output_info("--- %s ---\n", op_str);
	fla_test_output_info("\n");
	fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS,
							params, LIN, fla_test_getrf_experiment);

}


void fla_test_getrf_experiment(test_params_t *params,
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
	void* IPIV;
	void *A, *A_test;
	double time_min = 1e9;


	/* Determine the dimensions*/
	m = p_cur;
	n = q_cur;
	cs_A = m;
	/* Create the matrices for the current operation*/
	create_matrix(datatype, &A, m, n);

	create_vector(INTEGER, &IPIV, min(m, n));
	/* Initialize the test matrices*/
	rand_matrix(datatype, A, m, n, cs_A);

	/* Save the original matrix*/
	create_matrix(datatype, &A_test, m, n);
	copy_matrix(datatype, "full", m, n, A, cs_A, A_test, cs_A);

	/* call to API */
	prepare_getrf_run(m, n, A_test, IPIV, datatype, n_repeats, &time_min);

	/* execution time */
	*t = time_min;

	/* performance computation */
	/* 2mn^2 - (2/3)n^3 flops */
	*perf = (double)((2.0 * m * n * n) - ((2.0 / 3.0) * n * n * n)) / time_min / FLOPS_PER_UNIT_PERF;
	if (datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
		*perf *= 4.0;

	/* output validation */
	validate_getrf(m, n, A, A_test, IPIV, datatype, residual);

	/* Free up the buffers */
	free_matrix(A);
	free_matrix(A_test);
	free_vector(IPIV);
}


void prepare_getrf_run(integer m_A,
	integer n_A,
	void* A,
	integer* IPIV,
	integer datatype,
	integer n_repeats,
	double* time_min_)
{
	integer cs_A;
	integer i;
	void *A_save;
	integer info = 0;
	double time_min = 1e9, exe_time;

	/* Get column stride */
	cs_A = m_A;
	/* Save the original matrix */
	create_matrix(datatype, &A_save, m_A, n_A);
	copy_matrix(datatype, "full", m_A, n_A, A, cs_A, A_save, cs_A);

	for (i = 0; i < n_repeats; ++i)
	{

		/* Copy original input data */
		copy_matrix(datatype, "full", m_A, n_A, A, cs_A, A_save, cs_A);

		exe_time = fla_test_clock();

		/*  call to API */
		invoke_getrf(datatype, &m_A, &n_A, A_save, &cs_A, IPIV, &info);

		exe_time = fla_test_clock() - exe_time;

		/* Get the best execution time */
		time_min = min(time_min, exe_time);
	}

	*time_min_ = time_min;
	/*  Save the AFACT to matrix A */
	copy_matrix(datatype, "full", m_A, n_A, A_save, cs_A, A, cs_A);
	free_matrix(A_save);

}


/*
 *  GETRF_API calls LAPACK interface of
 *  Singular value decomposition - gesvd
 *  */
inline void invoke_getrf(integer datatype, integer *m, integer *n, void *a, integer *lda, integer *ipiv, integer *info)
{
	switch(datatype)
	{
		case FLOAT:
		{
			sgetrf_(m, n, a, lda, ipiv, info);
			break;
		}
		
		case DOUBLE:
		{
			dgetrf_(m, n, a, lda, ipiv, info);
			break;
		}

		case COMPLEX:
		{
			cgetrf_(m, n, a, lda, ipiv, info);
			break;
		}

		case DOUBLE_COMPLEX:
		{
			zgetrf_(m, n, a, lda, ipiv, info);
			break;
		}
	}
}

