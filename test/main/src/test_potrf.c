/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_libflame.h"
#include "test_common.h"
#include "test_prototype.h"

#define NUM_PARAM_COMBOS 2
#define NUM_MATRIX_ARGS  1

/* Local prototypes.*/
void fla_test_potrf_experiment(test_params_t *params, integer datatype, integer  p_cur, integer  q_cur, integer  pci, integer  n_repeats,double* perf, double* time_min, double* residual);
void prepare_potrf_run(char* uplo, integer m, void *A, integer datatype, integer n_repeats, double* time_min_);
inline void invoke_potrf(char* uplo, integer datatype, integer* m, void* a, integer* lda, integer* info);
void validate_potrf(char *uplo, integer m, void *A, void *A_test, integer datatype, double* residual);

void fla_test_potrf(test_params_t *params)
{
	char* op_str = "Cholesky factorization";
	char* front_str = "POTRF";
	char* lapack_str = "LAPACK";
	char* pc_str[NUM_PARAM_COMBOS] = {"U","L"};

	fla_test_output_info("--- %s ---\n", op_str);
	fla_test_output_info("\n");
	fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS,
							params, LIN, fla_test_potrf_experiment);
}

void fla_test_potrf_experiment(test_params_t *params,
	integer  datatype,
	integer  p_cur,
	integer  q_cur,
	integer  pci,
	integer  n_repeats,
	double* perf,
	double* time_min,
	double* residual)
{
	integer m;
	void *A = NULL, *A_test = NULL;
	char uplo = params->lin_solver_paramslist[pci].Uplo;

	/* Get input matrix dimensions */
	m = p_cur;

	/* Create input matrix parameters */
	create_matrix(datatype, &A, m, m);

	/* Initialize input symmetric positive definite matrix A */
	rand_spd_matrix(datatype, &uplo, &A, m, m);

	/* Make a copy of input matrix A. This is required to validate the API functionality */
	create_matrix(datatype, &A_test, m, m);
	copy_matrix(datatype, "full", m, m, A, m, A_test, m);

	prepare_potrf_run(&uplo, m, A_test, datatype, n_repeats, time_min);

	/* Compute the performance of the best experiment repeat */
	/* (1/3)m^3 for real and (4/3)m^3 for complex*/
        *perf = (double)(1.0 / 3.0 * m * m * m) / *time_min / FLOPS_PER_UNIT_PERF;
        if(datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
		*perf *= 4.0;

	validate_potrf(&uplo, m, A, A_test, datatype, residual);

	free_matrix(A);
	free_matrix(A_test);
}

void prepare_potrf_run(char* uplo, integer m,
	void *A,
	integer datatype,
	integer n_repeats,
	double* time_min_)
{
	void *A_save = NULL;
	double time_min = 1e9, exe_time;
	integer i, info = 0;

	/* Make a copy of the input matrix A. Same input values will be passed in
	   each itertaion.*/
	create_matrix(datatype, &A_save, m, m);
	copy_matrix(datatype, "full", m, m, A, m, A_save, m);

	for (i = 0; i < n_repeats; ++i)
	{
		/* Restore input matrix A value and allocate memory to output buffers
		for each iteration */
		copy_matrix(datatype, "full", m, m, A_save, m, A, m);
		exe_time = fla_test_clock();
		invoke_potrf(uplo, datatype, &m, A, &m, &info);
		exe_time = fla_test_clock() - exe_time;
		/* Get the best execution time */
		time_min = min(time_min, exe_time);
	}

	*time_min_ = time_min;
	free_matrix(A_save);
}

inline void invoke_potrf(char* uplo, integer datatype, integer* m, void* a, integer* lda, integer* info)
{
	switch(datatype)
	{
		case FLOAT:
		{
			spotrf_(uplo, m, a, lda, info);
			break;
		}
		case DOUBLE:
		{
			dpotrf_(uplo, m, a, lda, info);
			break;
		}
		case COMPLEX:
		{
			cpotrf_(uplo, m, a, lda, info);
			break;
		}
		case DOUBLE_COMPLEX:
		{
			zpotrf_(uplo, m, a, lda, info);
			break;
		}
	}
}

void validate_potrf(char *uplo, integer m, void *A, void *A_test, integer datatype, double* residual)
{
	void *b = NULL, *x = NULL;
	void *x_test = NULL, *b_test = NULL;
	void *work = NULL;
	void *A_save = NULL;
	void *buff_A = NULL, *buff_B = NULL;
	integer info;
	integer	nrhs = 1, incx = 1, incy = 1;
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
