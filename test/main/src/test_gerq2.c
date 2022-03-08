/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_libflame.h"
#include "test_common.h"
#include "test_prototype.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1


// Static variables.
static char* op_str = "RQ factorization with unblocked algorithm";
static char* front_str = "GERQ2";
static char* lapack_str = "LAPACK";
static char* pc_str[NUM_PARAM_COMBOS] = { "" };
static test_thresh_t thresh = { 50, 20,   // warn, pass for s
								50, 20,   // warn, pass for d
								50, 20,   // warn, pass for c
								50, 20 }; // warn, pass for z

// Local prototypes.
void fla_test_gerq2_experiment(test_params_t params, integer  datatype, integer  p_cur, integer  pci,
									integer  n_repeats, double* perf, double* t, double* residual);
void GERQ2_run(integer m_A, integer n_A, void *A, void *T, integer datatype, integer n_repeats, double* time_min_);
inline void GERQ2_API(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info);
void GERQ2_solve(integer m_A, integer n_A, void *A, void *A_test, void *T_test, int datatype, double* residual);


void fla_test_gerq2(test_params_t params)
{

	fla_test_output_info("--- %s ---\n", op_str);
	fla_test_output_info("\n");
	fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS,
							params, thresh, fla_test_gerq2_experiment);
}


void fla_test_gerq2_experiment(test_params_t params,
	integer  datatype,
	integer  p_cur,
	integer  pci,
	integer  n_repeats,
	double* perf,
	double* t,
	double* residual)
{
	integer m, n;
	void *A, *A_test, *T;
	double time_min = 1e9;

	// Determine the dimensions.
	m = p_cur;
	n = p_cur;

	// Create the matrices for the current operation.
	create_matrix(datatype, m, n, &A);
	create_matrix(datatype, m, n, &A_test);
	create_matrix(datatype, min(m, n), 1, &T);

	// Initialize the test matrices.
	rand_matrix(datatype, m, n, A);

	// Save the original matrix.
	copy_matrix(datatype, m, n, A, A_test);

	// call to API
	GERQ2_run(m, n, A_test, T, datatype, n_repeats, &time_min);

	// execution time
	*t = time_min;

	// TODO
	// performance computation will be added later
	*perf = 0;

	// output validation
	GERQ2_solve(m, n, A, A_test, T, datatype, residual);

	free_matrix(A);
	free_matrix(A_test);
	free_matrix(T);
}


void GERQ2_run(integer m_A, integer n_A,
	void *A,
	void *T,
	integer datatype,
	integer n_repeats,
	double* time_min_)
{
	integer cs_A;
	void *A_save, *T_test = NULL, *work = NULL;
	integer i;
	integer info = 0;
	double time_min = 1e9, exe_time;

	// Get column stride
	cs_A = m_A;

	// Save the original matrix.
	create_matrix(datatype, m_A, n_A, &A_save);
	copy_matrix(datatype, m_A, n_A, A, A_save);


	for (i = 0; i < n_repeats; ++i)
	{
		// reallocate memory
		free(T_test);
		free(work);
		create_matrix(datatype, min(m_A, n_A), 1, &T_test);
		create_matrix(datatype, cs_A, 1, &work);
			
		// Make a copy of the matrix
		copy_matrix(datatype, m_A, n_A, A_save, A);

		fla_start_timer();

		// call to API
		GERQ2_API(datatype, &m_A, &n_A, A, &cs_A, T_test, work, &info);

		exe_time = fla_end_timer();

		// Get the best execution time
		time_min = min(time_min, exe_time);
	}

	copy_matrix(datatype, min(m_A, n_A), 1, T_test, T);

	*time_min_ = time_min;

	free_matrix(A_save);
	free_matrix(T_test);
	free_matrix(work);
}


/*
 *  GERQ2_API calls LAPACK interface of
 *  RQ factorisation with unblocked algorithm - gerq2
 *  */
inline void GERQ2_API(integer datatype, integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info)
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


void GERQ2_solve(integer m_A,
	integer n_A,
	void *A,
	void *A_test,
	void *T_test,
	integer datatype,
	double* residual)
{

	void *b, *b_copy, *qbt, *work;
	integer cs_A, cs_b, mb, nb;
	integer  inc_qbt, inc_b;
	integer lwork, tinfo;
	integer min_m_n;

	// Compute the minimum dimension.
	min_m_n = min(m_A, n_A);

	// Get column stride.
	cs_A = m_A;

	// Create vectors to form a linear system.
	create_matrix(datatype, m_A, 1, &b);

	// Create a random vector b
	rand_matrix(datatype, m_A, 1, b);

	// copy the required dimensions of vector b.
	mb = m_A;
	nb = 1;
	cs_b = mb;
	inc_b = 1;

	// make a copy of vector b.
	create_matrix(datatype, m_A, 1, &b_copy);
	copy_matrix(datatype, m_A, 1, b, b_copy);

	// create a vector to hold the solution
	create_matrix(datatype, min_m_n, 1, &qbt);
	inc_qbt = 1;

	// Make a workspace query the first time through. This will provide us with
	// and ideal workspace size based on an internal block size.
	lwork = -1;
	create_matrix(datatype, 1, 1, &work);

	switch( datatype )
	{
		case FLOAT:
		{
			float one = 1, n_one = -1;
			float norm, norm_b;
			float eps;

			/* Solve R*x = QT*b linear equation*/
			// Rx' = b 
			strsm_( "left", "Upper", "No-Tran", "Non-unit", &mb, &nb, &one, A_test, &cs_A, b, &cs_b);

			// b' = QT*b
			sormrq_( "left", "Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );

			// allocate work buffer and calculate b'
			lwork = (integer) (*(float*)work);
			free(work);
			create_matrix(datatype, lwork, 1, &work);

			sormrq_( "left", "Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );

			free( work );

			for(integer i = 0; i < min_m_n; i++ )
				((float *)qbt)[i] = ((float *)b)[i];

			// compute norm of b
			norm_b = snrm2_(&m_A, b_copy, &inc_b);

			// compute output residue (b'-b)
			sgemv_("N", &m_A, &n_A, &one, A, &cs_A, qbt, &inc_qbt, &n_one, b_copy, &inc_b);
			
			// compute norm of norm(b'-b)
			norm = snrm2_(&m_A, b_copy, &inc_b);

			// get machine precision
			eps = slamch_("P");

			*residual = (double)(norm / (eps * norm_b * (float)m_A));

			break;
		}

		case DOUBLE:
		{
			double one = 1, n_one = -1;
			double norm, norm_b;
			double eps;

			/* Solve R*x = QT*b linear equation*/
			// Rx' = b 
			dtrsm_( "left", "Upper", "No-Tran", "Non-unit", &mb, &nb, &one, A_test, &cs_A, b, &cs_b);

			// b' = QT*b
			dormrq_( "left", "Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );

			// allocate work buffer and calculate b'
			lwork = (integer) (*(double*)work);
			free(work);
			create_matrix(datatype, lwork, 1, &work);

			dormrq_( "left", "Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );

			free( work );

			for(integer i = 0; i < min_m_n; i++ )
				((double *)qbt)[i] = ((double *)b)[i];

			// compute norm of b
			norm_b = dnrm2_(&m_A, b_copy, &inc_b);

			// compute output residue (b'-b)
			dgemv_("N", &m_A, &n_A, &one, A, &cs_A, qbt, &inc_qbt, &n_one, b_copy, &inc_b);
			
			// compute norm of norm(b'-b)
			norm = dnrm2_(&m_A, b_copy, &inc_b);

			// get machine precision
			eps = dlamch_("P");

			*residual = (double)(norm / (eps * norm_b * (double)m_A));
			
			break;
		}
		case COMPLEX:
		{
			scomplex one = {1, 0}, n_one = {-1, 0};
			float norm, norm_b;
			float eps;

			/* Solve R*x = QT*b linear equation*/
			// Rx' = b 
			ctrsm_( "left", "Upper", "No-Tran", "Non-unit", &mb, &nb, &one, A_test, &cs_A, b, &cs_b);

			// b' = QT*b
			cunmrq_( "left", "Conjugate Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );

			// allocate work buffer and calculate b'
			lwork = (integer) (*(float*)work);
			free(work);
			create_matrix(datatype, lwork, 1, &work);

			cunmrq_( "left", "Conjugate Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );

			free( work );

			for(integer i = 0; i < min_m_n; i++ )
			{
				((scomplex *)qbt)[i].real = ((scomplex *)b)[i].real;
				((scomplex *)qbt)[i].imag = ((scomplex *)b)[i].imag;
			}

			// compute norm of b
			norm_b = scnrm2_(&m_A, b_copy, &inc_b);

			// compute output residue (b'-b)
			cgemv_("N", &m_A, &n_A, &one, A, &cs_A, qbt, &inc_qbt, &n_one, b_copy, &inc_b);
			
			// compute norm of norm(b'-b)
			norm = scnrm2_(&m_A, b_copy, &inc_b);

			// get machine precision
			eps = slamch_("P");

			*residual = (double)(norm / (eps * norm_b * (float)m_A));

			break;
		}
		case DOUBLE_COMPLEX:
		{
			dcomplex one = {1, 0}, n_one = {-1, 0};
			double norm, norm_b;
			double eps;

			/* Solve R*x = QT*b linear equation*/
			// Rx' = b 
			ztrsm_( "left", "Upper", "No-Tran", "Non-unit", &mb, &nb, &one, A_test, &cs_A, b, &cs_b);

			// b' = QT*b
			zunmrq_( "left", "Conjugate Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );

			// allocate work buffer and calculate b'
			lwork = (integer) (*(double*)work);
			free(work);
			create_matrix(datatype, lwork, 1, &work);

			zunmrq_( "left", "Conjugate Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );
			free( work );

			for(integer i = 0; i < min_m_n; i++ )
			{
				((dcomplex *)qbt)[i].real = ((dcomplex *)b)[i].real;
				((dcomplex *)qbt)[i].imag = ((dcomplex *)b)[i].imag;
			}

			// compute norm of b
			norm_b = dznrm2_(&m_A, b_copy, &inc_b);

			// compute output residue (b'-b)
			zgemv_("N", &m_A, &n_A, &one, A, &cs_A, qbt, &inc_qbt, &n_one, b_copy, &inc_b);
			
			// compute norm of norm(b'-b)
			norm = dznrm2_(&m_A, b_copy, &inc_b);

			// get machine precision
			eps = dlamch_("P");

			*residual = (double)(norm / (eps * norm_b * (double)m_A));

			break;
		}
	}

	free_matrix( b );
	free_matrix( b_copy );
	free_matrix( qbt );
}