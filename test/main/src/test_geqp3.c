/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#include "test_libflame.h"
#include "test_common.h"
#include "test_prototype.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1


// Static variables
static char* op_str = "QR factorization with column pivoting";
static char* front_str = "GEQP3";
static char* lapack_str = "LAPACK";
static char* pc_str[NUM_PARAM_COMBOS] = { ""};
static test_thresh_t thresh = { 50, 20,   // warn, pass for s
								50, 20,   // warn, pass for d
								50, 20,   // warn, pass for c
								50, 20 }; // warn, pass for z


// Local prototypes.
void fla_test_geqp3_experiment(test_params_t *params, integer datatype, integer  p_cur, integer  pci, integer  n_repeats,
									double* perf, double* t,double* residual);
void GEQP3_run(integer m_A, integer n_A, void *A, void *j, void *T, integer datatype, integer n_repeats, double* time_min_);
inline void GEQP3_API(integer datatype, integer* m, integer* n, void* a, integer* lda, integer* jpvt, void* tau, void* work, void* rwork, integer* lwork, integer* info);
void GEQP3_solve(integer m_A, integer n_A, void *A, void *A_test, void *j_test, void *T_test, integer datatype, double* residual);


void fla_test_geqp3(test_params_t *params)
{
	fla_test_output_info("--- %s ---\n", op_str);
	fla_test_output_info("\n");
	fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS,
							params, thresh, LIN, fla_test_geqp3_experiment);
}


void fla_test_geqp3_experiment(test_params_t *params,
	int  datatype,
	integer  p_cur,
	integer  pci,
	integer  n_repeats,
	double* perf,
	double* t,
	double* residual)
{
	integer m, n;
	void *A, *T, *A_test, *j;
	double time_min = 1e9;

	// Determine the dimensions.
	m = p_cur;
	n = p_cur;

	// Create the matrices for the current operation.
	create_matrix(datatype, m, n, &A);
	create_matrix(datatype, m, n, &A_test);
	create_matrix(datatype, min(m, n), 1, &T);

	create_matrix(INT, n, 1, &j);

	// Initialize the test matrices.
	rand_matrix(datatype, m, n, A);

	// Save the original matrix.
	copy_matrix(datatype, m, n, A, A_test);

	// call to API
	GEQP3_run(m, n, A_test, j, T, datatype, n_repeats, &time_min);

	// execution time
	*t = time_min;

	// TODO
	// performance computation will be added later
	*perf = 0;

	// output validation
	GEQP3_solve(m, n, A, A_test, j, T, datatype, residual);

	free(A);
	free(A_test);
	free(T);
	free(j);
}


void GEQP3_run(integer m_A, integer n_A,
	void *A,
	void *j,
	void *T,
	integer datatype,
	integer n_repeats,
	double* time_min_)
{
	integer cs_A;
	integer lwork, lrwork;
	void *A_save, *j_test, *T_test, *work, *rwork;
	integer i;
	integer info = 0;
	double time_min = 1e9, exe_time;

	// Get column stride
	cs_A = m_A;

	// Save the original matrix.
	create_matrix(datatype, m_A, n_A, &A_save);
	copy_matrix(datatype, m_A, n_A, A, A_save);

	// Allocate the rwork array up front since its size is not dependent on internal block sizes.
	lrwork = 2 * n_A;
	if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
		create_realtype_matrix(datatype, lrwork, 1, &rwork);
	else
		rwork = NULL;

	// Make a workspace query the first time through. This will provide us with
	// and ideal workspace size based on an internal block size.
	lwork = -1;
	create_matrix(datatype, 1, 1, &work);

	create_matrix(INT, n_A, 1, &j_test);
	create_matrix(datatype, min(m_A, n_A), 1, &T_test);

	GEQP3_API(datatype, &m_A, &n_A, A, &cs_A, j_test, T_test, work, rwork, &lwork, &info);

	lwork = get_work_value( datatype, work );

	for (i = 0; i < n_repeats; ++i)
	{
		// reallocate memory
		free(j_test);
		free(T_test);
		free(work);
		free(rwork);
		create_matrix(INT, n_A, 1, &j_test);
		create_matrix(datatype, min(m_A, n_A), 1, &T_test);
		create_matrix(datatype, lwork, 1, &work);

		if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
			create_realtype_matrix(datatype, lrwork, 1, &rwork);
		else
			rwork = NULL;

		// Make a copy of the matrix
		copy_matrix(datatype, m_A, n_A, A_save, A);

		fla_start_timer();

		// call to API
		GEQP3_API(datatype, &m_A, &n_A, A, &cs_A, j_test, T_test, work, rwork, &lwork, &info);

		exe_time = fla_end_timer();

		// Get the best execution time
		time_min = min(time_min, exe_time);
	}

	copy_matrix(INT, n_A, 1, j_test, j);
	copy_matrix(datatype, min(m_A, n_A), 1, T_test, T);

	*time_min_ = time_min;

	free(A_save);
	free(j_test);
	free(T_test);
	free(work);
	if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
		free(rwork);
}


/*
 *  GEQP3 calls LAPACK interface of
 *  QR factorisation with column pivot - geqp3
 *  */
inline void GEQP3_API(integer datatype, integer* m, integer* n, void* a, integer* lda, integer* jpvt, void* tau, void* work, void* rwork, integer* lwork, integer* info)
{
	switch(datatype)
	{
		case FLOAT:
		{
			sgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
			break;
		}
		
		case DOUBLE:
		{
			dgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, info);
			break;
		}

		case COMPLEX:
		{
			cgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
			break;
		}

		case DOUBLE_COMPLEX:
		{
			zgeqp3_(m, n, a, lda, jpvt, tau, work, lwork, rwork, info);
			break;
		}
	}
}


void GEQP3_solve(integer m_A,
	integer n_A,
	void *A,
	void *A_test,
	void *j_test,
	void *T_test,
	int datatype,
	double* residual)
{
	
	void *qbt, *b, *y, *b_copy, *work;
	integer cs_A, cs_b, mb, nb;
	integer  inc_qbt, inc_b;
	integer lwork, tinfo;
	integer min_m_n;

	// Compute the minimum dimension.
	min_m_n = min(m_A, n_A);

	// Get column stride and row stride
	cs_A = m_A;

	// Create vectors to form a linear system.
	create_matrix(datatype, m_A, 1, &b);
	create_matrix(datatype, n_A, 1, &y);

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

			// Solve R*x = QT*b linear equation
			// b' = QT*b
			sormqr_( "Left", "Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );

			// allocate work buffer and calculate b'
			lwork = (integer) (*(float*)work);
			free(work);
			create_matrix(datatype, lwork, 1, &work);
			
			sormqr_( "Left", "Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );
			
			free( work );
			
			// Triangular Solve to find x in Ax=b
			strtrs_( "Upper", "No-Tran", "Non-unit", &min_m_n, &nb, A_test, &cs_A, b, &cs_b, &tinfo );
			
			for(int i = 0; i < min_m_n; i++ )
				((float *)qbt)[((integer *)j_test)[i]-1] = ((float *)b)[i];

			// compute norm of b
			norm_b = snrm2_(&m_A, b_copy, &inc_b);

			// compute output residue (b'-b)
			sgemv_("N", &m_A, &n_A, &one, A, &cs_A, qbt, &inc_qbt, &n_one, b_copy, &inc_b);
			
			// compute norm of norm(b'-b)
			norm = snrm2_(&m_A, b_copy, &inc_b);

			// get machine precision
			eps = slamch_("P");

			*residual = (double)(norm / (eps * norm_b * (double)m_A));

			break;
		}
		case DOUBLE:
		{
			double one = 1, n_one = -1;
			double norm, norm_b;
			double eps;

			// Solve R*x = QT*b linear equation
			// b' = QT*b
			dormqr_( "Left", "Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );
			
			// allocate work buffer and calculate b'
			lwork = (integer) (*(double*)work);
			free(work);
			create_matrix(datatype, lwork, 1, &work);
			
			dormqr_( "Left", "Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );
			
			free( work );
			
			// Triangular Solve to find x in Ax=b
			dtrtrs_( "Upper", "No-Tran", "Non-unit", &min_m_n, &nb, A_test, &cs_A, b, &cs_b, &tinfo );
			
			for(int i = 0; i < min_m_n; i++ )
				((double *)qbt)[((integer *)j_test)[i]-1] = ((double *)b)[i];

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

			// Solve R*x = QT*b linear equation
			// b' = QT*b	
			cunmqr_( "Left", "Conjugate Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );
			
			// allocate work buffer and calculate b'
			lwork = (integer) (*(float*)work);
			free(work);
			create_matrix(datatype, lwork, 1, &work);
			
			cunmqr_( "Left", "Conjugate Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );
			
			free( work );
			
			// Triangular Solve to find x in Ax=b
			ctrtrs_( "Upper", "No-Tran", "Non-unit", &min_m_n, &nb, A_test, &cs_A, b, &cs_b, &tinfo );
			
			for(int i = 0; i < min_m_n; i++ )
			{
				((scomplex *)qbt)[((integer *)j_test)[i]-1].real = ((scomplex *)b)[i].real;
				((scomplex *)qbt)[((integer *)j_test)[i]-1].imag = ((scomplex *)b)[i].imag;
			}

			// compute norm of b
			norm_b = scnrm2_(&m_A, b_copy, &inc_b);

			// compute output residue (b'-b)
			cgemv_("N", &m_A, &n_A, &one, A, &cs_A, qbt, &inc_qbt, &n_one, b_copy, &inc_b);
			
			// compute norm of norm(b'-b)
			norm = scnrm2_(&m_A, b_copy, &inc_b);

			// get machine precision
			eps = slamch_("P");

			*residual = (double)(norm / (eps * norm_b * (double)m_A));

			break;
		}
		case DOUBLE_COMPLEX:
		{
			dcomplex one = {1, 0}, n_one = {-1, 0};
			double norm, norm_b;
			double eps;

			// Solve R*x = QT*b linear equation
			// b' = QT*b
			zunmqr_( "Left", "Conjugate Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );
			
			// allocate work buffer and calculate b'
			lwork = (integer) (*(double*)work);
			free(work);
			create_matrix(datatype, lwork, 1, &work);
			
			zunmqr_( "Left", "Conjugate Transpose", &mb, &nb, &min_m_n, A_test, &cs_A, T_test, b, &cs_b, work, &lwork, &tinfo );
			
			free( work );
			// Triangular Solve to find x in Ax=b
			ztrtrs_( "Upper", "No-Tran", "Non-unit", &min_m_n, &nb, A_test, &cs_A, b, &cs_b, &tinfo );
			
			for(int i = 0; i < min_m_n; i++ )
			{
				((dcomplex *)qbt)[((integer *)j_test)[i]-1].real = ((dcomplex *)b)[i].real;
				((dcomplex *)qbt)[((integer *)j_test)[i]-1].imag = ((dcomplex *)b)[i].imag;
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
	
	free( b );
	free( b_copy );
	free( qbt );
	free( y );
}
