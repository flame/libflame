/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


#include "test_libflame.h"
#include "test_common.h"
#include "test_prototype.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1

// Svd_type
#define VECTORS_ALL           0
#define VECTORS_MIN_COPY      1
#define VECTORS_MIN_OVERWRITE 2
#define VECTORS_NONE          3


// Static variables
static char* op_str = "Singular value decomposition";
static char* front_str = "GESVD";
static char* lapack_str = "LAPACK";
static char* pc_str[NUM_PARAM_COMBOS] = { "" };
static test_thresh_t thresh = { 50, 20,   // warn, pass for s
								50, 20,   // warn, pass for d
								50, 20,   // warn, pass for c
								50, 20 }; // warn, pass for z

// Local prototypes.
void fla_test_gesvd_experiment(test_params_t *params, integer datatype, integer p_cur, integer pci,
									integer n_repeats, double* perf, double* t, double* residual);
void GESVD_run(int jobu, int jobv, integer m_A, integer n_A, void *A, void *s, void *U, void *V, integer datatype, integer n_repeats, double* time_min_);
void GESVD_solve(integer m_A, integer n_A, void *A, void *A_test, void *s_test, void *U_test, void *V_test, integer datatype, double* residual);
inline void GESVD_API(integer datatype, char* jobu, char* jobv, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* info);
void Param_map_svd_type( int svd_type, void* lapack_svd_type );

void fla_test_gesvd(test_params_t *params)
{
	fla_test_output_info("--- %s ---\n", op_str);
	fla_test_output_info("\n");
	fla_test_op_driver(front_str, lapack_str, NUM_PARAM_COMBOS, pc_str, NUM_MATRIX_ARGS,
							params, thresh, SVD, fla_test_gesvd_experiment);
}


void fla_test_gesvd_experiment(test_params_t *params,
	integer  datatype,
	integer  p_cur,
	integer pci,
	integer n_repeats,
	double* perf,
	double *t,
	double* residual)
{
	integer m, n;
	void *A, *U, *V, *s, *A_test;
	double time_min = 1e9;


	// Determine the dimensions.
	m = p_cur;
	n = p_cur;

	// Create the matrices for the current operation.
	create_matrix(datatype, m, n, &A);
	create_matrix(datatype, m, n, &A_test);
	create_matrix(datatype, m, m, &U);
	create_matrix(datatype, n, n, &V);

	// the datatype of s matrix must be real
	create_realtype_matrix(datatype, min(m, n), 1, &s);

	// Initialize the test matrices
	rand_matrix(datatype, m, n, A);

	// Save the original matrix.
	copy_matrix(datatype, m, n, A, A_test);

	// call to API
	GESVD_run(VECTORS_ALL, VECTORS_ALL, m, n, A_test, s, U, V, datatype, n_repeats, &time_min);

	// execution time
	*t = time_min;

	// TODO
	// performance computation will be added later
	*perf = 0;

	// output validation
	GESVD_solve(m, n, A, A_test, s, U, V, datatype, residual);

	// free all the matrixes
	free(A);
	free(A_test);
	free(U);
	free(V);
	free(s);
}



void GESVD_run(int jobu, int jobv,
	integer m_A, integer n_A,
	void *A,
	void *s,
	void *U,
	void *V,
	integer datatype,
	integer n_repeats,
	double* time_min_)
{
	char blas_jobu, blas_jobv;
	integer cs_A, cs_U, cs_V;
	integer min_m_n;
	void *A_save, *U_test, *V_test, *s_test, *work, *rwork;
	integer lwork, lrwork;
	integer i;
	integer info = 0;
	double time_min = 1e9, exe_time;


	cs_A = m_A;
	cs_U = m_A;
	cs_V = n_A;
	min_m_n = min(m_A, n_A);

	// Save the original matrix.
	create_matrix(datatype, m_A, n_A, &A_save);
	copy_matrix(datatype, m_A, n_A, A, A_save);

	// Allocate the rwork array up front since its size is not dependent on
	// internal block sizes.
	lrwork = 5 * min_m_n;
	if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
		create_realtype_matrix(datatype, lrwork, 1, &rwork);
	else
		rwork = NULL;

	Param_map_svd_type(jobu, &blas_jobu);
	Param_map_svd_type(jobv, &blas_jobv);

	// Make a workspace query the first time through. This will provide us with
	// and ideal workspace size based on an internal block size.
	lwork = -1;
	create_matrix(datatype, 1, 1, &work);

	create_matrix(datatype, m_A, m_A, &U_test);
	create_matrix(datatype, n_A, n_A, &V_test);
	create_realtype_matrix(datatype, min_m_n, 1, &s_test);

	GESVD_API(datatype, &blas_jobu, &blas_jobv, &m_A, &n_A, A, &cs_A, s_test, U_test, &cs_U, V_test, &cs_V, work, &lwork, rwork, &info);

	lwork = get_work_value( datatype, work );

	for (i = 0; i < n_repeats; ++i)
	{
		// reallocate memory
		free(U_test);
		free(V_test);
		free(s_test);
		free(work);
		free(rwork);
		create_matrix(datatype, m_A, m_A, &U_test);
		create_matrix(datatype, n_A, n_A, &V_test);
		create_realtype_matrix(datatype, min_m_n, 1, &s_test);
		create_matrix(datatype, lwork, 1, &work);

		if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
			create_realtype_matrix(datatype, lrwork, 1, &rwork);
		else
			rwork = NULL;

		// Copy original input data
		copy_matrix(datatype, m_A, n_A, A_save, A);

		fla_start_timer();

		// call to API
		GESVD_API(datatype, &blas_jobu, &blas_jobv, &m_A, &n_A, A, &cs_A, s_test, U_test, &cs_U, V_test, &cs_V, work, &lwork, rwork, &info);
		
		exe_time = fla_end_timer();

		// Get the best execution time
		time_min = min(time_min, exe_time);
	}

	copy_matrix(datatype, m_A, m_A, U_test, U);
	copy_matrix(datatype, n_A, n_A, V_test, V);
	copy_realtype_matrix(datatype, min_m_n, 1, s_test, s);

	*time_min_ = time_min;

	free(work);
	free(A_save);
	free(U_test);
	free(V_test);
	free(s_test);
	if ( datatype == COMPLEX || datatype == DOUBLE_COMPLEX)
		free(rwork);
}


/*
 *  GESVD_API calls LAPACK interface of
 *  Singular value decomposition - gesvd
 *  */
inline void GESVD_API(integer datatype, char* jobu, char* jobv, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* info)
{
	switch(datatype)
	{
		case FLOAT:
		{
			sgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
			break;
		}
		
		case DOUBLE:
		{
			dgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, info);
			break;
		}

		case COMPLEX:
		{
			cgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
			break;
		}

		case DOUBLE_COMPLEX:
		{
			zgesvd_(jobu, jobv, m, n, a, lda, s, u, ldu, vt, ldvt, work, lwork, rwork, info);
			break;
		}
	}
}


void GESVD_solve(integer m_A, integer n_A,
	void *A,
	void *A_test,
	void *s_test,
	void *U_test,
	void *V_test,
	integer datatype,
	double* residual)
{
	void *x, *b, *w;
	integer cs_A, cs_V, cs_U, cs_w, rs_w;
	integer m_U, n_U, m_V, n_V, m_s, n_s;
	integer inc_x, inc_b, inc_w, inc_s;

	// Get the matrix dimensions
	// matrix A
	cs_A = m_A;
	// matrix U
	m_U = m_A;
	n_U = m_A;
	cs_U = m_A;
	// matrix V
	m_V = n_A;
	n_V = n_A;
	cs_V = n_A;
	// matrix s
	m_s = n_A;
	n_s = 1;
	inc_s = 1;


	// Create vectors to form a linear system.
	create_matrix(datatype, n_A, 1, &x);
	inc_x = 1;
	create_matrix(datatype, m_A, 1, &b);
	inc_b = 1;
	create_matrix(datatype, min(m_A, n_A), 1, &w);
	cs_w = n_A;
	rs_w = 1;
	inc_w = 1;

	// Initialize the matrices
	rand_matrix(datatype, n_A, 1, x);

	switch( datatype )
	{
		case FLOAT:
		{
			float norm, norm_b;
			float eps;
			float one = 1, zero = 0, n_one = -1;

			// compute b = A*x
			sgemv_( "N", &m_A, &n_A, &one, A, &cs_A, x, &inc_x, &zero, b, &inc_b );

			// compute norm of b
			norm_b = snrm2_( &m_A, b, &inc_b);

			//compute b' = U*SIGMA*V**T
			sgemv_( "N", &m_V, &n_V, &one, V_test, &cs_V, x, &inc_x, &zero, w, &inc_w );
			diagmv( datatype, m_s, n_s, s_test, inc_s, w, rs_w, cs_w );

			// compute output residue (b'-b)
			sgemv_( "N", &m_U, &n_U, &one, U_test, &cs_U, w, &inc_w, &n_one, b, &inc_b );
			
			// compute norm(b'-b)
			norm = snrm2_( &m_A, b, &inc_b);

			// get machine precision
			eps = slamch_("P");

			*residual = (double)(norm / (eps * norm_b * (float)m_A));

			break;
		}

		case DOUBLE:
		{
			double norm, norm_b;
			double eps;
			double one = 1, zero = 0, n_one = -1;

			// compute b = A*x
			dgemv_( "N", &m_A, &n_A, &one, A, &cs_A, x, &inc_x, &zero, b, &inc_b );

			// compute norm of b
			norm_b = dnrm2_( &m_A, b, &inc_b);

			//compute b' = U*SIGMA*V**T
			dgemv_( "N", &m_V, &n_V, &one, V_test, &cs_V, x, &inc_x, &zero, w, &inc_w );
			diagmv( datatype, m_s, n_s, s_test, inc_s, w, rs_w, cs_w );
			
			// compute output residue (b'-b)
			dgemv_( "N", &m_U, &n_U, &one, U_test, &cs_U, w, &inc_w, &n_one, b, &inc_b );
			
			// compute norm(b'-b)
			norm = dnrm2_( &m_A, b, &inc_b);

			// get machine precision
			eps = dlamch_("P");

			*residual = (double)(norm / (eps * norm_b * (double)m_A));

			break;
		}

		case COMPLEX:
		{
			float norm, norm_b;
			float eps;
			scomplex zero = {0, 0}, one = {1, 0}, n_one = {-1, 0};

			// compute b = A*x
			cgemv_( "N", &m_A, &n_A, &one, A, &cs_A, x, &inc_x, &zero, b, &inc_b );

			// compute norm of b
			norm_b = scnrm2_( &m_A, b, &inc_b);

			//compute b' = U*SIGMA*V**T
			cgemv_( "N", &m_V, &n_V, &one, V_test, &cs_V, x, &inc_x, &zero, w, &inc_w );
			diagmv( datatype, m_s, n_s, s_test, inc_s, w, rs_w, cs_w );
			
			// compute output residue (b'-b)
			cgemv_( "N", &m_U, &n_U, &one, U_test, &cs_U, w, &inc_w, &n_one, b, &inc_b );
			
			// compute norm(b'-b)
			norm = scnrm2_( &m_A, b, &inc_b);

			// get machine precision
			eps = slamch_("P");

			*residual = (double)(norm / (eps * norm_b * (double)m_A));

			break;
		}

		case DOUBLE_COMPLEX:
		{

			double norm, norm_b;
			double eps;
			dcomplex zero = {0, 0}, one = {1, 0}, n_one = {-1, 0};

			// compute b = A*x
			zgemv_( "N", &m_A, &n_A, &one, A, &cs_A, x, &inc_x, &zero, b, &inc_b );

			// compute norm of b
			norm_b = dznrm2_( &m_A, b, &inc_b);

			//compute b' = U*SIGMA*V**T
			zgemv_( "N", &m_V, &n_V, &one, V_test, &cs_V, x, &inc_x, &zero, w, &inc_w );
			diagmv( datatype, m_s, n_s, s_test, inc_s, w, rs_w, cs_w );
			
			// compute output residue (b'-b)
			zgemv_( "N", &m_U, &n_U, &one, U_test, &cs_U, w, &inc_w, &n_one, b, &inc_b );
			
			// compute norm(b'-b)
			norm = dznrm2_( &m_A, b, &inc_b);

			// get machine precision
			eps = dlamch_("P");

			*residual = (double)(norm / (eps * norm_b * (double)m_A));

			break;
		}
	}

	free(x);
	free(b);
	free(w);
}


void Param_map_svd_type( int svd_type, void* lapack_svd_type )
{
	if( svd_type == VECTORS_ALL )
	{
		*( ( char* ) lapack_svd_type ) = 'A';
	}
	else if ( svd_type == VECTORS_MIN_COPY )
	{
		*( ( char* ) lapack_svd_type ) = 'S';
	}
	else if ( svd_type == VECTORS_MIN_OVERWRITE )
	{
		*( ( char* ) lapack_svd_type ) = 'O';
	}
	else if ( svd_type == VECTORS_NONE )
	{
		*( ( char* ) lapack_svd_type ) = 'N';
	}
	else
	{
		*( ( char* ) lapack_svd_type ) = 'A';
	}
}