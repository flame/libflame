
#include "FLAME.h"

FLA_Error FLA_Hevd_lv_var1_components( dim_t n_iter_max, FLA_Obj A, FLA_Obj l, dim_t k_accum, dim_t b_alg,
                                       double* dtime_tred, double* dtime_tevd, double* dtime_appq )
{
	FLA_Error    r_val = FLA_SUCCESS;
	FLA_Uplo     uplo = FLA_LOWER_TRIANGULAR;
	FLA_Datatype dt;
	FLA_Datatype dt_comp;
	FLA_Datatype dt_real;
	FLA_Obj      T, r, d, e, G;
	dim_t        mn_A;
	dim_t        n_G = k_accum;
	double       dtime_temp;

	mn_A    = FLA_Obj_length( A );
	dt      = FLA_Obj_datatype( A );
	dt_comp = FLA_Obj_datatype_proj_to_complex( A );
	dt_real = FLA_Obj_datatype_proj_to_real( A );

	// If the matrix is a scalar, then the EVD is easy.
	if ( mn_A == 1 )
	{
		FLA_Copy( A, l );
		FLA_Set( FLA_ONE, A );

		return FLA_SUCCESS;
	}

	// Create a matrix to hold block Householder transformations.
	FLA_Tridiag_UT_create_T( A, &T );

	// Create a vector to hold the realifying scalars.
	FLA_Obj_create( dt,      mn_A, 1, 0, 0, &r );

	// Create vectors to hold the diagonal and sub-diagonal.
	FLA_Obj_create( dt_real, mn_A,     1, 0, 0, &d );
	FLA_Obj_create( dt_real, mn_A-1,   1, 0, 0, &e );
	FLA_Obj_create( dt_comp, mn_A-1, n_G, 0, 0, &G );


  dtime_temp = FLA_Clock();
  {
	// Reduce the matrix to tridiagonal form.
	FLA_Tridiag_UT( uplo, A, T );
  }
  *dtime_tred = FLA_Clock() - dtime_temp;


	// Apply scalars to rotate elements on the sub-diagonal to the real domain.
	FLA_Tridiag_UT_realify( uplo, A, r );

	// Extract the diagonal and sub-diagonal from A.
	FLA_Tridiag_UT_extract_diagonals( uplo, A, d, e );


  dtime_temp = FLA_Clock();
  {
	// Form Q, overwriting A.
	FLA_Tridiag_UT_form_Q( uplo, A, T );
  }
  *dtime_appq = FLA_Clock() - dtime_temp;


	// Apply the scalars in r to Q.
	FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE, r, A );


  dtime_temp = FLA_Clock();
  {
	// Perform an eigenvalue decomposition on the tridiagonal matrix.
	r_val = FLA_Tevd_v_opt_var1( n_iter_max, d, e, G, A, b_alg );
  }
  *dtime_tevd = FLA_Clock() - dtime_temp;


	// Copy the converged eigenvalues to the output vector.
	FLA_Copy( d, l );

	// Sort the eigenvalues and eigenvectors in ascending order.
	FLA_Sort_evd( FLA_FORWARD, l, A );

//FLA_Obj_show( "var1: d", l, "%22.10e", "" );
//FLA_Obj_show( "var1: A", A, "%8.1e", "" );

	FLA_Obj_free( &T );
	FLA_Obj_free( &r );
	FLA_Obj_free( &d );
	FLA_Obj_free( &e );
	FLA_Obj_free( &G );

	return r_val;
}

