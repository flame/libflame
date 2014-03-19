/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Hevd_lv_unb_var1( dim_t n_iter_max, FLA_Obj A, FLA_Obj l, dim_t k_accum, dim_t b_alg )
{
	FLA_Uplo     uplo = FLA_LOWER_TRIANGULAR;
	FLA_Datatype dt;
	FLA_Datatype dt_real;
	FLA_Datatype dt_comp;
	FLA_Obj      scale, T, r, d, e, G;
	dim_t        mn_A;
	dim_t        n_G = k_accum;
	FLA_Error    r_val;

	mn_A    = FLA_Obj_length( A );
	dt      = FLA_Obj_datatype( A );
	dt_real = FLA_Obj_datatype_proj_to_real( A );
	dt_comp = FLA_Obj_datatype_proj_to_complex( A );

	// Make sure the matrix is column-stored.
	if ( FLA_Obj_row_stride( A ) != 1 )
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

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
	FLA_Obj_create( dt,      mn_A,     1, 0, 0, &r );

	// Create vectors to hold the diagonal and sub-diagonal.
	FLA_Obj_create( dt_real, mn_A,     1, 0, 0, &d );
	FLA_Obj_create( dt_real, mn_A-1,   1, 0, 0, &e );
	FLA_Obj_create( dt_comp, mn_A-1, n_G, 0, 0, &G );

	// Create a real scaling factor.
	FLA_Obj_create( dt_real, 1, 1, 0, 0, &scale );

	// Compute a scaling factor; If none is needed, sigma will be set to one.
	FLA_Hevd_compute_scaling( uplo, A, scale );

	// Scale the matrix if scale is non-unit.
	if ( !FLA_Obj_equals( scale, FLA_ONE ) )
		FLA_Scalr( uplo, scale, A );

	// Reduce the matrix to tridiagonal form.
	FLA_Tridiag_UT( uplo, A, T );

	// Apply scalars to rotate elements on the sub-diagonal to the real domain.
	FLA_Tridiag_UT_realify( uplo, A, r );

	// Extract the diagonal and sub-diagonal from A.
	FLA_Tridiag_UT_extract_real_diagonals( uplo, A, d, e );

	// Form Q, overwriting A.
	FLA_Tridiag_UT_form_Q( uplo, A, T, A );

	// Apply the scalars in r to Q.
	FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE, r, A );

	// Perform an eigenvalue decomposition on the tridiagonal matrix.
	r_val = FLA_Tevd_v_opt_var1( n_iter_max, d, e, G, A, b_alg );

	// Copy the converged eigenvalues to the output vector.
	FLA_Copy( d, l );

	// Sort the eigenvalues and eigenvectors in ascending order.
	FLA_Sort_evd( FLA_FORWARD, l, A );

	// If the matrix was scaled, rescale the eigenvalues.
	if ( !FLA_Obj_equals( scale, FLA_ONE ) )
		FLA_Inv_scal( scale, l );

	FLA_Obj_free( &scale );
	FLA_Obj_free( &T );
	FLA_Obj_free( &r );
	FLA_Obj_free( &d );
	FLA_Obj_free( &e );
	FLA_Obj_free( &G );

	return r_val;
}

