
#include "FLAME.h"

FLA_Error FLA_Svd_uv_unb_var2( dim_t n_iter_max, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V, dim_t k_accum, dim_t b_alg )
{
	FLA_Error    r_val = FLA_SUCCESS;
	FLA_Datatype dt;
	FLA_Datatype dt_real;
	FLA_Datatype dt_comp;
	FLA_Obj      scale, T, S, rL, rR, d, e, G, H, RG, RH, W;
	dim_t        m_A, n_A;
	dim_t        min_m_n;
	dim_t        n_GH;
	double       crossover_ratio = 17.0 / 9.0;

	n_GH    = k_accum;

	m_A     = FLA_Obj_length( A );
	n_A     = FLA_Obj_width( A );
	min_m_n = FLA_Obj_min_dim( A );
	dt      = FLA_Obj_datatype( A );
	dt_real = FLA_Obj_datatype_proj_to_real( A );
	dt_comp = FLA_Obj_datatype_proj_to_complex( A );

	// If the matrix is a scalar, then the SVD is easy.
	if ( min_m_n == 1 )
	{
		FLA_Copy( A, s );
		FLA_Set_to_identity( U );
		FLA_Set_to_identity( V );

		return FLA_SUCCESS;
	}

	// Create matrices to hold block Householder transformations.
	FLA_Bidiag_UT_create_T( A, &T, &S );

	// Create vectors to hold the realifying scalars.
	FLA_Obj_create( dt,      min_m_n,      1, 0, 0, &rL );
	FLA_Obj_create( dt,      min_m_n,      1, 0, 0, &rR );

	// Create vectors to hold the diagonal and sub-diagonal.
	FLA_Obj_create( dt_real, min_m_n,      1, 0, 0, &d );
	FLA_Obj_create( dt_real, min_m_n-1,    1, 0, 0, &e );

	// Create matrices to hold the left and right Givens scalars.
	FLA_Obj_create( dt_comp, min_m_n-1, n_GH, 0, 0, &G );
	FLA_Obj_create( dt_comp, min_m_n-1, n_GH, 0, 0, &H );

	// Create matrices to hold the left and right Givens matrices.
	FLA_Obj_create( dt_real, min_m_n, min_m_n, 0, 0, &RG );
	FLA_Obj_create( dt_real, min_m_n, min_m_n, 0, 0, &RH );
	FLA_Obj_create( dt,      m_A,     n_A,     0, 0, &W );

	// Create a real scaling factor.
	FLA_Obj_create( dt_real, 1, 1, 0, 0, &scale );

	// Compute a scaling factor; If none is needed, sigma will be set to one.
	FLA_Svd_compute_scaling( A, scale );

	// Scale the matrix if scale is non-unit.
	if ( !FLA_Obj_equals( scale, FLA_ONE ) )
		FLA_Scal( scale, A );

	if ( m_A >= n_A )
	{
		if ( m_A < crossover_ratio * n_A )
		{
			// Reduce the matrix to bidiagonal form.
			// Apply scalars to rotate elements on the sub-diagonal to the real domain.
			// Extract the diagonal and sub-diagonal from A.
			FLA_Bidiag_UT( A, T, S );
			FLA_Bidiag_UT_realify( A, rL, rR );
			FLA_Bidiag_UT_extract_real_diagonals( A, d, e );

			// Form U and V.
			FLA_Bidiag_UT_form_U( A, T, U );
			FLA_Bidiag_UT_form_V( A, S, V );

			// Apply the realifying scalars in rL and rR to U and V, respectively.
			{
				FLA_Obj UL, UR;
				FLA_Obj VL, VR;

				FLA_Part_1x2( U,   &UL, &UR,   min_m_n, FLA_LEFT );
				FLA_Part_1x2( V,   &VL, &VR,   min_m_n, FLA_LEFT );

				FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE,    rL, UL );
				FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, rR, VL );
			}

			// Perform a singular value decomposition on the bidiagonal matrix.
			r_val = FLA_Bsvd_v_opt_var2( n_iter_max, d, e, G, H, RG, RH, W, U, V, b_alg );
		}
		else // if ( crossover_ratio * n_A <= m_A )
		{
			FLA_Obj TQ, R;
			FLA_Obj AT,
			        AB;
			FLA_Obj UL, UR;

			// Perform a QR factorization on A and form Q in U.
			FLA_QR_UT_create_T( A, &TQ );
			FLA_QR_UT( A, TQ );
			FLA_QR_UT_form_Q( A, TQ, U );
			FLA_Obj_free( &TQ );

			// Set the lower triangle of R to zero and then copy the upper
			// triangle of A to R.
			FLA_Part_2x1( A,   &AT,
			                   &AB,   n_A, FLA_TOP );
			FLA_Obj_create( dt, n_A, n_A, 0, 0, &R );
			FLA_Setr( FLA_LOWER_TRIANGULAR, FLA_ZERO, R );
			FLA_Copyr( FLA_UPPER_TRIANGULAR, AT, R );

			// Reduce the matrix to bidiagonal form.
			// Apply scalars to rotate elements on the superdiagonal to the real domain.
			// Extract the diagonal and superdiagonal from A.
			FLA_Bidiag_UT( R, T, S );
			FLA_Bidiag_UT_realify( R, rL, rR );
			FLA_Bidiag_UT_extract_real_diagonals( R, d, e );

			// Form V from right Householder vectors in upper triangle of R.
			FLA_Bidiag_UT_form_V( R, S, V );

			// Form U in R.
			FLA_Bidiag_UT_form_U( R, T, R );

			// Apply the realifying scalars in rL and rR to U and V, respectively.
			FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE,    rL, R );
			FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, rR, V );

			// Perform a singular value decomposition on the bidiagonal matrix.
			r_val = FLA_Bsvd_v_opt_var2( n_iter_max, d, e, G, H, RG, RH, W, R, V, b_alg );

			// Multiply R into U, storing the result in A and then copying back
			// to U.
			FLA_Part_1x2( U,   &UL, &UR,   n_A, FLA_LEFT );
			FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
			          FLA_ONE, UL, R, FLA_ZERO, A );
			FLA_Copy( A, UL );

			FLA_Obj_free( &R );
		}
	}
	else // if ( m_A < n_A )
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	// Copy the converged eigenvalues to the output vector.
	FLA_Copy( d, s );

	// Sort the singular values and singular vectors in descending order.
	FLA_Sort_svd( FLA_BACKWARD, s, U, V );

	// If the matrix was scaled, rescale the singular values.
	if ( !FLA_Obj_equals( scale, FLA_ONE ) )
		FLA_Inv_scal( scale, s );

	FLA_Obj_free( &scale );
	FLA_Obj_free( &T );
	FLA_Obj_free( &S );
	FLA_Obj_free( &rL );
	FLA_Obj_free( &rR );
	FLA_Obj_free( &d );
	FLA_Obj_free( &e );
	FLA_Obj_free( &G );
	FLA_Obj_free( &H );
	FLA_Obj_free( &RG );
	FLA_Obj_free( &RH );
	FLA_Obj_free( &W );

	return r_val;
}

