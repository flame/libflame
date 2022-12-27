/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 16
#define NUM_MATRIX_ARGS  3
#define FIRST_VARIANT    1
#define LAST_VARIANT     1

// Static variables.
static char* op_str                   = "Apply Q via UT transform";
static char* flash_front_str          = "FLASH_Apply_Q_UT";
static char* fla_front_str            = "FLA_Apply_Q_UT";
//static char* fla_unb_var_str          = "";
//static char* fla_opt_var_str          = "";
//static char* fla_blk_var_str          = "";
static char* pc_str[NUM_PARAM_COMBOS] = { "lnfc", "lnfr", "lnbc", "lnbr",
                                          "lhfc", "lhfr", "lhbc", "lhbr",
                                          "rnfc", "rnfr", "rnbc", "rnbr",
                                          "rhfc", "rhfr", "rhbc", "rhbr" };
static test_thresh_t thresh           = { 1e-03, 1e-04,   // warn, pass for s
                                          1e-12, 1e-13,   // warn, pass for d
                                          1e-03, 1e-04,   // warn, pass for c
                                          1e-12, 1e-13 }; // warn, pass for z

// Local prototypes.
void libfla_test_apqut_experiment( test_params_t params,
                                   unsigned int  var,
                                   char*         sc_str,
                                   FLA_Datatype  datatype,
                                   uinteger  p,
                                   unsigned int  pci,
                                   unsigned int  n_repeats,
                                   signed int    impl,
                                   double*       perf,
                                   double*       t,
                                   double*       residual );
void libfla_test_apqut_impl( int        impl,
                             FLA_Side   side,
                             FLA_Trans  trans,
                             FLA_Direct direct,
                             FLA_Store  storev,
                             FLA_Obj    A,
                             FLA_Obj    T,
                             FLA_Obj    W,
                             FLA_Obj    B );

void libfla_test_apqut( FILE* output_stream, test_params_t params, test_op_t op )
{
	libfla_test_output_info( "--- %s ---\n", op_str );
	libfla_test_output_info( "\n" );

	if ( op.flash_front == ENABLE )
	{
		libfla_test_op_driver( flash_front_str, NULL,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_HIER_FRONT_END,
		                       params, thresh, libfla_test_apqut_experiment );
	}

	if ( op.fla_front == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, NULL,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_FRONT_END,
		                       params, thresh, libfla_test_apqut_experiment );
	}

}



void libfla_test_apqut_experiment( test_params_t params,
                                   unsigned int  var,
                                   char*         sc_str,
                                   FLA_Datatype  datatype,
                                   uinteger  p_cur,
                                   unsigned int  pci,
                                   unsigned int  n_repeats,
                                   signed int    impl,
                                   double*       perf,
                                   double*       t,
                                   double*       residual )
{
	dim_t        b_flash    = params.b_flash;
	dim_t        b_alg_flat = params.b_alg_flat;
	double       time_min   = 1e9;
	double       time;
	unsigned int i;
	uinteger m, n;
	uinteger min_m_n;
	integer   m_input;
	integer   n_input;
	FLA_Side     side;
	FLA_Trans    trans;
	FLA_Direct   direct;
	FLA_Store    storev;
	FLA_Obj      A, T, W, B, eye, norm;
	FLA_Obj      B_save;
	FLA_Obj      A_test, T_test, W_test, B_test;

	// Translate parameter characters to libflame constants.
	FLA_Param_map_char_to_flame_side( &pc_str[pci][0], &side );
	FLA_Param_map_char_to_flame_trans( &pc_str[pci][1], &trans );
	FLA_Param_map_char_to_flame_direct( &pc_str[pci][2], &direct );
	FLA_Param_map_char_to_flame_storev( &pc_str[pci][3], &storev );

	// We want to make sure the Apply_Q_UT routines work with rectangular
	// matrices. So we use m > n when testing with column-wise storage (via
	// QR factorization) and m < n when testing with row-wise storage (via
	// LQ factorization).
	if ( storev == FLA_COLUMNWISE )
	{
		m_input = -1;
		n_input = -1;
		//m_input = -1;
		//n_input = -1;
	}
	else // if ( storev == FLA_ROWWISE )
	{
		m_input = -1;
		n_input = -1;
		//m_input = -1;
		//n_input = -1;
	}

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur * abs(m_input);
	else               m = p_cur;
	if ( n_input < 0 ) n = p_cur * abs(n_input);
	else               n = p_cur;

	// Compute the minimum dimension.
	min_m_n = fla_min( m, n );

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, n, &A );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], b_alg_flat, min_m_n, &T );
	if ( storev == FLA_COLUMNWISE )
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], m, m, &B );
	else
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], n, n, &B );

	FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, B, &eye );

	FLA_Apply_Q_UT_create_workspace( T, B, &W );

	// Create a real scalar object to hold the norm of A.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Initialize the test matrices.
	FLA_Random_matrix( A );
	FLA_Set_to_identity( B );
	FLA_Set_to_identity( eye );

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, B, &B_save );

	// Use hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		if ( storev == FLA_COLUMNWISE )
			FLASH_QR_UT_create_hier_matrices( A, 1, &b_flash, &A_test, &T_test );
		else // if ( storev == FLA_ROWWISE )
			FLASH_LQ_UT_create_hier_matrices( A, 1, &b_flash, &A_test, &T_test );
		FLASH_Obj_create_hier_copy_of_flat( B, 1, &b_flash, &B_test );
		FLASH_Apply_Q_UT_create_workspace( T_test, B_test, &W_test );
	}
	else // if ( impl == FLA_TEST_FLAT_FRONT_END )
	{
		A_test = A;
		T_test = T;
		W_test = W;
		B_test = B;
	}

	// Compute a Householder factorization.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		if ( storev == FLA_COLUMNWISE ) FLASH_QR_UT( A_test, T_test );
		else                            FLASH_LQ_UT( A_test, T_test );
	}
	else // if ( impl == FLA_TEST_FLAT_FRONT_END )
	{
		if ( storev == FLA_COLUMNWISE ) FLA_QR_UT( A_test, T_test );
		else                            FLA_LQ_UT( A_test, T_test );
	}

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		if ( impl == FLA_TEST_HIER_FRONT_END )
			FLASH_Obj_hierarchify( B_save, B_test );
		else
			FLA_Copy_external( B_save, B_test );

		time = FLA_Clock();

		libfla_test_apqut_impl( impl, side, trans, direct, storev,
		                        A_test, T_test, W_test, B_test );
		
		time = FLA_Clock() - time;
		time_min = fla_min( time_min, time );
	}

	// Multiply by its conjugate-transpose to get what should be (near) identity
	// and then subtract from actual identity to get what should be (near) zero.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_flatten( B_test, B );
		FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
		                   FLA_ONE, B, B, FLA_MINUS_ONE, eye );
	}
	else // if ( impl == FLA_TEST_FLAT_FRONT_END )
	{
		FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
		                   FLA_ONE, B, B, FLA_MINUS_ONE, eye );
	}

	// Free the hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_free( &A_test );
		FLASH_Obj_free( &T_test );
		FLASH_Obj_free( &W_test );
		FLASH_Obj_free( &B_test );
	}

	// Compute the norm of eye, which contains I - Q * Q'.
	FLA_Norm1( eye, norm );
	FLA_Obj_extract_real_scalar( norm, residual );

	// Compute the performance of the best experiment repeat.
	*t = time_min;
  *perf = (  4.0 *       m * min_m_n * n -
	           2.0 * min_m_n * min_m_n * n ) / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

	// Free the supporting flat objects.
	FLA_Obj_free( &B_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
	FLA_Obj_free( &T );
	FLA_Obj_free( &W );
	FLA_Obj_free( &B );
	FLA_Obj_free( &eye );
	FLA_Obj_free( &norm );
}



void libfla_test_apqut_impl( int        impl,
                             FLA_Side   side,
                             FLA_Trans  trans,
                             FLA_Direct direct,
                             FLA_Store  storev,
                             FLA_Obj    A,
                             FLA_Obj    T,
                             FLA_Obj    W,
                             FLA_Obj    B )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_Apply_Q_UT( side, trans, direct, storev,
		                  A, T, W, B );
		break;

		case FLA_TEST_FLAT_FRONT_END:
		FLA_Apply_Q_UT( side, trans, direct, storev,
		                A, T, W, B );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

