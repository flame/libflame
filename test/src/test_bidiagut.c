/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 2
#define NUM_MATRIX_ARGS  3
#define FIRST_VARIANT    1
#define LAST_VARIANT     5

// Static variables.
static char* op_str                   = "Reduction to bidiagonal form via UT transform";
static char* fla_front_str            = "FLA_Bidiag_UT";
static char* fla_unb_var_str          = "unb_var";
static char* fla_opt_var_str          = "opt_var";
static char* fla_blk_var_str          = "blk_var";
static char* pc_str[NUM_PARAM_COMBOS] = { "l", "u" };
static test_thresh_t thresh           = { 1e-03, 1e-04,   // warn, pass for s
                                          1e-12, 1e-13,   // warn, pass for d
                                          1e-03, 1e-04,   // warn, pass for c
                                          1e-12, 1e-13 }; // warn, pass for z

static fla_bidiagut_t*  bidiagut_cntl_opt;
static fla_bidiagut_t*  bidiagut_cntl_unb;
static fla_bidiagut_t*  bidiagut_cntl_blk;
static fla_blocksize_t* bidiagut_cntl_bsize;

// Local prototypes.
void libfla_test_bidiagut_experiment( test_params_t params,
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
void libfla_test_bidiagut_impl( int      impl,
                                FLA_Obj  A,
                                FLA_Obj  TU,
                                FLA_Obj  TV );
void libfla_test_bidiagut_cntl_create( unsigned int var,
                                       dim_t        b_alg_flat );
void libfla_test_bidiagut_cntl_free( void );


void libfla_test_bidiagut( FILE* output_stream, test_params_t params, test_op_t op )
{
	libfla_test_output_info( "--- %s ---\n", op_str );
	libfla_test_output_info( "\n" );

	if ( op.fla_front == ENABLE )
	{
		//libfla_test_output_info( "%s() front-end...\n", fla_front_str );
		//libfla_test_output_info( "\n" );
		libfla_test_op_driver( fla_front_str, NULL,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_FRONT_END,
		                       params, thresh, libfla_test_bidiagut_experiment );
	}

	if ( op.fla_unb_vars == ENABLE )
	{
		//libfla_test_output_info( "%s() unblocked variants...\n", fla_front_str );
		//libfla_test_output_info( "\n" );
		libfla_test_op_driver( fla_front_str, fla_unb_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_UNB_VAR,
		                       params, thresh, libfla_test_bidiagut_experiment );
	}

	if ( op.fla_opt_vars == ENABLE )
	{
		//libfla_test_output_info( "%s() optimized unblocked variants...\n", fla_front_str );
		//libfla_test_output_info( "\n" );
		libfla_test_op_driver( fla_front_str, fla_opt_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_OPT_VAR,
		                       params, thresh, libfla_test_bidiagut_experiment );
	}

	if ( op.fla_blk_vars == ENABLE )
	{
		//libfla_test_output_info( "%s() blocked variants...\n", fla_front_str );
		//libfla_test_output_info( "\n" );
		libfla_test_op_driver( fla_front_str, fla_blk_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_BLK_VAR,
		                       params, thresh, libfla_test_bidiagut_experiment );
	}

}



void libfla_test_bidiagut_experiment( test_params_t params,
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
	dim_t        b_alg_flat = params.b_alg_flat;
	double       time_min   = 1e9;
	double       time;
	unsigned int i;
	uinteger m, n, min_m_n;
	integer   m_input    = -1;
	integer   n_input    = -2;
	FLA_Obj      A, TU, TV, WU, WV, QUh, QV, AQV, QUhAQV, norm;
	FLA_Obj      AL, AR;
	FLA_Obj      QVL, QVR;
	FLA_Obj      A_save;

	// Return early for 'l' parameter case since it does not yet exist.
	if ( pc_str[pci][0] == 'l' )
	{
		*perf     = 0.0;
		*residual = 0.0;
		return;
	}

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur / abs(m_input);
	else               m = p_cur;
	if ( n_input < 0 ) n = p_cur / abs(n_input);
	else               n = p_cur;

/// Uncomment when lower bidiagonal case is implemented.
/// This code swaps the m and n dimension because we want to use "tall"
/// rectangles for the upper bidiagonal case and "fat" rectangles for the
/// lower bidiagonal case. The ratios are hard-coded as m = -1 and n = -2
/// which only works for one or the other, depending on whether we divide
/// or multiply when computing m and n above.
/*
	// Swap the larger and smaller dimensions for the lower bidiagonal case.
	if ( pc_str[pci][0] == 'u' )
	{
		// Do nothing.
		;
	}
	else // if ( pc_str[pci][0] == 'l' )
	{
		uinteger m_temp = n;
		uinteger n_temp = m;
		m = m_temp;
		n = n_temp;
	}
*/

	// Compute the minimum dimension of A.
	min_m_n = min( m, n );

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, n, &A );

	if ( impl == FLA_TEST_FLAT_FRONT_END ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
	{
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], b_alg_flat, min_m_n, &TU );
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], b_alg_flat, min_m_n, &TV );
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], b_alg_flat,       m, &WU );
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], b_alg_flat,       n, &WV );
	}
	else
	{
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], min_m_n, min_m_n, &TU );
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], min_m_n, min_m_n, &TV );
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], min_m_n,       m, &WU );
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], min_m_n,       n, &WV );
	}

	// Initialize the test matrices.
	FLA_Random_matrix( A );

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_save );

	// Create auxiliary matrices to be used when checking the result.
	FLA_Obj_create( datatype, m, m, 0, 0, &QUh );
	FLA_Obj_create( datatype, n, n, 0, 0, &QV );
	FLA_Obj_create( datatype, m, n, 0, 0, &AQV );
	FLA_Obj_create( datatype, m, n, 0, 0, &QUhAQV );

	// Create a real scalar object to hold the norm of A.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Create a control tree for the individual variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_bidiagut_cntl_create( var, b_alg_flat );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		FLA_Copy_external( A_save, A );
		
		time = FLA_Clock();

		libfla_test_bidiagut_impl( impl, A, TU, TV );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Free the control trees if we're testing the variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_bidiagut_cntl_free();

	// Compute the performance of the best experiment repeat.
	*t = time_min;
  *perf = ( 4.0 * n * n * ( m - n / 3.0 ) ) / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

	// Check the result by computing R - Q' A_orig Q.
	if ( m >= n )
	{
		FLA_Set_to_identity( QUh );
		FLA_Set_to_identity( QV );
		FLA_Part_1x2( QV,   &QVL, &QVR,    1, FLA_LEFT );
		FLA_Part_1x2( A,    &AL,  &AR,     1, FLA_LEFT );
		FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
		                A, TU, WU, QUh );
		FLA_Apply_Q_UT( FLA_RIGHT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_ROWWISE,
		                AR, TV, WV, QVR );
		FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
		          FLA_ONE, A_save, QV, FLA_ZERO, AQV );
		FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
		          FLA_ONE, QUh, AQV, FLA_ZERO, QUhAQV );
		FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
		FLA_Triangularize( FLA_LOWER_TRIANGULAR, FLA_NONUNIT_DIAG, AR );
		*residual = FLA_Max_elemwise_diff( A, QUhAQV );
	}
	else // if ( m < n )
	{
		;
	}

	// Free the supporting flat objects.
	FLA_Obj_free( &WU );
	FLA_Obj_free( &WV );
	FLA_Obj_free( &QUh );
	FLA_Obj_free( &QV );
	FLA_Obj_free( &AQV );
	FLA_Obj_free( &QUhAQV );
	FLA_Obj_free( &norm );
	FLA_Obj_free( &A_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
	FLA_Obj_free( &TU );
	FLA_Obj_free( &TV );
}



void libfla_test_bidiagut_cntl_create( unsigned int var,
                                       dim_t        b_alg_flat )
{
	int var_unb  = FLA_UNB_VAR_OFFSET + var;
	int var_opt  = FLA_OPT_VAR_OFFSET + var;
	int var_blk  = FLA_BLK_VAR_OFFSET + var;

	bidiagut_cntl_bsize = FLA_Blocksize_create( b_alg_flat,
	                                            b_alg_flat,
	                                            b_alg_flat,
	                                            b_alg_flat );

	bidiagut_cntl_unb   = FLA_Cntl_bidiagut_obj_create( FLA_FLAT,
	                                                    var_unb,
	                                                    NULL );

	bidiagut_cntl_opt   = FLA_Cntl_bidiagut_obj_create( FLA_FLAT,
	                                                    var_opt,
	                                                    NULL );

	bidiagut_cntl_blk   = FLA_Cntl_bidiagut_obj_create( FLA_FLAT,
	                                                    var_blk,
	                                                    NULL );
}



void libfla_test_bidiagut_cntl_free( void )
{
	FLA_Blocksize_free( bidiagut_cntl_bsize );

	FLA_Cntl_obj_free( bidiagut_cntl_unb );
	FLA_Cntl_obj_free( bidiagut_cntl_opt );
	FLA_Cntl_obj_free( bidiagut_cntl_blk );
}



void libfla_test_bidiagut_impl( int      impl,
                                FLA_Obj  A,
                                FLA_Obj  TU,
                                FLA_Obj  TV )
{
	switch ( impl )
	{
		case FLA_TEST_FLAT_FRONT_END:
		FLA_Bidiag_UT( A, TU, TV );
		break;

		case FLA_TEST_FLAT_UNB_VAR:
		FLA_Bidiag_UT_internal( A, TU, TV, bidiagut_cntl_unb );
		break;

		case FLA_TEST_FLAT_OPT_VAR:
		FLA_Bidiag_UT_internal( A, TU, TV, bidiagut_cntl_opt );
		break;

		case FLA_TEST_FLAT_BLK_VAR:
		FLA_Bidiag_UT_internal( A, TU, TV, bidiagut_cntl_blk );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

