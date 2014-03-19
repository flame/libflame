/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 2
#define NUM_MATRIX_ARGS  2
#define FIRST_VARIANT    1
#define LAST_VARIANT     3

// Static variables.
static char* op_str                   = "Reduction to tridiagonal form via UT transform";
static char* fla_front_str            = "FLA_Tridiag_UT";
static char* fla_unb_var_str          = "unb_var";
static char* fla_opt_var_str          = "opt_var";
static char* fla_blk_var_str          = "blk_var";
static char* pc_str[NUM_PARAM_COMBOS] = { "l", "u" };
static test_thresh_t thresh           = { 1e-03, 1e-04,   // warn, pass for s
                                          1e-12, 1e-13,   // warn, pass for d
                                          1e-03, 1e-04,   // warn, pass for c
                                          1e-12, 1e-13 }; // warn, pass for z

static fla_tridiagut_t* tridiagut_cntl_opt;
static fla_tridiagut_t* tridiagut_cntl_unb;
static fla_tridiagut_t* tridiagut_cntl_blk;
static fla_blocksize_t* tridiagut_cntl_bsize;

// Local prototypes.
void libfla_test_tridiagut_experiment( test_params_t params,
                                       unsigned int  var,
                                       char*         sc_str,
                                       FLA_Datatype  datatype,
                                       unsigned int  p,
                                       unsigned int  pci,
                                       unsigned int  n_repeats,
                                       signed int    impl,
                                       double*       perf,
                                       double*       residual );
void libfla_test_tridiagut_impl( int      impl,
                                 FLA_Uplo uplo,
                                 FLA_Obj  A,
                                 FLA_Obj  T );
void libfla_test_tridiagut_cntl_create( unsigned int var,
                                        dim_t        b_alg_flat );
void libfla_test_tridiagut_cntl_free( void );


void libfla_test_tridiagut( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_tridiagut_experiment );
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
		                       params, thresh, libfla_test_tridiagut_experiment );
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
		                       params, thresh, libfla_test_tridiagut_experiment );
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
		                       params, thresh, libfla_test_tridiagut_experiment );
	}

}



void libfla_test_tridiagut_experiment( test_params_t params,
                                       unsigned int  var,
                                       char*         sc_str,
                                       FLA_Datatype  datatype,
                                       unsigned int  p_cur,
                                       unsigned int  pci,
                                       unsigned int  n_repeats,
                                       signed int    impl,
                                       double*       perf,
                                       double*       residual )
{
	dim_t        b_alg_flat = params.b_alg_flat;
	double       time_min   = 1e9;
	double       time;
	unsigned int i;
	unsigned int m;
	signed int   m_input    = -1;
	FLA_Uplo     uplo;
	FLA_Obj      A, T, W, Qh, AQ, QhAQ, norm;
	FLA_Obj      AT, AB;
	FLA_Obj      QhT, QhB;
	FLA_Obj      A_save;

	// Return early for 'u' parameter case since it does not yet exist.
	if ( pc_str[pci][0] == 'u' )
	{
		*perf     = 0.0;
		*residual = 0.0;
		return;
	}

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur * abs(m_input);
	else               m = p_cur;

	// Translate parameter characters to libflame constants.
	FLA_Param_map_char_to_flame_uplo( &pc_str[pci][0], &uplo );

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, m, &A );

	if ( impl == FLA_TEST_FLAT_FRONT_END ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
	{
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], b_alg_flat, m, &T );
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], b_alg_flat, m, &W );
	}
	else
	{
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], m, m, &T );
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], m, m, &W );
	}

	// Initialize the test matrices.
	FLA_Random_spd_matrix( uplo, A );
	FLA_Hermitianize( uplo, A );

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_save );

	// Create auxiliary matrices to be used when checking the result.
	FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Qh );
	FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &AQ );
	FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &QhAQ );

	// Create a real scalar object to hold the norm of A.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Create a control tree for the individual variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_tridiagut_cntl_create( var, b_alg_flat );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		FLA_Copy_external( A_save, A );
		
		time = FLA_Clock();

		libfla_test_tridiagut_impl( impl, uplo, A, T );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Free the control trees if we're testing the variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_tridiagut_cntl_free();

	// Compute the performance of the best experiment repeat.
	*perf = ( 4.0 / 3.0 * m * m * m ) / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

	// Check the result by computing R - Q' A_orig Q.
	if ( uplo == FLA_LOWER_TRIANGULAR )
	{
		FLA_Set_to_identity( Qh );
		FLA_Part_2x1( Qh,   &QhT,
		                    &QhB,   1, FLA_TOP );
		FLA_Part_2x1( A,    &AT,
		                    &AB,    1, FLA_TOP );
		FLA_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
		                AB, T, W, QhB );
		FLA_Gemm( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE,
		          FLA_ONE, A_save, Qh, FLA_ZERO, AQ );
		FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
		          FLA_ONE, Qh, AQ, FLA_ZERO, QhAQ );
		FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, AB );
		FLA_Hermitianize( FLA_LOWER_TRIANGULAR, A );
		*residual = FLA_Max_elemwise_diff( A, QhAQ );
	}
	else
	{
		;
	}

	// Free the supporting flat objects.
	FLA_Obj_free( &W );
	FLA_Obj_free( &Qh );
	FLA_Obj_free( &AQ );
	FLA_Obj_free( &QhAQ );
	FLA_Obj_free( &norm );
	FLA_Obj_free( &A_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
	FLA_Obj_free( &T );
}



void libfla_test_tridiagut_cntl_create( unsigned int var,
                                        dim_t        b_alg_flat )
{
	int var_unb  = FLA_UNB_VAR_OFFSET + var;
	int var_opt  = FLA_OPT_VAR_OFFSET + var;
	int var_blk  = FLA_BLK_VAR_OFFSET + var;

	tridiagut_cntl_bsize = FLA_Blocksize_create( b_alg_flat,
	                                             b_alg_flat,
	                                             b_alg_flat,
	                                             b_alg_flat );

	tridiagut_cntl_unb   = FLA_Cntl_tridiagut_obj_create( FLA_FLAT,
	                                                      var_unb,
	                                                      NULL );

	tridiagut_cntl_opt   = FLA_Cntl_tridiagut_obj_create( FLA_FLAT,
	                                                      var_opt,
	                                                      NULL );

	tridiagut_cntl_blk   = FLA_Cntl_tridiagut_obj_create( FLA_FLAT,
	                                                      var_blk,
	                                                      NULL );
}



void libfla_test_tridiagut_cntl_free( void )
{
	FLA_Blocksize_free( tridiagut_cntl_bsize );

	FLA_Cntl_obj_free( tridiagut_cntl_unb );
	FLA_Cntl_obj_free( tridiagut_cntl_opt );
	FLA_Cntl_obj_free( tridiagut_cntl_blk );
}



void libfla_test_tridiagut_impl( int      impl,
                                 FLA_Uplo uplo,
                                 FLA_Obj  A,
                                 FLA_Obj  T )
{
	switch ( impl )
	{
		case FLA_TEST_FLAT_FRONT_END:
		FLA_Tridiag_UT( uplo, A, T );
		break;

		case FLA_TEST_FLAT_UNB_VAR:
		FLA_Tridiag_UT_internal( uplo, A, T, tridiagut_cntl_unb );
		break;

		case FLA_TEST_FLAT_OPT_VAR:
		FLA_Tridiag_UT_internal( uplo, A, T, tridiagut_cntl_opt );
		break;

		case FLA_TEST_FLAT_BLK_VAR:
		FLA_Tridiag_UT_internal( uplo, A, T, tridiagut_cntl_blk );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

