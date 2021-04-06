/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 4
#define NUM_MATRIX_ARGS  3
#define FIRST_VARIANT    1
#define LAST_VARIANT     18

// Static variables.
static char* op_str                   = "Triangular Sylvester equation solve";
static char* flash_front_str          = "FLASH_Sylv";
static char* fla_front_str            = "FLA_Sylv";
//static char* fla_unb_var_str          = "unb_var";
static char* fla_opt_var_str          = "opt_var";
static char* fla_blk_var_str          = "blk_var";
static char* pc_str[NUM_PARAM_COMBOS] = { "nn", "nh", "hn", "hh" };
static test_thresh_t thresh           = { 1e-02, 1e-03,   // warn, pass for s
                                          1e-11, 1e-12,   // warn, pass for d
                                          1e-02, 1e-03,   // warn, pass for c
                                          1e-11, 1e-12 }; // warn, pass for z

static fla_sylv_t*      sylv_cntl_opt;
//static fla_sylv_t*      sylv_cntl_unb;
static fla_sylv_t*      sylv_cntl_opt1;
static fla_sylv_t*      sylv_cntl_blk;
static fla_blocksize_t* sylv_cntl_bsize;

// Local prototypes.
void libfla_test_sylv_experiment( test_params_t params,
                                  unsigned int  var,
                                  char*         sc_str,
                                  FLA_Datatype  datatype,
                                  uinteger  p_cur,
                                  unsigned int  pci,
                                  unsigned int  n_repeats,
                                  signed int    impl,
                                  double*       perf,
                                  double*       residual );
void libfla_test_sylv_impl( int         impl,
                            FLA_Trans   transa,
                            FLA_Trans   transb,
                            FLA_Obj     isgn,
                            FLA_Obj     A,
                            FLA_Obj     B,
                            FLA_Obj     C,
                            FLA_Obj     scale );
void libfla_test_sylv_cntl_create( unsigned int var,
                                   dim_t        b_alg_flat );
void libfla_test_sylv_cntl_free( void );


void libfla_test_sylv( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_sylv_experiment );
	}

	if ( op.fla_front == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, NULL,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_FRONT_END,
		                       params, thresh, libfla_test_sylv_experiment );
	}

/*
	if ( op.fla_unb_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_unb_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_UNB_VAR,
		                       params, thresh, libfla_test_sylv_experiment );
	}
*/

	if ( op.fla_opt_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_opt_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_OPT_VAR,
		                       params, thresh, libfla_test_sylv_experiment );
	}

	if ( op.fla_blk_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_blk_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_BLK_VAR,
		                       params, thresh, libfla_test_sylv_experiment );
	}
}



void libfla_test_sylv_experiment( test_params_t params,
                                  unsigned int  var,
                                  char*         sc_str,
                                  FLA_Datatype  datatype,
                                  uinteger  p_cur,
                                  unsigned int  pci,
                                  unsigned int  n_repeats,
                                  signed int    impl,
                                  double*       perf,
                                  double*       residual )
{
	dim_t        b_flash    = params.b_flash;
	dim_t        b_alg_flat = params.b_alg_flat;
	double       time_min   = 1e9;
	double       time;
	unsigned int i;
	uinteger m;
	integer   m_input    = -1;
	uinteger n;
	integer   n_input    = -1;
	integer   sign       = 1;
	FLA_Trans    transa;
	FLA_Trans    transb;
	FLA_Obj      A, B, C, X, isgn, scale, norm;
	FLA_Obj      AX, XB;
	FLA_Obj      C_save;
	FLA_Obj      A_test, B_test, C_test;

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur / abs(m_input);
	else               m = p_cur;
	if ( n_input < 0 ) n = p_cur / abs(n_input);
	else               n = p_cur;

	// Translate parameter characters to libflame constants.
	FLA_Param_map_char_to_flame_trans( &pc_str[pci][0], &transa );
	FLA_Param_map_char_to_flame_trans( &pc_str[pci][1], &transb );

	// Determine the sign.
	if ( 0 < sign ) isgn = FLA_ONE;
	else            isgn = FLA_MINUS_ONE;

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, m, &A );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], n, n, &B );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], m, n, &C );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], m, n, &X );

	// Create the scale object and a temporary norm object.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &scale );
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Initialize the test matrices.
	FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
	FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, B );
	FLA_Random_matrix( C );

	// Tweak A and B.
	FLA_Norm1( A, norm );
	FLA_Shift_diag( FLA_NO_CONJUGATE, norm, A );

	FLA_Norm1( B, norm );
	if ( FLA_Obj_is( isgn, FLA_MINUS_ONE ) ) FLA_Negate( norm );
	FLA_Shift_diag( FLA_NO_CONJUGATE, norm, B );

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, C, &C_save );

	// Create temporary matrices to hold intermediate products, used when
	// checking the solution.
	FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &AX );
	FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &XB );

	// Use hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_create_hier_copy_of_flat( A, 1, &b_flash, &A_test );
		FLASH_Obj_create_hier_copy_of_flat( B, 1, &b_flash, &B_test );
		FLASH_Obj_create_hier_copy_of_flat( C, 1, &b_flash, &C_test );
	}
	else
	{
		A_test = A;
		B_test = B;
		C_test = C;
	}

	// Create a control tree for the individual variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_sylv_cntl_create( var, b_alg_flat );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		if ( impl == FLA_TEST_HIER_FRONT_END )
			FLASH_Obj_hierarchify( C_save, C_test );
		else
			FLA_Copy_external( C_save, C_test );
		
		time = FLA_Clock();

		libfla_test_sylv_impl( impl, transa, transb, isgn, A_test, B_test, C_test, scale );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Copy the solution to flat matrix X.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_flatten( C_test, X );
	}
	else
    {
		FLA_Copy_external( C_test, X );
	}

	// Free the hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_free( &A_test );
		FLASH_Obj_free( &B_test );
		FLASH_Obj_free( &C_test );
	}

	// Free the control trees if we're testing the variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_sylv_cntl_free();

	// Compute the performance of the best experiment repeat.
	*perf = ( m * m * n + n * n * m ) / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

	// Compute the maximum element-wise diff instead of a residual.
	FLA_Trmmsx_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, transa, FLA_NONUNIT_DIAG,
	                     FLA_ONE, A, X, FLA_ZERO, AX );
	FLA_Trmmsx_external( FLA_RIGHT, FLA_UPPER_TRIANGULAR, transb, FLA_NONUNIT_DIAG,
	                     FLA_ONE, B, X, FLA_ZERO, XB );
	FLA_Axpy_external( isgn, XB, AX );
	*residual = FLA_Max_elemwise_diff( AX, C_save );

	// Free the supporting flat objects.
	FLA_Obj_free( &XB );
	FLA_Obj_free( &AX );
	FLA_Obj_free( &C_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
	FLA_Obj_free( &B );
	FLA_Obj_free( &C );
	FLA_Obj_free( &X );
	FLA_Obj_free( &norm );
	FLA_Obj_free( &scale );
}



extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_gemm_t* fla_gemm_cntl_blas;

void libfla_test_sylv_cntl_create( unsigned int var,
                                   dim_t        b_alg_flat )
{
	//int var_unb  = FLA_UNB_VAR_OFFSET + var;
	// Fix all optimized unblocked variants to var1 since internal back-end
	// doesn't handle other variants yet.
	//int var_opt  = FLA_OPT_VAR_OFFSET + var;
	int var_opt  = FLA_OPT_VAR_OFFSET + 1;
	int var1_opt = FLA_OPT_VAR_OFFSET + 1;
	int var_blk  = FLA_BLK_VAR_OFFSET + var;

	sylv_cntl_bsize = FLA_Blocksize_create( b_alg_flat, b_alg_flat, b_alg_flat, b_alg_flat );

/*
	sylv_cntl_unb   = FLA_Cntl_sylv_obj_create( FLA_FLAT,
	                                            var_unb,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL );
*/

	sylv_cntl_opt   = FLA_Cntl_sylv_obj_create( FLA_FLAT,
	                                            var_opt,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL );

	sylv_cntl_opt1  = FLA_Cntl_sylv_obj_create( FLA_FLAT,
	                                            var1_opt,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL );

	sylv_cntl_blk   = FLA_Cntl_sylv_obj_create( FLA_FLAT,
	                                            var_blk,
	                                            sylv_cntl_bsize,
	                                            sylv_cntl_opt1,
	                                            sylv_cntl_opt1,
	                                            sylv_cntl_opt1,
	                                            fla_gemm_cntl_blas,
	                                            fla_gemm_cntl_blas,
	                                            fla_gemm_cntl_blas,
	                                            fla_gemm_cntl_blas,
	                                            fla_gemm_cntl_blas,
	                                            fla_gemm_cntl_blas,
	                                            fla_gemm_cntl_blas,
	                                            fla_gemm_cntl_blas );
}



void libfla_test_sylv_cntl_free( void )
{
	FLA_Blocksize_free( sylv_cntl_bsize );

	//FLA_Cntl_obj_free( sylv_cntl_unb );
	FLA_Cntl_obj_free( sylv_cntl_opt );
	FLA_Cntl_obj_free( sylv_cntl_opt1 );
	FLA_Cntl_obj_free( sylv_cntl_blk );
}



void libfla_test_sylv_impl( int       impl,
                            FLA_Trans transa,
                            FLA_Trans transb,
                            FLA_Obj   isgn,
                            FLA_Obj   A,
                            FLA_Obj   B,
                            FLA_Obj   C,
                            FLA_Obj   scale )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_Sylv( transa, transb, isgn, A, B, C, scale );
		break;

		case FLA_TEST_FLAT_FRONT_END:
		FLA_Sylv( transa, transb, isgn, A, B, C, scale );
		break;

/*
		case FLA_TEST_FLAT_UNB_VAR:
		FLA_Sylv_internal( transa, transb, isgn, A, B, C, scale, sylv_cntl_unb );
		break;
*/

		case FLA_TEST_FLAT_OPT_VAR:
		FLA_Sylv_internal( transa, transb, isgn, A, B, C, scale, sylv_cntl_opt );
		break;

		case FLA_TEST_FLAT_BLK_VAR:
		FLA_Sylv_internal( transa, transb, isgn, A, B, C, scale, sylv_cntl_blk );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

