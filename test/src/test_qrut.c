/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  2
#define FIRST_VARIANT    1
#define LAST_VARIANT     2

// Static variables.
static char* op_str                   = "QR factorization via UT transform";
static char* flash_front_str          = "FLASH_QR_UT";
static char* fla_front_str            = "FLA_QR_UT";
static char* fla_unb_var_str          = "unb_var";
static char* fla_opt_var_str          = "opt_var";
static char* fla_blk_var_str          = "blk_var";
static char* pc_str[NUM_PARAM_COMBOS] = { "" };
static test_thresh_t thresh           = { 1e-02, 1e-03,   // warn, pass for s
                                          1e-11, 1e-12,   // warn, pass for d
                                          1e-02, 1e-03,   // warn, pass for c
                                          1e-11, 1e-12 }; // warn, pass for z

static fla_apqut_t*     apqut_cntl_blk;
static fla_qrut_t*      qrut_cntl_opt;
static fla_qrut_t*      qrut_cntl_unb;
static fla_qrut_t*      qrut_cntl_blk;
static fla_qrut_t*      qrut_cntl_blk_sub;
static fla_blocksize_t* qrut_cntl_bsize;

// Local prototypes.
void libfla_test_qrut_experiment( test_params_t params,
                                  unsigned int  var,
                                  char*         sc_str,
                                  FLA_Datatype  datatype,
                                  unsigned int  p,
                                  unsigned int  pci,
                                  unsigned int  n_repeats,
                                  signed int    impl,
                                  double*       perf,
                                  double*       residual );
void libfla_test_qrut_impl( int     impl,
                            FLA_Obj A,
                            FLA_Obj T );
void libfla_test_qrut_cntl_create( unsigned int var,
                                   dim_t        b_alg_flat );
void libfla_test_qrut_cntl_free( void );


void libfla_test_qrut( FILE* output_stream, test_params_t params, test_op_t op )
{
	libfla_test_output_info( "--- %s ---\n", op_str );
	libfla_test_output_info( "\n" );

	if ( op.flash_front == ENABLE )
	{
		//libfla_test_output_info( "%s() front-end...\n", flash_front_str );
		//libfla_test_output_info( "\n" );
		libfla_test_op_driver( flash_front_str, NULL,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_HIER_FRONT_END,
		                       params, thresh, libfla_test_qrut_experiment );
	}

	if ( op.fla_front == ENABLE )
	{
		//libfla_test_output_info( "%s() front-end...\n", fla_front_str );
		//libfla_test_output_info( "\n" );
		libfla_test_op_driver( fla_front_str, NULL,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_FRONT_END,
		                       params, thresh, libfla_test_qrut_experiment );
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
		                       params, thresh, libfla_test_qrut_experiment );
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
		                       params, thresh, libfla_test_qrut_experiment );
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
		                       params, thresh, libfla_test_qrut_experiment );
	}

}



void libfla_test_qrut_experiment( test_params_t params,
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
	dim_t        b_flash    = params.b_flash;
	dim_t        b_alg_flat = params.b_alg_flat;
	double       time_min   = 1e9;
	double       time;
	unsigned int i;
	unsigned int m, n;
	unsigned int min_m_n;
	signed int   m_input    = -2;
	signed int   n_input    = -1;
	FLA_Obj      A, T, x, b, y, norm;
	FLA_Obj      A_save;
	FLA_Obj      A_test, T_test, x_test, b_test;

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur * abs(m_input);
	else               m = p_cur;
	if ( n_input < 0 ) n = p_cur * abs(n_input);
	else               n = p_cur;

	// Compute the minimum dimension.
	min_m_n = min( m, n );

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, n, &A );

	if ( impl == FLA_TEST_FLAT_FRONT_END ||
	     ( impl == FLA_TEST_FLAT_BLK_VAR && var == 1 ) )
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], b_alg_flat, min_m_n, &T );
	else if ( var == 2 )
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], min_m_n, min_m_n, &T );
	else
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], 1, min_m_n, &T );

	// Initialize the test matrices.
	FLA_Random_matrix( A );

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_save );

	// Create vectors to form a linear system.
	FLA_Obj_create( datatype, n, 1, 0, 0, &x );
	FLA_Obj_create( datatype, m, 1, 0, 0, &b );
	FLA_Obj_create( datatype, n, 1, 0, 0, &y );

	// Create a real scalar object to hold the norm of A.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Create a random right-hand side vector.
	FLA_Random_matrix( b );

	// Use hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_QR_UT_create_hier_matrices( A, 1, &b_flash, &A_test, &T_test );
		FLASH_Obj_create_hier_copy_of_flat( b, 1, &b_flash, &b_test );
		FLASH_Obj_create_hier_copy_of_flat( x, 1, &b_flash, &x_test );
	}
	else
	{
		A_test = A;
		T_test = T;
	}

	// Create a control tree for the individual variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_qrut_cntl_create( var, b_alg_flat );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		if ( impl == FLA_TEST_HIER_FRONT_END )
			FLASH_Obj_hierarchify( A_save, A_test );
		else
			FLA_Copy_external( A_save, A_test );
		
		time = FLA_Clock();

		libfla_test_qrut_impl( impl, A_test, T_test );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Perform a linear solve with the result.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_QR_UT_solve( A_test, T_test, b_test, x_test );
		FLASH_Obj_flatten( x_test, x );
	}
	else
    {
		FLA_QR_UT_solve( A_test, T_test, b, x );
	}

	// Free the hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_free( &A_test );
		FLASH_Obj_free( &T_test );
		FLASH_Obj_free( &b_test );
		FLASH_Obj_free( &x_test );
	}

	// Free the control trees if we're testing the variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_qrut_cntl_free();

	// Compute the performance of the best experiment repeat.
	*perf = (         2.0   * m * n * n - 
	          ( 2.0 / 3.0 ) * n * n * n ) / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

	// Compute the residual.
	FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A_save, x, FLA_MINUS_ONE, b );
	FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A_save, b, FLA_ZERO, y );
	FLA_Nrm2_external( y, norm );
	FLA_Obj_extract_real_scalar( norm, residual );

	// Free the supporting flat objects.
	FLA_Obj_free( &x );
	FLA_Obj_free( &b );
	FLA_Obj_free( &y );
	FLA_Obj_free( &norm );
	FLA_Obj_free( &A_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
	FLA_Obj_free( &T );
}



extern __thread fla_axpyt_t* fla_axpyt_cntl_blas;
extern __thread fla_copyt_t* fla_copyt_cntl_blas;
extern __thread fla_gemm_t*  fla_gemm_cntl_blas;
extern __thread fla_trmm_t*  fla_trmm_cntl_blas;
extern __thread fla_trsm_t*  fla_trsm_cntl_blas;

void libfla_test_qrut_cntl_create( unsigned int var,
                                   dim_t        b_alg_flat )
{
	int var_unb  = FLA_UNB_VAR_OFFSET + var;
	int var_opt  = FLA_OPT_VAR_OFFSET + var;
	int var_blk  = FLA_BLK_VAR_OFFSET + var;
	int var_opt2 = FLA_OPT_VAR_OFFSET + 2;

	qrut_cntl_bsize = FLA_Blocksize_create( b_alg_flat, b_alg_flat, b_alg_flat, b_alg_flat );

	apqut_cntl_blk  = FLA_Cntl_apqut_obj_create( FLA_FLAT,
	                                             FLA_BLOCKED_VARIANT1,
	                                             qrut_cntl_bsize,
	                                             NULL,
	                                             fla_trmm_cntl_blas,
	                                             fla_trmm_cntl_blas,
	                                             fla_gemm_cntl_blas,
	                                             fla_gemm_cntl_blas,
	                                             fla_trsm_cntl_blas,
	                                             fla_copyt_cntl_blas,
	                                             fla_axpyt_cntl_blas );

	qrut_cntl_unb   = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                                var_unb,
                                                NULL,
                                                NULL,
                                                NULL );

	qrut_cntl_opt   = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                                var_opt,
                                                NULL,
                                                NULL,
                                                NULL );

	// Use unblocked variant 2 for blocked variant 1 subproblem.
	qrut_cntl_blk_sub = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                                  var_opt2,
                                                  NULL,
                                                  NULL,
                                                  NULL );

	qrut_cntl_blk   = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                                var_blk,
                                                qrut_cntl_bsize,
                                                qrut_cntl_blk_sub,
                                                apqut_cntl_blk );
}



void libfla_test_qrut_cntl_free( void )
{
	FLA_Blocksize_free( qrut_cntl_bsize );

	FLA_Cntl_obj_free( apqut_cntl_blk );
	FLA_Cntl_obj_free( qrut_cntl_unb );
	FLA_Cntl_obj_free( qrut_cntl_opt );
	FLA_Cntl_obj_free( qrut_cntl_blk );
	FLA_Cntl_obj_free( qrut_cntl_blk_sub );
}



void libfla_test_qrut_impl( int     impl,
                            FLA_Obj A,
                            FLA_Obj T )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_QR_UT( A, T );
		break;

		case FLA_TEST_FLAT_FRONT_END:
		FLA_QR_UT( A, T );
		break;

		case FLA_TEST_FLAT_UNB_VAR:
		FLA_QR_UT_internal( A, T, qrut_cntl_unb );
		break;

		case FLA_TEST_FLAT_OPT_VAR:
		FLA_QR_UT_internal( A, T, qrut_cntl_opt );
		break;

		case FLA_TEST_FLAT_BLK_VAR:
		FLA_QR_UT_internal( A, T, qrut_cntl_blk );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

