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
#define LAST_VARIANT     10

// Static variables.
static char* op_str                   = "Hermitian matrix-matrix multiply";
static char* flash_front_str          = "FLASH_Hemm";
static char* fla_front_str            = "FLA_Hemm";
static char* fla_unb_var_str          = "unb_var";
//static char* fla_opt_var_str          = "opt_var";
static char* fla_blk_var_str          = "blk_var";
static char* fla_unb_ext_str          = "external";
static char* pc_str[NUM_PARAM_COMBOS] = { "ll", "lu",
                                          "rl", "ru" };
static test_thresh_t thresh           = { 1e-02, 1e-03,   // warn, pass for s
                                          1e-11, 1e-12,   // warn, pass for d
                                          1e-02, 1e-03,   // warn, pass for c
                                          1e-11, 1e-12 }; // warn, pass for z

static fla_hemm_t*      hemm_cntl_unb;
static fla_hemm_t*      hemm_cntl_blk;
static fla_blocksize_t* hemm_cntl_bsize;

// Local prototypes.
void libfla_test_hemm_experiment( test_params_t params,
                                  unsigned int  var,
                                  char*         sc_str,
                                  FLA_Datatype  datatype,
                                  unsigned int  p_cur,
                                  unsigned int  pci,
                                  unsigned int  n_repeats,
                                  signed int    impl,
                                  double*       perf,
                                  double*       residual );
void libfla_test_hemm_impl( int         impl,
                            FLA_Side    side,
                            FLA_Uplo    uplo,
                            FLA_Obj     alpha,
                            FLA_Obj     A,
                            FLA_Obj     B,
                            FLA_Obj     beta,
                            FLA_Obj     C );
void libfla_test_hemm_cntl_create( unsigned int var,
                                   dim_t        b_alg_flat );
void libfla_test_hemm_cntl_free( void );


void libfla_test_hemm( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_hemm_experiment );
	}

	if ( op.fla_front == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, NULL,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_FRONT_END,
		                       params, thresh, libfla_test_hemm_experiment );
	}

	if ( op.fla_unb_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_unb_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_UNB_VAR,
		                       params, thresh, libfla_test_hemm_experiment );
	}

	if ( op.fla_blk_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_blk_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_BLK_VAR,
		                       params, thresh, libfla_test_hemm_experiment );
	}

	if ( op.fla_unb_ext  == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_unb_ext_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_UNB_EXT,
		                       params, thresh, libfla_test_hemm_experiment );
	}
}



void libfla_test_hemm_experiment( test_params_t params,
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
	unsigned int m;
	signed int   m_input    = -1;
	unsigned int n;
	signed int   n_input    = -1;
	FLA_Side     side;
	FLA_Uplo     uplo;
	FLA_Obj      A, B, C, x, y, z, w, norm;
	FLA_Obj      alpha, beta;
	FLA_Obj      C_save;
	FLA_Obj      A_test, B_test, C_test;

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur / abs(m_input);
	else               m = p_cur;
	if ( n_input < 0 ) n = p_cur / abs(n_input);
	else               n = p_cur;

	// Translate parameter characters to libflame constants.
	FLA_Param_map_char_to_flame_side( &pc_str[pci][0], &side );
	FLA_Param_map_char_to_flame_uplo( &pc_str[pci][1], &uplo );

	// Create the matrices for the current operation.
	if ( side == FLA_LEFT )
	{
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, m, &A );

		// Create vectors for use in test.
		FLA_Obj_create( datatype, n, 1, 0, 0, &x );
		FLA_Obj_create( datatype, m, 1, 0, 0, &y );
		FLA_Obj_create( datatype, m, 1, 0, 0, &z );
		FLA_Obj_create( datatype, m, 1, 0, 0, &w );
	}
	else
	{
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n, n, &A );

		// Create vectors for use in test.
		FLA_Obj_create( datatype, n, 1, 0, 0, &x );
		FLA_Obj_create( datatype, m, 1, 0, 0, &y );
		FLA_Obj_create( datatype, m, 1, 0, 0, &z );
		FLA_Obj_create( datatype, n, 1, 0, 0, &w );
	}
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], m, n, &B );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], m, n, &C );

	// Create a norm scalar.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Initialize the test matrices.
	FLA_Random_herm_matrix( uplo, A );
	FLA_Random_matrix( B );
	FLA_Random_matrix( C );

	// Initialize the test vectors.
	FLA_Random_matrix( x );
    FLA_Set( FLA_ZERO, y );
    FLA_Set( FLA_ZERO, z );
    FLA_Set( FLA_ZERO, w );

	// Set constants.
	alpha = FLA_TWO;
	beta  = FLA_MINUS_ONE;

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, C, &C_save );

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
	     impl == FLA_TEST_FLAT_BLK_VAR ||
	     impl == FLA_TEST_FLAT_UNB_EXT ||
	     impl == FLA_TEST_FLAT_BLK_EXT )
		libfla_test_hemm_cntl_create( var, b_alg_flat );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		if ( impl == FLA_TEST_HIER_FRONT_END )
			FLASH_Obj_hierarchify( C_save, C_test );
		else
			FLA_Copy_external( C_save, C_test );
		
		time = FLA_Clock();

		libfla_test_hemm_impl( impl, side, uplo, alpha, A_test, B_test, beta, C_test );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Copy the solution to flat matrix X.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_flatten( C_test, C );
	}
	else
    {
		// No action needed since C_test and C refer to the same object.
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
	     impl == FLA_TEST_FLAT_BLK_VAR ||
	     impl == FLA_TEST_FLAT_UNB_EXT ||
	     impl == FLA_TEST_FLAT_BLK_EXT )
		libfla_test_hemm_cntl_free();

	// Compute the performance of the best experiment repeat.
	if ( side == FLA_LEFT )
		*perf = ( 1 * m * m * n ) / time_min / FLOPS_PER_UNIT_PERF;
	else
		*perf = ( 1 * m * n * n ) / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

	// Compute:
	//   y = C * x
	// and compare to
	//   z = ( beta * C_orig + alpha * A * B ) x      (side = left)
	//   z = ( beta * C_orig + alpha * B * A ) x      (side = right)
	FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, C, x, FLA_ZERO, y );

	if ( side == FLA_LEFT )
	{
		FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, B, x, FLA_ZERO, w );
		FLA_Hemv_external( uplo,             alpha,   A, w, FLA_ZERO, z );
	}
	else
	{
		FLA_Hemv_external( uplo,             FLA_ONE, A, x, FLA_ZERO, w );
		FLA_Gemv_external( FLA_NO_TRANSPOSE, alpha,   B, w, FLA_ZERO, z );
	}
	FLA_Gemv_external( FLA_NO_TRANSPOSE, beta, C_save, x, FLA_ONE, z );

	// Compute || y - z ||.
	//FLA_Axpy_external( FLA_MINUS_ONE, y, z );
	//FLA_Nrm2_external( z, norm );
	//FLA_Obj_extract_real_scalar( norm, residual );
	*residual = FLA_Max_elemwise_diff( y, z );

	// Free the supporting flat objects.
	FLA_Obj_free( &C_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
	FLA_Obj_free( &B );
	FLA_Obj_free( &C );
	FLA_Obj_free( &x );
	FLA_Obj_free( &y );
	FLA_Obj_free( &z );
	FLA_Obj_free( &w );
	FLA_Obj_free( &norm );
}



extern fla_scal_t* fla_scal_cntl_blas;
extern fla_gemm_t* fla_gemm_cntl_blas;
extern fla_hemm_t* fla_hemm_cntl_blas;

void libfla_test_hemm_cntl_create( unsigned int var,
                                   dim_t        b_alg_flat )
{
	int var_unb = FLA_UNB_VAR_OFFSET + var;
	int var_blk = FLA_BLK_VAR_OFFSET + var;

	hemm_cntl_bsize = FLA_Blocksize_create( b_alg_flat,
	                                        b_alg_flat,
	                                        b_alg_flat,
	                                        b_alg_flat );

	hemm_cntl_unb   = FLA_Cntl_hemm_obj_create( FLA_FLAT,
	                                            var_unb,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL );

	hemm_cntl_blk   = FLA_Cntl_hemm_obj_create( FLA_FLAT,
	                                            var_blk,
	                                            hemm_cntl_bsize,
	                                            fla_scal_cntl_blas,
	                                            fla_hemm_cntl_blas,
	                                            fla_gemm_cntl_blas,
	                                            fla_gemm_cntl_blas );
}



void libfla_test_hemm_cntl_free( void )
{
	FLA_Blocksize_free( hemm_cntl_bsize );

	FLA_Cntl_obj_free( hemm_cntl_unb );
	FLA_Cntl_obj_free( hemm_cntl_blk );
}



void libfla_test_hemm_impl( int       impl,
                            FLA_Side  side,
                            FLA_Uplo  uplo,
                            FLA_Obj   alpha,
                            FLA_Obj   A,
                            FLA_Obj   B,
                            FLA_Obj   beta,
                            FLA_Obj   C )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_Hemm( side, uplo, alpha, A, B, beta, C );
		break;

		case FLA_TEST_FLAT_FRONT_END:
		FLA_Hemm( side, uplo, alpha, A, B, beta, C );
		break;

		case FLA_TEST_FLAT_UNB_VAR:
		FLA_Hemm_internal( side, uplo, alpha, A, B, beta, C, hemm_cntl_unb );
		break;

		case FLA_TEST_FLAT_BLK_VAR:
		FLA_Hemm_internal( side, uplo, alpha, A, B, beta, C, hemm_cntl_blk );
		break;

		case FLA_TEST_FLAT_UNB_EXT:
		FLA_Hemm_external( side, uplo, alpha, A, B, beta, C );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

