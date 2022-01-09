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
#define LAST_VARIANT     5

// Static variables.
static char* op_str                   = "Reduction of Hermitian-definite eigenproblem to std. form";
static char* flash_front_str          = "FLASH_Eig_gest";
static char* fla_front_str            = "FLA_Eig_gest";
static char* fla_unb_var_str          = "unb_var";
static char* fla_opt_var_str          = "opt_var";
static char* fla_blk_var_str          = "blk_var";
static char* pc_str[NUM_PARAM_COMBOS] = { "il", "iu", "nl", "nu" };
static test_thresh_t thresh           = { 1e-02, 1e-03,   // warn, pass for s
                                          1e-11, 1e-12,   // warn, pass for d
                                          1e-02, 1e-03,   // warn, pass for c
                                          1e-11, 1e-12 }; // warn, pass for z

static fla_eig_gest_t*  eig_gest_cntl_opt;
static fla_eig_gest_t*  eig_gest_cntl_unb;
static fla_eig_gest_t*  eig_gest_cntl_blk;
static fla_blocksize_t* eig_gest_cntl_bsize;

// Local prototypes.
void libfla_test_eig_gest_experiment( test_params_t params,
                                      unsigned int  var,
                                      char*         sc_str,
                                      FLA_Datatype  datatype,
                                      uinteger  p_cur,
                                      unsigned int  pci,
                                      unsigned int  n_repeats,
                                      signed int    impl,
                                      double*       perf,
                                      double*       t,
                                      double*       residual );
void libfla_test_eig_gest_impl( int         impl,
                                FLA_Uplo    inv,
                                FLA_Inv uplo,
                                FLA_Obj     A,
                                FLA_Obj     Y,
                                FLA_Obj     B );
void libfla_test_eig_gest_cntl_create( unsigned int var,
                                       dim_t        b_alg_flat );
void libfla_test_eig_gest_cntl_free( void );


void libfla_test_eig_gest( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_eig_gest_experiment );
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
		                       params, thresh, libfla_test_eig_gest_experiment );
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
		                       params, thresh, libfla_test_eig_gest_experiment );
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
		                       params, thresh, libfla_test_eig_gest_experiment );
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
		                       params, thresh, libfla_test_eig_gest_experiment );
	}
}



void libfla_test_eig_gest_experiment( test_params_t params,
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
	uinteger m;
	integer   m_input    = -1;
	FLA_Uplo     inv;
	FLA_Uplo     uplo;
	FLA_Obj      A, B, Y, norm;
	FLA_Obj      A_save, B_save;
	FLA_Obj      A_test, B_test, Y_test;

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur / abs(m_input);
	else               m = p_cur;

	// Translate parameter characters to libflame constants.
	FLA_Param_map_char_to_flame_inv( &pc_str[pci][0], &inv );
	FLA_Param_map_char_to_flame_uplo( &pc_str[pci][1], &uplo );

	if ( inv == FLA_NO_INVERSE &&
         ( ( impl == FLA_TEST_FLAT_UNB_VAR && var == 3 ) ||
	       ( impl == FLA_TEST_FLAT_OPT_VAR && var == 3 ) ||
	       ( impl == FLA_TEST_FLAT_BLK_VAR && var == 3 ) )
       )
	{
		*perf     = 0.0;
		*residual = 0.0;
		return;
	}

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, m, &A );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], m, m, &Y );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], m, m, &B );

	// Initialize the test matrices.
	FLA_Random_spd_matrix( uplo, A );
    FLA_Scalr( uplo, FLA_TWO, A );
	FLA_Hermitianize( uplo, A );

	FLA_Random_spd_matrix( uplo, B );
    FLA_Scalr( uplo, FLA_TWO, B );
	FLA_Chol( uplo, B );

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_save );
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, B, &B_save );

	// Create a real scalar object to hold the norm of A.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Use hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_create_hier_copy_of_flat( A, 1, &b_flash, &A_test );
		FLASH_Obj_create_hier_copy_of_flat( Y, 1, &b_flash, &Y_test );
		FLASH_Obj_create_hier_copy_of_flat( B, 1, &b_flash, &B_test );
	}
	else
	{
		A_test = A;
		Y_test = Y;
		B_test = B;
	}

	// Create a control tree for the individual variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_eig_gest_cntl_create( var, b_alg_flat );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		if ( impl == FLA_TEST_HIER_FRONT_END )
		{
			FLASH_Obj_hierarchify( A_save, A_test );
			FLASH_Obj_hierarchify( B_save, B_test );
		}
		else
		{
			FLA_Copy_external( A_save, A_test );
			FLA_Copy_external( B_save, B_test );
		}
		
		time = FLA_Clock();
		
		libfla_test_eig_gest_impl( impl, inv, uplo, A_test, Y_test, B_test );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Check our solution.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLA_Trans trans_left, trans_right;
		
		FLASH_Hermitianize( uplo, A_test );
		
		if ( ( inv == FLA_NO_INVERSE && uplo == FLA_LOWER_TRIANGULAR ) ||
		     ( inv == FLA_INVERSE    && uplo == FLA_UPPER_TRIANGULAR ) )
		{
			trans_left  = FLA_CONJ_TRANSPOSE;
			trans_right = FLA_NO_TRANSPOSE;
		}
		else
		{
			trans_left  = FLA_NO_TRANSPOSE;
			trans_right = FLA_CONJ_TRANSPOSE;
		}

		if ( inv == FLA_NO_INVERSE )
		{
			FLASH_Trsm( FLA_LEFT, uplo, trans_left, FLA_NONUNIT_DIAG,
			            FLA_ONE, B_test, A_test );
			FLASH_Trsm( FLA_RIGHT, uplo, trans_right, FLA_NONUNIT_DIAG,
			            FLA_ONE, B_test, A_test );
		}
		else // if ( inv == FLA_INVERSE )
		{
			FLASH_Trmm( FLA_LEFT, uplo, trans_left, FLA_NONUNIT_DIAG,
			            FLA_ONE, B_test, A_test );
			FLASH_Trmm( FLA_RIGHT, uplo, trans_right, FLA_NONUNIT_DIAG,
			            FLA_ONE, B_test, A_test );
		}
		FLASH_Obj_flatten( A_test, A );
	}
	else
	{
		FLA_Trans trans_left, trans_right;

		FLA_Hermitianize( uplo, A_test );

		if ( ( inv == FLA_NO_INVERSE && uplo == FLA_LOWER_TRIANGULAR ) ||
		     ( inv == FLA_INVERSE    && uplo == FLA_UPPER_TRIANGULAR ) )
		{
			trans_left  = FLA_CONJ_TRANSPOSE;
			trans_right = FLA_NO_TRANSPOSE;
		}
		else
		{
			trans_left  = FLA_NO_TRANSPOSE;
			trans_right = FLA_CONJ_TRANSPOSE;
		}

		if ( inv == FLA_NO_INVERSE )
		{
			FLA_Trsm( FLA_LEFT, uplo, trans_left, FLA_NONUNIT_DIAG,
			          FLA_ONE, B_test, A_test );
			FLA_Trsm( FLA_RIGHT, uplo, trans_right, FLA_NONUNIT_DIAG,
		              FLA_ONE, B_test, A_test );
		}
		else // if ( inv == FLA_INVERSE )
		{
			FLA_Trmm( FLA_LEFT, uplo, trans_left, FLA_NONUNIT_DIAG,
			          FLA_ONE, B_test, A_test );
			FLA_Trmm( FLA_RIGHT, uplo, trans_right, FLA_NONUNIT_DIAG,
		              FLA_ONE, B_test, A_test );
		}
	}

	// Free the hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_free( &A_test );
		FLASH_Obj_free( &Y_test );
		FLASH_Obj_free( &B_test );
	}

	// Free the control trees if we're testing the variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_eig_gest_cntl_free();

	// Compute the performance of the best experiment repeat.
	*t = time_min;
  *perf = 1.0 * m * m * m / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

	// Compute the residual.
	FLA_Axpy_external( FLA_MINUS_ONE, A_save, A );
	FLA_Norm1( A, norm );
	FLA_Obj_extract_real_scalar( norm, residual );

	// Free the supporting flat objects.
	FLA_Obj_free( &norm );
	FLA_Obj_free( &A_save );
	FLA_Obj_free( &B_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
	FLA_Obj_free( &Y );
	FLA_Obj_free( &B );
}



extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_axpy_t*  fla_axpy_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_gemm_t*  fla_gemm_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_hemm_t*  fla_hemm_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_her2k_t* fla_her2k_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_trmm_t*  fla_trmm_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_trsm_t*  fla_trsm_cntl_blas;

void libfla_test_eig_gest_cntl_create( unsigned int var,
                                       dim_t        b_alg_flat )
{
	int var_unb = FLA_UNB_VAR_OFFSET + var;
	int var_opt = FLA_OPT_VAR_OFFSET + var;
	int var_blk = FLA_BLK_VAR_OFFSET + var;

	eig_gest_cntl_bsize = FLA_Blocksize_create( b_alg_flat, b_alg_flat, b_alg_flat, b_alg_flat );

	eig_gest_cntl_unb   = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
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
	                                                    NULL,
	                                                    NULL );

	eig_gest_cntl_opt   = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
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
	                                                    NULL,
	                                                    NULL );

	eig_gest_cntl_blk   = FLA_Cntl_eig_gest_obj_create( FLA_FLAT,
	                                                    var_blk,
	                                                    eig_gest_cntl_bsize,
	                                                    eig_gest_cntl_opt,
	                                                    fla_axpy_cntl_blas,
	                                                    fla_axpy_cntl_blas,
	                                                    fla_gemm_cntl_blas,
	                                                    fla_gemm_cntl_blas,
	                                                    fla_gemm_cntl_blas,
	                                                    fla_hemm_cntl_blas,
	                                                    fla_her2k_cntl_blas,
	                                                    fla_trmm_cntl_blas,
	                                                    fla_trmm_cntl_blas,
	                                                    fla_trsm_cntl_blas,
	                                                    fla_trsm_cntl_blas );
}



void libfla_test_eig_gest_cntl_free( void )
{
	FLA_Blocksize_free( eig_gest_cntl_bsize );

	FLA_Cntl_obj_free( eig_gest_cntl_unb );
	FLA_Cntl_obj_free( eig_gest_cntl_opt );
	FLA_Cntl_obj_free( eig_gest_cntl_blk );
}



void libfla_test_eig_gest_impl( int impl,
                            FLA_Inv inv,
                            FLA_Uplo uplo,
                            FLA_Obj A,
                            FLA_Obj Y,
                            FLA_Obj B )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_Eig_gest( inv, uplo, A, B );
		break;

		case FLA_TEST_FLAT_FRONT_END:
		FLA_Eig_gest( inv, uplo, A, B );
		break;

		case FLA_TEST_FLAT_UNB_VAR:
		FLA_Eig_gest_internal( inv, uplo, A, Y, B, eig_gest_cntl_unb );
		break;

		case FLA_TEST_FLAT_OPT_VAR:
		FLA_Eig_gest_internal( inv, uplo, A, Y, B, eig_gest_cntl_opt );
		break;

		case FLA_TEST_FLAT_BLK_VAR:
		FLA_Eig_gest_internal( inv, uplo, A, Y, B, eig_gest_cntl_blk );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

