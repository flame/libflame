/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  4
#define FIRST_VARIANT    1
#define LAST_VARIANT     1

// Static variables.
static char* op_str                   = "Up/downdate via UD UT transform";
static char* fla_front_str            = "FLA_UDdate_UT";
static char* fla_unb_var_str          = "unb_var";
static char* fla_opt_var_str          = "opt_var";
static char* fla_blk_var_str          = "blk_var";
static char* pc_str[NUM_PARAM_COMBOS] = { "" };
static test_thresh_t thresh           = { 1e-02, 1e-03,   // warn, pass for s
                                          1e-11, 1e-12,   // warn, pass for d
                                          1e-02, 1e-03,   // warn, pass for c
                                          1e-11, 1e-12 }; // warn, pass for z

static fla_apqudut_t*   apqudut_cntl_blk;
static fla_uddateut_t*  uddateut_cntl_opt;
static fla_uddateut_t*  uddateut_cntl_unb;
static fla_uddateut_t*  uddateut_cntl_blk;
static fla_blocksize_t* uddateut_cntl_bsize;

// Local prototypes.
void libfla_test_uddateut_experiment( test_params_t params,
                                      unsigned int  var,
                                      char*         sc_str,
                                      FLA_Datatype  datatype,
                                      uinteger  p,
                                      unsigned int  pci,
                                      unsigned int  n_repeats,
                                      signed int    impl,
                                      double*       perf,
                                      double*       residual );
void libfla_test_uddateut_impl( int     impl,
                                FLA_Obj R,
                                FLA_Obj C,
                                FLA_Obj D,
                                FLA_Obj T );
void libfla_test_uddateut_cntl_create( unsigned int var,
                                       dim_t        b_alg_flat );
void libfla_test_uddateut_cntl_free( void );


void libfla_test_uddateut( FILE* output_stream, test_params_t params, test_op_t op )
{
	libfla_test_output_info( "--- %s ---\n", op_str );
	libfla_test_output_info( "\n" );

	if ( op.fla_front == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, NULL,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_FRONT_END,
		                       params, thresh, libfla_test_uddateut_experiment );
	}

	if ( op.fla_unb_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_unb_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_UNB_VAR,
		                       params, thresh, libfla_test_uddateut_experiment );
	}

	if ( op.fla_opt_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_opt_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_OPT_VAR,
		                       params, thresh, libfla_test_uddateut_experiment );
	}

	if ( op.fla_blk_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_blk_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_BLK_VAR,
		                       params, thresh, libfla_test_uddateut_experiment );
	}

}



void libfla_test_uddateut_experiment( test_params_t params,
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
	dim_t        b_alg_flat = params.b_alg_flat;
	double       time_min   = 1e9;
	double       time;
	unsigned int i;
	uinteger mB, mC, mD, n;
	integer   mB_input   = -1;
	integer   mC_input   = -4;
	integer   mD_input   = -4;
	integer   n_input    = -1;
	FLA_Obj      B, C, D, T, R, RR, E, EE;
	FLA_Obj      R_save, C_save, D_save;

	// Determine the dimensions.
	if ( mB_input < 0 ) mB = p_cur / abs(mB_input);
	else                mB = p_cur;
	if ( mC_input < 0 ) mC = p_cur / abs(mC_input);
	else                mC = p_cur;
	if ( mD_input < 0 ) mD = p_cur / abs(mD_input);
	else                mD = p_cur;
	if ( n_input  < 0 ) n  = p_cur / abs(n_input);
	else                n  = p_cur;

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], mB, n, &B );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], mC, n, &C );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], mD, n, &D );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n,  n, &R );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n,  n, &E );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n,  n, &RR );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n,  n, &EE );

	if ( impl == FLA_TEST_FLAT_FRONT_END ||
	     ( impl == FLA_TEST_FLAT_BLK_VAR && var == 1 ) )
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[3], b_alg_flat, n, &T );
	else
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[3], n, n, &T );

	// Initialize the test matrices.
	FLA_Random_matrix( B );
	FLA_Random_matrix( C );
	FLA_Random_matrix( D );

	// Intialize the test factorization.
	FLA_Set( FLA_ZERO, R );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B, FLA_ONE, R );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, D, FLA_ONE, R );
	FLA_Chol( FLA_UPPER_TRIANGULAR, R );

	// Initialize the solution factorization.
	FLA_Set( FLA_ZERO, E );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B, FLA_ONE, E );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, C, FLA_ONE, E );
	FLA_Chol( FLA_UPPER_TRIANGULAR, E );

	// Save the original test matrices to temporary objects.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, R, &R_save );
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, C, &C_save );
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, D, &D_save );

	// Create a control tree for the individual variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_uddateut_cntl_create( var, b_alg_flat );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		FLA_Copy_external( R_save, R );
		FLA_Copy_external( C_save, C );
		FLA_Copy_external( D_save, D );

		time = FLA_Clock();

		libfla_test_uddateut_impl( impl, R, C, D, T );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Compute R'R and E'E.
	FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, R, R, FLA_ZERO, RR );
	FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, E, E, FLA_ZERO, EE );

	// Free the control trees if we're testing the variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_uddateut_cntl_free();

	// Compute the performance of the best experiment repeat.
	*perf = 2.0 * ( ( mC + mD ) * n * n +
	                ( mC + mD ) * n * 6.0 ) / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( R ) ) *perf *= 4.0;

	// Compute the maximum element-wise difference between R'R and E'E and use
	// this instead of the residual.
	*residual = FLA_Max_elemwise_diff( RR, EE );

	// Free the supporting flat objects.
	FLA_Obj_free( &R_save );
	FLA_Obj_free( &C_save );
	FLA_Obj_free( &D_save );

	// Free the flat test matrices.
	FLA_Obj_free( &B );
	FLA_Obj_free( &C );
	FLA_Obj_free( &D );
	FLA_Obj_free( &T );
	FLA_Obj_free( &R );
	FLA_Obj_free( &RR );
	FLA_Obj_free( &E );
	FLA_Obj_free( &EE );
}



extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_axpyt_t* fla_axpyt_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_copyt_t* fla_copyt_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_gemm_t*  fla_gemm_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_trmm_t*  fla_trmm_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_trsm_t*  fla_trsm_cntl_blas;

void libfla_test_uddateut_cntl_create( unsigned int var,
                                       dim_t        b_alg_flat )
{
	int var_unb  = FLA_UNB_VAR_OFFSET + var;
	int var_opt  = FLA_OPT_VAR_OFFSET + var;
	int var_blk  = FLA_BLK_VAR_OFFSET + var;

	uddateut_cntl_bsize = FLA_Blocksize_create( b_alg_flat, b_alg_flat, b_alg_flat, b_alg_flat );

	apqudut_cntl_blk  = FLA_Cntl_apqudut_obj_create( FLA_FLAT,
	                                                 FLA_BLOCKED_VARIANT1,
	                                                 uddateut_cntl_bsize,
	                                                 NULL,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_gemm_cntl_blas,
	                                                 fla_trsm_cntl_blas,
	                                                 fla_copyt_cntl_blas,
	                                                 fla_axpyt_cntl_blas );

	uddateut_cntl_unb   = FLA_Cntl_uddateut_obj_create( FLA_FLAT,
	                                                    var_unb,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL );

	uddateut_cntl_opt   = FLA_Cntl_uddateut_obj_create( FLA_FLAT,
	                                                    var_opt,
	                                                    NULL,
	                                                    NULL,
	                                                    NULL );

	uddateut_cntl_blk   = FLA_Cntl_uddateut_obj_create( FLA_FLAT,
	                                                    var_blk,
	                                                    uddateut_cntl_bsize,
	                                                    uddateut_cntl_opt,
	                                                    apqudut_cntl_blk );
}



void libfla_test_uddateut_cntl_free( void )
{
	FLA_Blocksize_free( uddateut_cntl_bsize );

	FLA_Cntl_obj_free( apqudut_cntl_blk );
	FLA_Cntl_obj_free( uddateut_cntl_unb );
	FLA_Cntl_obj_free( uddateut_cntl_opt );
	FLA_Cntl_obj_free( uddateut_cntl_blk );
}



void libfla_test_uddateut_impl( int     impl,
                                FLA_Obj R,
                                FLA_Obj C,
                                FLA_Obj D,
                                FLA_Obj T )
{
	switch ( impl )
	{
		case FLA_TEST_FLAT_FRONT_END:
		FLA_UDdate_UT( R, C, D, T );
		break;

		case FLA_TEST_FLAT_UNB_VAR:
		FLA_UDdate_UT_internal( R, C, D, T, uddateut_cntl_unb );
		break;

		case FLA_TEST_FLAT_OPT_VAR:
		FLA_UDdate_UT_internal( R, C, D, T, uddateut_cntl_opt );
		break;

		case FLA_TEST_FLAT_BLK_VAR:
		FLA_UDdate_UT_internal( R, C, D, T, uddateut_cntl_blk );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

