/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  6
#define FIRST_VARIANT    1
#define LAST_VARIANT     1

// Static variables.
static char* op_str                   = "Apply up/downdating Q via UD UT transform";
static char* fla_front_str            = "FLA_Apply_QUD_UT";
//static char* fla_unb_var_str          = "";
//static char* fla_opt_var_str          = "";
//static char* fla_blk_var_str          = "";
static char* pc_str[NUM_PARAM_COMBOS] = { "" };
static test_thresh_t thresh           = { 1e-02, 1e-03,   // warn, pass for s
                                          1e-11, 1e-12,   // warn, pass for d
                                          1e-02, 1e-03,   // warn, pass for c
                                          1e-11, 1e-12 }; // warn, pass for z

// Local prototypes.
void libfla_test_apqudut_experiment( test_params_t params,
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
void libfla_test_apqudut_impl( int     impl,
                               FLA_Obj T, FLA_Obj W,
                                          FLA_Obj bR,
                               FLA_Obj C, FLA_Obj bC,
                               FLA_Obj D, FLA_Obj bD );

void libfla_test_apqudut( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_apqudut_experiment );
	}

}



void libfla_test_apqudut_experiment( test_params_t params,
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
	uinteger mB, mC, mD, n, n_rhs;
	integer   mB_input    = -1;
	integer   mC_input    = -4;
	integer   mD_input    = -4;
	integer   n_input     = -1;
	integer   n_rhs_input = -1;
	FLA_Obj      R_BD, R_BC, B, C, D, T, W;
	FLA_Obj      bR_BD, bR_BC, bB, bC, bD;
	FLA_Obj      bR_BD_save, bC_save, bD_save;

	// Determine the dimensions.
	if ( mB_input    < 0 ) mB    = p_cur / abs(mB_input);
	else                   mB    = p_cur;
	if ( mC_input    < 0 ) mC    = p_cur / abs(mC_input);
	else                   mC    = p_cur;
	if ( mD_input    < 0 ) mD    = p_cur / abs(mD_input);
	else                   mD    = p_cur;
	if ( n_input     < 0 ) n     = p_cur / abs(n_input);
	else                   n     = p_cur;
	if ( n_rhs_input < 0 ) n_rhs = p_cur / abs(n_rhs_input);
	else                   n_rhs = p_cur;

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], b_alg_flat, n, &T );

	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], mB, n, &B );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], mC, n, &C );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[2], mD, n, &D );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n,  n, &R_BC );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n,  n, &R_BD );

	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], mB, n_rhs, &bB );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[4], mC, n_rhs, &bC );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[5], mD, n_rhs, &bD );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[3], n,  n_rhs, &bR_BC );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[3], n,  n_rhs, &bR_BD );

	FLA_Apply_QUD_UT_create_workspace( T, bR_BD, &W );

	// Initialize the test matrices.
	FLA_Random_matrix( B );
	FLA_Random_matrix( C );
	FLA_Random_matrix( D );

	// Initialize the right-hand sides.
	FLA_Random_matrix( bB );
	FLA_Random_matrix( bC );
	FLA_Random_matrix( bD );

	// Intialize the test factorization.
	FLA_Set( FLA_ZERO, R_BD );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B, FLA_ONE, R_BD );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, D, FLA_ONE, R_BD );
	FLA_Chol( FLA_UPPER_TRIANGULAR, R_BD );

	// Initialize the solution factorization.
	FLA_Set( FLA_ZERO, R_BC );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B, FLA_ONE, R_BC );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, C, FLA_ONE, R_BC );
	FLA_Chol( FLA_UPPER_TRIANGULAR, R_BC );

	// Initialize the test right-hand side.
	FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, B, bB, FLA_ZERO, bR_BD );
	FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, D, bD, FLA_ONE,  bR_BD );
	FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, R_BD, bR_BD );

	// Initialize the solution right-hand side.
	FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, B, bB, FLA_ZERO, bR_BC );
	FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, C, bC, FLA_ONE,  bR_BC );
	FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, R_BC, bR_BC );

	// Perform the up/downdate on R_BD, C, D, and T.
	FLA_UDdate_UT( R_BD, C, D, T );

	// Save the original test right-hand sides to temporary objects.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, bR_BD, &bR_BD_save );
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, bC, &bC_save );
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, bD, &bD_save );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		FLA_Copy_external( bR_BD_save, bR_BD );
		FLA_Copy_external( bC_save, bC );
		FLA_Copy_external( bD_save, bD );

		time = FLA_Clock();

		libfla_test_apqudut_impl( impl, T, W,
		                                   bR_BD,
		                                C, bC,
		                                D, bD );
		
		time = FLA_Clock() - time;
		time_min = fla_min( time_min, time );
	}

	// Solve for the solutions of our two systems.
	FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
	                   FLA_ONE, R_BD, bR_BD );
	FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
	                   FLA_ONE, R_BC, bR_BC );

	// Compute the maximum element-wise difference between the solutions of
	// the two systems.
	*residual = FLA_Max_elemwise_diff( bR_BD, bR_BC );

	// Compute the performance of the best experiment repeat.
	*t = time_min;
  *perf = n * n_rhs * ( 2.0 * mC + 2.0 * mD + 0.5 * b_alg_flat + 0.5 ) /
	        time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( bR_BD ) ) *perf *= 4.0;

	// Free the supporting flat objects.
	FLA_Obj_free( &bR_BD_save );
	FLA_Obj_free( &bC_save );
	FLA_Obj_free( &bD_save );

	// Free the flat test matrices.
	FLA_Obj_free( &B );
	FLA_Obj_free( &C );
	FLA_Obj_free( &D );
	FLA_Obj_free( &R_BC );
	FLA_Obj_free( &R_BD );
	FLA_Obj_free( &bB );
	FLA_Obj_free( &bC );
	FLA_Obj_free( &bD );
	FLA_Obj_free( &bR_BC );
	FLA_Obj_free( &bR_BD );
	FLA_Obj_free( &T );
	FLA_Obj_free( &W );
}



void libfla_test_apqudut_impl( int     impl,
                               FLA_Obj T, FLA_Obj W,
                                          FLA_Obj bR_BD,
                               FLA_Obj C, FLA_Obj bC,
                               FLA_Obj D, FLA_Obj bD )
{
	switch ( impl )
	{
		case FLA_TEST_FLAT_FRONT_END:
		FLA_Apply_QUD_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
		                  T, W,
		                     bR_BD,
		                  C, bC,
		                  D, bD );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

