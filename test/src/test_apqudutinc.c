/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1
#define FIRST_VARIANT    1
#define LAST_VARIANT     1

// Static variables.
static char* op_str                   = "Apply up/downdating Q via UD UT transform (incremental)";
static char* flash_front_str          = "FLASH_Apply_QUD_UT_inc";
//static char* fla_front_str            = "";
//static char* fla_unb_var_str          = "";
//static char* fla_opt_var_str          = "";
//static char* fla_blk_var_str          = "";
static char* pc_str[NUM_PARAM_COMBOS] = { "" };
static test_thresh_t thresh           = { 1e-03, 1e-04,   // warn, pass for s
                                          1e-12, 1e-13,   // warn, pass for d
                                          1e-03, 1e-04,   // warn, pass for c
                                          1e-12, 1e-13 }; // warn, pass for z

// Local prototypes.
void libfla_test_apqudutinc_experiment( test_params_t params,
                                        unsigned int  var,
                                        char*         sc_str,
                                        FLA_Datatype  datatype,
                                        unsigned int  p,
                                        unsigned int  pci,
                                        unsigned int  n_repeats,
                                        signed int    impl,
                                        double*       perf,
                                        double*       residual );
void libfla_test_apqudutinc_impl( int     impl,
                                  FLA_Obj T, FLA_Obj W,
                                             FLA_Obj bR,
                                  FLA_Obj C, FLA_Obj bC,
                                  FLA_Obj D, FLA_Obj bD );


void libfla_test_apqudutinc( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_apqudutinc_experiment );
	}

}



void libfla_test_apqudutinc_experiment( test_params_t params,
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
	dim_t        b_alg_hier = params.b_alg_hier;
	double       time_min   = 1e9;
	double       time;
	unsigned int i;
	unsigned int mB, mC, mD, n, n_rhs;
	signed int   mB_input    = -1;
	signed int   mC_input    = -4;
	signed int   mD_input    = -4;
	signed int   n_input     = -1;
	signed int   n_rhs_input = -1;
	FLA_Obj      R_BD, R_BC, B, C, D, T, W, W2;
	FLA_Obj      bR_BD, bR_BC, bB, bC, bD;
	FLA_Obj      R_BD_flat, R_BC_flat, B_flat, C_flat, D_flat;
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
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], mB, n, &B_flat );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], mC, n, &C_flat );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], mD, n, &D_flat );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n,  n, &R_BC_flat );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n,  n, &R_BD_flat );

	// Initialize the test matrices.
	FLA_Random_matrix( B_flat );
	FLA_Random_matrix( C_flat );
	FLA_Random_matrix( D_flat );

	// Intialize the test factorization.
	FLA_Set( FLA_ZERO, R_BD_flat );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
	                   FLA_ONE, B_flat, FLA_ONE, R_BD_flat );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
	                   FLA_ONE, D_flat, FLA_ONE, R_BD_flat );
	FLA_Chol( FLA_UPPER_TRIANGULAR, R_BD_flat );

	// Initialize the solution factorization.
	FLA_Set( FLA_ZERO, R_BC_flat );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
	                   FLA_ONE, B_flat, FLA_ONE, R_BC_flat );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE,
	                   FLA_ONE, C_flat, FLA_ONE, R_BC_flat );
	FLA_Chol( FLA_UPPER_TRIANGULAR, R_BC_flat );

	// Create the hierarchical matrices for the UDdate_UT_inc operation.
	FLASH_UDdate_UT_inc_create_hier_matrices( R_BD_flat, C_flat, D_flat,
	                                          1, &b_flash, b_alg_hier,
	                                          &R_BD, &C, &D, &T, &W );
	FLASH_Obj_create_hier_copy_of_flat( B_flat, 1, &b_flash, &B );
	FLASH_Obj_create_hier_copy_of_flat( R_BC_flat, 1, &b_flash, &R_BC );

	// Create the hierarchical right-hand side matrices.
	FLASH_Obj_create( datatype, mB, n_rhs, 1, &b_flash, &bB );
	FLASH_Obj_create( datatype, mC, n_rhs, 1, &b_flash, &bC );
	FLASH_Obj_create( datatype, mD, n_rhs, 1, &b_flash, &bD );
	FLASH_Obj_create( datatype, n,  n_rhs, 1, &b_flash, &bR_BC );
	FLASH_Obj_create( datatype, n,  n_rhs, 1, &b_flash, &bR_BD );

	// Initialize the right-hand side matrices.
	FLASH_Random_matrix( bB );
	FLASH_Random_matrix( bC );
	FLASH_Random_matrix( bD );

	FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, B, bB, FLA_ZERO, bR_BD );
	FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, D, bD, FLA_ONE,  bR_BD );
	FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, R_BD, bR_BD );

	FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, B, bB, FLA_ZERO, bR_BC );
	FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, C, bC, FLA_ONE,  bR_BC );
	FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_NONUNIT_DIAG, FLA_ONE, R_BC, bR_BC );

	// Perform the updowndate on R_BD, C, and D.
	FLASH_UDdate_UT_inc( R_BD, C, D, T, W );

	// Allocate workspace for use when applying Q'.
	FLASH_Apply_QUD_UT_inc_create_workspace( T, bR_BD, &W2 );

	// Make copies of the right-hand side matrices so we can restore them
	// before each experiment repitition.
	FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, bR_BD, &bR_BD_save );
	FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, bC, &bC_save );
	FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, bD, &bD_save );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		FLASH_Copy( bR_BD_save, bR_BD );
		FLASH_Copy( bC_save, bC );
		FLASH_Copy( bD_save, bD );

		time = FLA_Clock();

		libfla_test_apqudutinc_impl( impl,
		                             T, W2,
		                                bR_BD,
		                             C, bC,
		                             D, bD );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Solve for the solutions of our two systems.
	FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
	            FLA_ONE, R_BD, bR_BD );
	FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
	            FLA_ONE, R_BC, bR_BC );

	// Compute the maximum element-wise difference between the solutions of
	// the two systems.
	*residual = FLASH_Max_elemwise_diff( bR_BD, bR_BC );

	// Compute the performance of the best experiment repeat.
	*perf = n * n_rhs * ( 2.0 * mC + 2.0 * mD + 0.5 * b_alg_hier + 0.5 ) /
	        time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( R_BD ) ) *perf *= 4.0;

	// Free the temporary hierarchical matrices.
	FLASH_Obj_free( &bR_BD_save );
	FLASH_Obj_free( &bC_save );
	FLASH_Obj_free( &bD_save );

	// Free the supporting flat objects.
	FLA_Obj_free( &B_flat );
	FLA_Obj_free( &C_flat );
	FLA_Obj_free( &D_flat );
	FLA_Obj_free( &R_BC_flat );
	FLA_Obj_free( &R_BD_flat );

	// Free the hierarchical test matrices.
	FLASH_Obj_free( &B );
	FLASH_Obj_free( &C );
	FLASH_Obj_free( &D );
	FLASH_Obj_free( &R_BC );
	FLASH_Obj_free( &R_BD );

	FLASH_Obj_free( &T );
	FLASH_Obj_free( &W );
	FLASH_Obj_free( &W2 );

	FLASH_Obj_free( &bB );
	FLASH_Obj_free( &bC );
	FLASH_Obj_free( &bD );
	FLASH_Obj_free( &bR_BC );
	FLASH_Obj_free( &bR_BD );
}



void libfla_test_apqudutinc_impl( int     impl,
                                  FLA_Obj T, FLA_Obj W,
                                             FLA_Obj bR,
                                  FLA_Obj C, FLA_Obj bC,
                                  FLA_Obj D, FLA_Obj bD )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_Apply_QUD_UT_inc( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
		                        T, W,
		                           bR,
		                        C, bC,
		                        D, bD );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

