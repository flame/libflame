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
static char* op_str                   = "CAQR factorization via UT transform (incremental)";
static char* flash_front_str          = "FLASH_CAQR_UT_inc";
//static char* fla_front_str            = "";
//static char* fla_unb_var_str          = "";
//static char* fla_opt_var_str          = "";
//static char* fla_blk_var_str          = "";
static char* pc_str[NUM_PARAM_COMBOS] = { "" };
static test_thresh_t thresh           = { 1e-02, 1e-03,   // warn, pass for s
                                          1e-11, 1e-12,   // warn, pass for d
                                          1e-02, 1e-03,   // warn, pass for c
                                          1e-11, 1e-12 }; // warn, pass for z

// Local prototypes.
void libfla_test_caqrutinc_experiment( test_params_t params,
                                       unsigned int  var,
                                       char*         sc_str,
                                       FLA_Datatype  datatype,
                                       unsigned int  p,
                                       unsigned int  pci,
                                       unsigned int  n_repeats,
                                       signed int    impl,
                                       double*       perf,
                                       double*       residual );
void libfla_test_caqrutinc_impl( int     impl,
                                 dim_t   p,
                                 FLA_Obj A,
                                 FLA_Obj ATW,
                                 FLA_Obj R,
                                 FLA_Obj RTW );


void libfla_test_caqrutinc( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_caqrutinc_experiment );
	}

}



void libfla_test_caqrutinc_experiment( test_params_t params,
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
	unsigned int m, n;
	signed int   m_input    = -8;
	signed int   n_input    = -1;
	unsigned int p          = 4;  // p <= abs(m_input) must hold!
	FLA_Obj      A, x, b, y, norm;
	FLA_Obj      A_save;
	FLA_Obj      A_test, ATW_test;
	FLA_Obj      R_test, RTW_test;
	FLA_Obj      x_test, b_test;

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur * abs(m_input);
	else               m = p_cur;
	if ( n_input < 0 ) n = p_cur * abs(n_input);
	else               n = p_cur;

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, n, &A );

	// Initialize the test matrices.
	FLA_Random_matrix( A );

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_save );

	// Create vectors to form a linear system.
	FLA_Obj_create( datatype, n, 1, 0, 0, &x );
	FLA_Obj_create( datatype, n, 1, 0, 0, &y );
	FLA_Obj_create( datatype, m, 1, 0, 0, &b );

	// Create a real scalar object to hold the norm of A.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Create a random right-hand side vector.
	FLA_Random_matrix( b );

	// Use hierarchical matrices since we're testing the FLASH front-end.
	FLASH_CAQR_UT_inc_create_hier_matrices( p, A, 1, &b_flash, b_alg_hier,
	                                        &A_test, &ATW_test, &R_test, &RTW_test );
	FLASH_Obj_create_hier_copy_of_flat( b, 1, &b_flash, &b_test );
	FLASH_Obj_create_hier_copy_of_flat( x, 1, &b_flash, &x_test );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		FLASH_Obj_hierarchify( A_save, A_test );
		
		time = FLA_Clock();

		libfla_test_caqrutinc_impl( impl, p, A_test, ATW_test, R_test, RTW_test );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Perform a linear solve with the result.
	FLASH_CAQR_UT_inc_solve( p, A_test, ATW_test, R_test, RTW_test, b_test, x_test );
	FLASH_Obj_flatten( x_test, x );

	// Free the FLASH objects.
	FLASH_Obj_free( &A_test );
	FLASH_Obj_free( &ATW_test );
	FLASH_Obj_free( &R_test );
	FLASH_Obj_free( &RTW_test );
	FLASH_Obj_free( &b_test );
	FLASH_Obj_free( &x_test );

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
	FLA_Obj_free( &y );
	FLA_Obj_free( &b );
	FLA_Obj_free( &norm );
	FLA_Obj_free( &A_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
}



void libfla_test_caqrutinc_impl( int     impl,
                                 dim_t   p,
                                 FLA_Obj A,
                                 FLA_Obj ATW,
                                 FLA_Obj R,
                                 FLA_Obj RTW )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_CAQR_UT_inc( p, A, ATW, R, RTW );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

