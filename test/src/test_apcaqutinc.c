/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "test_libflame.h"

//#define NUM_PARAM_COMBOS 16
#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  2
#define FIRST_VARIANT    1
#define LAST_VARIANT     1

// Static variables.
static char* op_str                   = "Apply CAQ via UT transform (incremental)";
static char* flash_front_str          = "FLASH_Apply_CAQ_UT_inc";
//static char* fla_front_str            = "FLA_Apply_Q_UT_inc";
//static char* fla_unb_var_str          = "";
//static char* fla_opt_var_str          = "";
//static char* fla_blk_var_str          = "";
//static char* pc_str[NUM_PARAM_COMBOS] = { "lnfc", "lnfr", "lnbc", "lnbr",
//                                          "lhfc", "lhfr", "lhbc", "lhbr",
//                                          "rnfc", "rnfr", "rnbc", "rnbr",
//                                          "rhfc", "rhfr", "rhbc", "rhbr" };
static char* pc_str[NUM_PARAM_COMBOS] = { "lhfc" };
static test_thresh_t thresh           = { 1e-02, 1e-03,   // warn, pass for s
                                          1e-11, 1e-12,   // warn, pass for d
                                          1e-02, 1e-03,   // warn, pass for c
                                          1e-11, 1e-12 }; // warn, pass for z

// Local prototypes.
void libfla_test_apcaqutinc_experiment( test_params_t params,
                                        unsigned int  var,
                                        char*         sc_str,
                                        FLA_Datatype  datatype,
                                        unsigned int  p,
                                        unsigned int  pci,
                                        unsigned int  n_repeats,
                                        signed int    impl,
                                        double*       perf,
                                        double*       residual );
void libfla_test_apcaqutinc_impl( int        impl,
                                  dim_t      p,
                                  FLA_Side   side,
                                  FLA_Trans  trans,
                                  FLA_Direct direct,
                                  FLA_Store  storev,
                                  FLA_Obj    A,
                                  FLA_Obj    ATW,
                                  FLA_Obj    R,
                                  FLA_Obj    RTW,
                                  FLA_Obj    W,
                                  FLA_Obj    B );

void libfla_test_apcaqutinc( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_apcaqutinc_experiment );
	}
}



void libfla_test_apcaqutinc_experiment( test_params_t params,
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
	unsigned int min_m_n, k;
	unsigned int p;
	signed int   m_input;
	signed int   n_input;
	FLA_Side     side;
	FLA_Trans    trans;
	FLA_Direct   direct;
	FLA_Store    storev;
	FLA_Obj      A, X, B, Y, norm;
	FLA_Obj      B_save;
	FLA_Obj      A_test, ATW_test;
	FLA_Obj      R_test, RTW_test;
	FLA_Obj      W_test, X_test, B_test;

	// Translate parameter characters to libflame constants.
	FLA_Param_map_char_to_flame_side( &pc_str[pci][0], &side );
	FLA_Param_map_char_to_flame_trans( &pc_str[pci][1], &trans );
	FLA_Param_map_char_to_flame_direct( &pc_str[pci][2], &direct );
	FLA_Param_map_char_to_flame_storev( &pc_str[pci][3], &storev );

	// We want to make sure the Apply_Q_UT_inc routines work with rectangular
	// matrices. So we use m > n when testing with column-wise storage (via
	// QR factorization) and m < n when testing with row-wise storage (via
	// LQ factorization).
	if ( storev == FLA_COLUMNWISE )
	{
		m_input = -8;
		n_input = -1;
		p       = 4;  // p <= abs(m_input) must hold!
		//m_input = -1;
		//n_input = -1;
	}
	else // if ( storev == FLA_ROWWISE )
	{
	}

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur * abs(m_input);
	else               m = p_cur;
	if ( n_input < 0 ) n = p_cur * abs(n_input);
	else               n = p_cur;

	// Compute the minimum dimension.
	min_m_n = min( m, n );

	// Choose the size of B based on the storev parameter.
	if ( storev == FLA_COLUMNWISE ) k = m;
	//else                            k = n;

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, n, &A );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], m, 1, &B );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], n, 1, &X );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], n, 1, &Y );

	// Initialize the test matrices.
	FLA_Random_matrix( A );
	FLA_Random_matrix( B );

	// Use hierarchical matrices since we're testing the FLASH front-end.
	if ( storev == FLA_COLUMNWISE )
		FLASH_CAQR_UT_inc_create_hier_matrices( p, A, 1, &b_flash, b_alg_hier,
		                                        &A_test, &ATW_test, &R_test, &RTW_test );
	//else // if ( storev == FLA_ROWWISE )
	//  FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	FLASH_Obj_create_hier_copy_of_flat( B, 1, &b_flash, &B_test );
	FLASH_Obj_create_hier_copy_of_flat( X, 1, &b_flash, &X_test );
	FLASH_Apply_CAQ_UT_inc_create_workspace( p, RTW_test, B_test, &W_test );

	// Create a real scalar object to hold the norm of A.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, B, &B_save );

	// Compute a Householder factorization.
	if ( storev == FLA_COLUMNWISE ) FLASH_CAQR_UT_inc( p, A_test, ATW_test, R_test, RTW_test );
	//else                            FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		FLASH_Obj_hierarchify( B_save, B_test );

		time = FLA_Clock();

		libfla_test_apcaqutinc_impl( impl, p, side, trans, direct, storev,
		                             A_test, ATW_test, R_test, RTW_test, W_test, B_test );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Multiply by its conjugate-transpose to get what should be (near) identity
	// and then subtract from actual identity to get what should be (near) zero.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLA_Obj  RT, RB;
		FLA_Obj  BT, BB;

		FLASH_Part_create_2x1( R_test,   &RT,
		                                 &RB,    FLASH_Obj_scalar_width( R_test ), FLA_TOP );
		FLASH_Part_create_2x1( B_test,   &BT,
		                                 &BB,    FLASH_Obj_scalar_width( R_test ), FLA_TOP );

		FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
		            FLA_ONE, RT, BT );
		FLASH_Copy( BT, X_test );

		FLASH_Part_free_2x1( &RT,
		                     &RB );
		FLASH_Part_free_2x1( &BT,
		                     &BB );

		FLASH_Obj_flatten( X_test, X );

		FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A, X, FLA_MINUS_ONE, B );
		FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A, B, FLA_ZERO, Y );
	}

	// Free the hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_free( &A_test );
		FLASH_Obj_free( &ATW_test );
		FLASH_Obj_free( &R_test );
		FLASH_Obj_free( &RTW_test );
		FLASH_Obj_free( &W_test );
		FLASH_Obj_free( &B_test );
		FLASH_Obj_free( &X_test );
	}

	// Compute the norm of Y.
	FLA_Nrm2_external( Y, norm );
	FLA_Obj_extract_real_scalar( norm, residual );

	// Compute the performance of the best experiment repeat.
	*perf = (  4.0 *       m * n * 1 -
	           2.0 *       n * n * 1 ) / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

	// Free the supporting flat objects.
	FLA_Obj_free( &B_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
	FLA_Obj_free( &B );
	FLA_Obj_free( &X );
	FLA_Obj_free( &Y );
	FLA_Obj_free( &norm );
}



void libfla_test_apcaqutinc_impl( int        impl,
                                  dim_t      p,
                                  FLA_Side   side,
                                  FLA_Trans  trans,
                                  FLA_Direct direct,
                                  FLA_Store  storev,
                                  FLA_Obj    A,
                                  FLA_Obj    ATW,
                                  FLA_Obj    R,
                                  FLA_Obj    RTW,
                                  FLA_Obj    W,
                                  FLA_Obj    B )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_Apply_CAQ_UT_inc( p, side, trans, direct, storev,
		                        A, ATW, R, RTW, W, B );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

