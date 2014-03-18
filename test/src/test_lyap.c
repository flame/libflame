
#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 2
#define NUM_MATRIX_ARGS  2
#define FIRST_VARIANT    1
#define LAST_VARIANT     4

// Static variables.
static char* op_str                   = "Triangular Lyapunov equation solve";
static char* flash_front_str          = "FLASH_Lyap";
static char* fla_front_str            = "FLA_Lyap";
static char* fla_unb_var_str          = "unb_var";
static char* fla_opt_var_str          = "opt_var";
static char* fla_blk_var_str          = "blk_var";
static char* pc_str[NUM_PARAM_COMBOS] = { "n", "h" };
static test_thresh_t thresh           = { 1e-02, 1e-03,   // warn, pass for s
                                          1e-11, 1e-12,   // warn, pass for d
                                          1e-02, 1e-03,   // warn, pass for c
                                          1e-11, 1e-12 }; // warn, pass for z

static fla_lyap_t*      lyap_cntl_opt;
static fla_lyap_t*      lyap_cntl_unb;
static fla_lyap_t*      lyap_cntl_blk;
static fla_blocksize_t* lyap_cntl_bsize;

// Local prototypes.
void libfla_test_lyap_experiment( test_params_t params,
                                  unsigned int  var,
                                  char*         sc_str,
                                  FLA_Datatype  datatype,
                                  unsigned int  p_cur,
                                  unsigned int  pci,
                                  unsigned int  n_repeats,
                                  signed int    impl,
                                  double*       perf,
                                  double*       residual );
void libfla_test_lyap_impl( int         impl,
                            FLA_Trans   trans,
                            FLA_Obj     isgn,
                            FLA_Obj     A,
                            FLA_Obj     C,
                            FLA_Obj     scale );
void libfla_test_lyap_cntl_create( unsigned int var,
                                   dim_t        b_alg_flat );
void libfla_test_lyap_cntl_free( void );


void libfla_test_lyap( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_lyap_experiment );
	}

	if ( op.fla_front == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, NULL,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_FRONT_END,
		                       params, thresh, libfla_test_lyap_experiment );
	}

	if ( op.fla_unb_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_unb_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_UNB_VAR,
		                       params, thresh, libfla_test_lyap_experiment );
	}

	if ( op.fla_opt_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_opt_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_OPT_VAR,
		                       params, thresh, libfla_test_lyap_experiment );
	}

	if ( op.fla_blk_vars == ENABLE )
	{
		libfla_test_op_driver( fla_front_str, fla_blk_var_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_BLK_VAR,
		                       params, thresh, libfla_test_lyap_experiment );
	}
}



void libfla_test_lyap_experiment( test_params_t params,
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
	signed int   sign       = -1;
	FLA_Trans    trans;
	FLA_Obj      A, C, X, isgn, scale, norm;
	FLA_Obj      AXpXA;
	FLA_Obj      C_save;
	FLA_Obj      A_test, C_test;

/*
	// Return early for optimized unblocked variants since they don't exist.
	if ( impl == FLA_TEST_FLAT_OPT_VAR && pc_str[pci][0] == 'n' )
	{
		*perf     = 0.0;
		*residual = 0.0;
		return;
	}
*/
	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur / abs(m_input);
	else               m = p_cur;

	// Translate parameter characters to libflame constants.
	FLA_Param_map_char_to_flame_trans( &pc_str[pci][0], &trans );

	// Determine the sign.
	if ( 0 < sign ) isgn = FLA_ONE;
	else            isgn = FLA_MINUS_ONE;

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, m, &A );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], m, m, &C );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], m, m, &X );

	// Create the scale object and a temporary norm object.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &scale );
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Initialize the test matrices.
	FLA_Random_tri_matrix( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
	FLA_Random_matrix( C );

	// Tweak A and C.
	FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, A );
	FLA_Norm1( A, norm );
	FLA_Shift_diag( FLA_NO_CONJUGATE, norm, A );

	FLA_Hermitianize( FLA_UPPER_TRIANGULAR, C );

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, C, &C_save );

	// Create temporary matrices to hold intermediate products, used when
	// checking the solution.
	FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, C, &AXpXA );

	// Use hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_Obj_create_hier_copy_of_flat( A, 1, &b_flash, &A_test );
		FLASH_Obj_create_hier_copy_of_flat( C, 1, &b_flash, &C_test );
	}
	else
	{
		A_test = A;
		C_test = C;
	}

	// Create a control tree for the individual variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_lyap_cntl_create( var, b_alg_flat );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		if ( impl == FLA_TEST_HIER_FRONT_END )
			FLASH_Obj_hierarchify( C_save, C_test );
		else
			FLA_Copy_external( C_save, C_test );
		
		time = FLA_Clock();

		libfla_test_lyap_impl( impl, trans, isgn, A_test, C_test, scale );
		
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
		FLASH_Obj_free( &C_test );
	}

	// Free the control trees if we're testing the variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_lyap_cntl_free();

	// Compute the performance of the best experiment repeat.
	*perf = 2.0 / 3.0 * m * m * m / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

	// Compute || +/-C - (AX - XB) ||.
	FLA_Hermitianize( FLA_UPPER_TRIANGULAR, X );
	if ( trans == FLA_NO_TRANSPOSE )
	{
		FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,   FLA_ONE, A, X, FLA_ZERO, AXpXA );
		FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_CONJ_TRANSPOSE, FLA_ONE, X, A, FLA_ONE,  AXpXA );
	}
	else // if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
	{
		FLA_Gemm_external( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, A, X, FLA_ZERO, AXpXA );
		FLA_Gemm_external( FLA_NO_TRANSPOSE,   FLA_NO_TRANSPOSE, FLA_ONE, X, A, FLA_ONE,  AXpXA );
	}
	FLA_Scal_external( isgn, AXpXA );
	FLA_Axpy_external( FLA_MINUS_ONE, C_save, AXpXA );
	FLA_Norm1( AXpXA, norm );
	FLA_Obj_extract_real_scalar( norm, residual );

	// Free the supporting flat objects.
	FLA_Obj_free( &AXpXA );
	FLA_Obj_free( &C_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
	FLA_Obj_free( &C );
	FLA_Obj_free( &X );
	FLA_Obj_free( &norm );
	FLA_Obj_free( &scale );
}



extern fla_scal_t*  fla_scal_cntl_blas;
extern fla_gemm_t*  fla_gemm_cntl_blas;
extern fla_hemm_t*  fla_hemm_cntl_blas;
extern fla_her2k_t* fla_her2k_cntl_blas;
extern fla_sylv_t*  fla_sylv_cntl;
extern fla_lyap_t*  fla_lyap_cntl_leaf;

void libfla_test_lyap_cntl_create( unsigned int var,
                                   dim_t        b_alg_flat )
{
	int var_unb  = FLA_UNB_VAR_OFFSET + var;
	int var_opt  = FLA_OPT_VAR_OFFSET + var;
	int var_blk  = FLA_BLK_VAR_OFFSET + var;

	lyap_cntl_bsize = FLA_Blocksize_create( b_alg_flat, b_alg_flat, b_alg_flat, b_alg_flat );

	lyap_cntl_unb   = FLA_Cntl_lyap_obj_create( FLA_FLAT,
	                                            var_unb,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL );

	lyap_cntl_opt   = FLA_Cntl_lyap_obj_create( FLA_FLAT,
	                                            var_opt,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL,
	                                            NULL );

	lyap_cntl_blk   = FLA_Cntl_lyap_obj_create( FLA_FLAT,
	                                            var_blk,
	                                            lyap_cntl_bsize,
	                                            fla_scal_cntl_blas,
	                                            fla_lyap_cntl_leaf,
	                                            fla_sylv_cntl,
	                                            fla_gemm_cntl_blas,
	                                            fla_gemm_cntl_blas,
	                                            fla_hemm_cntl_blas,
	                                            fla_her2k_cntl_blas );
}



void libfla_test_lyap_cntl_free( void )
{
	FLA_Blocksize_free( lyap_cntl_bsize );

	FLA_Cntl_obj_free( lyap_cntl_unb );
	FLA_Cntl_obj_free( lyap_cntl_opt );
	FLA_Cntl_obj_free( lyap_cntl_blk );
}



void libfla_test_lyap_impl( int       impl,
                            FLA_Trans trans,
                            FLA_Obj   isgn,
                            FLA_Obj   A,
                            FLA_Obj   C,
                            FLA_Obj   scale )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_Lyap( trans, isgn, A, C, scale );
		break;

		case FLA_TEST_FLAT_FRONT_END:
		FLA_Lyap( trans, isgn, A, C, scale );
		break;

		case FLA_TEST_FLAT_UNB_VAR:
		FLA_Lyap_internal( trans, isgn, A, C, scale, lyap_cntl_unb );
		break;

		case FLA_TEST_FLAT_OPT_VAR:
		FLA_Lyap_internal( trans, isgn, A, C, scale, lyap_cntl_opt );
		break;

		case FLA_TEST_FLAT_BLK_VAR:
		FLA_Lyap_internal( trans, isgn, A, C, scale, lyap_cntl_blk );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

