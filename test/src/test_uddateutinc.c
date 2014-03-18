
#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1
#define FIRST_VARIANT    1
#define LAST_VARIANT     1

// Static variables.
static char* op_str                   = "Up/downdate via UD UT transform (incremental)";
static char* flash_front_str          = "FLASH_UDdate_UT_inc";
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
void libfla_test_uddateutinc_experiment( test_params_t params,
                                         unsigned int  var,
                                         char*         sc_str,
                                         FLA_Datatype  datatype,
                                         unsigned int  p,
                                         unsigned int  pci,
                                         unsigned int  n_repeats,
                                         signed int    impl,
                                         double*       perf,
                                         double*       residual );
void libfla_test_uddateutinc_impl( int     impl,
                                   FLA_Obj R,
                                   FLA_Obj C,
                                   FLA_Obj D,
                                   FLA_Obj T,
                                   FLA_Obj W );


void libfla_test_uddateutinc( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_uddateutinc_experiment );
	}

}



void libfla_test_uddateutinc_experiment( test_params_t params,
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
	unsigned int mB, mC, mD, n;
	signed int   mB_input   = -1;
	signed int   mC_input   = -4;
	signed int   mD_input   = -4;
	signed int   n_input    = -1;
	FLA_Obj      C, D, T, W, R, E, RR, EE;
	FLA_Obj      B_flat, C_flat, D_flat, R_flat, E_flat;

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
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], mB, n, &B_flat );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], mC, n, &C_flat );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], mD, n, &D_flat );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n,  n, &R_flat );
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], n,  n, &E_flat );

	// Initialize the test matrices.
	FLA_Random_matrix( B_flat );
	FLA_Random_matrix( C_flat );
	FLA_Random_matrix( D_flat );

	// Intialize the test factorization.
	FLA_Set( FLA_ZERO, R_flat );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B_flat, FLA_ONE, R_flat );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, D_flat, FLA_ONE, R_flat );
	FLA_Chol( FLA_UPPER_TRIANGULAR, R_flat );

	// Initialize the solution factorization.
	FLA_Set( FLA_ZERO, E_flat );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, B_flat, FLA_ONE, E_flat );
	FLA_Herk_external( FLA_UPPER_TRIANGULAR, FLA_CONJ_TRANSPOSE, FLA_ONE, C_flat, FLA_ONE, E_flat );
	FLA_Chol( FLA_UPPER_TRIANGULAR, E_flat );

	// Create the hierarchical matrices for the UDdate_UT_inc operation.
	FLASH_UDdate_UT_inc_create_hier_matrices( R_flat, C_flat, D_flat, 1, &b_flash, b_alg_hier,
	                                          &R, &C, &D, &T, &W );
	FLASH_Obj_create_hier_copy_of_flat( E_flat, 1, &b_flash, &E );
	FLASH_Obj_create_hier_copy_of_flat( R_flat, 1, &b_flash, &RR );
	FLASH_Obj_create_hier_copy_of_flat( E_flat, 1, &b_flash, &EE );

	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
		FLASH_Obj_hierarchify( R_flat, R );
		FLASH_Obj_hierarchify( C_flat, C );
		FLASH_Obj_hierarchify( D_flat, D );

		time = FLA_Clock();

		libfla_test_uddateutinc_impl( impl, R, C, D, T, W );
		
		time = FLA_Clock() - time;
		time_min = min( time_min, time );
	}

	// Compute R'R and E'E.
	FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, R, R, FLA_ZERO, RR );
	FLASH_Gemm( FLA_CONJ_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE, E, E, FLA_ZERO, EE );

	// Compute the performance of the best experiment repeat.
	*perf = 2.0 * ( ( mC + mD ) * n * n +
	                ( mC + mD ) * n * 6.0 ) / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( R ) ) *perf *= 4.0;

	// Compute the maximum element-wise difference between R'R and E'E and use
	// this instead of the residual.
	*residual = FLASH_Max_elemwise_diff( RR, EE );

	// Free the supporting flat objects.
	FLA_Obj_free( &B_flat );
	FLA_Obj_free( &C_flat );
	FLA_Obj_free( &D_flat );
	FLA_Obj_free( &R_flat );
	FLA_Obj_free( &E_flat );

	// Free the hierarchical test matrices.
	FLASH_Obj_free( &C );
	FLASH_Obj_free( &D );
	FLASH_Obj_free( &T );
	FLASH_Obj_free( &W );
	FLASH_Obj_free( &R );
	FLASH_Obj_free( &E );
	FLASH_Obj_free( &RR );
	FLASH_Obj_free( &EE );
}



void libfla_test_uddateutinc_impl( int     impl,
                                   FLA_Obj R,
                                   FLA_Obj C,
                                   FLA_Obj D,
                                   FLA_Obj T,
                                   FLA_Obj W )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_UDdate_UT_inc( R, C, D, T, W );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

