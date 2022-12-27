/*
* Copyright (c) 2020 Advanced Micro Devices, Inc. All rights reserved.
* */
#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  1
#define FIRST_VARIANT    1
#define LAST_VARIANT     6

// Static variables.
static char* op_str                   = "LU factorization complete/incomplete w/o pivoting on nfact";
static char* fla_front_str            = "FLA_LU_nopiv_i";
static char* fla_blk_ext_str          = "external";
static char* pc_str[NUM_PARAM_COMBOS] = { "" };
static test_thresh_t thresh           = { 1e+01, 1e-01,   // warn, pass for s
                                          1e-07, 1e-10,   // warn, pass for d
                                          1e+01, 1e-01,   // warn, pass for c
                                          1e-07, 1e-10 }; // warn, pass for z

// Local prototypes.
void libfla_test_lu_nopiv_i_experiment( test_params_t params,
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
void libfla_test_lu_nopiv_i_impl( int         impl,
                                  FLA_Obj     A );
void FLA_GETRFNPI(  integer m,
		    integer n,
		    integer nfact,
		    FLA_Obj A_save,
		    FLA_Obj A,
		    integer lda,
		    FLA_Datatype datatype,
		    unsigned int n_repeats,
		    double* time_min_ );


void libfla_test_lu_nopiv_i( FILE* output_stream, test_params_t params, test_op_t op )
{
	libfla_test_output_info( "--- %s ---\n", op_str );
	libfla_test_output_info( "\n" );
        //op.fla_blk_ext = 1;
	if ( op.fla_blk_ext == ENABLE )
	{
		// libfla_test_output_info( "%s() blocked external variants...\n", fla_front_str );
		// libfla_test_output_info( "\n" );
		libfla_test_op_driver( fla_front_str, fla_blk_ext_str,
		                       FIRST_VARIANT, LAST_VARIANT,
		                       NUM_PARAM_COMBOS, pc_str,
		                       NUM_MATRIX_ARGS,
		                       FLA_TEST_FLAT_BLK_EXT,
		                       params, thresh, libfla_test_lu_nopiv_i_experiment );
	}
}



void libfla_test_lu_nopiv_i_experiment( test_params_t params,
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
        double time_min   = 1e9;
        double time;
	uinteger i;
	uinteger m, n, nfact;
	uinteger lda;
	integer   m_input    = -1;
	integer   n_input    = -1;
	FLA_Obj      A, x, b, norm;
	FLA_Obj      A_save;
	FLA_Obj      A_test;

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur / abs(m_input);
	else               m = p_cur;
	if ( n_input < 0 ) n = p_cur / abs(n_input);
	else               n = p_cur;

        // Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, n, &A );

        // Initialize the test matrices.
	FLA_Random_matrix( A );

        // Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_save );

        // Create vectors to form a linear system.
	FLA_Obj_create( datatype, m, 1, 0, 0, &x );
	FLA_Obj_create( datatype, m, 1, 0, 0, &b );

	// Create a real scalar object to hold the norm of A.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Create a random right-hand side vector.
	FLA_Random_matrix( b );

	    A_test = A;
        //For -2 , it calculates for all and prints min time result
      /*  if( params.p_nfact == -2 )
    {

        //TODO
		//for now

        for( int nfact=0; nfact <=(fla_min(m,n)); nfact++ )
		{

			 if ( impl == FLA_TEST_FLAT_BLK_EXT )
	        {
				time_min = 1e9;
		        lda     = FLA_Obj_col_stride( A );
		        A_test = A;


		        FLA_GETRFNPI( m, n, nfact, A_save, A_test, lda, datatype, n_repeats, &time_min );

	        }
        }

    }
*/
            if( params.p_nfact == -1 || params.p_nfact == -2 )
	    {
            nfact = rand() % fla_min(m,n);
        }
            else
		{
		nfact=params.p_nfact;
		}
        // Invoke FLA_GETRFNPI for external getrfnpi call
            if ( impl == FLA_TEST_FLAT_BLK_EXT )
	   {
		    lda     = FLA_Obj_col_stride( A );
		    A_test = A;
		    FLA_GETRFNPI( m, n, nfact, A_save, A_test, lda, datatype, n_repeats, &time_min );
	    }

	// Compute the performance of the best experiment repeat.
	*t = time_min;
  *perf = 2.0 / 3.0 * m * m * m / time_min / FLOPS_PER_UNIT_PERF;
	if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

        FLA_Obj ATL,ATR,ABL,ABR,AL,AR,L12,U12,AT,AB,A22;

        /* Partition to get (L1 / L2) & (U1 | U2) */

        /* First get L12 = (L1 / L2) */
        FLA_Part_1x2( A_test, &AL, &AR, nfact, FLA_LEFT );
        FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, AL, &L12 );
        FLA_Set_diag( FLA_ONE, L12 );

        /* Zero out the upper triangular part of L12 */
        for( i = 1; i < nfact; i++ )
        {
         FLA_Set_offdiag( i, FLA_ZERO, L12 );
        }

        /* Get U12 = (U1 | U2) */
        FLA_Part_2x1( A_test, &AT,
                              &AB, nfact, FLA_TOP );
        FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, AT, &U12 );

        /* Zero out the lower triangular part of U12 */
        for( i = 1; i < nfact; i++ )
        {
        FLA_Set_offdiag( -i, FLA_ZERO, U12 );
        }

        /* Create A22 by copying A_test and initialize the TL, TR, BL parts to zero */
        FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A_test, &A22 );
        FLA_Part_2x2( A22, &ATL, &ATR,
                           &ABL, &ABR, nfact, nfact, FLA_TL );

        /* Initialize TL, TR, BL parts of A22 to zero */
        FLA_Set( FLA_ZERO, ATL );
        FLA_Set( FLA_ZERO, ATR );
        FLA_Set( FLA_ZERO, ABL );

        /* Reconstruct A in A22 by A22=L12*U12+A22 */
        FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE, FLA_ONE,
                  L12, U12, FLA_ONE, A22 );

        //Calculate norm
        FLA_Axpy_external(FLA_MINUS_ONE,A22,A_save);
        FLA_Norm1(A_save,norm);
        FLA_Obj_extract_real_scalar( norm, residual );
        /* Free the created objects */
        FLA_Obj_free( &L12 );
        FLA_Obj_free( &U12 );
        FLA_Obj_free( &A22 );

        // Free the supporting flat objects.
         FLA_Obj_free( &x );
         FLA_Obj_free( &b );
         FLA_Obj_free( &norm );
         FLA_Obj_free( &A_save );

        // Free the flat test matrices.
        FLA_Obj_free( &A );
}

/*
 * FLA_GETRFNPI calls getrfnpi
 * */
void FLA_GETRFNPI(  integer m,
		    integer n,
		    integer nfact,
		    FLA_Obj A_save,
		    FLA_Obj A,
		    integer lda,
		    FLA_Datatype datatype,
		    unsigned int n_repeats,
		    double* time_min_ )
{
	integer      info;
	unsigned int i;
	double       time;
	double       time_min   = 1e9;

	switch( datatype )
	{
		case FLA_FLOAT:
		{
			for ( i = 0; i < n_repeats; ++i )
			{
			   FLA_Copy_external( A_save, A );
			   float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );

			   time = FLA_Clock();

			   sgetrfnpi_(&m, &n, &nfact, buff_A, &lda, &info);
			   time = FLA_Clock() - time;
			   time_min = fla_min( time_min, time );
			}
			break;
		}
		case FLA_DOUBLE:
		{
			for ( i = 0; i < n_repeats; ++i )
			{
			   FLA_Copy_external( A_save, A );
			   double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );

			   time = FLA_Clock();

			   dgetrfnpi_(&m, &n, &nfact, buff_A, &lda, &info);
			   time = FLA_Clock() - time;
			   time_min = fla_min( time_min, time );
			}
			break;
		}
		case FLA_COMPLEX:
		{
			for ( i = 0; i < n_repeats; ++i )
			{
			   FLA_Copy_external( A_save, A );
			   scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );

			   time = FLA_Clock();

			   cgetrfnpi_(&m, &n, &nfact,buff_A, &lda, &info);
			   time = FLA_Clock() - time;
			   time_min = fla_min( time_min, time );
			}
			break;
		}
		case FLA_DOUBLE_COMPLEX:
		{
			for ( i = 0; i < n_repeats; ++i )
			{
			   FLA_Copy_external( A_save, A );
			   dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

			   time = FLA_Clock();

			   zgetrfnpi_(&m, &n, &nfact, buff_A, &lda, &info);
			   time = FLA_Clock() - time;
			   time_min = fla_min( time_min, time );
			}
			break;
		}
	}
	*time_min_ = time_min;
}
