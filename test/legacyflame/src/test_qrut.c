/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"
#include "test_libflame.h"

#define NUM_PARAM_COMBOS 1
#define NUM_MATRIX_ARGS  2
#define FIRST_VARIANT    1
#define LAST_VARIANT     2

// Static variables.
static char* op_str                   = "QR factorization via UT transform";
static char* flash_front_str          = "FLASH_QR_UT";
static char* fla_front_str            = "FLA_QR_UT";
static char* fla_unb_var_str          = "unb_var";
static char* fla_opt_var_str          = "opt_var";
static char* fla_blk_var_str          = "blk_var";
static char* fla_blk_ext_str          = "external";
static char* pc_str[NUM_PARAM_COMBOS] = { "" };
static test_thresh_t thresh           = { 1e-02, 1e-03,   // warn, pass for s
                                          1e-11, 1e-12,   // warn, pass for d
                                          1e-02, 1e-03,   // warn, pass for c
                                          1e-11, 1e-12 }; // warn, pass for z

static fla_apqut_t*     apqut_cntl_blk;
static fla_qrut_t*      qrut_cntl_opt;
static fla_qrut_t*      qrut_cntl_unb;
static fla_qrut_t*      qrut_cntl_blk;
static fla_qrut_t*      qrut_cntl_blk_sub;
static fla_blocksize_t* qrut_cntl_bsize;

// Local prototypes.
void libfla_test_qrut_experiment( test_params_t params,
                                  unsigned int  var,
                                  char*         sc_str,
                                  FLA_Datatype  datatype,
                                  uinteger      p,
                                  unsigned int  pci,
                                  unsigned int  n_repeats,
                                  signed int    impl,
                                  double*       perf,
                                  double*       t,
                                  double*       residual );
void libfla_test_qrut_impl( int     impl,
                            FLA_Obj A,
                            FLA_Obj T );
void libfla_test_qrut_cntl_create( unsigned int var,
                                   dim_t        b_alg_flat );
void libfla_test_qrut_cntl_free( void );
void FLA_GEQRF( integer m,
                integer n,
                FLA_Obj A_save,
                FLA_Obj A,
                FLA_Obj T_obj,
                FLA_Datatype datatype,
                unsigned int n_repeats,
                double* time_min_ );


void libfla_test_qrut( FILE* output_stream, test_params_t params, test_op_t op )
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
		                       params, thresh, libfla_test_qrut_experiment );
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
		                       params, thresh, libfla_test_qrut_experiment );
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
		                       params, thresh, libfla_test_qrut_experiment );
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
		                       params, thresh, libfla_test_qrut_experiment );
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
		                       params, thresh, libfla_test_qrut_experiment );
	}
        if ( op.fla_blk_ext == ENABLE )
        {
                libfla_test_op_driver( fla_front_str, fla_blk_ext_str,
                                       FIRST_VARIANT, LAST_VARIANT,
                                       NUM_PARAM_COMBOS, pc_str,
                                       NUM_MATRIX_ARGS,
                                       FLA_TEST_FLAT_BLK_EXT,
                                       params, thresh, libfla_test_qrut_experiment );
        }
}



void libfla_test_qrut_experiment( test_params_t params,
                                  unsigned int  var,
                                  char*         sc_str,
                                  FLA_Datatype  datatype,
                                  uinteger      p_cur,
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
	uinteger     m, n;
	uinteger     min_m_n;
	integer      m_input    = -2;
	integer      n_input    = -1;
	FLA_Obj      A, T, x, b, y, norm, b_copy;
	FLA_Obj      A_save;
	FLA_Obj      A_test, T_test, x_test, b_test;

	// Determine the dimensions.
	if ( m_input < 0 ) m = p_cur * abs(m_input);
	else               m = p_cur;
	if ( n_input < 0 ) n = p_cur * abs(n_input);
	else               n = p_cur;

	// Compute the minimum dimension.
	min_m_n = fla_min( m, n );

	// Create the matrices for the current operation.
	libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[0], m, n, &A );

	if ( impl == FLA_TEST_FLAT_FRONT_END ||
	     ( impl == FLA_TEST_FLAT_BLK_EXT ) ||
	     ( impl == FLA_TEST_FLAT_BLK_VAR && var == 1 ) )
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], b_alg_flat, min_m_n, &T );
	else if ( var == 2 )
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], min_m_n, min_m_n, &T );
	else
		libfla_test_obj_create( datatype, FLA_NO_TRANSPOSE, sc_str[1], 1, min_m_n, &T );

	// Initialize the test matrices.
	FLA_Random_matrix( A );

	// Save the original object contents in a temporary object.
	FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, A, &A_save );

	// Create vectors to form a linear system.
	FLA_Obj_create( datatype, n, 1, 0, 0, &x );
	FLA_Obj_create( datatype, m, 1, 0, 0, &b );
	FLA_Obj_create( datatype, n, 1, 0, 0, &y );

	// Create a real scalar object to hold the norm of A.
	FLA_Obj_create( FLA_Obj_datatype_proj_to_real( A ), 1, 1, 0, 0, &norm );

	// Create a random right-hand side vector.
	FLA_Random_matrix( b );
	
	// Use hierarchical matrices if we're testing the FLASH front-end.
	if ( impl == FLA_TEST_HIER_FRONT_END )
	{
		FLASH_QR_UT_create_hier_matrices( A, 1, &b_flash, &A_test, &T_test );
		FLASH_Obj_create_hier_copy_of_flat( b, 1, &b_flash, &b_test );
		FLASH_Obj_create_hier_copy_of_flat( x, 1, &b_flash, &x_test );
	}
	else
	{
		A_test = A;
		T_test = T;
	}

	// Create a control tree for the individual variants.
	if ( impl == FLA_TEST_FLAT_UNB_VAR ||
	     impl == FLA_TEST_FLAT_OPT_VAR ||
	     impl == FLA_TEST_FLAT_BLK_VAR )
		libfla_test_qrut_cntl_create( var, b_alg_flat );
    
    // Invoke FLA_GETRF for external geqrf call
    if ( impl != FLA_TEST_FLAT_BLK_EXT )
    {
      
        // Repeat the experiment n_repeats times and record results.
		for ( i = 0; i < n_repeats; ++i )
		{
			if ( impl == FLA_TEST_HIER_FRONT_END )
				FLASH_Obj_hierarchify( A_save, A_test );
			else
				FLA_Copy_external( A_save, A_test );
		
			time = FLA_Clock();

			libfla_test_qrut_impl( impl, A_test, T_test );
		
			time = FLA_Clock() - time;
			time_min = fla_min( time_min, time );
		}
        
        // Compute the performance of the best experiment repeat.
        *t = time_min;
        *perf = (         2.0   * m * n * n - 
                ( 2.0 / 3.0 ) * n * n * n ) / time_min / FLOPS_PER_UNIT_PERF;
        if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;
	
        // Perform a linear solve with the result.
        if ( impl == FLA_TEST_HIER_FRONT_END )
        {
            FLASH_QR_UT_solve( A_test, T_test, b_test, x_test );
            FLASH_Obj_flatten( x_test, x );
        }
        else
        {
            FLA_QR_UT_solve( A_test, T_test, b, x );
        }

        // Free the hierarchical matrices if we're testing the FLASH front-end.
        if ( impl == FLA_TEST_HIER_FRONT_END )
        {
            FLASH_Obj_free( &A_test );
            FLASH_Obj_free( &T_test );
            FLASH_Obj_free( &b_test );
            FLASH_Obj_free( &x_test );
        }
        
        // Free the control trees if we're testing the variants.
        if ( impl == FLA_TEST_FLAT_UNB_VAR ||
             impl == FLA_TEST_FLAT_OPT_VAR ||
             impl == FLA_TEST_FLAT_BLK_VAR )
            libfla_test_qrut_cntl_free();

        // Compute the residual.
        FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A_save, x, FLA_MINUS_ONE, b );
        FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A_save, b, FLA_ZERO, y );
        FLA_Nrm2_external( y, norm );
        FLA_Obj_extract_real_scalar( norm, residual );
	}
    else
    {
	// Repeat the experiment n_repeats times and record results.
	for ( i = 0; i < n_repeats; ++i )
	{
            FLA_Copy_external( A_save, A_test );
			
            
            FLA_GEQRF( m, n, A_save, A_test, T_test, datatype, n_repeats, &time_min );
			
        }

        // Compute the performance of the best experiment repeat.
        *t = time_min;
        *perf = (         2.0   * m * n * n - 
                ( 2.0 / 3.0 ) * n * n * n ) / time_min / FLOPS_PER_UNIT_PERF;
        if ( FLA_Obj_is_complex( A ) ) *perf *= 4.0;

        integer cs_A;
        integer cs_b, lwork, tinfo;
        integer mb, nb;
        integer min_m_n;
        FLA_Obj qbt;

	min_m_n = fla_min( m, n );

        switch( datatype )
        {
            case FLA_FLOAT:
            {
                float* twork, wsize;
                FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &qbt );
                
                float* out_T   = FLA_FLOAT_PTR( T_test );
                float* mat_A   = FLA_FLOAT_PTR( A_test );
                float* ptr_b   = FLA_FLOAT_PTR( b );
                float* ptr_qbt = FLA_FLOAT_PTR( qbt );

                cs_A = FLA_Obj_col_stride( A_test );


                mb = b.m;
                nb = b.n;
                cs_b = FLA_Obj_col_stride( b );

                // Save the original object contents in a temporary object.
                FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, b, &b_copy );

                // b' = QT*b
                // Get size of work buffer
                lwork = -1;
                sormqr_( "Left", "Transpose", &mb, &nb, &min_m_n, mat_A, &cs_A, out_T,
                          ptr_b, &cs_b, &wsize, &lwork, &tinfo );

                // allocate work buffer and calculate b'
                lwork = (integer) wsize;
                twork = malloc( lwork * sizeof( float ) );
                sormqr_( "Left", "Transpose", &mb, &nb, &min_m_n, mat_A, &cs_A, out_T,
                          ptr_b, &cs_b, twork, &lwork, &tinfo );
                free( twork );

                // Triangular Solve to find x in Ax=b
                strtrs_( "Upper", "No-Tran", "Non-unit", &min_m_n, &nb, mat_A, &cs_A,
                          ptr_b, &cs_b, &tinfo );

                for( i = 0; i < min_m_n; i++ )
                {
                    ptr_qbt[i] = ptr_b[i];
                }
                break;
            }

            case FLA_DOUBLE:
            {
                double* twork, wsize;
                FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &qbt );
                
                double* out_T   = FLA_DOUBLE_PTR( T_test );
                double* mat_A   = FLA_DOUBLE_PTR( A_test );
                double* ptr_b   = FLA_DOUBLE_PTR( b );
                double* ptr_qbt = FLA_DOUBLE_PTR( qbt );

                cs_A = FLA_Obj_col_stride( A_test );

                mb = b.m;
                nb = b.n;
                cs_b = FLA_Obj_col_stride( b );

                // Save the original object contents in a temporary object.
                FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, b, &b_copy );

                // b' = QT*b
                // Get size of work buffer
                lwork = -1;
                dormqr_( "Left", "Transpose", &mb, &nb, &min_m_n, mat_A, &cs_A, out_T,
                          ptr_b, &cs_b, &wsize, &lwork, &tinfo );

                // allocate work buffer and calculate b'
                lwork = (integer) wsize;
                twork = malloc( lwork * sizeof( double ) );
                dormqr_( "Left", "Transpose", &mb, &nb, &min_m_n, mat_A, &cs_A, out_T,
                          ptr_b, &cs_b, twork, &lwork, &tinfo );
                free( twork );

                // Triangular Solve to find x in Ax=b
                dtrtrs_( "Upper", "No-Tran", "Non-unit", &min_m_n, &nb, mat_A, &cs_A,
                          ptr_b, &cs_b, &tinfo );

                for( i = 0; i < min_m_n; i++ )
                {
                    ptr_qbt[i] = ptr_b[i];
                }
                break;
            }
            case FLA_COMPLEX:
            {
                scomplex* twork, wsize;
                FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &qbt );
                scomplex* out_T   = FLA_COMPLEX_PTR( T_test );
                scomplex* mat_A   = FLA_COMPLEX_PTR( A_test );
                scomplex* ptr_b   = FLA_COMPLEX_PTR( b );
                scomplex* ptr_qbt = FLA_COMPLEX_PTR( qbt );

                cs_A = FLA_Obj_col_stride( A_test );

                mb = b.m;
                nb = b.n;
                cs_b = FLA_Obj_col_stride( b );

                // Save the original object contents in a temporary object.
                FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, b, &b_copy );

                // b' = QT*b
                // Get size of work buffer
                lwork = -1;
                cunmqr_( "Left", "Conjugate Transpose", &mb, &nb, &min_m_n, mat_A, &cs_A, out_T,
                          ptr_b, &cs_b, &wsize, &lwork, &tinfo );

                // allocate work buffer and calculate b'
                lwork = (integer) wsize.real;
                twork = malloc( lwork * sizeof( scomplex ) );
                cunmqr_( "Left", "Conjugate Transpose", &mb, &nb, &min_m_n, mat_A, &cs_A, out_T,
                          ptr_b, &cs_b, twork, &lwork, &tinfo );
                free( twork );
                // Triangular Solve to find x in Ax=b
                ctrtrs_( "Upper", "No-Tran", "Non-unit", &min_m_n, &nb, mat_A, &cs_A,
                          ptr_b, &cs_b, &tinfo );

                for( i = 0; i < min_m_n; i++ )
                {
                    ptr_qbt[i].real = ptr_b[i].real;
		            ptr_qbt[i].imag = ptr_b[i].imag;
                }
                break;
            }
            case FLA_DOUBLE_COMPLEX:
            {
                dcomplex* twork, wsize;
                FLA_Obj_create( datatype, min_m_n, 1, 0, 0, &qbt );
                dcomplex* out_T   = FLA_DOUBLE_COMPLEX_PTR( T_test );
                dcomplex* mat_A   = FLA_DOUBLE_COMPLEX_PTR( A_test );
                dcomplex* ptr_b   = FLA_DOUBLE_COMPLEX_PTR( b );
                dcomplex* ptr_qbt = FLA_DOUBLE_COMPLEX_PTR( qbt );

                cs_A = FLA_Obj_col_stride( A_test );

                mb = b.m;
                nb = b.n;
                cs_b = FLA_Obj_col_stride( b );

                // Save the original object contents in a temporary object.
                FLA_Obj_create_copy_of( FLA_NO_TRANSPOSE, b, &b_copy );

                // b' = QT*b
                // Get size of work buffer
                lwork = -1;
                zunmqr_( "Left", "Conjugate Transpose", &mb, &nb, &min_m_n, mat_A, &cs_A, out_T,
                          ptr_b, &cs_b, &wsize, &lwork, &tinfo );
                // allocate work buffer and calculate b'
                lwork = (integer) wsize.real;
                twork = malloc( lwork * sizeof( dcomplex ) );
                zunmqr_( "Left", "Conjugate Transpose", &mb, &nb, &min_m_n, mat_A, &cs_A, out_T,
                          ptr_b, &cs_b, twork, &lwork, &tinfo );
                free( twork );
                // Triangular Solve to find x in Ax=b
                ztrtrs_( "Upper", "No-Tran", "Non-unit", &min_m_n, &nb, mat_A, &cs_A,
                          ptr_b, &cs_b, &tinfo );

                for( i = 0; i < min_m_n; i++ )
                {
                    ptr_qbt[i].real = ptr_b[i].real;
        		    ptr_qbt[i].imag = ptr_b[i].imag;
                }
                break;
            }
        }
        
        // Compute the residual.
        FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_ONE, A_save, qbt, FLA_MINUS_ONE, b_copy );
        FLA_Gemv_external( FLA_CONJ_TRANSPOSE, FLA_ONE, A_save, b_copy, FLA_ZERO, y );
        FLA_Nrm2_external( y, norm );
        FLA_Obj_extract_real_scalar( norm, residual );

        FLA_Obj_free( &b_copy );
        FLA_Obj_free( &qbt );
    }

	// Free the supporting flat objects.
	FLA_Obj_free( &x );
	FLA_Obj_free( &b );
	FLA_Obj_free( &y );
	FLA_Obj_free( &norm );
	FLA_Obj_free( &A_save );

	// Free the flat test matrices.
	FLA_Obj_free( &A );
	FLA_Obj_free( &T );
}



extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_axpyt_t* fla_axpyt_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_copyt_t* fla_copyt_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_gemm_t*  fla_gemm_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_trmm_t*  fla_trmm_cntl_blas;
extern LIBFLAME_IMPORT TLS_CLASS_SPEC fla_trsm_t*  fla_trsm_cntl_blas;

void libfla_test_qrut_cntl_create( unsigned int var,
                                   dim_t        b_alg_flat )
{
	int var_unb  = FLA_UNB_VAR_OFFSET + var;
	int var_opt  = FLA_OPT_VAR_OFFSET + var;
	int var_blk  = FLA_BLK_VAR_OFFSET + var;
	int var_opt2 = FLA_OPT_VAR_OFFSET + 2;

	qrut_cntl_bsize = FLA_Blocksize_create( b_alg_flat, b_alg_flat, b_alg_flat, b_alg_flat );

	apqut_cntl_blk  = FLA_Cntl_apqut_obj_create( FLA_FLAT,
	                                             FLA_BLOCKED_VARIANT1,
	                                             qrut_cntl_bsize,
	                                             NULL,
	                                             fla_trmm_cntl_blas,
	                                             fla_trmm_cntl_blas,
	                                             fla_gemm_cntl_blas,
	                                             fla_gemm_cntl_blas,
	                                             fla_trsm_cntl_blas,
	                                             fla_copyt_cntl_blas,
	                                             fla_axpyt_cntl_blas );

	qrut_cntl_unb   = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                                var_unb,
                                                NULL,
                                                NULL,
                                                NULL );

	qrut_cntl_opt   = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                                var_opt,
                                                NULL,
                                                NULL,
                                                NULL );

	// Use unblocked variant 2 for blocked variant 1 subproblem.
	qrut_cntl_blk_sub = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                                  var_opt2,
                                                  NULL,
                                                  NULL,
                                                  NULL );

	qrut_cntl_blk   = FLA_Cntl_qrut_obj_create( FLA_FLAT,
                                                var_blk,
                                                qrut_cntl_bsize,
                                                qrut_cntl_blk_sub,
                                                apqut_cntl_blk );
}



void libfla_test_qrut_cntl_free( void )
{
	FLA_Blocksize_free( qrut_cntl_bsize );

	FLA_Cntl_obj_free( apqut_cntl_blk );
	FLA_Cntl_obj_free( qrut_cntl_unb );
	FLA_Cntl_obj_free( qrut_cntl_opt );
	FLA_Cntl_obj_free( qrut_cntl_blk );
	FLA_Cntl_obj_free( qrut_cntl_blk_sub );
}



void libfla_test_qrut_impl( int     impl,
                            FLA_Obj A,
                            FLA_Obj T )
{
	switch ( impl )
	{
		case FLA_TEST_HIER_FRONT_END:
		FLASH_QR_UT( A, T );
		break;

		case FLA_TEST_FLAT_FRONT_END:
		FLA_QR_UT( A, T );
		break;

		case FLA_TEST_FLAT_UNB_VAR:
		FLA_QR_UT_internal( A, T, qrut_cntl_unb );
		break;

		case FLA_TEST_FLAT_OPT_VAR:
		FLA_QR_UT_internal( A, T, qrut_cntl_opt );
		break;

		case FLA_TEST_FLAT_BLK_VAR:
		FLA_QR_UT_internal( A, T, qrut_cntl_blk );
		break;

		default:
		libfla_test_output_error( "Invalid implementation type.\n" );
	}
}

/*
 *  FLA_GEQRF calls LAPACK interface of
 *  QR Factorization  - geqrf
 *    *  */
void FLA_GEQRF( integer m,
                integer n,
                FLA_Obj A_save,
                FLA_Obj A,
                FLA_Obj T_obj,
                FLA_Datatype datatype,
                unsigned int n_repeats,
                double* time_min_ )
{
    integer          info;
    unsigned int i, j, k;
    double       time;
    double       time_min   = 1e9;
    integer lda;
    integer lwork;
	integer min_m_n = fla_min(m, n);
    integer rs_A, cs_A;

    lda     = (integer)FLA_Obj_col_stride( A );
    
    switch( datatype )
    {
        case FLA_FLOAT:
            {
                float* buff_T, *out_T;
                float* work;
                float  worksize;
                float* mat_A = ( float * ) FLA_FLOAT_PTR( A );
                float* buff_A = (float *) malloc( m * n * sizeof( float ) );
                float* aptr = buff_A;
			
        		buff_T    = ( float * )malloc( min_m_n * sizeof( float ) );
        		out_T     = FLA_FLOAT_PTR( T_obj );
                
                lwork = -1;
                sgeqrf_(&m, &n, buff_A, &lda, buff_T, &worksize, &lwork, &info);
                
                lwork = (int)worksize;
                work = (float * )malloc(lwork * sizeof(float));
 
                rs_A     = FLA_Obj_row_stride( A );
                cs_A     = FLA_Obj_col_stride( A );

                for ( i = 0; i < n_repeats; ++i )
                {
                    FLA_Copy_external( A_save, A );
                    
                    for( k = 0; k < m; k++ )
                    {
                        for( j = 0; j < n; j++ )
                        {
                            aptr[m*j+k] = mat_A[cs_A * j + rs_A * k];
                        }
                    }
                    time = FLA_Clock();
                    
                    sgeqrf_(&m, &n, buff_A, &lda, buff_T, work, &lwork, &info);

                    time = FLA_Clock() - time;
                    time_min = fla_min( time_min, time );
                }

                // copy tau values to output FLA obj
                for( i = 0; i < min_m_n; i++ )
                {
                    out_T[i] = buff_T[i];
                }

                // copy output matrix
                for( k = 0; k < m; k++ )
                {
                    for( j = 0; j < n; j++ )
                    {
                        mat_A[cs_A * j + rs_A * k] = aptr[m*j+k];
                    }
                }

                free(buff_T);
                free(work);
                free(buff_A);
                break;
            }
        case FLA_DOUBLE:
            {
                double* buff_T, *out_T;
                double* work;
                double  worksize;
                double* mat_A = ( double * ) FLA_DOUBLE_PTR( A );
                double* buff_A = ( double *) malloc( m * n * sizeof( double ) );
                double* aptr = buff_A;
			
        		buff_T    = ( double * )malloc( min_m_n * sizeof( double ) );
        		out_T     = FLA_DOUBLE_PTR( T_obj );
                
                lwork = -1;
                dgeqrf_(&m, &n, buff_A, &lda, buff_T, &worksize, &lwork, &info);
                
                lwork = (int)worksize;
                work = (double * )malloc(lwork * sizeof(double));

                rs_A     = FLA_Obj_row_stride( A );
                cs_A     = FLA_Obj_col_stride( A );
                
                for ( i = 0; i < n_repeats; ++i )
                {
                    FLA_Copy_external( A_save, A );
                    
                    for( k = 0; k < m; k++ )
                    {
                        for( j = 0; j < n; j++ )
                        {
                            aptr[m*j+k] = mat_A[cs_A * j + rs_A * k];
                        }
                    }
                    time = FLA_Clock();
                    
                    dgeqrf_(&m, &n, buff_A, &lda, buff_T, work, &lwork, &info);

                    time = FLA_Clock() - time;
                    time_min = fla_min( time_min, time );
                }

                // copy tau values to output FLA obj
                for( i = 0; i < min_m_n; i++ )
                {
                    out_T[i] = buff_T[i];
                }

                // copy output matrix
                for( k = 0; k < m; k++ )
                {
                    for( j = 0; j < n; j++ )
                    {
                        mat_A[cs_A * j + rs_A * k] = aptr[m*j+k];
                    }
                }

                free(buff_T);
                free(work);
                free(buff_A);
                break;
            }
        case FLA_COMPLEX:
                 {
                  scomplex* buff_T, *out_T;
                  scomplex* work;
                  scomplex worksize;
                  scomplex* mat_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
                  scomplex* buff_A = ( scomplex *) malloc( m * n * sizeof( scomplex ) );
                  scomplex* aptr = buff_A;
                  buff_T     = ( scomplex * )malloc(min_m_n * sizeof(scomplex));
                  out_T     = FLA_COMPLEX_PTR( T_obj );
                  lwork = -1;
                  cgeqrf_(&m, &n, buff_A, &lda, buff_T, &worksize, &lwork, &info);

                  lwork = (integer)worksize.real;
                  work = (scomplex * )malloc(lwork * sizeof(scomplex));
                  rs_A     = FLA_Obj_row_stride( A );
                  cs_A     = FLA_Obj_col_stride( A );

                 for ( i = 0; i < n_repeats; ++i )
                 {
                           FLA_Copy_external( A_save, A );
                           for( k = 0; k < m; k++ )
                           {
                              for( j = 0; j < n; j++ )
                             {
                                aptr[m*j+k].real = mat_A[cs_A * j + rs_A * k].real;
				aptr[m*j+k].imag = mat_A[cs_A * j + rs_A * k].imag;

                             }
                           }                    

                           time = FLA_Clock();

                           cgeqrf_(&m, &n, buff_A, &lda, buff_T, work, &lwork, &info);

                           time = FLA_Clock() - time;
                           time_min = fla_min( time_min, time );
                 }
                 // copy tau values to output FLA obj
                 for( i = 0; i < min_m_n; i++ )
                 {
                       out_T[i].real = buff_T[i].real;
		       out_T[i].imag = buff_T[i].imag;

                 }

                    // copy output matrix
                 for( k = 0; k < m; k++ )
                 {
                     for( j = 0; j < n; j++ )
                     {
                           mat_A[cs_A * j + rs_A * k].real = aptr[m*j+k].real;
			   mat_A[cs_A * j + rs_A * k].imag = aptr[m*j+k].imag;


                     }
                 }

                  free(buff_T);
                  free(work);
                  free(buff_A);
                  break;
                }
                case FLA_DOUBLE_COMPLEX:
                {
                dcomplex* buff_T, *out_T;
                dcomplex* work;
                dcomplex worksize;
                dcomplex* mat_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
                dcomplex* buff_A = ( dcomplex *) malloc( m * n * sizeof( dcomplex ) );
                dcomplex* aptr = buff_A;
                buff_T     = ( dcomplex * )malloc(min_m_n * sizeof(dcomplex));
                out_T     = FLA_DOUBLE_COMPLEX_PTR( T_obj );
                lwork = -1;
                zgeqrf_(&m, &n, buff_A, &lda, buff_T, &worksize, &lwork, &info);
                lwork = (integer)worksize.real;
                work = (dcomplex * )malloc(lwork * sizeof(dcomplex));

                rs_A     = FLA_Obj_row_stride( A );
                cs_A     = FLA_Obj_col_stride( A );

                for ( i = 0; i < n_repeats; ++i )
                {
                    FLA_Copy_external( A_save, A );
                    for( k = 0; k < m; k++ )
                    {
                       for( j = 0; j < n; j++ )
                       {
                                aptr[m*j+k].real = mat_A[cs_A * j + rs_A * k].real;
				aptr[m*j+k].imag = mat_A[cs_A * j + rs_A * k].imag;

                       }
                     }                    

                     time = FLA_Clock();

                     zgeqrf_(&m, &n, buff_A, &lda, buff_T, work, &lwork, &info);

                     time = FLA_Clock() - time;
                     time_min = fla_min( time_min, time );
                 }
                 // copy tau values to output FLA obj
                 for( i = 0; i < min_m_n; i++ )
                 {
                       out_T[i].real = buff_T[i].real;
		               out_T[i].imag = buff_T[i].imag;
                 }

                    // copy output matrix
                 for( k = 0; k < m; k++ )
                 {
                     for( j = 0; j < n; j++ )
                     {
                            mat_A[cs_A * j + rs_A * k].real = aptr[m*j+k].real;
                            mat_A[cs_A * j + rs_A * k].imag = aptr[m*j+k].imag;
                     }
                 }

                 free(buff_T);
                 free(work);
                 free(buff_A);
                 break;
                }
        }
        *time_min_ = time_min;
}

