/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#include "test_libflame.h"

// Operation modules.
#include "test_gemm.h"
#include "test_hemm.h"
#include "test_herk.h"
#include "test_her2k.h"
#include "test_symm.h"
#include "test_syrk.h"
#include "test_syr2k.h"
#include "test_trmm.h"
#include "test_trsm.h"
#include "test_chol.h"
#include "test_lu_nopiv_i.h"
#include "test_lu_nopiv.h"
#include "test_lu_piv.h"
#include "test_lu_incpiv.h"
#include "test_qrut.h"
#include "test_qrutinc.h"
#include "test_lqut.h"
#include "test_apqut.h"
#include "test_apqutinc.h"
#include "test_caqrutinc.h"
#include "test_apcaqutinc.h"
#include "test_uddateut.h"
#include "test_uddateutinc.h"
#include "test_apqudut.h"
#include "test_apqudutinc.h"
#include "test_hessut.h"
#include "test_tridiagut.h"
#include "test_bidiagut.h"
#include "test_eig_gest.h"
#include "test_trinv.h"
#include "test_spdinv.h"
#include "test_sylv.h"
#include "test_lyap.h"
#include "test_ldlt2_nopiv_ps.h"
#include "test_ldltx_nopiv_ps.h"


// Global variables.
char libfla_test_binary_name[ MAX_BINARY_NAME_LENGTH + 1 ];
char libfla_test_pass_string[ MAX_PASS_STRING_LENGTH + 1 ];
char libfla_test_warn_string[ MAX_PASS_STRING_LENGTH + 1 ];
char libfla_test_fail_string[ MAX_PASS_STRING_LENGTH + 1 ];
char libfla_test_storage_format_string[ 200 ];

char libfla_test_stor_chars[ NUM_STORAGE_CHARS + 1 ];
void libfla_test_read_tests_for_op_ext( FILE* input_stream, test_op_t* op );
void libfla_test_output_op_struct_ext( char* op_str, test_op_t op );
void libfla_test_read_tests_for_op_fla_ext( FILE* input_stream, test_op_t* op );
void libfla_test_output_op_struct_fla_ext( char* op_str, test_op_t op );

int main( int argc, char** argv )
{
	test_params_t params;
	test_ops_t    ops;

	// Initialize libflame.
	FLA_Init();

	printf(" LibFlame version: %s \n", FLA_Get_AOCL_Version() );

	// Initialize some strings.
	libfla_test_init_strings();

	// Parse the command line parameters.
	libfla_test_parse_command_line( argc, argv );

	// Read the main test suite parameters.
	libfla_test_read_parameter_file( PARAMETERS_FILENAME, &params );

	// Read which operations we're going to test.
	libfla_test_read_operation_file( OPERATIONS_FILENAME, &ops );

	// Test the BLAS level-3 operations.
	libfla_test_blas3_suite( stdout, params, ops );

	// Test the LAPACK-level operations.
	libfla_test_lapack_suite( stdout, params, ops );

	// Finalize libflame.
	FLA_Finalize();

	// Return peacefully.
	return 0;
}


void libfla_test_blas3_suite( FILE* output_stream, test_params_t params, test_ops_t ops )
{
	// Run the individual test modules.

	libfla_test_output_info( "\n" );
	libfla_test_output_info( "--- BLAS level-3 operation tests ---------------------\n" );
	libfla_test_output_info( "\n" );

	// General matrix-matrix multiply.
	libfla_test_gemm( output_stream, params, ops.gemm );

	// Hermitian matrix-matrix multiply.
	libfla_test_hemm( output_stream, params, ops.hemm );

	// Hermitian rank-k update.
	libfla_test_herk( output_stream, params, ops.herk );

	// Hermitian rank-2k update.
	libfla_test_her2k( output_stream, params, ops.her2k );

	// Symmetric matrix-matrix multiply.
	libfla_test_symm( output_stream, params, ops.symm );

	// Symmetric rank-k update.
	libfla_test_syrk( output_stream, params, ops.syrk );

	// Symmetric rank-2k update.
	libfla_test_syr2k( output_stream, params, ops.syr2k );

	// Triangular matrix-matrix multiply.
	libfla_test_trmm( output_stream, params, ops.trmm );

	// Triangular solve with multiple rhs.
	libfla_test_trsm( output_stream, params, ops.trsm );
}



void libfla_test_lapack_suite( FILE* output_stream, test_params_t params, test_ops_t ops )
{
	// Run the individual test modules.

	libfla_test_output_info( "\n" );
	libfla_test_output_info( "--- LAPACK-level operation tests ---------------------\n" );
	libfla_test_output_info( "\n" );

    // Cholesky factorization.
	libfla_test_chol( output_stream, params, ops.chol );

	// LU factorization without pivoting.
	libfla_test_lu_nopiv( output_stream, params, ops.lu_nopiv );

	// LU complete/incomplete factorization based on nfact without pivoting.
	libfla_test_lu_nopiv_i( output_stream, params, ops.lu_nopiv_i );

	// LU factorization with partial pivoting.
	libfla_test_lu_piv( output_stream, params, ops.lu_piv );

	// LU factorization with incremental pivoting.
	libfla_test_lu_incpiv( output_stream, params, ops.lu_incpiv );

	// QR factorization via the UT transform.
	libfla_test_qrut( output_stream, params, ops.qrut );

	// QR factorization via the UT transform (incremental).
	libfla_test_qrutinc( output_stream, params, ops.qrutinc );

	// LQ factorization via the UT transform.
	libfla_test_lqut( output_stream, params, ops.lqut );

	// Apply Q via the UT transform.
	libfla_test_apqut( output_stream, params, ops.apqut );

	// Apply Q via the UT transform (incremental).
	libfla_test_apqutinc( output_stream, params, ops.apqutinc );

	// Communication-avoiding QR factorization via the UT transform (incremental).
	libfla_test_caqrutinc( output_stream, params, ops.caqrutinc );

	// Apply communication-avoiding Q via the UT transform (incremental).
	libfla_test_apcaqutinc( output_stream, params, ops.apcaqutinc );

	// Up/downdating via the UT transform.
	libfla_test_uddateut( output_stream, params, ops.uddateut );

	// Up/downdating via the UT transform (incremental).
	libfla_test_uddateutinc( output_stream, params, ops.uddateutinc );

	// Apply up/downdating Q via the UD UT transform.
	libfla_test_apqudut( output_stream, params, ops.apqudut );

	// Apply up/downdating Q via the UD UT transform (incremental).
	libfla_test_apqudutinc( output_stream, params, ops.apqudutinc );

	// Reduction to upper Hessenberg form via the UT transform.
	libfla_test_hessut( output_stream, params, ops.hessut );

	// Reduction to tridiagonal form via the UT transform.
	libfla_test_tridiagut( output_stream, params, ops.tridiagut );

	// Reduction to bidiagonal form via the UT transform.
	libfla_test_bidiagut( output_stream, params, ops.bidiagut );

	// Reduction of Hermitian-definite eigenproblem to standard form.
	libfla_test_eig_gest( output_stream, params, ops.eig_gest );

	// Triangular matrix inversion.
	libfla_test_trinv( output_stream, params, ops.trinv );

	// Hermitian positive-definite matrix inversion.
	libfla_test_spdinv( output_stream, params, ops.spdinv );

	// Triangular Sylvester equation solve.
	libfla_test_sylv( output_stream, params, ops.sylv );

	// Triangular Lyapunov equation solve.
	libfla_test_lyap( output_stream, params, ops.lyap );

	// LDLT Transform incomplete
	libfla_test_ldlt2_nopiv_ps( &params, ops.ldlt_nopiv_part );
	libfla_test_ldltx_nopiv_ps( &params, ops.ldlt_nopiv_part );

}



void libfla_test_read_operation_file( char* input_filename, test_ops_t* ops )
{
	FILE* input_stream;

	// Attempt to open input file corresponding to input_filename as
	// read-only/binary.
	input_stream = fopen( input_filename, "rb" );

	// Check for success.
	if ( input_stream == NULL )
	{
		libfla_test_output_error( "Failed to open input file %s. Check existence and permissions.\n",
		                          input_filename );
	}

	libfla_test_output_info( "\n" );
	libfla_test_output_info( "--- operations to test -------------------------------\n" );
	libfla_test_output_info( "\n" );

	// Read the operation tests for general matrix-matrix multiply.
	libfla_test_read_tests_for_op_blas3( input_stream, &(ops->gemm) );
	libfla_test_output_op_struct_blas3( "gemm", ops->gemm );

	// Read the operation tests for Hermitian matrix-matrix multiply.
	libfla_test_read_tests_for_op_blas3( input_stream, &(ops->hemm) );
	libfla_test_output_op_struct_blas3( "hemm", ops->hemm );

	// Read the operation tests for Hermitian rank-k update.
	libfla_test_read_tests_for_op_blas3( input_stream, &(ops->herk) );
	libfla_test_output_op_struct_blas3( "herk", ops->herk );

	// Read the operation tests for Hermitian rank-2k update.
	libfla_test_read_tests_for_op_blas3( input_stream, &(ops->her2k) );
	libfla_test_output_op_struct_blas3( "her2k", ops->her2k );

	// Read the operation tests for Symmetric matrix-matrix multiply.
	libfla_test_read_tests_for_op_blas3( input_stream, &(ops->symm) );
	libfla_test_output_op_struct_blas3( "symm", ops->symm );

	// Read the operation tests for Symmetric rank-k update.
	libfla_test_read_tests_for_op_blas3( input_stream, &(ops->syrk) );
	libfla_test_output_op_struct_blas3( "syrk", ops->syrk );

	// Read the operation tests for Symmetric rank-2k update.
	libfla_test_read_tests_for_op_blas3( input_stream, &(ops->syr2k) );
	libfla_test_output_op_struct_blas3( "syr2k", ops->syr2k );

	// Read the operation tests for Triangular matrix-matrix multiply.
	libfla_test_read_tests_for_op_blas3( input_stream, &(ops->trmm) );
	libfla_test_output_op_struct_blas3( "trmm", ops->trmm );

	// Read the operation tests for Triangular solve with multiple rhs.
	libfla_test_read_tests_for_op_blas3( input_stream, &(ops->trsm) );
	libfla_test_output_op_struct_blas3( "trsm", ops->trsm );

	// Read the operation tests for Cholesky factorization.
	libfla_test_read_tests_for_op_ext( input_stream, &(ops->chol) );
	libfla_test_output_op_struct_ext( "chol", ops->chol );

	// Read the operation tests for LU_nopiv factorization.
	libfla_test_read_tests_for_op_ext( input_stream, &(ops->lu_nopiv) );
	libfla_test_output_op_struct_ext( "lu_nopiv", ops->lu_nopiv );

	//Read the operation tests for LU_nopiv_i factorization
	libfla_test_read_tests_for_op_fla_ext( input_stream, &(ops->lu_nopiv_i) );
	libfla_test_output_op_struct_fla_ext( "lu_nopiv_i", ops->lu_nopiv_i );

	// Read the operation tests for LU_piv factorization.
	libfla_test_read_tests_for_op_ext( input_stream, &(ops->lu_piv) );
	libfla_test_output_op_struct_ext( "lu_piv", ops->lu_piv );

	// Read the operation tests for LU_incpiv factorization.
	libfla_test_read_tests_for_op_flash_only( input_stream, &(ops->lu_incpiv) );
	libfla_test_output_op_struct_flash_only( "lu_incpiv", ops->lu_incpiv );

	// Read the operation tests for LDLT_nopiv_part factorization.
	libfla_test_read_tests_for_op_ext( input_stream, &(ops->ldlt_nopiv_part) );
	libfla_test_output_op_struct_ext( "ldlt_nopiv_part", ops->ldlt_nopiv_part );

	// Read the operation tests for QR_UT factorization.
	libfla_test_read_tests_for_op_ext( input_stream, &(ops->qrut) );
	libfla_test_output_op_struct_ext( "qrut", ops->qrut );

	// Read the operation tests for QR_UT_inc factorization.
	libfla_test_read_tests_for_op_flash_only( input_stream, &(ops->qrutinc) );
	libfla_test_output_op_struct_flash_only( "qrutinc", ops->qrutinc );

	// Read the operation tests for LQ_UT factorization.
	libfla_test_read_tests_for_op( input_stream, &(ops->lqut) );
	libfla_test_output_op_struct( "lqut", ops->lqut );

	// Read the operation tests for Apply_Q_UT.
	libfla_test_read_tests_for_op_front_only( input_stream, &(ops->apqut) );
	libfla_test_output_op_struct_front_only( "apqut", ops->apqut );

	// Read the operation tests for Apply_Q_UT_inc.
	libfla_test_read_tests_for_op_flash_only( input_stream, &(ops->apqutinc) );
	libfla_test_output_op_struct_flash_only( "apqutinc", ops->apqutinc );

	// Read the operation tests for CAQR_UT_inc factorization.
	libfla_test_read_tests_for_op_flash_only( input_stream, &(ops->caqrutinc) );
	libfla_test_output_op_struct_flash_only( "caqrutinc", ops->caqrutinc );

	// Read the operation tests for Apply_CAQ_UT_inc.
	libfla_test_read_tests_for_op_flash_only( input_stream, &(ops->apcaqutinc) );
	libfla_test_output_op_struct_flash_only( "apcaqutinc", ops->apcaqutinc );

	// Read the operation tests for UDdate_UT.
	libfla_test_read_tests_for_op_fla_only( input_stream, &(ops->uddateut) );
	libfla_test_output_op_struct_fla_only( "uddateut", ops->uddateut );

	// Read the operation tests for UDdate_UT_inc.
	libfla_test_read_tests_for_op_flash_only( input_stream, &(ops->uddateutinc) );
	libfla_test_output_op_struct_flash_only( "uddateutinc", ops->uddateutinc );

	// Read the operation tests for Apply_QUD_UT.
	libfla_test_read_tests_for_op_front_fla_only( input_stream, &(ops->apqudut) );
	libfla_test_output_op_struct_front_fla_only( "apqudut", ops->apqudut );

	// Read the operation tests for Apply_QUD_UT_inc.
	libfla_test_read_tests_for_op_flash_only( input_stream, &(ops->apqudutinc) );
	libfla_test_output_op_struct_flash_only( "apqudutinc", ops->apqudutinc );

	// Read the operation tests for Hess_UT reduction.
	libfla_test_read_tests_for_op_fla_only( input_stream, &(ops->hessut) );
	libfla_test_output_op_struct_fla_only( "hessut", ops->hessut );

	// Read the operation tests for Tridiag_UT reduction.
	libfla_test_read_tests_for_op_fla_only( input_stream, &(ops->tridiagut) );
	libfla_test_output_op_struct_fla_only( "tridiagut", ops->tridiagut );

	// Read the operation tests for Bidiag_UT reduction.
	libfla_test_read_tests_for_op_fla_only( input_stream, &(ops->bidiagut) );
	libfla_test_output_op_struct_fla_only( "bidiagut", ops->bidiagut );

	// Read the operation tests for Eig_gest.
	libfla_test_read_tests_for_op( input_stream, &(ops->eig_gest) );
	libfla_test_output_op_struct( "eig_gest", ops->eig_gest );

	// Read the operation tests for triangular matrix inversion.
	libfla_test_read_tests_for_op( input_stream, &(ops->trinv) );
	libfla_test_output_op_struct( "trinv", ops->trinv );

	// Read the operation tests for SPD/HPD matrix inversion.
	libfla_test_read_tests_for_op_front_only( input_stream, &(ops->spdinv) );
	libfla_test_output_op_struct_front_only( "spdinv", ops->spdinv );

	// Read the operation tests for triangular Sylvester equation solve.
	libfla_test_read_tests_for_op( input_stream, &(ops->sylv) );
	libfla_test_output_op_struct( "sylv", ops->sylv );

	// Read the operation tests for triangular Lyapunov equation solve.
	libfla_test_read_tests_for_op( input_stream, &(ops->lyap) );
	libfla_test_output_op_struct( "lyap", ops->lyap );

	// Close the file.
	fclose( input_stream );

}



void libfla_test_output_op_struct( char* op_str, test_op_t op )
{
	libfla_test_output_info( "%s flash_front  %d\n", op_str, op.flash_front );
	libfla_test_output_info( "%s fla_front    %d\n", op_str, op.fla_front );
	libfla_test_output_info( "%s fla_unb_vars %d\n", op_str, op.fla_unb_vars );
	libfla_test_output_info( "%s fla_opt_vars %d\n", op_str, op.fla_opt_vars );
	libfla_test_output_info( "%s fla_blk_vars %d\n", op_str, op.fla_blk_vars );
}

void libfla_test_output_op_struct_ext( char* op_str, test_op_t op )
{
	libfla_test_output_info( "%s flash_front  %d\n", op_str, op.flash_front );
	libfla_test_output_info( "%s fla_front    %d\n", op_str, op.fla_front );
	libfla_test_output_info( "%s fla_unb_vars %d\n", op_str, op.fla_unb_vars );
	libfla_test_output_info( "%s fla_opt_vars %d\n", op_str, op.fla_opt_vars );
	libfla_test_output_info( "%s fla_blk_vars %d\n", op_str, op.fla_blk_vars );
	libfla_test_output_info( "%s fla_blk_ext  %d\n", op_str, op.fla_blk_ext );
}

void libfla_test_output_op_struct_fla_ext( char* op_str, test_op_t op )
{
	libfla_test_output_info( "%s fla_blk_ext  %d\n", op_str, op.fla_blk_ext );
}

void libfla_test_output_op_struct_flash_only( char* op_str, test_op_t op )
{
	libfla_test_output_info( "%s flash_front  %d\n", op_str, op.flash_front );
}



void libfla_test_output_op_struct_front_only( char* op_str, test_op_t op )
{
	libfla_test_output_info( "%s flash_front  %d\n", op_str, op.flash_front );
	libfla_test_output_info( "%s fla_front    %d\n", op_str, op.fla_front );
}



void libfla_test_output_op_struct_front_fla_only( char* op_str, test_op_t op )
{
	libfla_test_output_info( "%s fla_front    %d\n", op_str, op.fla_front );
}



void libfla_test_output_op_struct_fla_only( char* op_str, test_op_t op )
{
	libfla_test_output_info( "%s fla_front    %d\n", op_str, op.fla_front );
	libfla_test_output_info( "%s fla_unb_vars %d\n", op_str, op.fla_unb_vars );
	libfla_test_output_info( "%s fla_opt_vars %d\n", op_str, op.fla_opt_vars );
	libfla_test_output_info( "%s fla_blk_vars %d\n", op_str, op.fla_blk_vars );
}



void libfla_test_output_op_struct_blas3( char* op_str, test_op_t op )
{
	libfla_test_output_info( "%s flash_front  %d\n", op_str, op.flash_front );
	libfla_test_output_info( "%s fla_front    %d\n", op_str, op.fla_front );
	libfla_test_output_info( "%s fla_unb_vars %d\n", op_str, op.fla_unb_vars );
	libfla_test_output_info( "%s fla_blk_vars %d\n", op_str, op.fla_blk_vars );
	libfla_test_output_info( "%s fla_unb_ext  %d\n", op_str, op.fla_unb_ext );
}



void libfla_test_read_tests_for_op( FILE* input_stream, test_op_t* op )
{
	char buffer[ INPUT_BUFFER_SIZE ];
	int  op_switch;
	int  flash_front;
	int  fla_front;
	int  fla_unb_vars;
	int  fla_opt_vars;
	int  fla_blk_vars;

	// Read the line for the overall operation switch.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &op_switch );

	// Read the line for the FLASH front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &flash_front );

	// Read the line for the FLA front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_front );

	// Read the line for the unblocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_unb_vars );

	// Read the line for the optimized unblocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_opt_vars );

	// Read the line for the blocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_blk_vars );

	if ( op_switch == DISABLE_ALL )
	{
		op->flash_front  = DISABLE;
		op->fla_front    = DISABLE;
		op->fla_unb_vars = DISABLE;
		op->fla_opt_vars = DISABLE;
		op->fla_blk_vars = DISABLE;
	}
	else
	{
		op->flash_front  = flash_front;
		op->fla_front    = fla_front;
		op->fla_unb_vars = fla_unb_vars;
		op->fla_opt_vars = fla_opt_vars;
		op->fla_blk_vars = fla_blk_vars;
	}
}

void libfla_test_read_tests_for_op_ext( FILE* input_stream, test_op_t* op )
{
	char buffer[ INPUT_BUFFER_SIZE ];
	int  op_switch;
	int  flash_front;
	int  fla_front;
	int  fla_unb_vars;
	int  fla_opt_vars;
	int  fla_blk_vars;
	int  fla_blk_ext;

	// Read the line for the overall operation switch.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &op_switch );

	// Read the line for the FLASH front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &flash_front );

	// Read the line for the FLA front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_front );

	// Read the line for the unblocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_unb_vars );

	// Read the line for the optimized unblocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_opt_vars );

	// Read the line for the blocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_blk_vars );

	// Read the line for the blocked external variant.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_blk_ext );

	if ( op_switch == DISABLE_ALL )
	{
		op->flash_front  = DISABLE;
		op->fla_front    = DISABLE;
		op->fla_unb_vars = DISABLE;
		op->fla_opt_vars = DISABLE;
		op->fla_blk_vars = DISABLE;
		op->fla_blk_ext  = DISABLE;
	}
	else
	{
		op->flash_front  = flash_front;
		op->fla_front    = fla_front;
		op->fla_unb_vars = fla_unb_vars;
		op->fla_opt_vars = fla_opt_vars;
		op->fla_blk_vars = fla_blk_vars;
		op->fla_blk_ext  = fla_blk_ext;
	}
}

void libfla_test_read_tests_for_op_fla_ext( FILE* input_stream, test_op_t* op )
{
	char buffer[ INPUT_BUFFER_SIZE ];
	int  op_switch;
	int  fla_blk_ext;

    // Read the line for the overall operation switch.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &op_switch );

    // Read the line for the blocked external variant.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_blk_ext );

	if ( op_switch == DISABLE_ALL )
	{
		op->fla_blk_ext  = DISABLE;
	}
	else
	{
		op->fla_blk_ext  = fla_blk_ext;
	}
}

void libfla_test_read_tests_for_op_flash_only( FILE* input_stream, test_op_t* op )
{
	char buffer[ INPUT_BUFFER_SIZE ];
	int  op_switch;
	int  flash_front;

	// Read the line for the overall operation switch.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &op_switch );

	// Read the line for the FLASH front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &flash_front );

	if ( op_switch == DISABLE_ALL )
	{
		op->flash_front  = DISABLE;
	}
	else
	{
		op->flash_front  = flash_front;
	}

	op->fla_front    = DISABLE; // not used
	op->fla_unb_vars = DISABLE; // not used
	op->fla_opt_vars = DISABLE; // not used
	op->fla_blk_vars = DISABLE; // not used
}



void libfla_test_read_tests_for_op_fla_only( FILE* input_stream, test_op_t* op )
{
	char buffer[ INPUT_BUFFER_SIZE ];
	int  op_switch;
	int  fla_front;
	int  fla_unb_vars;
	int  fla_opt_vars;
	int  fla_blk_vars;

	// Read the line for the overall operation switch.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &op_switch );

	// Read the line for the FLA front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_front );

	// Read the line for the unblocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_unb_vars );

	// Read the line for the optimized unblocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_opt_vars );

	// Read the line for the blocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_blk_vars );

	if ( op_switch == DISABLE_ALL )
	{
		op->fla_front    = DISABLE;
		op->fla_unb_vars = DISABLE;
		op->fla_opt_vars = DISABLE;
		op->fla_blk_vars = DISABLE;
	}
	else
	{
		op->fla_front    = fla_front;
		op->fla_unb_vars = fla_unb_vars;
		op->fla_opt_vars = fla_opt_vars;
		op->fla_blk_vars = fla_blk_vars;
	}

	op->flash_front  = DISABLE; // not used
}



void libfla_test_read_tests_for_op_front_only( FILE* input_stream, test_op_t* op )
{
	char buffer[ INPUT_BUFFER_SIZE ];
	int  op_switch;
	int  flash_front;
	int  fla_front;

	// Read the line for the overall operation switch.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &op_switch );

	// Read the line for the FLASH front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &flash_front );

	// Read the line for the FLA front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_front );

	if ( op_switch == DISABLE_ALL )
	{
		op->flash_front  = DISABLE;
		op->fla_front    = DISABLE;
	}
	else
	{
		op->flash_front  = flash_front;
		op->fla_front    = fla_front;
	}

	op->fla_unb_vars = DISABLE; // not used
	op->fla_opt_vars = DISABLE; // not used
	op->fla_blk_vars = DISABLE; // not used
}



void libfla_test_read_tests_for_op_front_fla_only( FILE* input_stream, test_op_t* op )
{
	char buffer[ INPUT_BUFFER_SIZE ];
	int  op_switch;
	int  fla_front;

	// Read the line for the overall operation switch.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &op_switch );

	// Read the line for the FLA front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_front );

	if ( op_switch == DISABLE_ALL )
	{
		op->fla_front    = DISABLE;
	}
	else
	{
		op->fla_front    = fla_front;
	}

	op->flash_front  = DISABLE; // not used
	op->fla_unb_vars = DISABLE; // not used
	op->fla_opt_vars = DISABLE; // not used
	op->fla_blk_vars = DISABLE; // not used
}



void libfla_test_read_tests_for_op_blas3( FILE* input_stream, test_op_t* op )
{
	char buffer[ INPUT_BUFFER_SIZE ];
	int  op_switch;
	int  flash_front;
	int  fla_front;
	int  fla_unb_vars;
	int  fla_blk_vars;
	int  fla_unb_ext;

	// Read the line for the overall operation switch.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &op_switch );

	// Read the line for the FLASH front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &flash_front );

	// Read the line for the FLA front-end.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_front );

	// Read the line for the unblocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_unb_vars );

	// Read the line for the blocked variants.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_blk_vars );

	// Read the line for the unblocked external implementation.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &fla_unb_ext );

	if ( op_switch == DISABLE_ALL )
	{
		op->flash_front  = DISABLE;
		op->fla_front    = DISABLE;
		op->fla_unb_vars = DISABLE;
		op->fla_blk_vars = DISABLE;
		op->fla_unb_ext  = DISABLE;
	}
	else
	{
		op->flash_front  = flash_front;
		op->fla_front    = fla_front;
		op->fla_unb_vars = fla_unb_vars;
		op->fla_blk_vars = fla_blk_vars;
		op->fla_unb_ext  = fla_unb_ext ;
	}

	op->fla_opt_vars = DISABLE; // not used
	op->fla_blk_ext  = DISABLE; // not used
}



void libfla_test_read_parameter_file( char* input_filename, test_params_t* params )
{
	FILE* input_stream;
	char  buffer[ INPUT_BUFFER_SIZE ];
	char  temp[ INPUT_BUFFER_SIZE ];
	int   i;

	// Attempt to open input file corresponding to input_filename as
	// read-only/binary.
	input_stream = fopen( input_filename, "rb" );

	// Check for success.
	if ( input_stream == NULL )
	{
		libfla_test_output_error( "Failed to open input file %s. Check existence and permissions.\n",
		                          input_filename );
	}

	// Read the number of repeats.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%u ", &(params->n_repeats) );

	// Read the storage schemes to test. We should have at most three: 'r' for
	// row-major, 'c' for column-major, and 'g' for general strides, OR just
	// 'm' for mixed storage.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%s ", temp );

	params->n_storage = strlen( temp );
	if ( params->n_storage > MAX_NUM_STORAGE )
	{
		libfla_test_output_error( "Detected too many storage schemes (%u) in input file.\n",
		                          params->n_storage );
	}
	strcpy( params->storage, temp );

	// If 'm' is in the string, then remove all other chars.
	for ( i = 0; i < params->n_storage; ++i )
	{
		if ( params->storage[i] == 'm' )
		{
			sprintf( params->storage, "m" );
			break;
		}
	}

	// Read the datatypes to test. We should have at most four: 's', 'd', 'c',
	// and 'z'.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%s ", temp );

	params->n_datatypes = strlen( temp );
	if ( params->n_datatypes > MAX_NUM_DATATYPES )
	{
		libfla_test_output_error( "Detected too many datatype requests (%u) in input file.\n",
		                          params->n_datatypes );
	}

	for( i = 0; i < params->n_datatypes; ++i )
	{
		if      ( temp[i] == 's' ) params->datatype[i] = FLA_FLOAT;
		else if ( temp[i] == 'd' ) params->datatype[i] = FLA_DOUBLE;
		else if ( temp[i] == 'c' ) params->datatype[i] = FLA_COMPLEX;
		else if ( temp[i] == 'z' ) params->datatype[i] = FLA_DOUBLE_COMPLEX;

		params->datatype_char[i] = temp[i];
	}

	// Read the blocksize to use for blocked algorithms on flat matrices.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%lu ", &(params->b_alg_flat) );

	// Read the algorithmic blocksize to use for algorithms-by-blocks.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%lu ", &(params->b_alg_hier) );

	// Read the storage (FLASH) blocksize to use for algorithms-by-blocks.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%lu ", &(params->b_flash) );

	// Read the initial problem size to test.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%lu ", &(params->p_first) );

	// Read the maximum problem size to test.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%lu ", &(params->p_max) );

	// Read the problem size increment to test.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%lu ", &(params->p_inc) );

	// Read the partial number of matrix size for incomplete factorization.
        libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%d ", &(params->p_nfact) );

        // Read the number of SuperMatrix threads to test with.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%u ", &(params->n_threads) );

	// Read the requested course of action if a test fails.
	libfla_test_read_next_line( buffer, input_stream );
	sscanf( buffer, "%c ", &(params->reaction_to_failure) );

	if ( params->reaction_to_failure != ON_FAILURE_IGNORE_CHAR &&
	     params->reaction_to_failure != ON_FAILURE_SLEEP_CHAR  &&
	     params->reaction_to_failure != ON_FAILURE_ABORT_CHAR  )
	{
		libfla_test_output_error( "Invalid reaction-to-failure character code (%c) in input file.\n",
		                          params->reaction_to_failure );
	}

	// Close the file.
	fclose( input_stream );

	libfla_test_output_info( "\n" );
	libfla_test_output_info( "--- test suite parameters ----------------------------\n" );
	libfla_test_output_info( "\n" );
	libfla_test_output_info( "n_repeats            %u\n", params->n_repeats );
	libfla_test_output_info( "n_storage            %u\n", params->n_storage );
	libfla_test_output_info( "storage              %s\n", params->storage );
	libfla_test_output_info( "n_datatypes          %u\n", params->n_datatypes );
	libfla_test_output_info( "datatype[0]          %d (%c)\n", params->datatype[0],
	                                                         params->datatype_char[0] );
	for( i = 1; i < params->n_datatypes; ++i )
	libfla_test_output_info( "        [%d]          %d (%c)\n", i, params->datatype[i],
	                                                             params->datatype_char[i] );
	libfla_test_output_info( "b_alg_flat           %u\n", params->b_alg_flat );
	libfla_test_output_info( "b_alg_hier           %u\n", params->b_alg_hier );
	libfla_test_output_info( "b_flash              %u\n", params->b_flash );
	libfla_test_output_info( "p_first              %u\n", params->p_first );
	libfla_test_output_info( "p_max                %u\n", params->p_max );
	libfla_test_output_info( "p_inc                %u\n", params->p_inc );
	libfla_test_output_info( "p_nfact              %d\n", params->p_nfact );
	libfla_test_output_info( "n_threads            %u\n", params->n_threads );
	libfla_test_output_info( "reaction_to_failure  %c\n", params->reaction_to_failure );
}



void libfla_test_read_next_line( char* buffer, FILE* input_stream )
{
	char temp[ INPUT_BUFFER_SIZE ];

	// We want to read at least one line, so we use a do-while loop.
	do
	{
		// Read the next line into a temporary buffer and check success.
		if ( fgets( temp, INPUT_BUFFER_SIZE-1, input_stream ) == NULL )
		{
			if ( feof( input_stream ) )
				libfla_test_output_error( "Error reading input file: encountered unexpected EOF." );
			else
				libfla_test_output_error( "Error (non-EOF) reading input file." );
		}
	}
    // We continue to read lines into buffer until the line is neither
	// commented nor blank.
	while ( temp[0] == COMMENT_CHAR || temp[0] == '\n' ||
	        temp[0] == ' '          || temp[0] == '\t' );


	// Save the string in temp, up to first white space character, into buffer.
	sscanf( temp, "%s ", buffer );
}



void libfla_test_output_info( char* message, ... )
{
	FILE* output_stream = stdout;
    va_list args;

    //fprintf( output_stream, "%s: ", libfla_test_binary_name );

    // Initialize variable argument environment.
    va_start( args, message );

    // Parse the received message and print its components.
    libfla_test_parse_message( output_stream, message, args );

    // Shutdown variable argument environment and clean up stack.
    va_end( args );

	// Flush the output stream.
    fflush( output_stream );
}



void libfla_test_output_error( char* message, ... )
{
	FILE*   output_stream = stderr;
    va_list args;

    fprintf( output_stream, "%s: *** error ***: ", libfla_test_binary_name );

    // Initialize variable argument environment.
    va_start( args, message );

    // Parse the received message and print its components.
    libfla_test_parse_message( output_stream, message, args );

    // Shutdown variable argument environment and clean up stack.
    va_end( args );

	// Flush the output stream.
    fflush( output_stream );

	// Exit.
	exit(1);
}



void libfla_test_parse_message( FILE* output_stream, char* message, va_list args )
{
	int           c, cf;
	char          format_spec[8];
	unsigned int  the_uint;
	int           the_int;
	double        the_double;
	char*         the_string;
	char          the_char;

	// Begin looping over message to insert variables wherever there are
	// format specifiers.
	for ( c = 0; message[c] != '\0'; )
	{
		if ( message[c] != '%' )
		{
			fprintf( output_stream, "%c", message[c] );
			c += 1;
		}
		else if ( message[c] == '%' && message[c+1] == '%' ) // handle escaped '%' chars.
		{
			fprintf( output_stream, "%c", message[c] );
			c += 2;
		}
		else
		{
			// Save the format string if there is one.
			format_spec[0] = '%';
			for ( c += 1, cf = 1; strchr( "udefsc", message[c] ) == NULL; ++c, ++cf )
			{
				format_spec[cf] = message[c];
			}

			// Add the final type specifier, and null-terminate the string.
			format_spec[cf] = message[c];
			format_spec[cf+1] = '\0';

			// Switch based on type, since we can't predict what will
			// va_args() will return.
			switch ( message[c] )
			{
				case 'u':
				the_uint = va_arg( args, unsigned int );
				fprintf( output_stream, format_spec, the_uint );
				break;

				case 'd':
				the_int = va_arg( args, int );
				fprintf( output_stream, format_spec, the_int );
				break;

				case 'e':
				the_double = va_arg( args, double );
				fprintf( output_stream, format_spec, the_double );
				break;

				case 'f':
				the_double = va_arg( args, double );
				fprintf( output_stream, format_spec, the_double );
				break;

				case 's':
				the_string = va_arg( args, char* );
				//fprintf( output_stream, "%s", the_string );
				fprintf( output_stream, format_spec, the_string );
				break;

				case 'c':
				the_char = va_arg( args, int );
				fprintf( output_stream, "%c", the_char );
				break;
			}

			// Move to next character past type specifier.
			c += 1;
		}
	}
}



void libfla_test_parse_command_line( int argc, char** argv )
{
	if ( argc > 1 )
	{
		fprintf( stderr, "Too many command line arguments.\n" );
		exit(1);
	}
	// Copy the binary name to a global string so we can use it later.
	strncpy( libfla_test_binary_name, argv[0], MAX_BINARY_NAME_LENGTH );
}



char* libfla_test_get_string_for_result( double         residual,
                                         FLA_Datatype   datatype,
                                         test_thresh_t* thresh )
{
	char* r_val;

	if      ( datatype == FLA_FLOAT )
	{
		if      ( residual > thresh->failwarn_s ) r_val = libfla_test_fail_string;
		else if ( residual > thresh->warnpass_s ) r_val = libfla_test_warn_string;
		else                                      r_val = libfla_test_pass_string;
	}
	else if ( datatype == FLA_DOUBLE )
	{
		if      ( residual > thresh->failwarn_d ) r_val = libfla_test_fail_string;
		else if ( residual > thresh->warnpass_d ) r_val = libfla_test_warn_string;
		else                                      r_val = libfla_test_pass_string;
	}
	else if ( datatype == FLA_COMPLEX )
	{
		if      ( residual > thresh->failwarn_c ) r_val = libfla_test_fail_string;
		else if ( residual > thresh->warnpass_c ) r_val = libfla_test_warn_string;
		else                                      r_val = libfla_test_pass_string;
	}
	else // if ( datatype == FLA_DOUBLE_COMPLEX )
	{
		if      ( residual > thresh->failwarn_z ) r_val = libfla_test_fail_string;
		else if ( residual > thresh->warnpass_z ) r_val = libfla_test_warn_string;
		else                                      r_val = libfla_test_pass_string;
	}

	return r_val;
}



void libfla_test_init_strings( void )
{
	sprintf( libfla_test_pass_string, "PASS" );
	sprintf( libfla_test_warn_string, "MARGINAL" );
	sprintf( libfla_test_fail_string, "FAILURE" );
	sprintf( libfla_test_storage_format_string, "Row(r) and General(g) storage format is not supported\n \ 
							\t\t  by External LAPACK interface" );

	sprintf( libfla_test_stor_chars, STORAGE_SCHEME_CHARS );
}



void libfla_test_fill_storage_strings( char** sc_str, unsigned int n_storage_run,
                                                      unsigned int n_matrices )
{
	unsigned int  sci, mi, i;
	unsigned int* c;

	// Allocate an array with one element per matrix argument. We will use
	// this array to keep track of our progress as we canonically move
	// though all possible storage combinations.
	c = ( unsigned int* ) malloc( n_matrices * sizeof( unsigned int ) );

	// Initialize all values in c to zero.
	for ( i = 0; i < n_matrices; ++i ) c[i] = 0;

	for ( sci = 0; sci < n_storage_run; ++sci )
	{
		// Iterate backwards since we want to form (for example):
		// (1) ccc, (2) ccr, (3) crc, (4) crr, etc.
		for ( i = 0, mi = n_matrices - 1; i < n_matrices; --mi, ++i )
		{
			// Map the current values in c to storage characters.
			sc_str[sci][mi] = libfla_test_stor_chars[ c[mi] ];
		}

		// Terminate the string.
		sc_str[sci][n_matrices] = '\0';

		// Only try to increment/carryover if this is NOT the last storage
		// run/combo.
		if ( sci < n_storage_run - 1 )
		{
			// Increment the least-most significant counter.
			c[ n_matrices - 1 ]++;

			// Perform "carryover" if needed.
			carryover( &c[ n_matrices - 1 ], n_matrices );
		}
	}


	// Free the array.
	free( c );
}


void carryover( unsigned int* c, unsigned int n_matrices )
{
	if ( n_matrices == 1 ) return;
	else
	{
		if ( *c == NUM_STORAGE_CHARS )
		{
			*c = 0;
			*(c-1) += 1;
			carryover( c-1, n_matrices-1 );
		}
	}
}


void libfla_test_op_driver( char*         func_str,
                            char*         impl_var_str,
                            unsigned int  first_var,
                            unsigned int  last_var,
                            unsigned int  n_pc,
                            char**        pc_str,
                            unsigned int  n_matrices,
                            signed int    impl,
                            test_params_t params,
                            test_thresh_t thresh,
                            void (*f_exp) (test_params_t, // params
                                           unsigned int,  // var
                                           char*,         // sc_cur_str (current storage string)
                                           FLA_Datatype,  // datatype
                                           uinteger,  // p_cur
                                           unsigned int,  // pci (param combo counter)
                                           unsigned int,  // n_repeats
                                           signed int,    // impl
                                           double*,       // perf
                                           double*,       //time
                                           double* ) )    // residual
{
	unsigned int n_threads           = params.n_threads;
	unsigned int n_storage           = params.n_storage;
	unsigned int n_datatypes         = params.n_datatypes;
	uinteger p_first             = params.p_first;
	uinteger p_max               = params.p_max;
	uinteger p_inc               = params.p_inc;
	uinteger n_repeats           = params.n_repeats;
	uinteger reaction_to_failure = params.reaction_to_failure;
	uinteger sci, dt, p_cur, mat, pci, var;
	char         datatype_char;
	FLA_Datatype datatype;
	double       perf, time, residual;
	char*        pass_str;
	char         blank_str[32];
	char         func_param_str[64];
	unsigned int n_spaces;
	unsigned int n_storage_run;
	char**       sc_str;

	// Set the number of threads and/or disable SuperMatrix.
	if ( n_threads == 0 ) FLASH_Queue_disable();
	else                  FLASH_Queue_set_num_threads( n_threads );

	// Execute the variant loop only once if we're runing a front-end test.
	if ( impl == FLA_TEST_HIER_FRONT_END ||
	     impl == FLA_TEST_FLAT_FRONT_END ||
	     impl == FLA_TEST_FLAT_UNB_EXT   ||
	     impl == FLA_TEST_FLAT_BLK_EXT )
	{
		first_var = 0;
		last_var  = 0;
	}

	// Determine the total number of storage schemes.
	if ( params.storage[0] == 'm' )
	{
		// Prepare to run all NUM_STORAGE_SCHEMES combinations for each
		// matrix argument.
		n_storage_run = ( unsigned int ) pow( ( double ) NUM_STORAGE_CHARS,
		                                      ( double ) n_matrices );

		sc_str = ( char** ) malloc( n_storage_run * sizeof( char* ) );
		for ( sci = 0; sci < n_storage_run; ++sci )
			sc_str[sci] = ( char* ) malloc( ( n_matrices + 1 ) * sizeof( char ) );

		libfla_test_fill_storage_strings( sc_str, n_storage_run, n_matrices );
	}
	else // if ( params.storage[0] == 'c' ||
	     //      params.storage[0] == 'r' ||
	     //      params.storage[0] == 'g' )
	{
		// Only run combinations where all matrices are stored in one
		// storage scheme or another (no mixed storage).
		n_storage_run = n_storage;

		sc_str = ( char** ) malloc( n_storage_run * sizeof( char* ) );
		for ( sci = 0; sci < n_storage_run; ++sci )
		{
			// Allocate a string for a storage combination.
			sc_str[sci] = ( char* ) malloc( ( n_matrices + 1 ) * sizeof( char ) );

			// Fill the string with the current storage scheme character
			// for each matrix operand.
			for ( mat = 0; mat < n_matrices; ++mat )
				sc_str[sci][mat] = params.storage[sci];
			sc_str[sci][n_matrices] = '\0';
		}
	}
  
  libfla_test_output_info( "%3sAPI%28s DATA_TYPE%4s SIZE%1s FLOPS%2s TIME(s)%6s ERROR%5s STATUS\n", "", "", "", "", "", "", "" );
  libfla_test_output_info( "%3s====%28s==========%4s====%1s=======%2s========%5s==========%2s========\n", "", "", "", "", "", "", "" );
	// Loop over variant, if applicable.
	for ( var = first_var; var <= last_var; ++var )
	{
		// Loop over the requested storage schemes.
		for ( sci = 0; sci < n_storage_run; ++sci )
		{
			// Loop over the requested datatypes.
			for ( dt = 0; dt < n_datatypes; ++dt )
			{
				datatype      = params.datatype[dt];
				datatype_char = params.datatype_char[dt];

				// Loop over the requested problem sizes.
				for ( p_cur = p_first; p_cur <= p_max; p_cur += p_inc )
				{
					// Loop over the operation's parameter combinations.
					for ( pci = 0; pci < n_pc; ++pci )
					{
						//If external interface is selected and row or general 
						//storage is set, then do not proceed. Row and General storage 
						//is not supported for external lapack interface
						if (impl == FLA_TEST_FLAT_BLK_EXT && 
							(sc_str[sci][0] == 'r' || sc_str[sci][0] == 'g'))
						{
						  pass_str = libfla_test_storage_format_string;
						  perf = time = residual = 0.0f;
						}
						else
						{
						  f_exp( params,
						       var,
						       sc_str[sci],
						       datatype,
						       p_cur, pci, n_repeats, impl,
						       &perf, &time, &residual );

						  pass_str = libfla_test_get_string_for_result( residual,
						                                              datatype,
						                                              &thresh );
						}

						// Output the results. Use different formats depending on
						// whether the results are from a front-end or variant.
						libfla_test_build_function_string( func_str, impl,
						                                   impl_var_str, var,
						                                   n_pc, pc_str[pci],
						                                   func_param_str );

						n_spaces = MAX_FUNC_STRING_LENGTH - strlen( func_param_str );
						fill_string_with_n_spaces( blank_str, n_spaces );
            
						libfla_test_output_info( "   %s%s  %c|%-6s  %5u  %6.3lf  %6.10lf  %9.2le   %s\n",
						                         func_param_str, blank_str,
						                         datatype_char, sc_str[sci],
						                         p_cur, perf, time, residual, pass_str );

						// If we need to check whether to do something on failure,
						// do so now.
						if ( reaction_to_failure == ON_FAILURE_SLEEP_CHAR )
						{
							if ( strstr( pass_str, FLA_TEST_FAIL_STRING ) == pass_str )
								libfla_test_sleep();
						}
						else if ( reaction_to_failure == ON_FAILURE_ABORT_CHAR )
						{
							if ( strstr( pass_str, FLA_TEST_FAIL_STRING ) == pass_str )
								libfla_test_abort();
						}
					}
				}

				libfla_test_output_info( "\n" );
			}
		}
	}

	for ( sci = 0; sci < n_storage_run; ++sci )
		free( sc_str[sci] );
	free( sc_str );
}



void libfla_test_print_result_info(char  *func_param_str,
                                   char  *datatype_char,
                                   char  *sc_str,
                                   integer    p_cur,
                                   double perf,
                                   double time,
                                   double residual,
                                   char  *pass_str,
                                   int    nfact )
{
	char blank_str[32];
	integer  n_spaces;

	n_spaces = MAX_FUNC_STRING_LENGTH - strlen( func_param_str );
	fill_string_with_n_spaces( blank_str, n_spaces );
  libfla_test_output_info( "   %s%s  %c|%-6s  %5u  %6.3lf  %6.10lf  %9.2le   %s for nfact=%d\n",
                                 func_param_str, blank_str,
                                 datatype_char, sc_str,
                                 p_cur, perf, time, residual, pass_str, nfact );
}



void libfla_test_build_function_string( char*        func_base_str,
                                        signed int   impl,
                                        char*        impl_var_str,
                                        unsigned int var,
                                        unsigned int n_pc,
                                        char*        pc_str,
                                        char*        func_str )
{

	sprintf( func_str, "%s", func_base_str );

	if ( impl == FLA_TEST_HIER_FRONT_END || impl == FLA_TEST_FLAT_FRONT_END )
	{
		//sprintf( &func_str[strlen(func_str)], "()" );
		if ( n_pc > 1 )
			sprintf( &func_str[strlen(func_str)], ":%s", pc_str );
	}
	else if ( impl == FLA_TEST_FLAT_UNB_EXT || impl == FLA_TEST_FLAT_BLK_EXT )
	{
		sprintf( &func_str[strlen(func_str)], "_%s", impl_var_str );

		if ( n_pc > 1 )
			sprintf( &func_str[strlen(func_str)], ":%s", pc_str );
	}
	else
	{
		if ( n_pc > 1 )
			sprintf( &func_str[strlen(func_str)], "_%s", pc_str );

		sprintf( &func_str[strlen(func_str)], "_%s%u", impl_var_str, var );
	}
}



void fill_string_with_n_spaces( char* str, unsigned int n_spaces )
{
	unsigned int i;

	for ( i = 0; i < n_spaces; ++i )
		sprintf( &str[i], " " );
}



void libfla_test_obj_create( FLA_Datatype dt, FLA_Trans trans, char storage, dim_t m, dim_t n, FLA_Obj* A )
{
	dim_t m_trans = m;
	dim_t n_trans = n;
	dim_t rs_g;
	dim_t cs_g;

	if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
	{
		m_trans = n;
		n_trans = m;
	}

	// In case of general strides, use strides with a column-major tilt.
	rs_g = 2 * 1;
	cs_g = 2 * m_trans;

	if      ( storage == 'c' ) FLA_Obj_create( dt, m_trans, n_trans, 0,       0, A );
	else if ( storage == 'r' ) FLA_Obj_create( dt, m_trans, n_trans, n_trans, 1, A );
	else if ( storage == 'g' ) FLA_Obj_create( dt, m_trans, n_trans, rs_g, cs_g, A );
	else                       FLA_Abort();
}



void libfla_test_sleep( void )
{
	int i;

	libfla_test_output_info( "Resuming in " );
	for ( i = SECONDS_TO_SLEEP; i > 0; --i )
	{
		libfla_test_output_info( "%d ", i );
#ifdef _WIN32
		Sleep(1);
#else
        sleep(1);
#endif
	}
	libfla_test_output_info( "\n" );
}


void libfla_test_abort( void )
{
	abort();
}

