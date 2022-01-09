/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

TLS_CLASS_SPEC fla_qr2ut_t*     flash_qr2ut_cntl_leaf = NULL;
TLS_CLASS_SPEC fla_qr2ut_t*     flash_qr2ut_cntl = NULL;
TLS_CLASS_SPEC fla_blocksize_t* flash_qr2ut_var2_bsize = NULL;

void FLASH_QR2_UT_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_qr2ut_var2_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to invoke variant 1.
	flash_qr2ut_cntl_leaf = FLA_Cntl_qr2ut_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM, 
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree to invoke variant 2.
	flash_qr2ut_cntl    = FLA_Cntl_qr2ut_obj_create( FLA_HIER,
	                                                 FLA_BLOCKED_VARIANT2, 
	                                                 flash_qr2ut_var2_bsize,
	                                                 flash_qr2ut_cntl_leaf,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL,
	                                                 NULL );
}

void FLASH_QR2_UT_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_qr2ut_cntl_leaf );
	FLA_Cntl_obj_free( flash_qr2ut_cntl );

	FLA_Blocksize_free( flash_qr2ut_var2_bsize );
}

