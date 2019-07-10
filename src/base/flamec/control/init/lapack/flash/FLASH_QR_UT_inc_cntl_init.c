/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern fla_qrut_t*   flash_qrut_cntl_leaf;
extern fla_apqut_t*  flash_apqut_cntl;
extern fla_qr2ut_t*  flash_qr2ut_cntl;
extern fla_apq2ut_t* flash_apq2ut_cntl;

fla_qrutinc_t*       flash_qrutinc_cntl = NULL;
fla_blocksize_t*     flash_qrutinc_var1_bsize = NULL;

void FLASH_QR_UT_inc_cntl_init()
{
	// Set blocksizes for hierarchical storage.
	flash_qrutinc_var1_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to invoke variant 1.
	flash_qrutinc_cntl = FLA_Cntl_qrutinc_obj_create( FLA_HIER,
	                                                  FLA_BLOCKED_VARIANT1, 
	                                                  flash_qrutinc_var1_bsize,
	                                                  flash_qrut_cntl_leaf,
	                                                  flash_apqut_cntl,
	                                                  flash_qr2ut_cntl,
	                                                  flash_apq2ut_cntl );
}

void FLASH_QR_UT_inc_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_qrutinc_cntl );

	FLA_Blocksize_free( flash_qrutinc_var1_bsize );
}

