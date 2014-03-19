/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLASH_QR_UT_inc( FLA_Obj A, FLA_Obj TW );

FLA_Error FLASH_QR_UT_inc_noopt( FLA_Obj A, FLA_Obj TW );
FLA_Error FLASH_QR_UT_inc_opt1( FLA_Obj A, FLA_Obj TW );

FLA_Error FLA_QR_UT_inc_blk_var1( FLA_Obj A, FLA_Obj TW, fla_qrutinc_t* cntl );
FLA_Error FLA_QR_UT_inc_blk_var2( FLA_Obj A, FLA_Obj TW, FLA_Obj U, fla_qrutinc_t* cntl );

FLA_Error FLASH_QR_UT_inc_create_hier_matrices( FLA_Obj A_flat, dim_t depth, dim_t* b_flash, dim_t b_alg, FLA_Obj* A, FLA_Obj* TW );
dim_t     FLASH_QR_UT_inc_determine_alg_blocksize( FLA_Obj A );

FLA_Error FLASH_QR_UT_inc_solve( FLA_Obj A, FLA_Obj TW, FLA_Obj B, FLA_Obj X );

