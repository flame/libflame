/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLASH_CAQR_UT_inc( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW );

FLA_Error FLASH_CAQR_UT_inc_noopt( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW );

FLA_Error FLASH_CAQR_UT_inc_create_hier_matrices( dim_t p, FLA_Obj A_flat, dim_t depth, dim_t* b_flash, dim_t b_alg, FLA_Obj* A, FLA_Obj* ATW, FLA_Obj* R, FLA_Obj* RTW );
dim_t     FLASH_CAQR_UT_inc_determine_alg_blocksize( FLA_Obj A );
FLA_Error FLASH_CAQR_UT_inc_adjust_views( FLA_Obj A, FLA_Obj TW );

void      FLA_CAQR_UT_inc_init_structure( dim_t p, dim_t nb_part, FLA_Obj R );

dim_t     FLA_CAQR_UT_inc_compute_blocks_per_part( dim_t p, FLA_Obj A );

FLA_Error FLA_CAQR_UT_inc_factorize_panels( dim_t nb_part, FLA_Obj A, FLA_Obj ATW );

FLA_Error FLA_CAQR_UT_inc_copy_triangles( dim_t nb_part, FLA_Obj A, FLA_Obj R );

FLA_Error FLA_CAQR_UT_inc_blk_var1( FLA_Obj R, FLA_Obj TW, fla_caqrutinc_t* cntl );

FLA_Error FLASH_CAQR_UT_inc_solve( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW, FLA_Obj B, FLA_Obj X );

