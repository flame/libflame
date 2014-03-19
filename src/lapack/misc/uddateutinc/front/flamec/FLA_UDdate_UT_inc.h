/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLASH_UDdate_UT_inc( FLA_Obj R,
                               FLA_Obj C,
                               FLA_Obj D, FLA_Obj T, FLA_Obj W );

FLA_Error FLA_UDdate_UT_inc_blk_var1( FLA_Obj R,
                                      FLA_Obj C,
                                      FLA_Obj D, FLA_Obj T, FLA_Obj W, fla_uddateutinc_t* cntl );

FLA_Error FLASH_UDdate_UT_inc_create_hier_matrices( FLA_Obj R_flat, FLA_Obj C_flat, FLA_Obj D_flat, dim_t depth, dim_t* b_flash, dim_t b_alg, FLA_Obj* R, FLA_Obj* C, FLA_Obj* D, FLA_Obj* T, FLA_Obj* W );
dim_t     FLASH_UDdate_UT_inc_determine_alg_blocksize( FLA_Obj R );

FLA_Error FLASH_UDdate_UT_inc_update_rhs( FLA_Obj T, FLA_Obj bR,
                                          FLA_Obj C, FLA_Obj bC,
                                          FLA_Obj D, FLA_Obj bD );
FLA_Error FLASH_UDdate_UT_inc_solve( FLA_Obj R, FLA_Obj bR, FLA_Obj x );
