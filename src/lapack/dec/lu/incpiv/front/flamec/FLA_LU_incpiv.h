/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_LU_incpiv_aux.h"

FLA_Error FLASH_LU_incpiv_create_hier_matrices( FLA_Obj A_flat, dim_t depth, dim_t* b_flash, dim_t b_alg, FLA_Obj* A, FLA_Obj* p, FLA_Obj* L );
dim_t     FLASH_LU_incpiv_determine_alg_blocksize( FLA_Obj A );

FLA_Error FLASH_LU_incpiv_noopt( FLA_Obj A, FLA_Obj p, FLA_Obj L );
FLA_Error FLASH_LU_incpiv_opt1( FLA_Obj A, FLA_Obj p, FLA_Obj L );

FLA_Error FLASH_LU_incpiv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj B, FLA_Obj X );
