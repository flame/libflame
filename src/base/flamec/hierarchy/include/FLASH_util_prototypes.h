/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLASH utility routine prototypes ----------------------------------------

double    FLASH_Max_elemwise_diff( FLA_Obj A, FLA_Obj B );

FLA_Error FLASH_Random_matrix( FLA_Obj H );
FLA_Error FLASH_Random_spd_matrix( FLA_Uplo uplo, FLA_Obj H );

FLA_Error FLASH_Norm1( FLA_Obj H, FLA_Obj norm );
FLA_Error FLASH_Obj_shift_diagonal( FLA_Conj conj, FLA_Obj sigma, FLA_Obj H );

FLA_Error FLASH_Set( FLA_Obj alpha, FLA_Obj H );

FLA_Error FLASH_Obj_create_diag_panel( FLA_Obj A, FLA_Obj* U );

FLA_Error FLASH_LU_find_zero_on_diagonal( FLA_Obj A );

FLA_Error FLASH_Triangularize( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLASH_Hermitianize( FLA_Uplo uplo, FLA_Obj A );

// --- FLASH utility check routine prototypes ----------------------------------

FLA_Error FLASH_LU_find_zero_on_diagonal_check( FLA_Obj A );

