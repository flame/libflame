/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- top-level wrapper prototypes --------------------------------------------

// Implemented:
FLA_Error FLASH_Chol( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLASH_Chol_solve( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_LU_nopiv( FLA_Obj A );
FLA_Error FLASH_LU_nopiv_solve( FLA_Obj A, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_LU_piv( FLA_Obj A, FLA_Obj p );
FLA_Error FLASH_LU_piv_solve( FLA_Obj A, FLA_Obj p, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_LU_incpiv( FLA_Obj A, FLA_Obj p, FLA_Obj L );
FLA_Error FLASH_FS_incpiv( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj b );
FLA_Error FLASH_Trinv( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLASH_Ttmm( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLASH_SPDinv( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLASH_Sylv( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLASH_Apply_Q_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B );
FLA_Error FLASH_Apply_Q2_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E );
FLA_Error FLASH_QR2_UT( FLA_Obj B, FLA_Obj D, FLA_Obj T );
FLA_Error FLASH_QR_UT( FLA_Obj A, FLA_Obj TW );
FLA_Error FLASH_QR_UT_solve( FLA_Obj A, FLA_Obj T, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_QR_UT_inc( FLA_Obj A, FLA_Obj TW );
FLA_Error FLASH_QR_UT_inc_solve( FLA_Obj A, FLA_Obj TW, FLA_Obj B, FLA_Obj X );
FLA_Error FLASH_Apply_Q_UT_inc( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj TW, FLA_Obj W1, FLA_Obj B );
FLA_Error FLASH_Apply_pivots( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A );
FLA_Error FLASH_Eig_gest( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B );

// Not yet implemented:
FLA_Error FLASH_LQ_UT_inv( FLA_Obj A, FLA_Obj TW );
FLA_Error FLASH_LQ2_UT( FLA_Obj B, FLA_Obj C, FLA_Obj T );
