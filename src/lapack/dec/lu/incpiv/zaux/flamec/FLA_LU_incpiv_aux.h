/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_SA_Apply_pivots( FLA_Obj C, FLA_Obj E, FLA_Obj p );
FLA_Error FLA_SA_LU_blk( FLA_Obj U,
                         FLA_Obj D, FLA_Obj p, FLA_Obj L, dim_t nb_alg );
FLA_Error FLA_SA_LU_unb( FLA_Obj U, 
                         FLA_Obj D, FLA_Obj p, FLA_Obj L );
FLA_Error FLA_SA_FS_blk( FLA_Obj L,
                         FLA_Obj D, FLA_Obj p, FLA_Obj C,
                                               FLA_Obj E, dim_t nb_alg );

FLA_Error FLASH_LU_incpiv_var1( FLA_Obj A, FLA_Obj p, FLA_Obj L, dim_t nb_alg, fla_lu_t* cntl );
FLA_Error FLASH_LU_incpiv_var2( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj U, dim_t nb_alg, fla_lu_t* cntl );
FLA_Error FLASH_Trsm_piv( FLA_Obj A, FLA_Obj B, FLA_Obj p, fla_trsm_t* cntl );
FLA_Error FLASH_SA_LU( FLA_Obj B, FLA_Obj C,
                       FLA_Obj D, FLA_Obj E, FLA_Obj p, FLA_Obj L, dim_t nb_alg, fla_lu_t* cntl );
FLA_Error FLASH_SA_FS( FLA_Obj L,
                       FLA_Obj D, FLA_Obj p, FLA_Obj C,
                                             FLA_Obj E, dim_t nb_alg, fla_gemm_t* cntl );

FLA_Error FLASH_FS_incpiv_aux1( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj b, dim_t nb_alg );
FLA_Error FLASH_FS_incpiv_aux2( FLA_Obj L,
                                FLA_Obj D, FLA_Obj p, FLA_Obj C,
                                                      FLA_Obj E, dim_t nb_alg );

