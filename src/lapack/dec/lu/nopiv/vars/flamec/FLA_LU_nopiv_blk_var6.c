/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.Â  All rights reserved.
* */

#include "FLAME.h"

FLA_Error FLA_LU_nopiv_blk_var6( FLA_Obj A , dim_t nfact)
{
  FLA_Obj ATL, ATR, ABL, ABR;

  dim_t b;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     nfact, nfact, FLA_TL );

  /*
   * the partition is done so that we can efficiently solve it
   * variant 6 is directly called to remove recursion, however FLA_LU_nopiv will use
   * recursion to solve nfactxnfact matrix independently
   * ATL size : nfact*nfact
   * ABL size : (m-nfact)*(nfact)
   * ATR size : nfact*(n-nfact)
   * ABR size : (m-nfact)*(n-nfact)
   * */

  /*------------------------------------------------------------*/

    // ATL = L1*U1
    FLA_LU_nopiv( ATL );

    // U2 or new ATR = L1^-1 * ATR
    FLA_Trsm( FLA_LEFT, FLA_LOWER_TRIANGULAR,
              FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
              FLA_ONE, ATL, ATR
             );

    // L2 or new ABL = ABL * U1^-1
    FLA_Trsm( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
                       FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
                       FLA_ONE, ATL, ABL
            );

    // ABR = ABR - new ABL * new ATR
    FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
               FLA_MINUS_ONE, ABL, ATR, FLA_ONE, ABR
             );
  return FLA_SUCCESS;
}
