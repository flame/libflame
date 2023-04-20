/*
 *  Copyright (c) 2021-2022 Advanced Micro Devices, Inc. All rights reserved.
 * */

#include "FLAME.h"

FLA_Error FLA_LU_nopiv_iz_blk_var1( integer m_A, integer n_A, FLA_Obj A, dcomplex* buff_A, integer nfact, integer rs_A, integer cs_A )
{

  void* FLA_memset( void* str, integer c, uinteger len );
  dcomplex* copy_A = (dcomplex*)FLA_malloc(m_A*n_A*sizeof(dcomplex));
  FLA_memset(copy_A,0,sizeof(copy_A));

  for(integer i=0;i<nfact;i++)
  {
    for(integer j=0;j<nfact;j++)
    {
      *(copy_A+i+j*cs_A) = *(buff_A+i+j*cs_A);
    }
  }

  FLA_Obj ATL, ATR, ABL, ABR;

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     nfact, nfact, FLA_TL );

  /*
   * the partition uses recursion based method to calculate LU incomplete factorization
   * ATL size : nfact*nfact
   * ABL size : (m-nfact)*(nfact)
   * ATR size : nfact*(n-nfact)
   * ABR size : (m-nfact)*(n-nfact)
   * */

  FLA_Error e_val = FLA_SUCCESS;

  // ATL = L1*U1
  e_val = FLA_LU_nopiv( ATL );                                             // Singular check, returns e_val = (i) where i is index on diagonal where value is 0
  if( e_val != FLA_SUCCESS )
  {
    for(integer i=0;i<nfact;i++)
    {
      for(integer j=0;j<nfact;j++)
      {
        *(buff_A+i+j*cs_A) = *(copy_A+i+j*cs_A);
      }
    }

    FLA_LU_nopiv_iz_unblk_var1(m_A,n_A,buff_A,e_val,rs_A,cs_A);            //restore original matrix and use unblocked variant 1  to recover last valid state$
    FLA_free(copy_A);
    return e_val;
  }

  FLA_free(copy_A);

  // U2 or new ATR = L1^-1 * ATR
  FLA_Trsm( FLA_LEFT, FLA_LOWER_TRIANGULAR,
            FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
            FLA_ONE, ATL, ATR );

  // L2 or new ABL = ABL * U1^-1
  FLA_Trsm( FLA_RIGHT, FLA_UPPER_TRIANGULAR,
            FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
            FLA_ONE, ATL, ABL );

  // ABR = ABR - new ABL * new ATR
  FLA_Gemm( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
            FLA_MINUS_ONE, ABL, ATR, FLA_ONE,
            ABR );

  return e_val;
}
