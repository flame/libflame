/*
 *  Copyright (c) 2021 Advanced Micro Devices, Inc. All rights reserved.
 * */

#include "FLAME.h"

/**************************************************************************************
 This method uses nonrecursive ger based method to calculate ATL,ABL,ATR,ABR implicitly
   * ATL size : nfact*nfact
   * ABL size : (m_A-nfact)*(nfact)
   * ATR size : nfact*(n_A-nfact)
   * ABR size : (m_A-nfact)*(n_A-nfact)
***************************************************************************************/

FLA_Error FLA_LU_nopiv_id_unblk_var1( int m_A, int n_A, double* A , int nfact, int rs_A, int cs_A )
{
  double rminusone = bl1_dm1();
  double *Minusone = &rminusone;
  double rone = bl1_d1();
  double *One = &rone;
  double rzero = bl1_d0();
  int inc_x, inc_y, i, diff, mnfactdiff, nnfactdiff;
  double alpha_inv;
  double *alpha;
  double dscalinv;
  mnfactdiff = m_A - nfact;
  nnfactdiff = n_A - nfact;
  FLA_Error e_val = FLA_SUCCESS;

  // LU calculation
  // ATL calculation if nfact!=1
  // ATL = L1*U1
  if( nfact != 1 )                         // If nfact == 1 ,no need to calculate ATL part
  {
     for( i = 0; i < nfact - 1; i++ )
     {
        alpha = A + i + i * cs_A;
        diff = nfact - i -1;
        if( *alpha != 1.0 )
        {
           alpha_inv = 1.0 / *alpha;
           dscal_( &diff, &alpha_inv, (A + i + 1 + i*cs_A), &rs_A );
        }
        inc_x = (i == nfact - 2) ? cs_A : rs_A;
        inc_y = cs_A;
        dger_( &diff, &diff, Minusone, (A + i + 1 + i * cs_A), &inc_x, (A + i + (i + 1) * cs_A), &inc_y, (A + i + 1 + (i + 1) * cs_A), &cs_A );
     }
  }

  // Singular check  and info population
  for( i = 0; i < nfact; i++ )
  {
     if( *(A + i + i * cs_A) == 0.0 )
     {
        e_val = i;
        return e_val;                               // Early return if 0 found on LU factorized matrix's main diagonal,e_val contains index where 0 exist
     }
  }
  // U2 or new ATR = L1^-1 * ATR
  if( mnfactdiff != 0 )                             // if size of ATR is 0 then trsm is invalid
  {
     if( nfact != 1 )                               // for nfact =1 , there is no need to calculate ATR part
     {
        dtrsm_( "L", "L", "N", "U", &nfact, &mnfactdiff, One, A, &cs_A, (A + cs_A * nfact), &cs_A );
     }
  }

  // L2 or new ABL calculation = U1^-1 * ABL
  if( nnfactdiff != 0 )                             // base invalid cases when size of  ABL is 0, nfact=1 case handled seperately
  {
     if( nfact != 1 )
     {
        dtrsm_( "R", "U", "N", "N", &nnfactdiff, &nfact, One, A, &cs_A, (A + nfact), &cs_A );
     }
     else                                           // if nfact==1
     {
        dscalinv = 1 / *A;
        dscal_( &nnfactdiff, &dscalinv, (A + 1), &rs_A );
     }
  }

  // ABR = ABR - new ABL * new ATR
  if( (mnfactdiff != 0) && (nnfactdiff != 0) )      // base invalid cases check
  {
     dgemm_( "N", "N", &nnfactdiff, &mnfactdiff,  &nfact, Minusone, (A + nfact), &cs_A, (A + nfact * cs_A), &cs_A, One, (A + nfact + nfact * cs_A), &cs_A );
  }
  return e_val;
}
