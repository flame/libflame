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

FLA_Error FLA_LU_nopiv_is_unblk_var1( integer m_A, integer n_A, float* A , integer nfact, integer rs_A, integer cs_A )
{
  float rminusone = bl1_sm1();
  float *Minusone = &rminusone;
  float rone = bl1_s1();
  float *One = &rone;
  integer inc_x, inc_y, i, diff, tr_m, tr_n, tr_nfe, tr_ne;
  float alpha_inv;
  float *alpha;
  float sscalinv;
  tr_m = m_A - nfact;
  tr_n = n_A - nfact;
  FLA_Error e_val = FLA_SUCCESS;

  // LU calculation
  // ATL calculation if nfact!=1
  // ATL = L1*U1
  for( i = 0; i < nfact - 1; i++ )
  {
     alpha = A + i + i * cs_A;
     diff = nfact - i -1;
     if( *alpha != 1.0 )
     {
        alpha_inv = 1.0 / *alpha;
        sscal_( &diff, &alpha_inv, (A + i + 1 + i*cs_A), &rs_A );
     }
     inc_x = (i == nfact - 2) ? cs_A : rs_A;
     inc_y = cs_A;
     sger_( &diff, &diff, Minusone, (A + i + 1 + i * cs_A), &inc_x, (A + i + (i + 1) * cs_A), &inc_y, (A + i + 1 + (i + 1) * cs_A), &cs_A );
     if( *alpha == 0.0 )
     {
        e_val = i;
        break;
     }
  }

  if( e_val == FLA_SUCCESS )
  {
    // U2 or new ATR = L1^-1 * ATR
    if( tr_m != 0 )                             // if size of ATR is 0 then trsm is invalid
    {
       if( nfact != 1 )                               // for nfact =1 , there is no need to calculate ATR part
       {
          strsm_( "L", "L", "N", "U", &nfact, &tr_m, One, A, &cs_A, (A + cs_A * nfact), &cs_A );
       }
    }

    // L2 or new ABL calculation = U1^-1 * ABL
    if( tr_n != 0 )                             // base invalid cases when size of  ABL is 0, nfact=1 case handled seperately
    {
       if( nfact != 1 )
       {
          strsm_( "R", "U", "N", "N", &tr_n, &nfact, One, A, &cs_A, (A + nfact), &cs_A );
       }
       else                                           // if nfact==1
       {
          sscalinv = 1 / *A;
          sscal_( &tr_n, &sscalinv, (A + 1), &rs_A );
       }
    }

    // ABR = ABR - new ABL * new ATR
    if( (tr_m != 0) && (tr_n != 0) )      // base invalid cases check
    {
       sgemm_( "N", "N", &tr_n, &tr_m,  &nfact, Minusone, (A + nfact), &cs_A, (A + nfact * cs_A), &cs_A, One, (A + nfact + nfact * cs_A), &cs_A );
    }
  }
  else                                              // setting matrix to last known valid state when encountering singular matrix
  {
    // U2 or new ATR = L1^-1 * ATR                    // ATR size to be updated is e_val x (m_A - nfact)
    if( tr_m != 0 )                             // if size of ATR is 0 then trsm is invalid
    {
      strsm_( "L", "L", "N", "U", (integer *) &e_val, &tr_m, One, A, &cs_A, (A + cs_A * nfact), &cs_A );
    }

    // L2 valid ABL calculation = U1^-1 * valid ABL   // ABL size to be updated is (n_A - nfact) * e_val
    if( tr_n != 0 )                             // base invalid cases when size of  ABL is 0, nfact=1 case handled seperately
    {
       strsm_( "R", "U", "N", "N", &tr_n, (integer *) &e_val, One, A, &cs_A, (A + nfact), &cs_A );
    }


    // new ABR1 = ABR1 - valid ABL * valid updated nfact  block             // ABR update
                                                    // The ABR is now a compound shape of 2 rectanges so 2 updates will happen using dgemm_ ABR = ABR1 + ABR2

    if( ( tr_n != 0 ) )                         // base invalid cases check
    {
       tr_nfe = nfact - e_val;
       sgemm_( "N", "N", &tr_n, &tr_nfe, (integer *) &e_val, Minusone, (A + nfact), &cs_A, (A + e_val * cs_A), &cs_A, One, (A + nfact + e_val * cs_A), &cs_A );
    }

  // new ABR2 = ABR2 - valid ATR * (valid updated ATL's ABL part + new ABL)
    if( ( tr_m != 0 ) )
    {
       tr_ne = n_A - e_val;
       sgemm_( "N", "N", &tr_ne, &tr_m, (integer *) &e_val, Minusone, (A + e_val), &cs_A, (A + nfact * cs_A), &cs_A, One, (A + e_val + nfact * cs_A), &cs_A );
    }
  }

  return e_val;
}
