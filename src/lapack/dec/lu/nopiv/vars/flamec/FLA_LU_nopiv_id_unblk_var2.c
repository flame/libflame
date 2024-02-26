/*
 *  Copyright (c) 2021-2023 Advanced Micro Devices, Inc. All rights reserved.
 * */

#include "FLAME.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

/*******************************************************************************************
 This algorithm uses noncrecursive nonpivot based LU factorization using same logic as getrf
 except that it stops at  nfact iteration. This algorithm works better than the other
 in very small cases  satisfying the following  constraints,
 for sizes <=611 such that nfact is less than 51% of (m or n)
 nfact*(m or n) <=200 for nfact<=4
 nfact*(m or n) <=200+45 for nfact>4 and nfact<=8
 nfact*(m or n) <= 200+45+45 nfact>8 and nfact<=12 etc
*******************************************************************************************/

FLA_Error FLA_LU_nopiv_id_unblk_var2( integer m_A, integer n_A, double* A, integer nfact, integer rs_A, integer cs_A )
{
  double rminusone = bl1_dm1();
  double *Minusone = &rminusone;
  integer inc_x, inc_y, i, mdiff, ndiff;
  double alpha_inv;
  double *alpha;
  FLA_Error e_val = FLA_SUCCESS;

  // LU calculation
  // Algorithm - M rank1 updates can solve an LU system . For Proof see Sherman Morrison identity
  // The below algortithm performs nfact rank-1 updates to solve LU system partially
  for( i = 0; (i < nfact) && (m_A - i -1); i++ )                                      // i goes till nfact times , there is no point of scaling 0 hence m_A-i-1 check is present
  {
     alpha = A + i + i * cs_A;
     mdiff = m_A - i - 1;
     ndiff = n_A - i - 1;
     if( *alpha == 0.0 ) break;
     if( *alpha != 1.0 )
     {
        alpha_inv = 1.0 / *alpha;
        dscal_( &mdiff, &alpha_inv, (A + i + 1 + i * cs_A), &rs_A );                  // rank 1 update
     }
     inc_x = (i == m_A - 2) ? cs_A : rs_A;                                            // the vector will be 1 and equal to column stride for m_A - 2 aka row vector
     inc_y = cs_A;                                                                    // always equal to column stride
     dger_( &ndiff, &mdiff, Minusone, (A + i + 1 + i * cs_A), &inc_x, (A + i + (i + 1) * cs_A), &inc_y, (A + i + 1 + (i + 1) * cs_A), &cs_A ); // rank 1 update
  }
  // Singular check  and info population
  for( i = 0; i < nfact; i++ )
  {
     if( *(A + i + i * cs_A) == 0.0 )
     {
        e_val = i;
        return e_val;                                                                 // Early return if 0 found on LU factorized matrix's main diagonal,e_val contains index where 0 exist
     }
  }
  return e_val;
}
