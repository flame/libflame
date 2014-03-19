/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error REF_Svdd_uv_components( FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V,
                                  double* dtime_bred, double* dtime_bsvd, double* dtime_appq,
                                  double* dtime_qrfa, double* dtime_gemm )
/*
{
  *dtime_bred = 1;
  *dtime_bsvd = 1;
  *dtime_appq = 1;
  *dtime_qrfa = 1;
  *dtime_gemm = 1;

  return FLA_Svdd_external( FLA_SVD_VECTORS_ALL, A, s, U, V );
}
*/

{
  FLA_Datatype dt_A;
  FLA_Datatype dt_A_real;
  dim_t        m_A, n_A;
  dim_t        min_m_n;
  FLA_Obj      tq, tu, tv, d, e, Ur, Vr, W;
  FLA_Obj      eT, epsilonB;
  FLA_Uplo     uplo = FLA_UPPER_TRIANGULAR;
  double       crossover_ratio = 16.0 / 10.0;
  double       dtime_temp;

  dt_A      = FLA_Obj_datatype( A );
  dt_A_real = FLA_Obj_datatype_proj_to_real( A );
  m_A       = FLA_Obj_length( A );
  n_A       = FLA_Obj_width( A );

  min_m_n   = FLA_Obj_min_dim( A );

  FLA_Obj_create( dt_A,      min_m_n, 1,   0, 0, &tq );
  FLA_Obj_create( dt_A,      min_m_n, 1,   0, 0, &tu );
  FLA_Obj_create( dt_A,      min_m_n, 1,   0, 0, &tv );
  FLA_Obj_create( dt_A_real, min_m_n, 1,   0, 0, &d );
  FLA_Obj_create( dt_A_real, min_m_n, 1,   0, 0, &e );
  FLA_Obj_create( dt_A_real, n_A,     n_A, 0, 0, &Ur );
  FLA_Obj_create( dt_A_real, n_A,     n_A, 0, 0, &Vr );


  FLA_Part_2x1( e,   &eT,
                     &epsilonB,    1, FLA_BOTTOM );

  if ( m_A >= n_A )
  {
    if ( m_A < crossover_ratio * n_A )
    {
      dtime_temp = FLA_Clock();
      {
        // Reduce to bidiagonal form.
        FLA_Bidiag_blk_external( A, tu, tv );
        FLA_Bidiag_UT_extract_diagonals( A, d, eT );
      }
      *dtime_bred = FLA_Clock() - dtime_temp;


      dtime_temp = FLA_Clock();
      {
        // Divide-and-conquor algorithm.
        FLA_Bsvdd_external( uplo, d, e, Ur, Vr );
      }
      *dtime_bsvd = FLA_Clock() - dtime_temp;


      dtime_temp = FLA_Clock();
      {
        // Form U.
        FLA_Copy_external( Ur, U );
        FLA_Bidiag_apply_U_external( FLA_LEFT, FLA_NO_TRANSPOSE, A, tu, U );

        // Form V.
        FLA_Copy_external( Vr, V );
        FLA_Bidiag_apply_V_external( FLA_RIGHT, FLA_CONJ_TRANSPOSE, A, tv, V );
      }
      *dtime_appq = FLA_Clock() - dtime_temp;


      *dtime_qrfa = 0.0;
      *dtime_gemm = 0.0;
    }
    else
    {
      FLA_Obj AT,
              AB;
      FLA_Obj UL, UR;

      FLA_Part_2x1( A,   &AT,
                         &AB,        n_A, FLA_TOP );
      FLA_Part_1x2( U,   &UL, &UR,   n_A, FLA_LEFT );

      // Create a temporary n-by-n matrix R.
      FLA_Obj_create( dt_A, n_A, n_A, 0, 0, &W );

      dtime_temp = FLA_Clock();
      {
        // Perform a QR factorization.
        FLA_QR_blk_external( A, tq );
        FLA_Copyr_external( FLA_LOWER_TRIANGULAR, A, UL );
        FLA_Setr( FLA_LOWER_TRIANGULAR, FLA_ZERO, A );
      }
      *dtime_qrfa = FLA_Clock() - dtime_temp;


      dtime_temp = FLA_Clock();
      {
        // Form Q.
        FLA_QR_form_Q_external( U, tq );
      }
      *dtime_appq = FLA_Clock() - dtime_temp;


      dtime_temp = FLA_Clock();
      {
        // Reduce R to bidiagonal form.
        FLA_Bidiag_blk_external( AT, tu, tv );
        FLA_Bidiag_UT_extract_diagonals( A, d, eT );
      }
      *dtime_bred = FLA_Clock() - dtime_temp;


      dtime_temp = FLA_Clock();
      {
        // Divide-and-conquor algorithm.
        FLA_Bsvdd_external( uplo, d, e, Ur, Vr );
      }
      *dtime_bsvd = FLA_Clock() - dtime_temp;


      dtime_temp = FLA_Clock();
      {
        // Form U in W.
        FLA_Copy_external( Ur, W );
        FLA_Bidiag_apply_U_external( FLA_LEFT, FLA_NO_TRANSPOSE, AT, tu, W );

        // Form V.
        FLA_Copy_external( Vr, V );
        FLA_Bidiag_apply_V_external( FLA_RIGHT, FLA_CONJ_TRANSPOSE, AT, tv, V );
      }
      *dtime_appq += FLA_Clock() - dtime_temp;


      dtime_temp = FLA_Clock();
      {
        // Multiply R into U, storing the result in A and then copying
        // back to U.
        FLA_Gemm_external( FLA_NO_TRANSPOSE, FLA_NO_TRANSPOSE,
                           FLA_ONE, UL, W, FLA_ZERO, A );
        FLA_Copy( A, UL );
      }
      *dtime_gemm = FLA_Clock() - dtime_temp;


      // Free R.
      FLA_Obj_free( &W );
    }
  }
  else
  {
    FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
  }

  // Copy singular values to output vector.
  FLA_Copy( d, s );

  // Sort singular values and vectors.
  FLA_Sort_svd( FLA_BACKWARD, s, U, V );

  FLA_Obj_free( &tq );
  FLA_Obj_free( &tu );
  FLA_Obj_free( &tv );
  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
  FLA_Obj_free( &Ur );
  FLA_Obj_free( &Vr );

  return FLA_SUCCESS;
}

