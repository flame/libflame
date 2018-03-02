/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error REF_Hevdd_lv_components( FLA_Obj A, FLA_Obj l,
                                   double* dtime_tred, double* dtime_tevd, double* dtime_appq )
/*
{
  return FLA_Hevdd_external( FLA_EVD_WITH_VECTORS, FLA_LOWER_TRIANGULAR, A, l );
}
*/
{
  FLA_Datatype dt_A;
  FLA_Datatype dt_A_real;
  FLA_Uplo     uplo;
  dim_t        m_A, n_A;
  FLA_Obj      t, d, e, W;
  double       dtime_temp;

  uplo      = FLA_LOWER_TRIANGULAR;
  dt_A      = FLA_Obj_datatype( A );
  dt_A_real = FLA_Obj_datatype_proj_to_real( A );
  m_A       = FLA_Obj_length( A );
  n_A       = FLA_Obj_width( A );

  FLA_Obj_create( dt_A,      m_A,   1, 0, 0, &t );
  FLA_Obj_create( dt_A_real, m_A,   1, 0, 0, &d );
  FLA_Obj_create( dt_A_real, m_A-1, 1, 0, 0, &e );
  FLA_Obj_create( dt_A,      m_A, n_A, 0, 0, &W );


  dtime_temp = FLA_Clock();
  {
    // Reduce to tridiagonal form.
    FLA_Tridiag_blk_external( uplo, A, t );
    FLA_Tridiag_UT_extract_real_diagonals( uplo, A, d, e );
  }
  *dtime_tred = FLA_Clock() - dtime_temp;


  dtime_temp = FLA_Clock();
  {
    // Divide-and-conquor algorithm.
    FLA_Tevdd_external( FLA_EVD_OF_TRIDIAG_WITH_VECTORS, d, e, W );
  }
  *dtime_tevd = FLA_Clock() - dtime_temp;


  dtime_temp = FLA_Clock();
  {
    // Apply Q.
    FLA_Tridiag_apply_Q_external( FLA_LEFT, uplo, FLA_NO_TRANSPOSE, A, t, W );
    FLA_Copy( W, A );
  }
  *dtime_appq = FLA_Clock() - dtime_temp;


  // Copy eigenvalues to output vector.
  FLA_Copy( d, l );

  // Sort eigenvalues and eigenvectors.
  FLA_Sort_evd( FLA_FORWARD, l, A );

  FLA_Obj_free( &t );
  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
  FLA_Obj_free( &W );

  return FLA_SUCCESS;
}

