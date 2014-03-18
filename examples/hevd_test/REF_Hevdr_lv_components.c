
#include "FLAME.h"

FLA_Error REF_Hevdr_lv_components( FLA_Obj A, FLA_Obj l, FLA_Obj Z,
                                   double* dtime_tred, double* dtime_tevd, double* dtime_appq )
/*
{
  return FLA_Hevdr_external( FLA_EVD_WITH_VECTORS, FLA_LOWER_TRIANGULAR, A, l, Z );
}
*/

{
  FLA_Datatype dt_A;
  FLA_Datatype dt_A_real;
  FLA_Uplo     uplo;
  dim_t        m_A, n_A;
  FLA_Obj      t, d, e, W, eT, eB;
  double       dtime_temp;

  uplo      = FLA_LOWER_TRIANGULAR;
  dt_A      = FLA_Obj_datatype( A );
  dt_A_real = FLA_Obj_datatype_proj_to_real( A );
  m_A       = FLA_Obj_length( A );
  n_A       = FLA_Obj_width( A );

  FLA_Obj_create( dt_A,      m_A,   1, 0, 0, &t );
  FLA_Obj_create( dt_A_real, m_A,   1, 0, 0, &d );
  FLA_Obj_create( dt_A_real, m_A,   1, 0, 0, &e );
  FLA_Obj_create( dt_A,      m_A, n_A, 0, 0, &W );

  FLA_Part_2x1( e,   &eT,
                     &eB,   m_A-1, FLA_TOP );

  dtime_temp = FLA_Clock();
  {
    // Reduce to tridiagonal form.
    FLA_Tridiag_blk_external( uplo, A, t );
    FLA_Tridiag_UT_extract_diagonals( uplo, A, d, eT );
  }
  *dtime_tred = FLA_Clock() - dtime_temp;


  dtime_temp = FLA_Clock();
  {
    // MRRR algorithm.
    FLA_Tevdr_external( FLA_EVD_WITH_VECTORS, d, e, l, W );
  }
  *dtime_tevd = FLA_Clock() - dtime_temp;


  dtime_temp = FLA_Clock();
  {
    // Apply Q.
    FLA_Tridiag_apply_Q_external( FLA_LEFT, uplo, FLA_NO_TRANSPOSE, A, t, W );
    FLA_Copy( W, Z );
  }
  *dtime_appq = FLA_Clock() - dtime_temp;


  // Sort eigenvalues and eigenvectors.
  FLA_Sort_evd( FLA_FORWARD, l, Z );

//FLA_Obj_show( "refr: d", l, "%22.10e", "" );
//FLA_Obj_show( "refr: A", A, "%8.1e", "" );

  FLA_Obj_free( &t );
  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
  FLA_Obj_free( &W );

  return FLA_SUCCESS;
}


