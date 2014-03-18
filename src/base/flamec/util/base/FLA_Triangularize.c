
#include "FLAME.h"

FLA_Error FLA_Triangularize( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  uplo1_t       blis_uplo;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Triangularize_check( uplo, diag, A );

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  // We have to toggle the uplo parameter because we will use it to specify
  // which triangle to zero out.
  if ( uplo == FLA_LOWER_TRIANGULAR ) uplo = FLA_UPPER_TRIANGULAR;
  else                                uplo = FLA_LOWER_TRIANGULAR;

  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_0 = ( float * ) FLA_FLOAT_PTR( FLA_ZERO );
    float *buff_1 = ( float * ) FLA_FLOAT_PTR( FLA_ONE );

    bl1_ssetmr( blis_uplo,
                m_A,
                n_A,
                buff_0,
                buff_A, rs_A, cs_A );

    if ( diag == FLA_UNIT_DIAG )
      bl1_ssetdiag( 0,
                    m_A,
                    n_A,
                    buff_1,
                    buff_A, rs_A, cs_A );
    else if ( diag == FLA_ZERO_DIAG )
      bl1_ssetdiag( 0,
                    m_A,
                    n_A,
                    buff_0,
                    buff_A, rs_A, cs_A );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_0 = ( double * ) FLA_DOUBLE_PTR( FLA_ZERO );
    double *buff_1 = ( double * ) FLA_DOUBLE_PTR( FLA_ONE );

    bl1_dsetmr( blis_uplo,
                m_A,
                n_A,
                buff_0,
                buff_A, rs_A, cs_A );

    if ( diag == FLA_UNIT_DIAG )
      bl1_dsetdiag( 0,
                    m_A,
                    n_A,
                    buff_1,
                    buff_A, rs_A, cs_A );
    else if ( diag == FLA_ZERO_DIAG )
      bl1_dsetdiag( 0,
                    m_A,
                    n_A,
                    buff_0,
                    buff_A, rs_A, cs_A );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_0 = ( scomplex * ) FLA_COMPLEX_PTR( FLA_ZERO );
    scomplex *buff_1 = ( scomplex * ) FLA_COMPLEX_PTR( FLA_ONE );

    bl1_csetmr( blis_uplo,
                m_A,
                n_A,
                buff_0,
                buff_A, rs_A, cs_A );

    if ( diag == FLA_UNIT_DIAG )
      bl1_csetdiag( 0,
                    m_A,
                    n_A,
                    buff_1,
                    buff_A, rs_A, cs_A );
    else if ( diag == FLA_ZERO_DIAG )
      bl1_csetdiag( 0,
                    m_A,
                    n_A,
                    buff_0,
                    buff_A, rs_A, cs_A );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_0 = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( FLA_ZERO );
    dcomplex *buff_1 = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );

    bl1_zsetmr( blis_uplo,
                m_A,
                n_A,
                buff_0,
                buff_A, rs_A, cs_A );

    if ( diag == FLA_UNIT_DIAG )
      bl1_zsetdiag( 0,
                    m_A,
                    n_A,
                    buff_1,
                    buff_A, rs_A, cs_A );
    else if ( diag == FLA_ZERO_DIAG )
      bl1_zsetdiag( 0,
                    m_A,
                    n_A,
                    buff_0,
                    buff_A, rs_A, cs_A );

    break;
  }

  }

  return FLA_SUCCESS;
}

