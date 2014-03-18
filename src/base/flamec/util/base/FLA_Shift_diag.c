
#include "FLAME.h"

FLA_Error FLA_Shift_diag( FLA_Conj conj, FLA_Obj sigma, FLA_Obj A )
{
  FLA_Datatype datatype_A;
  FLA_Datatype datatype_sigma;
  dim_t        m_A, n_A;
  dim_t        rs_A, cs_A;
  conj1_t       blis_conj;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Shift_diag_check( conj, sigma, A );

  datatype_A     = FLA_Obj_datatype( A );
  datatype_sigma = FLA_Obj_datatype( sigma );
  m_A            = FLA_Obj_length( A );
  n_A            = FLA_Obj_width( A );
  rs_A           = FLA_Obj_row_stride( A );
  cs_A           = FLA_Obj_col_stride( A );

  FLA_Param_map_flame_to_blis_conj( conj, &blis_conj );

  switch( datatype_A ){

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_sigma = ( float * ) FLA_FLOAT_PTR( sigma );

    bl1_sshiftdiag( blis_conj,
                    0,
                    m_A,
                    n_A,
                    buff_sigma,
                    buff_A, rs_A, cs_A );

    break;
  }
  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_sigma = ( double * ) FLA_DOUBLE_PTR( sigma );

    bl1_dshiftdiag( blis_conj,
                    0,
                    m_A,
                    n_A,
                    buff_sigma,
                    buff_A, rs_A, cs_A );

    break;
  }

  case FLA_COMPLEX:
  {
    if ( datatype_sigma == FLA_COMPLEX )
    {
      scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex *buff_sigma = ( scomplex * ) FLA_COMPLEX_PTR( sigma );

      bl1_cshiftdiag( blis_conj,
                      0,
                      m_A,
                      n_A,
                      buff_sigma,
                      buff_A, rs_A, cs_A );
    }
    else
    {
     scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
      float    *buff_sigma = ( float    * ) FLA_FLOAT_PTR( sigma );

      bl1_csshiftdiag( blis_conj,
                       0,
                       m_A,
                       n_A,
                       buff_sigma,
                       buff_A, rs_A, cs_A );
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    if ( datatype_sigma == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex *buff_sigma = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( sigma );

      bl1_zshiftdiag( blis_conj,
                      0,
                      m_A,
                      n_A,
                      buff_sigma,
                      buff_A, rs_A, cs_A );
    }
    else
    {
      dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      double   *buff_sigma = ( double   * ) FLA_DOUBLE_PTR( sigma );

      bl1_zdshiftdiag( blis_conj,
                       0,
                       m_A,
                       n_A,
                       buff_sigma,
                       buff_A, rs_A, cs_A );
    }

    break;
  }

  }

  return FLA_SUCCESS;
}

