
#include "FLAME.h"

FLA_Error FLA_Norm_frob( FLA_Obj A, FLA_Obj norm )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Norm_frob_check( A, norm );

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );
 
 
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_norm  = ( float * ) FLA_FLOAT_PTR( norm );

    bl1_sfnorm( m_A,
                n_A,
                buff_A, rs_A, cs_A,
                buff_norm );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_norm  = ( double * ) FLA_DOUBLE_PTR( norm );

    bl1_dfnorm( m_A,
                n_A,
                buff_A, rs_A, cs_A,
                buff_norm );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
    float    *buff_norm  = ( float    * ) FLA_FLOAT_PTR( norm );

    bl1_cfnorm( m_A,
                n_A,
                buff_A, rs_A, cs_A,
                buff_norm );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    double   *buff_norm  = ( double   * ) FLA_DOUBLE_PTR( norm );

    bl1_zfnorm( m_A,
                n_A,
                buff_A, rs_A, cs_A,
                buff_norm );

    break;
  }

  }

  return FLA_SUCCESS;
}

