
#include "FLAME.h"

FLA_Error FLA_Add_to_diag( void* diag_value, FLA_Obj A )
{
  FLA_Datatype datatype;
  dim_t        j, min_m_n;
  dim_t        rs, cs;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Add_to_diag_check( diag_value, A );

  datatype = FLA_Obj_datatype( A );
  min_m_n  = FLA_Obj_min_dim( A );
  rs       = FLA_Obj_row_stride( A );
  cs       = FLA_Obj_col_stride( A );

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A    = ( float * ) FLA_FLOAT_PTR( A );
    float *value_ptr = ( float * ) diag_value;

    for ( j = 0; j < min_m_n; j++ )
      buff_A[ j*cs + j*rs ] += *value_ptr;

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A    = ( double * ) FLA_DOUBLE_PTR( A );
    double *value_ptr = ( double * ) diag_value;

    for ( j = 0; j < min_m_n; j++ )
      buff_A[ j*cs + j*rs ] += *value_ptr;

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A    = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *value_ptr = ( scomplex * ) diag_value;

    for ( j = 0; j < min_m_n; j++ )
    {
      buff_A[ j*cs + j*rs ].real += value_ptr->real;
      buff_A[ j*cs + j*rs ].imag += value_ptr->imag;
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *value_ptr = ( dcomplex * ) diag_value;

    for ( j = 0; j < min_m_n; j++ )
    {
      buff_A[ j*cs + j*rs ].real += value_ptr->real;
      buff_A[ j*cs + j*rs ].imag += value_ptr->imag;
    }

    break;
  }

  }

  return FLA_SUCCESS;
}

