
#include "FLAME.h"

double FLA_Max_elemwise_diff( FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype datatype;
  dim_t        i, j;
  dim_t        m_A, n_A;
  dim_t        rs_A, cs_A;
  dim_t        rs_B, cs_B;
  double       diff;
  double       d_max = 0.0;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Max_elemwise_diff_check( A, B );

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );
  
  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_a = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_b = ( float * ) FLA_FLOAT_PTR( B );

    for( j = 0; j < n_A; j++ )
    {
      for( i = 0; i < m_A; i++ )
      {
        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ] - buff_b[ j*cs_B + i*rs_B ] );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);
      }
    }

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_a = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_b = ( double * ) FLA_DOUBLE_PTR( B );

    for( j = 0; j < n_A; j++ )
    {
      for( i = 0; i < m_A; i++ )
      {
        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ] - buff_b[ j*cs_B + i*rs_B ] );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);
      }
    }

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_a = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_b = ( scomplex * ) FLA_COMPLEX_PTR( B );

    for( j = 0; j < n_A; j++ )
    {
      for( i = 0; i < m_A; i++ )
      {
        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ].real - buff_b[ j*cs_B + i*rs_B ].real );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);

        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ].imag - buff_b[ j*cs_B + i*rs_B ].imag );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);
      }
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_a = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_b = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );

    for( j = 0; j < n_A; j++ )
    {
      for( i = 0; i < m_A; i++ )
      {
        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ].real - buff_b[ j*cs_B + i*rs_B ].real );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);

        diff = ( double ) ( buff_a[ j*cs_A + i*rs_A ].imag - buff_b[ j*cs_B + i*rs_B ].imag );

        if( fabs(diff) > d_max )
          d_max = fabs(diff);
      }
    }

    break;
  }

  }

  
  return d_max;
}

