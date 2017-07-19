/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"



FLA_Error FLA_Obj_copy_view( FLA_Obj A, FLA_Obj* B )
{
  FLA_Obj  A_view;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_copy_view_check( A, B );

  // Set the m_inner and n_inner fields of a temporary copy of A.
  A_view         = A;
  A_view.m_inner = FLASH_Obj_scalar_length( A );
  A_view.n_inner = FLASH_Obj_scalar_width( A ); 

  // Copy the modified view into B. 
  *B = A_view; 

  return FLA_SUCCESS;
}



void FLA_Obj_extract_real_scalar( FLA_Obj alpha, double* alpha_value )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_real_scalar_check( alpha, alpha_value );

  if ( FLA_Obj_is_single_precision( alpha ) )
    *alpha_value = ( double ) *FLA_FLOAT_PTR( alpha );
  else
    *alpha_value = *FLA_DOUBLE_PTR( alpha );
}



void FLA_Obj_extract_complex_scalar( FLA_Obj alpha, dcomplex* alpha_value )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_complex_scalar_check( alpha, alpha_value );

  if ( FLA_Obj_is_single_precision( alpha ) )
  {
    scomplex temp = *FLA_COMPLEX_PTR( alpha );
    alpha_value->real = ( double ) temp.real;
    alpha_value->imag = ( double ) temp.imag;
  }
  else
    *alpha_value = *FLA_DOUBLE_COMPLEX_PTR( alpha );
}



void FLA_Obj_extract_real_part( FLA_Obj a, FLA_Obj b )
{
  FLA_Datatype datatype;
  int          m, inc_a, inc_b;
  
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_real_part_check( a, b );

  datatype = FLA_Obj_datatype( a );

  m        = FLA_Obj_vector_dim( a );
  inc_a    = FLA_Obj_vector_inc( a );
  inc_b    = FLA_Obj_vector_inc( b );

  switch ( datatype )
  {
  case FLA_FLOAT:
    {
      float* buff_a = FLA_FLOAT_PTR( a );
      float* buff_b = FLA_FLOAT_PTR( b );
      bl1_scopy( m, 
                 buff_a, inc_a, 
                 buff_b, inc_b );        
      break;
    }
  case FLA_DOUBLE: 
    {
      double* buff_a = FLA_DOUBLE_PTR( a );
      double* buff_b = FLA_DOUBLE_PTR( b );
      bl1_dcopy( m, 
                 buff_a, inc_a, 
                 buff_b, inc_b );        
      break;
    }
  case FLA_COMPLEX: 
    {
      float* buff_a = FLA_FLOAT_PTR( a ); 
      float* buff_b = FLA_FLOAT_PTR( b );
      bl1_scopy( m, 
                 buff_a, inc_a*2, 
                 buff_b, inc_b );        
      break;
    }
  case FLA_DOUBLE_COMPLEX: 
    {
      double* buff_a = FLA_DOUBLE_PTR( a );
      double* buff_b = FLA_DOUBLE_PTR( b );
      bl1_dcopy( m, 
                 buff_a, inc_a*2, 
                 buff_b, inc_b );        
      break;
    }
  }
}

void FLA_Obj_extract_imag_part( FLA_Obj a, FLA_Obj b )
{
  FLA_Datatype datatype;
  int          m, inc_a, inc_b;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_extract_imag_part_check( a, b );

  datatype = FLA_Obj_datatype( a );

  m        = FLA_Obj_vector_dim( a );
  inc_a    = FLA_Obj_vector_inc( a );
  inc_b    = FLA_Obj_vector_inc( b );

  switch ( datatype )
  {
  case FLA_FLOAT:
    {
      float* buff_b = FLA_FLOAT_PTR( b );
      float* buff_0 = FLA_FLOAT_PTR( FLA_ZERO );
      bl1_ssetv( m,
                 buff_0,
                 buff_b, inc_b );
      break;
    }
  case FLA_DOUBLE: 
    {
      double* buff_b = FLA_DOUBLE_PTR( b );
      double* buff_0 = FLA_DOUBLE_PTR( FLA_ZERO );
      bl1_dsetv( m,
                 buff_0,
                 buff_b, inc_b );
      break;
    }
  case FLA_COMPLEX: 
    {
      float* buff_a = FLA_FLOAT_PTR( a ); 
      float* buff_b = FLA_FLOAT_PTR( b );
      bl1_scopy( m, 
                 ++buff_a, inc_a*2, 
                 buff_b, inc_b );        
      break;
    }
  case FLA_DOUBLE_COMPLEX: 
    {
      double* buff_a = FLA_DOUBLE_PTR( a );
      double* buff_b = FLA_DOUBLE_PTR( b );
      bl1_dcopy( m, 
                 ++buff_a, inc_a*2, 
                 buff_b, inc_b );        
      break;
    }
  }
}


void FLA_Obj_set_real_part( FLA_Obj alpha, FLA_Obj B )
{
  dim_t m_B;
  dim_t n_B;
  dim_t rs_B;
  dim_t cs_B;
  dim_t i, j;

  m_B  = FLA_Obj_length( B );
  n_B  = FLA_Obj_width( B );
  rs_B = FLA_Obj_row_stride( B );
  cs_B = FLA_Obj_col_stride( B );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_set_real_part_check( alpha, B );

  if ( FLA_Obj_is_complex( B ) )
  {
    if      ( FLA_Obj_datatype( B ) == FLA_COMPLEX )
    {
      float*    buff_alpha = FLA_FLOAT_PTR( alpha );
      scomplex* buff_B     = FLA_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          scomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->real = *buff_alpha;
        }
      }
    }
    else if ( FLA_Obj_datatype( B ) == FLA_DOUBLE_COMPLEX )
    {
      double*   buff_alpha = FLA_DOUBLE_PTR( alpha );
      dcomplex* buff_B     = FLA_DOUBLE_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          dcomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->real = *buff_alpha;
        }
      }
    }
  }
}



void FLA_Obj_set_imag_part( FLA_Obj alpha, FLA_Obj B )
{
  dim_t m_B;
  dim_t n_B;
  dim_t rs_B;
  dim_t cs_B;
  dim_t i, j;

  m_B  = FLA_Obj_length( B );
  n_B  = FLA_Obj_width( B );
  rs_B = FLA_Obj_row_stride( B );
  cs_B = FLA_Obj_col_stride( B );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_set_imag_part_check( alpha, B );

  if ( FLA_Obj_is_complex( B ) )
  {
    if      ( FLA_Obj_datatype( B ) == FLA_COMPLEX )
    {
      float*    buff_alpha = FLA_FLOAT_PTR( alpha );
      scomplex* buff_B     = FLA_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          scomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->imag = *buff_alpha;
        }
      }
    }
    else if ( FLA_Obj_datatype( B ) == FLA_DOUBLE_COMPLEX )
    {
      double*   buff_alpha = FLA_DOUBLE_PTR( alpha );
      dcomplex* buff_B     = FLA_DOUBLE_COMPLEX_PTR( B );

      for ( j = 0; j < n_B; ++j )
      {
        for ( i = 0; i < m_B; ++i )
        {
          dcomplex* beta11 = buff_B + rs_B * i + cs_B * j;

          beta11->imag = *buff_alpha;
        }
      }
    }
  }
}



FLA_Error FLA_Obj_fshow( FILE* file, char *s1, FLA_Obj A, char *format, char *s2 )
{
  FLA_Datatype datatype;
  dim_t        i, j, m, n;
  dim_t        rs, cs;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_fshow_check( file, s1, A, format, s2 );

  datatype = FLA_Obj_datatype( A );
  m        = FLA_Obj_length( A );
  n        = FLA_Obj_width( A );
  rs       = FLA_Obj_row_stride( A );
  cs       = FLA_Obj_col_stride( A );

  fprintf( file, "%s\n", s1 );

  switch ( datatype ){

  case FLA_CONSTANT:
  {
    int*      consti = FLA_INT_PTR( A );
    float*    consts = FLA_FLOAT_PTR( A );
    double*   constd = FLA_DOUBLE_PTR( A );
    scomplex* constc = FLA_COMPLEX_PTR( A );
    dcomplex* constz = FLA_DOUBLE_COMPLEX_PTR( A );

    fprintf( file, "int      = %d\n", *consti );
    fprintf( file, "float    = %e\n", *consts );
    fprintf( file, "double   = %e\n", *constd );
    fprintf( file, "scomplex = %e + %e\n", constc->real, constc->imag );
    fprintf( file, "dcomplex = %e + %e\n", constz->real, constc->imag );

    break;
  }

  case FLA_FLOAT:
  {
    float *buffer = ( float * ) FLA_FLOAT_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf( file, format, buffer[ j*cs + i*rs ] );
        fprintf( file, " " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_DOUBLE:
  {
    double *buffer = ( double * ) FLA_DOUBLE_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf( file, format, buffer[ j*cs + i*rs ] );
        fprintf( file, " " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buffer = ( scomplex * ) FLA_COMPLEX_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        //fprintf( file, format, buffer[ j*cs + i*rs ].real, buffer[ j*cs + i*rs ].imag );
        //fprintf( file, " " );
		fprintf( file, format, buffer[ j*cs + i*rs ].real );
		fprintf( file, " + " );
		fprintf( file, format, buffer[ j*cs + i*rs ].imag );
		fprintf( file, "  " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buffer = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        //fprintf( file, format, buffer[ j*cs + i*rs ].real, buffer[ j*cs + i*rs ].imag );
        //fprintf( file, " " );
		fprintf( file, format, buffer[ j*cs + i*rs ].real );
		fprintf( file, " + " );
		fprintf( file, format, buffer[ j*cs + i*rs ].imag );
		fprintf( file, "  " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  case FLA_INT:
  {
    int *buffer = ( int * ) FLA_INT_PTR( A );

    for ( i = 0; i < m; i++ )
    {
      for ( j = 0; j < n; j++ )
      {
        fprintf( file, format, buffer[ j*cs + i*rs ] );
        fprintf( file, " " );
      }
      fprintf( file, "\n" );
    }

    break;
  }

  }

  fprintf( file, "%s\n", s2 );

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_show( char *s1, FLA_Obj A, char *format, char *s2 )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_show_check( s1, A, format, s2 );

  return FLA_Obj_fshow( stdout, s1, A, format, s2 );
}

