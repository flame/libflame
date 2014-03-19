/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_extract_real_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Error r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Bidiag_UT_extract_real_diagonals_check( A, d, e );

  if ( FLA_Obj_length( A ) >= FLA_Obj_width( A ) )
    r_val = FLA_Bidiag_UT_u_extract_real_diagonals( A, d, e );
  else
    r_val = FLA_Bidiag_UT_l_extract_real_diagonals( A, d, e );

  return r_val;
}

FLA_Error FLA_Bidiag_UT_u_extract_real_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Datatype datatype;
  int          n_A;
  int          rs_A, cs_A;
  int          inc_d;
  int          inc_e;
  int          i;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_d    = FLA_Obj_vector_inc( d );

  if ( n_A != 1 )
    inc_e  = FLA_Obj_vector_inc( e );
  else
    inc_e  = 0;

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float*    buff_A = FLA_FLOAT_PTR( A );
      float*    buff_d = FLA_FLOAT_PTR( d );
      float*    buff_e = ( n_A != 1 ? FLA_FLOAT_PTR( e ) : NULL );

      for ( i = 0; i < n_A; ++i )
      {
        float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        float*    a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        float*    delta1   = buff_d + (i  )*inc_d;
        float*    epsilon1 = buff_e + (i  )*inc_e;

        int       n_ahead  = n_A - i - 1;

        // delta1 = alpha11;
        *delta1 = *alpha11;

        // epsilon1 = a12t_l;
        if ( n_ahead > 0 )
          *epsilon1 = *a12t_l;
      }

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A = FLA_DOUBLE_PTR( A );
      double*   buff_d = FLA_DOUBLE_PTR( d );
      double*   buff_e = ( n_A != 1 ? FLA_DOUBLE_PTR( e ) : NULL );

      for ( i = 0; i < n_A; ++i )
      {
        double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        double*   a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        double*   delta1   = buff_d + (i  )*inc_d;
        double*   epsilon1 = buff_e + (i  )*inc_e;

        int       n_ahead  = n_A - i - 1;

        // delta1 = alpha11;
        *delta1 = *alpha11;

        // epsilon1 = a12t_l;
        if ( n_ahead > 0 )
          *epsilon1 = *a12t_l;
      }

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      float*    buff_d = FLA_FLOAT_PTR( d );
      float*    buff_e = ( n_A != 1 ? FLA_FLOAT_PTR( e ) : NULL );

      for ( i = 0; i < n_A; ++i )
      {
        scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        scomplex* a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        float*    delta1   = buff_d + (i  )*inc_d;
        float*    epsilon1 = buff_e + (i  )*inc_e;

        int       n_ahead  = n_A - i - 1;

        // delta1 = alpha11;
        *delta1 = alpha11->real;

        // epsilon1 = a12t_l;
        if ( n_ahead > 0 )
          *epsilon1 = a12t_l->real;
      }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_d = FLA_DOUBLE_PTR( d );
      double*   buff_e = ( n_A != 1 ? FLA_DOUBLE_PTR( e ) : NULL );

      for ( i = 0; i < n_A; ++i )
      {
        dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        dcomplex* a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        double*   delta1   = buff_d + (i  )*inc_d;
        double*   epsilon1 = buff_e + (i  )*inc_e;

        int       n_ahead  = n_A - i - 1;

        // delta1 = alpha11;
        *delta1 = alpha11->real;

        // epsilon1 = a12t_l;
        if ( n_ahead > 0 )
          *epsilon1 = a12t_l->real;
      }

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_l_extract_real_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Datatype datatype;
  int          m_A;
  int          rs_A, cs_A;
  int          inc_d;
  int          inc_e;
  int          i;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_d    = FLA_Obj_vector_inc( d );
  
  if ( m_A != 1 )
    inc_e  = FLA_Obj_vector_inc( e );
  else
    inc_e  = 0;

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float*    buff_A = FLA_FLOAT_PTR( A );
      float*    buff_d = FLA_FLOAT_PTR( d );
      float*    buff_e = ( m_A != 1 ? FLA_FLOAT_PTR( e ) : NULL );

      for ( i = 0; i < m_A; ++i )
      {
        float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        float*    a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        float*    delta1   = buff_d + (i  )*inc_d;
        float*    epsilon1 = buff_e + (i  )*inc_e;

        int       m_ahead  = m_A - i - 1;

        // delta1 = alpha11;
        *delta1 = *alpha11;

        // epsilon1 = a21_t;
        if ( m_ahead > 0 )
          *epsilon1 = *a21_t;
      }

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A = FLA_DOUBLE_PTR( A );
      double*   buff_d = FLA_DOUBLE_PTR( d );
      double*   buff_e = ( m_A != 1 ? FLA_DOUBLE_PTR( e ) : NULL );

      for ( i = 0; i < m_A; ++i )
      {
        double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        double*   a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        double*   delta1   = buff_d + (i  )*inc_d;
        double*   epsilon1 = buff_e + (i  )*inc_e;

        int       m_ahead  = m_A - i - 1;

        // delta1 = alpha11;
        *delta1 = *alpha11;

        // epsilon1 = a21_t;
        if ( m_ahead > 0 )
          *epsilon1 = *a21_t;
      }

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      float*    buff_d = FLA_FLOAT_PTR( d );
      float*    buff_e = ( m_A != 1 ? FLA_FLOAT_PTR( e ) : NULL );

      for ( i = 0; i < m_A; ++i )
      {
        scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        scomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        float*    delta1   = buff_d + (i  )*inc_d;
        float*    epsilon1 = buff_e + (i  )*inc_e;

        int       m_ahead  = m_A - i - 1;

        // delta1 = alpha11;
        *delta1 = alpha11->real;

        // epsilon1 = a21_t;
        if ( m_ahead > 0 )
          *epsilon1 = a21_t->real;
      }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_d = FLA_DOUBLE_PTR( d );
      double*   buff_e = ( m_A != 1 ? FLA_DOUBLE_PTR( e ) : NULL );

      for ( i = 0; i < m_A; ++i )
      {
        dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        dcomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        double*   delta1   = buff_d + (i  )*inc_d;
        double*   epsilon1 = buff_e + (i  )*inc_e;

        int       m_ahead  = m_A - i - 1;

        // delta1 = alpha11;
        *delta1 = alpha11->real;

        // epsilon1 = a21_t;
        if ( m_ahead > 0 )
          *epsilon1 = a21_t->real;
      }

      break;
    }
  }

  return FLA_SUCCESS;
}

