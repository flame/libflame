/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_scale_diagonals( FLA_Obj alpha, FLA_Obj A )
{
  FLA_Error r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Bidiag_UT_scale_diagonals_check( alpha, A );

  if ( FLA_Obj_length( A ) >= FLA_Obj_width( A ) )
    r_val = FLA_Bidiag_UT_u_scale_diagonals( alpha, A );
  else
    r_val = FLA_Bidiag_UT_l_scale_diagonals( alpha, A );

  return r_val;
}



FLA_Error FLA_Bidiag_UT_u_scale_diagonals( FLA_Obj alpha, FLA_Obj A )
{
  FLA_Datatype datatype;
  integer          n_A;
  integer          rs_A, cs_A;
  integer          i;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float*    buff_A     =  FLA_FLOAT_PTR( A );
      float*    buff_alpha =  FLA_FLOAT_PTR( alpha );
      for ( i = 0; i < n_A; ++i )
      {
        float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        float*    a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        integer       n_ahead  = n_A - i - 1;

        bl1_sscals( buff_alpha, alpha11 );

        if ( n_ahead > 0 )
          bl1_sscals( buff_alpha, a12t_l );
      }

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A     = FLA_DOUBLE_PTR( A );
      double*   buff_alpha = FLA_DOUBLE_PTR( alpha );
      for ( i = 0; i < n_A; ++i )
      {
        double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        double*   a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        integer       n_ahead  = n_A - i - 1;

        bl1_dscals( buff_alpha, alpha11 );

        if ( n_ahead > 0 )
          bl1_dscals( buff_alpha, a12t_l );
      }

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A     = FLA_COMPLEX_PTR( A );
      float*    buff_alpha = FLA_FLOAT_PTR( alpha );
      for ( i = 0; i < n_A; ++i )
      {
        scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        scomplex* a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        integer       n_ahead  = n_A - i - 1;

        bl1_csscals( buff_alpha, alpha11 );

        if ( n_ahead > 0 )
          bl1_csscals( buff_alpha, a12t_l );
      }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A     = FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_alpha = FLA_DOUBLE_PTR( alpha );
      for ( i = 0; i < n_A; ++i )
      {
        dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        dcomplex* a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        integer       n_ahead  = n_A - i - 1;

        bl1_zdscals( buff_alpha, alpha11 );

        if ( n_ahead > 0 )
          bl1_zdscals( buff_alpha, a12t_l );
      }

      break;
    }
  }

  return FLA_SUCCESS;
}

FLA_Error FLA_Bidiag_UT_l_scale_diagonals( FLA_Obj alpha, FLA_Obj A )
{
  FLA_Datatype datatype;
  integer          m_A;
  integer          rs_A, cs_A;
  integer          i;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float*    buff_A     = FLA_FLOAT_PTR( A );
      float*    buff_alpha = FLA_FLOAT_PTR( alpha );
      for ( i = 0; i < m_A; ++i )
      {
        float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        float*    a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        integer       m_ahead  = m_A - i - 1;

        bl1_sscals( buff_alpha, alpha11 );

        if ( m_ahead > 0 )
          bl1_sscals( buff_alpha, a21_t );
      }

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A     = FLA_DOUBLE_PTR( A );
      double*   buff_alpha = FLA_DOUBLE_PTR( alpha );
      for ( i = 0; i < m_A; ++i )
      {
        double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        double*   a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        integer       m_ahead  = m_A - i - 1;

        bl1_dscals( buff_alpha, alpha11 );

        if ( m_ahead > 0 )
          bl1_dscals( buff_alpha, a21_t );
      }

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A     = FLA_COMPLEX_PTR( A );
      float*    buff_alpha = FLA_FLOAT_PTR( alpha );
      for ( i = 0; i < m_A; ++i )
      {
        scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        scomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        integer       m_ahead  = m_A - i - 1;

        bl1_csscals( buff_alpha, alpha11 );

        if ( m_ahead > 0 )
          bl1_csscals( buff_alpha, a21_t );
      }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A     = FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_alpha = FLA_DOUBLE_PTR( alpha );
      for ( i = 0; i < m_A; ++i )
      {
        dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        dcomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        integer       m_ahead  = m_A - i - 1;

        bl1_zdscals( buff_alpha, alpha11 );

        if ( m_ahead > 0 )
          bl1_zdscals( buff_alpha, a21_t );
      }

      break;
    }
  }

  return FLA_SUCCESS;
}

