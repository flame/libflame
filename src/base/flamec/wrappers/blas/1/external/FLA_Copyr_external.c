/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Copyr_external( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype dt_A;
  FLA_Datatype dt_B;
  integer          m_B, n_B;
  integer          rs_A, cs_A;
  integer          rs_B, cs_B;
  uplo1_t       blis_uplo;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Copyr_check( uplo, A, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  dt_A     = FLA_Obj_datatype( A );
  dt_B     = FLA_Obj_datatype( B );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );

  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );

  // If A is of type FLA_CONSTANT, then we have to proceed based on the
  // datatype of B.
  if      ( dt_A == FLA_CONSTANT )
  {
    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bl1_scopymr( blis_uplo,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bl1_dcopymr( blis_uplo,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bl1_ccopymr( blis_uplo,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bl1_zcopymr( blis_uplo,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
  }
/*
  else if ( dt_A == FLA_INT )
  {
    integer*      buff_A = ( integer * ) FLA_INT_PTR( A );
    integer*      buff_B = ( integer * ) FLA_INT_PTR( B );

    bl1_icopymr( blis_uplo,
                 m_B,
                 n_B,
                 buff_A, rs_A, cs_A,
                 buff_B, rs_B, cs_B );
  }
*/
  else if ( dt_A == FLA_FLOAT )
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );

    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bl1_scopymr( blis_uplo,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bl1_sdcopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bl1_sccopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bl1_szcopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
  }
  else if ( dt_A == FLA_DOUBLE )
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );

    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bl1_dscopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bl1_dcopymr( blis_uplo,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bl1_dccopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bl1_dzcopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
  }
  else if ( dt_A == FLA_COMPLEX )
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );

    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bl1_cscopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bl1_cdcopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bl1_ccopymr( blis_uplo,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bl1_czcopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
  }
  else if ( dt_A == FLA_DOUBLE_COMPLEX )
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bl1_zscopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bl1_zdcopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bl1_zccopymr( blis_uplo,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bl1_zcopymr( blis_uplo,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
  }
  
  return FLA_SUCCESS;
}

