/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Copy_external( FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype dt_A;
  FLA_Datatype dt_B;
  integer          m_B, n_B;
  integer          rs_A, cs_A;
  integer          rs_B, cs_B;
  trans1_t      blis_trans;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Copy_check( A, B );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  dt_A     = FLA_Obj_datatype( A );
  dt_B     = FLA_Obj_datatype( B );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );

  if ( FLA_Obj_is_conformal_to( FLA_NO_TRANSPOSE, A, B ) )
    FLA_Param_map_flame_to_blis_trans( FLA_NO_TRANSPOSE, &blis_trans );
  else // if ( FLA_Obj_is_conformal_to( FLA_TRANSPOSE, A, B ) )
    FLA_Param_map_flame_to_blis_trans( FLA_TRANSPOSE, &blis_trans );

  // If A is of type FLA_CONSTANT, then we have to proceed based on the
  // datatype of B.
  if      ( dt_A == FLA_CONSTANT )
  {
    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bl1_scopymt( blis_trans,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bl1_dcopymt( blis_trans,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bl1_ccopymt( blis_trans,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bl1_zcopymt( blis_trans,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
  }
  else if ( dt_A == FLA_INT )
  {
    integer*      buff_A = ( integer * ) FLA_INT_PTR( A );
    integer*      buff_B = ( integer * ) FLA_INT_PTR( B );

    bl1_icopymt( blis_trans,
                 m_B,
                 n_B,
                 buff_A, rs_A, cs_A,
                 buff_B, rs_B, cs_B );
  }
  else if ( dt_A == FLA_FLOAT )
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );

    if      ( dt_B == FLA_FLOAT )
    {
      float *buff_B = ( float * ) FLA_FLOAT_PTR( B );
      
      bl1_scopymt( blis_trans,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bl1_sdcopymt( blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bl1_sccopymt( blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bl1_szcopymt( blis_trans,
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
      
      bl1_dscopymt( blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bl1_dcopymt( blis_trans,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bl1_dccopymt( blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bl1_dzcopymt( blis_trans,
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
      
      bl1_cscopymt( blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bl1_cdcopymt( blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bl1_ccopymt( blis_trans,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bl1_czcopymt( blis_trans,
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
      
      bl1_zscopymt( blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE )
    {
      double *buff_B = ( double * ) FLA_DOUBLE_PTR( B );
      
      bl1_zdcopymt( blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_COMPLEX )
    {
      scomplex *buff_B = ( scomplex * ) FLA_COMPLEX_PTR( B );
      
      bl1_zccopymt( blis_trans,
                    m_B,
                    n_B,
                    buff_A, rs_A, cs_A,
                    buff_B, rs_B, cs_B );
    }
    else if ( dt_B == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_B = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
      
      bl1_zcopymt( blis_trans,
                   m_B,
                   n_B,
                   buff_A, rs_A, cs_A,
                   buff_B, rs_B, cs_B );
    }
  }
  
  return FLA_SUCCESS;
}

