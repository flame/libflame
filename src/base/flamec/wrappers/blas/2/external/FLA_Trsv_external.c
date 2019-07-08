/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLA_Trsv_external_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, FLA_Obj x ) 
{
  FLA_Datatype datatype;
  int          m_A;
  int          rs_A, cs_A;
  int          inc_x;
  uplo1_t       blis_uplo;
  trans1_t      blis_trans;
  diag1_t       blis_diag;

  if ( FLA_Check_error_level_ts(FLA_cntl_init_i) == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsv_check_ts( FLA_cntl_init_i, uplo, trans, diag, A, x );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype_ts( FLA_cntl_init_i, A );

  m_A      = FLA_Obj_length( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_x    = FLA_Obj_vector_inc( x );

  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );
  FLA_Param_map_flame_to_blis_trans( trans, &blis_trans );
  FLA_Param_map_flame_to_blis_diag( diag, &blis_diag );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_x = ( float * ) FLA_FLOAT_PTR( x );

    bl1_strsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_x = ( double * ) FLA_DOUBLE_PTR( x );

    bl1_dtrsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_x = ( scomplex * ) FLA_COMPLEX_PTR( x );

    bl1_ctrsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_x = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );

    bl1_ztrsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  }

  return FLA_SUCCESS;
}
#endif

FLA_Error FLA_Trsv_external( FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj A, FLA_Obj x ) 
{
  FLA_Datatype datatype;
  int          m_A;
  int          rs_A, cs_A;
  int          inc_x;
  uplo1_t       blis_uplo;
  trans1_t      blis_trans;
  diag1_t       blis_diag;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsv_check( uplo, trans, diag, A, x );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_x    = FLA_Obj_vector_inc( x );

  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );
  FLA_Param_map_flame_to_blis_trans( trans, &blis_trans );
  FLA_Param_map_flame_to_blis_diag( diag, &blis_diag );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_x = ( float * ) FLA_FLOAT_PTR( x );

    bl1_strsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_x = ( double * ) FLA_DOUBLE_PTR( x );

    bl1_dtrsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_x = ( scomplex * ) FLA_COMPLEX_PTR( x );

    bl1_ctrsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_x = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );

    bl1_ztrsv( blis_uplo,
               blis_trans,
               blis_diag,
               m_A,
               buff_A, rs_A, cs_A,
               buff_x, inc_x );

    break;
  }

  }

  return FLA_SUCCESS;
}

