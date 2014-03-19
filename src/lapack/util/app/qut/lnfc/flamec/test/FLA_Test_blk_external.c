/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Test_blk_external( FLA_Side side, FLA_Trans trans, FLA_Store storev, FLA_Obj A, FLA_Obj t, FLA_Obj B )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  // int          m_A, n_A;
  int          m_B, n_B;
  int          cs_A;
  int          cs_B;
  int          k_t;
  int          lwork;
  char         blas_side;
  char         blas_trans;
  FLA_Obj      work_obj, d, e, tu, tv;

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  // m_A      = FLA_Obj_length( A );
  // n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  cs_B     = FLA_Obj_col_stride( B );

  k_t      = FLA_Obj_vector_dim( t );

  FLA_Param_map_flame_to_netlib_side( side, &blas_side );
  FLA_Param_map_flame_to_netlib_trans( trans, &blas_trans );

  if ( side == FLA_LEFT )
    lwork  = ( m_B + n_B ) * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );
  else
    lwork  = m_B * FLA_Query_blocksize( datatype, FLA_DIMENSION_MIN );

  FLA_Obj_create( datatype, lwork, 1, 0, 0, &work_obj );
  FLA_Obj_create( FLA_DOUBLE, n_B, 1, 0, 0, &d );
  FLA_Obj_create( FLA_DOUBLE, n_B, 1, 0, 0, &e );
  FLA_Obj_create( datatype,   n_B, 1, 0, 0, &tu );
  FLA_Obj_create( datatype,   n_B, 1, 0, 0, &tv );
  

  switch( datatype ){

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    double*   buff_d    = ( double   * ) FLA_DOUBLE_PTR( d );
    double*   buff_e    = ( double   * ) FLA_DOUBLE_PTR( e );
    dcomplex *buff_tu   = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( tu );
    dcomplex *buff_tv   = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( tv );
    dcomplex *buff_t    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( t );
    dcomplex *buff_B    = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
    dcomplex *buff_work = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( work_obj );

    zgebrd_( &m_B,
             &n_B,
             buff_A, &cs_A,
             buff_d,
             buff_e,
             buff_tu,
             buff_tv,
             buff_work, &lwork,
             &info );
    zunmqr_( &blas_side,
             &blas_trans,
             &m_B,
             &n_B,
             &k_t,
             buff_A, &cs_A,
             buff_tu,
             buff_B, &cs_B,
             buff_work, &lwork,
             &info );

    break;
  }

  }

  FLA_Obj_free( &work_obj );
  FLA_Obj_free( &d );
  FLA_Obj_free( &e );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

