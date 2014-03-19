/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_pivots_unb_external( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A )
{
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  int          n_A, cs_A;
  int          m_p;
  int          inc_p;
  int*         buff_p;
  int          k1_1, k2_1;
  int*         pivots_lapack;
  int          i;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Apply_pivots_check( side, trans, p, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_p    = FLA_Obj_vector_inc( p );
  m_p      = FLA_Obj_vector_dim( p );

  buff_p   = FLA_INT_PTR( p );

  // Use one-based indices for LAPACK.
  k1_1     = 1;
  k2_1     = m_p;

  // Translate FLAME pivot indices to LAPACK-compatible indices. It is
  // important to note that this conversion, unlike the one done by
  // FLA_Shift_pivots_to(), is NOT in-place, but rather done separately
  // in a temporary buffer.
#ifdef FLA_ENABLE_WINDOWS_BUILD
  pivots_lapack = ( int * ) _alloca( m_p * sizeof( int ) );
#else
  pivots_lapack = ( int * )  alloca( m_p * sizeof( int ) );
#endif

  for ( i = 0; i < m_p; i++ )
  {
    pivots_lapack[ i ] = buff_p[ i ] + i + 1;
  }

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_A = ( float * ) FLA_FLOAT_PTR( A );

    F77_slaswp( &n_A,
                buff_A, &cs_A,
                &k1_1, 
                &k2_1,
                pivots_lapack,
                &inc_p );
    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_A = ( double * ) FLA_DOUBLE_PTR( A );

    F77_dlaswp( &n_A,
                buff_A, &cs_A,
                &k1_1, 
                &k2_1,
                pivots_lapack,
                &inc_p );
    break;
  }

  case FLA_COMPLEX:
  {
    scomplex* buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );

    F77_claswp( &n_A,
                buff_A, &cs_A,
                &k1_1, 
                &k2_1,
                pivots_lapack,
                &inc_p );
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex* buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

    F77_zlaswp( &n_A,
                buff_A, &cs_A,
                &k1_1, 
                &k2_1,
                pivots_lapack,
                &inc_p );
    break;
  }

  }
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return FLA_SUCCESS;
}

FLA_Error FLA_Apply_pivots_ln_unb_ext( FLA_Obj p, FLA_Obj A )
{
  return FLA_Apply_pivots_unb_external( FLA_LEFT, FLA_NO_TRANSPOSE, p, A );
}

