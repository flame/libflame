/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_LU_piv_unb_external( FLA_Obj A, FLA_Obj p )
{
  FLA_Error    r_val = FLA_SUCCESS;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  int          info;
  FLA_Datatype datatype;
  int          m_A, n_A, cs_A;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_LU_piv_check( A, p );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );
    int   *buff_p = ( int   * ) FLA_INT_PTR( p );

    F77_sgetf2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_p,
                &info );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );
    int    *buff_p = ( int    * ) FLA_INT_PTR( p );

    F77_dgetf2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_p,
                &info );

    break;
  } 

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );
    int      *buff_p = ( int      * ) FLA_INT_PTR( p );

    F77_cgetf2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_p,
                &info );

    break;
  } 

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    int      *buff_p = ( int      * ) FLA_INT_PTR( p );

    F77_zgetf2( &m_A,
                &n_A,
                buff_A, &cs_A,
                buff_p,
                &info );

    break;
  } 

  }

  FLA_Shift_pivots_to( FLA_NATIVE_PIVOTS, p );

  // Convert to zero-based indexing, if an index was reported.
  if ( info > 0 ) r_val = info - 1;
  else            r_val = FLA_SUCCESS;

#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return r_val;
}

FLA_Error FLA_LU_piv_unb_ext( FLA_Obj A, FLA_Obj p )
{
  return FLA_LU_piv_unb_external( A, p );
}

