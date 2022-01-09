/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Scalr_external( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A )
{
  FLA_Datatype datatype, dt_alpha;
  integer          m_A, n_A;
  integer          rs_A, cs_A;
  uplo1_t       blis_uplo;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Scalr_check( uplo, alpha, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  if ( FLA_Obj_is_constant( alpha ) )
    dt_alpha = datatype;
  else
    dt_alpha = FLA_Obj_datatype( alpha );

  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    bl1_sscalmr( blis_uplo,
                 m_A,
                 n_A,
                 buff_alpha,
                 buff_A, rs_A, cs_A );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    bl1_dscalmr( blis_uplo,
                 m_A,
                 n_A,
                 buff_alpha,
                 buff_A, rs_A, cs_A );

    break;
  }

  case FLA_COMPLEX:
  {
    if ( dt_alpha == FLA_COMPLEX )
    {
      scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );
  
      bl1_cscalmr( blis_uplo,
                   m_A,
                   n_A,
                   buff_alpha,
                   buff_A, rs_A, cs_A );
    }
    else if ( dt_alpha == FLA_FLOAT )
    {
      scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
      float    *buff_alpha = ( float    * ) FLA_FLOAT_PTR( alpha );

      bl1_csscalmr( blis_uplo,
                    m_A,
                    n_A,
                    buff_alpha,
                    buff_A, rs_A, cs_A );
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    if ( dt_alpha == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
  
      bl1_zscalmr( blis_uplo,
                   m_A,
                   n_A,
                   buff_alpha,
                   buff_A, rs_A, cs_A );
    }
    else if ( dt_alpha == FLA_DOUBLE )
    {
      dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      double   *buff_alpha = ( double   * ) FLA_DOUBLE_PTR( alpha );
  
      bl1_zdscalmr( blis_uplo,
                    m_A,
                    n_A,
                    buff_alpha,
                    buff_A, rs_A, cs_A );
    }

    break;
  }

  }

  return FLA_SUCCESS;
}

