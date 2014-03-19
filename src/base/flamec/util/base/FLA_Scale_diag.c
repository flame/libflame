/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Scale_diag( FLA_Conj conj, FLA_Obj alpha, FLA_Obj A )
{
  FLA_Datatype datatype_A;
  FLA_Datatype datatype_alpha;
  dim_t        m_A, n_A;
  dim_t        rs_A, cs_A;
  conj1_t       blis_conj;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Scale_diag_check( conj, alpha, A );

  datatype_A     = FLA_Obj_datatype( A );
  datatype_alpha = FLA_Obj_datatype( alpha );
  m_A            = FLA_Obj_length( A );
  n_A            = FLA_Obj_width( A );
  rs_A           = FLA_Obj_row_stride( A );
  cs_A           = FLA_Obj_col_stride( A );

  FLA_Param_map_flame_to_blis_conj( conj, &blis_conj );

  switch( datatype_A ){

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    bl1_sscalediag( blis_conj,
                    0,
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

    bl1_dscalediag( blis_conj,
                    0,
                    m_A,
                    n_A,
                    buff_alpha,
                    buff_A, rs_A, cs_A );

    break;
  }

  case FLA_COMPLEX:
  {
    if ( datatype_alpha == FLA_COMPLEX )
    {
      scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
      scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );

      bl1_cscalediag( blis_conj,
                      0,
                      m_A,
                      n_A,
                      buff_alpha,
                      buff_A, rs_A, cs_A );
    }
    else
    {
     scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
      float    *buff_alpha = ( float    * ) FLA_FLOAT_PTR( alpha );

      bl1_csscalediag( blis_conj,
                       0,
                       m_A,
                       n_A,
                       buff_alpha,
                       buff_A, rs_A, cs_A );
    }

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    if ( datatype_alpha == FLA_DOUBLE_COMPLEX )
    {
      dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );

      bl1_zscalediag( blis_conj,
                      0,
                      m_A,
                      n_A,
                      buff_alpha,
                      buff_A, rs_A, cs_A );
    }
    else
    {
      dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
      double   *buff_alpha = ( double   * ) FLA_DOUBLE_PTR( alpha );

      bl1_zdscalediag( blis_conj,
                       0,
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

