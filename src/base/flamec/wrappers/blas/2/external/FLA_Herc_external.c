/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Herc_external( FLA_Uplo uplo, FLA_Conj conj, FLA_Obj alpha, FLA_Obj x, FLA_Obj A )
{
  FLA_Datatype datatype;
  int          m_A;
  int          rs_A, cs_A;
  int          inc_x; 
  uplo1_t       blis_uplo;
  conj1_t       blis_conj;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Herc_check( uplo, conj, alpha, x, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_x    = FLA_Obj_vector_inc( x );

  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );
  FLA_Param_map_flame_to_blis_conj( conj, &blis_conj );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_x     = ( float * ) FLA_FLOAT_PTR( x );
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    bl1_ssyr( blis_uplo,
              m_A,
              buff_alpha,
              buff_x, inc_x,
              buff_A, rs_A, cs_A ); 
    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_x     = ( double * ) FLA_DOUBLE_PTR( x );
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    bl1_dsyr( blis_uplo,
              m_A,
              buff_alpha,
              buff_x, inc_x,
              buff_A, rs_A, cs_A ); 
    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_x     = ( scomplex * ) FLA_COMPLEX_PTR( x );
    float    *buff_alpha = ( float    * ) FLA_FLOAT_PTR( alpha );

    bl1_cher( blis_uplo,
              blis_conj,
              m_A,
              buff_alpha,
              buff_x, inc_x,
              buff_A, rs_A, cs_A ); 

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_x     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );
    double   *buff_alpha = ( double   * ) FLA_DOUBLE_PTR( alpha );

    bl1_zher( blis_uplo,
              blis_conj,
              m_A,
              buff_alpha,
              buff_x, inc_x,
              buff_A, rs_A, cs_A ); 

    break;
  }

  }

  return FLA_SUCCESS;
}

