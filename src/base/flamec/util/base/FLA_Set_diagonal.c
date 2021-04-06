/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


#define alpha_to_delta( alpha, delta ) delta = alpha; 
#define delta_to_alpha( alpha, delta ) alpha = delta; 

#define FLA_OBJ_SET_DIAGONAL_VECTOR( mytype, mybuffer, equal_to, mydim ) \
  {                                                                     \
    mytype*   buff_A = mybuffer( A );                                   \
    mytype*   buff_d = mybuffer( d );                                   \
                                                                        \
    for ( i = 0; i < mydim; ++i ) {                                     \
      equal_to( *( buff_A + i*cs_A + i*rs_A ), *(buff_d + i*inc_d) )    \
    }                                                                   \
  }

FLA_Error FLA_Set_diagonal_vector( FLA_Obj A, FLA_Obj d )
{
  FLA_Datatype datatype;
  integer          i, m;
  integer          rs_A, cs_A;
  integer          inc_d;

  datatype = FLA_Obj_datatype( A );

  m        = FLA_Obj_min_dim( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_d    = FLA_Obj_vector_inc( d );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      FLA_OBJ_SET_DIAGONAL_VECTOR( float, FLA_FLOAT_PTR, alpha_to_delta, m );
      break;
    }
    case FLA_DOUBLE:
    {
      FLA_OBJ_SET_DIAGONAL_VECTOR( double, FLA_DOUBLE_PTR, alpha_to_delta, m );
      break;
    }

    case FLA_COMPLEX:
    {
      FLA_OBJ_SET_DIAGONAL_VECTOR( scomplex, FLA_COMPLEX_PTR, alpha_to_delta, m );
      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      FLA_OBJ_SET_DIAGONAL_VECTOR( dcomplex, FLA_DOUBLE_COMPLEX_PTR, alpha_to_delta, m );
      break;
    }
  }
  return FLA_SUCCESS;
}

FLA_Error FLA_Set_diagonal_matrix( FLA_Obj d, FLA_Obj A )
{
  FLA_Datatype datatype;
  integer          i, m;
  integer          rs_A, cs_A;
  integer          inc_d;

  datatype = FLA_Obj_datatype( A );

  m        = FLA_Obj_min_dim( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_d    = FLA_Obj_vector_inc( d );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      FLA_OBJ_SET_DIAGONAL_VECTOR( float, FLA_FLOAT_PTR, delta_to_alpha, m );
      break;
    }
    case FLA_DOUBLE:
    {
      FLA_OBJ_SET_DIAGONAL_VECTOR( double, FLA_DOUBLE_PTR, delta_to_alpha, m );
      break;
    }

    case FLA_COMPLEX:
    {
      FLA_OBJ_SET_DIAGONAL_VECTOR( scomplex, FLA_COMPLEX_PTR, delta_to_alpha, m );
      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      FLA_OBJ_SET_DIAGONAL_VECTOR( dcomplex, FLA_DOUBLE_COMPLEX_PTR, delta_to_alpha, m );
      break;
    }
  }
  return FLA_SUCCESS;
}


