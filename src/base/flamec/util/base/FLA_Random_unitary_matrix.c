/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Random_unitary_matrix( FLA_Obj A )
{
  FLA_Obj B, T;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Random_unitary_matrix_check( A );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &B );

  FLA_Random_matrix( B );

  FLA_QR_UT_create_T( B, &T );

  FLA_QR_UT( B, T );

  FLA_QR_UT_form_Q( B, T, A );
  //FLA_Apply_Q_UT_create_workspace( A, T, &W );
  //FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE, B, T, W, A );

  //FLA_Obj_free( &W );
  FLA_Obj_free( &T );
  FLA_Obj_free( &B );

/*
  FLA_Datatype datatype;
  FLA_Obj      v, tau;
  FLA_Obj      aT,
               AB;
  int          i, mn;
  int          k;

  datatype = FLA_Obj_datatype( A );
  mn       = FLA_Obj_length( A );
  k        = 1;

  FLA_Obj_create( datatype, mn-1, 1, 0, 0, &v );
  FLA_Obj_create( datatype, 1,    1, 0, 0, &tau );

  FLA_Obj_set_to_identity( A );

  FLA_Part_2x1( A,   &aT,
                     &AB,    1, FLA_TOP );

  for ( i = 0; i < k; ++i )
  {
    FLA_Random_matrix( tau );
    FLA_Random_matrix( v );

    FLA_Apply_H2_UT( FLA_LEFT, tau, v, aT,
                                       AB );
  }

  FLA_Obj_free( &tau );
  FLA_Obj_free( &v );
*/

  return FLA_SUCCESS;
}

