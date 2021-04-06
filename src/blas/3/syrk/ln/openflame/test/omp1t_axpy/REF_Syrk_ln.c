/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error REF_Syrk_ln( FLA_Obj A, FLA_Obj C )
{
  FLA_Datatype datatype;
  integer          k, m, ldim_A, ldim_C;

  datatype = FLA_Obj_datatype( A );
  ldim_A   = FLA_Obj_ldim( A );
  ldim_C   = FLA_Obj_ldim( C );
  k        = FLA_Obj_width( A );
  m        = FLA_Obj_length( A );
  
  switch( datatype ){
    case FLA_DOUBLE:
    {
      double *buff_A, *buff_C, d_one=1.0;

      buff_A = ( double * ) FLA_Obj_buffer_at_view( A );
      buff_C = ( double * ) FLA_Obj_buffer_at_view( C );
    
      dsyrk_( "L", "N", &m, &k,
              &d_one, buff_A, &ldim_A, &d_one, buff_C, &ldim_C );
    } break;
  }
  
  return 0;
}

