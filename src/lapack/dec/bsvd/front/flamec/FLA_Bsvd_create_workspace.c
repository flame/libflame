/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bsvd_create_workspace( FLA_Obj d, FLA_Obj *G, FLA_Obj *H )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
  {  
    e_val = FLA_Check_real_object( d );
    FLA_Check_error_code( e_val );
    
    e_val = FLA_Check_nonconstant_object( d );
    FLA_Check_error_code( e_val );
  }

  // G and H stores the left and right Givens scalars.
  FLA_Datatype dt_comp  = FLA_Obj_datatype_proj_to_complex( d );
  dim_t        m_d      = FLA_Obj_vector_dim( d );
  dim_t        k_accum  = fla_min( 32, m_d );  

  if ( G != NULL ) FLA_Obj_create( dt_comp, m_d-1, k_accum, 0, 0, G );
  if ( H != NULL ) FLA_Obj_create( dt_comp, m_d-1, k_accum, 0, 0, H );

  return FLA_SUCCESS;
}
