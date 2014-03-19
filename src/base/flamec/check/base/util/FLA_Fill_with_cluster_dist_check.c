/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fill_with_cluster_dist_check( FLA_Obj n_clusters, FLA_Obj cluster_width, FLA_Obj x )
{
  FLA_Error e_val;

  e_val = FLA_Check_floating_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_nonconstant_object( x );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_int_object( n_clusters );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_real_object( cluster_width );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_identical_object_precision( cluster_width, x );
  FLA_Check_error_code( e_val );
  
  e_val = FLA_Check_if_scalar( n_clusters );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_scalar( cluster_width );
  FLA_Check_error_code( e_val );

  e_val = FLA_Check_if_vector( x );
  FLA_Check_error_code( e_val );

  return FLA_SUCCESS;
}

