/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_eig_gest_t* fla_eig_gest_ix_cntl;
extern __thread fla_eig_gest_t* fla_eig_gest_nx_cntl;

FLA_Error FLA_Eig_gest( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  FLA_Obj   Y;
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Eig_gest_check( inv, uplo, A, B );

  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Y );

  // Invoke FLA_Eig_gest_internal() with the appropriate control tree.
  if ( inv == FLA_INVERSE )
    r_val = FLA_Eig_gest_internal( inv, uplo, A, Y, B, fla_eig_gest_ix_cntl );
  else
    r_val = FLA_Eig_gest_internal( inv, uplo, A, Y, B, fla_eig_gest_nx_cntl );

  FLA_Obj_free( &Y );

  return r_val;
}

