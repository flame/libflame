/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_eig_gest_t* fla_eig_gest_ix_cntl_leaf;
extern __thread fla_eig_gest_t* fla_eig_gest_nx_cntl_leaf;

FLA_Error FLA_Eig_gest_task( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
  fla_eig_gest_t* cntl_leaf;

  if ( inv == FLA_INVERSE )
    cntl_leaf = fla_eig_gest_ix_cntl_leaf;
  else
    cntl_leaf = fla_eig_gest_nx_cntl_leaf;

  return FLA_Eig_gest_internal( inv, uplo, A, Y, B,
                                cntl_leaf );
}

FLA_Error FLA_Eig_gest_il_task( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
  //return FLA_Eig_gest_unb_external( FLA_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
  return FLA_Eig_gest_internal( FLA_INVERSE, FLA_LOWER_TRIANGULAR, A, Y, B,
                                fla_eig_gest_ix_cntl_leaf );
}

FLA_Error FLA_Eig_gest_iu_task( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
  //return FLA_Eig_gest_unb_external( FLA_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
  return FLA_Eig_gest_internal( FLA_INVERSE, FLA_UPPER_TRIANGULAR, A, Y, B,
                                fla_eig_gest_ix_cntl_leaf );
}

FLA_Error FLA_Eig_gest_nl_task( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
  //return FLA_Eig_gest_unb_external( FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, B );
  return FLA_Eig_gest_internal( FLA_NO_INVERSE, FLA_LOWER_TRIANGULAR, A, Y, B,
                                fla_eig_gest_nx_cntl_leaf );
}

FLA_Error FLA_Eig_gest_nu_task( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl )
{
  //return FLA_Eig_gest_unb_external( FLA_NO_INVERSE, FLA_UPPER_TRIANGULAR, A, B );
  return FLA_Eig_gest_internal( FLA_NO_INVERSE, FLA_UPPER_TRIANGULAR, A, Y, B,
                                fla_eig_gest_nx_cntl_leaf );
}

