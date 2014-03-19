/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_CAQR_UT_inc_solve( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW, FLA_Obj B, FLA_Obj X )
{
  FLA_Obj  W, Y;
  FLA_Obj  RT, RB;
  FLA_Obj  YT, YB;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_CAQR_UT_inc_solve_check( p, A, ATW, R, RTW, B, X );

  FLASH_Apply_CAQ_UT_inc_create_workspace( p, RTW, B, &W );

  FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, B, &Y );

  // Create a temporary hierarchical view of only the top n-by-n part of R in
  // case m > n so that RT captures the upper triangle. We do the same for Y
  // to ensure conformality.
  FLASH_Part_create_2x1( R,   &RT,    
                              &RB,    FLASH_Obj_scalar_width( R ), FLA_TOP );
  FLASH_Part_create_2x1( Y,   &YT,    
                              &YB,    FLASH_Obj_scalar_width( R ), FLA_TOP );

  FLASH_Apply_CAQ_UT_inc( p,
                          FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                          A, ATW, R, RTW, W, Y );

  FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, RT, YT );

  FLASH_Copy( YT, X );

  // Free the temporary hierarchical views.
  FLASH_Part_free_2x1( &RT,
                       &RB );
  FLASH_Part_free_2x1( &YT,
                       &YB );

  FLASH_Obj_free( &Y );
  FLASH_Obj_free( &W );

  return FLA_SUCCESS;
}

