/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_QR_UT_solve( FLA_Obj A, FLA_Obj TW, FLA_Obj B, FLA_Obj X )
{
  FLA_Obj W, Y;
  FLA_Obj AT, AB;
  FLA_Obj YT, YB;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_QR_UT_solve_check( A, TW, B, X );

  FLASH_Apply_Q_UT_create_workspace( TW, B, &W );

  FLASH_Obj_create_copy_of( FLA_NO_TRANSPOSE, B, &Y );

  FLASH_Apply_Q_UT( FLA_LEFT, FLA_CONJ_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                    A, TW, W, Y );

  FLA_Part_2x1( A,   &AT,
                     &AB,    FLA_Obj_width( A ), FLA_TOP );
  FLA_Part_2x1( Y,   &YT,
                     &YB,    FLA_Obj_width( A ), FLA_TOP );

  FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
              FLA_ONE, AT, YT );

  FLASH_Copy( YT, X );

  FLASH_Obj_free( &Y );
  FLASH_Obj_free( &W );

  return FLA_SUCCESS;
}

