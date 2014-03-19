/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_QR_UT_inc( FLA_Obj A, FLA_Obj TW )
{
  FLA_Error r_val;

  if ( FLASH_Queue_stack_depth() == 0 )
    r_val = FLASH_QR_UT_inc_opt1( A, TW );
  else
    r_val = FLASH_QR_UT_inc_noopt( A, TW );

  return r_val;
}

