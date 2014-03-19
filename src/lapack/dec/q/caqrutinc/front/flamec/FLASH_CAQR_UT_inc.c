/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_CAQR_UT_inc( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW )
{
  FLA_Error r_val;

  //if ( FLASH_Queue_stack_depth() == 0 )
  //  r_val = FLASH_CAQR_UT_inc_opt1( A, ATW, R, RTW );
  //else
    r_val = FLASH_CAQR_UT_inc_noopt( p, A, ATW, R, RTW );

  return r_val;
}

