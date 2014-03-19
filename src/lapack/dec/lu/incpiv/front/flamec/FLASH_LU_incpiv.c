/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_LU_incpiv( FLA_Obj A, FLA_Obj p, FLA_Obj L )
{
   FLA_Error r_val;

   // Check parameters.
   if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
      FLA_LU_incpiv_check( A, p, L );

   // *** The current LU_incpiv algorithm implemented assumes that
   // the matrix has a hierarchical depth of 1. We check for that here, because
   // we anticipate that we'll use a more general algorithm in the future, and
   // we don't want to forget to remove the constraint. ***
   if ( FLASH_Obj_depth( A ) != 1 )
   {
      FLA_Print_message( "FLASH_LU_incpiv() currently only supports matrices of depth 1",
                         __FILE__, __LINE__ );
      FLA_Abort();
   }

   if ( FLASH_Queue_stack_depth() == 0 )
      r_val = FLASH_LU_incpiv_opt1( A, p, L );
   else
      r_val = FLASH_LU_incpiv_noopt( A, p, L );

   return r_val;
}

