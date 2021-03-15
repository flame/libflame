/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern TLS_CLASS_SPEC fla_eig_gest_t* flash_eig_gest_cntl;

FLA_Error FLASH_Eig_gest( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B )
{
  FLA_Obj   Y;
  FLA_Error r_val;

  // Check parameters.
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Eig_gest_check( inv, uplo, A, B );

  // The temporary matrix object Y must exist when execution occurs, NOT just
  // when enqueuing occurs. So if the SuperMatrix stack depth is positive, then
  // it means the user has declared a "parallel region" in his code, and thus
  // execution won't occur until sometime after FLASH_Eig_gest() returns, at
  // which time Y will have been deallocated. Thus, we disallow this scenario
  // for now, until we can think of a more general solution.
  if ( FLASH_Queue_stack_depth() > 0 )
  {
     FLA_Print_message( "FLASH_Eig_gest() MUST be invoked with standalone parallelism, and may not be called from within a user-level parallel region",
                        __FILE__, __LINE__ );
     FLA_Abort();
  }


  FLASH_Obj_create_conf_to( FLA_NO_TRANSPOSE, A, &Y );

  // Begin a parallel region.
  FLASH_Queue_begin();
  
  // Enqueue tasks via a SuperMatrix-aware control tree.
  r_val = FLA_Eig_gest_internal( inv, uplo, A, Y, B, flash_eig_gest_cntl );
  
  // End the parallel region.
  FLASH_Queue_end();

  FLASH_Obj_free( &Y );

  return r_val;
}

