/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_QR_UT_copy_internal( FLA_Obj A, FLA_Obj T, FLA_Obj U, fla_qrut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_QR_UT_copy_internal_check( A, T, U, cntl );

	if ( FLASH_Queue_get_enabled() )
	{
		// Enqueue task.
		ENQUEUE_FLASH_QR_UT_copy( *FLASH_OBJ_PTR_AT( A ),
		                          *FLASH_OBJ_PTR_AT( T ),
		                          *FLASH_OBJ_PTR_AT( U ),
		                          NULL );
	}
	else
	{
		// Execute task immediately.
		FLA_QR_UT_copy_task( *FLASH_OBJ_PTR_AT( A ),
		                     *FLASH_OBJ_PTR_AT( T ),
		                     *FLASH_OBJ_PTR_AT( U ),
                             NULL );
	}

	return r_val;
}

