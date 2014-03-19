/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_UDdate_UT_inc_solve( FLA_Obj R, FLA_Obj bR, FLA_Obj x )
{
	// Check parameters.
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_UDdate_UT_inc_solve_check( R, bR, x );

	// Copy the contents of bR to x so that after the triangular solve, the
	// solution resides in x (and bR is preserved).
	FLASH_Copy( bR, x );
	
	// Perform a triangular solve with R the right-hand side.
	FLASH_Trsm( FLA_LEFT, FLA_UPPER_TRIANGULAR,
	            FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
	            FLA_ONE, R, x );

	return FLA_SUCCESS;
}

