
#include "FLAME.h"

FLA_Error FLA_UDdate_UT_solve( FLA_Obj R, FLA_Obj bR, FLA_Obj x )
{
	// Check parameters.
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLA_UDdate_UT_solve_check( R, bR, x );

	// Copy the contents of bR to x so that after the triangular solve, the
	// solution resides in x (and bR is preserved).
	FLA_Copy_external( bR, x );
	
	// Perform a triangular solve with R the right-hand side.
	FLA_Trsm_external( FLA_LEFT, FLA_UPPER_TRIANGULAR,
	                   FLA_NO_TRANSPOSE, FLA_NONUNIT_DIAG,
	                   FLA_ONE, R, x );

	return FLA_SUCCESS;
}

