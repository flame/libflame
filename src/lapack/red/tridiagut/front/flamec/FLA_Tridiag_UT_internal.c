/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_internal( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T, fla_tridiagut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Tridiag_UT_internal_check( uplo, A, T, cntl );

	if ( uplo == FLA_LOWER_TRIANGULAR )
	{
          r_val = FLA_Tridiag_UT_l( A, T, cntl );
	}
	else // if ( uplo == FLA_UPPER_TRIANGULAR )
	{
          FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
  	}

	return r_val;
}

