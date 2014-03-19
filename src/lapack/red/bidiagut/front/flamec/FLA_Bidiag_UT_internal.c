/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_internal( FLA_Obj A, FLA_Obj TU, FLA_Obj TV, fla_bidiagut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Bidiag_UT_internal_check( A, TU, TV, cntl );

	if ( FLA_Obj_length( A ) >= FLA_Obj_width( A ) )
	{
          r_val = FLA_Bidiag_UT_u( A, TU, TV, cntl );
	}
	else // if ( FLA_Obj_length( A ) < FLA_Obj_width( A ) )
	{
          // Flip A; swap(rs, cs), swap(m, n)
          FLA_Obj_flip_base( &A );
          FLA_Obj_flip_view( &A );
          
          r_val = FLA_Bidiag_UT_u( A, TV, TU, cntl );

          // Recover A; swap(rs, cs), swap(m, n)
          FLA_Obj_flip_base( &A );
          FLA_Obj_flip_view( &A );

          // According to the UT transform, the house-holder vectors are conjugated 
          // when they are applied from the right. 
          if ( FLA_Obj_is_complex( A ) ) 
          {
            FLA_Obj ATL, ATR,
                    ABL, ABR;
            dim_t   b;

            FLA_Conjugate( TU );
            FLA_Conjugate( TV );

            // U
            b = ( FLA_Obj_length( A ) - 1 );
            FLA_Part_2x2( A,    &ATL, &ATR,
                                &ABL, &ABR,    2, b, FLA_TL );
            FLA_Conjugate_r( FLA_LOWER_TRIANGULAR, ABL );
            
            // V
            b = ( FLA_Obj_width( A ) - 1 );
            FLA_Part_1x2( A,    &ATL, &ATR,    b, FLA_RIGHT );
            FLA_Conjugate_r( FLA_UPPER_TRIANGULAR, ATR );
          }
        }

	return r_val;
}

