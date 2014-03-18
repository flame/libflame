
#include "FLAME.h"

FLA_Error FLA_QR_UT_piv_internal( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p, fla_qrut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
          FLA_QR_UT_piv_internal_check( A, T, w, p, cntl );

        // Blocked blas2 oriented version
        // ------------ Variant 1 set -------------
        if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
          {
            r_val = FLA_QR_UT_piv_unb_var1( A, T, w, p );
          }
        else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
          {
            //r_val = FLA_QR_UT_piv_opt_var1( A, T );
            FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
          }
        else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
          {
            r_val = FLA_QR_UT_piv_blk_var1( A, T, w, p, cntl );
          }
        // ----------------------------------------
        // ------------ Variant 2 set -------------
        else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT2 )
          {
            r_val = FLA_QR_UT_piv_unb_var2( A, T, w, p );
          }

        else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT2 )
          {
            //r_val = FLA_QR_UT_opt_var2( A, T );
            FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
          }
        else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
          {
            r_val = FLA_QR_UT_piv_blk_var2( A, T, w, p, cntl );
          }
        // ----------------------------------------
        else
          {
            FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
          }

	return r_val;
}

