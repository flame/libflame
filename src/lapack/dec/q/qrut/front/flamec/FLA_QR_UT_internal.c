
#include "FLAME.h"

extern fla_qrut_t* flash_qrut_cntl;
extern fla_qrut_t* flash_qrut_cntl_leaf;
extern fla_qrut_t* fla_qrut_cntl_leaf;

FLA_Error FLA_QR_UT_internal( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl )
{
	FLA_Error r_val = FLA_SUCCESS;
	
	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_QR_UT_internal_check( A, T, cntl );

	if      ( FLA_Cntl_matrix_type( cntl ) == FLA_HIER &&
	          FLA_Cntl_variant( cntl ) == FLA_SUBPROBLEM )
	{
		if ( FLASH_Queue_get_enabled( ) )
		{
			// Enqueue
			ENQUEUE_FLASH_QR_UT_macro( A, *FLASH_OBJ_PTR_AT( T ), cntl );
		}
		else
		{
			// Execute
			r_val = FLA_QR_UT_macro_task( A, *FLASH_OBJ_PTR_AT( T ), cntl );
		}
	}
	else
	{
		if      ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT1 )
		{
			r_val = FLA_QR_UT_unb_var1( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT1 )
		{
			r_val = FLA_QR_UT_opt_var1( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT1 )
		{
			r_val = FLA_QR_UT_blk_var1( A, T, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNBLOCKED_VARIANT2 )
		{
			r_val = FLA_QR_UT_unb_var2( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_UNB_OPT_VARIANT2 )
		{
			r_val = FLA_QR_UT_opt_var2( A, T );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT2 )
		{
			r_val = FLA_QR_UT_blk_var2( A, T, cntl );
		}
		else if ( FLA_Cntl_variant( cntl ) == FLA_BLOCKED_VARIANT3 )
		{
			r_val = FLA_QR_UT_blk_var3( A, T, cntl );
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}

	return r_val;
}

