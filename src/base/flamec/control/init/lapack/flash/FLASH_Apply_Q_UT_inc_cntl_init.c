
#include "FLAME.h"

extern fla_apqut_t*  flash_apqut_cntl;
extern fla_apq2ut_t* flash_apq2ut_cntl;

fla_apqutinc_t*      flash_apqutinc_cntl;
fla_blocksize_t*     flash_apqutinc_var1_bsize;

void FLASH_Apply_Q_UT_inc_cntl_init()
{
	// Set blocksizes for hierarchical storage.
	flash_apqutinc_var1_bsize = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree to invoke variant 1.
	flash_apqutinc_cntl = FLA_Cntl_apqutinc_obj_create( FLA_HIER,
	                                                    FLA_BLOCKED_VARIANT1, 
	                                                    flash_apqutinc_var1_bsize,
	                                                    flash_apqut_cntl,
	                                                    flash_apq2ut_cntl );
}

void FLASH_Apply_Q_UT_inc_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_apqutinc_cntl );

	FLA_Blocksize_free( flash_apqutinc_var1_bsize );
}

