/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
void FLASH_Copyr_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_cntl_flash_init_i)
{
	// Set blocksize for hierarchical storage.
	FLA_cntl_flash_init_i->flash_copyr_bsize     = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are small.
	FLA_cntl_flash_init_i->flash_copyr_cntl_blas = FLA_Cntl_copyr_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that marches through A and B from TL to BR.
	FLA_cntl_flash_init_i->flash_copyr_cntl      = FLA_Cntl_copyr_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   FLA_cntl_flash_init_i->flash_copyr_bsize,
	                                                   FLA_cntl_flash_init_i->flash_copyr_cntl_blas,
	                                                   FLA_cntl_flash_init_i->flash_copy_cntl_tb );
}

void FLASH_Copyr_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_cntl_flash_init_i)
{
	FLA_Cntl_obj_free( FLA_cntl_flash_init_i->flash_copyr_cntl_blas );

	FLA_Cntl_obj_free( FLA_cntl_flash_init_i->flash_copyr_cntl );

	FLA_Blocksize_free( FLA_cntl_flash_init_i->flash_copyr_bsize );
}

#endif

extern fla_copy_t* flash_copy_cntl_tb;

fla_copyr_t*       flash_copyr_cntl_blas = NULL;
fla_copyr_t*       flash_copyr_cntl;
fla_blocksize_t*   flash_copyr_bsize;

void FLASH_Copyr_cntl_init()
{
	// Set blocksize for hierarchical storage.
	flash_copyr_bsize     = FLA_Blocksize_create( 1, 1, 1, 1 );

	// Create a control tree that assumes A and B are small.
	flash_copyr_cntl_blas = FLA_Cntl_copyr_obj_create( FLA_HIER,
	                                                   FLA_SUBPROBLEM,
	                                                   NULL,
	                                                   NULL,
	                                                   NULL );

	// Create a control tree that marches through A and B from TL to BR.
	flash_copyr_cntl      = FLA_Cntl_copyr_obj_create( FLA_HIER,
	                                                   FLA_BLOCKED_VARIANT3,
	                                                   flash_copyr_bsize,
	                                                   flash_copyr_cntl_blas,
	                                                   flash_copy_cntl_tb );
}

void FLASH_Copyr_cntl_finalize()
{
	FLA_Cntl_obj_free( flash_copyr_cntl_blas );

	FLA_Cntl_obj_free( flash_copyr_cntl );

	FLA_Blocksize_free( flash_copyr_bsize );
}

