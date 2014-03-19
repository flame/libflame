/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

//
// Level-1 BLAS
//

fla_axpy_t* FLA_Cntl_axpy_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_axpy_t*      sub_axpy )
{
	fla_axpy_t* cntl;
	
	cntl = ( fla_axpy_t* ) FLA_malloc( sizeof(fla_axpy_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_axpy    = sub_axpy;

	return cntl;
}

fla_axpyt_t* FLA_Cntl_axpyt_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_axpyt_t*     sub_axpyt )
{
	fla_axpyt_t* cntl;
	
	cntl = ( fla_axpyt_t* ) FLA_malloc( sizeof(fla_axpyt_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_axpyt   = sub_axpyt;

	return cntl;
}

fla_copy_t* FLA_Cntl_copy_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_copy_t*      sub_copy )
{
	fla_copy_t* cntl;
	
	cntl = ( fla_copy_t* ) FLA_malloc( sizeof(fla_copy_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_copy    = sub_copy;

	return cntl;
}

fla_copyt_t* FLA_Cntl_copyt_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_copyt_t*     sub_copyt )
{
	fla_copyt_t* cntl;
	
	cntl = ( fla_copyt_t* ) FLA_malloc( sizeof(fla_copyt_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_copyt   = sub_copyt;

	return cntl;
}

fla_copyr_t* FLA_Cntl_copyr_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_copyr_t*     sub_copyr,
                                        fla_copy_t*      sub_copy )
{
	fla_copyr_t* cntl;
	
	cntl = ( fla_copyr_t* ) FLA_malloc( sizeof(fla_copyr_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_copyr   = sub_copyr;
	cntl->sub_copy    = sub_copy;

	return cntl;
}

fla_scal_t* FLA_Cntl_scal_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal )
{
	fla_scal_t* cntl;
	
	cntl = ( fla_scal_t* ) FLA_malloc( sizeof(fla_scal_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scal    = sub_scal;

	return cntl;
}

fla_scalr_t* FLA_Cntl_scalr_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_scalr_t*     sub_scalr,
                                        fla_scal_t*      sub_scal )
{
	fla_scalr_t* cntl;
	
	cntl = ( fla_scalr_t* ) FLA_malloc( sizeof(fla_scalr_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scalr   = sub_scalr;
	cntl->sub_scal    = sub_scal;

	return cntl;
}

fla_swap_t* FLA_Cntl_swap_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_swap_t*      sub_swap )
{
	fla_swap_t* cntl;
	
	cntl = ( fla_swap_t* ) FLA_malloc( sizeof(fla_swap_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_swap    = sub_swap;

	return cntl;
}

fla_tpose_t* FLA_Cntl_tpose_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_tpose_t*     sub_trans,
                                        fla_swap_t*      sub_swap )
{
	fla_tpose_t* cntl;
	
	cntl = ( fla_tpose_t* ) FLA_malloc( sizeof(fla_tpose_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_trans   = sub_trans;
	cntl->sub_swap    = sub_swap;

	return cntl;
}

