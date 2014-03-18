
#include "FLAME.h"

//
// Level-2 BLAS
//

fla_gemv_t* FLA_Cntl_gemv_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize, 
                                      fla_scal_t*      sub_scal,
                                      fla_gemv_t*      sub_gemv )
{
	fla_gemv_t* cntl;
	
	cntl = ( fla_gemv_t* ) FLA_malloc( sizeof(fla_gemv_t) );
	
	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scal    = sub_scal;
	cntl->sub_gemv    = sub_gemv;

	return cntl;
}

fla_trsv_t* FLA_Cntl_trsv_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize, 
                                      fla_trsv_t*      sub_trsv,
                                      fla_gemv_t*      sub_gemv )
{
	fla_trsv_t* cntl;
	
	cntl = ( fla_trsv_t* ) FLA_malloc( sizeof(fla_trsv_t) );
	
	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_trsv    = sub_trsv;
	cntl->sub_gemv    = sub_gemv;

	return cntl;
}

