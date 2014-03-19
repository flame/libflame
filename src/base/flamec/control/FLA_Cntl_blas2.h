/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


//
// Level-2 BLAS
//

struct fla_gemv_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_scal_s* sub_scal;
	struct fla_gemv_s* sub_gemv;
};
typedef struct fla_gemv_s fla_gemv_t;

struct fla_trsv_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_trsv_s* sub_trsv;
	struct fla_gemv_s* sub_gemv;
};
typedef struct fla_trsv_s fla_trsv_t;


#define FLA_Cntl_sub_gemv( cntl )     cntl->sub_gemv
#define FLA_Cntl_sub_trsv( cntl )     cntl->sub_trsv


fla_gemv_t* FLA_Cntl_gemv_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_gemv_t*      sub_gemv );
fla_trsv_t* FLA_Cntl_trsv_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_trsv_t*      sub_trsv,
                                      fla_gemv_t*      sub_gemv );

