/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


//
// Level-3 BLAS
//

struct fla_gemm_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_scal_s* sub_scal;
	struct fla_gemm_s* sub_gemm;
};
typedef struct fla_gemm_s fla_gemm_t;


struct fla_hemm_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_scal_s* sub_scal;
	struct fla_hemm_s* sub_hemm;
	struct fla_gemm_s* sub_gemm1;
	struct fla_gemm_s* sub_gemm2;
};
typedef struct fla_hemm_s fla_hemm_t;


struct fla_herk_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_scalr_s* sub_scalr;
	struct fla_herk_s*  sub_herk;
	struct fla_gemm_s*  sub_gemm;
};
typedef struct fla_herk_s fla_herk_t;


struct fla_her2k_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_scalr_s* sub_scalr;
	struct fla_her2k_s* sub_her2k;
	struct fla_gemm_s*  sub_gemm1;
	struct fla_gemm_s*  sub_gemm2;
};
typedef struct fla_her2k_s fla_her2k_t;


struct fla_symm_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_scal_s* sub_scal;
	struct fla_symm_s* sub_symm;
	struct fla_gemm_s* sub_gemm1;
	struct fla_gemm_s* sub_gemm2;
};
typedef struct fla_symm_s fla_symm_t;


struct fla_syrk_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_scalr_s* sub_scalr;
	struct fla_syrk_s*  sub_syrk;
	struct fla_gemm_s*  sub_gemm;
};
typedef struct fla_syrk_s fla_syrk_t;


struct fla_syr2k_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_scalr_s* sub_scalr;
	struct fla_syr2k_s* sub_syr2k;
	struct fla_gemm_s*  sub_gemm1;
	struct fla_gemm_s*  sub_gemm2;
};
typedef struct fla_syr2k_s fla_syr2k_t;


struct fla_trmm_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_scal_s* sub_scal;
	struct fla_trmm_s* sub_trmm;
	struct fla_gemm_s* sub_gemm;
};
typedef struct fla_trmm_s fla_trmm_t;


struct fla_trsm_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_scal_s* sub_scal;
	struct fla_trsm_s* sub_trsm;
	struct fla_gemm_s* sub_gemm;
};
typedef struct fla_trsm_s fla_trsm_t;


#define FLA_Cntl_sub_gemm( cntl )     cntl->sub_gemm
#define FLA_Cntl_sub_gemm1( cntl )    cntl->sub_gemm1
#define FLA_Cntl_sub_gemm2( cntl )    cntl->sub_gemm2
#define FLA_Cntl_sub_gemm3( cntl )    cntl->sub_gemm3
#define FLA_Cntl_sub_gemm4( cntl )    cntl->sub_gemm4
#define FLA_Cntl_sub_gemm5( cntl )    cntl->sub_gemm5
#define FLA_Cntl_sub_gemm6( cntl )    cntl->sub_gemm6
#define FLA_Cntl_sub_gemm7( cntl )    cntl->sub_gemm7
#define FLA_Cntl_sub_gemm8( cntl )    cntl->sub_gemm8
#define FLA_Cntl_sub_hemm( cntl )     cntl->sub_hemm
#define FLA_Cntl_sub_hemm1( cntl )    cntl->sub_hemm1
#define FLA_Cntl_sub_hemm2( cntl )    cntl->sub_hemm2
#define FLA_Cntl_sub_herk( cntl )     cntl->sub_herk
#define FLA_Cntl_sub_herk1( cntl )    cntl->sub_herk1
#define FLA_Cntl_sub_herk2( cntl )    cntl->sub_herk2
#define FLA_Cntl_sub_her2k( cntl )    cntl->sub_her2k
#define FLA_Cntl_sub_symm( cntl )     cntl->sub_symm
#define FLA_Cntl_sub_syrk( cntl )     cntl->sub_syrk
#define FLA_Cntl_sub_syr2k( cntl )    cntl->sub_syr2k
#define FLA_Cntl_sub_trmm( cntl )     cntl->sub_trmm
#define FLA_Cntl_sub_trmm1( cntl )    cntl->sub_trmm1
#define FLA_Cntl_sub_trmm2( cntl )    cntl->sub_trmm2
#define FLA_Cntl_sub_trsm( cntl )     cntl->sub_trsm
#define FLA_Cntl_sub_trsm1( cntl )    cntl->sub_trsm1
#define FLA_Cntl_sub_trsm2( cntl )    cntl->sub_trsm2
#define FLA_Cntl_sub_trsm3( cntl )    cntl->sub_trsm3
#define FLA_Cntl_sub_trsm4( cntl )    cntl->sub_trsm4


fla_gemm_t* FLA_Cntl_gemm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_gemm_t*      sub_gemm );
fla_hemm_t* FLA_Cntl_hemm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_hemm_t*      sub_hemm,
                                      fla_gemm_t*      sub_gemm1,
                                      fla_gemm_t*      sub_gemm2 );
fla_herk_t* FLA_Cntl_herk_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scalr_t*     sub_scalr,
                                      fla_herk_t*      sub_herk,
                                      fla_gemm_t*      sub_gemm );
fla_her2k_t* FLA_Cntl_her2k_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_scalr_t*     sub_scalr,
                                        fla_her2k_t*     sub_her2k,
                                        fla_gemm_t*      sub_gemm1,
                                        fla_gemm_t*      sub_gemm2 );
fla_symm_t* FLA_Cntl_symm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_symm_t*      sub_symm,
                                      fla_gemm_t*      sub_gemm1,
                                      fla_gemm_t*      sub_gemm2 );
fla_syrk_t* FLA_Cntl_syrk_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scalr_t*     sub_scalr,
                                      fla_syrk_t*      sub_syrk,
                                      fla_gemm_t*      sub_gemm );
fla_syr2k_t* FLA_Cntl_syr2k_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_scalr_t*     sub_scalr,
                                        fla_syr2k_t*     sub_syr2k,
                                        fla_gemm_t*      sub_gemm1,
                                        fla_gemm_t*      sub_gemm2 );
fla_trmm_t* FLA_Cntl_trmm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_trmm_t*      sub_trmm,
                                      fla_gemm_t*      sub_gemm );
fla_trsm_t* FLA_Cntl_trsm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_trsm_t*      sub_trsm,
                                      fla_gemm_t*      sub_gemm );

