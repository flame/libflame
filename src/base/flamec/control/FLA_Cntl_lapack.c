/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

//
// LAPACK-level
//

fla_chol_t* FLA_Cntl_chol_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_chol_t*      sub_chol,
                                      fla_herk_t*      sub_herk,
                                      fla_trsm_t*      sub_trsm,
                                      fla_gemm_t*      sub_gemm )
{
	fla_chol_t* cntl;
	
	cntl = ( fla_chol_t* ) FLA_malloc( sizeof(fla_chol_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_chol    = sub_chol;
	cntl->sub_herk    = sub_herk;
	cntl->sub_trsm    = sub_trsm;
	cntl->sub_gemm    = sub_gemm;

	return cntl;
}

fla_lu_t* FLA_Cntl_lu_obj_create( FLA_Matrix_type  matrix_type,
                                  int              variant,
                                  fla_blocksize_t* blocksize,
                                  fla_lu_t*        sub_lu,
                                  fla_gemm_t*      sub_gemm1,
                                  fla_gemm_t*      sub_gemm2,
                                  fla_gemm_t*      sub_gemm3,
                                  fla_trsm_t*      sub_trsm1,
                                  fla_trsm_t*      sub_trsm2,
                                  fla_appiv_t*     sub_appiv1,
                                  fla_appiv_t*     sub_appiv2 )
{
	fla_lu_t* cntl;
	
	cntl = ( fla_lu_t* ) FLA_malloc( sizeof(fla_lu_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_lu      = sub_lu;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;
	cntl->sub_gemm3   = sub_gemm3;
	cntl->sub_trsm1   = sub_trsm1;
	cntl->sub_trsm2   = sub_trsm2;
	cntl->sub_appiv1  = sub_appiv1;
	cntl->sub_appiv2  = sub_appiv2;

	return cntl;
}


fla_appiv_t* FLA_Cntl_appiv_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_appiv_t*     sub_appiv )
{
	fla_appiv_t* cntl;
	
	cntl = ( fla_appiv_t* ) FLA_malloc( sizeof(fla_appiv_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_appiv   = sub_appiv;

	return cntl;
}


fla_qrut_t* FLA_Cntl_qrut_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_qrut_t*      sub_qrut,
                                      fla_apqut_t*     sub_apqut )
{
	fla_qrut_t* cntl;
	
	cntl = ( fla_qrut_t* ) FLA_malloc( sizeof(fla_qrut_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_qrut    = sub_qrut;
	cntl->sub_apqut   = sub_apqut;

	return cntl;
}

fla_qr2ut_t* FLA_Cntl_qr2ut_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_qr2ut_t*     sub_qr2ut,
                                        fla_gemm_t*      sub_gemm1,
                                        fla_gemm_t*      sub_gemm2,
                                        fla_trsm_t*      sub_trsm,
                                        fla_copy_t*      sub_copy,
                                        fla_axpy_t*      sub_axpy )
{
	fla_qr2ut_t* cntl;
	
	cntl = ( fla_qr2ut_t* ) FLA_malloc( sizeof(fla_qr2ut_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_qr2ut   = sub_qr2ut;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;
	cntl->sub_trsm    = sub_trsm ;
	cntl->sub_copy    = sub_copy ;
	cntl->sub_axpy    = sub_axpy ;

	return cntl;
}

fla_qrutinc_t* FLA_Cntl_qrutinc_obj_create( FLA_Matrix_type  matrix_type,
                                            int              variant,
                                            fla_blocksize_t* blocksize,
                                            fla_qrut_t*      sub_qrut,
                                            fla_apqut_t*     sub_apqut,
                                            fla_qr2ut_t*     sub_qr2ut,
                                            fla_apq2ut_t*    sub_apq2ut )
{
	fla_qrutinc_t* cntl;
	
	cntl = ( fla_qrutinc_t* ) FLA_malloc( sizeof(fla_qrutinc_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_qrut    = sub_qrut;
	cntl->sub_apqut   = sub_apqut;
	cntl->sub_qr2ut   = sub_qr2ut;
	cntl->sub_apq2ut  = sub_apq2ut;

	return cntl;
}

fla_caqrutinc_t* FLA_Cntl_caqrutinc_obj_create( FLA_Matrix_type  matrix_type,
                                                int              variant,
                                                fla_blocksize_t* blocksize,
                                                fla_caqr2ut_t*   sub_caqr2ut,
                                                fla_apcaq2ut_t*  sub_apcaq2ut )
{
	fla_caqrutinc_t* cntl;
	
	cntl = ( fla_caqrutinc_t* ) FLA_malloc( sizeof(fla_caqrutinc_t) );

	cntl->matrix_type  = matrix_type;
	cntl->variant      = variant;
	cntl->blocksize    = blocksize;
	cntl->sub_caqr2ut  = sub_caqr2ut;
	cntl->sub_apcaq2ut = sub_apcaq2ut;

	return cntl;
}

fla_lqut_t* FLA_Cntl_lqut_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_lqut_t*      sub_lqut,
                                      fla_apqut_t*     sub_apqut )
{
	fla_lqut_t* cntl;
	
	cntl = ( fla_lqut_t* ) FLA_malloc( sizeof(fla_lqut_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_lqut    = sub_lqut;
	cntl->sub_apqut   = sub_apqut;

	return cntl;
}

fla_caqr2ut_t* FLA_Cntl_caqr2ut_obj_create( FLA_Matrix_type  matrix_type,
                                            int              variant,
                                            fla_blocksize_t* blocksize,
                                            fla_caqr2ut_t*   sub_caqr2ut,
                                            fla_gemm_t*      sub_gemm1,
                                            fla_gemm_t*      sub_gemm2,
                                            fla_trmm_t*      sub_trmm1,
                                            fla_trmm_t*      sub_trmm2,
                                            fla_trsm_t*      sub_trsm,
                                            fla_axpy_t*      sub_axpy1,
                                            fla_axpy_t*      sub_axpy2,
                                            fla_axpy_t*      sub_axpy3,
                                            fla_copy_t*      sub_copy )
{
	fla_caqr2ut_t* cntl;
	
	cntl = ( fla_caqr2ut_t* ) FLA_malloc( sizeof(fla_caqr2ut_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_caqr2ut = sub_caqr2ut;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;
	cntl->sub_trmm1   = sub_trmm1;
	cntl->sub_trmm2   = sub_trmm2;
	cntl->sub_trsm    = sub_trsm ;
	cntl->sub_axpy1   = sub_axpy1;
	cntl->sub_axpy2   = sub_axpy2;
	cntl->sub_axpy3   = sub_axpy3;
	cntl->sub_copy    = sub_copy ;

	return cntl;
}

fla_hessut_t* FLA_Cntl_hessut_obj_create( FLA_Matrix_type  matrix_type,
                                          int              variant,
                                          fla_blocksize_t* blocksize )
{
	fla_hessut_t* cntl;
	
	cntl = ( fla_hessut_t* ) FLA_malloc( sizeof(fla_hessut_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;

	return cntl;
}

fla_tridiagut_t* FLA_Cntl_tridiagut_obj_create( FLA_Matrix_type  matrix_type,
                                                int              variant,
                                                fla_blocksize_t* blocksize )
{
	fla_tridiagut_t* cntl;
	
	cntl = ( fla_tridiagut_t* ) FLA_malloc( sizeof(fla_tridiagut_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;

	return cntl;
}

fla_bidiagut_t* FLA_Cntl_bidiagut_obj_create( FLA_Matrix_type  matrix_type,
                                              int              variant,
                                              fla_blocksize_t* blocksize )
{
	fla_bidiagut_t* cntl;
	
	cntl = ( fla_bidiagut_t* ) FLA_malloc( sizeof(fla_bidiagut_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;

	return cntl;
}

fla_trinv_t* FLA_Cntl_trinv_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_trinv_t*     sub_trinv,
                                        fla_trmm_t*      sub_trmm,
                                        fla_trsm_t*      sub_trsm1,
                                        fla_trsm_t*      sub_trsm2,
                                        fla_gemm_t*      sub_gemm )
{
	fla_trinv_t* cntl;
	
	cntl = ( fla_trinv_t* ) FLA_malloc( sizeof(fla_trinv_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_trinv   = sub_trinv;
	cntl->sub_trmm    = sub_trmm;
	cntl->sub_trsm1   = sub_trsm1;
	cntl->sub_trsm2   = sub_trsm2;
	cntl->sub_gemm    = sub_gemm;

	return cntl;
}

fla_ttmm_t* FLA_Cntl_ttmm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_ttmm_t*      sub_ttmm,
                                      fla_herk_t*      sub_herk,
                                      fla_trmm_t*      sub_trmm,
                                      fla_gemm_t*      sub_gemm )
{
	fla_ttmm_t* cntl;
	
	cntl = ( fla_ttmm_t* ) FLA_malloc( sizeof(fla_ttmm_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_ttmm    = sub_ttmm;
	cntl->sub_herk    = sub_herk;
	cntl->sub_trmm    = sub_trmm;
	cntl->sub_gemm    = sub_gemm;

	return cntl;
}

fla_uddateut_t* FLA_Cntl_uddateut_obj_create( FLA_Matrix_type  matrix_type,
                                              int              variant,
                                              fla_blocksize_t* blocksize,
                                              fla_uddateut_t*  sub_uddateut,
                                              fla_apqudut_t*   sub_apqudut )
{
	fla_uddateut_t* cntl;
	
	cntl = ( fla_uddateut_t* ) FLA_malloc( sizeof(fla_uddateut_t) );

	cntl->matrix_type  = matrix_type;
	cntl->variant      = variant;
	cntl->blocksize    = blocksize;
	cntl->sub_uddateut = sub_uddateut;
	cntl->sub_apqudut  = sub_apqudut;

	return cntl;
}

fla_uddateutinc_t* FLA_Cntl_uddateutinc_obj_create( FLA_Matrix_type  matrix_type,
                                                    int              variant,
                                                    fla_blocksize_t* blocksize,
                                                    fla_uddateut_t*  sub_uddateut,
                                                    fla_apqudut_t*   sub_apqudut )
{
	fla_uddateutinc_t* cntl;
	
	cntl = ( fla_uddateutinc_t* ) FLA_malloc( sizeof(fla_uddateutinc_t) );

	cntl->matrix_type  = matrix_type;
	cntl->variant      = variant;
	cntl->blocksize    = blocksize;
	cntl->sub_uddateut = sub_uddateut;
	cntl->sub_apqudut  = sub_apqudut;

	return cntl;
}

fla_apqudutinc_t* FLA_Cntl_apqudutinc_obj_create( FLA_Matrix_type  matrix_type,
                                                  int              variant,
                                                  fla_blocksize_t* blocksize,
                                                  fla_apqudut_t*   sub_apqudut )
{
	fla_apqudutinc_t* cntl;
	
	cntl = ( fla_apqudutinc_t* ) FLA_malloc( sizeof(fla_apqudutinc_t) );

	cntl->matrix_type  = matrix_type;
	cntl->variant      = variant;
	cntl->blocksize    = blocksize;
	cntl->sub_apqudut  = sub_apqudut;

	return cntl;
}

fla_sylv_t* FLA_Cntl_sylv_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_sylv_t*      sub_sylv1,
                                      fla_sylv_t*      sub_sylv2,
                                      fla_sylv_t*      sub_sylv3,
                                      fla_gemm_t*      sub_gemm1,
                                      fla_gemm_t*      sub_gemm2,
                                      fla_gemm_t*      sub_gemm3,
                                      fla_gemm_t*      sub_gemm4,
                                      fla_gemm_t*      sub_gemm5,
                                      fla_gemm_t*      sub_gemm6,
                                      fla_gemm_t*      sub_gemm7,
                                      fla_gemm_t*      sub_gemm8 )
{
	fla_sylv_t* cntl;
	
	cntl = ( fla_sylv_t* ) FLA_malloc( sizeof(fla_sylv_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_sylv1   = sub_sylv1;
	cntl->sub_sylv2   = sub_sylv2;
	cntl->sub_sylv3   = sub_sylv3;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;
	cntl->sub_gemm3   = sub_gemm3;
	cntl->sub_gemm4   = sub_gemm4;
	cntl->sub_gemm5   = sub_gemm5;
	cntl->sub_gemm6   = sub_gemm6;
	cntl->sub_gemm7   = sub_gemm7;
	cntl->sub_gemm8   = sub_gemm8;

	return cntl;
}

fla_lyap_t* FLA_Cntl_lyap_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_lyap_t*      sub_lyap,
                                      fla_sylv_t*      sub_sylv,
                                      fla_gemm_t*      sub_gemm1,
                                      fla_gemm_t*      sub_gemm2,
                                      fla_hemm_t*      sub_hemm,
                                      fla_her2k_t*     sub_her2k )
{
	fla_lyap_t* cntl;
	
	cntl = ( fla_lyap_t* ) FLA_malloc( sizeof(fla_lyap_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_scal    = sub_scal;
	cntl->sub_lyap    = sub_lyap;
	cntl->sub_sylv    = sub_sylv;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;
	cntl->sub_hemm    = sub_hemm;
	cntl->sub_her2k   = sub_her2k;

	return cntl;
}

fla_spdinv_t* FLA_Cntl_spdinv_obj_create( FLA_Matrix_type  matrix_type,
                                          int              variant,
                                          fla_blocksize_t* blocksize,
                                          fla_chol_t*      sub_chol,
                                          fla_trinv_t*     sub_trinv,
                                          fla_ttmm_t*      sub_ttmm )
{
	fla_spdinv_t* cntl;
	
	cntl = ( fla_spdinv_t* ) FLA_malloc( sizeof(fla_spdinv_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_chol    = sub_chol;
	cntl->sub_trinv   = sub_trinv;
	cntl->sub_ttmm    = sub_ttmm;

	return cntl;
}

fla_apqut_t* FLA_Cntl_apqut_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_apqut_t*     sub_apqut,
                                        fla_trmm_t*      sub_trmm1,
                                        fla_trmm_t*      sub_trmm2,
                                        fla_gemm_t*      sub_gemm1,
                                        fla_gemm_t*      sub_gemm2,
                                        fla_trsm_t*      sub_trsm,
                                        fla_copyt_t*     sub_copyt,
                                        fla_axpyt_t*     sub_axpyt )
{
	fla_apqut_t* cntl;
	
	cntl = ( fla_apqut_t* ) FLA_malloc( sizeof(fla_apqut_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_apqut   = sub_apqut;
	cntl->sub_trmm1   = sub_trmm1;
	cntl->sub_trmm2   = sub_trmm2;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;
	cntl->sub_trsm    = sub_trsm;
	cntl->sub_copyt   = sub_copyt;
	cntl->sub_axpyt   = sub_axpyt;

	return cntl;
}

fla_apq2ut_t* FLA_Cntl_apq2ut_obj_create( FLA_Matrix_type  matrix_type,
                                          int              variant,
                                          fla_blocksize_t* blocksize,
                                          fla_apq2ut_t*    sub_apq2ut,
                                          fla_gemm_t*      sub_gemm1,
                                          fla_gemm_t*      sub_gemm2,
                                          fla_trsm_t*      sub_trsm,
                                          fla_copyt_t*     sub_copyt,
                                          fla_axpyt_t*     sub_axpyt )
{
	fla_apq2ut_t* cntl;
	
	cntl = ( fla_apq2ut_t* ) FLA_malloc( sizeof(fla_apq2ut_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_apq2ut  = sub_apq2ut;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;
	cntl->sub_trsm    = sub_trsm;
	cntl->sub_copyt   = sub_copyt;
	cntl->sub_axpyt   = sub_axpyt;

	return cntl;
}

fla_apcaq2ut_t* FLA_Cntl_apcaq2ut_obj_create( FLA_Matrix_type  matrix_type,
                                              int              variant,
                                              fla_blocksize_t* blocksize,
                                              fla_apcaq2ut_t*  sub_apcaq2ut,
                                              fla_gemm_t*      sub_gemm1,
                                              fla_gemm_t*      sub_gemm2,
                                              fla_trmm_t*      sub_trmm1,
                                              fla_trmm_t*      sub_trmm2,
                                              fla_trsm_t*      sub_trsm,
                                              fla_axpy_t*      sub_axpy1,
                                              fla_axpy_t*      sub_axpy2,
                                              fla_axpy_t*      sub_axpy3,
                                              fla_copy_t*      sub_copy )
{
	fla_apcaq2ut_t* cntl;
	
	cntl = ( fla_apcaq2ut_t* ) FLA_malloc( sizeof(fla_apcaq2ut_t) );

	cntl->matrix_type  = matrix_type;
	cntl->variant      = variant;
	cntl->blocksize    = blocksize;
	cntl->sub_apcaq2ut = sub_apcaq2ut;
	cntl->sub_gemm1    = sub_gemm1;
	cntl->sub_gemm2    = sub_gemm2;
	cntl->sub_trmm1    = sub_trmm1;
	cntl->sub_trmm2    = sub_trmm2;
	cntl->sub_trsm     = sub_trsm;
	cntl->sub_axpy1    = sub_axpy1;
	cntl->sub_axpy2    = sub_axpy2;
	cntl->sub_axpy3    = sub_axpy3;
	cntl->sub_copy     = sub_copy;

	return cntl;
}


fla_apqutinc_t* FLA_Cntl_apqutinc_obj_create( FLA_Matrix_type  matrix_type,
                                              int              variant,
                                              fla_blocksize_t* blocksize,
                                              fla_apqut_t*     sub_apqut,
                                              fla_apq2ut_t*    sub_apq2ut )
{
	fla_apqutinc_t* cntl;
	
	cntl = ( fla_apqutinc_t* ) FLA_malloc( sizeof(fla_apqutinc_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_apqut   = sub_apqut;
	cntl->sub_apq2ut  = sub_apq2ut;

	return cntl;
}

fla_apcaqutinc_t* FLA_Cntl_apcaqutinc_obj_create( FLA_Matrix_type  matrix_type,
                                                  int              variant,
                                                  fla_blocksize_t* blocksize,
                                                  fla_apcaq2ut_t*  sub_apcaq2ut )
{
	fla_apcaqutinc_t* cntl;
	
	cntl = ( fla_apcaqutinc_t* ) FLA_malloc( sizeof(fla_apcaqutinc_t) );

	cntl->matrix_type  = matrix_type;
	cntl->variant      = variant;
	cntl->blocksize    = blocksize;
	cntl->sub_apcaq2ut = sub_apcaq2ut;

	return cntl;
}

fla_apqudut_t* FLA_Cntl_apqudut_obj_create( FLA_Matrix_type  matrix_type,
                                            int              variant,
                                            fla_blocksize_t* blocksize,
                                            fla_apqudut_t*   sub_apqudut,
                                            fla_gemm_t*      sub_gemm1,
                                            fla_gemm_t*      sub_gemm2,
                                            fla_gemm_t*      sub_gemm3,
                                            fla_gemm_t*      sub_gemm4,
                                            fla_trsm_t*      sub_trsm,
                                            fla_copyt_t*     sub_copyt,
                                            fla_axpyt_t*     sub_axpyt )
{
	fla_apqudut_t* cntl;
	
	cntl = ( fla_apqudut_t* ) FLA_malloc( sizeof(fla_apqudut_t) );

	cntl->matrix_type = matrix_type;
	cntl->variant     = variant;
	cntl->blocksize   = blocksize;
	cntl->sub_apqudut = sub_apqudut;
	cntl->sub_gemm1   = sub_gemm1;
	cntl->sub_gemm2   = sub_gemm2;
	cntl->sub_gemm3   = sub_gemm3;
	cntl->sub_gemm4   = sub_gemm4;
	cntl->sub_trsm    = sub_trsm;
	cntl->sub_copyt   = sub_copyt;
	cntl->sub_axpyt   = sub_axpyt;

	return cntl;
}

fla_eig_gest_t* FLA_Cntl_eig_gest_obj_create( FLA_Matrix_type  matrix_type,
                                              int              variant,
                                              fla_blocksize_t* blocksize,
                                              fla_eig_gest_t*  sub_eig_gest,
                                              fla_axpy_t*      sub_axpy1,
                                              fla_axpy_t*      sub_axpy2,
                                              fla_gemm_t*      sub_gemm1,
                                              fla_gemm_t*      sub_gemm2,
                                              fla_gemm_t*      sub_gemm3,
                                              fla_hemm_t*      sub_hemm,
                                              fla_her2k_t*     sub_her2k,
                                              fla_trmm_t*      sub_trmm1,
                                              fla_trmm_t*      sub_trmm2,
                                              fla_trsm_t*      sub_trsm1,
                                              fla_trsm_t*      sub_trsm2 )
{
	fla_eig_gest_t* cntl;
	
	cntl = ( fla_eig_gest_t* ) FLA_malloc( sizeof(fla_eig_gest_t) );

	cntl->matrix_type  = matrix_type;
	cntl->variant      = variant;
	cntl->blocksize    = blocksize;
	cntl->sub_eig_gest = sub_eig_gest;
	cntl->sub_axpy1    = sub_axpy1;
	cntl->sub_axpy2    = sub_axpy2;
	cntl->sub_gemm1    = sub_gemm1;
	cntl->sub_gemm2    = sub_gemm2;
	cntl->sub_gemm3    = sub_gemm3;
	cntl->sub_hemm     = sub_hemm;
	cntl->sub_her2k    = sub_her2k;
	cntl->sub_trmm1    = sub_trmm1;
	cntl->sub_trmm2    = sub_trmm2;
	cntl->sub_trsm1    = sub_trsm1;
	cntl->sub_trsm2    = sub_trsm2;

	return cntl;
}

