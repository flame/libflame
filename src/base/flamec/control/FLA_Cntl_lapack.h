/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/


//
// LAPACK-level
//

struct fla_chol_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_chol_s* sub_chol;
	struct fla_herk_s* sub_herk;
	struct fla_trsm_s* sub_trsm;
	struct fla_gemm_s* sub_gemm;
};
typedef struct fla_chol_s fla_chol_t;


struct fla_ttmm_s
{
	FLA_Matrix_type    matrix_type;
	int                variant;
	fla_blocksize_t*   blocksize;
	struct fla_ttmm_s* sub_ttmm;
	struct fla_herk_s* sub_herk;
	struct fla_trmm_s* sub_trmm;
	struct fla_gemm_s* sub_gemm;
};
typedef struct fla_ttmm_s fla_ttmm_t;


struct fla_appiv_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_appiv_s* sub_appiv;
};
typedef struct fla_appiv_s fla_appiv_t;


struct fla_lu_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_lu_s*    sub_lu;
	struct fla_gemm_s*  sub_gemm1;
	struct fla_gemm_s*  sub_gemm2;
	struct fla_gemm_s*  sub_gemm3;
	struct fla_trsm_s*  sub_trsm1;
	struct fla_trsm_s*  sub_trsm2;
	struct fla_appiv_s* sub_appiv1;
	struct fla_appiv_s* sub_appiv2;
};
typedef struct fla_lu_s fla_lu_t;


struct fla_qr_ut_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_qr_ut_s* sub_qrut;
	struct fla_apqut_s* sub_apqut;
};
typedef struct fla_qr_ut_s fla_qrut_t;


struct fla_qr2_ut_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_qr2_ut_s*   sub_qr2ut;
	struct fla_gemm_s*     sub_gemm1;
	struct fla_gemm_s*     sub_gemm2;
	struct fla_trsm_s*     sub_trsm;
	struct fla_copy_s*     sub_copy;
	struct fla_axpy_s*     sub_axpy;
};
typedef struct fla_qr2_ut_s fla_qr2ut_t;


struct fla_lq_ut_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_lq_ut_s* sub_lqut;
	struct fla_apqut_s* sub_apqut;
};
typedef struct fla_lq_ut_s fla_lqut_t;

struct fla_caqr2ut_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_caqr2ut_s*  sub_caqr2ut;
	struct fla_gemm_s*     sub_gemm1;
	struct fla_gemm_s*     sub_gemm2;
	struct fla_trmm_s*     sub_trmm1;
	struct fla_trmm_s*     sub_trmm2;
	struct fla_trsm_s*     sub_trsm;
	struct fla_axpy_s*     sub_axpy1;
	struct fla_axpy_s*     sub_axpy2;
	struct fla_axpy_s*     sub_axpy3;
	struct fla_copy_s*     sub_copy;
};
typedef struct fla_caqr2ut_s fla_caqr2ut_t;


struct fla_hess_ut_s
{
	FLA_Matrix_type       matrix_type;
	int                   variant;
	fla_blocksize_t*      blocksize;
};
typedef struct fla_hess_ut_s fla_hessut_t;

struct fla_tridiag_ut_s
{
	FLA_Matrix_type       matrix_type;
	int                   variant;
	fla_blocksize_t*      blocksize;
};
typedef struct fla_tridiag_ut_s fla_tridiagut_t;

struct fla_bidiag_ut_s
{
	FLA_Matrix_type       matrix_type;
	int                   variant;
	fla_blocksize_t*      blocksize;
};
typedef struct fla_bidiag_ut_s fla_bidiagut_t;

struct fla_trinv_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_trinv_s* sub_trinv;
	struct fla_gemm_s*  sub_gemm;
	struct fla_trmm_s*  sub_trmm;
	struct fla_trsm_s*  sub_trsm1;
	struct fla_trsm_s*  sub_trsm2;
};
typedef struct fla_trinv_s fla_trinv_t;


struct fla_sylv_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_sylv_s*  sub_sylv1;
	struct fla_sylv_s*  sub_sylv2;
	struct fla_sylv_s*  sub_sylv3;
	struct fla_gemm_s*  sub_gemm1;
	struct fla_gemm_s*  sub_gemm2;
	struct fla_gemm_s*  sub_gemm3;
	struct fla_gemm_s*  sub_gemm4;
	struct fla_gemm_s*  sub_gemm5;
	struct fla_gemm_s*  sub_gemm6;
	struct fla_gemm_s*  sub_gemm7;
	struct fla_gemm_s*  sub_gemm8;
};
typedef struct fla_sylv_s fla_sylv_t;


struct fla_lyap_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_scal_s*  sub_scal;
	struct fla_lyap_s*  sub_lyap;
	struct fla_sylv_s*  sub_sylv;
	struct fla_gemm_s*  sub_gemm1;
	struct fla_gemm_s*  sub_gemm2;
	struct fla_hemm_s*  sub_hemm;
	struct fla_her2k_s* sub_her2k;
};
typedef struct fla_lyap_s fla_lyap_t;


struct fla_spdinv_s
{
	FLA_Matrix_type     matrix_type;
	int                 variant;
	fla_blocksize_t*    blocksize;
	struct fla_chol_s*  sub_chol;
	struct fla_trinv_s* sub_trinv;
	struct fla_ttmm_s*  sub_ttmm;
};
typedef struct fla_spdinv_s fla_spdinv_t;


struct fla_apqut_s
{
	FLA_Matrix_type      matrix_type;
	int                  variant;
	fla_blocksize_t*     blocksize;
	struct fla_apqut_s*  sub_apqut;
	struct fla_trmm_s*   sub_trmm1;
	struct fla_trmm_s*   sub_trmm2;
	struct fla_gemm_s*   sub_gemm1;
	struct fla_gemm_s*   sub_gemm2;
	struct fla_trsm_s*   sub_trsm;
	struct fla_copyt_s*  sub_copyt;
	struct fla_axpyt_s*  sub_axpyt;
};
typedef struct fla_apqut_s fla_apqut_t;


struct fla_apq2ut_s
{
	FLA_Matrix_type       matrix_type;
	int                   variant;
	fla_blocksize_t*      blocksize;
	struct fla_apq2ut_s*  sub_apq2ut;
	struct fla_gemm_s*    sub_gemm1;
	struct fla_gemm_s*    sub_gemm2;
	struct fla_trsm_s*    sub_trsm;
	struct fla_copyt_s*   sub_copyt;
	struct fla_axpyt_s*   sub_axpyt;
};
typedef struct fla_apq2ut_s fla_apq2ut_t;


struct fla_caqrutinc_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_caqr2ut_s*  sub_caqr2ut;
	struct fla_apcaq2ut_s* sub_apcaq2ut;
};
typedef struct fla_caqrutinc_s fla_caqrutinc_t;


struct fla_apcaqutinc_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_apcaq2ut_s* sub_apcaq2ut;
};
typedef struct fla_apcaqutinc_s fla_apcaqutinc_t;


struct fla_apcaq2ut_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_apcaq2ut_s* sub_apcaq2ut;
	struct fla_gemm_s*     sub_gemm1;
	struct fla_gemm_s*     sub_gemm2;
	struct fla_trmm_s*     sub_trmm1;
	struct fla_trmm_s*     sub_trmm2;
	struct fla_trsm_s*     sub_trsm;
	struct fla_axpy_s*     sub_axpy1;
	struct fla_axpy_s*     sub_axpy2;
	struct fla_axpy_s*     sub_axpy3;
	struct fla_copy_s*     sub_copy;
};
typedef struct fla_apcaq2ut_s fla_apcaq2ut_t;


struct fla_qr_ut_inc_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_qr_ut_s*    sub_qrut;
	struct fla_qr2_ut_s*   sub_qr2ut;
	struct fla_apqut_s*    sub_apqut;
	struct fla_apq2ut_s*   sub_apq2ut;
};
typedef struct fla_qr_ut_inc_s fla_qrutinc_t;


struct fla_apqutinc_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_apqut_s*    sub_apqut;
	struct fla_apq2ut_s*   sub_apq2ut;
};
typedef struct fla_apqutinc_s fla_apqutinc_t;


struct fla_uddateut_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_uddateut_s* sub_uddateut;
	struct fla_apqudut_s*  sub_apqudut;
};
typedef struct fla_uddateut_s fla_uddateut_t;


struct fla_apqudut_s
{
	FLA_Matrix_type       matrix_type;
	int                   variant;
	fla_blocksize_t*      blocksize;
	struct fla_apqudut_s* sub_apqudut;
	struct fla_gemm_s*    sub_gemm1;
	struct fla_gemm_s*    sub_gemm2;
	struct fla_gemm_s*    sub_gemm3;
	struct fla_gemm_s*    sub_gemm4;
	struct fla_trsm_s*    sub_trsm;
	struct fla_copyt_s*   sub_copyt;
	struct fla_axpyt_s*   sub_axpyt;
};
typedef struct fla_apqudut_s fla_apqudut_t;


struct fla_uddateutinc_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_uddateut_s* sub_uddateut;
	struct fla_apqudut_s*  sub_apqudut;
};
typedef struct fla_uddateutinc_s fla_uddateutinc_t;


struct fla_apqudutinc_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_apqudut_s*  sub_apqudut;
};
typedef struct fla_apqudutinc_s fla_apqudutinc_t;


struct fla_eig_gest_s
{
	FLA_Matrix_type        matrix_type;
	int                    variant;
	fla_blocksize_t*       blocksize;
	struct fla_eig_gest_s* sub_eig_gest;
	struct fla_axpy_s*     sub_axpy1;
	struct fla_axpy_s*     sub_axpy2;
	struct fla_gemm_s*     sub_gemm1;
	struct fla_gemm_s*     sub_gemm2;
	struct fla_gemm_s*     sub_gemm3;
	struct fla_hemm_s*     sub_hemm;
	struct fla_her2k_s*    sub_her2k;
	struct fla_trmm_s*     sub_trmm1;
	struct fla_trmm_s*     sub_trmm2;
	struct fla_trsm_s*     sub_trsm1;
	struct fla_trsm_s*     sub_trsm2;
};
typedef struct fla_eig_gest_s fla_eig_gest_t;


#define FLA_Cntl_sub_chol( cntl )      cntl->sub_chol
#define FLA_Cntl_sub_lu( cntl )        cntl->sub_lu
#define FLA_Cntl_sub_qr( cntl )        cntl->sub_qr
#define FLA_Cntl_sub_qrut( cntl )      cntl->sub_qrut
#define FLA_Cntl_sub_qr2ut( cntl )     cntl->sub_qr2ut
#define FLA_Cntl_sub_lq( cntl )        cntl->sub_lq
#define FLA_Cntl_sub_lqut( cntl )      cntl->sub_lqut
#define FLA_Cntl_sub_caqr2ut( cntl )   cntl->sub_caqr2ut
#define FLA_Cntl_sub_trinv( cntl )     cntl->sub_trinv
#define FLA_Cntl_sub_ttmm( cntl )      cntl->sub_ttmm
#define FLA_Cntl_sub_sylv( cntl )      cntl->sub_sylv
#define FLA_Cntl_sub_sylv1( cntl )     cntl->sub_sylv1
#define FLA_Cntl_sub_sylv2( cntl )     cntl->sub_sylv2
#define FLA_Cntl_sub_sylv3( cntl )     cntl->sub_sylv3
#define FLA_Cntl_sub_lyap( cntl )      cntl->sub_lyap
#define FLA_Cntl_sub_appiv( cntl )     cntl->sub_appiv
#define FLA_Cntl_sub_appiv1( cntl )    cntl->sub_appiv1
#define FLA_Cntl_sub_appiv2( cntl )    cntl->sub_appiv2
#define FLA_Cntl_sub_apqut( cntl )     cntl->sub_apqut
#define FLA_Cntl_sub_apq2ut( cntl )    cntl->sub_apq2ut
#define FLA_Cntl_sub_apcaq2ut( cntl )  cntl->sub_apcaq2ut
#define FLA_Cntl_sub_uddateut( cntl )  cntl->sub_uddateut
#define FLA_Cntl_sub_apqudut( cntl )   cntl->sub_apqudut
#define FLA_Cntl_sub_hessut( cntl )    cntl->sub_hessut
#define FLA_Cntl_sub_tridiagut( cntl ) cntl->sub_tridiagut
#define FLA_Cntl_sub_bidiagut( cntl )  cntl->sub_bidiagut
#define FLA_Cntl_sub_eig_gest( cntl )  cntl->sub_eig_gest


fla_chol_t* FLA_Cntl_chol_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_chol_t*      sub_chol,
                                      fla_herk_t*      sub_herk,
                                      fla_trsm_t*      sub_trsm,
                                      fla_gemm_t*      sub_gemm );
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
                                  fla_appiv_t*     sub_appiv2 );
fla_appiv_t* FLA_Cntl_appiv_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_appiv_t*     sub_appiv );
fla_qrut_t* FLA_Cntl_qrut_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_qrut_t*      sub_qrut,
                                      fla_apqut_t*     sub_apqut );
fla_qr2ut_t* FLA_Cntl_qr2ut_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_qr2ut_t*     sub_qr2ut,
                                        fla_gemm_t*      sub_gemm1,
                                        fla_gemm_t*      sub_gemm2,
                                        fla_trsm_t*      sub_trsm,
                                        fla_copy_t*      sub_copy,
                                        fla_axpy_t*      sub_axpy );
fla_lqut_t* FLA_Cntl_lqut_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_lqut_t*      sub_lqut,
                                      fla_apqut_t*     sub_apqut );
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
                                            fla_copy_t*      sub_copy );
fla_hessut_t* FLA_Cntl_hessut_obj_create( FLA_Matrix_type  matrix_type,
                                          int              variant,
                                          fla_blocksize_t* blocksize );
fla_tridiagut_t* FLA_Cntl_tridiagut_obj_create( FLA_Matrix_type  matrix_type,
                                                int              variant,
                                                fla_blocksize_t* blocksize );
fla_bidiagut_t* FLA_Cntl_bidiagut_obj_create( FLA_Matrix_type  matrix_type,
                                              int              variant,
                                              fla_blocksize_t* blocksize );
fla_trinv_t* FLA_Cntl_trinv_obj_create( FLA_Matrix_type  matrix_type,
                                        int              variant,
                                        fla_blocksize_t* blocksize,
                                        fla_trinv_t*     sub_trinv,
                                        fla_trmm_t*      sub_trmm,
                                        fla_trsm_t*      sub_trsm1,
                                        fla_trsm_t*      sub_trsm2,
                                        fla_gemm_t*      sub_gemm );
fla_ttmm_t* FLA_Cntl_ttmm_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_ttmm_t*      sub_ttmm,
                                      fla_herk_t*      sub_herk,
                                      fla_trmm_t*      sub_trmm,
                                      fla_gemm_t*      sub_gemm );
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
                                      fla_gemm_t*      sub_gemm8 );
fla_lyap_t* FLA_Cntl_lyap_obj_create( FLA_Matrix_type  matrix_type,
                                      int              variant,
                                      fla_blocksize_t* blocksize,
                                      fla_scal_t*      sub_scal,
                                      fla_lyap_t*      sub_lyap,
                                      fla_sylv_t*      sub_sylv,
                                      fla_gemm_t*      sub_gemm1,
                                      fla_gemm_t*      sub_gemm2,
                                      fla_hemm_t*      sub_hemm,
                                      fla_her2k_t*     sub_her2k );
fla_spdinv_t* FLA_Cntl_spdinv_obj_create( FLA_Matrix_type  matrix_type,
                                          int              variant,
                                          fla_blocksize_t* blocksize,
                                          fla_chol_t*      sub_chol,
                                          fla_trinv_t*     sub_trinv,
                                          fla_ttmm_t*      sub_ttmm );
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
                                        fla_axpyt_t*     sub_axpyt );
fla_apq2ut_t* FLA_Cntl_apq2ut_obj_create( FLA_Matrix_type  matrix_type,
                                          int              variant,
                                          fla_blocksize_t* blocksize,
                                          fla_apq2ut_t*    sub_apq2ut,
                                          fla_gemm_t*      sub_gemm1,
                                          fla_gemm_t*      sub_gemm2,
                                          fla_trsm_t*      sub_trsm,
                                          fla_copyt_t*     sub_copyt,
                                          fla_axpyt_t*     sub_axpyt );
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
                                              fla_copy_t*      sub_copy );
fla_qrutinc_t* FLA_Cntl_qrutinc_obj_create( FLA_Matrix_type  matrix_type,
                                            int              variant,
                                            fla_blocksize_t* blocksize,
                                            fla_qrut_t*      sub_qrut,
                                            fla_apqut_t*     sub_apqut,
                                            fla_qr2ut_t*     sub_qr2ut,
                                            fla_apq2ut_t*    sub_apq2ut );
fla_apqutinc_t* FLA_Cntl_apqutinc_obj_create( FLA_Matrix_type  matrix_type,
                                              int              variant,
                                              fla_blocksize_t* blocksize,
                                              fla_apqut_t*     sub_apqut,
                                              fla_apq2ut_t*    sub_apq2ut );
fla_caqrutinc_t* FLA_Cntl_caqrutinc_obj_create( FLA_Matrix_type  matrix_type,
                                                int              variant,
                                                fla_blocksize_t* blocksize,
                                                fla_caqr2ut_t*   sub_caqr2ut,
                                                fla_apcaq2ut_t*  sub_apcaq2ut );
fla_apcaqutinc_t* FLA_Cntl_apcaqutinc_obj_create( FLA_Matrix_type  matrix_type,
                                                  int              variant,
                                                  fla_blocksize_t* blocksize,
                                                  fla_apcaq2ut_t*  sub_apcaq2ut );
fla_uddateut_t* FLA_Cntl_uddateut_obj_create( FLA_Matrix_type  matrix_type,
                                              int              variant,
                                              fla_blocksize_t* blocksize,
                                              fla_uddateut_t*  sub_uddateut,
                                              fla_apqudut_t*   sub_apqudut );
fla_apqudut_t* FLA_Cntl_apqudut_obj_create( FLA_Matrix_type  matrix_type,
                                            int              variant,
                                            fla_blocksize_t* blocksize,
                                            fla_apqudut_t*   sub_apq2ut,
                                            fla_gemm_t*      sub_gemm1,
                                            fla_gemm_t*      sub_gemm2,
                                            fla_gemm_t*      sub_gemm3,
                                            fla_gemm_t*      sub_gemm4,
                                            fla_trsm_t*      sub_trsm,
                                            fla_copyt_t*     sub_copyt,
                                            fla_axpyt_t*     sub_axpyt );
fla_uddateutinc_t* FLA_Cntl_uddateutinc_obj_create( FLA_Matrix_type  matrix_type,
                                                    int              variant,
                                                    fla_blocksize_t* blocksize,
                                                    fla_uddateut_t*  sub_uddateut,
                                                    fla_apqudut_t*   sub_apqudut );
fla_apqudutinc_t* FLA_Cntl_apqudutinc_obj_create( FLA_Matrix_type  matrix_type,
                                                  int              variant,
                                                  fla_blocksize_t* blocksize,
                                                  fla_apqudut_t*   sub_apqudut );
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
                                              fla_trsm_t*      sub_trsm2 );

