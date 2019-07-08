/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES

typedef struct FLA_Cntl_init_flash_s 
{
  /* FLASH_Axpy_cntl_init */
  fla_axpy_t*        flash_axpy_cntl_blas;
  fla_axpy_t*        flash_axpy_cntl_tb;
  fla_axpy_t*        flash_axpy_cntl;
  fla_blocksize_t*   flash_axpy_bsize;
  
  /* FLASH_Axpyt_cntl_init */
  fla_axpyt_t*       flash_axpyt_cntl_blas;
  fla_axpyt_t*       flash_axpyt_cntl_tb;
  fla_axpyt_t*       flash_axpyt_cntl_lr;
  fla_axpyt_t*       flash_axpyt_cntl;
  fla_blocksize_t*   flash_axpyt_bsize;
  
  /* FLASH_Copy_cntl_init */
  fla_copy_t*        flash_copy_cntl_blas;
  fla_copy_t*        flash_copy_cntl_tb;
  fla_copy_t*        flash_copy_cntl;
  fla_blocksize_t*   flash_copy_bsize;
  
  /* FLASH_Copyt_cntl_init */
  fla_copyt_t*       flash_copyt_cntl_blas;
  fla_copyt_t*       flash_copyt_cntl_tb;
  fla_copyt_t*       flash_copyt_cntl_lr;
  fla_copyt_t*       flash_copyt_cntl;
  fla_blocksize_t*   flash_copyt_bsize;
  
  /* FLASH_Copyr_cntl_init */
  fla_copyr_t*       flash_copyr_cntl_blas;
  fla_copyr_t*       flash_copyr_cntl;
  fla_blocksize_t*   flash_copyr_bsize;
  
  /* FLASH_Scal_cntl_init */
  fla_scal_t*        flash_scal_cntl_blas;
  fla_scal_t*        flash_scal_cntl_tb;
  fla_scal_t*        flash_scal_cntl_lr;
  fla_scal_t*        flash_scal_cntl;
  fla_blocksize_t*   flash_scal_bsize;
  
  /* FLASH_Scalr_cntl_init */
  fla_scalr_t*       flash_scalr_cntl_blas;
  fla_scalr_t*       flash_scalr_cntl;
  fla_blocksize_t*   flash_scalr_bsize;
  
  /* FLASH_Gemv_cntl_init */
  fla_gemv_t*      flash_gemv_cntl_blas;
  fla_gemv_t*      flash_gemv_cntl_fm_rp;
  fla_gemv_t*      flash_gemv_cntl_fm_cp;
  fla_gemv_t*      flash_gemv_cntl_rp_bv;
  fla_gemv_t*      flash_gemv_cntl_cp_bv;
  fla_blocksize_t* flash_gemv_bsize;
  
  /* FLASH_Trsv_cntl_init */
  fla_trsv_t*        flash_trsv_cntl_blas;
  fla_trsv_t*        flash_trsv_cntl;
  fla_blocksize_t*   flash_trsv_bsize;
  
  /* FLASH_Gemm_cntl_init */
  fla_gemm_t*      flash_gemm_cntl_blas;
  fla_gemm_t*      flash_gemm_cntl_mm_mp;
  fla_gemm_t*      flash_gemm_cntl_mm_pm;
  fla_gemm_t*      flash_gemm_cntl_mm_op;
  fla_gemm_t*      flash_gemm_cntl_mp_pb;
  fla_gemm_t*      flash_gemm_cntl_mp_ip;
  fla_gemm_t*      flash_gemm_cntl_pm_bp;
  fla_gemm_t*      flash_gemm_cntl_pm_ip;
  fla_gemm_t*      flash_gemm_cntl_op_bp;
  fla_gemm_t*      flash_gemm_cntl_op_pb;
  fla_gemm_t*      flash_gemm_cntl_pb_bb;
  fla_gemm_t*      flash_gemm_cntl_bp_bb;
  fla_gemm_t*      flash_gemm_cntl_ip_bb;
  
  fla_gemm_t*      flash_gemm_cntl_mm;
  fla_gemm_t*      flash_gemm_cntl_mp;
  fla_gemm_t*      flash_gemm_cntl_pm;
  fla_gemm_t*      flash_gemm_cntl_op;
  fla_gemm_t*      flash_gemm_cntl_pb;
  fla_gemm_t*      flash_gemm_cntl_bp;
  fla_gemm_t*      flash_gemm_cntl_ip;
  
  fla_blocksize_t* flash_gemm_bsize;
  
  /* FLASH_Hemm_cntl_init */
  fla_hemm_t*        flash_hemm_cntl_blas;
  fla_hemm_t*        flash_hemm_cntl_bp;
  fla_hemm_t*        flash_hemm_cntl_mp;
  fla_hemm_t*        flash_hemm_cntl_mm;
  fla_blocksize_t*   flash_hemm_bsize;
  
  /* FLASH_Herk_cntl_init */
  fla_herk_t*         flash_herk_cntl_blas;
  fla_herk_t*         flash_herk_cntl_ip;
  fla_herk_t*         flash_herk_cntl_op;
  fla_herk_t*         flash_herk_cntl_mm;
  fla_blocksize_t*    flash_herk_bsize;
  
  /* FLASH_Her2k_cntl_init */
  fla_her2k_t*        flash_her2k_cntl_blas;
  fla_her2k_t*        flash_her2k_cntl_ip;
  fla_her2k_t*        flash_her2k_cntl_op;
  fla_her2k_t*        flash_her2k_cntl_mm;
  fla_blocksize_t*    flash_her2k_bsize;
  
  /* FLASH_Symm_cntl_init */
  fla_symm_t*        flash_symm_cntl_blas;
  fla_symm_t*        flash_symm_cntl_bp;
  fla_symm_t*        flash_symm_cntl_mp;
  fla_symm_t*        flash_symm_cntl_mm;
  fla_blocksize_t*   flash_symm_bsize;
  
  /* FLASH_Syrk_cntl_init */
  fla_syrk_t*         flash_syrk_cntl_blas;
  fla_syrk_t*         flash_syrk_cntl_ip;
  fla_syrk_t*         flash_syrk_cntl_op;
  fla_syrk_t*         flash_syrk_cntl_mm;
  fla_blocksize_t*    flash_syrk_bsize;
  
  /* FLASH_Syr2k_cntl_init */
  fla_syr2k_t*        flash_syr2k_cntl_blas;
  fla_syr2k_t*        flash_syr2k_cntl_ip;
  fla_syr2k_t*        flash_syr2k_cntl_op;
  fla_syr2k_t*        flash_syr2k_cntl_mm;
  fla_blocksize_t*    flash_syr2k_bsize;
  
  /* FLASH_Trmm_cntl_init */
  fla_trmm_t*        flash_trmm_cntl_blas;
  fla_trmm_t*        flash_trmm_cntl_bp;
  fla_trmm_t*        flash_trmm_cntl_mp;
  fla_trmm_t*        flash_trmm_cntl_mm;
  fla_blocksize_t*   flash_trmm_bsize;
  
  /* FLASH_Trsm_cntl_init */
  fla_trsm_t*        flash_trsm_cntl_blas;
  fla_trsm_t*        flash_trsm_cntl_bp;
  fla_trsm_t*        flash_trsm_cntl_mp;
  fla_trsm_t*        flash_trsm_cntl_mm;
  fla_blocksize_t*   flash_trsm_bsize;
  
  /* FLASH_Apply_pivots_cntl_init */
  fla_appiv_t*       flash_appiv_cntl_leaf;
  fla_appiv_t*       flash_appiv_cntl_bp;
  fla_appiv_t*       flash_appiv_cntl;
  fla_blocksize_t*   flash_appiv_bsize;
  
  /* FLASH_Chol_cntl_init */
  fla_chol_t*        flash_chol_cntl_leaf;
  fla_chol_t*        flash_chol_cntl;
  fla_blocksize_t*   flash_chol_bsize;
  
  /* FLASH_LU_nopiv_cntl_init */
  fla_lu_t*          flash_lu_nopiv_cntl_leaf;
  fla_lu_t*          flash_lu_nopiv_cntl;
  fla_blocksize_t*   flash_lu_nopiv_bsize;
  
  /* FLASH_LU_piv_cntl_init */
  fla_lu_t*           flash_lu_piv_cntl_leaf;
  fla_lu_t*           flash_lu_piv_cntl;
  fla_blocksize_t*    flash_lu_piv_bsize;
  
  /* FLASH_LU_incpiv_cntl_init */
  fla_lu_t*           flash_lu_incpiv_cntl_leaf;
  fla_lu_t*           flash_lu_incpiv_cntl;
  fla_blocksize_t*    flash_lu_incpiv_bsize;
  
  /* FLASH_Trinv_cntl_init */
  fla_trinv_t*       flash_trinv_cntl_leaf;
  fla_trinv_t*       flash_trinv_cntl;
  fla_blocksize_t*   flash_trinv_bsize;
  
  /* FLASH_Ttmm_cntl_init */
  fla_ttmm_t*        flash_ttmm_cntl_leaf;
  fla_ttmm_t*        flash_ttmm_cntl;
  fla_blocksize_t*   flash_ttmm_bsize;
  
  /* FLASH_Sylv_cntl_init */
  fla_sylv_t*        flash_sylv_cntl_leaf;
  fla_sylv_t*        flash_sylv_cntl_mb;
  fla_sylv_t*        flash_sylv_cntl;
  fla_blocksize_t*   flash_sylv_bsize;
  
  /* FLASH_QR2_UT_cntl_init */
  fla_qr2ut_t*     flash_qr2ut_cntl_leaf;
  fla_qr2ut_t*     flash_qr2ut_cntl;
  fla_blocksize_t* flash_qr2ut_var2_bsize;
  
  /* FLASH_CAQR2_UT_cntl_init */
  fla_caqr2ut_t*   flash_caqr2ut_cntl_leaf;
  fla_caqr2ut_t*   flash_caqr2ut_cntl;
  fla_blocksize_t* flash_caqr2ut_var2_bsize;
  
  /* FLASH_Apply_Q_UT_cntl_init */
  fla_apqut_t*        flash_apqut_cntl_leaf;
  fla_apqut_t*        flash_apqut_cntl;
  fla_apqut_t*        flash_apqut_cntl_blas;
  fla_blocksize_t*    flash_apqut_var1_bsize;
  fla_blocksize_t*    flash_apqut_var2_bsize;
  
  /* FLASH_Apply_Q2_UT_cntl_init */
  fla_apq2ut_t*    flash_apq2ut_cntl_leaf;
  fla_apq2ut_t*    flash_apq2ut_cntl_mid;
  fla_apq2ut_t*    flash_apq2ut_cntl;
  fla_blocksize_t* flash_apq2ut_var2_bsize;
  fla_blocksize_t* flash_apq2ut_var3_bsize;
  
  /* FLASH_Apply_CAQ2_UT_cntl_init */
  fla_apcaq2ut_t*  flash_apcaq2ut_cntl_leaf;
  fla_apcaq2ut_t*  flash_apcaq2ut_cntl_mid;
  fla_apcaq2ut_t*  flash_apcaq2ut_cntl;
  fla_blocksize_t* flash_apcaq2ut_var2_bsize;
  fla_blocksize_t* flash_apcaq2ut_var3_bsize;
  
  /* FLASH_Apply_QUD_UT_cntl_init */
  fla_apqudut_t*   flash_apqudut_cntl_leaf;
  fla_apqudut_t*   flash_apqudut_cntl_mid;
  fla_apqudut_t*   flash_apqudut_cntl;
  fla_blocksize_t* flash_apqudut_var2_bsize;
  fla_blocksize_t* flash_apqudut_var3_bsize;
  
  /* FLASH_Eig_gest_cntl_init */
  fla_eig_gest_t*     flash_eig_gest_cntl_leaf;
  fla_eig_gest_t*     flash_eig_gest_cntl;
  fla_blocksize_t*    flash_eig_gest_bsize;
  
  /* FLASH_Lyap_cntl_init */
  fla_lyap_t*         flash_lyap_cntl_leaf;
  fla_lyap_t*         flash_lyap_cntl;
  fla_blocksize_t*    flash_lyap_bsize;
  
  /* FLASH_SPDinv_cntl_init */
  fla_spdinv_t*       flash_spdinv_cntl;
  fla_blocksize_t*    flash_spdinv_size_cutoff;
  
  /* FLASH_QR_UT_cntl_init */
  fla_qrut_t*          flash_qrut_cntl_leaf;
  fla_qrut_t*          flash_qrut_cntl;
  
  fla_blocksize_t*     flash_qrut_var3_bsize;
  
  /* FLASH_QR_UT_inc_cntl_init */
  fla_qrutinc_t*       flash_qrutinc_cntl;
  fla_blocksize_t*     flash_qrutinc_var1_bsize;
  
  /* FLASH_LQ_UT_cntl_init */
  fla_lqut_t*          flash_lqut_cntl_leaf;
  fla_lqut_t*          flash_lqut_cntl;
  
  fla_blocksize_t*     flash_lqut_var3_bsize;
  
  /* FLASH_CAQR_UT_inc_cntl_init */
  fla_caqrutinc_t*     flash_caqrutinc_cntl;
  fla_blocksize_t*     flash_caqrutinc_var1_bsize;
  
  /* FLASH_Apply_Q_UT_inc_cntl_init */
  fla_apqutinc_t*      flash_apqutinc_cntl;
  fla_blocksize_t*     flash_apqutinc_var1_bsize;
  
  /* FLASH_Apply_CAQ_UT_inc_cntl_init */
  fla_apcaqutinc_t*      flash_apcaqutinc_cntl;
  fla_blocksize_t*       flash_apcaqutinc_var1_bsize;
  
  /* FLASH_UDdate_UT_cntl_init */
  fla_uddateut_t*  flash_uddateut_cntl_leaf;
  fla_uddateut_t*  flash_uddateut_cntl;
  fla_blocksize_t* flash_uddateut_var2_bsize;
  
  /* FLASH_UDdate_UT_inc_cntl_init */
  fla_uddateutinc_t*     flash_uddateutinc_cntl;
  fla_blocksize_t*       flash_uddateutinc_var1_bsize;
  
  /* FLASH_Apply_QUD_UT_inc_cntl_init */
  fla_apqudutinc_t*     flash_apqudutinc_cntl;
  fla_blocksize_t*      flash_apqudutinc_var1_bsize;
}FLA_Cntl_init_flash_s;

void FLA_Cntl_init_flash_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLA_Cntl_finalize_flash_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);


// --- Base library prototypes -------------------------------------------------
void FLASH_Transpose_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);

void FLASH_Transpose_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);


// --- Level-1 BLAS prototypes -------------------------------------------------
void FLASH_Axpy_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Axpyt_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Copy_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Copyt_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Copyr_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Scal_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Scalr_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);

void FLASH_Axpy_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Axpyt_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Copy_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Copyt_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Copyr_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Scal_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Scalr_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);


// --- Level-2 BLAS prototypes -------------------------------------------------
void FLASH_Gemv_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Trsv_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);

void FLASH_Gemv_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Trsv_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);


// --- Level-3 BLAS prototypes -------------------------------------------------
void FLASH_Gemm_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Hemm_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Herk_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Her2k_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Symm_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Syrk_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Syr2k_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Trmm_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Trsm_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);

void FLASH_Gemm_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Hemm_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Herk_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Her2k_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Symm_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Syrk_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Syr2k_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Trmm_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Trsm_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);


// --- LAPACK-level prototypes -------------------------------------------------
void FLASH_Apply_pivots_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Chol_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_LU_nopiv_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_LU_piv_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_LU_incpiv_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Trinv_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Ttmm_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_SPDinv_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Sylv_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Lyap_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_QR_UT_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_QR2_UT_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_LQ_UT_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_CAQR2_UT_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_UDdate_UT_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_QR_UT_inc_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_CAQR_UT_inc_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_UDdate_UT_inc_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_Q_UT_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_Q2_UT_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_CAQ2_UT_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_QUD_UT_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_Q_UT_inc_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_CAQ_UT_inc_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_QUD_UT_inc_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Eig_gest_cntl_init_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);

void FLASH_Apply_pivots_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Chol_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_LU_nopiv_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_LU_piv_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_LU_incpiv_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Trinv_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Ttmm_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_SPDinv_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Sylv_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Lyap_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_QR_UT_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_QR2_UT_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_LQ_UT_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_CAQR2_UT_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_UDdate_UT_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_QR_UT_inc_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_CAQR_UT_inc_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_UDdate_UT_inc_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_Q_UT_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_Q2_UT_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_CAQ2_UT_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_QUD_UT_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_Q_UT_inc_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_CAQ_UT_inc_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Apply_QUD_UT_inc_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);
void FLASH_Eig_gest_cntl_finalize_ts(FLA_Cntl_init_flash_s *FLA_Cntl_init_flash_i);

#endif

void FLA_Cntl_init_flash( void );
void FLA_Cntl_finalize_flash( void );


// --- Base library prototypes -------------------------------------------------
void FLASH_Transpose_cntl_init( void );

void FLASH_Transpose_cntl_finalize( void );


// --- Level-1 BLAS prototypes -------------------------------------------------
void FLASH_Axpy_cntl_init( void );
void FLASH_Axpyt_cntl_init( void );
void FLASH_Copy_cntl_init( void );
void FLASH_Copyt_cntl_init( void );
void FLASH_Copyr_cntl_init( void );
void FLASH_Scal_cntl_init( void );
void FLASH_Scalr_cntl_init( void );

void FLASH_Axpy_cntl_finalize( void );
void FLASH_Axpyt_cntl_finalize( void );
void FLASH_Copy_cntl_finalize( void );
void FLASH_Copyt_cntl_finalize( void );
void FLASH_Copyr_cntl_finalize( void );
void FLASH_Scal_cntl_finalize( void );
void FLASH_Scalr_cntl_finalize( void );


// --- Level-2 BLAS prototypes -------------------------------------------------
void FLASH_Gemv_cntl_init( void );
void FLASH_Trsv_cntl_init( void );

void FLASH_Gemv_cntl_finalize( void );
void FLASH_Trsv_cntl_finalize( void );


// --- Level-3 BLAS prototypes -------------------------------------------------
void FLASH_Gemm_cntl_init( void );
void FLASH_Hemm_cntl_init( void );
void FLASH_Herk_cntl_init( void );
void FLASH_Her2k_cntl_init( void );
void FLASH_Symm_cntl_init( void );
void FLASH_Syrk_cntl_init( void );
void FLASH_Syr2k_cntl_init( void );
void FLASH_Trmm_cntl_init( void );
void FLASH_Trsm_cntl_init( void );

void FLASH_Gemm_cntl_finalize( void );
void FLASH_Hemm_cntl_finalize( void );
void FLASH_Herk_cntl_finalize( void );
void FLASH_Her2k_cntl_finalize( void );
void FLASH_Symm_cntl_finalize( void );
void FLASH_Syrk_cntl_finalize( void );
void FLASH_Syr2k_cntl_finalize( void );
void FLASH_Trmm_cntl_finalize( void );
void FLASH_Trsm_cntl_finalize( void );


// --- LAPACK-level prototypes -------------------------------------------------
void FLASH_Apply_pivots_cntl_init( void );
void FLASH_Chol_cntl_init( void );
void FLASH_LU_nopiv_cntl_init( void );
void FLASH_LU_piv_cntl_init( void );
void FLASH_LU_incpiv_cntl_init( void );
void FLASH_Trinv_cntl_init( void );
void FLASH_Ttmm_cntl_init( void );
void FLASH_SPDinv_cntl_init( void );
void FLASH_Sylv_cntl_init( void );
void FLASH_Lyap_cntl_init( void );
void FLASH_QR_UT_cntl_init( void );
void FLASH_QR2_UT_cntl_init( void );
void FLASH_LQ_UT_cntl_init( void );
void FLASH_CAQR2_UT_cntl_init( void );
void FLASH_UDdate_UT_cntl_init( void );
void FLASH_QR_UT_inc_cntl_init( void );
void FLASH_CAQR_UT_inc_cntl_init( void );
void FLASH_UDdate_UT_inc_cntl_init( void );
void FLASH_Apply_Q_UT_cntl_init( void );
void FLASH_Apply_Q2_UT_cntl_init( void );
void FLASH_Apply_CAQ2_UT_cntl_init( void );
void FLASH_Apply_QUD_UT_cntl_init( void );
void FLASH_Apply_Q_UT_inc_cntl_init( void );
void FLASH_Apply_CAQ_UT_inc_cntl_init( void );
void FLASH_Apply_QUD_UT_inc_cntl_init( void );
void FLASH_Eig_gest_cntl_init( void );

void FLASH_Apply_pivots_cntl_finalize( void );
void FLASH_Chol_cntl_finalize( void );
void FLASH_LU_nopiv_cntl_finalize( void );
void FLASH_LU_piv_cntl_finalize( void );
void FLASH_LU_incpiv_cntl_finalize( void );
void FLASH_Trinv_cntl_finalize( void );
void FLASH_Ttmm_cntl_finalize( void );
void FLASH_SPDinv_cntl_finalize( void );
void FLASH_Sylv_cntl_finalize( void );
void FLASH_Lyap_cntl_finalize( void );
void FLASH_QR_UT_cntl_finalize( void );
void FLASH_QR2_UT_cntl_finalize( void );
void FLASH_LQ_UT_cntl_finalize( void );
void FLASH_CAQR2_UT_cntl_finalize( void );
void FLASH_UDdate_UT_cntl_finalize( void );
void FLASH_QR_UT_inc_cntl_finalize( void );
void FLASH_CAQR_UT_inc_cntl_finalize( void );
void FLASH_UDdate_UT_inc_cntl_finalize( void );
void FLASH_Apply_Q_UT_cntl_finalize( void );
void FLASH_Apply_Q2_UT_cntl_finalize( void );
void FLASH_Apply_CAQ2_UT_cntl_finalize( void );
void FLASH_Apply_QUD_UT_cntl_finalize( void );
void FLASH_Apply_Q_UT_inc_cntl_finalize( void );
void FLASH_Apply_CAQ_UT_inc_cntl_finalize( void );
void FLASH_Apply_QUD_UT_inc_cntl_finalize( void );
void FLASH_Eig_gest_cntl_finalize( void );

