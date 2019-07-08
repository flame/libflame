/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES

typedef struct FLA_Transpose_cntl_init_s
{
  fla_swap_t*      fla_swap_cntl_panel;
  fla_swap_t*      fla_swap_cntl_blas;
  
  fla_tpose_t*     fla_tpose_cntl;
  fla_tpose_t*     fla_tpose_cntl_unb;
  fla_blocksize_t* fla_tpose_bsize;
  fla_blocksize_t* fla_tpose_swap_bsize;
}FLA_Transpose_cntl_init_s;

typedef struct FLA_Cntl_init_flamec_s
{
  FLA_Transpose_cntl_init_s *FLA_Transpose_cntl_init_i;

  /* FLA_Axpy_cntl_init */
  fla_axpy_t* fla_axpy_cntl_blas;
  
  /* FLA_Axpyt_cntl_init */
  fla_axpyt_t* fla_axpyt_cntl_blas;
  
  /* FLA_Copy_cntl_init */
  fla_copy_t* fla_copy_cntl_blas;
  
  /* FLA_Copyt_cntl_init */
  fla_copyt_t* fla_copyt_cntl_blas;
  
  /* FLA_Copyr_cntl_init */
  fla_copyr_t* fla_copyr_cntl_blas;
  
  /* FLA_Scal_cntl_init */
  fla_scal_t* fla_scal_cntl_blas;
  
  /* FLA_Scalr_cntl_init */
  fla_scalr_t* fla_scalr_cntl_blas;
  
  /* FLA_Gemv_cntl_init */
  fla_gemv_t* fla_gemv_cntl_blas;
  
  /* FLA_Trsv_cntl_init */
  fla_trsv_t* fla_trsv_cntl_blas;
  
  /* FLA_Gemm_cntl_init */
  fla_gemm_t*      fla_gemm_cntl_blas;
  
  fla_gemm_t*      fla_gemm_cntl_pb_bb;
  fla_gemm_t*      fla_gemm_cntl_bp_bb;
  fla_gemm_t*      fla_gemm_cntl_ip_bb;
  
  fla_gemm_t*      fla_gemm_cntl_mp_ip;
  fla_gemm_t*      fla_gemm_cntl_mp_ip_bb;
  fla_gemm_t*      fla_gemm_cntl_op_bp;
  fla_gemm_t*      fla_gemm_cntl_op_bp_bb;
  fla_gemm_t*      fla_gemm_cntl_pm_ip;
  fla_gemm_t*      fla_gemm_cntl_pm_ip_bb;
  fla_gemm_t*      fla_gemm_cntl_op_pb;
  fla_gemm_t*      fla_gemm_cntl_op_pb_bb;
  fla_gemm_t*      fla_gemm_cntl_mp_pb;
  fla_gemm_t*      fla_gemm_cntl_mp_pb_bb;
  fla_gemm_t*      fla_gemm_cntl_pm_bp;
  fla_gemm_t*      fla_gemm_cntl_pm_bp_bb;
  
  fla_gemm_t*      fla_gemm_cntl_mm_pm;
  fla_gemm_t*      fla_gemm_cntl_mm_pm_ip;
  fla_gemm_t*      fla_gemm_cntl_mm_pm_ip_bb;
  fla_gemm_t*      fla_gemm_cntl_mm_mp;
  fla_gemm_t*      fla_gemm_cntl_mm_mp_ip;
  fla_gemm_t*      fla_gemm_cntl_mm_mp_ip_bb;
  fla_gemm_t*      fla_gemm_cntl_mm_op;
  fla_gemm_t*      fla_gemm_cntl_mm_op_bp;
  fla_gemm_t*      fla_gemm_cntl_mm_op_bp_bb;
  
  fla_blocksize_t* fla_gemm_var1_bsize;
  fla_blocksize_t* fla_gemm_var3_bsize;
  fla_blocksize_t* fla_gemm_var5_bsize;
  
  /* FLA_Hemm_cntl_init */
  
  fla_hemm_t*        fla_hemm_cntl_blas;
  fla_hemm_t*        fla_hemm_cntl_bp;
  fla_hemm_t*        fla_hemm_cntl_mp;
  fla_hemm_t*        fla_hemm_cntl_mm;
  fla_blocksize_t*   fla_hemm_var1_bsize;
  fla_blocksize_t*   fla_hemm_var9_bsize;
  
  /* FLA_Herk_cntl_init */
  fla_herk_t*         fla_herk_cntl_blas;
  fla_herk_t*         fla_herk_cntl_ip;
  fla_herk_t*         fla_herk_cntl_op;
  fla_herk_t*         fla_herk_cntl_mm;
  fla_blocksize_t*    fla_herk_var2_bsize;
  fla_blocksize_t*    fla_herk_var5_bsize;
  
  /* FLA_Her2k_cntl_init */
  fla_her2k_t*        fla_her2k_cntl_blas;
  fla_her2k_t*        fla_her2k_cntl_ip;
  fla_her2k_t*        fla_her2k_cntl_op;
  fla_her2k_t*        fla_her2k_cntl_mm;
  fla_blocksize_t*    fla_her2k_var3_bsize;
  fla_blocksize_t*    fla_her2k_var9_bsize;
  
  /* FLA_Symm_cntl_init */
  fla_symm_t*        fla_symm_cntl_blas;
  fla_symm_t*        fla_symm_cntl_bp;
  fla_symm_t*        fla_symm_cntl_mp;
  fla_symm_t*        fla_symm_cntl_mm;
  fla_blocksize_t*   fla_symm_var1_bsize;
  fla_blocksize_t*   fla_symm_var9_bsize;
  
  /* FLA_Syrk_cntl_init */
  fla_syrk_t*         fla_syrk_cntl_blas;
  fla_syrk_t*         fla_syrk_cntl_ip;
  fla_syrk_t*         fla_syrk_cntl_op;
  fla_syrk_t*         fla_syrk_cntl_mm;
  fla_blocksize_t*    fla_syrk_var2_bsize;
  fla_blocksize_t*    fla_syrk_var5_bsize;
  
  /* FLA_Syr2k_cntl_init */
  fla_syr2k_t*        fla_syr2k_cntl_blas;
  fla_syr2k_t*        fla_syr2k_cntl_ip;
  fla_syr2k_t*        fla_syr2k_cntl_op;
  fla_syr2k_t*        fla_syr2k_cntl_mm;
  fla_blocksize_t*    fla_syr2k_var3_bsize;
  fla_blocksize_t*    fla_syr2k_var9_bsize;
  
  /* FLA_Trmm_cntl_init */
  fla_trmm_t*        fla_trmm_cntl_blas;
  fla_trmm_t*        fla_trmm_cntl_bp;
  fla_trmm_t*        fla_trmm_cntl_mp;
  fla_trmm_t*        fla_trmm_cntl_mm;
  fla_blocksize_t*   fla_trmm_var1_bsize;
  fla_blocksize_t*   fla_trmm_var3_bsize;
  
  /* FLA_Trsm_cntl_init */
  fla_trsm_t*        fla_trsm_cntl_blas;
  fla_trsm_t*        fla_trsm_cntl_bp;
  fla_trsm_t*        fla_trsm_cntl_mp;
  fla_trsm_t*        fla_trsm_cntl_mm;
  fla_blocksize_t*   fla_trsm_var2_bsize;
  fla_blocksize_t*   fla_trsm_var3_bsize;
  
  /* FLA_Apply_pivots_cntl_init */
  fla_appiv_t* fla_appiv_cntl_leaf;
  
  /* FLA_Chol_cntl_init */
  fla_chol_t*        fla_chol_cntl;
  fla_chol_t*        fla_chol_cntl2;
  
  fla_chol_t*        fla_chol_cntl_in;
  fla_chol_t*        fla_chol_cntl_leaf;
  fla_blocksize_t*   fla_chol_var3_bsize;
  fla_blocksize_t*   fla_chol_var3_bsize_in;
  
  /* FLA_LU_nopiv_cntl_init */
  fla_lu_t*          fla_lu_nopiv_cntl;
  fla_lu_t*          fla_lu_nopiv_cntl2;
  
  fla_lu_t*          fla_lu_nopiv_cntl_in;
  fla_lu_t*          fla_lu_nopiv_cntl_leaf;
  fla_blocksize_t*   fla_lu_nopiv_var5_bsize;
  fla_blocksize_t*   fla_lu_nopiv_var5_bsize_in;
  
  /* FLA_LU_piv_cntl_init */
  fla_lu_t*           fla_lu_piv_cntl;
  fla_lu_t*           fla_lu_piv_cntl2;
  
  fla_lu_t*           fla_lu_piv_cntl_in;
  fla_lu_t*           fla_lu_piv_cntl_leaf;
  fla_blocksize_t*    fla_lu_piv_var5_bsize;
  fla_blocksize_t*    fla_lu_piv_var5_bsize_in;
  
  /* FLA_Trinv_cntl_init */
  fla_trinv_t*       fla_trinv_cntl_leaf;
  fla_trinv_t*       fla_trinv_cntl;
  fla_blocksize_t*   fla_trinv_var3_bsize;
  
  /* FLA_Ttmm_cntl_init */
  fla_ttmm_t*        fla_ttmm_cntl_leaf;
  fla_ttmm_t*        fla_ttmm_cntl;
  fla_blocksize_t*   fla_ttmm_var1_bsize;
  
  /* FLA_Sylv_cntl_init */
  fla_sylv_t*        fla_sylv_cntl_leaf;
  fla_sylv_t*        fla_sylv_cntl_mb;
  fla_sylv_t*        fla_sylv_cntl;
  fla_blocksize_t*   fla_sylv_bsize;
  
  /* FLA_QR2_UT_cntl_init */
  fla_qr2ut_t*       fla_qr2ut_cntl_unb;
  fla_qr2ut_t*       fla_qr2ut_cntl_leaf;
  fla_blocksize_t*   fla_qr2ut_var1_bsize;
  
  /* FLA_CAQR2_UT_cntl_init */
  fla_caqr2ut_t*     fla_caqr2ut_cntl_unb;
  fla_caqr2ut_t*     fla_caqr2ut_cntl_leaf;
  fla_blocksize_t*   fla_caqr2ut_var1_bsize;
  
  /* FLA_Apply_Q_UT_cntl_init */
  fla_apqut_t*        fla_apqut_cntl_leaf;
  fla_apqut_t*        fla_apqut_cntl;
  fla_blocksize_t*    fla_apqut_var1_bsize;
  fla_blocksize_t*    fla_apqut_var2_bsize;
  
  /* FLA_Apply_Q2_UT_cntl_init */
  fla_apq2ut_t*       fla_apq2ut_cntl_leaf;
  fla_blocksize_t*    fla_apq2ut_var1_bsize;
  
  /* FLA_Apply_CAQ2_UT_cntl_init */
  fla_apcaq2ut_t*     fla_apcaq2ut_cntl_leaf;
  fla_blocksize_t*    fla_apcaq2ut_var1_bsize;
  
  /* FLA_Apply_QUD_UT_cntl_init */
  fla_apqudut_t*      fla_apqudut_cntl_leaf;
  fla_blocksize_t*    fla_apqudut_var1_bsize;
  
  /* FLA_Eig_gest_cntl_init */
  fla_eig_gest_t*     fla_eig_gest_ix_cntl;
  fla_eig_gest_t*     fla_eig_gest_nx_cntl;
  fla_eig_gest_t*     fla_eig_gest_ix_cntl_leaf;
  fla_eig_gest_t*     fla_eig_gest_nx_cntl_leaf;
  fla_blocksize_t*    fla_eig_gest_var1_bsize;
  
  /* FLA_Lyap_cntl_init */
  fla_lyap_t*         fla_lyap_cntl_leaf;
  fla_lyap_t*         fla_lyap_cntl;
  fla_blocksize_t*    fla_lyap_bsize;
  
  /* FLA_SPDinv_cntl_init */
  fla_spdinv_t*       fla_spdinv_cntl;
  fla_blocksize_t*    fla_spdinv_size_cutoff;
  
  /* FLA_QR_UT_cntl_init */
  fla_qrut_t*         fla_qrut_cntl_unb;
  fla_qrut_t*         fla_qrut_cntl_leaf;
  
  fla_qrut_t*         fla_qrut_piv_cntl_unb;
  fla_qrut_t*         fla_qrut_piv_cntl_leaf;
  
  fla_blocksize_t*    fla_qrut_var1_bsize_leaf;
  
  /* FLA_LQ_UT_cntl_init */
  fla_lqut_t*         fla_lqut_cntl_unb;
  fla_lqut_t*         fla_lqut_cntl_leaf;
  
  fla_blocksize_t*    fla_lqut_var1_bsize_leaf;
  
  /* FLA_UDdate_UT_cntl_init */
  fla_uddateut_t*       fla_uddateut_cntl_unb;
  fla_uddateut_t*       fla_uddateut_cntl_leaf;
  fla_blocksize_t*      fla_uddateut_var1_bsize;
  
  /* FLA_Hess_UT_cntl_init */
  fla_hessut_t*       fla_hessut_cntl_leaf;
  
  fla_blocksize_t*    fla_hessut_bsize_leaf;
  
  /* FLA_Tridiag_UT_cntl_init */
  fla_tridiagut_t*    fla_tridiagut_cntl_fused;
  fla_tridiagut_t*    fla_tridiagut_cntl_nofus;
  fla_tridiagut_t*    fla_tridiagut_cntl_plain;
  
  fla_blocksize_t*    fla_tridiagut_bsize_leaf;
  
  /* FLA_Bidiag_UT_cntl_init */
  fla_bidiagut_t*    fla_bidiagut_cntl_fused;
  fla_bidiagut_t*    fla_bidiagut_cntl_nofus;
  fla_bidiagut_t*    fla_bidiagut_cntl_plain;
  
  fla_blocksize_t*   fla_bidiagut_bsize_leaf;

  // unsigned int *fla_error_checking_level;
}FLA_Cntl_init_flamec_s;

void FLA_Cntl_init_flamec_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Cntl_finalize_flamec_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);


// --- Base library prototypes -------------------------------------------------
void FLA_Transpose_cntl_init_ts(FLA_Transpose_cntl_init_s* FLA_Transpose_cntl_init_i);

void FLA_Transpose_cntl_finalize_ts(FLA_Transpose_cntl_init_s* FLA_Transpose_cntl_init_i);


// --- Level-1 BLAS prototypes -------------------------------------------------
void FLA_Axpy_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Axpyt_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Copy_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Copyt_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Copyr_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Scal_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Scalr_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);

void FLA_Axpy_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Axpyt_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Copy_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Copyt_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Copyr_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Scal_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Scalr_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);


// --- Level-2 BLAS prototypes -------------------------------------------------
void FLA_Gemv_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Trsv_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);

void FLA_Gemv_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Trsv_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);


// --- Level-3 BLAS prototypes -------------------------------------------------
void FLA_Gemm_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Hemm_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Herk_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Her2k_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Symm_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Syrk_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Syr2k_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Trmm_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Trsm_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);

void FLA_Gemm_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Hemm_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Herk_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Her2k_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Symm_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Syrk_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Syr2k_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Trmm_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Trsm_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);


// --- LAPACK-level prototypes -------------------------------------------------
void FLA_Apply_pivots_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Chol_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_LU_piv_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_LU_nopiv_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_QR_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_QR2_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_LQ_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_CAQR2_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_UDdate_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Hess_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Tridiag_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Bidiag_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Trinv_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Ttmm_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Sylv_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Lyap_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_SPDinv_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Apply_Q_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Apply_Q2_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Apply_CAQ2_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Apply_QUD_UT_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Eig_gest_cntl_init_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);

void FLA_Apply_pivots_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Chol_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_LU_piv_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_LU_nopiv_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_QR_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_QR2_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_LQ_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_CAQR2_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_UDdate_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Hess_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Tridiag_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Bidiag_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Trinv_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Ttmm_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Sylv_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Lyap_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_SPDinv_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Apply_Q_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Apply_Q2_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Apply_CAQ2_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Apply_QUD_UT_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);
void FLA_Eig_gest_cntl_finalize_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i);

#endif

void FLA_Cntl_init_flamec( void );
void FLA_Cntl_finalize_flamec( void );


// --- Base library prototypes -------------------------------------------------
void FLA_Transpose_cntl_init( void );

void FLA_Transpose_cntl_finalize( void );


// --- Level-1 BLAS prototypes -------------------------------------------------
void FLA_Axpy_cntl_init( void );
void FLA_Axpyt_cntl_init( void );
void FLA_Copy_cntl_init( void );
void FLA_Copyt_cntl_init( void );
void FLA_Copyr_cntl_init( void );
void FLA_Scal_cntl_init( void );
void FLA_Scalr_cntl_init( void );

void FLA_Axpy_cntl_finalize( void );
void FLA_Axpyt_cntl_finalize( void );
void FLA_Copy_cntl_finalize( void );
void FLA_Copyt_cntl_finalize( void );
void FLA_Copyr_cntl_finalize( void );
void FLA_Scal_cntl_finalize( void );
void FLA_Scalr_cntl_finalize( void );


// --- Level-2 BLAS prototypes -------------------------------------------------
void FLA_Gemv_cntl_init( void );
void FLA_Trsv_cntl_init( void );

void FLA_Gemv_cntl_finalize( void );
void FLA_Trsv_cntl_finalize( void );


// --- Level-3 BLAS prototypes -------------------------------------------------
void FLA_Gemm_cntl_init( void );
void FLA_Hemm_cntl_init( void );
void FLA_Herk_cntl_init( void );
void FLA_Her2k_cntl_init( void );
void FLA_Symm_cntl_init( void );
void FLA_Syrk_cntl_init( void );
void FLA_Syr2k_cntl_init( void );
void FLA_Trmm_cntl_init( void );
void FLA_Trsm_cntl_init( void );

void FLA_Gemm_cntl_finalize( void );
void FLA_Hemm_cntl_finalize( void );
void FLA_Herk_cntl_finalize( void );
void FLA_Her2k_cntl_finalize( void );
void FLA_Symm_cntl_finalize( void );
void FLA_Syrk_cntl_finalize( void );
void FLA_Syr2k_cntl_finalize( void );
void FLA_Trmm_cntl_finalize( void );
void FLA_Trsm_cntl_finalize( void );


// --- LAPACK-level prototypes -------------------------------------------------
void FLA_Apply_pivots_cntl_init( void );
void FLA_Chol_cntl_init( void );
void FLA_LU_piv_cntl_init( void );
void FLA_LU_nopiv_cntl_init( void );
void FLA_QR_UT_cntl_init( void );
void FLA_QR2_UT_cntl_init( void );
void FLA_LQ_UT_cntl_init( void );
void FLA_CAQR2_UT_cntl_init( void );
void FLA_UDdate_UT_cntl_init( void );
void FLA_Hess_UT_cntl_init( void );
void FLA_Tridiag_UT_cntl_init( void );
void FLA_Bidiag_UT_cntl_init( void );
void FLA_Trinv_cntl_init( void );
void FLA_Ttmm_cntl_init( void );
void FLA_Sylv_cntl_init( void );
void FLA_Lyap_cntl_init( void );
void FLA_SPDinv_cntl_init( void );
void FLA_Apply_Q_UT_cntl_init( void );
void FLA_Apply_Q2_UT_cntl_init( void );
void FLA_Apply_CAQ2_UT_cntl_init( void );
void FLA_Apply_QUD_UT_cntl_init( void );
void FLA_Eig_gest_cntl_init( void );

void FLA_Apply_pivots_cntl_finalize( void );
void FLA_Chol_cntl_finalize( void );
void FLA_LU_piv_cntl_finalize( void );
void FLA_LU_nopiv_cntl_finalize( void );
void FLA_QR_UT_cntl_finalize( void );
void FLA_QR2_UT_cntl_finalize( void );
void FLA_LQ_UT_cntl_finalize( void );
void FLA_CAQR2_UT_cntl_finalize( void );
void FLA_UDdate_UT_cntl_finalize( void );
void FLA_Hess_UT_cntl_finalize( void );
void FLA_Tridiag_UT_cntl_finalize( void );
void FLA_Bidiag_UT_cntl_finalize( void );
void FLA_Trinv_cntl_finalize( void );
void FLA_Ttmm_cntl_finalize( void );
void FLA_Sylv_cntl_finalize( void );
void FLA_Lyap_cntl_finalize( void );
void FLA_SPDinv_cntl_finalize( void );
void FLA_Apply_Q_UT_cntl_finalize( void );
void FLA_Apply_Q2_UT_cntl_finalize( void );
void FLA_Apply_CAQ2_UT_cntl_finalize( void );
void FLA_Apply_QUD_UT_cntl_finalize( void );
void FLA_Eig_gest_cntl_finalize( void );

