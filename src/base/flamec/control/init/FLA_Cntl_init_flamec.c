/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
void FLA_Cntl_init_flamec_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i)
{
  FLA_cntl_flamec_init_i->FLA_Transpose_cntl_init_i = (FLA_Transpose_cntl_init_s *) FLA_malloc(sizeof(FLA_Transpose_cntl_init_s));

  FLA_Transpose_cntl_init_ts(FLA_cntl_flamec_init_i->FLA_Transpose_cntl_init_i);

  // Level-1 BLAS
  FLA_Axpy_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Axpyt_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Copy_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Copyt_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Copyr_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Scal_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Scalr_cntl_init_ts(FLA_cntl_flamec_init_i);

  // Level-2 BLAS
  FLA_Gemv_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Trsv_cntl_init_ts(FLA_cntl_flamec_init_i);

  // Level-3 BLAS
  // Note gemm must be first since it is used by all other level-3 BLAS.
  FLA_Gemm_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Hemm_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Herk_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Her2k_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Symm_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Syrk_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Syr2k_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Trmm_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Trsm_cntl_init_ts(FLA_cntl_flamec_init_i);

  // LAPACK-level
  // These require level-3 BLAS operations to be initialized.
  FLA_Apply_pivots_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Chol_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_LU_nopiv_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_LU_piv_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Trinv_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Ttmm_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Sylv_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_QR2_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_CAQR2_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Apply_Q_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Apply_Q2_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Apply_CAQ2_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Apply_QUD_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Eig_gest_cntl_init_ts(FLA_cntl_flamec_init_i);

  // Compound LAPACK operations
  // These require previous LAPACK operations to already be initialized.
  FLA_Lyap_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_SPDinv_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_QR_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_LQ_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_UDdate_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Hess_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Tridiag_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
  FLA_Bidiag_UT_cntl_init_ts(FLA_cntl_flamec_init_i);
}

void FLA_Cntl_finalize_flamec_ts(FLA_Cntl_init_flamec_s *FLA_cntl_flamec_init_i)
{
  FLA_Transpose_cntl_finalize_ts(FLA_cntl_flamec_init_i->FLA_Transpose_cntl_init_i);
  FLA_free(FLA_cntl_flamec_init_i->FLA_Transpose_cntl_init_i);

  // Level-1 BLAS
  FLA_Axpy_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Axpyt_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Copy_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Copyt_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Copyr_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Scal_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Scalr_cntl_finalize_ts(FLA_cntl_flamec_init_i);

  // Level-2 BLAS
  FLA_Gemv_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Trsv_cntl_finalize_ts(FLA_cntl_flamec_init_i);

  // Level-3 BLAS
  FLA_Gemm_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Hemm_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Herk_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Her2k_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Symm_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Syrk_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Syr2k_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Trmm_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Trsm_cntl_finalize_ts(FLA_cntl_flamec_init_i);

  // LAPACK-level
  FLA_Apply_pivots_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Chol_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_LU_nopiv_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_LU_piv_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Trinv_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Ttmm_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Sylv_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_QR2_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_CAQR2_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Apply_Q_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Apply_Q2_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Apply_CAQ2_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Apply_QUD_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Eig_gest_cntl_finalize_ts(FLA_cntl_flamec_init_i);

  // Compound LAPACK operations
  FLA_Lyap_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_SPDinv_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_QR_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_LQ_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_UDdate_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Hess_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Tridiag_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
  FLA_Bidiag_UT_cntl_finalize_ts(FLA_cntl_flamec_init_i);
}
#endif

void FLA_Cntl_init_flamec()
{
  FLA_Transpose_cntl_init();

  // Level-1 BLAS
  FLA_Axpy_cntl_init();
  FLA_Axpyt_cntl_init();
  FLA_Copy_cntl_init();
  FLA_Copyt_cntl_init();
  FLA_Copyr_cntl_init();
  FLA_Scal_cntl_init();
  FLA_Scalr_cntl_init();

  // Level-2 BLAS
  FLA_Gemv_cntl_init();
  FLA_Trsv_cntl_init();

  // Level-3 BLAS
  // Note gemm must be first since it is used by all other level-3 BLAS.
  FLA_Gemm_cntl_init();
  FLA_Hemm_cntl_init();
  FLA_Herk_cntl_init();
  FLA_Her2k_cntl_init();
  FLA_Symm_cntl_init();
  FLA_Syrk_cntl_init();
  FLA_Syr2k_cntl_init();
  FLA_Trmm_cntl_init();
  FLA_Trsm_cntl_init();

  // LAPACK-level
  // These require level-3 BLAS operations to be initialized.
  FLA_Apply_pivots_cntl_init();
  FLA_Chol_cntl_init();
  FLA_LU_nopiv_cntl_init();
  FLA_LU_piv_cntl_init();
  FLA_Trinv_cntl_init();
  FLA_Ttmm_cntl_init();
  FLA_Sylv_cntl_init();
  FLA_QR2_UT_cntl_init();
  FLA_CAQR2_UT_cntl_init();
  FLA_Apply_Q_UT_cntl_init();
  FLA_Apply_Q2_UT_cntl_init();
  FLA_Apply_CAQ2_UT_cntl_init();
  FLA_Apply_QUD_UT_cntl_init();
  FLA_Eig_gest_cntl_init();

  // Compound LAPACK operations
  // These require previous LAPACK operations to already be initialized.
  FLA_Lyap_cntl_init();
  FLA_SPDinv_cntl_init();
  FLA_QR_UT_cntl_init();
  FLA_LQ_UT_cntl_init();
  FLA_UDdate_UT_cntl_init();
  FLA_Hess_UT_cntl_init();
  FLA_Tridiag_UT_cntl_init();
  FLA_Bidiag_UT_cntl_init();
}

void FLA_Cntl_finalize_flamec()
{
  FLA_Transpose_cntl_finalize();

  // Level-1 BLAS
  FLA_Axpy_cntl_finalize();
  FLA_Axpyt_cntl_finalize();
  FLA_Copy_cntl_finalize();
  FLA_Copyt_cntl_finalize();
  FLA_Copyr_cntl_finalize();
  FLA_Scal_cntl_finalize();
  FLA_Scalr_cntl_finalize();

  // Level-2 BLAS
  FLA_Gemv_cntl_finalize();
  FLA_Trsv_cntl_finalize();

  // Level-3 BLAS
  FLA_Gemm_cntl_finalize();
  FLA_Hemm_cntl_finalize();
  FLA_Herk_cntl_finalize();
  FLA_Her2k_cntl_finalize();
  FLA_Symm_cntl_finalize();
  FLA_Syrk_cntl_finalize();
  FLA_Syr2k_cntl_finalize();
  FLA_Trmm_cntl_finalize();
  FLA_Trsm_cntl_finalize();

  // LAPACK-level
  FLA_Apply_pivots_cntl_finalize();
  FLA_Chol_cntl_finalize();
  FLA_LU_nopiv_cntl_finalize();
  FLA_LU_piv_cntl_finalize();
  FLA_Trinv_cntl_finalize();
  FLA_Ttmm_cntl_finalize();
  FLA_Sylv_cntl_finalize();
  FLA_QR2_UT_cntl_finalize();
  FLA_CAQR2_UT_cntl_finalize();
  FLA_Apply_Q_UT_cntl_finalize();
  FLA_Apply_Q2_UT_cntl_finalize();
  FLA_Apply_CAQ2_UT_cntl_finalize();
  FLA_Apply_QUD_UT_cntl_finalize();
  FLA_Eig_gest_cntl_finalize();

  // Compound LAPACK operations
  FLA_Lyap_cntl_finalize();
  FLA_SPDinv_cntl_finalize();
  FLA_QR_UT_cntl_finalize();
  FLA_LQ_UT_cntl_finalize();
  FLA_UDdate_UT_cntl_finalize();
  FLA_Hess_UT_cntl_finalize();
  FLA_Tridiag_UT_cntl_finalize();
  FLA_Bidiag_UT_cntl_finalize();
}

