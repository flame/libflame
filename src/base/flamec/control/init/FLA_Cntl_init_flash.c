/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

void FLA_Cntl_init_flash()
{
  // Level-1 BLAS
  FLASH_Axpy_cntl_init();
  FLASH_Axpyt_cntl_init();
  FLASH_Copy_cntl_init();
  FLASH_Copyt_cntl_init();
  FLASH_Copyr_cntl_init();
  FLASH_Scal_cntl_init();
  FLASH_Scalr_cntl_init();

  // Level-2 BLAS
  FLASH_Gemv_cntl_init();
  FLASH_Trsv_cntl_init();

  // Level-3 BLAS
  // Note gemm must be first since it is used by all other level-3 BLAS.
  FLASH_Gemm_cntl_init();
  FLASH_Hemm_cntl_init();
  FLASH_Herk_cntl_init();
  FLASH_Her2k_cntl_init();
  FLASH_Symm_cntl_init();
  FLASH_Syrk_cntl_init();
  FLASH_Syr2k_cntl_init();
  FLASH_Trmm_cntl_init();
  FLASH_Trsm_cntl_init();

  // LAPACK-level
  // These require level-3 BLAS operations to be initialized.
  FLASH_Apply_pivots_cntl_init();
  FLASH_Chol_cntl_init();
  FLASH_LU_nopiv_cntl_init();
  FLASH_LU_piv_cntl_init();
  FLASH_LU_incpiv_cntl_init();
  FLASH_Trinv_cntl_init();
  FLASH_Ttmm_cntl_init();
  FLASH_Sylv_cntl_init();
  FLASH_QR2_UT_cntl_init();
  FLASH_CAQR2_UT_cntl_init();
  FLASH_Apply_Q_UT_cntl_init();
  FLASH_Apply_Q2_UT_cntl_init();
  FLASH_Apply_CAQ2_UT_cntl_init();
  FLASH_Apply_QUD_UT_cntl_init();
  FLASH_Eig_gest_cntl_init();

  // Compound LAPACK operations
  // These require previous LAPACK operations to already be initialized.
  FLASH_Lyap_cntl_init();
  FLASH_SPDinv_cntl_init();
  FLASH_QR_UT_cntl_init();
  FLASH_QR_UT_inc_cntl_init();
  FLASH_LQ_UT_cntl_init();
  FLASH_CAQR_UT_inc_cntl_init();
  FLASH_Apply_Q_UT_inc_cntl_init();
  FLASH_Apply_CAQ_UT_inc_cntl_init();
  FLASH_UDdate_UT_cntl_init();
  FLASH_UDdate_UT_inc_cntl_init();
  FLASH_Apply_QUD_UT_inc_cntl_init();
}

void FLA_Cntl_finalize_flash()
{
  // Level-1 BLAS
  FLASH_Axpy_cntl_finalize();
  FLASH_Axpyt_cntl_finalize();
  FLASH_Copy_cntl_finalize();
  FLASH_Copyt_cntl_finalize();
  FLASH_Copyr_cntl_finalize();
  FLASH_Scal_cntl_finalize();
  FLASH_Scalr_cntl_finalize();

  // Level-2 BLAS
  FLASH_Gemv_cntl_finalize();
  FLASH_Trsv_cntl_finalize();

  // Level-3 BLAS
  FLASH_Gemm_cntl_finalize();
  FLASH_Hemm_cntl_finalize();
  FLASH_Herk_cntl_finalize();
  FLASH_Her2k_cntl_finalize();
  FLASH_Symm_cntl_finalize();
  FLASH_Syrk_cntl_finalize();
  FLASH_Syr2k_cntl_finalize();
  FLASH_Trmm_cntl_finalize();
  FLASH_Trsm_cntl_finalize();

  // LAPACK-level
  FLASH_Apply_pivots_cntl_finalize();
  FLASH_Chol_cntl_finalize();
  FLASH_LU_nopiv_cntl_finalize();
  FLASH_LU_piv_cntl_finalize();
  FLASH_LU_incpiv_cntl_finalize();
  FLASH_Trinv_cntl_finalize();
  FLASH_Ttmm_cntl_finalize();
  FLASH_Sylv_cntl_finalize();
  FLASH_QR2_UT_cntl_finalize();
  FLASH_CAQR2_UT_cntl_finalize();
  FLASH_Apply_Q_UT_cntl_finalize();
  FLASH_Apply_Q2_UT_cntl_finalize();
  FLASH_Apply_CAQ2_UT_cntl_finalize();
  FLASH_Apply_QUD_UT_cntl_finalize();
  FLASH_Eig_gest_cntl_finalize();

  // Compound LAPACK operations
  FLASH_Lyap_cntl_finalize();
  FLASH_SPDinv_cntl_finalize();
  FLASH_QR_UT_cntl_finalize();
  FLASH_QR_UT_inc_cntl_finalize();
  FLASH_LQ_UT_cntl_finalize();
  FLASH_CAQR_UT_inc_cntl_finalize();
  FLASH_Apply_Q_UT_inc_cntl_finalize();
  FLASH_Apply_CAQ_UT_inc_cntl_finalize();
  FLASH_UDdate_UT_cntl_finalize();
  FLASH_UDdate_UT_inc_cntl_finalize();
  FLASH_Apply_QUD_UT_inc_cntl_finalize();
}

