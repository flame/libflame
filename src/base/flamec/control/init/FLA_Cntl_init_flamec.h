/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

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

