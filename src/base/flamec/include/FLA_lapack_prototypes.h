/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- top-level wrapper prototypes --------------------------------------------

FLA_Error FLA_Chol( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_LU_nopiv( FLA_Obj A );
FLA_Error FLA_LU_piv( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_QR_UT( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_QR_UT_piv( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p );
FLA_Error FLA_LQ_UT( FLA_Obj A, FLA_Obj S );
FLA_Error FLA_Trinv( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLA_Ttmm( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Sylv( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_SPDinv( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Hess_UT( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Eig_gest( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B );

FLA_Error FLA_Accum_T_UT( FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj tau, FLA_Obj T );
FLA_Error FLA_Apply_H2_UT( FLA_Side side, FLA_Obj tau, FLA_Obj u2, FLA_Obj a1, FLA_Obj A2 );
FLA_Error FLA_Apply_HUD_UT( FLA_Side side, FLA_Obj tau, FLA_Obj w12t, FLA_Obj u2, FLA_Obj v2, FLA_Obj r12t, FLA_Obj C2, FLA_Obj D2 );
FLA_Error FLA_Apply_Q_UT( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B );
FLA_Error FLA_Apply_pivots( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A );

// --- task wrapper prototypes -------------------------------------------------

FLA_Error FLA_Chol_task( FLA_Uplo uplo, FLA_Obj A, fla_chol_t* cntl );
FLA_Error FLA_Chol_l_task( FLA_Obj A, fla_chol_t* cntl );
FLA_Error FLA_Chol_u_task( FLA_Obj A, fla_chol_t* cntl );
FLA_Error FLA_LU_piv_macro_task( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl );
FLA_Error FLA_Apply_pivots_task( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
FLA_Error FLA_Apply_pivots_ln_task( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
FLA_Error FLA_Apply_pivots_macro_task( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
FLA_Error FLA_LU_nopiv_task( FLA_Obj A, fla_lu_t* cntl );
FLA_Error FLA_LU_piv_task( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl );
FLA_Error FLA_LU_piv_copy_task( FLA_Obj A, FLA_Obj p, FLA_Obj U, fla_lu_t* cntl );
FLA_Error FLA_Trsm_piv_task( FLA_Obj A, FLA_Obj B, FLA_Obj p, fla_trsm_t* cntl );
FLA_Error FLA_SA_LU_task( FLA_Obj U, FLA_Obj D, FLA_Obj p, FLA_Obj L, dim_t nb_alg, fla_lu_t* cntl );
FLA_Error FLA_SA_FS_task( FLA_Obj L, FLA_Obj D, FLA_Obj p, FLA_Obj C, FLA_Obj E, dim_t nb_alg, fla_gemm_t* cntl );
FLA_Error FLA_Trinv_task( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_ln_task( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_lu_task( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_un_task( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Trinv_uu_task( FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Ttmm_task( FLA_Uplo uplo, FLA_Obj A, fla_ttmm_t* cntl );
FLA_Error FLA_Ttmm_l_task( FLA_Obj A, fla_ttmm_t* cntl );
FLA_Error FLA_Ttmm_u_task( FLA_Obj A, fla_ttmm_t* cntl );
FLA_Error FLA_Sylv_task( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nh_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hn_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Lyap_task( FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl );
FLA_Error FLA_Lyap_n_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl );
FLA_Error FLA_Lyap_h_task( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl );
FLA_Error FLA_Apply_Q_UT_task( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhbc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhbr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhfc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lhfr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lnbc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lnbr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lnfc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_lnfr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rhbc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rhbr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rhfc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rhfr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnbc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnbr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnfc_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q_UT_rnfr_task( FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q2_UT_task( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E, fla_apq2ut_t* cntl );
FLA_Error FLA_Apply_Q2_UT_lhfc_task( FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E, fla_apq2ut_t* cntl );
FLA_Error FLA_Apply_CAQ2_UT_task( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E, fla_apcaq2ut_t* cntl );
FLA_Error FLA_Apply_CAQ2_UT_lhfc_task( FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E, fla_apcaq2ut_t* cntl );
FLA_Error FLA_QR2_UT_task( FLA_Obj B, FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl );
FLA_Error FLA_CAQR2_UT_task( FLA_Obj B, FLA_Obj D, FLA_Obj T, fla_caqr2ut_t* cntl );
FLA_Error FLA_QR_UT_macro_task( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_task( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_copy_task( FLA_Obj A, FLA_Obj T, FLA_Obj U, fla_qrut_t* cntl );
FLA_Error FLA_LQ_UT_macro_task( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl );
FLA_Error FLA_LQ_UT_task( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl );
FLA_Error FLA_UDdate_UT_task( FLA_Obj R, FLA_Obj C, FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl );
FLA_Error FLA_Apply_QUD_UT_task( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj U, FLA_Obj C, FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl );
FLA_Error FLA_Apply_QUD_UT_lhfc_task( FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj U, FLA_Obj C, FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl );
FLA_Error FLA_Eig_gest_task( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_il_task( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_iu_task( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nl_task( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nu_task( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );

// --- external wrapper prototypes ---------------------------------------------

FLA_Error FLA_Apply_Q_blk_external( FLA_Side side, FLA_Trans trans, FLA_Store storev, FLA_Obj A, FLA_Obj t, FLA_Obj B );

FLA_Error FLA_Apply_pivots_unb_external( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A );
FLA_Error FLA_Apply_pivots_ln_unb_ext( FLA_Obj p, FLA_Obj A );

FLA_Error FLA_Apply_pivots_macro_external( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A );

FLA_Error FLA_Chol_blk_external( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Chol_l_blk_ext( FLA_Obj A );
FLA_Error FLA_Chol_u_blk_ext( FLA_Obj A );
FLA_Error FLA_Chol_unb_external( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Chol_l_unb_ext( FLA_Obj A );
FLA_Error FLA_Chol_u_unb_ext( FLA_Obj A );

FLA_Error FLA_LU_piv_blk_external( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_LU_piv_blk_ext( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_LU_piv_unb_external( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_LU_piv_unb_ext( FLA_Obj A, FLA_Obj p );

FLA_Error FLA_QR_blk_external( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_QR_unb_external( FLA_Obj A, FLA_Obj t );

FLA_Error FLA_LQ_blk_external( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_LQ_unb_external( FLA_Obj A, FLA_Obj t );

FLA_Error FLA_Hess_blk_external( FLA_Obj A, FLA_Obj t, int ilo, int ihi );
FLA_Error FLA_Hess_unb_external( FLA_Obj A, FLA_Obj t, int ilo, int ihi );

FLA_Error FLA_Tridiag_blk_external( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t );
FLA_Error FLA_Tridiag_unb_external( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t );

FLA_Error FLA_Bidiag_blk_external( FLA_Obj A, FLA_Obj tu, FLA_Obj tv );
FLA_Error FLA_Bidiag_unb_external( FLA_Obj A, FLA_Obj tu, FLA_Obj tv );

FLA_Error FLA_QR_form_Q_external( FLA_Obj A, FLA_Obj t );

FLA_Error FLA_Tridiag_form_Q_external( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t );
FLA_Error FLA_Tridiag_apply_Q_external( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Obj A, FLA_Obj t, FLA_Obj B );

FLA_Error FLA_Bidiag_form_U_external( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_Bidiag_form_V_external( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_Bidiag_apply_U_external( FLA_Side side, FLA_Trans trans, FLA_Obj A, FLA_Obj t, FLA_Obj B );
FLA_Error FLA_Bidiag_apply_V_external( FLA_Side side, FLA_Trans trans, FLA_Obj A, FLA_Obj t, FLA_Obj B );

FLA_Error FLA_Trinv_blk_external( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLA_Trinv_ln_blk_ext( FLA_Obj A );
FLA_Error FLA_Trinv_lu_blk_ext( FLA_Obj A );
FLA_Error FLA_Trinv_un_blk_ext( FLA_Obj A );
FLA_Error FLA_Trinv_uu_blk_ext( FLA_Obj A );
FLA_Error FLA_Trinv_unb_external( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLA_Trinv_ln_unb_ext( FLA_Obj A );
FLA_Error FLA_Trinv_lu_unb_ext( FLA_Obj A );
FLA_Error FLA_Trinv_un_unb_ext( FLA_Obj A );
FLA_Error FLA_Trinv_uu_unb_ext( FLA_Obj A );

FLA_Error FLA_Ttmm_blk_external( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Ttmm_l_blk_ext( FLA_Obj A );
FLA_Error FLA_Ttmm_u_blk_ext( FLA_Obj A );
FLA_Error FLA_Ttmm_unb_external( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Ttmm_l_unb_ext( FLA_Obj A );
FLA_Error FLA_Ttmm_u_unb_ext( FLA_Obj A );

FLA_Error FLA_Sylv_blk_external( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_blk_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nh_blk_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_blk_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_blk_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_unb_external( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_unb_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nh_unb_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hn_unb_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_unb_ext( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );

FLA_Error FLA_SPDinv_blk_external( FLA_Uplo uplo, FLA_Obj A );

FLA_Error FLA_Eig_gest_blk_external( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Eig_gest_il_blk_ext( FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Eig_gest_iu_blk_ext( FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_blk_ext( FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_blk_ext( FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Eig_gest_unb_external( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Eig_gest_il_unb_ext( FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Eig_gest_iu_unb_ext( FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Eig_gest_nl_unb_ext( FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_unb_ext( FLA_Obj A, FLA_Obj B );

FLA_Error FLA_Tevd_external( FLA_Evd_type jobz, FLA_Obj d, FLA_Obj e, FLA_Obj A );
FLA_Error FLA_Tevdd_external( FLA_Evd_type jobz, FLA_Obj d, FLA_Obj e, FLA_Obj A );
FLA_Error FLA_Tevdr_external( FLA_Evd_type jobz, FLA_Obj d, FLA_Obj e, FLA_Obj l, FLA_Obj A );
FLA_Error FLA_Hevd_external( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l );
FLA_Error FLA_Hevdd_external( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l );
FLA_Error FLA_Hevdr_external( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l, FLA_Obj Z );
FLA_Error FLA_Bsvd_external( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, FLA_Obj U, FLA_Obj V );
FLA_Error FLA_Bsvdd_external( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, FLA_Obj U, FLA_Obj V );
FLA_Error FLA_Svd_external( FLA_Svd_type jobu, FLA_Svd_type jobv, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V );
FLA_Error FLA_Svdd_external( FLA_Svd_type jobz, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V );

// --- external HIP prototypes -------------------------------------------------
#ifdef FLA_ENABLE_HIP
FLA_Error FLA_Apply_pivots_unb_external_hip( rocblas_handle handle, FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, void* A_hip );
FLA_Error FLA_Apply_pivots_ln_unb_ext_hip( rocblas_handle handle, FLA_Obj p, FLA_Obj A, void* A_hip );
FLA_Error FLA_Apply_Q_blk_external_hip( rocblas_handle handle, FLA_Side side, FLA_Trans trans, FLA_Store storev, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Bidiag_apply_U_external_hip( rocblas_handle handle, FLA_Side side, FLA_Trans trans, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Bidiag_apply_V_external_hip( rocblas_handle handle, FLA_Side side, FLA_Trans trans, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Bidiag_blk_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj tu, void* tu_hip, FLA_Obj tv, void* tv_hip );
FLA_Error FLA_Bidiag_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj tu, void* tu_hip, FLA_Obj tv, void* tv_hip );
FLA_Error FLA_Bidiag_unb_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj tu, void* tu_hip, FLA_Obj tv, void* tv_hip );
FLA_Error FLA_Bidiag_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj tu, void* tu_hip, FLA_Obj tv, void* tv_hip );
FLA_Error FLA_Bidiag_form_U_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_Bidiag_form_V_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_Bsvd_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj d, void* d_hip, FLA_Obj e, void* e_hip, FLA_Obj U, void* U_hip, FLA_Obj V, void* V_hip );
FLA_Error FLA_Chol_blk_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip );
FLA_Error FLA_Chol_l_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip );
FLA_Error FLA_Chol_u_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip );
FLA_Error FLA_Chol_unb_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip );
FLA_Error FLA_Chol_l_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip );
FLA_Error FLA_Chol_u_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip );
FLA_Error FLA_Eig_gest_blk_external_hip( rocblas_handle handle, FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Eig_gest_il_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Eig_gest_iu_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Eig_gest_nl_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Eig_gest_nu_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Eig_gest_unb_external_hip( rocblas_handle handle, FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Eig_gest_il_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Eig_gest_iu_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Eig_gest_nl_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Eig_gest_nu_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Hevdd_external_hip( rocblas_handle handle, FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj e, void* e_hip );
FLA_Error FLA_Hevd_external_hip( rocblas_handle handle, FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj e, void* e_hip );
FLA_Error FLA_LU_piv_blk_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj p );
FLA_Error FLA_LU_piv_copy_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj p, FLA_Obj U, void* U_hip );
FLA_Error FLA_LQ_blk_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_LQ_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_LQ_unb_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_LQ_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_QR_form_Q_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_QR_unb_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_QR_unb_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_SA_Apply_pivots_hip( rocblas_handle handle, FLA_Obj C, void* C_hip, FLA_Obj E, void* E_hip, FLA_Obj p );
FLA_Error FLA_SA_FS_blk_hip( rocblas_handle handle, FLA_Obj L, FLA_Obj D, void* D_hip, FLA_Obj p, FLA_Obj C, void* C_hip,FLA_Obj E, void* E_hip, dim_t nb_alg );
FLA_Error FLA_Svd_external_hip( rocblas_handle handle, FLA_Svd_type jobu, FLA_Svd_type jobv, FLA_Obj A, void* A_hip, FLA_Obj s, void* s_hip, FLA_Obj U, void* U_hip, FLA_Obj V, void* V_hip );
FLA_Error FLA_Tevdd_external_hip( rocblas_handle handle, FLA_Evd_type jobz, FLA_Obj d, void* d_hip, FLA_Obj e, void* e_hip, FLA_Obj A, void* A_hip );
FLA_Error FLA_Tevd_external_hip( rocblas_handle handle, FLA_Evd_type jobz, FLA_Obj d, void* d_hip, FLA_Obj e, void* e_hip, FLA_Obj A, void* A_hip );
FLA_Error FLA_Tridiag_apply_Q_external_hip( rocblas_handle handle, FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip, FLA_Obj B, void* B_hip );
FLA_Error FLA_Tridiag_blk_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* a_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_Tridiag_blk_ext_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* a_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_Tridiag_form_Q_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_Tridiag_unb_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_Tridiag_unb_ext_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Obj A, void* A_hip, FLA_Obj t, void* t_hip );
FLA_Error FLA_Trinv_blk_external_hip( rocblas_handle handle, FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A, void* A_hip );
FLA_Error FLA_Trinv_ln_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip );
FLA_Error FLA_Trinv_lu_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip );
FLA_Error FLA_Trinv_un_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip );
FLA_Error FLA_Trinv_uu_blk_ext_hip( rocblas_handle handle, FLA_Obj A, void* A_hip );
FLA_Error FLA_Trsm_piv_external_hip( rocblas_handle handle, FLA_Obj A, void* A_hip, FLA_Obj B, void* B_hip, FLA_Obj p );
#endif

// --- check routine prototypes ------------------------------------------------

FLA_Error FLA_Chol_check( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Chol_solve_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj B, FLA_Obj X );
FLA_Error FLA_LU_nopiv_check( FLA_Obj A );
FLA_Error FLA_LU_nopiv_solve_check( FLA_Obj A, FLA_Obj B, FLA_Obj X );
FLA_Error FLA_LU_piv_check( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_LU_piv_solve_check( FLA_Obj A, FLA_Obj p, FLA_Obj B, FLA_Obj X );
FLA_Error FLA_LU_incpiv_check( FLA_Obj A, FLA_Obj p, FLA_Obj L );
FLA_Error FLA_LU_incpiv_solve_check( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj B, FLA_Obj X );
FLA_Error FLA_FS_incpiv_check( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj b );
FLA_Error FLA_QR_check( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_QR_UT_check( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_QR_UT_solve_check( FLA_Obj A, FLA_Obj T, FLA_Obj B, FLA_Obj X );
FLA_Error FLA_QR_UT_recover_tau_check( FLA_Obj T, FLA_Obj tau );
FLA_Error FLA_QR_UT_form_Q_check( FLA_Obj A, FLA_Obj T, FLA_Obj Q );
FLA_Error FLA_LQ_check( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_LQ_UT_check( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_LQ_UT_solve_check( FLA_Obj A, FLA_Obj T, FLA_Obj B, FLA_Obj X );
FLA_Error FLA_LQ_UT_recover_tau_check( FLA_Obj T, FLA_Obj tau );
FLA_Error FLA_LQ_UT_form_Q_check( FLA_Obj A, FLA_Obj T, FLA_Obj Q );
FLA_Error FLA_Hess_check( FLA_Obj A, FLA_Obj t, int ilo, int ihi );
FLA_Error FLA_Hess_UT_check( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_recover_tau_check( FLA_Obj T, FLA_Obj tau );
FLA_Error FLA_Tridiag_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t );
FLA_Error FLA_Tridiag_UT_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_recover_tau_check( FLA_Obj T, FLA_Obj tau );
FLA_Error FLA_Tridiag_UT_scale_diagonals_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Tridiag_UT_extract_diagonals_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Tridiag_UT_extract_real_diagonals_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Tridiag_UT_realify_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_realify_subdiagonal_check( FLA_Obj b, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_shift_U_check( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Tridiag_UT_form_Q_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T, FLA_Obj Q );
FLA_Error FLA_Trinv_check( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLA_Bidiag_check( FLA_Obj A, FLA_Obj tu, FLA_Obj tv );
FLA_Error FLA_Bidiag_UT_check( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_recover_tau_check( FLA_Obj TU, FLA_Obj TV, FLA_Obj tu, FLA_Obj tv );
FLA_Error FLA_Bidiag_UT_extract_diagonals_check( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_extract_real_diagonals_check( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_scale_diagonals_check( FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Bidiag_UT_realify_check( FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_realify_diagonals_check( FLA_Uplo uplo, FLA_Obj a, FLA_Obj b, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Bidiag_UT_form_U_check( FLA_Obj A, FLA_Obj T, FLA_Obj U );
FLA_Error FLA_Bidiag_UT_form_V_check( FLA_Obj A, FLA_Obj S, FLA_Obj V );
FLA_Error FLA_Ttmm_check( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Sylv_check( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Lyap_check( FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_SPDinv_check( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Eig_gest_check( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj B );

FLA_Error FLA_Apply_Q_check( FLA_Side side, FLA_Trans trans, FLA_Store storev, FLA_Obj A, FLA_Obj t, FLA_Obj B );

FLA_Error FLA_QR_form_Q_check( FLA_Obj A, FLA_Obj t );

FLA_Error FLA_Tridiag_form_Q_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj t );
FLA_Error FLA_Tridiag_apply_Q_check( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Obj A, FLA_Obj t, FLA_Obj B );

FLA_Error FLA_Bidiag_form_U_check( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_Bidiag_form_V_check( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_Bidiag_apply_U_check( FLA_Side side, FLA_Trans trans, FLA_Obj A, FLA_Obj t, FLA_Obj B );
FLA_Error FLA_Bidiag_apply_V_check( FLA_Side side, FLA_Trans trans, FLA_Obj A, FLA_Obj t, FLA_Obj B );

FLA_Error FLA_Apply_Q_UT_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B );
FLA_Error FLA_Apply_Q2_UT_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E );
FLA_Error FLA_Apply_QUD_UT_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj U, FLA_Obj C, FLA_Obj V, FLA_Obj D );
FLA_Error FLA_Apply_pivots_check( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A );
FLA_Error FLA_QR2_UT_check( FLA_Obj B, FLA_Obj D, FLA_Obj T );
FLA_Error FLA_CAQR2_UT_check( FLA_Obj B, FLA_Obj D, FLA_Obj T );
FLA_Error FLA_QR_UT_inc_check( FLA_Obj A, FLA_Obj TW );
FLA_Error FLA_Apply_Q_UT_inc_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj TW, FLA_Obj W1, FLA_Obj B );
FLA_Error FLA_Apply_CAQ_UT_inc_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW, FLA_Obj W1, FLA_Obj B );

FLA_Error FLA_QR_UT_inc_solve_check( FLA_Obj A, FLA_Obj TW, FLA_Obj B, FLA_Obj X );
FLA_Error FLA_CAQR_UT_inc_solve_check( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW, FLA_Obj B, FLA_Obj X );

FLA_Error FLA_UDdate_UT_check( FLA_Obj R, FLA_Obj C, FLA_Obj D, FLA_Obj T );
FLA_Error FLA_UDdate_UT_update_rhs_check( FLA_Obj T, FLA_Obj bR, FLA_Obj C, FLA_Obj bC, FLA_Obj D, FLA_Obj bD );
FLA_Error FLA_UDdate_UT_solve_check( FLA_Obj R, FLA_Obj bR, FLA_Obj x );

FLA_Error FLA_UDdate_UT_inc_check( FLA_Obj R, FLA_Obj C, FLA_Obj D, FLA_Obj T, FLA_Obj W );
FLA_Error FLA_UDdate_UT_inc_update_rhs_check( FLA_Obj T, FLA_Obj bR, FLA_Obj C, FLA_Obj bC, FLA_Obj D, FLA_Obj bD );
FLA_Error FLA_UDdate_UT_inc_solve_check( FLA_Obj R, FLA_Obj bR, FLA_Obj x );

FLA_Error FLA_CAQR_UT_inc_check( dim_t p, FLA_Obj A, FLA_Obj ATW, FLA_Obj R, FLA_Obj RTW );

FLA_Error FLA_Apply_QUD_UT_inc_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj U, FLA_Obj C, FLA_Obj V, FLA_Obj D );

FLA_Error FLA_Apply_H2_UT_check( FLA_Side side, FLA_Obj tau, FLA_Obj u2, FLA_Obj a1t, FLA_Obj A2 );
FLA_Error FLA_Apply_HUD_UT_check( FLA_Side side, FLA_Obj tau, FLA_Obj w12t, FLA_Obj u2, FLA_Obj v2, FLA_Obj r12t, FLA_Obj C2, FLA_Obj D2 );
FLA_Error FLA_Accum_T_UT_check( FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj tau, FLA_Obj T );

FLA_Error FLA_Tevd_compute_scaling_check( FLA_Obj d, FLA_Obj e, FLA_Obj sigma );
FLA_Error FLA_Hevd_compute_scaling_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj sigma );
FLA_Error FLA_Hevd_check( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l );
FLA_Error FLA_Hevdd_check( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l );
FLA_Error FLA_Hevdr_check( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l, FLA_Obj Z );

FLA_Error FLA_Bsvd_check( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e,
                          FLA_Obj G, FLA_Obj H,
                          FLA_Svd_type jobu, FLA_Obj U,
                          FLA_Svd_type jobv, FLA_Obj V );
FLA_Error FLA_Bsvd_ext_check( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e,
                              FLA_Obj G, FLA_Obj H,
                              FLA_Svd_type jobu, FLA_Obj U,
                              FLA_Svd_type jobv, FLA_Obj V,
                              FLA_Bool apply_Uh2C, FLA_Obj C );
FLA_Error FLA_Bsvd_compute_scaling_check( FLA_Obj d, FLA_Obj e, FLA_Obj sigma );
FLA_Error FLA_Svd_compute_scaling_check( FLA_Obj A, FLA_Obj sigma );
FLA_Error FLA_Svd_check( FLA_Svd_type jobu, FLA_Svd_type jobv, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V );
FLA_Error FLA_Svd_ext_check( FLA_Svd_type jobu, FLA_Trans transu, FLA_Svd_type jobv, FLA_Trans transv,
                             FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V );
FLA_Error FLA_Svdd_check( FLA_Svd_type jobz, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V );

FLA_Error FLA_Chol_internal_check( FLA_Uplo uplo, FLA_Obj A, fla_chol_t* cntl );
FLA_Error FLA_LU_nopiv_internal_check( FLA_Obj A, fla_lu_t* cntl );
FLA_Error FLA_Trinv_internal_check( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A, fla_trinv_t* cntl );
FLA_Error FLA_Ttmm_internal_check( FLA_Uplo uplo, FLA_Obj A, fla_ttmm_t* cntl );
FLA_Error FLA_SPDinv_internal_check( FLA_Uplo uplo, FLA_Obj A, fla_spdinv_t* cntl );
FLA_Error FLA_Sylv_internal_check( FLA_Trans transa, FLA_Trans transb, FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Lyap_internal_check( FLA_Trans trans, FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl );
FLA_Error FLA_QR_UT_internal_check( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_copy_internal_check( FLA_Obj A, FLA_Obj T, FLA_Obj U, fla_qrut_t* cntl );
FLA_Error FLA_QR2_UT_internal_check( FLA_Obj B, FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl );
FLA_Error FLA_CAQR2_UT_internal_check( FLA_Obj B, FLA_Obj D, FLA_Obj T, fla_caqr2ut_t* cntl );
FLA_Error FLA_LQ_UT_internal_check( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl );
FLA_Error FLA_Hess_UT_internal_check( FLA_Obj A, FLA_Obj T, fla_hessut_t* cntl );
FLA_Error FLA_Tridiag_UT_internal_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T, fla_tridiagut_t* cntl );
FLA_Error FLA_Bidiag_UT_internal_check( FLA_Obj A, FLA_Obj TU, FLA_Obj TV, fla_bidiagut_t* cntl );

FLA_Error FLA_UDdate_UT_internal_check( FLA_Obj R, FLA_Obj C, FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl );

FLA_Error FLA_Apply_Q_UT_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj T, FLA_Obj W, FLA_Obj B, fla_apqut_t* cntl );
FLA_Error FLA_Apply_Q2_UT_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E, fla_apq2ut_t* cntl );
FLA_Error FLA_Apply_CAQ2_UT_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj D, FLA_Obj T, FLA_Obj W, FLA_Obj C, FLA_Obj E, fla_apcaq2ut_t* cntl );
FLA_Error FLA_Apply_QUD_UT_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj U, FLA_Obj C, FLA_Obj V, FLA_Obj D, fla_apqudut_t* cntl );

FLA_Error FLA_Apply_Q_UT_inc_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj A, FLA_Obj TW, FLA_Obj W1, FLA_Obj B, fla_apqutinc_t* cntl );
FLA_Error FLA_Apply_CAQ_UT_inc_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj R, FLA_Obj TW, FLA_Obj W, FLA_Obj B, fla_apcaqutinc_t* cntl );
FLA_Error FLA_Apply_QUD_UT_inc_internal_check( FLA_Side side, FLA_Trans trans, FLA_Direct direct, FLA_Store storev, FLA_Obj T, FLA_Obj W, FLA_Obj R, FLA_Obj U, FLA_Obj C, FLA_Obj V, FLA_Obj D, fla_apqudutinc_t* cntl );

FLA_Error FLA_Eig_gest_internal_check( FLA_Inv inv, FLA_Uplo uplo, FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
