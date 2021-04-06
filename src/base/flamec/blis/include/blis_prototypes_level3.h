/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Level-3 BLAS-like prototypes --------------------------------------------

// --- gemm ---

void bl1_sgemm( trans1_t transa, trans1_t transb, integer m, integer k, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs, float*    beta, float*    c, integer c_rs, integer c_cs );
void bl1_dgemm( trans1_t transa, trans1_t transb, integer m, integer k, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs, double*   beta, double*   c, integer c_rs, integer c_cs );
void bl1_cgemm( trans1_t transa, trans1_t transb, integer m, integer k, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs, scomplex* beta, scomplex* c, integer c_rs, integer c_cs );
void bl1_zgemm( trans1_t transa, trans1_t transb, integer m, integer k, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs, dcomplex* beta, dcomplex* c, integer c_rs, integer c_cs );

void bl1_sgemm_blas( trans1_t transa, trans1_t transb, integer m, integer n, integer k, float*    alpha, float*    a, integer lda, float*    b, integer ldb, float*    beta, float*    c, integer ldc );
void bl1_dgemm_blas( trans1_t transa, trans1_t transb, integer m, integer n, integer k, double*   alpha, double*   a, integer lda, double*   b, integer ldb, double*   beta, double*   c, integer ldc );
void bl1_cgemm_blas( trans1_t transa, trans1_t transb, integer m, integer n, integer k, scomplex* alpha, scomplex* a, integer lda, scomplex* b, integer ldb, scomplex* beta, scomplex* c, integer ldc );
void bl1_zgemm_blas( trans1_t transa, trans1_t transb, integer m, integer n, integer k, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* b, integer ldb, dcomplex* beta, dcomplex* c, integer ldc );

// --- hemm ---

void bl1_shemm( side1_t side, uplo1_t uplo, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs, float*    beta, float*    c, integer c_rs, integer c_cs );
void bl1_dhemm( side1_t side, uplo1_t uplo, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs, double*   beta, double*   c, integer c_rs, integer c_cs );
void bl1_chemm( side1_t side, uplo1_t uplo, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs, scomplex* beta, scomplex* c, integer c_rs, integer c_cs );
void bl1_zhemm( side1_t side, uplo1_t uplo, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs, dcomplex* beta, dcomplex* c, integer c_rs, integer c_cs );

void bl1_chemm_blas( side1_t side, uplo1_t uplo, integer m, integer n, scomplex* alpha, scomplex* a, integer lda, scomplex* b, integer ldb, scomplex* beta, scomplex* c, integer ldc );
void bl1_zhemm_blas( side1_t side, uplo1_t uplo, integer m, integer n, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* b, integer ldb, dcomplex* beta, dcomplex* c, integer ldc );

// --- herk ---

void bl1_sherk( uplo1_t uplo, trans1_t trans, integer m, integer k, float*  alpha, float*    a, integer a_rs, integer a_cs, float*  beta, float*    c, integer c_rs, integer c_cs );
void bl1_dherk( uplo1_t uplo, trans1_t trans, integer m, integer k, double* alpha, double*   a, integer a_rs, integer a_cs, double* beta, double*   c, integer c_rs, integer c_cs );
void bl1_cherk( uplo1_t uplo, trans1_t trans, integer m, integer k, float*  alpha, scomplex* a, integer a_rs, integer a_cs, float*  beta, scomplex* c, integer c_rs, integer c_cs );
void bl1_zherk( uplo1_t uplo, trans1_t trans, integer m, integer k, double* alpha, dcomplex* a, integer a_rs, integer a_cs, double* beta, dcomplex* c, integer c_rs, integer c_cs );

void bl1_cherk_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, float*  alpha, scomplex* a, integer lda, float*  beta, scomplex* c, integer ldc );
void bl1_zherk_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, double* alpha, dcomplex* a, integer lda, double* beta, dcomplex* c, integer ldc );

// --- her2k ---

void bl1_sher2k( uplo1_t uplo, trans1_t trans, integer m, integer k, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs, float*  beta, float*    c, integer c_rs, integer c_cs );
void bl1_dher2k( uplo1_t uplo, trans1_t trans, integer m, integer k, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs, double* beta, double*   c, integer c_rs, integer c_cs );
void bl1_cher2k( uplo1_t uplo, trans1_t trans, integer m, integer k, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs, float*  beta, scomplex* c, integer c_rs, integer c_cs );
void bl1_zher2k( uplo1_t uplo, trans1_t trans, integer m, integer k, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs, double* beta, dcomplex* c, integer c_rs, integer c_cs );

void bl1_cher2k_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, scomplex* alpha, scomplex* a, integer lda, scomplex* b, integer ldb, float*  beta, scomplex* c, integer ldc );
void bl1_zher2k_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* b, integer ldb, double* beta, dcomplex* c, integer ldc );

// --- symm ---

void bl1_ssymm( side1_t side, uplo1_t uplo, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs, float*    beta, float*    c, integer c_rs, integer c_cs );
void bl1_dsymm( side1_t side, uplo1_t uplo, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs, double*   beta, double*   c, integer c_rs, integer c_cs );
void bl1_csymm( side1_t side, uplo1_t uplo, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs, scomplex* beta, scomplex* c, integer c_rs, integer c_cs );
void bl1_zsymm( side1_t side, uplo1_t uplo, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs, dcomplex* beta, dcomplex* c, integer c_rs, integer c_cs );

void bl1_ssymm_blas( side1_t side, uplo1_t uplo, integer m, integer n, float*    alpha, float*    a, integer lda, float*    b, integer ldb, float*    beta, float*    c, integer ldc );
void bl1_dsymm_blas( side1_t side, uplo1_t uplo, integer m, integer n, double*   alpha, double*   a, integer lda, double*   b, integer ldb, double*   beta, double*   c, integer ldc );
void bl1_csymm_blas( side1_t side, uplo1_t uplo, integer m, integer n, scomplex* alpha, scomplex* a, integer lda, scomplex* b, integer ldb, scomplex* beta, scomplex* c, integer ldc );
void bl1_zsymm_blas( side1_t side, uplo1_t uplo, integer m, integer n, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* b, integer ldb, dcomplex* beta, dcomplex* c, integer ldc );

// --- syrk ---

void bl1_ssyrk( uplo1_t uplo, trans1_t trans, integer m, integer k, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    beta, float*    c, integer c_rs, integer c_cs );
void bl1_dsyrk( uplo1_t uplo, trans1_t trans, integer m, integer k, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   beta, double*   c, integer c_rs, integer c_cs );
void bl1_csyrk( uplo1_t uplo, trans1_t trans, integer m, integer k, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* beta, scomplex* c, integer c_rs, integer c_cs );
void bl1_zsyrk( uplo1_t uplo, trans1_t trans, integer m, integer k, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* beta, dcomplex* c, integer c_rs, integer c_cs );

void bl1_ssyrk_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, float*    alpha, float*    a, integer lda, float*    beta, float*    c, integer ldc );
void bl1_dsyrk_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, double*   alpha, double*   a, integer lda, double*   beta, double*   c, integer ldc );
void bl1_csyrk_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, scomplex* alpha, scomplex* a, integer lda, scomplex* beta, scomplex* c, integer ldc );
void bl1_zsyrk_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* beta, dcomplex* c, integer ldc );

// --- syr2k ---

void bl1_ssyr2k( uplo1_t uplo, trans1_t trans, integer m, integer k, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs, float*    beta, float*    c, integer c_rs, integer c_cs );
void bl1_dsyr2k( uplo1_t uplo, trans1_t trans, integer m, integer k, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs, double*   beta, double*   c, integer c_rs, integer c_cs );
void bl1_csyr2k( uplo1_t uplo, trans1_t trans, integer m, integer k, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs, scomplex* beta, scomplex* c, integer c_rs, integer c_cs );
void bl1_zsyr2k( uplo1_t uplo, trans1_t trans, integer m, integer k, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs, dcomplex* beta, dcomplex* c, integer c_rs, integer c_cs );

void bl1_ssyr2k_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, float*    alpha, float*    a, integer lda, float*    b, integer ldb, float*    beta, float*    c, integer ldc );
void bl1_dsyr2k_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, double*   alpha, double*   a, integer lda, double*   b, integer ldb, double*   beta, double*   c, integer ldc );
void bl1_csyr2k_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, scomplex* alpha, scomplex* a, integer lda, scomplex* b, integer ldb, scomplex* beta, scomplex* c, integer ldc );
void bl1_zsyr2k_blas( uplo1_t uplo, trans1_t trans, integer m, integer k, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* b, integer ldb, dcomplex* beta, dcomplex* c, integer ldc );

// --- trmm ---

void bl1_strmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_dtrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_ctrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_ztrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

void bl1_strmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, float*    alpha, float*    a, integer lda, float*    b, integer ldb );
void bl1_dtrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, double*   alpha, double*   a, integer lda, double*   b, integer ldb );
void bl1_ctrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, scomplex* alpha, scomplex* a, integer lda, scomplex* b, integer ldb );
void bl1_ztrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* b, integer ldb );

// --- trsm ---

void bl1_strsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_dtrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_ctrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_ztrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

void bl1_strsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, float*    alpha, float*    a, integer lda, float*    b, integer ldb );
void bl1_dtrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, double*   alpha, double*   a, integer lda, double*   b, integer ldb );
void bl1_ctrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, scomplex* alpha, scomplex* a, integer lda, scomplex* b, integer ldb );
void bl1_ztrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* b, integer ldb );

// --- trmmsx ---

void bl1_strmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs, float*    beta, float*    c, integer c_rs, integer c_cs );
void bl1_dtrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs, double*   beta, double*   c, integer c_rs, integer c_cs );
void bl1_ctrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs, scomplex* beta, scomplex* c, integer c_rs, integer c_cs );
void bl1_ztrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs, dcomplex* beta, dcomplex* c, integer c_rs, integer c_cs );

// --- trsmsx ---

void bl1_strsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs, float*    beta, float*    c, integer c_rs, integer c_cs );
void bl1_dtrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs, double*   beta, double*   c, integer c_rs, integer c_cs );
void bl1_ctrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs, scomplex* beta, scomplex* c, integer c_rs, integer c_cs );
void bl1_ztrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs, dcomplex* beta, dcomplex* c, integer c_rs, integer c_cs );

