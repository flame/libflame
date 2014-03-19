/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Level-3 BLAS-like prototypes --------------------------------------------

// --- gemm ---

void bl1_sgemm( trans1_t transa, trans1_t transb, int m, int k, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dgemm( trans1_t transa, trans1_t transb, int m, int k, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_cgemm( trans1_t transa, trans1_t transb, int m, int k, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_zgemm( trans1_t transa, trans1_t transb, int m, int k, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_sgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, float*    alpha, float*    a, int lda, float*    b, int ldb, float*    beta, float*    c, int ldc );
void bl1_dgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, double*   alpha, double*   a, int lda, double*   b, int ldb, double*   beta, double*   c, int ldc );
void bl1_cgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bl1_zgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- hemm ---

void bl1_shemm( side1_t side, uplo1_t uplo, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dhemm( side1_t side, uplo1_t uplo, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_chemm( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_zhemm( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_chemm_blas( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bl1_zhemm_blas( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- herk ---

void bl1_sherk( uplo1_t uplo, trans1_t trans, int m, int k, float*  alpha, float*    a, int a_rs, int a_cs, float*  beta, float*    c, int c_rs, int c_cs );
void bl1_dherk( uplo1_t uplo, trans1_t trans, int m, int k, double* alpha, double*   a, int a_rs, int a_cs, double* beta, double*   c, int c_rs, int c_cs );
void bl1_cherk( uplo1_t uplo, trans1_t trans, int m, int k, float*  alpha, scomplex* a, int a_rs, int a_cs, float*  beta, scomplex* c, int c_rs, int c_cs );
void bl1_zherk( uplo1_t uplo, trans1_t trans, int m, int k, double* alpha, dcomplex* a, int a_rs, int a_cs, double* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_cherk_blas( uplo1_t uplo, trans1_t trans, int m, int k, float*  alpha, scomplex* a, int lda, float*  beta, scomplex* c, int ldc );
void bl1_zherk_blas( uplo1_t uplo, trans1_t trans, int m, int k, double* alpha, dcomplex* a, int lda, double* beta, dcomplex* c, int ldc );

// --- her2k ---

void bl1_sher2k( uplo1_t uplo, trans1_t trans, int m, int k, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*  beta, float*    c, int c_rs, int c_cs );
void bl1_dher2k( uplo1_t uplo, trans1_t trans, int m, int k, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double* beta, double*   c, int c_rs, int c_cs );
void bl1_cher2k( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, float*  beta, scomplex* c, int c_rs, int c_cs );
void bl1_zher2k( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, double* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_cher2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, float*  beta, scomplex* c, int ldc );
void bl1_zher2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, double* beta, dcomplex* c, int ldc );

// --- symm ---

void bl1_ssymm( side1_t side, uplo1_t uplo, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dsymm( side1_t side, uplo1_t uplo, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_csymm( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_zsymm( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_ssymm_blas( side1_t side, uplo1_t uplo, int m, int n, float*    alpha, float*    a, int lda, float*    b, int ldb, float*    beta, float*    c, int ldc );
void bl1_dsymm_blas( side1_t side, uplo1_t uplo, int m, int n, double*   alpha, double*   a, int lda, double*   b, int ldb, double*   beta, double*   c, int ldc );
void bl1_csymm_blas( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bl1_zsymm_blas( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- syrk ---

void bl1_ssyrk( uplo1_t uplo, trans1_t trans, int m, int k, float*    alpha, float*    a, int a_rs, int a_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dsyrk( uplo1_t uplo, trans1_t trans, int m, int k, double*   alpha, double*   a, int a_rs, int a_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_csyrk( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_zsyrk( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_ssyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, float*    alpha, float*    a, int lda, float*    beta, float*    c, int ldc );
void bl1_dsyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, double*   alpha, double*   a, int lda, double*   beta, double*   c, int ldc );
void bl1_csyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int lda, scomplex* beta, scomplex* c, int ldc );
void bl1_zsyrk_blas( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* beta, dcomplex* c, int ldc );

// --- syr2k ---

void bl1_ssyr2k( uplo1_t uplo, trans1_t trans, int m, int k, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dsyr2k( uplo1_t uplo, trans1_t trans, int m, int k, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_csyr2k( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_zsyr2k( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

void bl1_ssyr2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, float*    alpha, float*    a, int lda, float*    b, int ldb, float*    beta, float*    c, int ldc );
void bl1_dsyr2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, double*   alpha, double*   a, int lda, double*   b, int ldb, double*   beta, double*   c, int ldc );
void bl1_csyr2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc );
void bl1_zsyr2k_blas( uplo1_t uplo, trans1_t trans, int m, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc );

// --- trmm ---

void bl1_strmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dtrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_ctrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_ztrmm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bl1_strmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int lda, float*    b, int ldb );
void bl1_dtrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int lda, double*   b, int ldb );
void bl1_ctrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb );
void bl1_ztrmm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb );

// --- trsm ---

void bl1_strsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dtrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_ctrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_ztrsm( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bl1_strsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int lda, float*    b, int ldb );
void bl1_dtrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int lda, double*   b, int ldb );
void bl1_ctrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb );
void bl1_ztrsm_blas( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb );

// --- trmmsx ---

void bl1_strmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dtrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_ctrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_ztrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

// --- trsmsx ---

void bl1_strsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs, float*    beta, float*    c, int c_rs, int c_cs );
void bl1_dtrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs, double*   beta, double*   c, int c_rs, int c_cs );
void bl1_ctrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs );
void bl1_ztrsmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs );

