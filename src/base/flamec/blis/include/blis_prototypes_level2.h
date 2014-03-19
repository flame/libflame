/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Level-2 BLAS-like prototypes --------------------------------------------

// --- gemv ---

void bl1_sgemv( trans1_t transa, conj1_t conjx, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_dgemv( trans1_t transa, conj1_t conjx, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_cgemv( trans1_t transa, conj1_t conjx, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zgemv( trans1_t transa, conj1_t conjx, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

void bl1_sgemv_blas( trans1_t transa, int m, int n, float*    alpha, float*    a, int lda, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_dgemv_blas( trans1_t transa, int m, int n, double*   alpha, double*   a, int lda, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_cgemv_blas( trans1_t transa, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zgemv_blas( trans1_t transa, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- ger ---

void bl1_sger( conj1_t conjx, conj1_t conjy, int m, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int a_rs, int a_cs );
void bl1_dger( conj1_t conjx, conj1_t conjy, int m, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int a_rs, int a_cs );
void bl1_cger( conj1_t conjx, conj1_t conjy, int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs );
void bl1_zger( conj1_t conjx, conj1_t conjy, int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs );

void bl1_sger_blas(  int m, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int lda );
void bl1_dger_blas(  int m, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int lda );
void bl1_cgerc_blas( int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bl1_cgeru_blas( int m, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bl1_zgerc_blas( int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );
void bl1_zgeru_blas( int m, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );

// --- hemv ---

void bl1_shemv( uplo1_t uplo, conj1_t conj, int m, float*    alpha, float*    a, int a_rs, int a_cs, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_dhemv( uplo1_t uplo, conj1_t conj, int m, double*   alpha, double*   a, int a_rs, int a_cs, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_chemv( uplo1_t uplo, conj1_t conj, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zhemv( uplo1_t uplo, conj1_t conj, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

void bl1_chemv_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zhemv_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- her ---

void bl1_sher( uplo1_t uplo, conj1_t conj, int m, float*  alpha, float*    x, int incx, float*    a, int a_rs, int a_cs );
void bl1_dher( uplo1_t uplo, conj1_t conj, int m, double* alpha, double*   x, int incx, double*   a, int a_rs, int a_cs );
void bl1_cher( uplo1_t uplo, conj1_t conj, int m, float*  alpha, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs );
void bl1_zher( uplo1_t uplo, conj1_t conj, int m, double* alpha, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs );

void bl1_cher_blas( uplo1_t uplo, int m, float*  alpha, scomplex* x, int incx, scomplex* a, int lda );
void bl1_zher_blas( uplo1_t uplo, int m, double* alpha, dcomplex* x, int incx, dcomplex* a, int lda );

// --- her2 ---

void bl1_sher2( uplo1_t uplo, conj1_t conj, int m, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int a_rs, int a_cs );
void bl1_dher2( uplo1_t uplo, conj1_t conj, int m, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int a_rs, int a_cs );
void bl1_cher2( uplo1_t uplo, conj1_t conj, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs );
void bl1_zher2( uplo1_t uplo, conj1_t conj, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs );

void bl1_cher2_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bl1_zher2_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );

// --- symv ---

void bl1_ssymv( uplo1_t uplo, int m, float*    alpha, float*    a, int a_rs, int a_cs, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_dsymv( uplo1_t uplo, int m, double*   alpha, double*   a, int a_rs, int a_cs, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_csymv( uplo1_t uplo, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zsymv( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

void bl1_ssymv_blas( uplo1_t uplo, int m, float*    alpha, float*    a, int lda, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_dsymv_blas( uplo1_t uplo, int m, double*   alpha, double*   a, int lda, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_csymv_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* a, int lda, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zsymv_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* a, int lda, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- syr ---

void bl1_ssyr( uplo1_t uplo, int m, float*    alpha, float*    x, int incx, float*    a, int a_rs, int a_cs );
void bl1_dsyr( uplo1_t uplo, int m, double*   alpha, double*   x, int incx, double*   a, int a_rs, int a_cs );
void bl1_csyr( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* a, int a_rs, int a_cs );
void bl1_zsyr( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* a, int a_rs, int a_cs );

void bl1_ssyr_blas( uplo1_t uplo, int m, float*    alpha, float*    x, int incx, float*    a, int lda );
void bl1_dsyr_blas( uplo1_t uplo, int m, double*   alpha, double*   x, int incx, double*   a, int lda );
void bl1_csyr_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* a, int lda );
void bl1_zsyr_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* a, int lda );

// --- syr2 ---

void bl1_ssyr2( uplo1_t uplo, int m, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int a_rs, int a_cs );
void bl1_dsyr2( uplo1_t uplo, int m, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int a_rs, int a_cs );
void bl1_csyr2( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int a_rs, int a_cs );
void bl1_zsyr2( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int a_rs, int a_cs );

void bl1_ssyr2_blas( uplo1_t uplo, int m, float*    alpha, float*    x, int incx, float*    y, int incy, float*    a, int lda );
void bl1_dsyr2_blas( uplo1_t uplo, int m, double*   alpha, double*   x, int incx, double*   y, int incy, double*   a, int lda );
void bl1_csyr2_blas( uplo1_t uplo, int m, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* a, int lda );
void bl1_zsyr2_blas( uplo1_t uplo, int m, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* a, int lda );

// --- trmv ---

void bl1_strmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float*    a, int a_rs, int a_cs, float*    x, int incx );
void bl1_dtrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double*   a, int a_rs, int a_cs, double*   x, int incx );
void bl1_ctrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx );
void bl1_ztrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx );

void bl1_strmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float*    a, int lda, float*    x, int incx );
void bl1_dtrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double*   a, int lda, double*   x, int incx );
void bl1_ctrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* a, int lda, scomplex* x, int incx );
void bl1_ztrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* a, int lda, dcomplex* x, int incx );

// --- trsv ---

void bl1_strsv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float*    a, int a_rs, int a_cs, float*    x, int incx );
void bl1_dtrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double*   a, int a_rs, int a_cs, double*   x, int incx );
void bl1_ctrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx );
void bl1_ztrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx );

void bl1_strsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float*    a, int lda, float*    x, int incx );
void bl1_dtrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double*   a, int lda, double*   x, int incx );
void bl1_ctrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* a, int lda, scomplex* x, int incx );
void bl1_ztrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* a, int lda, dcomplex* x, int incx );

// --- trmvsx ---

void bl1_strmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy );
void bl1_dtrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy );
void bl1_ctrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_ztrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- trsvsx ---

void bl1_strsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, float* alpha, float* a, int a_rs, int a_cs, float* x, int incx, float* beta, float* y, int incy );
void bl1_dtrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, double* alpha, double* a, int a_rs, int a_cs, double* x, int incx, double* beta, double* y, int incy );
void bl1_ctrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_ztrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, int m, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

