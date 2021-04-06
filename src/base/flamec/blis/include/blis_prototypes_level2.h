/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Level-2 BLAS-like prototypes --------------------------------------------

// --- gemv ---

void bl1_sgemv( trans1_t transa, conj1_t conjx, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    x, integer incx, float*    beta, float*    y, integer incy );
void bl1_dgemv( trans1_t transa, conj1_t conjx, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   x, integer incx, double*   beta, double*   y, integer incy );
void bl1_cgemv( trans1_t transa, conj1_t conjx, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy );
void bl1_zgemv( trans1_t transa, conj1_t conjx, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy );

void bl1_sgemv_blas( trans1_t transa, integer m, integer n, float*    alpha, float*    a, integer lda, float*    x, integer incx, float*    beta, float*    y, integer incy );
void bl1_dgemv_blas( trans1_t transa, integer m, integer n, double*   alpha, double*   a, integer lda, double*   x, integer incx, double*   beta, double*   y, integer incy );
void bl1_cgemv_blas( trans1_t transa, integer m, integer n, scomplex* alpha, scomplex* a, integer lda, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy );
void bl1_zgemv_blas( trans1_t transa, integer m, integer n, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy );

// --- ger ---

void bl1_sger( conj1_t conjx, conj1_t conjy, integer m, integer n, float*    alpha, float*    x, integer incx, float*    y, integer incy, float*    a, integer a_rs, integer a_cs );
void bl1_dger( conj1_t conjx, conj1_t conjy, integer m, integer n, double*   alpha, double*   x, integer incx, double*   y, integer incy, double*   a, integer a_rs, integer a_cs );
void bl1_cger( conj1_t conjx, conj1_t conjy, integer m, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer a_rs, integer a_cs );
void bl1_zger( conj1_t conjx, conj1_t conjy, integer m, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer a_rs, integer a_cs );

void bl1_sger_blas(  integer m, integer n, float*    alpha, float*    x, integer incx, float*    y, integer incy, float*    a, integer lda );
void bl1_dger_blas(  integer m, integer n, double*   alpha, double*   x, integer incx, double*   y, integer incy, double*   a, integer lda );
void bl1_cgerc_blas( integer m, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer lda );
void bl1_cgeru_blas( integer m, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer lda );
void bl1_zgerc_blas( integer m, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer lda );
void bl1_zgeru_blas( integer m, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer lda );

// --- hemv ---

void bl1_shemv( uplo1_t uplo, conj1_t conj, integer m, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    x, integer incx, float*    beta, float*    y, integer incy );
void bl1_dhemv( uplo1_t uplo, conj1_t conj, integer m, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   x, integer incx, double*   beta, double*   y, integer incy );
void bl1_chemv( uplo1_t uplo, conj1_t conj, integer m, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy );
void bl1_zhemv( uplo1_t uplo, conj1_t conj, integer m, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy );

void bl1_chemv_blas( uplo1_t uplo, integer m, scomplex* alpha, scomplex* a, integer lda, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy );
void bl1_zhemv_blas( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy );

// --- her ---

void bl1_sher( uplo1_t uplo, conj1_t conj, integer m, float*  alpha, float*    x, integer incx, float*    a, integer a_rs, integer a_cs );
void bl1_dher( uplo1_t uplo, conj1_t conj, integer m, double* alpha, double*   x, integer incx, double*   a, integer a_rs, integer a_cs );
void bl1_cher( uplo1_t uplo, conj1_t conj, integer m, float*  alpha, scomplex* x, integer incx, scomplex* a, integer a_rs, integer a_cs );
void bl1_zher( uplo1_t uplo, conj1_t conj, integer m, double* alpha, dcomplex* x, integer incx, dcomplex* a, integer a_rs, integer a_cs );

void bl1_cher_blas( uplo1_t uplo, integer m, float*  alpha, scomplex* x, integer incx, scomplex* a, integer lda );
void bl1_zher_blas( uplo1_t uplo, integer m, double* alpha, dcomplex* x, integer incx, dcomplex* a, integer lda );

// --- her2 ---

void bl1_sher2( uplo1_t uplo, conj1_t conj, integer m, float*    alpha, float*    x, integer incx, float*    y, integer incy, float*    a, integer a_rs, integer a_cs );
void bl1_dher2( uplo1_t uplo, conj1_t conj, integer m, double*   alpha, double*   x, integer incx, double*   y, integer incy, double*   a, integer a_rs, integer a_cs );
void bl1_cher2( uplo1_t uplo, conj1_t conj, integer m, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer a_rs, integer a_cs );
void bl1_zher2( uplo1_t uplo, conj1_t conj, integer m, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer a_rs, integer a_cs );

void bl1_cher2_blas( uplo1_t uplo, integer m, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer lda );
void bl1_zher2_blas( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer lda );

// --- symv ---

void bl1_ssymv( uplo1_t uplo, integer m, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    x, integer incx, float*    beta, float*    y, integer incy );
void bl1_dsymv( uplo1_t uplo, integer m, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   x, integer incx, double*   beta, double*   y, integer incy );
void bl1_csymv( uplo1_t uplo, integer m, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy );
void bl1_zsymv( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy );

void bl1_ssymv_blas( uplo1_t uplo, integer m, float*    alpha, float*    a, integer lda, float*    x, integer incx, float*    beta, float*    y, integer incy );
void bl1_dsymv_blas( uplo1_t uplo, integer m, double*   alpha, double*   a, integer lda, double*   x, integer incx, double*   beta, double*   y, integer incy );
void bl1_csymv_blas( uplo1_t uplo, integer m, scomplex* alpha, scomplex* a, integer lda, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy );
void bl1_zsymv_blas( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* a, integer lda, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy );

// --- syr ---

void bl1_ssyr( uplo1_t uplo, integer m, float*    alpha, float*    x, integer incx, float*    a, integer a_rs, integer a_cs );
void bl1_dsyr( uplo1_t uplo, integer m, double*   alpha, double*   x, integer incx, double*   a, integer a_rs, integer a_cs );
void bl1_csyr( uplo1_t uplo, integer m, scomplex* alpha, scomplex* x, integer incx, scomplex* a, integer a_rs, integer a_cs );
void bl1_zsyr( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* a, integer a_rs, integer a_cs );

void bl1_ssyr_blas( uplo1_t uplo, integer m, float*    alpha, float*    x, integer incx, float*    a, integer lda );
void bl1_dsyr_blas( uplo1_t uplo, integer m, double*   alpha, double*   x, integer incx, double*   a, integer lda );
void bl1_csyr_blas( uplo1_t uplo, integer m, scomplex* alpha, scomplex* x, integer incx, scomplex* a, integer lda );
void bl1_zsyr_blas( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* a, integer lda );

// --- syr2 ---

void bl1_ssyr2( uplo1_t uplo, integer m, float*    alpha, float*    x, integer incx, float*    y, integer incy, float*    a, integer a_rs, integer a_cs );
void bl1_dsyr2( uplo1_t uplo, integer m, double*   alpha, double*   x, integer incx, double*   y, integer incy, double*   a, integer a_rs, integer a_cs );
void bl1_csyr2( uplo1_t uplo, integer m, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer a_rs, integer a_cs );
void bl1_zsyr2( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer a_rs, integer a_cs );

void bl1_ssyr2_blas( uplo1_t uplo, integer m, float*    alpha, float*    x, integer incx, float*    y, integer incy, float*    a, integer lda );
void bl1_dsyr2_blas( uplo1_t uplo, integer m, double*   alpha, double*   x, integer incx, double*   y, integer incy, double*   a, integer lda );
void bl1_csyr2_blas( uplo1_t uplo, integer m, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* a, integer lda );
void bl1_zsyr2_blas( uplo1_t uplo, integer m, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* a, integer lda );

// --- trmv ---

void bl1_strmv( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, float*    a, integer a_rs, integer a_cs, float*    x, integer incx );
void bl1_dtrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, double*   a, integer a_rs, integer a_cs, double*   x, integer incx );
void bl1_ctrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, scomplex* a, integer a_rs, integer a_cs, scomplex* x, integer incx );
void bl1_ztrmv( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, dcomplex* a, integer a_rs, integer a_cs, dcomplex* x, integer incx );

void bl1_strmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, float*    a, integer lda, float*    x, integer incx );
void bl1_dtrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, double*   a, integer lda, double*   x, integer incx );
void bl1_ctrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, scomplex* a, integer lda, scomplex* x, integer incx );
void bl1_ztrmv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, dcomplex* a, integer lda, dcomplex* x, integer incx );

// --- trsv ---

void bl1_strsv( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, float*    a, integer a_rs, integer a_cs, float*    x, integer incx );
void bl1_dtrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, double*   a, integer a_rs, integer a_cs, double*   x, integer incx );
void bl1_ctrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, scomplex* a, integer a_rs, integer a_cs, scomplex* x, integer incx );
void bl1_ztrsv( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, dcomplex* a, integer a_rs, integer a_cs, dcomplex* x, integer incx );

void bl1_strsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, float*    a, integer lda, float*    x, integer incx );
void bl1_dtrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, double*   a, integer lda, double*   x, integer incx );
void bl1_ctrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, scomplex* a, integer lda, scomplex* x, integer incx );
void bl1_ztrsv_blas( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, dcomplex* a, integer lda, dcomplex* x, integer incx );

// --- trmvsx ---

void bl1_strmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, float* alpha, float* a, integer a_rs, integer a_cs, float* x, integer incx, float* beta, float* y, integer incy );
void bl1_dtrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, double* alpha, double* a, integer a_rs, integer a_cs, double* x, integer incx, double* beta, double* y, integer incy );
void bl1_ctrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy );
void bl1_ztrmvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy );

// --- trsvsx ---

void bl1_strsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, float* alpha, float* a, integer a_rs, integer a_cs, float* x, integer incx, float* beta, float* y, integer incy );
void bl1_dtrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, double* alpha, double* a, integer a_rs, integer a_cs, double* x, integer incx, double* beta, double* y, integer incy );
void bl1_ctrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy );
void bl1_ztrsvsx( uplo1_t uplo, trans1_t trans, diag1_t diag, integer m, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy );

