/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Name-mangling macro definitions -----------------------------------------

// --- Name-mangle level-1 BLAS routines ---------------------------

#define F77_isamax F77_FUNC( isamax , ISAMAX )
#define F77_idamax F77_FUNC( idamax , IDAMAX )
#define F77_icamax F77_FUNC( icamax , ICAMAX )
#define F77_izamax F77_FUNC( izamax , IZAMAX )
#define F77_sasum  F77_FUNC( sasum  , SASUM  )
#define F77_dasum  F77_FUNC( dasum  , DASUM  )
#define F77_scasum F77_FUNC( scasum , SCASUM )
#define F77_dzasum F77_FUNC( dzasum , DZASUM )
#define F77_saxpy  F77_FUNC( saxpy  , SAXPY  )
#define F77_daxpy  F77_FUNC( daxpy  , DAXPY  )
#define F77_caxpy  F77_FUNC( caxpy  , CAXPY  )
#define F77_zaxpy  F77_FUNC( zaxpy  , ZAXPY  )
#define F77_scopy  F77_FUNC( scopy  , SCOPY  )
#define F77_dcopy  F77_FUNC( dcopy  , DCOPY  )
#define F77_ccopy  F77_FUNC( ccopy  , CCOPY  )
#define F77_zcopy  F77_FUNC( zcopy  , ZCOPY  )
#define F77_sdot   F77_FUNC( sdot   , SDOT   )
#define F77_ddot   F77_FUNC( ddot   , DDOT   )
#define F77_cdotu  F77_FUNC( cdotu  , CDOTU  )
#define F77_cdotc  F77_FUNC( cdotc  , CDOTC  )
#define F77_zdotu  F77_FUNC( zdotu  , ZDOTU  )
#define F77_zdotc  F77_FUNC( zdotc  , ZDOTC  )
#define F77_snrm2  F77_FUNC( snrm2  , SNRM2  )
#define F77_dnrm2  F77_FUNC( dnrm2  , DNRM2  )
#define F77_scnrm2 F77_FUNC( scnrm2 , SCNRM2 )
#define F77_dznrm2 F77_FUNC( dznrm2 , DZNRM2 )
#define F77_sscal  F77_FUNC( sscal  , SSCAL  )
#define F77_dscal  F77_FUNC( dscal  , DSCAL  )
#define F77_cscal  F77_FUNC( cscal  , CSCAL  )
#define F77_csscal F77_FUNC( csscal , CSSCAL )
#define F77_zscal  F77_FUNC( zscal  , ZSCAL  )
#define F77_zdscal F77_FUNC( zdscal , ZDSCAL )
#define F77_sswap  F77_FUNC( sswap  , SSWAP  )
#define F77_dswap  F77_FUNC( dswap  , DSWAP  )
#define F77_cswap  F77_FUNC( cswap  , CSWAP  )
#define F77_zswap  F77_FUNC( zswap  , ZSWAP  )

// --- Name-mangle level-2 BLAS routines ---------------------------

#define F77_sgemv  F77_FUNC( sgemv  , SGEMV  )
#define F77_dgemv  F77_FUNC( dgemv  , DGEMV  )
#define F77_cgemv  F77_FUNC( cgemv  , CGEMV  )
#define F77_zgemv  F77_FUNC( zgemv  , ZGEMV  )
#define F77_sger   F77_FUNC( sger   , SGER   )
#define F77_dger   F77_FUNC( dger   , DGER   )
#define F77_cgerc  F77_FUNC( cgerc  , CGERC  )
#define F77_cgeru  F77_FUNC( cgeru  , CGERU  )
#define F77_zgerc  F77_FUNC( zgerc  , ZGERC  )
#define F77_zgeru  F77_FUNC( zgeru  , ZGERU  )
#define F77_chemv  F77_FUNC( chemv  , CHEMV  )
#define F77_zhemv  F77_FUNC( zhemv  , ZHEMV  )
#define F77_cher   F77_FUNC( cher   , CHER   )
#define F77_zher   F77_FUNC( zher   , ZHER   )
#define F77_cher2  F77_FUNC( cher2  , CHER2  )
#define F77_zher2  F77_FUNC( zher2  , ZHER2  )
#define F77_ssymv  F77_FUNC( ssymv  , SSYMV  )
#define F77_dsymv  F77_FUNC( dsymv  , DSYMV  )
#define F77_ssyr   F77_FUNC( ssyr   , SSYR   )
#define F77_dsyr   F77_FUNC( dsyr   , DSYR   )
#define F77_ssyr2  F77_FUNC( ssyr2  , SSYR2  )
#define F77_dsyr2  F77_FUNC( dsyr2  , DSYR2  )
#define F77_strmv  F77_FUNC( strmv  , STRMV  )
#define F77_dtrmv  F77_FUNC( dtrmv  , DTRMV  )
#define F77_ctrmv  F77_FUNC( ctrmv  , CTRMV  )
#define F77_ztrmv  F77_FUNC( ztrmv  , ZTRMV  )
#define F77_strsv  F77_FUNC( strsv  , STRSV  )
#define F77_dtrsv  F77_FUNC( dtrsv  , DTRSV  )
#define F77_ctrsv  F77_FUNC( ctrsv  , CTRSV  )
#define F77_ztrsv  F77_FUNC( ztrsv  , ZTRSV  )

// --- Name-mangle level-3 BLAS routines ---------------------------

#define F77_sgemm  F77_FUNC( sgemm  , SGEMM  )
#define F77_dgemm  F77_FUNC( dgemm  , DGEMM  )
#define F77_cgemm  F77_FUNC( cgemm  , CGEMM  )
#define F77_zgemm  F77_FUNC( zgemm  , ZGEMM  )
#define F77_chemm  F77_FUNC( chemm  , CHEMM  )
#define F77_zhemm  F77_FUNC( zhemm  , ZHEMM  )
#define F77_cherk  F77_FUNC( cherk  , CHERK  )
#define F77_zherk  F77_FUNC( zherk  , ZHERK  )
#define F77_cher2k F77_FUNC( cher2k , CHER2K )
#define F77_zher2k F77_FUNC( zher2k , ZHER2K )
#define F77_ssymm  F77_FUNC( ssymm  , SSYMM  )
#define F77_dsymm  F77_FUNC( dsymm  , DSYMM  )
#define F77_csymm  F77_FUNC( csymm  , CSYMM  )
#define F77_zsymm  F77_FUNC( zsymm  , ZSYMM  )
#define F77_ssyrk  F77_FUNC( ssyrk  , SSYRK  )
#define F77_dsyrk  F77_FUNC( dsyrk  , DSYRK  )
#define F77_csyrk  F77_FUNC( csyrk  , CSYRK  )
#define F77_zsyrk  F77_FUNC( zsyrk  , ZSYRK  )
#define F77_ssyr2k F77_FUNC( ssyr2k , SSYR2K )
#define F77_dsyr2k F77_FUNC( dsyr2k , DSYR2K )
#define F77_csyr2k F77_FUNC( csyr2k , CSYR2K )
#define F77_zsyr2k F77_FUNC( zsyr2k , ZSYR2K )
#define F77_strmm  F77_FUNC( strmm  , STRMM  )
#define F77_dtrmm  F77_FUNC( dtrmm  , DTRMM  )
#define F77_ctrmm  F77_FUNC( ctrmm  , CTRMM  )
#define F77_ztrmm  F77_FUNC( ztrmm  , ZTRMM  )
#define F77_strsm  F77_FUNC( strsm  , STRSM  )
#define F77_dtrsm  F77_FUNC( dtrsm  , DTRSM  )
#define F77_ctrsm  F77_FUNC( ctrsm  , CTRSM  )
#define F77_ztrsm  F77_FUNC( ztrsm  , ZTRSM  )

#ifdef BLIS1_FROM_LIBFLAME
// --- Prototypes --------------------------------------------------------------

// --- Level-1 BLAS prototypes -------------------

// --- amax ---
integer  F77_isamax ( integer* n, float*    x, integer* incx );
integer  F77_idamax ( integer* n, double*   x, integer* incx );
integer  F77_icamax ( integer* n, scomplex* x, integer* incx );
integer  F77_izamax ( integer* n, dcomplex* x, integer* incx );
// --- asum ---
float    F77_sasum  ( integer* n, float*    x, integer* incx );
double   F77_dasum  ( integer* n, double*   x, integer* incx );
float    F77_scasum ( integer* n, scomplex* x, integer* incx );
double   F77_dzasum ( integer* n, dcomplex* x, integer* incx );
// --- axpy ---
void     F77_saxpy  ( integer* n, float*    alpha, float*    x, integer* incx,  float*    y, integer* incy );
void     F77_daxpy  ( integer* n, double*   alpha, double*   x, integer* incx,  double*   y, integer* incy );
void     F77_caxpy  ( integer* n, scomplex* alpha, scomplex* x, integer* incx,  scomplex* y, integer* incy );
void     F77_zaxpy  ( integer* n, dcomplex* alpha, dcomplex* x, integer* incx,  dcomplex* y, integer* incy );
// --- copy ---
void     F77_scopy  ( integer* n, float*    x, integer* incx, float*    y, integer* incy );
void     F77_dcopy  ( integer* n, double*   x, integer* incx, double*   y, integer* incy );
void     F77_ccopy  ( integer* n, scomplex* x, integer* incx, scomplex* y, integer* incy );
void     F77_zcopy  ( integer* n, dcomplex* x, integer* incx, dcomplex* y, integer* incy );
// --- dot ---
float    F77_sdot   ( integer* n, float*    x, integer* incx, float*    y, integer* incy );
double   F77_ddot   ( integer* n, double*   x, integer* incx, double*   y, integer* incy );
scomplex F77_cdotu  ( integer* n, scomplex* x, integer* incx, scomplex* y, integer* incy );
scomplex F77_cdotc  ( integer* n, scomplex* x, integer* incx, scomplex* y, integer* incy );
dcomplex F77_zdotu  ( integer* n, dcomplex* x, integer* incx, dcomplex* y, integer* incy );
dcomplex F77_zdotc  ( integer* n, dcomplex* x, integer* incx, dcomplex* y, integer* incy );
// --- nrm2 ---
float    F77_snrm2  ( integer* n, float*    x, integer* incx );
double   F77_dnrm2  ( integer* n, double*   x, integer* incx );
float    F77_scnrm2 ( integer* n, scomplex* x, integer* incx );
double   F77_dznrm2 ( integer* n, dcomplex* x, integer* incx );
// --- scal ---
void     F77_sscal  ( integer* n, float*    alpha, float*    y, integer* incy );
void     F77_dscal  ( integer* n, double*   alpha, double*   y, integer* incy );
void     F77_cscal  ( integer* n, scomplex* alpha, scomplex* y, integer* incy );
void     F77_csscal ( integer* n, float*    alpha, scomplex* y, integer* incy );
void     F77_zscal  ( integer* n, dcomplex* alpha, dcomplex* y, integer* incy );
void     F77_zdscal ( integer* n, double*   alpha, dcomplex* y, integer* incy );
// --- swap ---
void     F77_sswap  ( integer* n, float*    x, integer* incx, float*    y, integer* incy );
void     F77_dswap  ( integer* n, double*   x, integer* incx, double*   y, integer* incy );
void     F77_cswap  ( integer* n, scomplex* x, integer* incx, scomplex* y, integer* incy );
void     F77_zswap  ( integer* n, dcomplex* x, integer* incx, dcomplex* y, integer* incy );

// --- Level-2 BLAS prototypes -------------------

// --- gemv ---
void     F77_sgemv  ( char* transa, integer* m, integer* n, float*    alpha, float*    a, integer* lda, float*    x, integer* incx, float*    beta, float*    y, integer* incy );
void     F77_dgemv  ( char* transa, integer* m, integer* n, double*   alpha, double*   a, integer* lda, double*   x, integer* incx, double*   beta, double*   y, integer* incy );
void     F77_cgemv  ( char* transa, integer* m, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* x, integer* incx, scomplex* beta, scomplex* y, integer* incy );
void     F77_zgemv  ( char* transa, integer* m, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* x, integer* incx, dcomplex* beta, dcomplex* y, integer* incy );
// --- ger ---
void     F77_sger   ( integer* m, integer* n, float*    alpha, float*    x, integer* incx, float*    y, integer* incy, float*    a, integer* lda );
void     F77_dger   ( integer* m, integer* n, double*   alpha, double*   x, integer* incx, double*   y, integer* incy, double*   a, integer* lda );
void     F77_cgerc  ( integer* m, integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* y, integer* incy, scomplex* a, integer* lda );
void     F77_cgeru  ( integer* m, integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* y, integer* incy, scomplex* a, integer* lda );
void     F77_zgerc  ( integer* m, integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* y, integer* incy, dcomplex* a, integer* lda );
void     F77_zgeru  ( integer* m, integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* y, integer* incy, dcomplex* a, integer* lda );
// --- hemv ---
void     F77_chemv  ( char* uplo, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* x, integer* incx, scomplex* beta, scomplex* y, integer* incy );
void     F77_zhemv  ( char* uplo, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* x, integer* incx, dcomplex* beta, dcomplex* y, integer* incy );
// --- her ---
void     F77_cher   ( char* uplo, integer* n, float*    alpha, scomplex* x, integer* incx, scomplex* a, integer* lda );
void     F77_zher   ( char* uplo, integer* n, double*   alpha, dcomplex* x, integer* incx, dcomplex* a, integer* lda );
// --- her2 ---
void     F77_cher2  ( char* uplo, integer* n, scomplex* alpha, scomplex* x, integer* incx, scomplex* y, integer* incy, scomplex* a, integer* lda );
void     F77_zher2  ( char* uplo, integer* n, dcomplex* alpha, dcomplex* x, integer* incx, dcomplex* y, integer* incy, dcomplex* a, integer* lda );
// --- symv ---
void     F77_ssymv  ( char* uplo, integer* n, float*    alpha, float*    a, integer* lda, float*    x, integer* incx, float*    beta, float*    y, integer* incy );
void     F77_dsymv  ( char* uplo, integer* n, double*   alpha, double*   a, integer* lda, double*   x, integer* incx, double*   beta, double*   y, integer* incy );
// --- syr ---
void     F77_ssyr   ( char* uplo, integer* n, float*    alpha, float*    x, integer* incx, float*    a, integer* lda );
void     F77_dsyr   ( char* uplo, integer* n, double*   alpha, double*   x, integer* incx, double*   a, integer* lda );
// --- syr2 ---
void     F77_ssyr2  ( char* uplo, integer* n, float*    alpha, float*    x, integer* incx, float*    y, integer* incy, float*    a, integer* lda );
void     F77_dsyr2  ( char* uplo, integer* n, double*   alpha, double*   x, integer* incx, double*   y, integer* incy, double*   a, integer* lda );
// --- trmv ---
void     F77_strmv  ( char* uplo, char* transa, char* diag, integer* n,  float*    a, integer* lda, float*    y, integer* incy );
void     F77_dtrmv  ( char* uplo, char* transa, char* diag, integer* n,  double*   a, integer* lda, double*   y, integer* incy );
void     F77_ctrmv  ( char* uplo, char* transa, char* diag, integer* n,  scomplex* a, integer* lda, scomplex* y, integer* incy );
void     F77_ztrmv  ( char* uplo, char* transa, char* diag, integer* n,  dcomplex* a, integer* lda, dcomplex* y, integer* incy );
// --- trsv ---
void     F77_strsv  ( char* uplo, char* transa, char* diag, integer* n,  float*    a, integer* lda, float*    y, integer* incy );
void     F77_dtrsv  ( char* uplo, char* transa, char* diag, integer* n,  double*   a, integer* lda, double*   y, integer* incy );
void     F77_ctrsv  ( char* uplo, char* transa, char* diag, integer* n,  scomplex* a, integer* lda, scomplex* y, integer* incy );
void     F77_ztrsv  ( char* uplo, char* transa, char* diag, integer* n,  dcomplex* a, integer* lda, dcomplex* y, integer* incy );

// --- Level-3 BLAS prototypes -------------------

// --- gemm ---
void     F77_sgemm  ( char* transa, char* transb, integer* m, integer* n, integer* k, float*    alpha, float*    a, integer* lda, float*    b, integer* ldb, float*    beta, float*    c, integer* ldc );
void     F77_dgemm  ( char* transa, char* transb, integer* m, integer* n, integer* k, double*   alpha, double*   a, integer* lda, double*   b, integer* ldb, double*   beta, double*   c, integer* ldc );
void     F77_cgemm  ( char* transa, char* transb, integer* m, integer* n, integer* k, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* beta, scomplex* c, integer* ldc );
void     F77_zgemm  ( char* transa, char* transb, integer* m, integer* n, integer* k, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* beta, dcomplex* c, integer* ldc );
// --- hemm ---
void     F77_chemm  ( char* side, char* uplo, integer* m, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* beta, scomplex* c, integer* ldc );
void     F77_zhemm  ( char* side, char* uplo, integer* m, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* beta, dcomplex* c, integer* ldc );
// --- herk ---
void     F77_cherk  ( char* uplo, char* transa, integer* n, integer* k, float*  alpha, scomplex* a, integer* lda, float*  beta, scomplex* c, integer* ldc );
void     F77_zherk  ( char* uplo, char* transa, integer* n, integer* k, double* alpha, dcomplex* a, integer* lda, double* beta, dcomplex* c, integer* ldc );
// --- her2k ---
void     F77_cher2k ( char* uplo, char* transa, integer* n, integer* k, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, float*  beta, scomplex* c, integer* ldc );
void     F77_zher2k ( char* uplo, char* transa, integer* n, integer* k, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, double* beta, dcomplex* c, integer* ldc );
// --- symm ---
void     F77_ssymm  ( char* side, char* uplo, integer* m, integer* n, float*    alpha, float*    a, integer* lda, float*    b, integer* ldb, float*    beta, float*    c, integer* ldc );
void     F77_dsymm  ( char* side, char* uplo, integer* m, integer* n, double*   alpha, double*   a, integer* lda, double*   b, integer* ldb, double*   beta, double*   c, integer* ldc );
void     F77_csymm  ( char* side, char* uplo, integer* m, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* beta, scomplex* c, integer* ldc );
void     F77_zsymm  ( char* side, char* uplo, integer* m, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* beta, dcomplex* c, integer* ldc );
// --- syrk ---
void     F77_ssyrk  ( char* uplo, char* transa, integer* n, integer* k, float*    alpha, float*    a, integer* lda, float*    beta, float*    c, integer* ldc );
void     F77_dsyrk  ( char* uplo, char* transa, integer* n, integer* k, double*   alpha, double*   a, integer* lda, double*   beta, double*   c, integer* ldc );
void     F77_csyrk  ( char* uplo, char* transa, integer* n, integer* k, scomplex* alpha, scomplex* a, integer* lda, scomplex* beta, scomplex* c, integer* ldc );
void     F77_zsyrk  ( char* uplo, char* transa, integer* n, integer* k, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* beta, dcomplex* c, integer* ldc );
// --- syr2k ---
void     F77_ssyr2k ( char* uplo, char* transa, integer* n, integer* k, float*    alpha, float*    a, integer* lda, float*    b, integer* ldb, float*    beta, float*    c, integer* ldc );
void     F77_dsyr2k ( char* uplo, char* transa, integer* n, integer* k, double*   alpha, double*   a, integer* lda, double*   b, integer* ldb, double*   beta, double*   c, integer* ldc );
void     F77_csyr2k ( char* uplo, char* transa, integer* n, integer* k, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb, scomplex* beta, scomplex* c, integer* ldc );
void     F77_zsyr2k ( char* uplo, char* transa, integer* n, integer* k, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb, dcomplex* beta, dcomplex* c, integer* ldc );
// --- trmm ---
void     F77_strmm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, float*    alpha, float*    a, integer* lda, float*    b, integer* ldb );
void     F77_dtrmm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, double*   alpha, double*   a, integer* lda, double*   b, integer* ldb );
void     F77_ctrmm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb );
void     F77_ztrmm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb );
// --- trsm ---
void     F77_strsm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, float*    alpha, float*    a, integer* lda, float*    b, integer* ldb );
void     F77_dtrsm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, double*   alpha, double*   a, integer* lda, double*   b, integer* ldb );
void     F77_ctrsm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, scomplex* alpha, scomplex* a, integer* lda, scomplex* b, integer* ldb );
void     F77_ztrsm  ( char* side, char* uplo, char* transa, char* diag, integer* m, integer* n, dcomplex* alpha, dcomplex* a, integer* lda, dcomplex* b, integer* ldb );

#endif
