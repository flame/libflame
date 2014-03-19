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


// --- Prototypes --------------------------------------------------------------

// --- Level-1 BLAS prototypes -------------------

// --- amax ---
int      F77_isamax ( int* n, float*    x, int* incx );
int      F77_idamax ( int* n, double*   x, int* incx );
int      F77_icamax ( int* n, scomplex* x, int* incx );
int      F77_izamax ( int* n, dcomplex* x, int* incx );
// --- asum ---
float    F77_sasum  ( int* n, float*    x, int* incx );
double   F77_dasum  ( int* n, double*   x, int* incx );
float    F77_scasum ( int* n, scomplex* x, int* incx );
double   F77_dzasum ( int* n, dcomplex* x, int* incx );
// --- axpy ---
void     F77_saxpy  ( int* n, float*    alpha, float*    x, int* incx,  float*    y, int* incy );
void     F77_daxpy  ( int* n, double*   alpha, double*   x, int* incx,  double*   y, int* incy );
void     F77_caxpy  ( int* n, scomplex* alpha, scomplex* x, int* incx,  scomplex* y, int* incy );
void     F77_zaxpy  ( int* n, dcomplex* alpha, dcomplex* x, int* incx,  dcomplex* y, int* incy );
// --- copy ---
void     F77_scopy  ( int* n, float*    x, int* incx, float*    y, int* incy );
void     F77_dcopy  ( int* n, double*   x, int* incx, double*   y, int* incy );
void     F77_ccopy  ( int* n, scomplex* x, int* incx, scomplex* y, int* incy );
void     F77_zcopy  ( int* n, dcomplex* x, int* incx, dcomplex* y, int* incy );
// --- dot ---
float    F77_sdot   ( int* n, float*    x, int* incx, float*    y, int* incy );
double   F77_ddot   ( int* n, double*   x, int* incx, double*   y, int* incy );
scomplex F77_cdotu  ( int* n, scomplex* x, int* incx, scomplex* y, int* incy );
scomplex F77_cdotc  ( int* n, scomplex* x, int* incx, scomplex* y, int* incy );
dcomplex F77_zdotu  ( int* n, dcomplex* x, int* incx, dcomplex* y, int* incy );
dcomplex F77_zdotc  ( int* n, dcomplex* x, int* incx, dcomplex* y, int* incy );
// --- nrm2 ---
float    F77_snrm2  ( int* n, float*    x, int* incx );
double   F77_dnrm2  ( int* n, double*   x, int* incx );
float    F77_scnrm2 ( int* n, scomplex* x, int* incx );
double   F77_dznrm2 ( int* n, dcomplex* x, int* incx );
// --- scal ---
void     F77_sscal  ( int* n, float*    alpha, float*    y, int* incy );
void     F77_dscal  ( int* n, double*   alpha, double*   y, int* incy );
void     F77_cscal  ( int* n, scomplex* alpha, scomplex* y, int* incy );
void     F77_csscal ( int* n, float*    alpha, scomplex* y, int* incy );
void     F77_zscal  ( int* n, dcomplex* alpha, dcomplex* y, int* incy );
void     F77_zdscal ( int* n, double*   alpha, dcomplex* y, int* incy );
// --- swap ---
void     F77_sswap  ( int* n, float*    x, int* incx, float*    y, int* incy );
void     F77_dswap  ( int* n, double*   x, int* incx, double*   y, int* incy );
void     F77_cswap  ( int* n, scomplex* x, int* incx, scomplex* y, int* incy );
void     F77_zswap  ( int* n, dcomplex* x, int* incx, dcomplex* y, int* incy );

// --- Level-2 BLAS prototypes -------------------

// --- gemv ---
void     F77_sgemv  ( char* transa, int* m, int* n, float*    alpha, float*    a, int* lda, float*    x, int* incx, float*    beta, float*    y, int* incy );
void     F77_dgemv  ( char* transa, int* m, int* n, double*   alpha, double*   a, int* lda, double*   x, int* incx, double*   beta, double*   y, int* incy );
void     F77_cgemv  ( char* transa, int* m, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* x, int* incx, scomplex* beta, scomplex* y, int* incy );
void     F77_zgemv  ( char* transa, int* m, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* x, int* incx, dcomplex* beta, dcomplex* y, int* incy );
// --- ger ---
void     F77_sger   ( int* m, int* n, float*    alpha, float*    x, int* incx, float*    y, int* incy, float*    a, int* lda );
void     F77_dger   ( int* m, int* n, double*   alpha, double*   x, int* incx, double*   y, int* incy, double*   a, int* lda );
void     F77_cgerc  ( int* m, int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* y, int* incy, scomplex* a, int* lda );
void     F77_cgeru  ( int* m, int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* y, int* incy, scomplex* a, int* lda );
void     F77_zgerc  ( int* m, int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* y, int* incy, dcomplex* a, int* lda );
void     F77_zgeru  ( int* m, int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* y, int* incy, dcomplex* a, int* lda );
// --- hemv ---
void     F77_chemv  ( char* uplo, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* x, int* incx, scomplex* beta, scomplex* y, int* incy );
void     F77_zhemv  ( char* uplo, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* x, int* incx, dcomplex* beta, dcomplex* y, int* incy );
// --- her ---
void     F77_cher   ( char* uplo, int* n, float*    alpha, scomplex* x, int* incx, scomplex* a, int* lda );
void     F77_zher   ( char* uplo, int* n, double*   alpha, dcomplex* x, int* incx, dcomplex* a, int* lda );
// --- her2 ---
void     F77_cher2  ( char* uplo, int* n, scomplex* alpha, scomplex* x, int* incx, scomplex* y, int* incy, scomplex* a, int* lda );
void     F77_zher2  ( char* uplo, int* n, dcomplex* alpha, dcomplex* x, int* incx, dcomplex* y, int* incy, dcomplex* a, int* lda );
// --- symv ---
void     F77_ssymv  ( char* uplo, int* n, float*    alpha, float*    a, int* lda, float*    x, int* incx, float*    beta, float*    y, int* incy );
void     F77_dsymv  ( char* uplo, int* n, double*   alpha, double*   a, int* lda, double*   x, int* incx, double*   beta, double*   y, int* incy );
// --- syr ---
void     F77_ssyr   ( char* uplo, int* n, float*    alpha, float*    x, int* incx, float*    a, int* lda );
void     F77_dsyr   ( char* uplo, int* n, double*   alpha, double*   x, int* incx, double*   a, int* lda );
// --- syr2 ---
void     F77_ssyr2  ( char* uplo, int* n, float*    alpha, float*    x, int* incx, float*    y, int* incy, float*    a, int* lda );
void     F77_dsyr2  ( char* uplo, int* n, double*   alpha, double*   x, int* incx, double*   y, int* incy, double*   a, int* lda );
// --- trmv ---
void     F77_strmv  ( char* uplo, char* transa, char* diag, int* n,  float*    a, int* lda, float*    y, int* incy );
void     F77_dtrmv  ( char* uplo, char* transa, char* diag, int* n,  double*   a, int* lda, double*   y, int* incy );
void     F77_ctrmv  ( char* uplo, char* transa, char* diag, int* n,  scomplex* a, int* lda, scomplex* y, int* incy );
void     F77_ztrmv  ( char* uplo, char* transa, char* diag, int* n,  dcomplex* a, int* lda, dcomplex* y, int* incy );
// --- trsv ---
void     F77_strsv  ( char* uplo, char* transa, char* diag, int* n,  float*    a, int* lda, float*    y, int* incy );
void     F77_dtrsv  ( char* uplo, char* transa, char* diag, int* n,  double*   a, int* lda, double*   y, int* incy );
void     F77_ctrsv  ( char* uplo, char* transa, char* diag, int* n,  scomplex* a, int* lda, scomplex* y, int* incy );
void     F77_ztrsv  ( char* uplo, char* transa, char* diag, int* n,  dcomplex* a, int* lda, dcomplex* y, int* incy );

// --- Level-3 BLAS prototypes -------------------

// --- gemm ---
void     F77_sgemm  ( char* transa, char* transb, int* m, int* n, int* k, float*    alpha, float*    a, int* lda, float*    b, int* ldb, float*    beta, float*    c, int* ldc );
void     F77_dgemm  ( char* transa, char* transb, int* m, int* n, int* k, double*   alpha, double*   a, int* lda, double*   b, int* ldb, double*   beta, double*   c, int* ldc );
void     F77_cgemm  ( char* transa, char* transb, int* m, int* n, int* k, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* beta, scomplex* c, int* ldc );
void     F77_zgemm  ( char* transa, char* transb, int* m, int* n, int* k, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* beta, dcomplex* c, int* ldc );
// --- hemm ---
void     F77_chemm  ( char* side, char* uplo, int* m, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* beta, scomplex* c, int* ldc );
void     F77_zhemm  ( char* side, char* uplo, int* m, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* beta, dcomplex* c, int* ldc );
// --- herk ---
void     F77_cherk  ( char* uplo, char* transa, int* n, int* k, float*  alpha, scomplex* a, int* lda, float*  beta, scomplex* c, int* ldc );
void     F77_zherk  ( char* uplo, char* transa, int* n, int* k, double* alpha, dcomplex* a, int* lda, double* beta, dcomplex* c, int* ldc );
// --- her2k ---
void     F77_cher2k ( char* uplo, char* transa, int* n, int* k, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb, float*  beta, scomplex* c, int* ldc );
void     F77_zher2k ( char* uplo, char* transa, int* n, int* k, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb, double* beta, dcomplex* c, int* ldc );
// --- symm ---
void     F77_ssymm  ( char* side, char* uplo, int* m, int* n, float*    alpha, float*    a, int* lda, float*    b, int* ldb, float*    beta, float*    c, int* ldc );
void     F77_dsymm  ( char* side, char* uplo, int* m, int* n, double*   alpha, double*   a, int* lda, double*   b, int* ldb, double*   beta, double*   c, int* ldc );
void     F77_csymm  ( char* side, char* uplo, int* m, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* beta, scomplex* c, int* ldc );
void     F77_zsymm  ( char* side, char* uplo, int* m, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* beta, dcomplex* c, int* ldc );
// --- syrk ---
void     F77_ssyrk  ( char* uplo, char* transa, int* n, int* k, float*    alpha, float*    a, int* lda, float*    beta, float*    c, int* ldc );
void     F77_dsyrk  ( char* uplo, char* transa, int* n, int* k, double*   alpha, double*   a, int* lda, double*   beta, double*   c, int* ldc );
void     F77_csyrk  ( char* uplo, char* transa, int* n, int* k, scomplex* alpha, scomplex* a, int* lda, scomplex* beta, scomplex* c, int* ldc );
void     F77_zsyrk  ( char* uplo, char* transa, int* n, int* k, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* beta, dcomplex* c, int* ldc );
// --- syr2k ---
void     F77_ssyr2k ( char* uplo, char* transa, int* n, int* k, float*    alpha, float*    a, int* lda, float*    b, int* ldb, float*    beta, float*    c, int* ldc );
void     F77_dsyr2k ( char* uplo, char* transa, int* n, int* k, double*   alpha, double*   a, int* lda, double*   b, int* ldb, double*   beta, double*   c, int* ldc );
void     F77_csyr2k ( char* uplo, char* transa, int* n, int* k, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb, scomplex* beta, scomplex* c, int* ldc );
void     F77_zsyr2k ( char* uplo, char* transa, int* n, int* k, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb, dcomplex* beta, dcomplex* c, int* ldc );
// --- trmm ---
void     F77_strmm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, float*    alpha, float*    a, int* lda, float*    b, int* ldb );
void     F77_dtrmm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, double*   alpha, double*   a, int* lda, double*   b, int* ldb );
void     F77_ctrmm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb );
void     F77_ztrmm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb );
// --- trsm ---
void     F77_strsm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, float*    alpha, float*    a, int* lda, float*    b, int* ldb );
void     F77_dtrsm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, double*   alpha, double*   a, int* lda, double*   b, int* ldb );
void     F77_ctrsm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, scomplex* alpha, scomplex* a, int* lda, scomplex* b, int* ldb );
void     F77_ztrsm  ( char* side, char* uplo, char* transa, char* diag, int* m, int* n, dcomplex* alpha, dcomplex* a, int* lda, dcomplex* b, int* ldb );

