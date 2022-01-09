/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include <stddef.h>

/*
 * Enumerated and derived types
 */
#define CBLAS_INDEX size_t  /* this may vary between platforms */
enum CBLAS_ORDER     {CblasRowMajor=101, CblasColMajor=102};
enum CBLAS_TRANSPOSE {CblasNoTrans=111, CblasTrans=112, CblasConjTrans=113};
enum CBLAS_UPLO      {CblasUpper=121, CblasLower=122};
enum CBLAS_DIAG      {CblasNonUnit=131, CblasUnit=132};
enum CBLAS_SIDE      {CblasLeft=141, CblasRight=142};

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS functions (complex are recast as routines)
 * ===========================================================================
 */
float  cblas_sdsdot(const integer N, const float alpha, const float *X,
                    const integer incX, const float *Y, const integer incY);
double cblas_dsdot(const integer N, const float *X, const integer incX, const float *Y,
                   const integer incY);
float  cblas_sdot(const integer N, const float  *X, const integer incX,
                  const float  *Y, const integer incY);
double cblas_ddot(const integer N, const double *X, const integer incX,
                  const double *Y, const integer incY);

/*
 * Functions having prefixes Z and C only
 */
void   cblas_cdotu_sub(const integer N, const void *X, const integer incX,
                       const void *Y, const integer incY, void *dotu);
void   cblas_cdotc_sub(const integer N, const void *X, const integer incX,
                       const void *Y, const integer incY, void *dotc);

void   cblas_zdotu_sub(const integer N, const void *X, const integer incX,
                       const void *Y, const integer incY, void *dotu);
void   cblas_zdotc_sub(const integer N, const void *X, const integer incX,
                       const void *Y, const integer incY, void *dotc);


/*
 * Functions having prefixes S D SC DZ
 */
float  cblas_snrm2(const integer N, const float *X, const integer incX);
float  cblas_sasum(const integer N, const float *X, const integer incX);

double cblas_dnrm2(const integer N, const double *X, const integer incX);
double cblas_dasum(const integer N, const double *X, const integer incX);

float  cblas_scnrm2(const integer N, const void *X, const integer incX);
float  cblas_scasum(const integer N, const void *X, const integer incX);

double cblas_dznrm2(const integer N, const void *X, const integer incX);
double cblas_dzasum(const integer N, const void *X, const integer incX);


/*
 * Functions having standard 4 prefixes (S D C Z)
 */
CBLAS_INDEX cblas_isamax(const integer N, const float  *X, const integer incX);
CBLAS_INDEX cblas_idamax(const integer N, const double *X, const integer incX);
CBLAS_INDEX cblas_icamax(const integer N, const void   *X, const integer incX);
CBLAS_INDEX cblas_izamax(const integer N, const void   *X, const integer incX);

/*
 * ===========================================================================
 * Prototypes for level 1 BLAS routines
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (s, d, c, z)
 */
void cblas_sswap(const integer N, float *X, const integer incX, 
                 float *Y, const integer incY);
void cblas_scopy(const integer N, const float *X, const integer incX, 
                 float *Y, const integer incY);
void cblas_saxpy(const integer N, const float alpha, const float *X,
                 const integer incX, float *Y, const integer incY);

void cblas_dswap(const integer N, double *X, const integer incX, 
                 double *Y, const integer incY);
void cblas_dcopy(const integer N, const double *X, const integer incX, 
                 double *Y, const integer incY);
void cblas_daxpy(const integer N, const double alpha, const double *X,
                 const integer incX, double *Y, const integer incY);

void cblas_cswap(const integer N, void *X, const integer incX, 
                 void *Y, const integer incY);
void cblas_ccopy(const integer N, const void *X, const integer incX, 
                 void *Y, const integer incY);
void cblas_caxpy(const integer N, const void *alpha, const void *X,
                 const integer incX, void *Y, const integer incY);

void cblas_zswap(const integer N, void *X, const integer incX, 
                 void *Y, const integer incY);
void cblas_zcopy(const integer N, const void *X, const integer incX, 
                 void *Y, const integer incY);
void cblas_zaxpy(const integer N, const void *alpha, const void *X,
                 const integer incX, void *Y, const integer incY);


/* 
 * Routines with S and D prefix only
 */
void cblas_srotg(float *a, float *b, float *c, float *s);
void cblas_srotmg(float *d1, float *d2, float *b1, const float b2, float *P);
void cblas_srot(const integer N, float *X, const integer incX,
                float *Y, const integer incY, const float c, const float s);
void cblas_srotm(const integer N, float *X, const integer incX,
                float *Y, const integer incY, const float *P);

void cblas_drotg(double *a, double *b, double *c, double *s);
void cblas_drotmg(double *d1, double *d2, double *b1, const double b2, double *P);
void cblas_drot(const integer N, double *X, const integer incX,
                double *Y, const integer incY, const double c, const double s);
void cblas_drotm(const integer N, double *X, const integer incX,
                double *Y, const integer incY, const double *P);


/* 
 * Routines with S D C Z CS and ZD prefixes
 */
void cblas_sscal(const integer N, const float alpha, float *X, const integer incX);
void cblas_dscal(const integer N, const double alpha, double *X, const integer incX);
void cblas_cscal(const integer N, const void *alpha, void *X, const integer incX);
void cblas_zscal(const integer N, const void *alpha, void *X, const integer incX);
void cblas_csscal(const integer N, const float alpha, void *X, const integer incX);
void cblas_zdscal(const integer N, const double alpha, void *X, const integer incX);

/*
 * ===========================================================================
 * Prototypes for level 2 BLAS
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const integer M, const integer N,
                 const float alpha, const float *A, const integer lda,
                 const float *X, const integer incX, const float beta,
                 float *Y, const integer incY);
void cblas_sgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const integer M, const integer N,
                 const integer KL, const integer KU, const float alpha,
                 const float *A, const integer lda, const float *X,
                 const integer incX, const float beta, float *Y, const integer incY);
void cblas_strmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const float *A, const integer lda, 
                 float *X, const integer incX);
void cblas_stbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const integer K, const float *A, const integer lda, 
                 float *X, const integer incX);
void cblas_stpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const float *Ap, float *X, const integer incX);
void cblas_strsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const float *A, const integer lda, float *X,
                 const integer incX);
void cblas_stbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const integer K, const float *A, const integer lda,
                 float *X, const integer incX);
void cblas_stpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const float *Ap, float *X, const integer incX);

void cblas_dgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const integer M, const integer N,
                 const double alpha, const double *A, const integer lda,
                 const double *X, const integer incX, const double beta,
                 double *Y, const integer incY);
void cblas_dgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const integer M, const integer N,
                 const integer KL, const integer KU, const double alpha,
                 const double *A, const integer lda, const double *X,
                 const integer incX, const double beta, double *Y, const integer incY);
void cblas_dtrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const double *A, const integer lda, 
                 double *X, const integer incX);
void cblas_dtbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const integer K, const double *A, const integer lda, 
                 double *X, const integer incX);
void cblas_dtpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const double *Ap, double *X, const integer incX);
void cblas_dtrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const double *A, const integer lda, double *X,
                 const integer incX);
void cblas_dtbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const integer K, const double *A, const integer lda,
                 double *X, const integer incX);
void cblas_dtpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const double *Ap, double *X, const integer incX);

void cblas_cgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const integer M, const integer N,
                 const void *alpha, const void *A, const integer lda,
                 const void *X, const integer incX, const void *beta,
                 void *Y, const integer incY);
void cblas_cgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const integer M, const integer N,
                 const integer KL, const integer KU, const void *alpha,
                 const void *A, const integer lda, const void *X,
                 const integer incX, const void *beta, void *Y, const integer incY);
void cblas_ctrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const void *A, const integer lda, 
                 void *X, const integer incX);
void cblas_ctbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const integer K, const void *A, const integer lda, 
                 void *X, const integer incX);
void cblas_ctpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const void *Ap, void *X, const integer incX);
void cblas_ctrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const void *A, const integer lda, void *X,
                 const integer incX);
void cblas_ctbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const integer K, const void *A, const integer lda,
                 void *X, const integer incX);
void cblas_ctpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const void *Ap, void *X, const integer incX);

void cblas_zgemv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const integer M, const integer N,
                 const void *alpha, const void *A, const integer lda,
                 const void *X, const integer incX, const void *beta,
                 void *Y, const integer incY);
void cblas_zgbmv(const enum CBLAS_ORDER order,
                 const enum CBLAS_TRANSPOSE TransA, const integer M, const integer N,
                 const integer KL, const integer KU, const void *alpha,
                 const void *A, const integer lda, const void *X,
                 const integer incX, const void *beta, void *Y, const integer incY);
void cblas_ztrmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const void *A, const integer lda, 
                 void *X, const integer incX);
void cblas_ztbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const integer K, const void *A, const integer lda, 
                 void *X, const integer incX);
void cblas_ztpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const void *Ap, void *X, const integer incX);
void cblas_ztrsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const void *A, const integer lda, void *X,
                 const integer incX);
void cblas_ztbsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const integer K, const void *A, const integer lda,
                 void *X, const integer incX);
void cblas_ztpsv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE TransA, const enum CBLAS_DIAG Diag,
                 const integer N, const void *Ap, void *X, const integer incX);


/* 
 * Routines with S and D prefixes only
 */
void cblas_ssymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const float alpha, const float *A,
                 const integer lda, const float *X, const integer incX,
                 const float beta, float *Y, const integer incY);
void cblas_ssbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const integer K, const float alpha, const float *A,
                 const integer lda, const float *X, const integer incX,
                 const float beta, float *Y, const integer incY);
void cblas_sspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const float alpha, const float *Ap,
                 const float *X, const integer incX,
                 const float beta, float *Y, const integer incY);
void cblas_sger(const enum CBLAS_ORDER order, const integer M, const integer N,
                const float alpha, const float *X, const integer incX,
                const float *Y, const integer incY, float *A, const integer lda);
void cblas_ssyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const float alpha, const float *X,
                const integer incX, float *A, const integer lda);
void cblas_sspr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const float alpha, const float *X,
                const integer incX, float *Ap);
void cblas_ssyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const float alpha, const float *X,
                const integer incX, const float *Y, const integer incY, float *A,
                const integer lda);
void cblas_sspr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const float alpha, const float *X,
                const integer incX, const float *Y, const integer incY, float *A);

void cblas_dsymv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const double alpha, const double *A,
                 const integer lda, const double *X, const integer incX,
                 const double beta, double *Y, const integer incY);
void cblas_dsbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const integer K, const double alpha, const double *A,
                 const integer lda, const double *X, const integer incX,
                 const double beta, double *Y, const integer incY);
void cblas_dspmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const double alpha, const double *Ap,
                 const double *X, const integer incX,
                 const double beta, double *Y, const integer incY);
void cblas_dger(const enum CBLAS_ORDER order, const integer M, const integer N,
                const double alpha, const double *X, const integer incX,
                const double *Y, const integer incY, double *A, const integer lda);
void cblas_dsyr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const double alpha, const double *X,
                const integer incX, double *A, const integer lda);
void cblas_dspr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const double alpha, const double *X,
                const integer incX, double *Ap);
void cblas_dsyr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const double alpha, const double *X,
                const integer incX, const double *Y, const integer incY, double *A,
                const integer lda);
void cblas_dspr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const double alpha, const double *X,
                const integer incX, const double *Y, const integer incY, double *A);


/* 
 * Routines with C and Z prefixes only
 */
void cblas_chemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const void *alpha, const void *A,
                 const integer lda, const void *X, const integer incX,
                 const void *beta, void *Y, const integer incY);
void cblas_chbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const integer K, const void *alpha, const void *A,
                 const integer lda, const void *X, const integer incX,
                 const void *beta, void *Y, const integer incY);
void cblas_chpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const void *alpha, const void *Ap,
                 const void *X, const integer incX,
                 const void *beta, void *Y, const integer incY);
void cblas_cgeru(const enum CBLAS_ORDER order, const integer M, const integer N,
                 const void *alpha, const void *X, const integer incX,
                 const void *Y, const integer incY, void *A, const integer lda);
void cblas_cgerc(const enum CBLAS_ORDER order, const integer M, const integer N,
                 const void *alpha, const void *X, const integer incX,
                 const void *Y, const integer incY, void *A, const integer lda);
void cblas_cher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const float alpha, const void *X, const integer incX,
                void *A, const integer lda);
void cblas_chpr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const float *alpha, const void *X,
                const integer incX, void *A);
void cblas_cher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const integer N,
                const void *alpha, const void *X, const integer incX,
                const void *Y, const integer incY, void *A, const integer lda);
void cblas_chpr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const integer N,
                const void *alpha, const void *X, const integer incX,
                const void *Y, const integer incY, void *Ap);

void cblas_zhemv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const void *alpha, const void *A,
                 const integer lda, const void *X, const integer incX,
                 const void *beta, void *Y, const integer incY);
void cblas_zhbmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const integer K, const void *alpha, const void *A,
                 const integer lda, const void *X, const integer incX,
                 const void *beta, void *Y, const integer incY);
void cblas_zhpmv(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                 const integer N, const void *alpha, const void *Ap,
                 const void *X, const integer incX,
                 const void *beta, void *Y, const integer incY);
void cblas_zgeru(const enum CBLAS_ORDER order, const integer M, const integer N,
                 const void *alpha, const void *X, const integer incX,
                 const void *Y, const integer incY, void *A, const integer lda);
void cblas_zgerc(const enum CBLAS_ORDER order, const integer M, const integer N,
                 const void *alpha, const void *X, const integer incX,
                 const void *Y, const integer incY, void *A, const integer lda);
void cblas_zher(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const double alpha, const void *X, const integer incX,
                void *A, const integer lda);
void cblas_zhpr(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo,
                const integer N, const double *alpha, const void *X,
                const integer incX, void *A);
void cblas_zher2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const integer N,
                const void *alpha, const void *X, const integer incX,
                const void *Y, const integer incY, void *A, const integer lda);
void cblas_zhpr2(const enum CBLAS_ORDER order, const enum CBLAS_UPLO Uplo, const integer N,
                const void *alpha, const void *X, const integer incX,
                const void *Y, const integer incY, void *Ap);

/*
 * ===========================================================================
 * Prototypes for level 3 BLAS
 * ===========================================================================
 */

/* 
 * Routines with standard 4 prefixes (S, D, C, Z)
 */
void cblas_sgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const integer M, const integer N,
                 const integer K, const float alpha, const float *A,
                 const integer lda, const float *B, const integer ldb,
                 const float beta, float *C, const integer ldc);
void cblas_ssymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const integer M, const integer N,
                 const float alpha, const float *A, const integer lda,
                 const float *B, const integer ldb, const float beta,
                 float *C, const integer ldc);
void cblas_ssyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                 const float alpha, const float *A, const integer lda,
                 const float beta, float *C, const integer ldc);
void cblas_ssyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                  const float alpha, const float *A, const integer lda,
                  const float *B, const integer ldb, const float beta,
                  float *C, const integer ldc);
void cblas_strmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const integer M, const integer N,
                 const float alpha, const float *A, const integer lda,
                 float *B, const integer ldb);
void cblas_strsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const integer M, const integer N,
                 const float alpha, const float *A, const integer lda,
                 float *B, const integer ldb);

void cblas_dgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const integer M, const integer N,
                 const integer K, const double alpha, const double *A,
                 const integer lda, const double *B, const integer ldb,
                 const double beta, double *C, const integer ldc);
void cblas_dsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const integer M, const integer N,
                 const double alpha, const double *A, const integer lda,
                 const double *B, const integer ldb, const double beta,
                 double *C, const integer ldc);
void cblas_dsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                 const double alpha, const double *A, const integer lda,
                 const double beta, double *C, const integer ldc);
void cblas_dsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                  const double alpha, const double *A, const integer lda,
                  const double *B, const integer ldb, const double beta,
                  double *C, const integer ldc);
void cblas_dtrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const integer M, const integer N,
                 const double alpha, const double *A, const integer lda,
                 double *B, const integer ldb);
void cblas_dtrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const integer M, const integer N,
                 const double alpha, const double *A, const integer lda,
                 double *B, const integer ldb);

void cblas_cgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const integer M, const integer N,
                 const integer K, const void *alpha, const void *A,
                 const integer lda, const void *B, const integer ldb,
                 const void *beta, void *C, const integer ldc);
void cblas_csymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const integer M, const integer N,
                 const void *alpha, const void *A, const integer lda,
                 const void *B, const integer ldb, const void *beta,
                 void *C, const integer ldc);
void cblas_csyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                 const void *alpha, const void *A, const integer lda,
                 const void *beta, void *C, const integer ldc);
void cblas_csyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                  const void *alpha, const void *A, const integer lda,
                  const void *B, const integer ldb, const void *beta,
                  void *C, const integer ldc);
void cblas_ctrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const integer M, const integer N,
                 const void *alpha, const void *A, const integer lda,
                 void *B, const integer ldb);
void cblas_ctrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const integer M, const integer N,
                 const void *alpha, const void *A, const integer lda,
                 void *B, const integer ldb);

void cblas_zgemm(const enum CBLAS_ORDER Order, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_TRANSPOSE TransB, const integer M, const integer N,
                 const integer K, const void *alpha, const void *A,
                 const integer lda, const void *B, const integer ldb,
                 const void *beta, void *C, const integer ldc);
void cblas_zsymm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const integer M, const integer N,
                 const void *alpha, const void *A, const integer lda,
                 const void *B, const integer ldb, const void *beta,
                 void *C, const integer ldc);
void cblas_zsyrk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                 const void *alpha, const void *A, const integer lda,
                 const void *beta, void *C, const integer ldc);
void cblas_zsyr2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                  const void *alpha, const void *A, const integer lda,
                  const void *B, const integer ldb, const void *beta,
                  void *C, const integer ldc);
void cblas_ztrmm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const integer M, const integer N,
                 const void *alpha, const void *A, const integer lda,
                 void *B, const integer ldb);
void cblas_ztrsm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const enum CBLAS_TRANSPOSE TransA,
                 const enum CBLAS_DIAG Diag, const integer M, const integer N,
                 const void *alpha, const void *A, const integer lda,
                 void *B, const integer ldb);


/* 
 * Routines with prefixes C and Z only
 */
void cblas_chemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const integer M, const integer N,
                 const void *alpha, const void *A, const integer lda,
                 const void *B, const integer ldb, const void *beta,
                 void *C, const integer ldc);
void cblas_cherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                 const float alpha, const void *A, const integer lda,
                 const float beta, void *C, const integer ldc);
void cblas_cher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                  const void *alpha, const void *A, const integer lda,
                  const void *B, const integer ldb, const float beta,
                  void *C, const integer ldc);

void cblas_zhemm(const enum CBLAS_ORDER Order, const enum CBLAS_SIDE Side,
                 const enum CBLAS_UPLO Uplo, const integer M, const integer N,
                 const void *alpha, const void *A, const integer lda,
                 const void *B, const integer ldb, const void *beta,
                 void *C, const integer ldc);
void cblas_zherk(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                 const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                 const double alpha, const void *A, const integer lda,
                 const double beta, void *C, const integer ldc);
void cblas_zher2k(const enum CBLAS_ORDER Order, const enum CBLAS_UPLO Uplo,
                  const enum CBLAS_TRANSPOSE Trans, const integer N, const integer K,
                  const void *alpha, const void *A, const integer lda,
                  const void *B, const integer ldb, const double beta,
                  void *C, const integer ldc);
