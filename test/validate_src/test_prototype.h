/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/

#ifndef TEST_PROTOTYPE_H
#define TEST_PROTOTYPE_H

/* These functions are API invoking functions used in other API test codes */
extern void invoke_getrf(integer datatype, integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);
extern void invoke_potrf(char* uplo, integer datatype, integer* m, void* a, integer* lda, integer* info);
extern void invoke_geqrf(integer datatype, integer* m, integer* n, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
/* --------BLAS APIs -------- */

extern int saxpy_(integer *n, void *a, void *x, integer *incx, void *y, integer *incy);
extern int daxpy_(integer *n, void *a, void *x, integer *incx, void *y, integer *incy);
extern int caxpy_(integer *n, void *a, void *x, integer *incx, void *y, integer *incy);
extern int zaxpy_(integer *n, void *a, void *x, integer *incx, void *y, integer *incy);

extern float snrm2_(integer *n, void *x, integer *incx);
extern double dnrm2_(integer *n, void *x, integer *incx);
extern float scnrm2_(integer *n, void *x, integer *incx);
extern double dznrm2_(integer *n, void *x, integer *incx);

extern int scopy_(integer *n, void *x, integer *incx, void *y, integer *incy);
extern int dcopy_(integer *n, void *x, integer *incx, void *y, integer *incy);
extern int ccopy_(integer *n, void *x, integer *incx, void *y, integer *incy);
extern int zcopy_(integer *n, void *x, integer *incx, void *y, integer *incy);

extern int sgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);
extern int dgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);
extern int cgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);
extern int zgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);

extern int sgemm_(char *transa, char *transb, integer *m, integer * n, integer *k, void *alpha, void *a, integer *lda, void *b, integer * ldb, void *beta, void *c__, integer *ldc);
extern int dgemm_(char *transa, char *transb, integer *m, integer * n, integer *k, void *alpha, void *a, integer *lda, void *b, integer * ldb, void *beta, void *c__, integer *ldc);
extern int cgemm_(char *transa, char *transb, integer *m, integer * n, integer *k, void *alpha, void *a, integer *lda, void *b, integer * ldb, void *beta, void *c__, integer *ldc);
extern int zgemm_(char *transa, char *transb, integer *m, integer * n, integer *k, void *alpha, void *a, integer *lda, void *b, integer * ldb, void *beta, void *c__, integer *ldc);

extern int spotrs_(char* uplo, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);
extern int dpotrs_(char* uplo, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);
extern int cpotrs_(char* uplo, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);
extern int zpotrs_(char* uplo, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);

/* --------LAPACK APIs ---------*/

extern float slamch_(char *);
extern double dlamch_(char *);

extern int claswp_(integer* n, scomplex* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx);
extern int dlaswp_(integer* n, double* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx);
extern int slaswp_(integer* n, float* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx);
extern int zlaswp_(integer* n, dcomplex* a, integer* lda, integer* k1, integer* k2, integer* ipiv, integer* incx);


extern float slange_(char* norm, integer* m, integer* n, void* a, integer* lda, void* work);
extern double dlange_(char* norm, integer* m, integer* n, void* a, integer* lda, void* work);
extern float clange_(char* norm, integer* m, integer* n, void* a, integer* lda, void* work);
extern double zlange_(char* norm, integer* m, integer* n, void* a, integer* lda, void* work);

extern int slaset_(char *uplo, integer *m, integer *n, void *alpha, void *beta, void *a, integer *lda);
extern int dlaset_(char *uplo, integer *m, integer *n, void *alpha, void *beta, void *a, integer *lda);
extern int claset_(char *uplo, integer *m, integer *n, void *alpha, void *beta, void *a, integer *lda);
extern int zlaset_(char *uplo, integer *m, integer *n, void *alpha, void *beta, void *a, integer *lda);

extern int slacpy_(char* uplo, integer* m, integer* n, void* a, integer* lda, void* b, integer* ldb);
extern int dlacpy_(char* uplo, integer* m, integer* n, void* a, integer* lda, void* b, integer* ldb);
extern int clacpy_(char* uplo, integer* m, integer* n, void* a, integer* lda, void* b, integer* ldb);
extern int zlacpy_(char* uplo, integer* m, integer* n, void* a, integer* lda, void* b, integer* ldb);

extern int sorgrq_(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int dorgrq_(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int cungrq_(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int zungrq_(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);

extern int sorgqr_(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int dorgqr_(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int cungqr_(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);
extern int zungqr_(integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* work, integer* lwork, integer* info);

/* Eigen value and Eigen vectors*/
extern int sgeevx_(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* wr, void* wi, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* scale, void* abnrm, void* rconde, void* rcondv, void* work, integer* lwork, integer* iwork, integer* info);
extern int dgeevx_(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* wr, void* wi, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* scale, void* abnrm, void* rconde, void* rcondv, void* work, integer* lwork, integer* iwork, integer* info);
extern int cgeevx_(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* w, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* scale, void* abnrm, void* rconde, void* rcondv, void* work, integer* lwork, void* rwork, integer* info);
extern int zgeevx_(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* w, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* scale, void* abnrm, void* rconde, void* rcondv, void* work, integer* lwork, void* rwork, integer* info);

int sgeev_(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* wr, void* wi, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, integer* info);
int dgeev_(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* wr, void* wi, void* vl, integer* ldvl, void* vr, void* ldvr, void* work, integer* lwork, integer* info);
int cgeev_(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* w, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, void* rwork, integer* info);
int zgeev_(char* jobvl, char* jobvr, integer* n, void* a, integer* lda, void* w, void* vl, integer* ldvl, void* vr, integer* ldvr, void* work, integer* lwork, void* rwork, integer* info);

/* Singular Value Decomposition */
extern int sgesdd_(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* iwork, integer* info);
extern int dgesdd_(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* iwork, integer* info);
extern int cgesdd_(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* iwork, integer* info);
extern int zgesdd_(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* iwork, integer* info);

extern int sgesvd_(char* jobu, char* jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* info);
extern int dgesvd_(char* jobu, char *jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* info);
extern int cgesvd_(char* jobu, char* jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* rwork, integer* info);
extern int zgesvd_(char* jobu, char* jobvt, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* rwork, integer* info);

/* QR factorization */
extern int sgeqrf_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int dgeqrf_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int cgeqrf_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int zgeqrf_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);

/* RQ factorization APIs */
extern int sgerqf_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int dgerqf_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int cgerqf_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);
extern int zgerqf_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *lwork, integer *info);

/* RQ factorization with unblocked algorithm APIs*/
extern int sgerq2_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info);
extern int dgerq2_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info);
extern int cgerq2_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info);
extern int zgerq2_(integer *m, integer *n, void *a, integer *lda, void *tau, void *work, integer *info);

/* QR factorization with column pivoting */
extern int sgeqp3_(integer *m, integer *n, void *a, integer *lda, integer *jpvt, void *tau, void *work, integer *lwork, integer *info);
extern int dgeqp3_(integer *m, integer *n, void *a, integer *lda, integer *jpvt, void *tau, void *work, integer *lwork, integer *info);
extern int cgeqp3_(integer *m, integer *n, void *a, integer *lda, integer *jpvt, void *tau, void *work, integer *lwork, void *rwork, integer *info);
extern int zgeqp3_(integer *m, integer *n, void *a, integer *lda, integer *jpvt, void *tau, void *work, integer *lwork, void *rwork, integer *info);

/* Cholesky factorization APIS*/
extern int spotrf_(char* uplo, integer* n, void* a, integer* lda, integer* info);
extern int dpotrf_(char* uplo, integer* n, void* a, integer* lda, integer* info);
extern int cpotrf_(char* uplo, integer* n, void* a, integer* lda, integer* info);
extern int zpotrf_(char* uplo, integer* n, void* a, integer* lda, integer* info);

/* LU factorization of a general m by n matrix */
extern int dgetrf_(integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);
extern int sgetrf_(integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);
extern int cgetrf_(integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);
extern int zgetrf_(integer* m, integer* n, void* a, integer* lda, integer* ipiv, integer* info);

/* LU factorization of a general m by n matrix */
extern int dgetri_(integer* n, void* a, integer* lda, integer* ipiv, void *work, integer *lwork, integer* info);
extern int sgetri_(integer* n, void* a, integer* lda, integer* ipiv, void* work, integer* lwork, integer* info);
extern int cgetri_(integer* n, void* a, integer* lda, integer* ipiv, void* work, integer* lwork, integer* info);
extern int zgetri_(integer* n, void* a, integer* lda, integer* ipiv, void* work, integer* lwork, integer* info);

/* LU factorization of a general m by n matrix */
extern int dgetrs_(char* trans, integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int sgetrs_(char* trans, integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int cgetrs_(char* trans, integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int zgetrs_(char* trans, integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);

/* Computation of Eigen Values and Eigen Vectors */
extern int ssyevd_(char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, integer* iwork, integer* liwork, integer* info);
extern int dsyevd_(char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, integer* iwork, integer* liwork, integer* info);
extern int cheevd_(char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, void* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info);
extern int zheevd_(char* jobz, char* uplo, integer* n, void* a, integer* lda, void* w, void* work, integer* lwork, void* rwork, integer* lrwork, integer* iwork, integer* liwork, integer* info);
/* Computation of Eigen Values and Eigen Vectors */
extern int sggevx_(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alphar, void* alphai, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, void* rcondv, void* work, integer* lwork, integer* iwork, integer* bwork, integer* info);
extern int dggevx_(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alphar, void* alphai, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, void* rcondv, void* work, integer* lwork, integer* iwork, integer* bwork, integer* info);
extern int cggevx_(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alpha, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, void* rcondv, void* work, integer* lwork, void* rwork, integer* iwork, integer* bwork, integer* info);
extern int zggevx_(char* balanc, char* jobvl, char* jobvr, char* sense, integer* n, void* a, integer* lda, void* b, integer* ldb, void* alpha, void* beta, void* vl, integer* ldvl, void* vr, integer* ldvr, integer* ilo, integer* ihi, void* lscale, void* rscale, void* abnrm, void* bbnrm, void* rconde, void* rcondv, void* work, integer* lwork, void* rwork, integer* iwork, integer* bwork, integer* info);

extern int sgesv_(integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int dgesv_(integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int cgesv_(integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);
extern int zgesv_(integer* n, integer* nrhs, void* a, integer* lda, integer* ipiv, void* b, integer* ldb, integer* info);

#endif  // TEST_PROTOTYPE_H
