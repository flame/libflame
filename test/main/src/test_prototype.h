/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


/* --------BLAS APIs -------- */

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


/* --------LAPACK APIs ---------*/

extern float slamch_(char *);
extern double dlamch_(char *);

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

/* Singular Value Decomposition */
extern int sgesdd_(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* iwork, integer* info);
extern int dgesdd_(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* iwork, integer* info);
extern int cgesdd_(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* iwork, integer* info);
extern int zgesdd_(char* jobz, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* iwork, integer* info);

/* QR factorization*/
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