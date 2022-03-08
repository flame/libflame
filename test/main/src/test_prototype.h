/*
	Copyright (C) 2022, Advanced Micro Devices, Inc. All rights reserved.
*/


/* --------BLIS APIs -------- */

extern float snrm2_(integer *n, void *x, integer *incx);
extern double dnrm2_(integer *n, void *x, integer *incx);
extern float scnrm2_(integer *n, void *x, integer *incx);
extern double dznrm2_(integer *n, void *x, integer *incx);

extern int sgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);
extern int dgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);
extern int cgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);
extern int zgemv_(char *trans, integer *m, integer *n, void *alpha, void *a, integer *lda, void *x, integer *incx, void *beta, void *y, integer *incy);

extern int strsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, void *alpha, void *a, integer *lda, void *b, integer *ldb);
extern int dtrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, void *alpha, void *a, integer *lda, void *b, integer *ldb);
extern int ctrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, void *alpha, void *a, integer *lda, void *b, integer *ldb);
extern int ztrsm_(char *side, char *uplo, char *transa, char *diag, integer *m, integer *n, void *alpha, void *a, integer *lda, void *b, integer *ldb);

extern int strtrs_(char* uplo, char* trans, char* diag, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);
extern int dtrtrs_(char* uplo, char* trans, char* diag, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);
extern int ctrtrs_(char* uplo, char* trans, char* diag, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);
extern int ztrtrs_(char* uplo, char* trans, char* diag, integer* n, integer* nrhs, void* a, integer* lda, void* b, integer* ldb, integer* info);


/* --------LAPACK APIs ---------*/

extern float slamch_(char *);
extern double dlamch_(char *);

extern int sormrq_(char *side, char *trans, integer *m, integer *n, integer *k, void *a, integer *lda, void *tau, void *c, integer *ldc, void *work, integer *lwork, integer *info);
extern int dormrq_(char *side, char *trans, integer *m, integer *n, integer *k, void *a, integer *lda, void *tau, void *c, integer *ldc, void *work, integer *lwork, integer *info);
extern int cunmrq_(char *side, char *trans, integer *m, integer *n, integer *k, void *a, integer *lda, void *tau, void *c, integer *ldc, void *work, integer *lwork, integer *info);
extern int zunmrq_(char *side, char *trans, integer *m, integer *n, integer *k, void *a, integer *lda, void *tau, void *c, integer *ldc, void *work, integer *lwork, integer *info);

extern int sormqr_(char* side, char* trans, integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* c, integer* ldc, void* work, integer* lwork, integer* info);
extern int dormqr_(char* side, char* trans, integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* c, integer* ldc, void* work, integer* lwork, integer* info);
extern int cunmqr_(char* side, char* trans, integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* c, integer* ldc, void* work, integer* lwork, integer* info);
extern int zunmqr_(char* side, char* trans, integer* m, integer* n, integer* k, void* a, integer* lda, void* tau, void* c, integer* ldc, void* work, integer* lwork, integer* info);

/* SVD APIs*/
extern int sgesvd_(char* jobu, char* jobv, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* info);
extern int dgesvd_(char* jobu, char* jobv, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, integer* info);
extern int cgesvd_(char* jobu, char* jobv, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* info);
extern int zgesvd_(char* jobu, char* jobv, integer* m, integer* n, void* a, integer* lda, void* s, void* u, integer* ldu, void* vt, integer* ldvt, void* work, integer* lwork, void* rwork, integer* info);

/* QR factorization with column pivoting*/
extern int sgeqp3_(integer* m, integer* n, void* a, integer* lda, integer* jpvt, void* tau, void* work, integer* lwork, integer* info);
extern int dgeqp3_(integer* m, integer* n, void* a, integer* lda, integer* jpvt, void* tau, void* work, integer* lwork, integer* info);
extern int cgeqp3_(integer* m, integer* n, void* a, integer* lda, integer* jpvt, void* tau, void* work, integer* lwork, void* rwork, integer* info);
extern int zgeqp3_(integer* m, integer* n, void* a, integer* lda, integer* jpvt, void* tau, void* work, integer* lwork, void* rwork, integer* info);

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