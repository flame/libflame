/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_Svd_ext_u_unb_var1( FLA_Svd_type jobu, FLA_Svd_type jobv, 
                                  dim_t n_iter_max,
                                  FLA_Obj A, FLA_Obj s, FLA_Obj V, FLA_Obj U,
                                  dim_t k_accum,
                                  dim_t b_alg );
int lapack_dbdsqr(char *uplo, integer *n, integer *ncvt, integer *
	          nru, integer *ncc, doublereal *d__, doublereal *e, doublereal *vt, 
	          integer *ldvt, doublereal *u, integer *ldu, doublereal *c__, integer *
	          ldc, doublereal *work, integer *info);
int lapack_dgebd2(integer *m, integer *n, doublereal *a, integer *
	          lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	          taup, doublereal *work, integer *info);
int lapack_dgebrd(integer *m, integer *n, doublereal *a, integer *
	          lda, doublereal *d__, doublereal *e, doublereal *tauq, doublereal *
	          taup, doublereal *work, integer *lwork, integer *info);
int lapack_dgelqf(integer *m, integer *n, doublereal *a, integer *
	          lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
int lapack_dgelq2(integer *m, integer *n, doublereal *a, integer * lda, 
                  doublereal *tau, doublereal *work, integer *info);
int lapack_dgesvd(char *jobu, char *jobvt, integer *m, integer *n, 
	          doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
	          ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork, 
	          integer *info);
int lapack_dorg2r(integer *m, integer *n, integer *k, doublereal *
	          a, integer *lda, doublereal *tau, doublereal *work, integer *info);
int lapack_dorgbr(char *vect, integer *m, integer *n, integer *k, 
	          doublereal *a, integer *lda, doublereal *tau, doublereal *work, 
	          integer *lwork, integer *info);
int lapack_dorgl2(integer *m, integer *n, integer *k, doublereal *
	          a, integer *lda, doublereal *tau, doublereal *work, integer *info);
int lapack_dorglq(integer *m, integer *n, integer *k, doublereal *
	          a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	          integer *info);
int lapack_dorgqr(integer *m, integer *n, integer *k, doublereal *
	          a, integer *lda, doublereal *tau, doublereal *work, integer *lwork, 
	          integer *info);
int lapack_dorm2r(char *side, char *trans, integer *m, integer *n, integer *k, 
                  doublereal *a, integer *lda, doublereal *tau, doublereal * c__, 
                  integer *ldc, doublereal *work, integer *info);
int lapack_dormbr(char *vect, char *side, char *trans, integer *m, 
	          integer *n, integer *k, doublereal *a, integer *lda, doublereal *tau, 
	          doublereal *c__, integer *ldc, doublereal *work, integer *lwork, 
	          integer *info);
int lapack_dormlq(char *side, char *trans, integer *m, integer *n, 
                  integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	          c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
int lapack_dorml2(char *side, char *trans, integer *m, integer *n, integer *k, 
                  doublereal *a, integer *lda, doublereal *tau, doublereal * c__, 
                  integer *ldc, doublereal *work, integer *info); 
int lapack_dormqr(char *side, char *trans, integer *m, integer *n, 
	          integer *k, doublereal *a, integer *lda, doublereal *tau, doublereal *
	          c__, integer *ldc, doublereal *work, integer *lwork, integer *info);
int  dgesvd2x2(   char *jobu, char *jobvt, integer *m, integer *n,
                  doublereal *a, integer *lda, doublereal *s, doublereal *u, integer *
                  ldu, doublereal *vt, integer *ldvt, doublereal *work, integer *lwork,
                  integer *info);



