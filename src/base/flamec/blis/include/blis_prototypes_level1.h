/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Level-1 BLAS-like prototypes --------------------------------------------

// --- amax ---

void bl1_samax( int n, float*    x, int incx, int* index );
void bl1_damax( int n, double*   x, int incx, int* index );
void bl1_camax( int n, scomplex* x, int incx, int* index );
void bl1_zamax( int n, dcomplex* x, int incx, int* index );

// --- asum ---

void bl1_sasum( int n, float*    x, int incx, float*  norm );
void bl1_dasum( int n, double*   x, int incx, double* norm );
void bl1_casum( int n, scomplex* x, int incx, float*  norm );
void bl1_zasum( int n, dcomplex* x, int incx, double* norm );

// --- axpy ---

void bl1_saxpy( int n, float*    alpha, float*    x, int incx, float*    y, int incy );
void bl1_daxpy( int n, double*   alpha, double*   x, int incx, double*   y, int incy );
void bl1_caxpy( int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zaxpy( int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy );

// --- axpyv ---

void bl1_saxpyv( conj1_t conj, int n, float*    alpha, float*    x, int incx, float*    y, int incy );
void bl1_daxpyv( conj1_t conj, int n, double*   alpha, double*   x, int incx, double*   y, int incy );
void bl1_caxpyv( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zaxpyv( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy );

// --- axpymt ---

void bl1_saxpymt( trans1_t trans, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_daxpymt( trans1_t trans, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_caxpymt( trans1_t trans, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zaxpymt( trans1_t trans, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- axpymrt ---

void bl1_saxpymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_daxpymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_caxpymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zaxpymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- axpysv ---

void bl1_saxpysv( int n, float*    alpha0, float*    alpha1, float*    x, int incx, float*    beta, float*    y, int incy );
void bl1_daxpysv( int n, double*   alpha0, double*   alpha1, double*   x, int incx, double*   beta, double*   y, int incy );
void bl1_caxpysv( int n, scomplex* alpha0, scomplex* alpha1, scomplex* x, int incx, scomplex* beta, scomplex* y, int incy );
void bl1_zaxpysv( int n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* x, int incx, dcomplex* beta, dcomplex* y, int incy );

// --- axpysmt ---

void bl1_saxpysmt( trans1_t trans, int m, int n, float*    alpha0, float*    alpha1, float*    a, int a_rs, int a_cs, float*    beta, float*    b, int b_rs, int b_cs );
void bl1_daxpysmt( trans1_t trans, int m, int n, double*   alpha0, double*   alpha1, double*   a, int a_rs, int a_cs, double*   beta, double*   b, int b_rs, int b_cs );
void bl1_caxpysmt( trans1_t trans, int m, int n, scomplex* alpha0, scomplex* alpha1, scomplex* a, int a_rs, int a_cs, scomplex* beta, scomplex* b, int b_rs, int b_cs );
void bl1_zaxpysmt( trans1_t trans, int m, int n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* a, int a_rs, int a_cs, dcomplex* beta, dcomplex* b, int b_rs, int b_cs );

// --- conjv ---

void bl1_sconjv( int m, float* x, int incx );
void bl1_dconjv( int m, double* x, int incx );
void bl1_cconjv( int m, scomplex* x, int incx );
void bl1_zconjv( int m, dcomplex* x, int incx );

// --- conjm ---

void bl1_sconjm( int m, int n, float*    a, int a_rs, int a_cs );
void bl1_dconjm( int m, int n, double*   a, int a_rs, int a_cs );
void bl1_cconjm( int m, int n, scomplex* a, int a_rs, int a_cs );
void bl1_zconjm( int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- conjmr ---

void bl1_sconjmr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs );
void bl1_dconjmr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs );
void bl1_cconjmr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs );
void bl1_zconjmr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs );

// --- copy ---

void bl1_scopy( int m, float*    x, int incx, float*    y, int incy );
void bl1_dcopy( int m, double*   x, int incx, double*   y, int incy );
void bl1_ccopy( int m, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zcopy( int m, dcomplex* x, int incx, dcomplex* y, int incy );

// --- copyv ---

void bl1_icopyv( conj1_t conj, int m, int*      x, int incx, int*      y, int incy );
void bl1_scopyv( conj1_t conj, int m, float*    x, int incx, float*    y, int incy );
void bl1_dcopyv( conj1_t conj, int m, double*   x, int incx, double*   y, int incy );
void bl1_ccopyv( conj1_t conj, int m, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zcopyv( conj1_t conj, int m, dcomplex* x, int incx, dcomplex* y, int incy );

void bl1_sdcopyv( conj1_t conj, int m, float*    x, int incx, double*   y, int incy );
void bl1_dscopyv( conj1_t conj, int m, double*   x, int incx, float*    y, int incy );
void bl1_sccopyv( conj1_t conj, int m, float*    x, int incx, scomplex* y, int incy );
void bl1_cscopyv( conj1_t conj, int m, scomplex* x, int incx, float*    y, int incy );
void bl1_szcopyv( conj1_t conj, int m, float*    x, int incx, dcomplex* y, int incy );
void bl1_zscopyv( conj1_t conj, int m, dcomplex* x, int incx, float*    y, int incy );
void bl1_dccopyv( conj1_t conj, int m, double*   x, int incx, scomplex* y, int incy );
void bl1_cdcopyv( conj1_t conj, int m, scomplex* x, int incx, double*   y, int incy );
void bl1_dzcopyv( conj1_t conj, int m, double*   x, int incx, dcomplex* y, int incy );
void bl1_zdcopyv( conj1_t conj, int m, dcomplex* x, int incx, double*   y, int incy );
void bl1_czcopyv( conj1_t conj, int m, scomplex* x, int incx, dcomplex* y, int incy );
void bl1_zccopyv( conj1_t conj, int m, dcomplex* x, int incx, scomplex* y, int incy );

// --- copymr ---

void bl1_scopymr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dcopymr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_ccopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zcopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bl1_sscopymr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_sdcopymr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dscopymr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_sccopymr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cscopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_szcopymr( uplo1_t uplo, int m, int n, float*    a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zscopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_ddcopymr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dccopymr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cdcopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dzcopymr( uplo1_t uplo, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zdcopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_cccopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_czcopymr( uplo1_t uplo, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zccopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zzcopymr( uplo1_t uplo, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- copymrt ---

void bl1_scopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_ccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bl1_sscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_sdcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_sccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_szcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_dscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_ddcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_dzcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_cscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_cdcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_cccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_czcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zscopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_zdcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_zccopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zzcopymrt( uplo1_t uplo, trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- copymt ---

void bl1_icopymt( trans1_t trans, int m, int n, int*      a, int a_rs, int a_cs, int*      b, int b_rs, int b_cs );
void bl1_scopymt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dcopymt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_ccopymt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zcopymt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

void bl1_sscopymt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_sdcopymt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dscopymt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_sccopymt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cscopymt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_szcopymt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zscopymt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_ddcopymt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dccopymt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_cdcopymt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_dzcopymt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zdcopymt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_cccopymt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_czcopymt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );
void bl1_zccopymt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zzcopymt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

// --- dot ---

void bl1_cdot_in( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy, scomplex* rho );
void bl1_zdot_in( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* rho );

void bl1_sdot( conj1_t conj, int n, float*    x, int incx, float*    y, int incy, float*    rho );
void bl1_ddot( conj1_t conj, int n, double*   x, int incx, double*   y, int incy, double*   rho );
void bl1_cdot( conj1_t conj, int n, scomplex* x, int incx, scomplex* y, int incy, scomplex* rho );
void bl1_zdot( conj1_t conj, int n, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* rho );

// --- dots ---

void bl1_sdots( conj1_t conj, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    beta, float*    rho );
void bl1_ddots( conj1_t conj, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   beta, double*   rho );
void bl1_cdots( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* beta, scomplex* rho );
void bl1_zdots( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* beta, dcomplex* rho );

// --- dot2s ---

void bl1_sdot2s( conj1_t conj, int n, float*    alpha, float*    x, int incx, float*    y, int incy, float*    beta, float*    rho );
void bl1_ddot2s( conj1_t conj, int n, double*   alpha, double*   x, int incx, double*   y, int incy, double*   beta, double*   rho );
void bl1_cdot2s( conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx, scomplex* y, int incy, scomplex* beta, scomplex* rho );
void bl1_zdot2s( conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx, dcomplex* y, int incy, dcomplex* beta, dcomplex* rho );

// --- fnorm ---

void bl1_sfnorm( int m, int n, float*    a, int a_rs, int a_cs, float*  norm );
void bl1_dfnorm( int m, int n, double*   a, int a_rs, int a_cs, double* norm );
void bl1_cfnorm( int m, int n, scomplex* a, int a_rs, int a_cs, float*  norm );
void bl1_zfnorm( int m, int n, dcomplex* a, int a_rs, int a_cs, double* norm );

// --- invscalv ---

void bl1_sinvscalv(  conj1_t conj, int n, float*    alpha, float*    x, int incx );
void bl1_dinvscalv(  conj1_t conj, int n, double*   alpha, double*   x, int incx );
void bl1_csinvscalv( conj1_t conj, int n, float*    alpha, scomplex* x, int incx );
void bl1_cinvscalv(  conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx );
void bl1_zdinvscalv( conj1_t conj, int n, double*   alpha, dcomplex* x, int incx );
void bl1_zinvscalv(  conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx );

// --- invscalm ---

void bl1_sinvscalm(  conj1_t conj, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs );
void bl1_dinvscalm(  conj1_t conj, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs );
void bl1_csinvscalm( conj1_t conj, int m, int n, float*    alpha, scomplex* a, int a_rs, int a_cs );
void bl1_cinvscalm(  conj1_t conj, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs );
void bl1_zdinvscalm( conj1_t conj, int m, int n, double*   alpha, dcomplex* a, int a_rs, int a_cs );
void bl1_zinvscalm(  conj1_t conj, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs );

// --- nrm2 ---

void bl1_snrm2( int n, float*    x, int incx, float*  norm );
void bl1_dnrm2( int n, double*   x, int incx, double* norm );
void bl1_cnrm2( int n, scomplex* x, int incx, float*  norm );
void bl1_znrm2( int n, dcomplex* x, int incx, double* norm );

// --- scal ---

void bl1_sscal(  int n, float*    alpha, float*    x, int incx );
void bl1_dscal(  int n, double*   alpha, double*   x, int incx );
void bl1_csscal( int n, float*    alpha, scomplex* x, int incx );
void bl1_cscal(  int n, scomplex* alpha, scomplex* x, int incx );
void bl1_zdscal( int n, double*   alpha, dcomplex* x, int incx );
void bl1_zscal(  int n, dcomplex* alpha, dcomplex* x, int incx );

// --- scalv ---

void bl1_sscalv(  conj1_t conj, int n, float*    alpha, float*    x, int incx );
void bl1_dscalv(  conj1_t conj, int n, double*   alpha, double*   x, int incx );
void bl1_csscalv( conj1_t conj, int n, float*    alpha, scomplex* x, int incx );
void bl1_cscalv(  conj1_t conj, int n, scomplex* alpha, scomplex* x, int incx );
void bl1_zdscalv( conj1_t conj, int n, double*   alpha, dcomplex* x, int incx );
void bl1_zscalv(  conj1_t conj, int n, dcomplex* alpha, dcomplex* x, int incx );

// --- scalm ---

void bl1_sscalm(  conj1_t conj, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs );
void bl1_dscalm(  conj1_t conj, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs );
void bl1_csscalm( conj1_t conj, int m, int n, float*    alpha, scomplex* a, int a_rs, int a_cs );
void bl1_cscalm(  conj1_t conj, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs );
void bl1_zdscalm( conj1_t conj, int m, int n, double*   alpha, dcomplex* a, int a_rs, int a_cs );
void bl1_zscalm(  conj1_t conj, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs );

// --- scalmr ---

void bl1_sscalmr(  uplo1_t uplo, int m, int n, float*    alpha, float*    a, int a_rs, int a_cs );
void bl1_dscalmr(  uplo1_t uplo, int m, int n, double*   alpha, double*   a, int a_rs, int a_cs );
void bl1_csscalmr( uplo1_t uplo, int m, int n, float*    alpha, scomplex* a, int a_rs, int a_cs );
void bl1_cscalmr(  uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs );
void bl1_zdscalmr( uplo1_t uplo, int m, int n, double*   alpha, dcomplex* a, int a_rs, int a_cs );
void bl1_zscalmr(  uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs );

// --- swap ---

void bl1_sswap( int n, float*    x, int incx, float*    y, int incy );
void bl1_dswap( int n, double*   x, int incx, double*   y, int incy );
void bl1_cswap( int n, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zswap( int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- swapv ---

void bl1_sswapv( int n, float*    x, int incx, float*    y, int incy );
void bl1_dswapv( int n, double*   x, int incx, double*   y, int incy );
void bl1_cswapv( int n, scomplex* x, int incx, scomplex* y, int incy );
void bl1_zswapv( int n, dcomplex* x, int incx, dcomplex* y, int incy );

// --- swapmt ---

void bl1_sswapmt( trans1_t trans, int m, int n, float*    a, int a_rs, int a_cs, float*    b, int b_rs, int b_cs );
void bl1_dswapmt( trans1_t trans, int m, int n, double*   a, int a_rs, int a_cs, double*   b, int b_rs, int b_cs );
void bl1_cswapmt( trans1_t trans, int m, int n, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs );
void bl1_zswapmt( trans1_t trans, int m, int n, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs );

