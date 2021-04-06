/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Level-1 BLAS-like prototypes --------------------------------------------

// --- amax ---

void bl1_samax( integer n, float*    x, integer incx, integer* index );
void bl1_damax( integer n, double*   x, integer incx, integer* index );
void bl1_camax( integer n, scomplex* x, integer incx, integer* index );
void bl1_zamax( integer n, dcomplex* x, integer incx, integer* index );

// --- asum ---

void bl1_sasum( integer n, float*    x, integer incx, float*  norm );
void bl1_dasum( integer n, double*   x, integer incx, double* norm );
void bl1_casum( integer n, scomplex* x, integer incx, float*  norm );
void bl1_zasum( integer n, dcomplex* x, integer incx, double* norm );

// --- axpy ---

void bl1_saxpy( integer n, float*    alpha, float*    x, integer incx, float*    y, integer incy );
void bl1_daxpy( integer n, double*   alpha, double*   x, integer incx, double*   y, integer incy );
void bl1_caxpy( integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy );
void bl1_zaxpy( integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy );

// --- axpyv ---

void bl1_saxpyv( conj1_t conj, integer n, float*    alpha, float*    x, integer incx, float*    y, integer incy );
void bl1_daxpyv( conj1_t conj, integer n, double*   alpha, double*   x, integer incx, double*   y, integer incy );
void bl1_caxpyv( conj1_t conj, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy );
void bl1_zaxpyv( conj1_t conj, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy );

// --- axpymt ---

void bl1_saxpymt( trans1_t trans, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_daxpymt( trans1_t trans, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_caxpymt( trans1_t trans, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zaxpymt( trans1_t trans, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

// --- axpymrt ---

void bl1_saxpymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_daxpymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_caxpymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zaxpymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

// --- axpysv ---

void bl1_saxpysv( integer n, float*    alpha0, float*    alpha1, float*    x, integer incx, float*    beta, float*    y, integer incy );
void bl1_daxpysv( integer n, double*   alpha0, double*   alpha1, double*   x, integer incx, double*   beta, double*   y, integer incy );
void bl1_caxpysv( integer n, scomplex* alpha0, scomplex* alpha1, scomplex* x, integer incx, scomplex* beta, scomplex* y, integer incy );
void bl1_zaxpysv( integer n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* x, integer incx, dcomplex* beta, dcomplex* y, integer incy );

// --- axpysmt ---

void bl1_saxpysmt( trans1_t trans, integer m, integer n, float*    alpha0, float*    alpha1, float*    a, integer a_rs, integer a_cs, float*    beta, float*    b, integer b_rs, integer b_cs );
void bl1_daxpysmt( trans1_t trans, integer m, integer n, double*   alpha0, double*   alpha1, double*   a, integer a_rs, integer a_cs, double*   beta, double*   b, integer b_rs, integer b_cs );
void bl1_caxpysmt( trans1_t trans, integer m, integer n, scomplex* alpha0, scomplex* alpha1, scomplex* a, integer a_rs, integer a_cs, scomplex* beta, scomplex* b, integer b_rs, integer b_cs );
void bl1_zaxpysmt( trans1_t trans, integer m, integer n, dcomplex* alpha0, dcomplex* alpha1, dcomplex* a, integer a_rs, integer a_cs, dcomplex* beta, dcomplex* b, integer b_rs, integer b_cs );

// --- conjv ---

void bl1_sconjv( integer m, float* x, integer incx );
void bl1_dconjv( integer m, double* x, integer incx );
void bl1_cconjv( integer m, scomplex* x, integer incx );
void bl1_zconjv( integer m, dcomplex* x, integer incx );

// --- conjm ---

void bl1_sconjm( integer m, integer n, float*    a, integer a_rs, integer a_cs );
void bl1_dconjm( integer m, integer n, double*   a, integer a_rs, integer a_cs );
void bl1_cconjm( integer m, integer n, scomplex* a, integer a_rs, integer a_cs );
void bl1_zconjm( integer m, integer n, dcomplex* a, integer a_rs, integer a_cs );

// --- conjmr ---

void bl1_sconjmr( uplo1_t uplo, integer m, integer n, float*    a, integer a_rs, integer a_cs );
void bl1_dconjmr( uplo1_t uplo, integer m, integer n, double*   a, integer a_rs, integer a_cs );
void bl1_cconjmr( uplo1_t uplo, integer m, integer n, scomplex* a, integer a_rs, integer a_cs );
void bl1_zconjmr( uplo1_t uplo, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs );

// --- copy ---

void bl1_scopy( integer m, float*    x, integer incx, float*    y, integer incy );
void bl1_dcopy( integer m, double*   x, integer incx, double*   y, integer incy );
void bl1_ccopy( integer m, scomplex* x, integer incx, scomplex* y, integer incy );
void bl1_zcopy( integer m, dcomplex* x, integer incx, dcomplex* y, integer incy );

// --- copyv ---

void bl1_icopyv( conj1_t conj, integer m, integer*      x, integer incx, integer*      y, integer incy );
void bl1_scopyv( conj1_t conj, integer m, float*    x, integer incx, float*    y, integer incy );
void bl1_dcopyv( conj1_t conj, integer m, double*   x, integer incx, double*   y, integer incy );
void bl1_ccopyv( conj1_t conj, integer m, scomplex* x, integer incx, scomplex* y, integer incy );
void bl1_zcopyv( conj1_t conj, integer m, dcomplex* x, integer incx, dcomplex* y, integer incy );

void bl1_sdcopyv( conj1_t conj, integer m, float*    x, integer incx, double*   y, integer incy );
void bl1_dscopyv( conj1_t conj, integer m, double*   x, integer incx, float*    y, integer incy );
void bl1_sccopyv( conj1_t conj, integer m, float*    x, integer incx, scomplex* y, integer incy );
void bl1_cscopyv( conj1_t conj, integer m, scomplex* x, integer incx, float*    y, integer incy );
void bl1_szcopyv( conj1_t conj, integer m, float*    x, integer incx, dcomplex* y, integer incy );
void bl1_zscopyv( conj1_t conj, integer m, dcomplex* x, integer incx, float*    y, integer incy );
void bl1_dccopyv( conj1_t conj, integer m, double*   x, integer incx, scomplex* y, integer incy );
void bl1_cdcopyv( conj1_t conj, integer m, scomplex* x, integer incx, double*   y, integer incy );
void bl1_dzcopyv( conj1_t conj, integer m, double*   x, integer incx, dcomplex* y, integer incy );
void bl1_zdcopyv( conj1_t conj, integer m, dcomplex* x, integer incx, double*   y, integer incy );
void bl1_czcopyv( conj1_t conj, integer m, scomplex* x, integer incx, dcomplex* y, integer incy );
void bl1_zccopyv( conj1_t conj, integer m, dcomplex* x, integer incx, scomplex* y, integer incy );

// --- copymr ---

void bl1_scopymr( uplo1_t uplo, integer m, integer n, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_dcopymr( uplo1_t uplo, integer m, integer n, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_ccopymr( uplo1_t uplo, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zcopymr( uplo1_t uplo, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

void bl1_sscopymr( uplo1_t uplo, integer m, integer n, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_sdcopymr( uplo1_t uplo, integer m, integer n, float*    a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_dscopymr( uplo1_t uplo, integer m, integer n, double*   a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_sccopymr( uplo1_t uplo, integer m, integer n, float*    a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_cscopymr( uplo1_t uplo, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_szcopymr( uplo1_t uplo, integer m, integer n, float*    a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_zscopymr( uplo1_t uplo, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_ddcopymr( uplo1_t uplo, integer m, integer n, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_dccopymr( uplo1_t uplo, integer m, integer n, double*   a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_cdcopymr( uplo1_t uplo, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_dzcopymr( uplo1_t uplo, integer m, integer n, double*   a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_zdcopymr( uplo1_t uplo, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_cccopymr( uplo1_t uplo, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_czcopymr( uplo1_t uplo, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_zccopymr( uplo1_t uplo, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zzcopymr( uplo1_t uplo, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

// --- copymrt ---

void bl1_scopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_dcopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_ccopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zcopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

void bl1_sscopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_sdcopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_sccopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_szcopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_dscopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_ddcopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_dccopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_dzcopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_cscopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_cdcopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_cccopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_czcopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_zscopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_zdcopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_zccopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zzcopymrt( uplo1_t uplo, trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

// --- copymt ---

void bl1_icopymt( trans1_t trans, integer m, integer n, integer*      a, integer a_rs, integer a_cs, integer*      b, integer b_rs, integer b_cs );
void bl1_scopymt( trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_dcopymt( trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_ccopymt( trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zcopymt( trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

void bl1_sscopymt( trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_sdcopymt( trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_dscopymt( trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_sccopymt( trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_cscopymt( trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_szcopymt( trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_zscopymt( trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_ddcopymt( trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_dccopymt( trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_cdcopymt( trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_dzcopymt( trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_zdcopymt( trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_cccopymt( trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_czcopymt( trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );
void bl1_zccopymt( trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zzcopymt( trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

// --- dot ---

void bl1_cdot_in( conj1_t conj, integer n, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* rho );
void bl1_zdot_in( conj1_t conj, integer n, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* rho );

void bl1_sdot( conj1_t conj, integer n, float*    x, integer incx, float*    y, integer incy, float*    rho );
void bl1_ddot( conj1_t conj, integer n, double*   x, integer incx, double*   y, integer incy, double*   rho );
void bl1_cdot( conj1_t conj, integer n, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* rho );
void bl1_zdot( conj1_t conj, integer n, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* rho );

// --- dots ---

void bl1_sdots( conj1_t conj, integer n, float*    alpha, float*    x, integer incx, float*    y, integer incy, float*    beta, float*    rho );
void bl1_ddots( conj1_t conj, integer n, double*   alpha, double*   x, integer incx, double*   y, integer incy, double*   beta, double*   rho );
void bl1_cdots( conj1_t conj, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* beta, scomplex* rho );
void bl1_zdots( conj1_t conj, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* beta, dcomplex* rho );

// --- dot2s ---

void bl1_sdot2s( conj1_t conj, integer n, float*    alpha, float*    x, integer incx, float*    y, integer incy, float*    beta, float*    rho );
void bl1_ddot2s( conj1_t conj, integer n, double*   alpha, double*   x, integer incx, double*   y, integer incy, double*   beta, double*   rho );
void bl1_cdot2s( conj1_t conj, integer n, scomplex* alpha, scomplex* x, integer incx, scomplex* y, integer incy, scomplex* beta, scomplex* rho );
void bl1_zdot2s( conj1_t conj, integer n, dcomplex* alpha, dcomplex* x, integer incx, dcomplex* y, integer incy, dcomplex* beta, dcomplex* rho );

// --- fnorm ---

void bl1_sfnorm( integer m, integer n, float*    a, integer a_rs, integer a_cs, float*  norm );
void bl1_dfnorm( integer m, integer n, double*   a, integer a_rs, integer a_cs, double* norm );
void bl1_cfnorm( integer m, integer n, scomplex* a, integer a_rs, integer a_cs, float*  norm );
void bl1_zfnorm( integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, double* norm );

// --- invscalv ---

void bl1_sinvscalv(  conj1_t conj, integer n, float*    alpha, float*    x, integer incx );
void bl1_dinvscalv(  conj1_t conj, integer n, double*   alpha, double*   x, integer incx );
void bl1_csinvscalv( conj1_t conj, integer n, float*    alpha, scomplex* x, integer incx );
void bl1_cinvscalv(  conj1_t conj, integer n, scomplex* alpha, scomplex* x, integer incx );
void bl1_zdinvscalv( conj1_t conj, integer n, double*   alpha, dcomplex* x, integer incx );
void bl1_zinvscalv(  conj1_t conj, integer n, dcomplex* alpha, dcomplex* x, integer incx );

// --- invscalm ---

void bl1_sinvscalm(  conj1_t conj, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs );
void bl1_dinvscalm(  conj1_t conj, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs );
void bl1_csinvscalm( conj1_t conj, integer m, integer n, float*    alpha, scomplex* a, integer a_rs, integer a_cs );
void bl1_cinvscalm(  conj1_t conj, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs );
void bl1_zdinvscalm( conj1_t conj, integer m, integer n, double*   alpha, dcomplex* a, integer a_rs, integer a_cs );
void bl1_zinvscalm(  conj1_t conj, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs );

// --- nrm2 ---

void bl1_snrm2( integer n, float*    x, integer incx, float*  norm );
void bl1_dnrm2( integer n, double*   x, integer incx, double* norm );
void bl1_cnrm2( integer n, scomplex* x, integer incx, float*  norm );
void bl1_znrm2( integer n, dcomplex* x, integer incx, double* norm );

// --- scal ---

void bl1_sscal(  integer n, float*    alpha, float*    x, integer incx );
void bl1_dscal(  integer n, double*   alpha, double*   x, integer incx );
void bl1_csscal( integer n, float*    alpha, scomplex* x, integer incx );
void bl1_cscal(  integer n, scomplex* alpha, scomplex* x, integer incx );
void bl1_zdscal( integer n, double*   alpha, dcomplex* x, integer incx );
void bl1_zscal(  integer n, dcomplex* alpha, dcomplex* x, integer incx );

// --- scalv ---

void bl1_sscalv(  conj1_t conj, integer n, float*    alpha, float*    x, integer incx );
void bl1_dscalv(  conj1_t conj, integer n, double*   alpha, double*   x, integer incx );
void bl1_csscalv( conj1_t conj, integer n, float*    alpha, scomplex* x, integer incx );
void bl1_cscalv(  conj1_t conj, integer n, scomplex* alpha, scomplex* x, integer incx );
void bl1_zdscalv( conj1_t conj, integer n, double*   alpha, dcomplex* x, integer incx );
void bl1_zscalv(  conj1_t conj, integer n, dcomplex* alpha, dcomplex* x, integer incx );

// --- scalm ---

void bl1_sscalm(  conj1_t conj, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs );
void bl1_dscalm(  conj1_t conj, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs );
void bl1_csscalm( conj1_t conj, integer m, integer n, float*    alpha, scomplex* a, integer a_rs, integer a_cs );
void bl1_cscalm(  conj1_t conj, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs );
void bl1_zdscalm( conj1_t conj, integer m, integer n, double*   alpha, dcomplex* a, integer a_rs, integer a_cs );
void bl1_zscalm(  conj1_t conj, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs );

// --- scalmr ---

void bl1_sscalmr(  uplo1_t uplo, integer m, integer n, float*    alpha, float*    a, integer a_rs, integer a_cs );
void bl1_dscalmr(  uplo1_t uplo, integer m, integer n, double*   alpha, double*   a, integer a_rs, integer a_cs );
void bl1_csscalmr( uplo1_t uplo, integer m, integer n, float*    alpha, scomplex* a, integer a_rs, integer a_cs );
void bl1_cscalmr(  uplo1_t uplo, integer m, integer n, scomplex* alpha, scomplex* a, integer a_rs, integer a_cs );
void bl1_zdscalmr( uplo1_t uplo, integer m, integer n, double*   alpha, dcomplex* a, integer a_rs, integer a_cs );
void bl1_zscalmr(  uplo1_t uplo, integer m, integer n, dcomplex* alpha, dcomplex* a, integer a_rs, integer a_cs );

// --- swap ---

void bl1_sswap( integer n, float*    x, integer incx, float*    y, integer incy );
void bl1_dswap( integer n, double*   x, integer incx, double*   y, integer incy );
void bl1_cswap( integer n, scomplex* x, integer incx, scomplex* y, integer incy );
void bl1_zswap( integer n, dcomplex* x, integer incx, dcomplex* y, integer incy );

// --- swapv ---

void bl1_sswapv( integer n, float*    x, integer incx, float*    y, integer incy );
void bl1_dswapv( integer n, double*   x, integer incx, double*   y, integer incy );
void bl1_cswapv( integer n, scomplex* x, integer incx, scomplex* y, integer incy );
void bl1_zswapv( integer n, dcomplex* x, integer incx, dcomplex* y, integer incy );

// --- swapmt ---

void bl1_sswapmt( trans1_t trans, integer m, integer n, float*    a, integer a_rs, integer a_cs, float*    b, integer b_rs, integer b_cs );
void bl1_dswapmt( trans1_t trans, integer m, integer n, double*   a, integer a_rs, integer a_cs, double*   b, integer b_rs, integer b_cs );
void bl1_cswapmt( trans1_t trans, integer m, integer n, scomplex* a, integer a_rs, integer a_cs, scomplex* b, integer b_rs, integer b_cs );
void bl1_zswapmt( trans1_t trans, integer m, integer n, dcomplex* a, integer a_rs, integer a_cs, dcomplex* b, integer b_rs, integer b_cs );

