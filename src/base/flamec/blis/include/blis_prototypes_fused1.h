/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Fused Level-1 BLAS-like prototypes --------------------------------------

// --- axmyv2 ---

void bl1_saxmyv2( conj1_t conjx, int n, float*    alpha, float*    beta, float*    x, int inc_x, float*    y, int inc_y, float*    z, int inc_z );
void bl1_daxmyv2( conj1_t conjx, int n, double*   alpha, double*   beta, double*   x, int inc_x, double*   y, int inc_y, double*   z, int inc_z );
void bl1_caxmyv2( conj1_t conjx, int n, scomplex* alpha, scomplex* beta, scomplex* x, int inc_x, scomplex* y, int inc_y, scomplex* z, int inc_z );
void bl1_zaxmyv2( conj1_t conjx, int n, dcomplex* alpha, dcomplex* beta, dcomplex* x, int inc_x, dcomplex* y, int inc_y, dcomplex* z, int inc_z );

// --- axpyv2b ---

void bl1_saxpyv2b( int n, float*    beta1, float*    beta2, float*    a1, int inc_a1, float*    a2, int inc_a2, float*    w, int inc_w );
void bl1_daxpyv2b( int n, double*   beta1, double*   beta2, double*   a1, int inc_a1, double*   a2, int inc_a2, double*   w, int inc_w );
void bl1_caxpyv2b( int n, scomplex* beta1, scomplex* beta2, scomplex* a1, int inc_a1, scomplex* a2, int inc_a2, scomplex* w, int inc_w );
void bl1_zaxpyv2b( int n, dcomplex* beta1, dcomplex* beta2, dcomplex* a1, int inc_a1, dcomplex* a2, int inc_a2, dcomplex* w, int inc_w );

// --- axpyv3b ---

void bl1_saxpyv3b( int n, float*    beta1, float*    beta2, float*    beta3, float*    a1, int inc_a1, float*    a2, int inc_a2, float*    a3, int inc_a3, float*    w, int inc_w );
void bl1_daxpyv3b( int n, double*   beta1, double*   beta2, double*   beta3, double*   a1, int inc_a1, double*   a2, int inc_a2, double*   a3, int inc_a3, double*   w, int inc_w );
void bl1_caxpyv3b( int n, scomplex* beta1, scomplex* beta2, scomplex* beta3, scomplex* a1, int inc_a1, scomplex* a2, int inc_a2, scomplex* a3, int inc_a3, scomplex* w, int inc_w );
void bl1_zaxpyv3b( int n, dcomplex* beta1, dcomplex* beta2, dcomplex* beta3, dcomplex* a1, int inc_a1, dcomplex* a2, int inc_a2, dcomplex* a3, int inc_a3, dcomplex* w, int inc_w );

// --- axpyv2bdotaxpy ---

void bl1_saxpyv2bdotaxpy( int n, float*    beta, float*    u, int inc_u, float*    gamma, float*    z, int inc_z, float*    a, int inc_a, float*    x, int inc_x, float*    kappa, float*    rho, float*    w, int inc_w );
void bl1_daxpyv2bdotaxpy( int n, double*   beta, double*   u, int inc_u, double*   gamma, double*   z, int inc_z, double*   a, int inc_a, double*   x, int inc_x, double*   kappa, double*   rho, double*   w, int inc_w );
void bl1_caxpyv2bdotaxpy( int n, scomplex* beta, scomplex* u, int inc_u, scomplex* gamma, scomplex* z, int inc_z, scomplex* a, int inc_a, scomplex* x, int inc_x, scomplex* kappa, scomplex* rho, scomplex* w, int inc_w );
void bl1_zaxpyv2bdotaxpy( int n, dcomplex* beta, dcomplex* u, int inc_u, dcomplex* gamma, dcomplex* z, int inc_z, dcomplex* a, int inc_a, dcomplex* x, int inc_x, dcomplex* kappa, dcomplex* rho, dcomplex* w, int inc_w );

// --- dotsv2 ---

void bl1_sdotsv2( conj1_t conjxy, int n, float*    x, int inc_x, float*    y, int inc_y, float*    z, int inc_z, float*    beta, float*    rho_xz, float*    rho_yz );
void bl1_ddotsv2( conj1_t conjxy, int n, double*   x, int inc_x, double*   y, int inc_y, double*   z, int inc_z, double*   beta, double*   rho_xz, double*   rho_yz );
void bl1_cdotsv2( conj1_t conjxy, int n, scomplex* x, int inc_x, scomplex* y, int inc_y, scomplex* z, int inc_z, scomplex* beta, scomplex* rho_xz, scomplex* rho_yz );
void bl1_zdotsv2( conj1_t conjxy, int n, dcomplex* x, int inc_x, dcomplex* y, int inc_y, dcomplex* z, int inc_z, dcomplex* beta, dcomplex* rho_xz, dcomplex* rho_yz );

// --- dotsv3 ---

void bl1_sdotsv3( conj1_t conjxyw, int n, float*    x, int inc_x, float*    y, int inc_y, float*    w, int inc_w, float*    z, int inc_z, float*    beta, float*    rho_xz, float*    rho_yz, float*    rho_wz );
void bl1_ddotsv3( conj1_t conjxyw, int n, double*   x, int inc_x, double*   y, int inc_y, double*   w, int inc_w, double*   z, int inc_z, double*   beta, double*   rho_xz, double*   rho_yz, double*   rho_wz );
void bl1_cdotsv3( conj1_t conjxyw, int n, scomplex* x, int inc_x, scomplex* y, int inc_y, scomplex* w, int inc_w, scomplex* z, int inc_z, scomplex* beta, scomplex* rho_xz, scomplex* rho_yz, scomplex* rho_wz );
void bl1_zdotsv3( conj1_t conjxyw, int n, dcomplex* x, int inc_x, dcomplex* y, int inc_y, dcomplex* w, int inc_w, dcomplex* z, int inc_z, dcomplex* beta, dcomplex* rho_xz, dcomplex* rho_yz, dcomplex* rho_wz );

// --- dotaxpy ---

void bl1_sdotaxpy( int n, float*    a, int inc_a, float*    x, int inc_x, float*    kappa, float*    rho, float*    w, int inc_w );
void bl1_ddotaxpy( int n, double*   a, int inc_a, double*   x, int inc_x, double*   kappa, double*   rho, double*   w, int inc_w );
void bl1_cdotaxpy( int n, scomplex* a, int inc_a, scomplex* x, int inc_x, scomplex* kappa, scomplex* rho, scomplex* w, int inc_w );
void bl1_zdotaxpy( int n, dcomplex* a, int inc_a, dcomplex* x, int inc_x, dcomplex* kappa, dcomplex* rho, dcomplex* w, int inc_w );

// --- dotaxmyv2 ---

void bl1_sdotaxmyv2( int n, float*    alpha, float*    beta, float*    x, int inc_x, float*    u, int inc_u, float*    rho, float*    y, int inc_y, float*    z, int inc_z );
void bl1_ddotaxmyv2( int n, double*   alpha, double*   beta, double*   x, int inc_x, double*   u, int inc_u, double*   rho, double*   y, int inc_y, double*   z, int inc_z );
void bl1_cdotaxmyv2( int n, scomplex* alpha, scomplex* beta, scomplex* x, int inc_x, scomplex* u, int inc_u, scomplex* rho, scomplex* y, int inc_y, scomplex* z, int inc_z );
void bl1_zdotaxmyv2( int n, dcomplex* alpha, dcomplex* beta, dcomplex* x, int inc_x, dcomplex* u, int inc_u, dcomplex* rho, dcomplex* y, int inc_y, dcomplex* z, int inc_z );

// --- dotv2axpyv2b ---

void bl1_sdotv2axpyv2b( int n, float*    a1, int inc_a1, float*    a2, int inc_a2, float*    x,  int inc_x, float*    kappa1, float*    kappa2, float*    rho1, float*    rho2, float*    w, int inc_w );
void bl1_ddotv2axpyv2b( int n, double*   a1, int inc_a1, double*   a2, int inc_a2, double*   x,  int inc_x, double*   kappa1, double*   kappa2, double*   rho1, double*   rho2, double*   w, int inc_w );
void bl1_cdotv2axpyv2b( int n, scomplex* a1, int inc_a1, scomplex* a2, int inc_a2, scomplex* x,  int inc_x, scomplex* kappa1, scomplex* kappa2, scomplex* rho1, scomplex* rho2, scomplex* w, int inc_w );
void bl1_zdotv2axpyv2b( int n, dcomplex* a1, int inc_a1, dcomplex* a2, int inc_a2, dcomplex* x,  int inc_x, dcomplex* kappa1, dcomplex* kappa2, dcomplex* rho1, dcomplex* rho2, dcomplex* w, int inc_w );

// --- axpyv2bdots ---

void bl1_zaxpyv2bdots( int       n,
                       dcomplex* alpha1,
                       dcomplex* alpha2,
                       dcomplex* x1, int inc_x1,
                       dcomplex* x2, int inc_x2,
                       dcomplex* y,  int inc_y,
                       dcomplex* u,  int inc_u,
                       dcomplex* beta,
                       dcomplex* rho );
