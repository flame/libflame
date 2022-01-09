/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Fused Level-1 BLAS-like prototypes --------------------------------------

// --- axmyv2 ---

void bl1_saxmyv2( conj1_t conjx, integer n, float*    alpha, float*    beta, float*    x, integer inc_x, float*    y, integer inc_y, float*    z, integer inc_z );
void bl1_daxmyv2( conj1_t conjx, integer n, double*   alpha, double*   beta, double*   x, integer inc_x, double*   y, integer inc_y, double*   z, integer inc_z );
void bl1_caxmyv2( conj1_t conjx, integer n, scomplex* alpha, scomplex* beta, scomplex* x, integer inc_x, scomplex* y, integer inc_y, scomplex* z, integer inc_z );
void bl1_zaxmyv2( conj1_t conjx, integer n, dcomplex* alpha, dcomplex* beta, dcomplex* x, integer inc_x, dcomplex* y, integer inc_y, dcomplex* z, integer inc_z );

// --- axpyv2b ---

void bl1_saxpyv2b( integer n, float*    beta1, float*    beta2, float*    a1, integer inc_a1, float*    a2, integer inc_a2, float*    w, integer inc_w );
void bl1_daxpyv2b( integer n, double*   beta1, double*   beta2, double*   a1, integer inc_a1, double*   a2, integer inc_a2, double*   w, integer inc_w );
void bl1_caxpyv2b( integer n, scomplex* beta1, scomplex* beta2, scomplex* a1, integer inc_a1, scomplex* a2, integer inc_a2, scomplex* w, integer inc_w );
void bl1_zaxpyv2b( integer n, dcomplex* beta1, dcomplex* beta2, dcomplex* a1, integer inc_a1, dcomplex* a2, integer inc_a2, dcomplex* w, integer inc_w );

// --- axpyv3b ---

void bl1_saxpyv3b( integer n, float*    beta1, float*    beta2, float*    beta3, float*    a1, integer inc_a1, float*    a2, integer inc_a2, float*    a3, integer inc_a3, float*    w, integer inc_w );
void bl1_daxpyv3b( integer n, double*   beta1, double*   beta2, double*   beta3, double*   a1, integer inc_a1, double*   a2, integer inc_a2, double*   a3, integer inc_a3, double*   w, integer inc_w );
void bl1_caxpyv3b( integer n, scomplex* beta1, scomplex* beta2, scomplex* beta3, scomplex* a1, integer inc_a1, scomplex* a2, integer inc_a2, scomplex* a3, integer inc_a3, scomplex* w, integer inc_w );
void bl1_zaxpyv3b( integer n, dcomplex* beta1, dcomplex* beta2, dcomplex* beta3, dcomplex* a1, integer inc_a1, dcomplex* a2, integer inc_a2, dcomplex* a3, integer inc_a3, dcomplex* w, integer inc_w );

// --- axpyv2bdotaxpy ---

void bl1_saxpyv2bdotaxpy( integer n, float*    beta, float*    u, integer inc_u, float*    gamma, float*    z, integer inc_z, float*    a, integer inc_a, float*    x, integer inc_x, float*    kappa, float*    rho, float*    w, integer inc_w );
void bl1_daxpyv2bdotaxpy( integer n, double*   beta, double*   u, integer inc_u, double*   gamma, double*   z, integer inc_z, double*   a, integer inc_a, double*   x, integer inc_x, double*   kappa, double*   rho, double*   w, integer inc_w );
void bl1_caxpyv2bdotaxpy( integer n, scomplex* beta, scomplex* u, integer inc_u, scomplex* gamma, scomplex* z, integer inc_z, scomplex* a, integer inc_a, scomplex* x, integer inc_x, scomplex* kappa, scomplex* rho, scomplex* w, integer inc_w );
void bl1_zaxpyv2bdotaxpy( integer n, dcomplex* beta, dcomplex* u, integer inc_u, dcomplex* gamma, dcomplex* z, integer inc_z, dcomplex* a, integer inc_a, dcomplex* x, integer inc_x, dcomplex* kappa, dcomplex* rho, dcomplex* w, integer inc_w );

// --- dotsv2 ---

void bl1_sdotsv2( conj1_t conjxy, integer n, float*    x, integer inc_x, float*    y, integer inc_y, float*    z, integer inc_z, float*    beta, float*    rho_xz, float*    rho_yz );
void bl1_ddotsv2( conj1_t conjxy, integer n, double*   x, integer inc_x, double*   y, integer inc_y, double*   z, integer inc_z, double*   beta, double*   rho_xz, double*   rho_yz );
void bl1_cdotsv2( conj1_t conjxy, integer n, scomplex* x, integer inc_x, scomplex* y, integer inc_y, scomplex* z, integer inc_z, scomplex* beta, scomplex* rho_xz, scomplex* rho_yz );
void bl1_zdotsv2( conj1_t conjxy, integer n, dcomplex* x, integer inc_x, dcomplex* y, integer inc_y, dcomplex* z, integer inc_z, dcomplex* beta, dcomplex* rho_xz, dcomplex* rho_yz );

// --- dotsv3 ---

void bl1_sdotsv3( conj1_t conjxyw, integer n, float*    x, integer inc_x, float*    y, integer inc_y, float*    w, integer inc_w, float*    z, integer inc_z, float*    beta, float*    rho_xz, float*    rho_yz, float*    rho_wz );
void bl1_ddotsv3( conj1_t conjxyw, integer n, double*   x, integer inc_x, double*   y, integer inc_y, double*   w, integer inc_w, double*   z, integer inc_z, double*   beta, double*   rho_xz, double*   rho_yz, double*   rho_wz );
void bl1_cdotsv3( conj1_t conjxyw, integer n, scomplex* x, integer inc_x, scomplex* y, integer inc_y, scomplex* w, integer inc_w, scomplex* z, integer inc_z, scomplex* beta, scomplex* rho_xz, scomplex* rho_yz, scomplex* rho_wz );
void bl1_zdotsv3( conj1_t conjxyw, integer n, dcomplex* x, integer inc_x, dcomplex* y, integer inc_y, dcomplex* w, integer inc_w, dcomplex* z, integer inc_z, dcomplex* beta, dcomplex* rho_xz, dcomplex* rho_yz, dcomplex* rho_wz );

// --- dotaxpy ---

void bl1_sdotaxpy( integer n, float*    a, integer inc_a, float*    x, integer inc_x, float*    kappa, float*    rho, float*    w, integer inc_w );
void bl1_ddotaxpy( integer n, double*   a, integer inc_a, double*   x, integer inc_x, double*   kappa, double*   rho, double*   w, integer inc_w );
void bl1_cdotaxpy( integer n, scomplex* a, integer inc_a, scomplex* x, integer inc_x, scomplex* kappa, scomplex* rho, scomplex* w, integer inc_w );
void bl1_zdotaxpy( integer n, dcomplex* a, integer inc_a, dcomplex* x, integer inc_x, dcomplex* kappa, dcomplex* rho, dcomplex* w, integer inc_w );

// --- dotaxmyv2 ---

void bl1_sdotaxmyv2( integer n, float*    alpha, float*    beta, float*    x, integer inc_x, float*    u, integer inc_u, float*    rho, float*    y, integer inc_y, float*    z, integer inc_z );
void bl1_ddotaxmyv2( integer n, double*   alpha, double*   beta, double*   x, integer inc_x, double*   u, integer inc_u, double*   rho, double*   y, integer inc_y, double*   z, integer inc_z );
void bl1_cdotaxmyv2( integer n, scomplex* alpha, scomplex* beta, scomplex* x, integer inc_x, scomplex* u, integer inc_u, scomplex* rho, scomplex* y, integer inc_y, scomplex* z, integer inc_z );
void bl1_zdotaxmyv2( integer n, dcomplex* alpha, dcomplex* beta, dcomplex* x, integer inc_x, dcomplex* u, integer inc_u, dcomplex* rho, dcomplex* y, integer inc_y, dcomplex* z, integer inc_z );

// --- dotv2axpyv2b ---

void bl1_sdotv2axpyv2b( integer n, float*    a1, integer inc_a1, float*    a2, integer inc_a2, float*    x,  integer inc_x, float*    kappa1, float*    kappa2, float*    rho1, float*    rho2, float*    w, integer inc_w );
void bl1_ddotv2axpyv2b( integer n, double*   a1, integer inc_a1, double*   a2, integer inc_a2, double*   x,  integer inc_x, double*   kappa1, double*   kappa2, double*   rho1, double*   rho2, double*   w, integer inc_w );
void bl1_cdotv2axpyv2b( integer n, scomplex* a1, integer inc_a1, scomplex* a2, integer inc_a2, scomplex* x,  integer inc_x, scomplex* kappa1, scomplex* kappa2, scomplex* rho1, scomplex* rho2, scomplex* w, integer inc_w );
void bl1_zdotv2axpyv2b( integer n, dcomplex* a1, integer inc_a1, dcomplex* a2, integer inc_a2, dcomplex* x,  integer inc_x, dcomplex* kappa1, dcomplex* kappa2, dcomplex* rho1, dcomplex* rho2, dcomplex* w, integer inc_w );

// --- axpyv2bdots ---

void bl1_zaxpyv2bdots( integer       n,
                       dcomplex* alpha1,
                       dcomplex* alpha2,
                       dcomplex* x1, integer inc_x1,
                       dcomplex* x2, integer inc_x2,
                       dcomplex* y,  integer inc_y,
                       dcomplex* u,  integer inc_u,
                       dcomplex* beta,
                       dcomplex* rho );
