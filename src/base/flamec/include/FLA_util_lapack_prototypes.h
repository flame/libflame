/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- LAPACK-related utility prototypes ---------------------------------------

FLA_Error FLA_Househ2_UT( FLA_Side side, FLA_Obj chi_1, FLA_Obj x2, FLA_Obj tau );
FLA_Error FLA_Househ2_UT_l_ops( int m_x2,
                                float* chi_1,
                                float* x2, int inc_x2,
                                float* tau );
FLA_Error FLA_Househ2_UT_l_opd( int m_x2,
                                double* chi_1,
                                double* x2, int inc_x2,
                                double* tau );
FLA_Error FLA_Househ2_UT_l_opc( int m_x2,
                                scomplex* chi_1,
                                scomplex* x2, int inc_x2,
                                scomplex* tau );
FLA_Error FLA_Househ2_UT_l_opz( int m_x2,
                                dcomplex* chi_1,
                                dcomplex* x2, int inc_x2,
                                dcomplex* tau );
FLA_Error FLA_Househ2_UT_r_ops( int m_x2,
                                float* chi_1,
                                float* x2, int inc_x2,
                                float* tau );
FLA_Error FLA_Househ2_UT_r_opd( int m_x2,
                                double* chi_1,
                                double* x2, int inc_x2,
                                double* tau );
FLA_Error FLA_Househ2_UT_r_opc( int m_x2,
                                scomplex* chi_1,
                                scomplex* x2, int inc_x2,
                                scomplex* tau );
FLA_Error FLA_Househ2_UT_r_opz( int m_x2,
                                dcomplex* chi_1,
                                dcomplex* x2, int inc_x2,
                                dcomplex* tau );

FLA_Error FLA_Househ3UD_UT( FLA_Obj chi_1, FLA_Obj x2, FLA_Obj y2, FLA_Obj tau );
FLA_Error FLA_Househ3UD_UT_ops( int m_x2,
                                int m_y2,
                                float* chi_1,
                                float* x2, int inc_x2,
                                float* y2, int inc_y2,
                                float* tau );
FLA_Error FLA_Househ3UD_UT_opd( int m_x2,
                                int m_y2,
                                double* chi_1,
                                double* x2, int inc_x2,
                                double* y2, int inc_y2,
                                double* tau );
FLA_Error FLA_Househ3UD_UT_opc( int m_x2,
                                int m_y2,
                                scomplex* chi_1,
                                scomplex* x2, int inc_x2,
                                scomplex* y2, int inc_y2,
                                scomplex* tau );
FLA_Error FLA_Househ3UD_UT_opz( int m_x2,
                                int m_y2,
                                dcomplex* chi_1,
                                dcomplex* x2, int inc_x2,
                                dcomplex* y2, int inc_y2,
                                dcomplex* tau );

FLA_Error FLA_Househ2s_UT( FLA_Side side, FLA_Obj chi_1, FLA_Obj x2, FLA_Obj alpha, FLA_Obj chi_1_minus_alpha, FLA_Obj tau );
FLA_Error FLA_Househ2s_UT_l_ops( int    m_x2,
                                 float* chi_1,
                                 float* x2, int inc_x2,
                                 float* alpha,
                                 float* chi_1_minus_alpha,
                                 float* tau );
FLA_Error FLA_Househ2s_UT_l_opd( int     m_x2,
                                 double* chi_1,
                                 double* x2, int inc_x2,
                                 double* alpha,
                                 double* chi_1_minus_alpha,
                                 double* tau );
FLA_Error FLA_Househ2s_UT_l_opc( int       m_x2,
                                 scomplex* chi_1,
                                 scomplex* x2, int inc_x2,
                                 scomplex* alpha,
                                 scomplex* chi_1_minus_alpha,
                                 scomplex* tau );
FLA_Error FLA_Househ2s_UT_l_opz( int       m_x2,
                                 dcomplex* chi_1,
                                 dcomplex* x2, int inc_x2,
                                 dcomplex* alpha,
                                 dcomplex* chi_1_minus_alpha,
                                 dcomplex* tau );
FLA_Error FLA_Househ2s_UT_r_ops( int    m_x2,
                                 float* chi_1,
                                 float* x2, int inc_x2,
                                 float* alpha,
                                 float* chi_1_minus_alpha,
                                 float* tau );
FLA_Error FLA_Househ2s_UT_r_opd( int     m_x2,
                                 double* chi_1,
                                 double* x2, int inc_x2,
                                 double* alpha,
                                 double* chi_1_minus_alpha,
                                 double* tau );
FLA_Error FLA_Househ2s_UT_r_opc( int       m_x2,
                                 scomplex* chi_1,
                                 scomplex* x2, int inc_x2,
                                 scomplex* alpha,
                                 scomplex* chi_1_minus_alpha,
                                 scomplex* tau );
FLA_Error FLA_Househ2s_UT_r_opz( int       m_x2,
                                 dcomplex* chi_1,
                                 dcomplex* x2, int inc_x2,
                                 dcomplex* alpha,
                                 dcomplex* chi_1_minus_alpha,
                                 dcomplex* tau );

FLA_Error FLA_Hev_2x2( FLA_Obj alpha11, FLA_Obj alpha21, FLA_Obj alpha22,
                       FLA_Obj lambda1, FLA_Obj lambda2 );
FLA_Error FLA_Hev_2x2_ops( float*    buff_alpha11,
                           float*    buff_alpha21,
                           float*    buff_alpha22,
                           float*    buff_lambda1,
                           float*    buff_lambda2 );
FLA_Error FLA_Hev_2x2_opd( double*   buff_alpha11,
                           double*   buff_alpha21,
                           double*   buff_alpha22,
                           double*   buff_lambda1,
                           double*   buff_lambda2 );

FLA_Error FLA_Hevv_2x2( FLA_Obj alpha11, FLA_Obj alpha21, FLA_Obj alpha22,
                        FLA_Obj lambda1, FLA_Obj lambda2,
                        FLA_Obj gamma1,  FLA_Obj sigma1 );
FLA_Error FLA_Hevv_2x2_ops( float*    alpha11,
                            float*    alpha21,
                            float*    alpha22,
                            float*    lambda1,
                            float*    lambda2,
                            float*    gamma1,
                            float*    sigma1 );
FLA_Error FLA_Hevv_2x2_opd( double*   alpha11,
                            double*   alpha21,
                            double*   alpha22,
                            double*   lambda1,
                            double*   lambda2,
                            double*   gamma1,
                            double*   sigma1 );
FLA_Error FLA_Hevv_2x2_opc( scomplex* alpha11,
                            scomplex* alpha21,
                            scomplex* alpha22,
                            float*    lambda1,
                            float*    lambda2,
                            float*    gamma1,
                            scomplex* sigma1 );
FLA_Error FLA_Hevv_2x2_opz( dcomplex* alpha11,
                            dcomplex* alpha21,
                            dcomplex* alpha22,
                            double*   lambda1,
                            double*   lambda2,
                            double*   gamma1,
                            dcomplex* sigma1 );

FLA_Error FLA_Wilkshift_tridiag( FLA_Obj delta1, FLA_Obj epsilon, FLA_Obj delta2, FLA_Obj kappa );
FLA_Error FLA_Wilkshift_tridiag_ops( float   delta1,
                                     float   epsilon,
                                     float   delta2,
                                     float*  kappa );
FLA_Error FLA_Wilkshift_tridiag_opd( double  delta1,
                                     double  epsilon,
                                     double  delta2,
                                     double* kappa );

FLA_Error FLA_Pythag2( FLA_Obj chi, FLA_Obj psi, FLA_Obj rho );
FLA_Error FLA_Pythag2_ops( float*    chi,
                           float*    psi,
                           float*    rho );
FLA_Error FLA_Pythag2_opd( double*   chi,
                           double*   psi,
                           double*   rho );

FLA_Error FLA_Pythag3( FLA_Obj chi, FLA_Obj psi, FLA_Obj zeta, FLA_Obj rho );
FLA_Error FLA_Pythag3_ops( float*    chi,
                           float*    psi,
                           float*    zeta,
                           float*    rho );
FLA_Error FLA_Pythag3_opd( double*   chi,
                           double*   psi,
                           double*   zeta,
                           double*   rho );

FLA_Error FLA_Sort_evd( FLA_Direct direct, FLA_Obj l, FLA_Obj V );
FLA_Error FLA_Sort_evd_f_ops( int       m_A,
                              float*    l, int inc_l,
                              float*    V, int rs_V, int cs_V );
FLA_Error FLA_Sort_evd_b_ops( int       m_A,
                              float*    l, int inc_l,
                              float*    V, int rs_V, int cs_V );
FLA_Error FLA_Sort_evd_f_opd( int       m_A,
                              double*   l, int inc_l,
                              double*   V, int rs_V, int cs_V );
FLA_Error FLA_Sort_evd_b_opd( int       m_A,
                              double*   l, int inc_l,
                              double*   V, int rs_V, int cs_V );
FLA_Error FLA_Sort_evd_f_opc( int       m_A,
                              float*    l, int inc_l,
                              scomplex* V, int rs_V, int cs_V );
FLA_Error FLA_Sort_evd_b_opc( int       m_A,
                              float*    l, int inc_l,
                              scomplex* V, int rs_V, int cs_V );
FLA_Error FLA_Sort_evd_f_opz( int       m_A,
                              double*   l, int inc_l,
                              dcomplex* V, int rs_V, int cs_V );
FLA_Error FLA_Sort_evd_b_opz( int       m_A,
                              double*   l, int inc_l,
                              dcomplex* V, int rs_V, int cs_V );

FLA_Error FLA_Sort_bsvd_ext( FLA_Direct direct, FLA_Obj s,
                             FLA_Bool apply_U, FLA_Obj U,
                             FLA_Bool apply_V, FLA_Obj V,
                             FLA_Bool apply_C, FLA_Obj C );
FLA_Error FLA_Sort_bsvd_ext_f_ops( int m_s, float* s, int inc_s,
                                   int m_U, float* U, int rs_U, int cs_U,
                                   int m_V, float* V, int rs_V, int cs_V,
                                   int n_C, float* C, int rs_C, int cs_C );
FLA_Error FLA_Sort_bsvd_ext_b_ops( int m_s, float* s, int inc_s,
                                   int m_U, float* U, int rs_U, int cs_U,
                                   int m_V, float* V, int rs_V, int cs_V,
                                   int n_C, float* C, int rs_C, int cs_C );
FLA_Error FLA_Sort_bsvd_ext_f_opd( int m_s, double* s, int inc_s,
                                   int m_U, double* U, int rs_U, int cs_U,
                                   int m_V, double* V, int rs_V, int cs_V,
                                   int n_C, double* C, int rs_C, int cs_C );
FLA_Error FLA_Sort_bsvd_ext_b_opd( int m_s, double* s, int inc_s,
                                   int m_U, double* U, int rs_U, int cs_U,
                                   int m_V, double* V, int rs_V, int cs_V,
                                   int n_C, double* C, int rs_C, int cs_C );
FLA_Error FLA_Sort_bsvd_ext_f_opc( int m_s, float*    s, int inc_s,
                                   int m_U, scomplex* U, int rs_U, int cs_U,
                                   int m_V, scomplex* V, int rs_V, int cs_V,
                                   int n_C, scomplex* C, int rs_C, int cs_C );
FLA_Error FLA_Sort_bsvd_ext_b_opc( int m_s, float*    s, int inc_s,
                                   int m_U, scomplex* U, int rs_U, int cs_U,
                                   int m_V, scomplex* V, int rs_V, int cs_V,
                                   int n_C, scomplex* C, int rs_C, int cs_C );
FLA_Error FLA_Sort_bsvd_ext_f_opz( int m_s, double*   s, int inc_s,
                                   int m_U, dcomplex* U, int rs_U, int cs_U,
                                   int m_V, dcomplex* V, int rs_V, int cs_V,
                                   int n_C, dcomplex* C, int rs_C, int cs_C );
FLA_Error FLA_Sort_bsvd_ext_b_opz( int m_s, double*   s, int inc_s,
                                   int m_U, dcomplex* U, int rs_U, int cs_U,
                                   int m_V, dcomplex* V, int rs_V, int cs_V,
                                   int n_C, dcomplex* C, int rs_C, int cs_C );

FLA_Error FLA_Sort_svd( FLA_Direct direct, FLA_Obj s, FLA_Obj U, FLA_Obj V );
FLA_Error FLA_Sort_svd_f_ops( int       m_U,
                              int       n_V,
                              float*    s, int inc_s,
                              float*    U, int rs_U, int cs_U,
                              float*    V, int rs_V, int cs_V );
FLA_Error FLA_Sort_svd_b_ops( int       m_U,
                              int       n_V,
                              float*    s, int inc_s,
                              float*    U, int rs_U, int cs_U,
                              float*    V, int rs_V, int cs_V );
FLA_Error FLA_Sort_svd_f_opd( int       m_U,
                              int       n_V,
                              double*   s, int inc_s,
                              double*   U, int rs_U, int cs_U,
                              double*   V, int rs_V, int cs_V );
FLA_Error FLA_Sort_svd_b_opd( int       m_U,
                              int       n_V,
                              double*   s, int inc_s,
                              double*   U, int rs_U, int cs_U,
                              double*   V, int rs_V, int cs_V );
FLA_Error FLA_Sort_svd_f_opc( int       m_U,
                              int       n_V,
                              float*    s, int inc_s,
                              scomplex* U, int rs_U, int cs_U,
                              scomplex* V, int rs_V, int cs_V );
FLA_Error FLA_Sort_svd_b_opc( int       m_U,
                              int       n_V,
                              float*    s, int inc_s,
                              scomplex* U, int rs_U, int cs_U,
                              scomplex* V, int rs_V, int cs_V );
FLA_Error FLA_Sort_svd_f_opz( int       m_U,
                              int       n_V,
                              double*   s, int inc_s,
                              dcomplex* U, int rs_U, int cs_U,
                              dcomplex* V, int rs_V, int cs_V );
FLA_Error FLA_Sort_svd_b_opz( int       m_U,
                              int       n_V,
                              double*   s, int inc_s,
                              dcomplex* U, int rs_U, int cs_U,
                              dcomplex* V, int rs_V, int cs_V );

FLA_Error FLA_Sv_2x2( FLA_Obj alpha11, FLA_Obj alpha12, FLA_Obj alpha22,
                      FLA_Obj sigma1, FLA_Obj sigma2 );
FLA_Error FLA_Sv_2x2_ops( float*    alpha11,
                          float*    alpha12,
                          float*    alpha22,
                          float*    sigma1,
                          float*    sigma2 );
FLA_Error FLA_Sv_2x2_opd( double*   alpha11,
                          double*   alpha12,
                          double*   alpha22,
                          double*   sigma1,
                          double*   sigma2 );

FLA_Error FLA_Svv_2x2( FLA_Obj alpha11, FLA_Obj alpha12, FLA_Obj alpha22,
                       FLA_Obj sigma1, FLA_Obj sigma2,
                       FLA_Obj gammaL, FLA_Obj sigmaL,
                       FLA_Obj gammaR, FLA_Obj sigmaR );
FLA_Error FLA_Svv_2x2_ops( float*    alpha11,
                           float*    alpha12,
                           float*    alpha22,
                           float*    sigma1,
                           float*    sigma2,
                           float*    gammaL,
                           float*    sigmaL,
                           float*    gammaR,
                           float*    sigmaR );
FLA_Error FLA_Svv_2x2_opd( double*   alpha11,
                           double*   alpha12,
                           double*   alpha22,
                           double*   sigma1,
                           double*   sigma2,
                           double*   gammaL,
                           double*   sigmaL,
                           double*   gammaR,
                           double*   sigmaR );

FLA_Error FLA_Mach_params( FLA_Machval machval, FLA_Obj val );
float     FLA_Mach_params_ops( FLA_Machval machval );
double    FLA_Mach_params_opd( FLA_Machval machval );

FLA_Error FLA_Apply_diag_matrix( FLA_Side side, FLA_Conj conj, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Shift_pivots_to( FLA_Pivot_type ptype, FLA_Obj p );
FLA_Error FLA_Form_perm_matrix( FLA_Obj p, FLA_Obj A );
FLA_Error FLA_LU_find_zero_on_diagonal( FLA_Obj A );

// --- f2c-converted routine prototypes ----------------------------------------

doublereal fla_dlamch( char* cmach, ftnlen cmach_len );
real       fla_slamch( char* cmach, ftnlen cmach_len );
logical    fla_lsame( char* ca, char* cb, ftnlen ca_len, ftnlen cb_len );
double     fla_pow_di( doublereal* a, integer* n );
real       fla_pow_ri( real* a, integer* n );

// --- LAPACK-related utility check routine prototypes -------------------------

FLA_Error FLA_Househ2_UT_check( FLA_Side side, FLA_Obj chi_1, FLA_Obj x2, FLA_Obj tau );
FLA_Error FLA_Househ3UD_UT_check( FLA_Obj chi_1, FLA_Obj x2, FLA_Obj y2, FLA_Obj tau );
FLA_Error FLA_Househ2s_UT_check( FLA_Side side, FLA_Obj chi_1, FLA_Obj x2, FLA_Obj alpha, FLA_Obj chi_1_minus_alpha, FLA_Obj tau );

FLA_Error FLA_Givens2_check( FLA_Obj chi_1, FLA_Obj chi_2, FLA_Obj gamma, FLA_Obj sigma, FLA_Obj chi_1_new );
FLA_Error FLA_Apply_GTG_check( FLA_Obj gamma, FLA_Obj sigma, FLA_Obj delta1, FLA_Obj epsilon1, FLA_Obj delta2 );
FLA_Error FLA_Apply_G_1x2_check( FLA_Obj gamma, FLA_Obj sigma, FLA_Obj beta, FLA_Obj epsilon );
FLA_Error FLA_Apply_G_mx2_check( FLA_Obj gamma, FLA_Obj sigma, FLA_Obj a1, FLA_Obj a2 );
FLA_Error FLA_Apply_G_check( FLA_Side side, FLA_Direct direct, FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Wilkshift_tridiag_check( FLA_Obj delta1, FLA_Obj epsilon, FLA_Obj delta2, FLA_Obj kappa );
FLA_Error FLA_Wilkshift_bidiag_check( FLA_Obj epsilon1, FLA_Obj delta1, FLA_Obj epsilon2, FLA_Obj delta2, FLA_Obj kappa );
FLA_Error FLA_Introduce_bulge_check( FLA_Obj shift, FLA_Obj gamma, FLA_Obj sigma, FLA_Obj delta1, FLA_Obj epsilon1, FLA_Obj delta2, FLA_Obj beta, FLA_Obj epsilon2 );
FLA_Error FLA_Mach_params_check( FLA_Machval machval, FLA_Obj val );

FLA_Error FLA_Sort_evd_check( FLA_Direct direct, FLA_Obj l, FLA_Obj V );
FLA_Error FLA_Sort_svd_check( FLA_Direct direct, FLA_Obj s, FLA_Obj U, FLA_Obj V );

FLA_Error FLA_Apply_diag_matrix_check( FLA_Side side, FLA_Conj conj, FLA_Obj x, FLA_Obj A );
FLA_Error FLA_Shift_pivots_to_check( FLA_Pivot_type ptype, FLA_Obj p );
FLA_Error FLA_Form_perm_matrix_check( FLA_Obj p, FLA_Obj A );
FLA_Error FLA_LU_find_zero_on_diagonal_check( FLA_Obj A );

