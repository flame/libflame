/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Tevd_iteracc_v.h"
#include "FLA_Tevd_eigval_v.h"
#include "FLA_Tevd_francis_v.h"

// --- FLA_Tevd_compute_scaling() ----------------------------------------------

FLA_Error FLA_Tevd_compute_scaling_ops( int       m_A,
                                        float*    buff_d, int inc_d, 
                                        float*    buff_e, int inc_e,
                                        float*    sigma );
FLA_Error FLA_Tevd_compute_scaling_opd( int       m_A,
                                        double*   buff_d, int inc_d, 
                                        double*   buff_e, int inc_e,
                                        double*   sigma );

// --- FLA_Tevd_find_submatrix() -----------------------------------------------

FLA_Error FLA_Tevd_find_submatrix_ops( int       m_A,
                                       int       ij_begin,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR );
FLA_Error FLA_Tevd_find_submatrix_opd( int       m_A,
                                       int       ij_begin,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR );

// --- FLA_Tevd_find_perfshift() -----------------------------------------------

FLA_Error FLA_Tevd_find_perfshift_ops( int       m_d,
                                       int       m_l,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e, 
                                       float*    buff_l, int inc_l, 
                                       int*      buff_lstat, int inc_lstat, 
                                       float*    buff_pu, int inc_pu, 
                                       int*      ij_shift );
FLA_Error FLA_Tevd_find_perfshift_opd( int       m_d,
                                       int       m_l,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e, 
                                       double*   buff_l, int inc_l, 
                                       int*      buff_lstat, int inc_lstat, 
                                       double*   buff_pu, int inc_pu, 
                                       int*      ij_shift );

// --- FLA_Norm1_tridiag() -----------------------------------------------------

FLA_Error FLA_Norm1_tridiag( FLA_Obj d, FLA_Obj e, FLA_Obj norm );
FLA_Error FLA_Norm1_tridiag_ops( int       m_A,
                                 float*    buff_d, int inc_d, 
                                 float*    buff_e, int inc_e,
                                 float*    norm );
FLA_Error FLA_Norm1_tridiag_opd( int       m_A,
                                 double*   buff_d, int inc_d, 
                                 double*   buff_e, int inc_e,
                                 double*   norm );

// --- FLA_Tevd_v_opt_var1() ---------------------------------------------------

FLA_Error FLA_Tevd_v_opt_var1( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U, dim_t b_alg );
FLA_Error FLA_Tevd_v_ops_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opd_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opc_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               scomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opz_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               dcomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );

// --- FLA_Tevd_v_opt_var2() ---------------------------------------------------

FLA_Error FLA_Tevd_v_opt_var2( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj R, FLA_Obj W, FLA_Obj U, dim_t b_alg );
FLA_Error FLA_Tevd_v_ops_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_G_extra,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_R, int rs_R, int cs_R,
                               float*    buff_W, int rs_W, int cs_W,
                               float*    buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opd_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_G_extra,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_R, int rs_R, int cs_R,
                               double*   buff_W, int rs_W, int cs_W,
                               double*   buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opc_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_G_extra,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G,
                               float*    buff_R, int rs_R, int cs_R,
                               scomplex* buff_W, int rs_W, int cs_W,
                               scomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );
FLA_Error FLA_Tevd_v_opz_var2( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_G_extra,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G,
                               double*   buff_R, int rs_R, int cs_R,
                               dcomplex* buff_W, int rs_W, int cs_W,
                               dcomplex* buff_U, int rs_U, int cs_U,
                               int       b_alg );

