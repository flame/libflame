/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Bsvd_iteracc_v.h"
#include "FLA_Bsvd_sinval_v.h"
#include "FLA_Bsvd_francis_v.h"

// --- FLA_Bsvd_compute_shift() ------------------------------------------------

FLA_Error FLA_Bsvd_compute_shift( FLA_Obj tol, FLA_Obj sminl, FLA_Obj smax, FLA_Obj d, FLA_Obj e, FLA_Obj shift );
FLA_Error FLA_Bsvd_compute_shift_ops( integer       m_A,
                                      float     tol,
                                      float     sminl,
                                      float     smax,
                                      float*    buff_d, integer inc_d,
                                      float*    buff_e, integer inc_e,
                                      float*    shift );
FLA_Error FLA_Bsvd_compute_shift_opd( integer       m_A,
                                      double    tol,
                                      double    sminl,
                                      double    smax,
                                      double*   buff_d, integer inc_d,
                                      double*   buff_e, integer inc_e,
                                      double*   shift );

// --- FLA_Bsvd_compute_tol_thresh() -------------------------------------------

FLA_Error FLA_Bsvd_compute_tol_thresh( FLA_Obj tolmul, FLA_Obj maxit, FLA_Obj d, FLA_Obj e, FLA_Obj tol, FLA_Obj thresh );
FLA_Error FLA_Bsvd_compute_tol_thresh_ops( integer       m_A,
                                           float     tolmul,
                                           float     maxit,
                                           float*    buff_d, integer inc_d, 
                                           float*    buff_e, integer inc_e, 
                                           float*    tol,
                                           float*    thresh );
FLA_Error FLA_Bsvd_compute_tol_thresh_opd( integer       m_A,
                                           double    tolmul,
                                           double    maxit,
                                           double*   buff_d, integer inc_d, 
                                           double*   buff_e, integer inc_e, 
                                           double*   tol,
                                           double*   thresh );

// --- FLA_Bsvd_find_converged() -----------------------------------------------

FLA_Error FLA_Bsvd_find_converged( FLA_Obj tol, FLA_Obj d, FLA_Obj e, FLA_Obj sminl );
FLA_Error FLA_Bsvd_find_converged_ops( integer       m_A,
                                       float     tol, 
                                       float*    buff_d, integer inc_d, 
                                       float*    buff_e, integer inc_e,
                                       float*    sminl );
FLA_Error FLA_Bsvd_find_converged_opd( integer       m_A,
                                       double    tol, 
                                       double*   buff_d, integer inc_d, 
                                       double*   buff_e, integer inc_e,
                                       double*   sminl );

// --- FLA_Bsvd_find_max_min() -------------------------------------------------

FLA_Error FLA_Bsvd_find_max_min( FLA_Obj d, FLA_Obj e, FLA_Obj smax, FLA_Obj smin );
FLA_Error FLA_Bsvd_find_max_min_ops( integer       m_A,
                                     float*    buff_d, integer inc_d, 
                                     float*    buff_e, integer inc_e, 
                                     float*    smax,
                                     float*    smin );
FLA_Error FLA_Bsvd_find_max_min_opd( integer       m_A,
                                     double*   buff_d, integer inc_d, 
                                     double*   buff_e, integer inc_e, 
                                     double*   smax,
                                     double*   smin );

// --- FLA_Bsvd_find_submatrix() -----------------------------------------------

FLA_Error FLA_Bsvd_find_submatrix_ops( integer       mn_A,
                                       integer       ij_begin,
                                       float*    buff_d, integer inc_d,
                                       float*    buff_e, integer inc_e,
                                       integer*      ijTL,
                                       integer*      ijBR );
FLA_Error FLA_Bsvd_find_submatrix_opd( integer       mn_A,
                                       integer       ij_begin,
                                       double*   buff_d, integer inc_d,
                                       double*   buff_e, integer inc_e,
                                       integer*      ijTL,
                                       integer*      ijBR );

// --- FLA_Bsvd_v_opt_var1() ---------------------------------------------------

FLA_Error FLA_Bsvd_v_opt_var1( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H, FLA_Obj U, FLA_Obj V, dim_t b_alg );
FLA_Error FLA_Bsvd_v_ops_var1( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               scomplex* buff_H, integer rs_H, integer cs_H,
                               float*    buff_U, integer rs_U, integer cs_U,
                               float*    buff_V, integer rs_V, integer cs_V,
                               integer       b_alg );
FLA_Error FLA_Bsvd_v_opd_var1( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               dcomplex* buff_H, integer rs_H, integer cs_H,
                               double*   buff_U, integer rs_U, integer cs_U,
                               double*   buff_V, integer rs_V, integer cs_V,
                               integer       b_alg );
FLA_Error FLA_Bsvd_v_opc_var1( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               scomplex* buff_H, integer rs_H, integer cs_H,
                               scomplex* buff_U, integer rs_U, integer cs_U,
                               scomplex* buff_V, integer rs_V, integer cs_V,
                               integer       b_alg );
FLA_Error FLA_Bsvd_v_opz_var1( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               dcomplex* buff_H, integer rs_H, integer cs_H,
                               dcomplex* buff_U, integer rs_U, integer cs_U,
                               dcomplex* buff_V, integer rs_V, integer cs_V,
                               integer       b_alg );

// --- FLA_Bsvd_v_opt_var2() ---------------------------------------------------

FLA_Error FLA_Bsvd_v_opt_var2( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H, FLA_Obj RG, FLA_Obj RH, FLA_Obj W, FLA_Obj U, FLA_Obj V, dim_t b_alg );
FLA_Error FLA_Bsvd_v_ops_var2( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               scomplex* buff_H, integer rs_H, integer cs_H,
                               float*    buff_RG, integer rs_RG, integer cs_RG,
                               float*    buff_RH, integer rs_RH, integer cs_RH,
                               float*    buff_W, integer rs_W, integer cs_W,
                               float*    buff_U, integer rs_U, integer cs_U,
                               float*    buff_V, integer rs_V, integer cs_V,
                               integer       b_alg );
FLA_Error FLA_Bsvd_v_opd_var2( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               dcomplex* buff_H, integer rs_H, integer cs_H,
                               double*   buff_RG, integer rs_RG, integer cs_RG,
                               double*   buff_RH, integer rs_RH, integer cs_RH,
                               double*   buff_W, integer rs_W, integer cs_W,
                               double*   buff_U, integer rs_U, integer cs_U,
                               double*   buff_V, integer rs_V, integer cs_V,
                               integer       b_alg );
FLA_Error FLA_Bsvd_v_opc_var2( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               scomplex* buff_H, integer rs_H, integer cs_H,
                               float*    buff_RG, integer rs_RG, integer cs_RG,
                               float*    buff_RH, integer rs_RH, integer cs_RH,
                               scomplex* buff_W, integer rs_W, integer cs_W,
                               scomplex* buff_U, integer rs_U, integer cs_U,
                               scomplex* buff_V, integer rs_V, integer cs_V,
                               integer       b_alg );
FLA_Error FLA_Bsvd_v_opz_var2( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               dcomplex* buff_H, integer rs_H, integer cs_H,
                               double*   buff_RG, integer rs_RG, integer cs_RG,
                               double*   buff_RH, integer rs_RH, integer cs_RH,
                               dcomplex* buff_W, integer rs_W, integer cs_W,
                               dcomplex* buff_U, integer rs_U, integer cs_U,
                               dcomplex* buff_V, integer rs_V, integer cs_V,
                               integer       b_alg );

