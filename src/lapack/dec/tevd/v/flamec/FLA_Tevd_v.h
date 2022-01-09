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

FLA_Error FLA_Tevd_compute_scaling_ops( integer       m_A,
                                        float*    buff_d, integer inc_d, 
                                        float*    buff_e, integer inc_e,
                                        float*    sigma );
FLA_Error FLA_Tevd_compute_scaling_opd( integer       m_A,
                                        double*   buff_d, integer inc_d, 
                                        double*   buff_e, integer inc_e,
                                        double*   sigma );

// --- FLA_Tevd_find_submatrix() -----------------------------------------------

FLA_Error FLA_Tevd_find_submatrix_ops( integer       m_A,
                                       integer       ij_begin,
                                       float*    buff_d, integer inc_d, 
                                       float*    buff_e, integer inc_e,
                                       integer*      ijTL,
                                       integer*      ijBR );
FLA_Error FLA_Tevd_find_submatrix_opd( integer       m_A,
                                       integer       ij_begin,
                                       double*   buff_d, integer inc_d, 
                                       double*   buff_e, integer inc_e,
                                       integer*      ijTL,
                                       integer*      ijBR );

// --- FLA_Tevd_find_perfshift() -----------------------------------------------

FLA_Error FLA_Tevd_find_perfshift_ops( integer       m_d,
                                       integer       m_l,
                                       float*    buff_d, integer inc_d, 
                                       float*    buff_e, integer inc_e, 
                                       float*    buff_l, integer inc_l, 
                                       integer*      buff_lstat, integer inc_lstat, 
                                       float*    buff_pu, integer inc_pu, 
                                       integer*      ij_shift );
FLA_Error FLA_Tevd_find_perfshift_opd( integer       m_d,
                                       integer       m_l,
                                       double*   buff_d, integer inc_d, 
                                       double*   buff_e, integer inc_e, 
                                       double*   buff_l, integer inc_l, 
                                       integer*      buff_lstat, integer inc_lstat, 
                                       double*   buff_pu, integer inc_pu, 
                                       integer*      ij_shift );

// --- FLA_Norm1_tridiag() -----------------------------------------------------

FLA_Error FLA_Norm1_tridiag( FLA_Obj d, FLA_Obj e, FLA_Obj norm );
FLA_Error FLA_Norm1_tridiag_ops( integer       m_A,
                                 float*    buff_d, integer inc_d, 
                                 float*    buff_e, integer inc_e,
                                 float*    norm );
FLA_Error FLA_Norm1_tridiag_opd( integer       m_A,
                                 double*   buff_d, integer inc_d, 
                                 double*   buff_e, integer inc_e,
                                 double*   norm );

// --- FLA_Tevd_v_opt_var1() ---------------------------------------------------

FLA_Error FLA_Tevd_v_opt_var1( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U, dim_t b_alg );
FLA_Error FLA_Tevd_v_ops_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               float*    buff_U, integer rs_U, integer cs_U,
                               integer       b_alg );
FLA_Error FLA_Tevd_v_opd_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               double*   buff_U, integer rs_U, integer cs_U,
                               integer       b_alg );
FLA_Error FLA_Tevd_v_opc_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               scomplex* buff_U, integer rs_U, integer cs_U,
                               integer       b_alg );
FLA_Error FLA_Tevd_v_opz_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               dcomplex* buff_U, integer rs_U, integer cs_U,
                               integer       b_alg );

// --- FLA_Tevd_v_opt_var2() ---------------------------------------------------

FLA_Error FLA_Tevd_v_opt_var2( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj R, FLA_Obj W, FLA_Obj U, dim_t b_alg );
FLA_Error FLA_Tevd_v_ops_var2( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_G_extra,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               float*    buff_R, integer rs_R, integer cs_R,
                               float*    buff_W, integer rs_W, integer cs_W,
                               float*    buff_U, integer rs_U, integer cs_U,
                               integer       b_alg );
FLA_Error FLA_Tevd_v_opd_var2( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_G_extra,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               double*   buff_R, integer rs_R, integer cs_R,
                               double*   buff_W, integer rs_W, integer cs_W,
                               double*   buff_U, integer rs_U, integer cs_U,
                               integer       b_alg );
FLA_Error FLA_Tevd_v_opc_var2( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_G_extra,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               float*    buff_R, integer rs_R, integer cs_R,
                               scomplex* buff_W, integer rs_W, integer cs_W,
                               scomplex* buff_U, integer rs_U, integer cs_U,
                               integer       b_alg );
FLA_Error FLA_Tevd_v_opz_var2( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_G_extra,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               double*   buff_R, integer rs_R, integer cs_R,
                               dcomplex* buff_W, integer rs_W, integer cs_W,
                               dcomplex* buff_U, integer rs_U, integer cs_U,
                               integer       b_alg );

