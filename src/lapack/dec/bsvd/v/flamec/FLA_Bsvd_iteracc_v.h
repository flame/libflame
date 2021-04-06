/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Bsvd_iteracc_v_opt_var1() -------------------------------------------

FLA_Error FLA_Bsvd_iteracc_v_ops_var1( integer       m_A,
                                       integer       n_GH,
                                       integer       ijTL,
                                       float     tol,
                                       float     thresh,
                                       float*    buff_d, integer inc_d, 
                                       float*    buff_e, integer inc_e,
                                       scomplex* buff_G, integer rs_G, integer cs_G,
                                       scomplex* buff_H, integer rs_H, integer cs_H,
                                       integer*      n_iter_perf );
FLA_Error FLA_Bsvd_iteracc_v_opd_var1( integer       m_A,
                                       integer       n_GH,
                                       integer       ijTL,
                                       double    tol,
                                       double    thresh,
                                       double*   buff_d, integer inc_d, 
                                       double*   buff_e, integer inc_e,
                                       dcomplex* buff_G, integer rs_G, integer cs_G,
                                       dcomplex* buff_H, integer rs_H, integer cs_H,
                                       integer*      n_iter_perf );
