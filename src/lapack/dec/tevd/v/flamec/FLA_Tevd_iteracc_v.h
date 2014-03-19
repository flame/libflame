/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Tevd_iteracc_v_opt_var1() -------------------------------------------

FLA_Error FLA_Tevd_iteracc_v_ops_var1( int       m_A,
                                       int       n_G,
                                       int       ijTL,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       scomplex* buff_G, int rs_G, int cs_G,
                                       int*      n_iter_perf );
FLA_Error FLA_Tevd_iteracc_v_opd_var1( int       m_A,
                                       int       n_G,
                                       int       ijTL,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       dcomplex* buff_G, int rs_G, int cs_G,
                                       int*      n_iter_perf );

FLA_Error FLA_Tevd_iteracc_v_ops_var3( int       m_A,
                                       int       m_U,
                                       int       n_G,
                                       int       ijTL,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       float*    buff_l, int inc_l,
                                       int*      buff_ls, int inc_ls,
                                       float*    buff_pu, int inc_pu,
                                       scomplex* buff_G, int rs_G, int cs_G,
                                       int*      n_iter_perf );
FLA_Error FLA_Tevd_iteracc_v_opd_var3( int       m_A,
                                       int       m_U,
                                       int       n_G,
                                       int       ijTL,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       double*   buff_l, int inc_l,
                                       int*      buff_ls, int inc_ls,
                                       double*   buff_pu, int inc_pu,
                                       dcomplex* buff_G, int rs_G, int cs_G,
                                       int*      n_iter_perf );

