
// --- FLA_Bsvd_iteracc_v_opt_var1() -------------------------------------------

FLA_Error FLA_Bsvd_iteracc_v_ops_var1( int       m_A,
                                       int       n_GH,
                                       int       ijTL,
                                       float     tol,
                                       float     thresh,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       scomplex* buff_G, int rs_G, int cs_G,
                                       scomplex* buff_H, int rs_H, int cs_H,
                                       int*      n_iter_perf );
FLA_Error FLA_Bsvd_iteracc_v_opd_var1( int       m_A,
                                       int       n_GH,
                                       int       ijTL,
                                       double    tol,
                                       double    thresh,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       dcomplex* buff_G, int rs_G, int cs_G,
                                       dcomplex* buff_H, int rs_H, int cs_H,
                                       int*      n_iter_perf );
