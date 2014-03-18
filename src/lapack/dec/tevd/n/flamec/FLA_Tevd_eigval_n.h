
// --- FLA_Tevd_eigval_n_opt_var1() --------------------------------------------

FLA_Error FLA_Tevd_eigval_n_opt_var1( FLA_Obj G, FLA_Obj d, FLA_Obj e, FLA_Obj n_iter );
FLA_Error FLA_Tevd_eigval_n_ops_var1( int       m_A,
                                      int       n_G,
                                      float*    buff_d, int inc_d, 
                                      float*    buff_e, int inc_e,
                                      int*      n_iter );
FLA_Error FLA_Tevd_eigval_n_opd_var1( int       m_A,
                                      int       n_G,
                                      double*   buff_d, int inc_d, 
                                      double*   buff_e, int inc_e,
                                      int*      n_iter );

