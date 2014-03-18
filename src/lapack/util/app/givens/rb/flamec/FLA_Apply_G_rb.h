
FLA_Error FLA_Apply_G_rb_opt_var1( FLA_Obj c, FLA_Obj s, FLA_Obj A );
FLA_Error FLA_Apply_G_rb_ops_var1( int       m_A,
                                   int       n_A,
                                   float*    buff_c, int inc_c,
                                   float*    buff_s, int inc_s,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rb_opd_var1( int       m_A,
                                   int       n_A,
                                   double*   buff_c, int inc_c,
                                   double*   buff_s, int inc_s,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rb_opc_var1( int       m_A,
                                   int       n_A,
                                   float*    buff_c, int inc_c,
                                   float*    buff_s, int inc_s,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rb_opz_var1( int       m_A,
                                   int       n_A,
                                   double*   buff_c, int inc_c,
                                   double*   buff_s, int inc_s,
                                   dcomplex* buff_A, int rs_A, int cs_A );

