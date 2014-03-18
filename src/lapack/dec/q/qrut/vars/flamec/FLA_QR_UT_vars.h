
FLA_Error FLA_QR_UT_unb_var1( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_QR_UT_blk_var1( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_opt_var1( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_QR_UT_ops_var1( int m_A,
                              int n_A,
                              float* A, int rs_A, int cs_A,
                              float* t, int inc_t );
FLA_Error FLA_QR_UT_opd_var1( int m_A,
                              int n_A,
                              double* A, int rs_A, int cs_A,
                              double* t, int inc_t );
FLA_Error FLA_QR_UT_opc_var1( int m_A,
                              int n_A,
                              scomplex* A, int rs_A, int cs_A,
                              scomplex* t, int inc_t );
FLA_Error FLA_QR_UT_opz_var1( int m_A,
                              int n_A,
                              dcomplex* A, int rs_A, int cs_A,
                              dcomplex* t, int inc_t );

FLA_Error FLA_QR_UT_unb_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_QR_UT_blk_var2( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_opt_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_QR_UT_ops_var2( int m_A,
                              int n_A,
                              float* A, int rs_A, int cs_A,
                              float* T, int rs_T, int cs_T );
FLA_Error FLA_QR_UT_opd_var2( int m_A,
                              int n_A,
                              double* A, int rs_A, int cs_A,
                              double* T, int rs_T, int cs_T );
FLA_Error FLA_QR_UT_opc_var2( int m_A,
                              int n_A,
                              scomplex* A, int rs_A, int cs_A,
                              scomplex* T, int rs_T, int cs_T );
FLA_Error FLA_QR_UT_opz_var2( int m_A,
                              int n_A,
                              dcomplex* A, int rs_A, int cs_A,
                              dcomplex* T, int rs_T, int cs_T );

FLA_Error FLA_QR_UT_blk_var3( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );

