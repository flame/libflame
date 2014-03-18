

// BLAS 2 version
FLA_Error FLA_QR_UT_piv_unb_var1( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p );
FLA_Error FLA_QR_UT_piv_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p, fla_qrut_t* cntl );

// BLAS 3 version
FLA_Error FLA_QR_UT_piv_unb_var2( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p );
FLA_Error FLA_QR_UT_piv_blk_var2( FLA_Obj A, FLA_Obj T, FLA_Obj w, FLA_Obj p, fla_qrut_t* cntl );
FLA_Error FLA_Apply_H2_UT_piv_row( FLA_Obj tau, FLA_Obj a1t, FLA_Obj u1t, FLA_Obj W,
                                   FLA_Obj u2,  FLA_Obj A2,  FLA_Obj U2,  FLA_Obj w1t,
                                   FLA_Obj vt );

