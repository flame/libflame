/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_QR_UT_unb_var1( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_QR_UT_blk_var1( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_opt_var1( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_QR_UT_ops_var1( integer m_A,
                              integer n_A,
                              float* A, integer rs_A, integer cs_A,
                              float* t, integer inc_t );
FLA_Error FLA_QR_UT_opd_var1( integer m_A,
                              integer n_A,
                              double* A, integer rs_A, integer cs_A,
                              double* t, integer inc_t );
FLA_Error FLA_QR_UT_opc_var1( integer m_A,
                              integer n_A,
                              scomplex* A, integer rs_A, integer cs_A,
                              scomplex* t, integer inc_t );
FLA_Error FLA_QR_UT_opz_var1( integer m_A,
                              integer n_A,
                              dcomplex* A, integer rs_A, integer cs_A,
                              dcomplex* t, integer inc_t );

FLA_Error FLA_QR_UT_unb_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_QR_UT_blk_var2( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_opt_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_QR_UT_ops_var2( integer m_A,
                              integer n_A,
                              float* A, integer rs_A, integer cs_A,
                              float* T, integer rs_T, integer cs_T );
FLA_Error FLA_QR_UT_opd_var2( integer m_A,
                              integer n_A,
                              double* A, integer rs_A, integer cs_A,
                              double* T, integer rs_T, integer cs_T );
FLA_Error FLA_QR_UT_opc_var2( integer m_A,
                              integer n_A,
                              scomplex* A, integer rs_A, integer cs_A,
                              scomplex* T, integer rs_T, integer cs_T );
FLA_Error FLA_QR_UT_opz_var2( integer m_A,
                              integer n_A,
                              dcomplex* A, integer rs_A, integer cs_A,
                              dcomplex* T, integer rs_T, integer cs_T );

FLA_Error FLA_QR_UT_blk_var3( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );

