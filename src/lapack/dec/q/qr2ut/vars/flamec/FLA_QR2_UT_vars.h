/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_QR2_UT_blk_var1( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl );
FLA_Error FLA_QR2_UT_blk_var2( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T, fla_qr2ut_t* cntl );

FLA_Error FLA_QR2_UT_unb_var1( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T );

FLA_Error FLA_QR2_UT_opt_var1( FLA_Obj U,
                               FLA_Obj D, FLA_Obj T );

FLA_Error FLA_QR2_UT_ops_var1( int m_UT,
                               int m_D,
                               float* U, int rs_U, int cs_U,
                               float* D, int rs_D, int cs_D,
                               float* T, int rs_T, int cs_T );
FLA_Error FLA_QR2_UT_opd_var1( int m_UT,
                               int m_D,
                               double* U, int rs_U, int cs_U,
                               double* D, int rs_D, int cs_D,
                               double* T, int rs_T, int cs_T );
FLA_Error FLA_QR2_UT_opc_var1( int m_UT,
                               int m_D,
                               scomplex* U, int rs_U, int cs_U,
                               scomplex* D, int rs_D, int cs_D,
                               scomplex* T, int rs_T, int cs_T );
FLA_Error FLA_QR2_UT_opz_var1( int m_UT,
                               int m_D,
                               dcomplex* U, int rs_U, int cs_U,
                               dcomplex* D, int rs_D, int cs_D,
                               dcomplex* T, int rs_T, int cs_T );
