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

FLA_Error FLA_QR2_UT_ops_var1( integer m_UT,
                               integer m_D,
                               float* U, integer rs_U, integer cs_U,
                               float* D, integer rs_D, integer cs_D,
                               float* T, integer rs_T, integer cs_T );
FLA_Error FLA_QR2_UT_opd_var1( integer m_UT,
                               integer m_D,
                               double* U, integer rs_U, integer cs_U,
                               double* D, integer rs_D, integer cs_D,
                               double* T, integer rs_T, integer cs_T );
FLA_Error FLA_QR2_UT_opc_var1( integer m_UT,
                               integer m_D,
                               scomplex* U, integer rs_U, integer cs_U,
                               scomplex* D, integer rs_D, integer cs_D,
                               scomplex* T, integer rs_T, integer cs_T );
FLA_Error FLA_QR2_UT_opz_var1( integer m_UT,
                               integer m_D,
                               dcomplex* U, integer rs_U, integer cs_U,
                               dcomplex* D, integer rs_D, integer cs_D,
                               dcomplex* T, integer rs_T, integer cs_T );
