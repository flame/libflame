/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_UDdate_UT_blk_var1( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl );
FLA_Error FLA_UDdate_UT_blk_var2( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T, fla_uddateut_t* cntl );

FLA_Error FLA_UDdate_UT_unb_var1( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T );

FLA_Error FLA_UDdate_UT_opt_var1( FLA_Obj R,
                                  FLA_Obj C,
                                  FLA_Obj D, FLA_Obj T );
FLA_Error FLA_UDdate_UT_ops_var1( integer mn_RT,
                                  integer m_C,
                                  integer m_D,
                                  float* R, integer rs_R, integer cs_R,
                                  float* C, integer rs_C, integer cs_C,
                                  float* D, integer rs_D, integer cs_D,
                                  float* T, integer rs_T, integer cs_T );
FLA_Error FLA_UDdate_UT_opd_var1( integer mn_RT,
                                  integer m_C,
                                  integer m_D,
                                  double* R, integer rs_R, integer cs_R,
                                  double* C, integer rs_C, integer cs_C,
                                  double* D, integer rs_D, integer cs_D,
                                  double* T, integer rs_T, integer cs_T );
FLA_Error FLA_UDdate_UT_opc_var1( integer mn_RT,
                                  integer m_C,
                                  integer m_D,
                                  scomplex* R, integer rs_R, integer cs_R,
                                  scomplex* C, integer rs_C, integer cs_C,
                                  scomplex* D, integer rs_D, integer cs_D,
                                  scomplex* T, integer rs_T, integer cs_T );
FLA_Error FLA_UDdate_UT_opz_var1( integer mn_RT,
                                  integer m_C,
                                  integer m_D,
                                  dcomplex* R, integer rs_R, integer cs_R,
                                  dcomplex* C, integer rs_C, integer cs_C,
                                  dcomplex* D, integer rs_D, integer cs_D,
                                  dcomplex* T, integer rs_T, integer cs_T );

