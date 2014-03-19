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
FLA_Error FLA_UDdate_UT_ops_var1( int mn_RT,
                                  int m_C,
                                  int m_D,
                                  float* R, int rs_R, int cs_R,
                                  float* C, int rs_C, int cs_C,
                                  float* D, int rs_D, int cs_D,
                                  float* T, int rs_T, int cs_T );
FLA_Error FLA_UDdate_UT_opd_var1( int mn_RT,
                                  int m_C,
                                  int m_D,
                                  double* R, int rs_R, int cs_R,
                                  double* C, int rs_C, int cs_C,
                                  double* D, int rs_D, int cs_D,
                                  double* T, int rs_T, int cs_T );
FLA_Error FLA_UDdate_UT_opc_var1( int mn_RT,
                                  int m_C,
                                  int m_D,
                                  scomplex* R, int rs_R, int cs_R,
                                  scomplex* C, int rs_C, int cs_C,
                                  scomplex* D, int rs_D, int cs_D,
                                  scomplex* T, int rs_T, int cs_T );
FLA_Error FLA_UDdate_UT_opz_var1( int mn_RT,
                                  int m_C,
                                  int m_D,
                                  dcomplex* R, int rs_R, int cs_R,
                                  dcomplex* C, int rs_C, int cs_C,
                                  dcomplex* D, int rs_D, int cs_D,
                                  dcomplex* T, int rs_T, int cs_T );

