/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_H2_UT_r_unb_var1( FLA_Obj tau, FLA_Obj u2h,
                                      FLA_Obj a1, FLA_Obj A2 );

FLA_Error FLA_Apply_H2_UT_r_opt_var1( FLA_Obj tau, FLA_Obj u2h,
                                      FLA_Obj a1, FLA_Obj A2 );

FLA_Error FLA_Apply_H2_UT_r_ops_var1( int n_u2h_A2,
                                      int m_a1,
                                      float* tau,
                                      float* u2h, int inc_u2h,
                                      float* a1, int inc_a1,
                                      float* A2, int rs_A2, int cs_A2 );

FLA_Error FLA_Apply_H2_UT_r_opd_var1( int n_u2h_A2,
                                      int m_a1,
                                      double* tau,
                                      double* u2h, int inc_u2h,
                                      double* a1, int inc_a1,
                                      double* A2, int rs_A2, int cs_A2 );

FLA_Error FLA_Apply_H2_UT_r_opc_var1( int n_u2h_A2,
                                      int m_a1,
                                      scomplex* tau,
                                      scomplex* u2h, int inc_u2h,
                                      scomplex* a1, int inc_a1,
                                      scomplex* A2, int rs_A2, int cs_A2 );

FLA_Error FLA_Apply_H2_UT_r_opz_var1( int n_u2h_A2,
                                      int m_a1,
                                      dcomplex* tau,
                                      dcomplex* u2h, int inc_u2h,
                                      dcomplex* a1, int inc_a1,
                                      dcomplex* A2, int rs_A2, int cs_A2 );

