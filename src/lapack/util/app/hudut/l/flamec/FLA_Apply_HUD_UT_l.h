/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_HUD_UT_l_unb_var1( FLA_Obj tau, FLA_Obj w12t,
                                                    FLA_Obj r12t,
                                       FLA_Obj u1,  FLA_Obj C2,
                                       FLA_Obj v1,  FLA_Obj D2 );

FLA_Error FLA_Apply_HUD_UT_l_opt_var1( FLA_Obj tau, FLA_Obj w12t,
                                                    FLA_Obj r12t,
                                       FLA_Obj u1,  FLA_Obj C2,
                                       FLA_Obj v1,  FLA_Obj D2 );

FLA_Error FLA_Apply_HUD_UT_l_ops_var1( int m_u1_C2,
                                       int m_v1_D2,
                                       int n_r12t,
                                       float* tau,
                                       float* w12t, int inc_w12t,
                                       float* r12t, int inc_r12t,
                                       float* u1, int inc_u1,
                                       float* C2, int rs_C2, int cs_C2,
                                       float* v1, int inc_v1,
                                       float* D2, int rs_D2, int cs_D2 );

FLA_Error FLA_Apply_HUD_UT_l_opd_var1( int m_u1_C2,
                                       int m_v1_D2,
                                       int n_r12t,
                                       double* tau,
                                       double* w12t, int inc_w12t,
                                       double* r12t, int inc_r12t,
                                       double* u1, int inc_u1,
                                       double* C2, int rs_C2, int cs_C2,
                                       double* v1, int inc_v1,
                                       double* D2, int rs_D2, int cs_D2 );

FLA_Error FLA_Apply_HUD_UT_l_opc_var1( int m_u1_C2,
                                       int m_v1_D2,
                                       int n_r12t,
                                       scomplex* tau,
                                       scomplex* w12t, int inc_w12t,
                                       scomplex* r12t, int inc_r12t,
                                       scomplex* u1, int inc_u1,
                                       scomplex* C2, int rs_C2, int cs_C2,
                                       scomplex* v1, int inc_v1,
                                       scomplex* D2, int rs_D2, int cs_D2 );

FLA_Error FLA_Apply_HUD_UT_l_opz_var1( int m_u1_C2,
                                       int m_v1_D2,
                                       int n_r12t,
                                       dcomplex* tau,
                                       dcomplex* w12t, int inc_w12t,
                                       dcomplex* r12t, int inc_r12t,
                                       dcomplex* u1, int inc_u1,
                                       dcomplex* C2, int rs_C2, int cs_C2,
                                       dcomplex* v1, int inc_v1,
                                       dcomplex* D2, int rs_D2, int cs_D2 );

