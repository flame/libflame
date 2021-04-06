/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_H2_UT_l_unb_var1( FLA_Obj tau, FLA_Obj u2, FLA_Obj a1t,
                                                               FLA_Obj A2 );

FLA_Error FLA_Apply_H2_UT_l_opt_var1( FLA_Obj tau, FLA_Obj u2, FLA_Obj a1t,
                                                               FLA_Obj A2 );

FLA_Error FLA_Apply_H2_UT_l_ops_var1( integer m_u2_A2,
                                      integer n_a1t,
                                      float* tau,
                                      float* u2, integer inc_u2,
                                      float* a1t, integer inc_a1t,
                                      float* A2, integer rs_A2, integer cs_A2 );

FLA_Error FLA_Apply_H2_UT_l_opd_var1( integer m_u2_A2,
                                      integer n_a1t,
                                      double* tau,
                                      double* u2, integer inc_u2,
                                      double* a1t, integer inc_a1t,
                                      double* A2, integer rs_A2, integer cs_A2 );

FLA_Error FLA_Apply_H2_UT_l_opc_var1( integer m_u2_A2,
                                      integer n_a1t,
                                      scomplex* tau,
                                      scomplex* u2, integer inc_u2,
                                      scomplex* a1t, integer inc_a1t,
                                      scomplex* A2, integer rs_A2, integer cs_A2 );

FLA_Error FLA_Apply_H2_UT_l_opz_var1( integer m_u2_A2,
                                      integer n_a1t,
                                      dcomplex* tau,
                                      dcomplex* u2, integer inc_u2,
                                      dcomplex* a1t, integer inc_a1t,
                                      dcomplex* A2, integer rs_A2, integer cs_A2 );

