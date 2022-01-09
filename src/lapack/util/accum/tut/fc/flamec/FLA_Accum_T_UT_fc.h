/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Accum_T_UT_fc_unb_var1( FLA_Obj A, FLA_Obj t, FLA_Obj T );
FLA_Error FLA_Accum_T_UT_fc_blk_var2( FLA_Obj A, FLA_Obj t, FLA_Obj T );

FLA_Error FLA_Accum_T_UT_fc_opt_var1( FLA_Obj A, FLA_Obj t, FLA_Obj T );

FLA_Error FLA_Accum_T_UT_fc_ops_var1( integer m_A,
                                      integer n_AT,
                                      float* A, integer rs_A, integer cs_A,
                                      integer m_t, 
                                      float* t, integer inc_t,
                                      float* T, integer rs_T, integer cs_T );
FLA_Error FLA_Accum_T_UT_fc_opd_var1( integer m_A,
                                      integer n_AT,
                                      double* A, integer rs_A, integer cs_A,
                                      integer m_t, 
                                      double* t, integer inc_t,
                                      double* T, integer rs_T, integer cs_T );
FLA_Error FLA_Accum_T_UT_fc_opc_var1( integer m_A,
                                      integer n_AT,
                                      scomplex* A, integer rs_A, integer cs_A,
                                      integer m_t, 
                                      scomplex* t, integer inc_t,
                                      scomplex* T, integer rs_T, integer cs_T );
FLA_Error FLA_Accum_T_UT_fc_opz_var1( integer m_A,
                                      integer n_AT,
                                      dcomplex* A, integer rs_A, integer cs_A,
                                      integer m_t, 
                                      dcomplex* t, integer inc_t,
                                      dcomplex* T, integer rs_T, integer cs_T );
