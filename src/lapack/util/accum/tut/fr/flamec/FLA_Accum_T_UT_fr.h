/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Accum_T_UT_fr_unb_var1( FLA_Obj A, FLA_Obj t, FLA_Obj T );
FLA_Error FLA_Accum_T_UT_fr_blk_var2( FLA_Obj A, FLA_Obj t, FLA_Obj T );

FLA_Error FLA_Accum_T_UT_fr_opt_var1( FLA_Obj A, FLA_Obj t, FLA_Obj T );

FLA_Error FLA_Accum_T_UT_fr_ops_var1( int m_A,
                                      int n_A,
                                      float* A, int rs_A, int cs_A,
                                      int m_t,
                                      float* t, int inc_t,
                                      float* T, int rs_T, int cs_T );
FLA_Error FLA_Accum_T_UT_fr_opd_var1( int m_A,
                                      int n_A,
                                      double* A, int rs_A, int cs_A,
                                      int m_t,
                                      double* t, int inc_t,
                                      double* T, int rs_T, int cs_T );
FLA_Error FLA_Accum_T_UT_fr_opc_var1( int m_A,
                                      int n_A,
                                      scomplex* A, int rs_A, int cs_A,
                                      int m_t,
                                      scomplex* t, int inc_t,
                                      scomplex* T, int rs_T, int cs_T );
FLA_Error FLA_Accum_T_UT_fr_opz_var1( int m_A,
                                      int n_A,
                                      dcomplex* A, int rs_A, int cs_A,
                                      int m_t,
                                      dcomplex* t, int inc_t,
                                      dcomplex* T, int rs_T, int cs_T );
