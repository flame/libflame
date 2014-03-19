/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_LQ_UT_unb_var1( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_LQ_UT_blk_var1( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl );
FLA_Error FLA_LQ_UT_opt_var1( FLA_Obj A, FLA_Obj t );
FLA_Error FLA_LQ_UT_ops_var1( int m_A,
                              int n_A,
                              float* A, int rs_A, int cs_A,
                              float* t, int inc_t );
FLA_Error FLA_LQ_UT_opd_var1( int m_A,
                              int n_A,
                              double* A, int rs_A, int cs_A,
                              double* t, int inc_t );
FLA_Error FLA_LQ_UT_opc_var1( int m_A,
                              int n_A,
                              scomplex* A, int rs_A, int cs_A,
                              scomplex* t, int inc_t );
FLA_Error FLA_LQ_UT_opz_var1( int m_A,
                              int n_A,
                              dcomplex* A, int rs_A, int cs_A,
                              dcomplex* t, int inc_t );

FLA_Error FLA_LQ_UT_unb_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_LQ_UT_blk_var2( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl );
FLA_Error FLA_LQ_UT_opt_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_LQ_UT_ops_var2( int m_A,
                              int n_A,
                              float* A, int rs_A, int cs_A,
                              float* T, int rs_T, int cs_T );
FLA_Error FLA_LQ_UT_opd_var2( int m_A,
                              int n_A,
                              double* A, int rs_A, int cs_A,
                              double* T, int rs_T, int cs_T );
FLA_Error FLA_LQ_UT_opc_var2( int m_A,
                              int n_A,
                              scomplex* A, int rs_A, int cs_A,
                              scomplex* T, int rs_T, int cs_T );
FLA_Error FLA_LQ_UT_opz_var2( int m_A,
                              int n_A,
                              dcomplex* A, int rs_A, int cs_A,
                              dcomplex* T, int rs_T, int cs_T );

FLA_Error FLA_LQ_UT_blk_var3( FLA_Obj A, FLA_Obj T, fla_lqut_t* cntl );

