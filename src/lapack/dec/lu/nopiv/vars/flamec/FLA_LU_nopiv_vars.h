/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_LU_nopiv_blk_var1( FLA_Obj A, fla_lu_t* cntl );
FLA_Error FLA_LU_nopiv_blk_var2( FLA_Obj A, fla_lu_t* cntl );
FLA_Error FLA_LU_nopiv_blk_var3( FLA_Obj A, fla_lu_t* cntl );
FLA_Error FLA_LU_nopiv_blk_var4( FLA_Obj A, fla_lu_t* cntl );
FLA_Error FLA_LU_nopiv_blk_var5( FLA_Obj A, fla_lu_t* cntl );

FLA_Error FLA_LU_nopiv_unb_var1( FLA_Obj A );
FLA_Error FLA_LU_nopiv_unb_var2( FLA_Obj A );
FLA_Error FLA_LU_nopiv_unb_var3( FLA_Obj A );
FLA_Error FLA_LU_nopiv_unb_var4( FLA_Obj A );
FLA_Error FLA_LU_nopiv_unb_var5( FLA_Obj A );

FLA_Error FLA_LU_nopiv_opt_var1( FLA_Obj A );
FLA_Error FLA_LU_nopiv_ops_var1( int m_A,
                                 int n_A,
                                 float* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opd_var1( int m_A,
                                 int n_A,
                                 double* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opc_var1( int m_A,
                                 int n_A,
                                 scomplex* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opz_var1( int m_A,
                                 int n_A,
                                 dcomplex* A, int rs_A, int cs_A );

FLA_Error FLA_LU_nopiv_opt_var2( FLA_Obj A );
FLA_Error FLA_LU_nopiv_ops_var2( int m_A,
                                 int n_A,
                                 float* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opd_var2( int m_A,
                                 int n_A,
                                 double* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opc_var2( int m_A,
                                 int n_A,
                                 scomplex* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opz_var2( int m_A,
                                 int n_A,
                                 dcomplex* A, int rs_A, int cs_A );

FLA_Error FLA_LU_nopiv_opt_var3( FLA_Obj A );
FLA_Error FLA_LU_nopiv_ops_var3( int m_A,
                                 int n_A,
                                 float* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opd_var3( int m_A,
                                 int n_A,
                                 double* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opc_var3( int m_A,
                                 int n_A,
                                 scomplex* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opz_var3( int m_A,
                                 int n_A,
                                 dcomplex* A, int rs_A, int cs_A );

FLA_Error FLA_LU_nopiv_opt_var4( FLA_Obj A );
FLA_Error FLA_LU_nopiv_ops_var4( int m_A,
                                 int n_A,
                                 float* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opd_var4( int m_A,
                                 int n_A,
                                 double* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opc_var4( int m_A,
                                 int n_A,
                                 scomplex* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opz_var4( int m_A,
                                 int n_A,
                                 dcomplex* A, int rs_A, int cs_A );

FLA_Error FLA_LU_nopiv_opt_var5( FLA_Obj A );
FLA_Error FLA_LU_nopiv_ops_var5( int m_A,
                                 int n_A,
                                 float* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opd_var5( int m_A,
                                 int n_A,
                                 double* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opc_var5( int m_A,
                                 int n_A,
                                 scomplex* A, int rs_A, int cs_A );
FLA_Error FLA_LU_nopiv_opz_var5( int m_A,
                                 int n_A,
                                 dcomplex* A, int rs_A, int cs_A );

