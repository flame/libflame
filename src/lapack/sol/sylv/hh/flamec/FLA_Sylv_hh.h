/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Sylv_hh_blk_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var2( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var3( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var4( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var5( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var6( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var7( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var8( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var9( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var10( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var11( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var12( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var13( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var14( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var15( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var16( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var17( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_hh_blk_var18( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );

FLA_Error FLA_Sylv_hh_opt_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var2( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var3( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var4( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var5( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var6( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var7( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var8( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var9( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var10( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var11( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var12( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var13( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var14( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var15( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var16( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var17( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_hh_opt_var18( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );

FLA_Error FLA_Sylv_hh_ops_var1( float sgn,
                                int m_C,
                                int n_C,
                                float* buff_A, int rs_A, int cs_A,
                                float* buff_B, int rs_B, int cs_B,
                                float* buff_C, int rs_C, int cs_C,
                                float* buff_scale,
                                int* info );
FLA_Error FLA_Sylv_hh_opd_var1( double sgn,
                                int m_C,
                                int n_C,
                                double* buff_A, int rs_A, int cs_A,
                                double* buff_B, int rs_B, int cs_B,
                                double* buff_C, int rs_C, int cs_C,
                                double* buff_scale,
                                int* info );
FLA_Error FLA_Sylv_hh_opc_var1( float sgn,
                                int m_C,
                                int n_C,
                                scomplex* buff_A, int rs_A, int cs_A,
                                scomplex* buff_B, int rs_B, int cs_B,
                                scomplex* buff_C, int rs_C, int cs_C,
                                scomplex* buff_scale,
                                int* info );
FLA_Error FLA_Sylv_hh_opz_var1( double sgn,
                                int m_C,
                                int n_C,
                                dcomplex* buff_A, int rs_A, int cs_A,
                                dcomplex* buff_B, int rs_B, int cs_B,
                                dcomplex* buff_C, int rs_C, int cs_C,
                                dcomplex* buff_scale,
                                int* info );
