/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Sylv_nn_blk_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var2( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var3( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var4( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var5( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var6( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var7( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var8( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var9( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var10( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var11( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var12( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var13( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var14( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var15( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var16( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var17( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );
FLA_Error FLA_Sylv_nn_blk_var18( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale, fla_sylv_t* cntl );

FLA_Error FLA_Sylv_nn_opt_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var2( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var3( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var4( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var5( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var6( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var7( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var8( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var9( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var10( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var11( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var12( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var13( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var14( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var15( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var16( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var17( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );
FLA_Error FLA_Sylv_nn_opt_var18( FLA_Obj isgn, FLA_Obj A, FLA_Obj B, FLA_Obj C, FLA_Obj scale );

FLA_Error FLA_Sylv_nn_ops_var1( float sgn,
                                integer m_C,
                                integer n_C,
                                float* buff_A, integer rs_A, integer cs_A,
                                float* buff_B, integer rs_B, integer cs_B,
                                float* buff_C, integer rs_C, integer cs_C,
                                float* buff_scale,
                                integer* info );
FLA_Error FLA_Sylv_nn_opd_var1( double sgn,
                                integer m_C,
                                integer n_C,
                                double* buff_A, integer rs_A, integer cs_A,
                                double* buff_B, integer rs_B, integer cs_B,
                                double* buff_C, integer rs_C, integer cs_C,
                                double* buff_scale,
                                integer* info );
FLA_Error FLA_Sylv_nn_opc_var1( float sgn,
                                integer m_C,
                                integer n_C,
                                scomplex* buff_A, integer rs_A, integer cs_A,
                                scomplex* buff_B, integer rs_B, integer cs_B,
                                scomplex* buff_C, integer rs_C, integer cs_C,
                                scomplex* buff_scale,
                                integer* info );
FLA_Error FLA_Sylv_nn_opz_var1( double sgn,
                                integer m_C,
                                integer n_C,
                                dcomplex* buff_A, integer rs_A, integer cs_A,
                                dcomplex* buff_B, integer rs_B, integer cs_B,
                                dcomplex* buff_C, integer rs_C, integer cs_C,
                                dcomplex* buff_scale,
                                integer* info );
