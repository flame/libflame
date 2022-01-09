/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
   Copyright (c) 2021 Advanced Micro Devices, Inc. All rights reserved.
* */

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
FLA_Error FLA_LU_nopiv_ops_var1( integer m_A,
                                 integer n_A,
                                 float* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opd_var1( integer m_A,
                                 integer n_A,
                                 double* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opc_var1( integer m_A,
                                 integer n_A,
                                 scomplex* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opz_var1( integer m_A,
                                 integer n_A,
                                 dcomplex* A, integer rs_A, integer cs_A );

FLA_Error FLA_LU_nopiv_opt_var2( FLA_Obj A );
FLA_Error FLA_LU_nopiv_ops_var2( integer m_A,
                                 integer n_A,
                                 float* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opd_var2( integer m_A,
                                 integer n_A,
                                 double* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opc_var2( integer m_A,
                                 integer n_A,
                                 scomplex* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opz_var2( integer m_A,
                                 integer n_A,
                                 dcomplex* A, integer rs_A, integer cs_A );

FLA_Error FLA_LU_nopiv_opt_var3( FLA_Obj A );
FLA_Error FLA_LU_nopiv_ops_var3( integer m_A,
                                 integer n_A,
                                 float* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opd_var3( integer m_A,
                                 integer n_A,
                                 double* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opc_var3( integer m_A,
                                 integer n_A,
                                 scomplex* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opz_var3( integer m_A,
                                 integer n_A,
                                 dcomplex* A, integer rs_A, integer cs_A );

FLA_Error FLA_LU_nopiv_opt_var4( FLA_Obj A );
FLA_Error FLA_LU_nopiv_ops_var4( integer m_A,
                                 integer n_A,
                                 float* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opd_var4( integer m_A,
                                 integer n_A,
                                 double* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opc_var4( integer m_A,
                                 integer n_A,
                                 scomplex* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opz_var4( integer m_A,
                                 integer n_A,
                                 dcomplex* A, integer rs_A, integer cs_A );

FLA_Error FLA_LU_nopiv_opt_var5( FLA_Obj A );
FLA_Error FLA_LU_nopiv_ops_var5( integer m_A,
                                 integer n_A,
                                 float* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opd_var5( integer m_A,
                                 integer n_A,
                                 double* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opc_var5( integer m_A,
                                 integer n_A,
                                 scomplex* A, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_opz_var5( integer m_A,
                                 integer n_A,
                                 dcomplex* A, integer rs_A, integer cs_A );

FLA_Error FLA_LU_nopiv_is_blk_var1( integer m_A, integer n_A,FLA_Obj A, float* buff_A, integer nfact, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_id_blk_var1( integer m_A, integer n_A,FLA_Obj A, double* buff_A, integer nfact, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_ic_blk_var1( integer m_A, integer n_A,FLA_Obj A, scomplex* buff_A, integer nfact, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_iz_blk_var1( integer m_A, integer n_A,FLA_Obj A, dcomplex* buff_A, integer nfact, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_is_unblk_var1( integer m_A,
                                      integer n_A,
                                      float* A, integer nfact, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_id_unblk_var1( integer m_A,
                                      integer n_A,
                                      double* A, integer nfact, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_ic_unblk_var1( integer m_A,
                                      integer n_A,
                                      scomplex* A, integer nfact, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_iz_unblk_var1( integer m_A,
                                      integer n_A,
                                      dcomplex* A, integer nfact, integer rs_A, integer cs_A );
FLA_Error FLA_LU_nopiv_id_unblk_var2( integer m_A,
                                      integer n_A,
                                      double* A, integer nfact, integer rs_A, integer cs_A );
