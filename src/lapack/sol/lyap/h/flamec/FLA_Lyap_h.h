/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Lyap_h_unb_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj C );
FLA_Error FLA_Lyap_h_unb_var2( FLA_Obj isgn, FLA_Obj A, FLA_Obj C );
FLA_Error FLA_Lyap_h_unb_var3( FLA_Obj isgn, FLA_Obj A, FLA_Obj C );
FLA_Error FLA_Lyap_h_unb_var4( FLA_Obj isgn, FLA_Obj A, FLA_Obj C );

FLA_Error FLA_Lyap_h_blk_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl );
FLA_Error FLA_Lyap_h_blk_var2( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl );
FLA_Error FLA_Lyap_h_blk_var3( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl );
FLA_Error FLA_Lyap_h_blk_var4( FLA_Obj isgn, FLA_Obj A, FLA_Obj C, FLA_Obj scale, fla_lyap_t* cntl );

FLA_Error FLA_Lyap_h_opt_var1( FLA_Obj isgn, FLA_Obj A, FLA_Obj C );
FLA_Error FLA_Lyap_h_ops_var1( integer m_AC,
                               float* buff_sgn,
                               float* buff_A, integer rs_A, integer cs_A, 
                               float* buff_W, integer rs_W, integer cs_W, 
                               float* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opd_var1( integer m_AC,
                               double* buff_sgn,
                               double* buff_A, integer rs_A, integer cs_A, 
                               double* buff_W, integer rs_W, integer cs_W, 
                               double* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opc_var1( integer m_AC,
                               scomplex* buff_sgn,
                               scomplex* buff_A, integer rs_A, integer cs_A, 
                               scomplex* buff_W, integer rs_W, integer cs_W, 
                               scomplex* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opz_var1( integer m_AC,
                               dcomplex* buff_sgn,
                               dcomplex* buff_A, integer rs_A, integer cs_A, 
                               dcomplex* buff_W, integer rs_W, integer cs_W, 
                               dcomplex* buff_C, integer rs_C, integer cs_C );

FLA_Error FLA_Lyap_h_opt_var2( FLA_Obj isgn, FLA_Obj A, FLA_Obj C );
FLA_Error FLA_Lyap_h_ops_var2( integer m_AC,
                               float* buff_sgn,
                               float* buff_A, integer rs_A, integer cs_A, 
                               float* buff_W, integer rs_W, integer cs_W, 
                               float* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opd_var2( integer m_AC,
                               double* buff_sgn,
                               double* buff_A, integer rs_A, integer cs_A, 
                               double* buff_W, integer rs_W, integer cs_W, 
                               double* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opc_var2( integer m_AC,
                               scomplex* buff_sgn,
                               scomplex* buff_A, integer rs_A, integer cs_A, 
                               scomplex* buff_W, integer rs_W, integer cs_W, 
                               scomplex* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opz_var2( integer m_AC,
                               dcomplex* buff_sgn,
                               dcomplex* buff_A, integer rs_A, integer cs_A, 
                               dcomplex* buff_W, integer rs_W, integer cs_W, 
                               dcomplex* buff_C, integer rs_C, integer cs_C );

FLA_Error FLA_Lyap_h_opt_var3( FLA_Obj isgn, FLA_Obj A, FLA_Obj C );
FLA_Error FLA_Lyap_h_ops_var3( integer m_AC,
                               float* buff_sgn,
                               float* buff_A, integer rs_A, integer cs_A, 
                               float* buff_W, integer rs_W, integer cs_W, 
                               float* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opd_var3( integer m_AC,
                               double* buff_sgn,
                               double* buff_A, integer rs_A, integer cs_A, 
                               double* buff_W, integer rs_W, integer cs_W, 
                               double* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opc_var3( integer m_AC,
                               scomplex* buff_sgn,
                               scomplex* buff_A, integer rs_A, integer cs_A, 
                               scomplex* buff_W, integer rs_W, integer cs_W, 
                               scomplex* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opz_var3( integer m_AC,
                               dcomplex* buff_sgn,
                               dcomplex* buff_A, integer rs_A, integer cs_A, 
                               dcomplex* buff_W, integer rs_W, integer cs_W, 
                               dcomplex* buff_C, integer rs_C, integer cs_C );

FLA_Error FLA_Lyap_h_opt_var4( FLA_Obj isgn, FLA_Obj A, FLA_Obj C );
FLA_Error FLA_Lyap_h_ops_var4( integer m_AC,
                               float* buff_sgn,
                               float* buff_A, integer rs_A, integer cs_A, 
                               float* buff_W, integer rs_W, integer cs_W, 
                               float* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opd_var4( integer m_AC,
                               double* buff_sgn,
                               double* buff_A, integer rs_A, integer cs_A, 
                               double* buff_W, integer rs_W, integer cs_W, 
                               double* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opc_var4( integer m_AC,
                               scomplex* buff_sgn,
                               scomplex* buff_A, integer rs_A, integer cs_A, 
                               scomplex* buff_W, integer rs_W, integer cs_W, 
                               scomplex* buff_C, integer rs_C, integer cs_C );
FLA_Error FLA_Lyap_h_opz_var4( integer m_AC,
                               dcomplex* buff_sgn,
                               dcomplex* buff_A, integer rs_A, integer cs_A, 
                               dcomplex* buff_W, integer rs_W, integer cs_W, 
                               dcomplex* buff_C, integer rs_C, integer cs_C );
