/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_LU_piv_blk_var3( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl );
FLA_Error FLA_LU_piv_blk_var4( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl );
FLA_Error FLA_LU_piv_blk_var5( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl );

FLA_Error FLA_LU_piv_unb_var3( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_LU_piv_unb_var3b( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_LU_piv_unb_var4( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_LU_piv_unb_var5( FLA_Obj A, FLA_Obj p );

FLA_Error FLA_LU_piv_opt_var3( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_LU_piv_ops_var3( integer m_A,
                               integer n_A,
                               float*    buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );
FLA_Error FLA_LU_piv_opd_var3( integer m_A,
                               integer n_A,
                               double*   buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );
FLA_Error FLA_LU_piv_opc_var3( integer m_A,
                               integer n_A,
                               scomplex* buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );
FLA_Error FLA_LU_piv_opz_var3( integer m_A,
                               integer n_A,
                               dcomplex* buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );

FLA_Error FLA_LU_piv_opt_var4( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_LU_piv_ops_var4( integer m_A,
                               integer n_A,
                               float*    buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );
FLA_Error FLA_LU_piv_opd_var4( integer m_A,
                               integer n_A,
                               double*   buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );
FLA_Error FLA_LU_piv_opc_var4( integer m_A,
                               integer n_A,
                               scomplex* buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );
FLA_Error FLA_LU_piv_opz_var4( integer m_A,
                               integer n_A,
                               dcomplex* buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );

FLA_Error FLA_LU_piv_opt_var5( FLA_Obj A, FLA_Obj p );
FLA_Error FLA_LU_piv_ops_var5( integer m_A,
                               integer n_A,
                               float*    buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );
FLA_Error FLA_LU_piv_opd_var5( integer m_A,
                               integer n_A,
                               double*   buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );
FLA_Error FLA_LU_piv_opc_var5( integer m_A,
                               integer n_A,
                               scomplex* buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );
FLA_Error FLA_LU_piv_opz_var5( integer m_A,
                               integer n_A,
                               dcomplex* buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p );
