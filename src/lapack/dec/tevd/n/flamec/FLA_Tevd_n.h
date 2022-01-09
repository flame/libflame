/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Tevd_iteracc_n.h"
#include "FLA_Tevd_eigval_n.h"
#include "FLA_Tevd_francis_n.h"

// --- FLA_Tevd_find_submatrix() -----------------------------------------------

FLA_Error FLA_Tevd_find_submatrix_ops( integer       m_A,
                                       integer       ij_begin,
                                       float*    buff_d, integer inc_d, 
                                       float*    buff_e, integer inc_e,
                                       integer*      ijTL,
                                       integer*      ijBR );
FLA_Error FLA_Tevd_find_submatrix_opd( integer       m_A,
                                       integer       ij_begin,
                                       double*   buff_d, integer inc_d, 
                                       double*   buff_e, integer inc_e,
                                       integer*      ijTL,
                                       integer*      ijBR );

// --- FLA_Norm1_tridiag() -----------------------------------------------------

FLA_Error FLA_Norm1_tridiag( FLA_Obj d, FLA_Obj e, FLA_Obj norm );
FLA_Error FLA_Norm1_tridiag_ops( integer       m_A,
                                 float*    buff_d, integer inc_d, 
                                 float*    buff_e, integer inc_e,
                                 float*    norm );
FLA_Error FLA_Norm1_tridiag_opd( integer       m_A,
                                 double*   buff_d, integer inc_d, 
                                 double*   buff_e, integer inc_e,
                                 double*   norm );

// --- FLA_Tevd_n_opt_var1() ---------------------------------------------------

FLA_Error FLA_Tevd_n_opt_var1( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U );
FLA_Error FLA_Tevd_n_ops_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G );
FLA_Error FLA_Tevd_n_opd_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G );
FLA_Error FLA_Tevd_n_opc_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d, 
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G );
FLA_Error FLA_Tevd_n_opz_var1( integer       m_A,
                               integer       m_U,
                               integer       n_G,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d, 
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G );


