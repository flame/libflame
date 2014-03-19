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

FLA_Error FLA_Tevd_find_submatrix_ops( int       m_A,
                                       int       ij_begin,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR );
FLA_Error FLA_Tevd_find_submatrix_opd( int       m_A,
                                       int       ij_begin,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR );

// --- FLA_Norm1_tridiag() -----------------------------------------------------

FLA_Error FLA_Norm1_tridiag( FLA_Obj d, FLA_Obj e, FLA_Obj norm );
FLA_Error FLA_Norm1_tridiag_ops( int       m_A,
                                 float*    buff_d, int inc_d, 
                                 float*    buff_e, int inc_e,
                                 float*    norm );
FLA_Error FLA_Norm1_tridiag_opd( int       m_A,
                                 double*   buff_d, int inc_d, 
                                 double*   buff_e, int inc_e,
                                 double*   norm );

// --- FLA_Tevd_n_opt_var1() ---------------------------------------------------

FLA_Error FLA_Tevd_n_opt_var1( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj U );
FLA_Error FLA_Tevd_n_ops_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G );
FLA_Error FLA_Tevd_n_opd_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G );
FLA_Error FLA_Tevd_n_opc_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               float*    buff_d, int inc_d, 
                               float*    buff_e, int inc_e,
                               scomplex* buff_G, int rs_G, int cs_G );
FLA_Error FLA_Tevd_n_opz_var1( int       m_A,
                               int       m_U,
                               int       n_G,
                               int       n_iter_max,
                               double*   buff_d, int inc_d, 
                               double*   buff_e, int inc_e,
                               dcomplex* buff_G, int rs_G, int cs_G );


