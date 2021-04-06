/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Tevd_eigval_v_opt_var1() --------------------------------------------

FLA_Error FLA_Tevd_eigval_v_opt_var1( FLA_Obj G, FLA_Obj d, FLA_Obj e, FLA_Obj n_iter );
FLA_Error FLA_Tevd_eigval_v_ops_var1( integer       m_A,
                                      integer       n_G,
                                      scomplex* buff_G, integer rs_G, integer cs_G,
                                      float*    buff_d, integer inc_d, 
                                      float*    buff_e, integer inc_e,
                                      integer*      n_iter );
FLA_Error FLA_Tevd_eigval_v_opd_var1( integer       m_A,
                                      integer       n_G,
                                      dcomplex* buff_G, integer rs_G, integer cs_G,
                                      double*   buff_d, integer inc_d, 
                                      double*   buff_e, integer inc_e,
                                      integer*      n_iter );

FLA_Error FLA_Tevd_eigval_v_ops_var3( integer       m_A,
                                      integer       m_U,
                                      integer       n_G,
                                      scomplex* buff_G, integer rs_G, integer cs_G,
                                      float*    buff_d, integer inc_d, 
                                      float*    buff_e, integer inc_e,
                                      float*    buff_l, integer inc_l,
                                      integer*      buff_ls, integer inc_ls,
                                      float*    buff_pu, integer inc_pu,
                                      integer*      n_iter );
FLA_Error FLA_Tevd_eigval_v_opd_var3( integer       m_A,
                                      integer       m_U,
                                      integer       n_G,
                                      dcomplex* buff_G, integer rs_G, integer cs_G,
                                      double*   buff_d, integer inc_d, 
                                      double*   buff_e, integer inc_e,
                                      double*   buff_l, integer inc_l,
                                      integer*      buff_ls, integer inc_ls,
                                      double*   buff_pu, integer inc_pu,
                                      integer*      n_iter );

