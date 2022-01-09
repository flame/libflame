/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Tevd_eigval_n_opt_var1() --------------------------------------------

FLA_Error FLA_Tevd_eigval_n_opt_var1( FLA_Obj G, FLA_Obj d, FLA_Obj e, FLA_Obj n_iter );
FLA_Error FLA_Tevd_eigval_n_ops_var1( integer       m_A,
                                      integer       n_G,
                                      float*    buff_d, integer inc_d, 
                                      float*    buff_e, integer inc_e,
                                      integer*      n_iter );
FLA_Error FLA_Tevd_eigval_n_opd_var1( integer       m_A,
                                      integer       n_G,
                                      double*   buff_d, integer inc_d, 
                                      double*   buff_e, integer inc_e,
                                      integer*      n_iter );

