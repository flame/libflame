/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Tevd_iteracc_n_opt_var1() -------------------------------------------

FLA_Error FLA_Tevd_iteracc_n_ops_var1( int       m_A,
                                       int       n_G,
                                       int       ijTL,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e,
                                       int*      n_iter_perf );
FLA_Error FLA_Tevd_iteracc_n_opd_var1( int       m_A,
                                       int       n_G,
                                       int       ijTL,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e,
                                       int*      n_iter_perf );

