/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Tevd_francis_n_opt_var1() -------------------------------------------

FLA_Error FLA_Tevd_francis_n_opt_var1( FLA_Obj shift, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Tevd_francis_n_ops_var1( int       m_A,
                                       float*    buff_shift,
                                       float*    buff_d, int inc_d, 
                                       float*    buff_e, int inc_e ); 
FLA_Error FLA_Tevd_francis_n_opd_var1( int       m_A,
                                       double*   buff_shift,
                                       double*   buff_d, int inc_d, 
                                       double*   buff_e, int inc_e ); 

