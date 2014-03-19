/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_Apply_G_rb_opt_var1( FLA_Obj c, FLA_Obj s, FLA_Obj A );
FLA_Error FLA_Apply_G_rb_ops_var1( int       m_A,
                                   int       n_A,
                                   float*    buff_c, int inc_c,
                                   float*    buff_s, int inc_s,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rb_opd_var1( int       m_A,
                                   int       n_A,
                                   double*   buff_c, int inc_c,
                                   double*   buff_s, int inc_s,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rb_opc_var1( int       m_A,
                                   int       n_A,
                                   float*    buff_c, int inc_c,
                                   float*    buff_s, int inc_s,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rb_opz_var1( int       m_A,
                                   int       n_A,
                                   double*   buff_c, int inc_c,
                                   double*   buff_s, int inc_s,
                                   dcomplex* buff_A, int rs_A, int cs_A );

