/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_Apply_G_lb_opt_var1( FLA_Obj c, FLA_Obj s, FLA_Obj A );
FLA_Error FLA_Apply_G_lb_ops_var1( integer       m_A,
                                   integer       n_A,
                                   float*    buff_c, integer inc_c,
                                   float*    buff_s, integer inc_s,
                                   float*    buff_A, integer rs_A, integer cs_A );
FLA_Error FLA_Apply_G_lb_opd_var1( integer       m_A,
                                   integer       n_A,
                                   double*   buff_c, integer inc_c,
                                   double*   buff_s, integer inc_s,
                                   double*   buff_A, integer rs_A, integer cs_A );
FLA_Error FLA_Apply_G_lb_opc_var1( integer       m_A,
                                   integer       n_A,
                                   float*    buff_c, integer inc_c,
                                   float*    buff_s, integer inc_s,
                                   scomplex* buff_A, integer rs_A, integer cs_A );
FLA_Error FLA_Apply_G_lb_opz_var1( integer       m_A,
                                   integer       n_A,
                                   double*   buff_c, integer inc_c,
                                   double*   buff_s, integer inc_s,
                                   dcomplex* buff_A, integer rs_A, integer cs_A );

