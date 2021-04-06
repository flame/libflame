/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- FLA_Bsvd_ext_opt_var1() ---------------------------------------------------

FLA_Error FLA_Bsvd_ext_opt_var1( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H, 
                                 FLA_Svd_type jobu, FLA_Obj U, 
                                 FLA_Svd_type jobv, FLA_Obj V, 
                                 FLA_Bool apply_Uh2C, FLA_Obj C,
                                 dim_t b_alg );
FLA_Error FLA_Bsvd_ext_ops_var1( integer       m_d,
                                 integer       m_U,
                                 integer       m_V,
                                 integer       m_C,
                                 integer       n_C,
                                 integer       n_GH,
                                 integer       n_iter_max,
                                 float*    buff_d, integer inc_d, 
                                 float*    buff_e, integer inc_e,
                                 scomplex* buff_G, integer rs_G, integer cs_G,
                                 scomplex* buff_H, integer rs_H, integer cs_H,
                                 float*    buff_U, integer rs_U, integer cs_U,
                                 float*    buff_V, integer rs_V, integer cs_V,
                                 float*    buff_C, integer rs_C, integer cs_C,
                                 integer       b_alg );
FLA_Error FLA_Bsvd_ext_opd_var1( integer       m_d,
                                 integer       m_U,
                                 integer       m_V,
                                 integer       m_C,
                                 integer       n_C,
                                 integer       n_GH,
                                 integer       n_iter_max,
                                 double*   buff_d, integer inc_d, 
                                 double*   buff_e, integer inc_e,
                                 dcomplex* buff_G, integer rs_G, integer cs_G,
                                 dcomplex* buff_H, integer rs_H, integer cs_H,
                                 double*   buff_U, integer rs_U, integer cs_U,
                                 double*   buff_V, integer rs_V, integer cs_V,
                                 double*   buff_C, integer rs_C, integer cs_C,
                                 integer       b_alg );
FLA_Error FLA_Bsvd_ext_opc_var1( integer       m_d,
                                 integer       m_U,
                                 integer       m_V,
                                 integer       m_C,
                                 integer       n_C,
                                 integer       n_GH,
                                 integer       n_iter_max,
                                 float*    buff_d, integer inc_d, 
                                 float*    buff_e, integer inc_e,
                                 scomplex* buff_G, integer rs_G, integer cs_G,
                                 scomplex* buff_H, integer rs_H, integer cs_H,
                                 scomplex* buff_U, integer rs_U, integer cs_U,
                                 scomplex* buff_V, integer rs_V, integer cs_V,
                                 scomplex* buff_C, integer rs_C, integer cs_C,
                                 integer       b_alg );
FLA_Error FLA_Bsvd_ext_opz_var1( integer       m_d,
                                 integer       m_U,
                                 integer       m_V,
                                 integer       m_C,
                                 integer       n_C,
                                 integer       n_GH,
                                 integer       n_iter_max,
                                 double*   buff_d, integer inc_d, 
                                 double*   buff_e, integer inc_e,
                                 dcomplex* buff_G, integer rs_G, integer cs_G,
                                 dcomplex* buff_H, integer rs_H, integer cs_H,
                                 dcomplex* buff_U, integer rs_U, integer cs_U,
                                 dcomplex* buff_V, integer rs_V, integer cs_V,
                                 dcomplex* buff_C, integer rs_C, integer cs_C,
                                 integer       b_alg );

