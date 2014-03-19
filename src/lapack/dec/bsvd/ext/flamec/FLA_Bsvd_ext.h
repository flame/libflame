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
FLA_Error FLA_Bsvd_ext_ops_var1( int       m_d,
                                 int       m_U,
                                 int       m_V,
                                 int       m_C,
                                 int       n_C,
                                 int       n_GH,
                                 int       n_iter_max,
                                 float*    buff_d, int inc_d, 
                                 float*    buff_e, int inc_e,
                                 scomplex* buff_G, int rs_G, int cs_G,
                                 scomplex* buff_H, int rs_H, int cs_H,
                                 float*    buff_U, int rs_U, int cs_U,
                                 float*    buff_V, int rs_V, int cs_V,
                                 float*    buff_C, int rs_C, int cs_C,
                                 int       b_alg );
FLA_Error FLA_Bsvd_ext_opd_var1( int       m_d,
                                 int       m_U,
                                 int       m_V,
                                 int       m_C,
                                 int       n_C,
                                 int       n_GH,
                                 int       n_iter_max,
                                 double*   buff_d, int inc_d, 
                                 double*   buff_e, int inc_e,
                                 dcomplex* buff_G, int rs_G, int cs_G,
                                 dcomplex* buff_H, int rs_H, int cs_H,
                                 double*   buff_U, int rs_U, int cs_U,
                                 double*   buff_V, int rs_V, int cs_V,
                                 double*   buff_C, int rs_C, int cs_C,
                                 int       b_alg );
FLA_Error FLA_Bsvd_ext_opc_var1( int       m_d,
                                 int       m_U,
                                 int       m_V,
                                 int       m_C,
                                 int       n_C,
                                 int       n_GH,
                                 int       n_iter_max,
                                 float*    buff_d, int inc_d, 
                                 float*    buff_e, int inc_e,
                                 scomplex* buff_G, int rs_G, int cs_G,
                                 scomplex* buff_H, int rs_H, int cs_H,
                                 scomplex* buff_U, int rs_U, int cs_U,
                                 scomplex* buff_V, int rs_V, int cs_V,
                                 scomplex* buff_C, int rs_C, int cs_C,
                                 int       b_alg );
FLA_Error FLA_Bsvd_ext_opz_var1( int       m_d,
                                 int       m_U,
                                 int       m_V,
                                 int       m_C,
                                 int       n_C,
                                 int       n_GH,
                                 int       n_iter_max,
                                 double*   buff_d, int inc_d, 
                                 double*   buff_e, int inc_e,
                                 dcomplex* buff_G, int rs_G, int cs_G,
                                 dcomplex* buff_H, int rs_H, int cs_H,
                                 dcomplex* buff_U, int rs_U, int cs_U,
                                 dcomplex* buff_V, int rs_V, int cs_V,
                                 dcomplex* buff_C, int rs_C, int cs_C,
                                 int       b_alg );

