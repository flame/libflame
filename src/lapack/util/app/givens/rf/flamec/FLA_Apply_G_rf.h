/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// Variant 1

FLA_Error FLA_Apply_G_rf_opt_var1( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ops_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opd_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opc_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opz_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_asm_var1( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var1( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var1( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );

// Variant 2

FLA_Error FLA_Apply_G_rf_opt_var2( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ops_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opd_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opc_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opz_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_asm_var2( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var2( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var2( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );

// Variant 3

FLA_Error FLA_Apply_G_rf_opt_var3( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ops_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opd_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opc_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opz_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_asm_var3( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var3( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );

// Variant 4

FLA_Error FLA_Apply_G_rf_opt_var4( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ops_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opd_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opc_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opz_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_asm_var4( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var4( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var4( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );

// Variant 5

FLA_Error FLA_Apply_G_rf_opt_var5( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ops_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opd_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opc_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opz_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_asm_var5( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var5( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var5( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );

// Variant 6

FLA_Error FLA_Apply_G_rf_opt_var6( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ops_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opd_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opc_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opz_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_asm_var6( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var6( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var6( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );

// Variant 7

FLA_Error FLA_Apply_G_rf_opt_var7( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ops_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opd_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opc_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opz_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_asm_var7( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var7( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var7( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );

// Variant 8

FLA_Error FLA_Apply_G_rf_opt_var8( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ops_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opd_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opc_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opz_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_asm_var8( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var8( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var8( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );

// Variant 9

FLA_Error FLA_Apply_G_rf_opt_var9( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ops_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opd_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opc_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_opz_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_asm_var9( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var9( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var9( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );







// Variant 3b

FLA_Error FLA_Apply_G_rf_asm_var3b( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var3b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var3b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var3b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var3b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var3b( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var3b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var3b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var3b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var3b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );


// Variant 5b

FLA_Error FLA_Apply_G_rf_asm_var5b( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var5b( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var5b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );


// Variant 6b

FLA_Error FLA_Apply_G_rf_asm_var6b( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var6b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var6b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var6b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var6b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var6b( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var6b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var6b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var6b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var6b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );


// Variant 8b

FLA_Error FLA_Apply_G_rf_asm_var8b( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var8b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var8b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var8b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var8b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var8b( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var8b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var8b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var8b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var8b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );

FLA_Error FLA_Apply_G_rf_bhs_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bhd_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bhc_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bhz_var3( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   FLA_Obj*  buff_A, int rs_A, int cs_A,
                                   int       b_alg );


// Variant 9b

FLA_Error FLA_Apply_G_rf_asm_var9b( FLA_Obj G, FLA_Obj A );
FLA_Error FLA_Apply_G_rf_ass_var9b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asd_var9b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asc_var9b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Apply_G_rf_asz_var9b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   int       iTL,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Apply_G_rf_blk_var9b( FLA_Obj G, FLA_Obj A, dim_t b_alg );
FLA_Error FLA_Apply_G_rf_bls_var9b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   float*    buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_bld_var9b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   double*   buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blc_var9b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   scomplex* buff_G, int rs_G, int cs_G,
                                   scomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );
FLA_Error FLA_Apply_G_rf_blz_var9b( int       k_G,
                                   int       m_A,
                                   int       n_A,
                                   int       i_k,
                                   dcomplex* buff_G, int rs_G, int cs_G,
                                   dcomplex* buff_A, int rs_A, int cs_A,
                                   int       b_alg );


