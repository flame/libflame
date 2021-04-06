/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_Hess_UT_blk_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_unb_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_unb_var1( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Hess_UT_blk_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_blf_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_unb_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_unb_var2( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Hess_UT_blk_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_blf_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_unb_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_unb_var3( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Hess_UT_blk_var4( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_blf_var4( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_unb_var4( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_unb_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T );

FLA_Error FLA_Hess_UT_blk_var5( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_unb_var5( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_unb_var5( FLA_Obj A, FLA_Obj U, FLA_Obj Z, FLA_Obj T );


FLA_Error FLA_Hess_UT_opt_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_opt_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ops_var1( integer m_A,
                                     integer m_T,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opd_var1( integer m_A,
                                     integer m_T,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opc_var1( integer m_A,
                                     integer m_T,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opz_var1( integer m_A,
                                     integer m_T,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_T, integer rs_T, integer cs_T );


FLA_Error FLA_Hess_UT_opt_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_opt_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ops_var2( integer m_A,
                                     integer m_T,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opd_var2( integer m_A,
                                     integer m_T,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opc_var2( integer m_A,
                                     integer m_T,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opz_var2( integer m_A,
                                     integer m_T,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_T, integer rs_T, integer cs_T );


FLA_Error FLA_Hess_UT_opt_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_opt_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ops_var3( integer m_A,
                                     integer m_T,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opd_var3( integer m_A,
                                     integer m_T,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opc_var3( integer m_A,
                                     integer m_T,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opz_var3( integer m_A,
                                     integer m_T,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_T, integer rs_T, integer cs_T );


FLA_Error FLA_Hess_UT_opt_var4( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_opt_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ops_var4( integer m_A,
                                     integer m_T,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_Y, integer rs_Y, integer cs_Y, 
                                     float* buff_Z, integer rs_Z, integer cs_Z, 
                                     float* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opd_var4( integer m_A,
                                     integer m_T,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_Y, integer rs_Y, integer cs_Y, 
                                     double* buff_Z, integer rs_Z, integer cs_Z, 
                                     double* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opc_var4( integer m_A,
                                     integer m_T,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_Y, integer rs_Y, integer cs_Y, 
                                     scomplex* buff_Z, integer rs_Z, integer cs_Z, 
                                     scomplex* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opz_var4( integer m_A,
                                     integer m_T,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_Y, integer rs_Y, integer cs_Y, 
                                     dcomplex* buff_Z, integer rs_Z, integer cs_Z, 
                                     dcomplex* buff_T, integer rs_T, integer cs_T );


FLA_Error FLA_Hess_UT_opt_var5( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_opt_var5( FLA_Obj A, FLA_Obj U, FLA_Obj Z, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ops_var5( integer m_A,
                                     integer m_T,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_U, integer rs_U, integer cs_U, 
                                     float* buff_Z, integer rs_Z, integer cs_Z, 
                                     float* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opd_var5( integer m_A,
                                     integer m_T,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_U, integer rs_U, integer cs_U, 
                                     double* buff_Z, integer rs_Z, integer cs_Z, 
                                     double* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opc_var5( integer m_A,
                                     integer m_T,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_U, integer rs_U, integer cs_U, 
                                     scomplex* buff_Z, integer rs_Z, integer cs_Z, 
                                     scomplex* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_opz_var5( integer m_A,
                                     integer m_T,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_U, integer rs_U, integer cs_U, 
                                     dcomplex* buff_Z, integer rs_Z, integer cs_Z, 
                                     dcomplex* buff_T, integer rs_T, integer cs_T );


FLA_Error FLA_Hess_UT_ofu_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ofu_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ofs_var1( integer m_A,
                                     integer m_T,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofd_var1( integer m_A,
                                     integer m_T,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofc_var1( integer m_A,
                                     integer m_T,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofz_var1( integer m_A,
                                     integer m_T,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_T, integer rs_T, integer cs_T );


FLA_Error FLA_Hess_UT_ofu_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ofu_var2( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ofs_var2( integer m_A,
                                     integer m_T,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofd_var2( integer m_A,
                                     integer m_T,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofc_var2( integer m_A,
                                     integer m_T,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofz_var2( integer m_A,
                                     integer m_T,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_T, integer rs_T, integer cs_T );


FLA_Error FLA_Hess_UT_ofu_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ofu_var3( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ofs_var3( integer m_A,
                                     integer m_T,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofd_var3( integer m_A,
                                     integer m_T,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofc_var3( integer m_A,
                                     integer m_T,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofz_var3( integer m_A,
                                     integer m_T,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_T, integer rs_T, integer cs_T );


FLA_Error FLA_Hess_UT_ofu_var4( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ofu_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T );
FLA_Error FLA_Hess_UT_step_ofs_var4( integer m_A,
                                     integer m_T,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_Y, integer rs_Y, integer cs_Y,
                                     float* buff_Z, integer rs_Z, integer cs_Z,
                                     float* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofd_var4( integer m_A,
                                     integer m_T,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_Y, integer rs_Y, integer cs_Y,
                                     double* buff_Z, integer rs_Z, integer cs_Z,
                                     double* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofc_var4( integer m_A,
                                     integer m_T,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_Y, integer rs_Y, integer cs_Y,
                                     scomplex* buff_Z, integer rs_Z, integer cs_Z,
                                     scomplex* buff_T, integer rs_T, integer cs_T );
FLA_Error FLA_Hess_UT_step_ofz_var4( integer m_A,
                                     integer m_T,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_Y, integer rs_Y, integer cs_Y,
                                     dcomplex* buff_Z, integer rs_Z, integer cs_Z,
                                     dcomplex* buff_T, integer rs_T, integer cs_T );


// --- Fused operations --------------------------------------------------------

FLA_Error FLA_Fused_Ahx_Ax_ops_var1( integer m_A,
                                     integer n_A,
                                     float* buff_A, integer rs_A, integer cs_A, 
                                     float* buff_x, integer inc_x, 
                                     float* buff_v, integer inc_v, 
                                     float* buff_w, integer inc_w );
FLA_Error FLA_Fused_Ahx_Ax_opd_var1( integer m_A,
                                     integer n_A,
                                     double* buff_A, integer rs_A, integer cs_A, 
                                     double* buff_x, integer inc_x, 
                                     double* buff_v, integer inc_v, 
                                     double* buff_w, integer inc_w );
FLA_Error FLA_Fused_Ahx_Ax_opc_var1( integer m_A,
                                     integer n_A,
                                     scomplex* buff_A, integer rs_A, integer cs_A, 
                                     scomplex* buff_x, integer inc_x, 
                                     scomplex* buff_v, integer inc_v, 
                                     scomplex* buff_w, integer inc_w );
FLA_Error FLA_Fused_Ahx_Ax_opz_var1( integer m_A,
                                     integer n_A,
                                     dcomplex* buff_A, integer rs_A, integer cs_A, 
                                     dcomplex* buff_x, integer inc_x, 
                                     dcomplex* buff_v, integer inc_v, 
                                     dcomplex* buff_w, integer inc_w );


FLA_Error FLA_Fused_Gerc2_Ahx_Ax_ops_var1( integer m_A,
                                           integer n_A,
                                           float* buff_alpha, 
                                           float* buff_u, integer inc_u, 
                                           float* buff_y, integer inc_y, 
                                           float* buff_z, integer inc_z, 
                                           float* buff_A, integer rs_A, integer cs_A, 
                                           float* buff_x, integer inc_x, 
                                           float* buff_v, integer inc_v, 
                                           float* buff_w, integer inc_w );
FLA_Error FLA_Fused_Gerc2_Ahx_Ax_opd_var1( integer m_A,
                                           integer n_A,
                                           double* buff_alpha, 
                                           double* buff_u, integer inc_u, 
                                           double* buff_y, integer inc_y, 
                                           double* buff_z, integer inc_z, 
                                           double* buff_A, integer rs_A, integer cs_A, 
                                           double* buff_x, integer inc_x, 
                                           double* buff_v, integer inc_v, 
                                           double* buff_w, integer inc_w );
FLA_Error FLA_Fused_Gerc2_Ahx_Ax_opc_var1( integer m_A,
                                           integer n_A,
                                           scomplex* buff_alpha, 
                                           scomplex* buff_u, integer inc_u, 
                                           scomplex* buff_y, integer inc_y, 
                                           scomplex* buff_z, integer inc_z, 
                                           scomplex* buff_A, integer rs_A, integer cs_A, 
                                           scomplex* buff_x, integer inc_x, 
                                           scomplex* buff_v, integer inc_v, 
                                           scomplex* buff_w, integer inc_w );
FLA_Error FLA_Fused_Gerc2_Ahx_Ax_opz_var1( integer m_A,
                                           integer n_A,
                                           dcomplex* buff_alpha, 
                                           dcomplex* buff_u, integer inc_u, 
                                           dcomplex* buff_y, integer inc_y, 
                                           dcomplex* buff_z, integer inc_z, 
                                           dcomplex* buff_A, integer rs_A, integer cs_A, 
                                           dcomplex* buff_x, integer inc_x, 
                                           dcomplex* buff_v, integer inc_v, 
                                           dcomplex* buff_w, integer inc_w );


FLA_Error FLA_Fused_Uhu_Yhu_Zhu_ops_var1( integer m_U,
                                          integer n_U,
                                          float* buff_delta,
                                          float* buff_U, integer rs_U, integer cs_U,
                                          float* buff_Y, integer rs_Y, integer cs_Y,
                                          float* buff_Z, integer rs_Z, integer cs_Z,
                                          float* buff_t, integer inc_t,
                                          float* buff_u, integer inc_u,
                                          float* buff_y, integer inc_y,
                                          float* buff_z, integer inc_z );
FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opd_var1( integer m_U,
                                          integer n_U,
                                          double* buff_delta,
                                          double* buff_U, integer rs_U, integer cs_U,
                                          double* buff_Y, integer rs_Y, integer cs_Y,
                                          double* buff_Z, integer rs_Z, integer cs_Z,
                                          double* buff_t, integer inc_t,
                                          double* buff_u, integer inc_u,
                                          double* buff_y, integer inc_y,
                                          double* buff_z, integer inc_z );
FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opc_var1( integer m_U,
                                          integer n_U,
                                          scomplex* buff_delta,
                                          scomplex* buff_U, integer rs_U, integer cs_U,
                                          scomplex* buff_Y, integer rs_Y, integer cs_Y,
                                          scomplex* buff_Z, integer rs_Z, integer cs_Z,
                                          scomplex* buff_t, integer inc_t,
                                          scomplex* buff_u, integer inc_u,
                                          scomplex* buff_y, integer inc_y,
                                          scomplex* buff_z, integer inc_z );
FLA_Error FLA_Fused_Uhu_Yhu_Zhu_opz_var1( integer m_U,
                                          integer n_U,
                                          dcomplex* buff_delta,
                                          dcomplex* buff_U, integer rs_U, integer cs_U,
                                          dcomplex* buff_Y, integer rs_Y, integer cs_Y,
                                          dcomplex* buff_Z, integer rs_Z, integer cs_Z,
                                          dcomplex* buff_t, integer inc_t,
                                          dcomplex* buff_u, integer inc_u,
                                          dcomplex* buff_y, integer inc_y,
                                          dcomplex* buff_z, integer inc_z );

