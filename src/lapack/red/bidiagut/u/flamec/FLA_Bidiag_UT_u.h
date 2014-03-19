/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

FLA_Error FLA_Bidiag_UT_u_unb_var1( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blk_var1( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_step_unb_var1( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_u_unb_var2( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blk_var2( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blf_var2( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_step_unb_var2( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_u_unb_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blk_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blf_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_step_unb_var3( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_u_unb_var4( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blk_var4( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blf_var4( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_step_unb_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_u_unb_var5( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_blk_var5( FLA_Obj A, FLA_Obj TU, FLA_Obj TV );
FLA_Error FLA_Bidiag_UT_u_step_unb_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj TU, FLA_Obj TV );

FLA_Error FLA_Bidiag_UT_u_opt_var1( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_opt_var1( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ops_var1( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opd_var1( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opc_var1( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opz_var1( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S );

FLA_Error FLA_Bidiag_UT_u_opt_var2( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_opt_var2( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ops_var2( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opd_var2( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opc_var2( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opz_var2( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S );

FLA_Error FLA_Bidiag_UT_u_opt_var3( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_opt_var3( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ops_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opd_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opc_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opz_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S );

FLA_Error FLA_Bidiag_UT_u_opt_var4( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_opt_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ops_var4( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_Y, int rs_Y, int cs_Y, 
                                         float* buff_Z, int rs_Z, int cs_Z, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opd_var4( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_Y, int rs_Y, int cs_Y, 
                                         double* buff_Z, int rs_Z, int cs_Z, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opc_var4( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_Y, int rs_Y, int cs_Y, 
                                         scomplex* buff_Z, int rs_Z, int cs_Z, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opz_var4( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_Y, int rs_Y, int cs_Y, 
                                         dcomplex* buff_Z, int rs_Z, int cs_Z, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S );

FLA_Error FLA_Bidiag_UT_u_opt_var5( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_opt_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ops_var5( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_Y, int rs_Y, int cs_Y, 
                                         float* buff_Z, int rs_Z, int cs_Z, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opd_var5( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_Y, int rs_Y, int cs_Y, 
                                         double* buff_Z, int rs_Z, int cs_Z, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opc_var5( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_Y, int rs_Y, int cs_Y, 
                                         scomplex* buff_Z, int rs_Z, int cs_Z, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_opz_var5( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_Y, int rs_Y, int cs_Y, 
                                         dcomplex* buff_Z, int rs_Z, int cs_Z, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S );


FLA_Error FLA_Bidiag_UT_u_ofu_var2( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofu_var2( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofs_var2( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofd_var2( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofc_var2( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofz_var2( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S );

FLA_Error FLA_Bidiag_UT_u_ofu_var3( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofu_var3( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofs_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofd_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofc_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofz_var3( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S );

FLA_Error FLA_Bidiag_UT_u_ofu_var4( FLA_Obj A, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofu_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj Z, FLA_Obj T, FLA_Obj S );
FLA_Error FLA_Bidiag_UT_u_step_ofs_var4( int m_A,
                                         int n_A,
                                         int m_TS,
                                         float* buff_A, int rs_A, int cs_A, 
                                         float* buff_Y, int rs_Y, int cs_Y, 
                                         float* buff_Z, int rs_Z, int cs_Z, 
                                         float* buff_T, int rs_T, int cs_T, 
                                         float* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofd_var4( int m_A,
                                         int n_A,
                                         int m_TS,
                                         double* buff_A, int rs_A, int cs_A, 
                                         double* buff_Y, int rs_Y, int cs_Y, 
                                         double* buff_Z, int rs_Z, int cs_Z, 
                                         double* buff_T, int rs_T, int cs_T, 
                                         double* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofc_var4( int m_A,
                                         int n_A,
                                         int m_TS,
                                         scomplex* buff_A, int rs_A, int cs_A, 
                                         scomplex* buff_Y, int rs_Y, int cs_Y, 
                                         scomplex* buff_Z, int rs_Z, int cs_Z, 
                                         scomplex* buff_T, int rs_T, int cs_T, 
                                         scomplex* buff_S, int rs_S, int cs_S );
FLA_Error FLA_Bidiag_UT_u_step_ofz_var4( int m_A,
                                         int n_A,
                                         int m_TS,
                                         dcomplex* buff_A, int rs_A, int cs_A, 
                                         dcomplex* buff_Y, int rs_Y, int cs_Y, 
                                         dcomplex* buff_Z, int rs_Z, int cs_Z, 
                                         dcomplex* buff_T, int rs_T, int cs_T, 
                                         dcomplex* buff_S, int rs_S, int cs_S );

// --- Fused operations ---

FLA_Error FLA_Fused_Gerc2_opt_var1( FLA_Obj alpha, FLA_Obj u, FLA_Obj y, FLA_Obj z, FLA_Obj v, FLA_Obj A );
FLA_Error FLA_Fused_Gerc2_ops_var1( int m_A,
                                    int n_A,
                                    float* buff_alpha, 
                                    float* buff_u, int inc_u, 
                                    float* buff_y, int inc_y, 
                                    float* buff_z, int inc_z, 
                                    float* buff_v, int inc_v, 
                                    float* buff_A, int rs_A, int cs_A ); 
FLA_Error FLA_Fused_Gerc2_opd_var1( int m_A,
                                    int n_A,
                                    double* buff_alpha, 
                                    double* buff_u, int inc_u, 
                                    double* buff_y, int inc_y, 
                                    double* buff_z, int inc_z, 
                                    double* buff_v, int inc_v, 
                                    double* buff_A, int rs_A, int cs_A ); 
FLA_Error FLA_Fused_Gerc2_opc_var1( int m_A,
                                    int n_A,
                                    scomplex* buff_alpha, 
                                    scomplex* buff_u, int inc_u, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_z, int inc_z, 
                                    scomplex* buff_v, int inc_v, 
                                    scomplex* buff_A, int rs_A, int cs_A ); 
FLA_Error FLA_Fused_Gerc2_opz_var1( int m_A,
                                    int n_A,
                                    dcomplex* buff_alpha, 
                                    dcomplex* buff_u, int inc_u, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_z, int inc_z, 
                                    dcomplex* buff_v, int inc_v, 
                                    dcomplex* buff_A, int rs_A, int cs_A ); 


FLA_Error FLA_Fused_Ahx_Axpy_Ax_opt_var1( FLA_Obj A, FLA_Obj u, FLA_Obj tau, FLA_Obj a, FLA_Obj beta, FLA_Obj y, FLA_Obj w );
FLA_Error FLA_Fused_Ahx_Axpy_Ax_ops_var1( int m_A,
                                          int n_A,
                                          float* buff_tau, 
                                          float* buff_beta, 
                                          float* buff_A, int rs_A, int cs_A, 
                                          float* buff_u, int inc_u, 
                                          float* buff_a, int inc_a, 
                                          float* buff_y, int inc_y, 
                                          float* buff_w, int inc_w );
FLA_Error FLA_Fused_Ahx_Axpy_Ax_opd_var1( int m_A,
                                          int n_A,
                                          double* buff_tau, 
                                          double* buff_beta, 
                                          double* buff_A, int rs_A, int cs_A, 
                                          double* buff_u, int inc_u, 
                                          double* buff_a, int inc_a, 
                                          double* buff_y, int inc_y, 
                                          double* buff_w, int inc_w );
FLA_Error FLA_Fused_Ahx_Axpy_Ax_opc_var1( int m_A,
                                          int n_A,
                                          scomplex* buff_tau, 
                                          scomplex* buff_beta, 
                                          scomplex* buff_A, int rs_A, int cs_A, 
                                          scomplex* buff_u, int inc_u, 
                                          scomplex* buff_a, int inc_a, 
                                          scomplex* buff_y, int inc_y, 
                                          scomplex* buff_w, int inc_w );
FLA_Error FLA_Fused_Ahx_Axpy_Ax_opz_var1( int m_A,
                                          int n_A,
                                          dcomplex* buff_tau, 
                                          dcomplex* buff_beta, 
                                          dcomplex* buff_A, int rs_A, int cs_A, 
                                          dcomplex* buff_u, int inc_u, 
                                          dcomplex* buff_a, int inc_a, 
                                          dcomplex* buff_y, int inc_y, 
                                          dcomplex* buff_w, int inc_w );

FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opt_var1( FLA_Obj alpha, FLA_Obj tau, FLA_Obj u, FLA_Obj y, FLA_Obj z, FLA_Obj v, FLA_Obj A, FLA_Obj up, FLA_Obj a, FLA_Obj w );
FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_ops_var1( int m_A,
                                                int n_A,
                                                float* buff_tau, 
                                                float* buff_alpha, 
                                                float* buff_u, int inc_u, 
                                                float* buff_y, int inc_y, 
                                                float* buff_z, int inc_z, 
                                                float* buff_v, int inc_v, 
                                                float* buff_A, int rs_A, int cs_A, 
                                                float* buff_up, int inc_up, 
                                                float* buff_a, int inc_a, 
                                                float* buff_w, int inc_w );
FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opd_var1( int m_A,
                                                int n_A,
                                                double* buff_tau, 
                                                double* buff_alpha, 
                                                double* buff_u, int inc_u, 
                                                double* buff_y, int inc_y, 
                                                double* buff_z, int inc_z, 
                                                double* buff_v, int inc_v, 
                                                double* buff_A, int rs_A, int cs_A, 
                                                double* buff_up, int inc_up, 
                                                double* buff_a, int inc_a, 
                                                double* buff_w, int inc_w );
FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opc_var1( int m_A,
                                                int n_A,
                                                scomplex* buff_tau, 
                                                scomplex* buff_alpha, 
                                                scomplex* buff_u, int inc_u, 
                                                scomplex* buff_y, int inc_y, 
                                                scomplex* buff_z, int inc_z, 
                                                scomplex* buff_v, int inc_v, 
                                                scomplex* buff_A, int rs_A, int cs_A, 
                                                scomplex* buff_up, int inc_up, 
                                                scomplex* buff_a, int inc_a, 
                                                scomplex* buff_w, int inc_w );
FLA_Error FLA_Fused_Gerc2_Ahx_Axpy_Ax_opz_var1( int m_A,
                                                int n_A,
                                                dcomplex* buff_tau, 
                                                dcomplex* buff_alpha, 
                                                dcomplex* buff_u, int inc_u, 
                                                dcomplex* buff_y, int inc_y, 
                                                dcomplex* buff_z, int inc_z, 
                                                dcomplex* buff_v, int inc_v, 
                                                dcomplex* buff_A, int rs_A, int cs_A, 
                                                dcomplex* buff_up, int inc_up, 
                                                dcomplex* buff_a, int inc_a, 
                                                dcomplex* buff_w, int inc_w );

FLA_Error FLA_Fused_UYx_ZVx_opt_var1( FLA_Obj delta, FLA_Obj a, FLA_Obj U, FLA_Obj Y, FLA_Obj Z, FLA_Obj V, FLA_Obj A, FLA_Obj temp, FLA_Obj t, FLA_Obj w, FLA_Obj al );
FLA_Error FLA_Fused_UYx_ZVx_ops_var1( int m_U,
                                      int n_U,
                                      int m_V,
                                      int n_V,
                                      float* buff_delta, 
                                      float* buff_U, int rs_U, int cs_U, 
                                      float* buff_Y, int rs_Y, int cs_Y, 
                                      float* buff_Z, int rs_Z, int cs_Z, 
                                      float* buff_V, int rs_V, int cs_V, 
                                      float* buff_A, int rs_A, int cs_A, 
                                      float* buff_temp, int inc_temp, 
                                      float* buff_t, int inc_t, 
                                      float* buff_a, int inc_a, 
                                      float* buff_w, int inc_w, 
                                      float* buff_al, int inc_al );
FLA_Error FLA_Fused_UYx_ZVx_opd_var1( int m_U,
                                      int n_U,
                                      int m_V,
                                      int n_V,
                                      double* buff_delta, 
                                      double* buff_U, int rs_U, int cs_U, 
                                      double* buff_Y, int rs_Y, int cs_Y, 
                                      double* buff_Z, int rs_Z, int cs_Z, 
                                      double* buff_V, int rs_V, int cs_V, 
                                      double* buff_A, int rs_A, int cs_A, 
                                      double* buff_temp, int inc_temp, 
                                      double* buff_t, int inc_t, 
                                      double* buff_a, int inc_a, 
                                      double* buff_w, int inc_w, 
                                      double* buff_al, int inc_al );
FLA_Error FLA_Fused_UYx_ZVx_opc_var1( int m_U,
                                      int n_U,
                                      int m_V,
                                      int n_V,
                                      scomplex* buff_delta, 
                                      scomplex* buff_U, int rs_U, int cs_U, 
                                      scomplex* buff_Y, int rs_Y, int cs_Y, 
                                      scomplex* buff_Z, int rs_Z, int cs_Z, 
                                      scomplex* buff_V, int rs_V, int cs_V, 
                                      scomplex* buff_A, int rs_A, int cs_A, 
                                      scomplex* buff_temp, int inc_temp, 
                                      scomplex* buff_t, int inc_t, 
                                      scomplex* buff_a, int inc_a, 
                                      scomplex* buff_w, int inc_w, 
                                      scomplex* buff_al, int inc_al );
FLA_Error FLA_Fused_UYx_ZVx_opz_var1( int m_U,
                                      int n_U,
                                      int m_V,
                                      int n_V,
                                      dcomplex* buff_delta, 
                                      dcomplex* buff_U, int rs_U, int cs_U, 
                                      dcomplex* buff_Y, int rs_Y, int cs_Y, 
                                      dcomplex* buff_Z, int rs_Z, int cs_Z, 
                                      dcomplex* buff_V, int rs_V, int cs_V, 
                                      dcomplex* buff_A, int rs_A, int cs_A, 
                                      dcomplex* buff_temp, int inc_temp, 
                                      dcomplex* buff_t, int inc_t, 
                                      dcomplex* buff_a, int inc_a, 
                                      dcomplex* buff_w, int inc_w, 
                                      dcomplex* buff_al, int inc_al );
