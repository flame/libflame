/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_Tridiag_UT_l.h"
//#include "FLA_Tridiag_UT_u.h"

FLA_Error FLA_Tridiag_UT( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T );

FLA_Error FLA_Tridiag_UT_internal( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T, fla_tridiagut_t* cntl );

FLA_Error FLA_Tridiag_UT_l( FLA_Obj A, FLA_Obj T, fla_tridiagut_t* cntl );
FLA_Error FLA_Tridiag_UT_u( FLA_Obj A, FLA_Obj T, fla_tridiagut_t* cntl );

FLA_Error FLA_Tridiag_UT_create_T( FLA_Obj A, FLA_Obj* T );
FLA_Error FLA_Tridiag_UT_recover_tau( FLA_Obj T, FLA_Obj t );

FLA_Error FLA_Tridiag_UT_scale_diagonals( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A );

FLA_Error FLA_Tridiag_UT_extract_diagonals( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d, FLA_Obj e );
FLA_Error FLA_Tridiag_UT_extract_real_diagonals( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d, FLA_Obj e );
//// FLA_Error FLA_Tridiag_UT_l_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e );
//// FLA_Error FLA_Tridiag_UT_u_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e );

FLA_Error FLA_Tridiag_UT_realify( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_l_realify_unb( FLA_Obj A, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_l_realify_opt( FLA_Obj A, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_u_realify_unb( FLA_Obj A, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_u_realify_opt( FLA_Obj A, FLA_Obj d );

FLA_Error FLA_Tridiag_UT_realify_subdiagonal( FLA_Obj b, FLA_Obj d );
FLA_Error FLA_Tridiag_UT_realify_subdiagonal_opt( FLA_Obj b, FLA_Obj d );

FLA_Error FLA_Tridiag_UT_shift_U( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Tridiag_UT_shift_U_l_ops( int       m_A,
                                        float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_u_ops( int       m_A,
                                        float*    buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_l_opd( int       m_A,
                                        double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_u_opd( int       m_A,
                                        double*   buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_l_opc( int       m_A,
                                        scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_u_opc( int       m_A,
                                        scomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_l_opz( int       m_A,
                                        dcomplex* buff_A, int rs_A, int cs_A );
FLA_Error FLA_Tridiag_UT_shift_U_u_opz( int       m_A,
                                        dcomplex* buff_A, int rs_A, int cs_A );

FLA_Error FLA_Tridiag_UT_form_Q( FLA_Uplo uplo, FLA_Obj A, FLA_Obj T, FLA_Obj Q );
FLA_Error FLA_Tridiag_UT_form_Q_l_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W );
FLA_Error FLA_Tridiag_UT_form_Q_u_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W );
FLA_Error FLA_Tridiag_UT_form_Q_l_opt_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_Tridiag_UT_form_Q_l_ops_var1( int       m_A,
                                            int       n_AT,
                                            float*    buff_A, int rs_A, int cs_A,
                                            float*    buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_form_Q_l_opd_var1( int       m_A,
                                            int       n_AT,
                                            double*   buff_A, int rs_A, int cs_A,
                                            double*   buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_form_Q_l_opc_var1( int       m_A,
                                            int       n_AT,
                                            scomplex* buff_A, int rs_A, int cs_A,
                                            scomplex* buff_T, int rs_T, int cs_T );
FLA_Error FLA_Tridiag_UT_form_Q_l_opz_var1( int       m_A,
                                            int       n_AT,
                                            dcomplex* buff_A, int rs_A, int cs_A,
                                            dcomplex* buff_T, int rs_T, int cs_T );
