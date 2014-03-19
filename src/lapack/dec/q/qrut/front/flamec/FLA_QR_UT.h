/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLA_QR_UT_vars.h"

FLA_Error FLA_QR_UT( FLA_Obj A, FLA_Obj T );

FLA_Error FLA_QR_UT_internal( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl );
FLA_Error FLA_QR_UT_copy_internal( FLA_Obj A, FLA_Obj T, FLA_Obj U, fla_qrut_t* cntl );

FLA_Error FLA_QR_UT_create_T( FLA_Obj A, FLA_Obj* T );

FLA_Error FLA_QR_UT_recover_tau( FLA_Obj T, FLA_Obj tau );

FLA_Error FLA_QR_UT_solve( FLA_Obj A, FLA_Obj T, FLA_Obj B, FLA_Obj X );

FLA_Error FLASH_QR_UT( FLA_Obj A, FLA_Obj TW );
FLA_Error FLASH_QR_UT_create_hier_matrices( FLA_Obj A_flat, dim_t depth, dim_t* b_flash, FLA_Obj* A, FLA_Obj* TW );
FLA_Error FLASH_QR_UT_solve( FLA_Obj A, FLA_Obj T, FLA_Obj B, FLA_Obj X );


FLA_Error FLA_QR_UT_form_Q( FLA_Obj A, FLA_Obj T, FLA_Obj Q );
FLA_Error FLA_QR_UT_form_Q_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W );
FLA_Error FLA_QR_UT_form_Q_opt_var1( FLA_Obj A, FLA_Obj T );
FLA_Error FLA_QR_UT_form_Q_ops_var1( int       m_A,
                                     int       n_AT,
                                     float*    buff_A, int rs_A, int cs_A,
                                     float*    buff_T, int rs_T, int cs_T );
FLA_Error FLA_QR_UT_form_Q_opd_var1( int       m_A,
                                     int       n_AT,
                                     double*   buff_A, int rs_A, int cs_A,
                                     double*   buff_T, int rs_T, int cs_T );
FLA_Error FLA_QR_UT_form_Q_opc_var1( int       m_A,
                                     int       n_AT,
                                     scomplex* buff_A, int rs_A, int cs_A,
                                     scomplex* buff_T, int rs_T, int cs_T );
FLA_Error FLA_QR_UT_form_Q_opz_var1( int       m_A,
                                     int       n_AT,
                                     dcomplex* buff_A, int rs_A, int cs_A,
                                     dcomplex* buff_T, int rs_T, int cs_T );
