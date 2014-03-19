/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Eig_gest_nu_blk_var1( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nu_blk_var2( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nu_blk_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nu_blk_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );
FLA_Error FLA_Eig_gest_nu_blk_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B, fla_eig_gest_t* cntl );

FLA_Error FLA_Eig_gest_nu_unb_var1( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_unb_var2( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_unb_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_unb_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_unb_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B );

FLA_Error FLA_Eig_gest_nu_opt_var1( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_ops_var1( int m_AB,
                                    float*    buff_A, int rs_A, int cs_A, 
                                    float*    buff_y, int inc_y,
                                    float*    buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opd_var1( int m_AB,
                                    double*   buff_A, int rs_A, int cs_A, 
                                    double*   buff_y, int inc_y,
                                    double*   buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opc_var1( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y,
                                    scomplex* buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opz_var1( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y,
                                    dcomplex* buff_B, int rs_B, int cs_B );

FLA_Error FLA_Eig_gest_nu_opt_var2( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_ops_var2( int m_AB,
                                    float*    buff_A, int rs_A, int cs_A, 
                                    float*    buff_y, int inc_y,
                                    float*    buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opd_var2( int m_AB,
                                    double*   buff_A, int rs_A, int cs_A, 
                                    double*   buff_y, int inc_y,
                                    double*   buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opc_var2( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y,
                                    scomplex* buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opz_var2( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y,
                                    dcomplex* buff_B, int rs_B, int cs_B );

FLA_Error FLA_Eig_gest_nu_opt_var3( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_ops_var3( int m_AB,
                                    float*    buff_A, int rs_A, int cs_A, 
                                    float*    buff_y, int inc_y, 
                                    float*    buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opd_var3( int m_AB,
                                    double*   buff_A, int rs_A, int cs_A, 
                                    double*   buff_y, int inc_y, 
                                    double*   buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opc_var3( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opz_var3( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B );

FLA_Error FLA_Eig_gest_nu_opt_var4( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_ops_var4( int m_AB,
                                    float*    buff_A, int rs_A, int cs_A, 
                                    float*    buff_y, int inc_y, 
                                    float*    buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opd_var4( int m_AB,
                                    double*   buff_A, int rs_A, int cs_A, 
                                    double*   buff_y, int inc_y, 
                                    double*   buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opc_var4( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opz_var4( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B );

FLA_Error FLA_Eig_gest_nu_opt_var5( FLA_Obj A, FLA_Obj Y, FLA_Obj B );
FLA_Error FLA_Eig_gest_nu_ops_var5( int m_AB,
                                    float*    buff_A, int rs_A, int cs_A, 
                                    float*    buff_y, int inc_y, 
                                    float*    buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opd_var5( int m_AB,
                                    double*   buff_A, int rs_A, int cs_A, 
                                    double*   buff_y, int inc_y, 
                                    double*   buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opc_var5( int m_AB,
                                    scomplex* buff_A, int rs_A, int cs_A, 
                                    scomplex* buff_y, int inc_y, 
                                    scomplex* buff_B, int rs_B, int cs_B );
FLA_Error FLA_Eig_gest_nu_opz_var5( int m_AB,
                                    dcomplex* buff_A, int rs_A, int cs_A, 
                                    dcomplex* buff_y, int inc_y, 
                                    dcomplex* buff_B, int rs_B, int cs_B );

