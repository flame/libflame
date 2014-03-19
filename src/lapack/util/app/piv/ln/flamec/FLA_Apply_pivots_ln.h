/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_pivots_ln_blk_var1( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );
FLA_Error FLA_Apply_pivots_ln_blk_var2( FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl );

FLA_Error FLA_Apply_pivots_ln_opt_var1( FLA_Obj p, FLA_Obj A );
FLA_Error FLA_Apply_pivots_ln_opi_var1( int n, 
                                        int*      a, int a_rs, int a_cs, 
                                        int k1, 
                                        int k2, 
                                        int* p, int incp );
FLA_Error FLA_Apply_pivots_ln_ops_var1( int n, 
                                        float*    a, int a_rs, int a_cs, 
                                        int k1, 
                                        int k2, 
                                        int* p, int incp );
FLA_Error FLA_Apply_pivots_ln_opd_var1( int n, 
                                        double*   a, int a_rs, int a_cs, 
                                        int k1, 
                                        int k2, 
                                        int* p, int incp );
FLA_Error FLA_Apply_pivots_ln_opc_var1( int n, 
                                        scomplex* a, int a_rs, int a_cs, 
                                        int k1, 
                                        int k2, 
                                        int* p, int incp );
FLA_Error FLA_Apply_pivots_ln_opz_var1( int n, 
                                        dcomplex* a, int a_rs, int a_cs, 
                                        int k1, 
                                        int k2, 
                                        int* p, int incp );
