/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_pivots_rn_opt_var1( FLA_Obj p, FLA_Obj A );
FLA_Error FLA_Apply_pivots_rn_ops_var1( integer n, 
                                        float*    a, integer a_rs, integer a_cs, 
                                        integer k1, 
                                        integer k2, 
                                        integer* p, integer incp );
FLA_Error FLA_Apply_pivots_rn_opd_var1( integer n, 
                                        double*   a, integer a_rs, integer a_cs, 
                                        integer k1, 
                                        integer k2, 
                                        integer* p, integer incp );
FLA_Error FLA_Apply_pivots_rn_opc_var1( integer n, 
                                        scomplex* a, integer a_rs, integer a_cs, 
                                        integer k1, 
                                        integer k2, 
                                        integer* p, integer incp );
FLA_Error FLA_Apply_pivots_rn_opz_var1( integer n, 
                                        dcomplex* a, integer a_rs, integer a_cs, 
                                        integer k1, 
                                        integer k2, 
                                        integer* p, integer incp );
