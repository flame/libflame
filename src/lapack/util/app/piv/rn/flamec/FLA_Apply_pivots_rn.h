
#include "FLAME.h"

FLA_Error FLA_Apply_pivots_rn_opt_var1( FLA_Obj p, FLA_Obj A );
FLA_Error FLA_Apply_pivots_rn_ops_var1( int n, 
                                        float*    a, int a_rs, int a_cs, 
                                        int k1, 
                                        int k2, 
                                        int* p, int incp );
FLA_Error FLA_Apply_pivots_rn_opd_var1( int n, 
                                        double*   a, int a_rs, int a_cs, 
                                        int k1, 
                                        int k2, 
                                        int* p, int incp );
FLA_Error FLA_Apply_pivots_rn_opc_var1( int n, 
                                        scomplex* a, int a_rs, int a_cs, 
                                        int k1, 
                                        int k2, 
                                        int* p, int incp );
FLA_Error FLA_Apply_pivots_rn_opz_var1( int n, 
                                        dcomplex* a, int a_rs, int a_cs, 
                                        int k1, 
                                        int k2, 
                                        int* p, int incp );
