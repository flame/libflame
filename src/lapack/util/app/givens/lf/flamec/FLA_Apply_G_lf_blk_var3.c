/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Apply_G_lf_blk_var3( FLA_Obj G, FLA_Obj A, dim_t b_alg )
{
    FLA_Datatype datatype;
    int          k_G, m_A, n_A;
    int          rs_G, cs_G;
    int          rs_A, cs_A;

    datatype = FLA_Obj_datatype( A );

    k_G      = FLA_Obj_width( G );
    rs_G     = FLA_Obj_row_stride( G );
    cs_G     = FLA_Obj_col_stride( G );

    n_A      = FLA_Obj_length( A );
    m_A      = FLA_Obj_width( A );
    cs_A     = FLA_Obj_row_stride( A );
    rs_A     = FLA_Obj_col_stride( A );

    switch ( datatype )
    {
    case FLA_FLOAT:
    {
        scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
        float*    buff_A = ( float*    ) FLA_FLOAT_PTR( A );

        FLA_Apply_G_rf_bls_var3( k_G,
                                 m_A,
                                 n_A,
                                 buff_G, rs_G, cs_G,
                                 buff_A, rs_A, cs_A,
                                 b_alg );

        break;
    }

    case FLA_DOUBLE:
    {
        dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
        double*   buff_A = ( double*   ) FLA_DOUBLE_PTR( A );

        FLA_Apply_G_rf_bld_var3( k_G,
                                 m_A,
                                 n_A,
                                 buff_G, rs_G, cs_G,
                                 buff_A, rs_A, cs_A,
                                 b_alg );

        break;
    }

    case FLA_COMPLEX:
    {
        scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
        scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );

        FLA_Apply_G_rf_blc_var3( k_G,
                                 m_A,
                                 n_A,
                                 buff_G, rs_G, cs_G,
                                 buff_A, rs_A, cs_A,
                                 b_alg );

        break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
        dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
        dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );

        FLA_Apply_G_rf_blz_var3( k_G,
                                 m_A,
                                 n_A,
                                 buff_G, rs_G, cs_G,
                                 buff_A, rs_A, cs_A,
                                 b_alg );

        break;
    }
    }

    return FLA_SUCCESS;
}

