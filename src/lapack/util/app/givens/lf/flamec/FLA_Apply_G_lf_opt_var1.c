
#include "FLAME.h"

FLA_Error FLA_Apply_G_lf_opt_var1( FLA_Obj G, FLA_Obj A )
/*
  Apply k sets of Givens rotations to a matrix A from the left,
  where each set takes the form:

    A := ( G(n-1,k) ... G(1,k) G(0,k) ) A

  where Gik is the ith Givens rotation formed from the kth set,
  stored in the (i,k) entries of of G:

    Gik  =  / gamma_ik  -sigma_ik \
            \ sigma_ik   gamma_ik /

  This variant iterates naively and applies rotations to two columns
  at a time.

  -FGVZ
*/
{
    FLA_Datatype datatype;
    int          k_G, m_A, n_A;
    int          rs_G, cs_G;
    int          rs_A, cs_A;

    datatype = FLA_Obj_datatype( A );

    k_G      = FLA_Obj_width( G );
    rs_G     = FLA_Obj_row_stride( G );
    cs_G     = FLA_Obj_col_stride( G );

    // Swap dimensions of A.
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

        FLA_Apply_G_rf_ops_var1( k_G,
                                 m_A,
                                 n_A,
                                 buff_G, rs_G, cs_G,
                                 buff_A, rs_A, cs_A );

        break;
    }

    case FLA_DOUBLE:
    {
        dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
        double*   buff_A = ( double*   ) FLA_DOUBLE_PTR( A );

        FLA_Apply_G_rf_opd_var1( k_G,
                                 m_A,
                                 n_A,
                                 buff_G, rs_G, cs_G,
                                 buff_A, rs_A, cs_A );

        break;
    }

    case FLA_COMPLEX:
    {
        scomplex* buff_G = ( scomplex* ) FLA_COMPLEX_PTR( G );
        scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );

        FLA_Apply_G_rf_opc_var1( k_G,
                                 m_A,
                                 n_A,
                                 buff_G, rs_G, cs_G,
                                 buff_A, rs_A, cs_A );

        break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
        dcomplex* buff_G = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( G );
        dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );

        FLA_Apply_G_rf_opz_var1( k_G,
                                 m_A,
                                 n_A,
                                 buff_G, rs_G, cs_G,
                                 buff_A, rs_A, cs_A );

        break;
    }
    }

    return FLA_SUCCESS;
}
