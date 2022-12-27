/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_QR_UT_form_Q( FLA_Obj A, FLA_Obj T, FLA_Obj Q )
{
    FLA_Error r_val = FLA_SUCCESS;
    FLA_Obj   QTL, QTR,
              QBL, QBR;
    FLA_Obj   W;
    dim_t     b;

    if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
        FLA_QR_UT_form_Q_check( A, T, Q );

    if ( FLA_Obj_is_overlapped( A, Q ) == FALSE )
    {
        // If A and Q are different objects, Q is explicitly formed with A.

        // Set Q identify
        FLA_Set_to_identity( Q );

        // Q = H_{0} H_{1} ... H_{k-1}
        FLA_Apply_Q_UT_create_workspace_side( FLA_LEFT, T, Q, &W );
        r_val = FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE,
                                FLA_FORWARD, FLA_COLUMNWISE,
                                A, T, W, Q );
        FLA_Obj_free( &W );

    }
    else
    {
        // If A and Q are the same objects, Q is formed in-place.
        // - even if A and Q has the same base, they may have different
        //   dimensions.
        // - width of T controls the loop in FLA_QR_UT_form_Q_blk_var1.

        // Zero out the upper triangle of Q.
        FLA_Setr( FLA_UPPER_TRIANGULAR, FLA_ZERO, Q );

        // Adjust T w.r.t A; W is a place holder.
        if ( FLA_Obj_width( T ) > FLA_Obj_width( A ) )
            FLA_Part_1x2( T,    &T,   &W,
                          FLA_Obj_width( A ),
                          FLA_LEFT );

        // Zero out the lower triangle of QBR
        if ( FLA_Obj_width( Q ) > FLA_Obj_width( T ) )
        {
            b = FLA_Obj_width( T );
            FLA_Part_2x2( Q, &QTL, &QTR,
                             &QBL, &QBR, b, b, FLA_TL );
            FLA_Setr( FLA_LOWER_TRIANGULAR, FLA_ZERO, QBR );
        }

        // Set the digaonal to one.
        FLA_Set_diag( FLA_ONE, Q );

        // Create workspace for applying the block Householder transforms.
        FLA_Apply_Q_UT_create_workspace_side( FLA_LEFT, T, Q, &W );

        // Overwrite Q, which currently contains Householder vectors in the
        // strictly lower triangle and identity in the upper triangle, with
        // the unitary matrix associated with those Householder transforms.
        r_val = FLA_QR_UT_form_Q_blk_var1( Q, T, W );

        // Free the temporary workspace.
        FLA_Obj_free( &W );
    }
    /*
      FLA_Apply_Q_UT_create_workspace( T, Q, &W );
      FLA_Set_to_identity( Q );
      FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
      A, T, W, Q );
      FLA_Obj_free( &W );
      FLA_Obj_show( "Q", Q, "%8.1e %8.1e ", "" );
    */

    return r_val;
}

FLA_Error FLA_QR_UT_form_Q_blk_var1( FLA_Obj A, FLA_Obj T, FLA_Obj W )
{
    FLA_Obj ATL,   ATR,      A00, A01, A02,
            ABL,   ABR,      A10, A11, A12,
            A20, A21, A22;

    FLA_Obj TL,    TR,       T0,  T1,  T2;

    FLA_Obj T1T,
            T2B;

    FLA_Obj WTL,  WTR,
            WBL,  WBR;

    FLA_Obj AB1,   AB2;

    dim_t   b, b_alg;
    dim_t   m_BR, n_BR;

    b_alg = FLA_Obj_length( T );


    // If A is wider than T, then we need to position ourseves carefully
    // within the matrix for the initial partitioning.
    if ( FLA_Obj_width( A ) > FLA_Obj_width( T ) )
    {
        m_BR = FLA_Obj_length( A ) - FLA_Obj_width( T );
        n_BR = FLA_Obj_width( A ) - FLA_Obj_width( T );
    }
    else
    {
        m_BR = FLA_Obj_length( A ) - FLA_Obj_width( A );
        n_BR = 0;
    }

    FLA_Part_2x2( A,    &ATL, &ATR,
                        &ABL, &ABR,     m_BR, n_BR, FLA_BR );

    FLA_Part_1x2( T,    &TL,  &TR,      0, FLA_RIGHT );

    while ( /* FLA_Obj_min_dim( ATL ) > 0 && */ FLA_Obj_width( TL ) > 0 )
    {
        b = fla_min( b_alg, FLA_Obj_min_dim( ATL ) );

        // Since T was filled from left to right, and since we need to access them
        // in reverse order, we need to handle the case where the last block is
        // smaller than the other b x b blocks.
        if ( FLA_Obj_width( TR ) == 0 && FLA_Obj_width( T ) % b_alg > 0 )
            b = FLA_Obj_width( T ) % b_alg;

        FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, &A01, /**/ &A02,
                                                    &A10, &A11, /**/ &A12,
                               /* ************* */   /* ******************** */
                               ABL, /**/ ABR,       &A20, &A21, /**/ &A22,
                               b, b, FLA_TL );

        FLA_Repart_1x2_to_1x3( TL,  /**/ TR,        &T0, &T1, /**/ &T2,
                               b, FLA_LEFT );

        /*------------------------------------------------------------*/

        FLA_Part_2x1( T1,    &T1T,
                             &T2B,     b, FLA_TOP );

        FLA_Part_2x2( W,     &WTL, &WTR,
                             &WBL, &WBR,     b, FLA_Obj_width( A12 ), FLA_TL );

        // Use an unblocked algorithm for the first (or only) block.
        if ( FLA_Obj_length( ABR ) == 0 )
        {
            FLA_QR_UT_form_Q_opt_var1( A11, T1T );
        }
        else
        {
            FLA_Merge_2x1( A11,
                           A21,   &AB1 );
            FLA_Merge_2x1( A12,
                           A22,   &AB2 );

            // Apply the block Householder transforms to A12 and A22.
            FLA_Apply_Q_UT( FLA_LEFT, FLA_NO_TRANSPOSE, FLA_FORWARD, FLA_COLUMNWISE,
                            AB1, T1T, WTL, AB2 );

            // Apply H to the current block panel consisting of A11 and A21.
            FLA_QR_UT_form_Q_opt_var1( AB1, T1T );
        }

        /*------------------------------------------------------------*/

        FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, /**/ A01, A02,
                                  /* ************** */  /* ****************** */
                                                         A10, /**/ A11, A12,
                                  &ABL, /**/ &ABR,       A20, /**/ A21, A22,
                                  FLA_BR );

        FLA_Cont_with_1x3_to_1x2( &TL,  /**/ &TR,        T0, /**/ T1, T2,
                                  FLA_RIGHT );
    }

    return FLA_SUCCESS;
}


FLA_Error FLA_QR_UT_form_Q_opt_var1( FLA_Obj A, FLA_Obj T )
{
    FLA_Datatype datatype;
    integer          m_A, n_A;
    integer          rs_A, cs_A;
    integer          rs_T, cs_T;

    datatype = FLA_Obj_datatype( A );

    m_A      = FLA_Obj_length( A );
    n_A      = FLA_Obj_width( A );
    rs_A     = FLA_Obj_row_stride( A );
    cs_A     = FLA_Obj_col_stride( A );

    rs_T     = FLA_Obj_row_stride( T );
    cs_T     = FLA_Obj_col_stride( T );

    switch ( datatype )
    {
    case FLA_FLOAT:
    {
        float*    buff_A = ( float* ) FLA_FLOAT_PTR( A );
        float*    buff_T = ( float* ) FLA_FLOAT_PTR( T );

        FLA_QR_UT_form_Q_ops_var1( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_T, rs_T, cs_T );

        break;
    }

    case FLA_DOUBLE:
    {
        double*   buff_A = ( double* ) FLA_DOUBLE_PTR( A );
        double*   buff_T = ( double* ) FLA_DOUBLE_PTR( T );

        FLA_QR_UT_form_Q_opd_var1( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_T, rs_T, cs_T );

        break;
    }

    case FLA_COMPLEX:
    {
        scomplex* buff_A = ( scomplex* ) FLA_COMPLEX_PTR( A );
        scomplex* buff_T = ( scomplex* ) FLA_COMPLEX_PTR( T );

        FLA_QR_UT_form_Q_opc_var1( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_T, rs_T, cs_T );

        break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
        dcomplex* buff_A = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );
        dcomplex* buff_T = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( T );

        FLA_QR_UT_form_Q_opz_var1( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_T, rs_T, cs_T );

        break;
    }
    }

    return FLA_SUCCESS;
}

FLA_Error FLA_QR_UT_form_Q_ops_var1( integer       m_A,
                                     integer       n_A,
                                     float*    buff_A, integer rs_A, integer cs_A,
                                     float*    buff_T, integer rs_T, integer cs_T )
{
    float    one     = bl1_d1();
    integer      min_m_n = fla_min( m_A, n_A );
    integer      i;

    for ( i = min_m_n - 1; i >= 0; --i )
    {
        //float*    a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
        float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
        float*    a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
        float*    A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

        float*    tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

        float     minus_inv_tau11;

        //integer       m_behind = i;
        integer       n_ahead  = n_A - i - 1;
        integer       m_ahead  = m_A - i - 1;

        FLA_Apply_H2_UT_l_ops_var1( m_ahead,
                                    n_ahead,
                                    tau11,
                                    a21,  rs_A,
                                    a12t, cs_A,
                                    A22,  rs_A, cs_A );

        minus_inv_tau11 = -one / *tau11;

        *alpha11 = one + minus_inv_tau11;

        bl1_sscalv( BLIS1_NO_CONJUGATE,
                    m_ahead,
                    &minus_inv_tau11,
                    a21, rs_A );

        // Not necessary if upper triangle of A is initialized to identity.
        //bl1_ssetv( m_behind,
        //           &zero,
        //           a01, rs_A );
    }

    return FLA_SUCCESS;
}

FLA_Error FLA_QR_UT_form_Q_opd_var1( integer       m_A,
                                     integer       n_A,
                                     double*   buff_A, integer rs_A, integer cs_A,
                                     double*   buff_T, integer rs_T, integer cs_T )
{
    double   one     = bl1_d1();
    integer      min_m_n = fla_min( m_A, n_A );
    integer      i;

    for ( i = min_m_n - 1; i >= 0; --i )
    {
        //double*   a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
        double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
        double*   a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
        double*   A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

        double*   tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

        double    minus_inv_tau11;

        //integer       m_behind = i;
        integer       n_ahead  = n_A - i - 1;
        integer       m_ahead  = m_A - i - 1;

        FLA_Apply_H2_UT_l_opd_var1( m_ahead,
                                    n_ahead,
                                    tau11,
                                    a21,  rs_A,
                                    a12t, cs_A,
                                    A22,  rs_A, cs_A );

        minus_inv_tau11 = -one / *tau11;

        *alpha11 = one + minus_inv_tau11;

        bl1_dscalv( BLIS1_NO_CONJUGATE,
                    m_ahead,
                    &minus_inv_tau11,
                    a21, rs_A );

        // Not necessary if upper triangle of A is initialized to identity.
        //bl1_dsetv( m_behind,
        //           &zero,
        //           a01, rs_A );
    }

    return FLA_SUCCESS;
}


FLA_Error FLA_QR_UT_form_Q_opc_var1( integer       m_A,
                                     integer       n_A,
                                     scomplex* buff_A, integer rs_A, integer cs_A,
                                     scomplex* buff_T, integer rs_T, integer cs_T )
{
    scomplex zero    = bl1_c0();
    scomplex one     = bl1_c1();
    integer      min_m_n = fla_min( m_A, n_A );
    integer      i;

    for ( i = min_m_n - 1; i >= 0; --i )
    {
        //scomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
        scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
        scomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
        scomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

        scomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

        scomplex  minus_inv_tau11;

        //integer       m_behind = i;
        integer       n_ahead  = n_A - i - 1;
        integer       m_ahead  = m_A - i - 1;

        FLA_Apply_H2_UT_l_opc_var1( m_ahead,
                                    n_ahead,
                                    tau11,
                                    a21,  rs_A,
                                    a12t, cs_A,
                                    A22,  rs_A, cs_A );

        minus_inv_tau11.real = -one.real / tau11->real;
        minus_inv_tau11.imag = zero.imag;

        alpha11->real = one.real + minus_inv_tau11.real;
        alpha11->imag = zero.imag;

        bl1_cscalv( BLIS1_NO_CONJUGATE,
                    m_ahead,
                    &minus_inv_tau11,
                    a21, rs_A );

        // Not necessary if upper triangle of A is initialized to identity.
        //bl1_csetv( m_behind,
        //           &zero,
        //           a01, rs_A );
    }

    return FLA_SUCCESS;
}

FLA_Error FLA_QR_UT_form_Q_opz_var1( integer       m_A,
                                     integer       n_A,
                                     dcomplex* buff_A, integer rs_A, integer cs_A,
                                     dcomplex* buff_T, integer rs_T, integer cs_T )
{
    dcomplex zero    = bl1_z0();
    dcomplex one     = bl1_z1();
    integer      min_m_n = fla_min( m_A, n_A );
    integer      i;

    for ( i = min_m_n - 1; i >= 0; --i )
    {
        //dcomplex* a01      = buff_A + (i  )*cs_A + (0  )*rs_A;
        dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        dcomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
        dcomplex* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;
        dcomplex* A22      = buff_A + (i+1)*cs_A + (i+1)*rs_A;

        dcomplex* tau11    = buff_T + (i  )*cs_T + (i  )*rs_T;

        dcomplex  minus_inv_tau11;

        //integer       m_behind = i;
        integer       n_ahead  = n_A - i - 1;
        integer       m_ahead  = m_A - i - 1;

        FLA_Apply_H2_UT_l_opz_var1( m_ahead,
                                    n_ahead,
                                    tau11,
                                    a21,  rs_A,
                                    a12t, cs_A,
                                    A22,  rs_A, cs_A );

        minus_inv_tau11.real = -one.real / tau11->real;
        minus_inv_tau11.imag = zero.imag;

        alpha11->real = one.real + minus_inv_tau11.real;
        alpha11->imag = zero.imag;

        bl1_zscalv( BLIS1_NO_CONJUGATE,
                    m_ahead,
                    &minus_inv_tau11,
                    a21, rs_A );

        // Not necessary if upper triangle of A is initialized to identity.
        //bl1_zsetv( m_behind,
        //           &zero,
        //           a01, rs_A );
    }

    return FLA_SUCCESS;
}

