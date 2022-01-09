/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bsvd_sinval_v_opt_var1( FLA_Obj tol, FLA_Obj thresh, 
                                      FLA_Obj G, FLA_Obj H, 
                                      FLA_Obj d, FLA_Obj e, 
                                      FLA_Obj k )
{
    FLA_Datatype datatype;
    integer          m_A, n_GH;
    integer          rs_G, cs_G;
    integer          rs_H, cs_H;
    integer          inc_d;
    integer          inc_e;

    datatype = FLA_Obj_datatype( d );

    m_A      = FLA_Obj_vector_dim( d );
    n_GH     = FLA_Obj_width( G );

    rs_G     = FLA_Obj_row_stride( G );
    cs_G     = FLA_Obj_col_stride( G );

    rs_H     = FLA_Obj_row_stride( H );
    cs_H     = FLA_Obj_col_stride( H );

    inc_d    = FLA_Obj_vector_inc( d );
    inc_e    = FLA_Obj_vector_inc( e );


    switch ( datatype )
    {
    case FLA_FLOAT:
    {
        float*    buff_tol    = FLA_FLOAT_PTR( tol );
        float*    buff_thresh = FLA_FLOAT_PTR( thresh );
        scomplex* buff_G      = FLA_COMPLEX_PTR( G );
        scomplex* buff_H      = FLA_COMPLEX_PTR( H );
        float*    buff_d      = FLA_FLOAT_PTR( d );
        float*    buff_e      = FLA_FLOAT_PTR( e );
        integer*      buff_k      = FLA_INT_PTR( k );

        FLA_Bsvd_sinval_v_ops_var1( m_A,
                                    n_GH,
                                    9,
                                    *buff_tol,
                                    *buff_thresh,
                                    buff_G, rs_G, cs_G,
                                    buff_H, rs_H, cs_H,
                                    buff_d, inc_d,
                                    buff_e, inc_e,
                                    buff_k );

        break;
    }

    case FLA_DOUBLE:
    {
        double*   buff_tol    = FLA_DOUBLE_PTR( tol );
        double*   buff_thresh = FLA_DOUBLE_PTR( thresh );
        dcomplex* buff_G      = FLA_DOUBLE_COMPLEX_PTR( G );
        dcomplex* buff_H      = FLA_DOUBLE_COMPLEX_PTR( H );
        double*   buff_d      = FLA_DOUBLE_PTR( d );
        double*   buff_e      = FLA_DOUBLE_PTR( e );
        integer*      buff_k      = FLA_INT_PTR( k );

        FLA_Bsvd_sinval_v_opd_var1( m_A,
                                    n_GH,
                                    9,
                                    *buff_tol,
                                    *buff_thresh,
                                    buff_G, rs_G, cs_G,
                                    buff_H, rs_H, cs_H,
                                    buff_d, inc_d,
                                    buff_e, inc_e,
                                    buff_k );

        break;
    }
    }

    return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_sinval_v_ops_var1( integer       m_A,
                                      integer       n_GH,
                                      integer       n_iter_allowed,
                                      float     tol,
                                      float     thresh,
                                      scomplex* buff_G, integer rs_G, integer cs_G,
                                      scomplex* buff_H, integer rs_H, integer cs_H,
                                      float*    buff_d, integer inc_d,
                                      float*    buff_e, integer inc_e,
                                      integer*      n_iter )
{
    FLA_Error r_val;
    float     one = bl1_s1();
    //float*   d_first;
    //float*   d_last_m1;
    float*    e_last;
    //float*   d_last;
    float     smax;
    float     smin;
    float     sminl;
    float     shift;
    integer       k;

    // Initialize pointers to some diagonal and superdiagonal elements
    // that we will refer to later.
    e_last    = buff_e + (m_A-2)*inc_e;
    //d_last_m1 = buff_d + (m_A-2)*inc_d;
    //d_last    = buff_d + (m_A-1)*inc_d;
    //d_first   = buff_d + (0    )*inc_d;

    // Find the largest element of the diagonal or superdiagonal.
    // This is used later when checking the shift.
    FLA_Bsvd_find_max_min_ops( m_A,
                               buff_d, inc_d,
                               buff_e, inc_e,
                               &smax,
                               &smin );

    // Perform some iterations.
    for ( k = 0; k < n_iter_allowed; ++k )
    {
        scomplex* g1 = buff_G + (k  )*cs_G;
        scomplex* h1 = buff_H + (k  )*cs_H;

        /*------------------------------------------------------------*/

        // Before we perform any rotations, check for pre-existing deflation.
        r_val = FLA_Bsvd_find_converged_ops( m_A,
                                             tol,
                                             buff_d, inc_d,
                                             buff_e, inc_e,
                                             &sminl );

        // If r_val is positive, then deflation was found.
        if ( 0 <= r_val )
        {
            // Set the off-diagonal element to zero.
            buff_e[ (r_val)*inc_e ] = 0.0F;

            *n_iter = k;
            return r_val;
        }


        // Compute a shift with the last 2x2 matrix.
        FLA_Bsvd_compute_shift_ops( m_A,
                                    tol,
                                    sminl,
                                    smax,
                                    buff_d, inc_d,
                                    buff_e, inc_e,
                                    &shift );

        // Perform a Francis step.
        r_val = FLA_Bsvd_francis_v_ops_var1( m_A,
                                             shift,
                                             g1,     rs_G,
                                             h1,     rs_H,
                                             buff_d, inc_d,
                                             buff_e, inc_e );

        // Check for convergence using thresh.
        if ( MAC_Bsvd_sinval_is_converged_ops( thresh, one, *e_last ) )
        {
            *e_last = 0.0F;
            *n_iter = k + 1;
            return m_A - 1;
        }

        /*------------------------------------------------------------*/
    }

    *n_iter = n_iter_allowed;
    return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Bsvd_sinval_v_opd_var1( integer       m_A,
                                      integer       n_GH,
                                      integer       n_iter_allowed,
                                      double    tol,
                                      double    thresh,
                                      dcomplex* buff_G, integer rs_G, integer cs_G,
                                      dcomplex* buff_H, integer rs_H, integer cs_H,
                                      double*   buff_d, integer inc_d,
                                      double*   buff_e, integer inc_e,
                                      integer*      n_iter )
{
    FLA_Error r_val;
    double    one = bl1_d1();
    //double*   d_first;
    //double*   d_last_m1;
    double*   e_last;
    //double*   d_last;
    double    smax;
    double    smin;
    double    sminl;
    double    shift;
    integer       k;

    // Initialize pointers to some diagonal and superdiagonal elements
    // that we will refer to later.
    e_last    = buff_e + (m_A-2)*inc_e;
    //d_last_m1 = buff_d + (m_A-2)*inc_d;
    //d_last    = buff_d + (m_A-1)*inc_d;
    //d_first   = buff_d + (0    )*inc_d;

    // Find the largest element of the diagonal or superdiagonal.
    // This is used later when checking the shift.
    FLA_Bsvd_find_max_min_opd( m_A,
                               buff_d, inc_d,
                               buff_e, inc_e,
                               &smax,
                               &smin );

    // Perform some iterations.
    for ( k = 0; k < n_iter_allowed; ++k )
    {
        dcomplex* g1 = buff_G + (k  )*cs_G;
        dcomplex* h1 = buff_H + (k  )*cs_H;

        /*------------------------------------------------------------*/

        // Before we perform any rotations, check for pre-existing deflation.
        r_val = FLA_Bsvd_find_converged_opd( m_A,
                                             tol,
                                             buff_d, inc_d,
                                             buff_e, inc_e,
                                             &sminl );

        // If r_val is positive, then deflation was found.
        if ( 0 <= r_val )
        {
#ifdef PRINTF
            printf( "FLA_Bsvd_sinval_v_opt_var1: Deflation detected in col %d, sval %d\n", r_val, m_A - 1 );
            printf( "FLA_Bsvd_sinval_v_opt_var1: alpha11 alpha12 = %23.19e %23.19e\n", buff_d[r_val*inc_d], buff_e[r_val*inc_e] );
            printf( "FLA_Bsvd_sinval_v_opt_var1:         alpha22 =         %43.19e\n", buff_d[(r_val+1)*inc_d] );
#endif

            // Set the off-diagonal element to zero.
            buff_e[ (r_val)*inc_e ] = 0.0;

            *n_iter = k;
            return r_val;
        }


        // Compute a shift with the last 2x2 matrix.
        FLA_Bsvd_compute_shift_opd( m_A,
                                    tol,
                                    sminl,
                                    smax,
                                    buff_d, inc_d,
                                    buff_e, inc_e,
                                    &shift );

        // Perform a Francis step.
        r_val = FLA_Bsvd_francis_v_opd_var1( m_A,
                                             shift,
                                             g1,     rs_G,
                                             h1,     rs_H,
                                             buff_d, inc_d,
                                             buff_e, inc_e );

        // Check for convergence using thresh.
        if ( MAC_Bsvd_sinval_is_converged_opd( thresh, one, *e_last ) )
        {
            *e_last = 0.0;
            *n_iter = k + 1;
            return m_A - 1;
        }

        /*------------------------------------------------------------*/
    }

    *n_iter = n_iter_allowed;
    return FLA_FAILURE;
}

