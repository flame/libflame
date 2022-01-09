/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

// Note that this operation is designed for tall rectangular matrix A.
// If m_A < n_A then, U and V should be swapped before entering this function.
FLA_Error FLA_Bsvd_v_opt_var1( dim_t n_iter_max, FLA_Obj d, FLA_Obj e, FLA_Obj G, FLA_Obj H, FLA_Obj U, FLA_Obj V, dim_t b_alg )
{
    FLA_Error    r_val = FLA_SUCCESS;
    FLA_Datatype datatype;
    integer          m_U, m_V, n_GH;
    integer          inc_d;
    integer          inc_e;
    integer          rs_G, cs_G;
    integer          rs_H, cs_H;
    integer          rs_U, cs_U;
    integer          rs_V, cs_V;

    datatype = FLA_Obj_datatype( U );

    m_U        = FLA_Obj_length( U );
    m_V        = FLA_Obj_length( V );
    n_GH       = FLA_Obj_width( G );

    inc_d      = FLA_Obj_vector_inc( d );
    inc_e      = FLA_Obj_vector_inc( e );

    rs_G       = FLA_Obj_row_stride( G );
    cs_G       = FLA_Obj_col_stride( G );

    rs_H       = FLA_Obj_row_stride( H );
    cs_H       = FLA_Obj_col_stride( H );

    rs_U       = FLA_Obj_row_stride( U );
    cs_U       = FLA_Obj_col_stride( U );

    rs_V       = FLA_Obj_row_stride( V );
    cs_V       = FLA_Obj_col_stride( V );


    switch ( datatype )
    {
    case FLA_FLOAT:
    {
        float*    buff_d = FLA_FLOAT_PTR( d );
        float*    buff_e = FLA_FLOAT_PTR( e );
        scomplex* buff_G = FLA_COMPLEX_PTR( G );
        scomplex* buff_H = FLA_COMPLEX_PTR( H );
        float*    buff_U = FLA_FLOAT_PTR( U );
        float*    buff_V = FLA_FLOAT_PTR( V );

        r_val = FLA_Bsvd_v_ops_var1( min( m_U, m_V ),
                                     m_U,
                                     m_V,
                                     n_GH,
                                     n_iter_max,
                                     buff_d, inc_d,
                                     buff_e, inc_e,
                                     buff_G, rs_G, cs_G,
                                     buff_H, rs_H, cs_H,
                                     buff_U, rs_U, cs_U,
                                     buff_V, rs_V, cs_V,
                                     b_alg );

        break;
    }

    case FLA_DOUBLE:
    {
        double*   buff_d = FLA_DOUBLE_PTR( d );
        double*   buff_e = FLA_DOUBLE_PTR( e );
        dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
        dcomplex* buff_H = FLA_DOUBLE_COMPLEX_PTR( H );
        double*   buff_U = FLA_DOUBLE_PTR( U );
        double*   buff_V = FLA_DOUBLE_PTR( V );

        r_val = FLA_Bsvd_v_opd_var1( min( m_U, m_V ),
                                     m_U,
                                     m_V,
                                     n_GH,
                                     n_iter_max,
                                     buff_d, inc_d,
                                     buff_e, inc_e,
                                     buff_G, rs_G, cs_G,
                                     buff_H, rs_H, cs_H,
                                     buff_U, rs_U, cs_U,
                                     buff_V, rs_V, cs_V,
                                     b_alg );

        break;
    }

    case FLA_COMPLEX:
    {
        float*    buff_d = FLA_FLOAT_PTR( d );
        float*    buff_e = FLA_FLOAT_PTR( e );
        scomplex* buff_G = FLA_COMPLEX_PTR( G );
        scomplex* buff_H = FLA_COMPLEX_PTR( H );
        scomplex* buff_U = FLA_COMPLEX_PTR( U );
        scomplex* buff_V = FLA_COMPLEX_PTR( V );

        r_val = FLA_Bsvd_v_opc_var1( min( m_U, m_V ),
                                     m_U,
                                     m_V,
                                     n_GH,
                                     n_iter_max,
                                     buff_d, inc_d,
                                     buff_e, inc_e,
                                     buff_G, rs_G, cs_G,
                                     buff_H, rs_H, cs_H,
                                     buff_U, rs_U, cs_U,
                                     buff_V, rs_V, cs_V,
                                     b_alg );

        break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
        double*   buff_d = FLA_DOUBLE_PTR( d );
        double*   buff_e = FLA_DOUBLE_PTR( e );
        dcomplex* buff_G = FLA_DOUBLE_COMPLEX_PTR( G );
        dcomplex* buff_H = FLA_DOUBLE_COMPLEX_PTR( H );
        dcomplex* buff_U = FLA_DOUBLE_COMPLEX_PTR( U );
        dcomplex* buff_V = FLA_DOUBLE_COMPLEX_PTR( V );

        r_val = FLA_Bsvd_v_opz_var1( min( m_U, m_V ),
                                     m_U,
                                     m_V,
                                     n_GH,
                                     n_iter_max,
                                     buff_d, inc_d,
                                     buff_e, inc_e,
                                     buff_G, rs_G, cs_G,
                                     buff_H, rs_H, cs_H,
                                     buff_U, rs_U, cs_U,
                                     buff_V, rs_V, cs_V,
                                     b_alg );

        break;
    }
    }

    return r_val;
}



FLA_Error FLA_Bsvd_v_ops_var1( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d,
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               scomplex* buff_H, integer rs_H, integer cs_H,
                               float*    buff_U, integer rs_U, integer cs_U,
                               float*    buff_V, integer rs_V, integer cs_V,
                               integer       b_alg )
{
    scomplex  one        = bl1_c1();
    float     rzero      = bl1_s0();

    integer       maxitr     = 6;

    float     eps;
    float     tolmul;
    float     tol;
    float     thresh;

    scomplex* G;
    scomplex* H;
    float*    d1;
    float*    e1;
    integer       r_val;
    integer       done;
    integer       m_GH_sweep_max;
    integer       ij_begin;
    integer       ijTL, ijBR;
    integer       m_A11;
    integer       n_iter_perf;
    integer       n_UV_apply;
    integer       total_deflations;
    integer       n_deflations;
    integer       n_iter_prev;
    integer       n_iter_perf_sweep_max;

    // Compute some convergence constants.
    eps    = FLA_Mach_params_ops( FLA_MACH_EPS );
    tolmul = max( 10.0F, min( 100.0F, powf( eps, -0.125F ) ) );
    FLA_Bsvd_compute_tol_thresh_ops( min_m_n,
                                     tolmul,
                                     maxitr,
                                     buff_d, inc_d,
                                     buff_e, inc_e,
                                     &tol,
                                     &thresh );

    // Initialize our completion flag.
    done = FALSE;

    // Initialize a counter that holds the maximum number of rows of G
    // and H that we would need to initialize for the next sweep.
    m_GH_sweep_max = min_m_n - 1;

    // Initialize a counter for the total number of iterations performed.
    n_iter_prev = 0;

    // Iterate until the matrix has completely deflated.
    for ( total_deflations = 0; done != TRUE; )
    {

        // Initialize G and H to contain only identity rotations.
        bl1_csetm( m_GH_sweep_max,
                   n_GH,
                   &one,
                   buff_G, rs_G, cs_G );
        bl1_csetm( m_GH_sweep_max,
                   n_GH,
                   &one,
                   buff_H, rs_H, cs_H );

        // Keep track of the maximum number of iterations performed in the
        // current sweep. This is used when applying the sweep's Givens
        // rotations.
        n_iter_perf_sweep_max = 0;

        // Perform a sweep: Move through the matrix and perform a bidiagonal
        // SVD on each non-zero submatrix that is encountered. During the
        // first time through, ijTL will be 0 and ijBR will be min_m_n - 1.
        for ( ij_begin = 0; ij_begin < min_m_n;  )
        {

            // Search for the first submatrix along the diagonal that is
            // bounded by zeroes (or endpoints of the matrix). If no
            // submatrix is found (ie: if the entire superdiagonal is zero
            // then FLA_FAILURE is returned. This function also inspects
            // superdiagonal elements for proximity to zero. If a given
            // element is close enough to zero, then it is deemed
            // converged and manually set to zero.
            r_val = FLA_Bsvd_find_submatrix_ops( min_m_n,
                                                 ij_begin,
                                                 buff_d, inc_d,
                                                 buff_e, inc_e,
                                                 &ijTL,
                                                 &ijBR );

            // Verify that a submatrix was found. If one was not found,
            // then we are done with the current sweep. Furthermore, if
            // a submatrix was not found AND we began our search at the
            // beginning of the matrix (ie: ij_begin == 0), then the
            // matrix has completely deflated and so we are done with
            // Francis step iteration.
            if ( r_val == FLA_FAILURE )
            {
                if ( ij_begin == 0 )
                {
                    done = TRUE;
                }

                // Break out of the current sweep so we can apply the last
                // remaining Givens rotations.
                break;
            }

            // If we got this far, then:
            //   (a) ijTL refers to the index of the first non-zero
            //       superdiagonal along the diagonal, and
            //   (b) ijBR refers to either:
            //       - the first zero element that occurs after ijTL, or
            //       - the the last diagonal element.
            // Note that ijTL and ijBR also correspond to the first and
            // last diagonal elements of the submatrix of interest. Thus,
            // we may compute the dimension of this submatrix as:
            m_A11 = ijBR - ijTL + 1;

            // Adjust ij_begin, which gets us ready for the next submatrix
            // search in the current sweep.
            ij_begin = ijBR + 1;

            // Index to the submatrices upon which we will operate.
            d1 = buff_d + ijTL * inc_d;
            e1 = buff_e + ijTL * inc_e;
            G  = buff_G + ijTL * rs_G;
            H  = buff_H + ijTL * rs_H;

            // Search for a batch of singular values, recursing on deflated
            // subproblems whenever a split occurs. Iteration continues as
            // long as
            //   (a) there is still matrix left to operate on, and
            //   (b) the number of iterations performed in this batch is
            //       less than n_G.
            // If/when either of the two above conditions fails to hold,
            // the function returns.
            n_deflations = FLA_Bsvd_iteracc_v_ops_var1( m_A11,
                                                        n_GH,
                                                        ijTL,
                                                        tol,
                                                        thresh,
                                                        d1, inc_d,
                                                        e1, inc_e,
                                                        G,  rs_G, cs_G,
                                                        H,  rs_H, cs_H,
                                                        &n_iter_perf );

            // Record the number of deflations that were observed.
            total_deflations += n_deflations;

            // Update the maximum number of iterations performed in the
            // current sweep.
            n_iter_perf_sweep_max = max( n_iter_perf_sweep_max, n_iter_perf );

            // Store the most recent value of ijBR in m_G_sweep_max.
            // When the sweep is done, this value will contain the minimum
            // number of rows of G and H we can apply and safely include all
            // non-identity rotations that were computed during the
            // singular value searches.
            m_GH_sweep_max = ijBR;

            // Make sure we haven't exceeded our maximum iteration count.
            if ( n_iter_prev >= n_iter_max * min_m_n )
            {
                FLA_Abort();
                //return FLA_FAILURE;
            }
        }

        // The sweep is complete. Now we must apply the Givens rotations
        // that were accumulated during the sweep.

        // Recall that the number of columns of U and V to which we apply
        // rotations is one more than the number of rotations.
        n_UV_apply = m_GH_sweep_max + 1;

        // Apply the Givens rotations. Note that we only apply k sets of
        // rotations, where k = n_iter_perf_sweep_max. Also note that we only
        // apply to n_UV_apply columns of U and V since this is the most we
        // need to touch given the most recent value stored to m_GH_sweep_max.
        //FLA_Apply_G_rf_bls_var5( n_iter_perf_sweep_max,
        FLA_Apply_G_rf_bls_var3( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_bld_var9( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_bld_var6( n_iter_perf_sweep_max,
                                 m_U,
                                 n_UV_apply,
                                 buff_G, rs_G, cs_G,
                                 buff_U, rs_U, cs_U,
                                 b_alg );
        //FLA_Apply_G_rf_blc_var5( n_iter_perf_sweep_max,
        FLA_Apply_G_rf_bls_var3( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_bld_var9( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_bld_var6( n_iter_perf_sweep_max,
                                 m_V,
                                 n_UV_apply,
                                 buff_H, rs_H, cs_H,
                                 buff_V, rs_V, cs_V,
                                 b_alg );

        // Increment the total number of iterations previously performed.
        n_iter_prev += n_iter_perf_sweep_max;

    }

    // Make all the singular values positive.
    {
        integer    i;
        float  minus_one = bl1_sm1();

        for ( i = 0; i < min_m_n; ++i )
        {
            if ( buff_d[ (i  )*inc_d ] < rzero )
            {
                buff_d[ (i  )*inc_d ] = -buff_d[ (i  )*inc_d ];

                // Scale the right singular vectors.
                bl1_sscalv( BLIS1_NO_CONJUGATE,
                            m_V,
                            &minus_one,
                            buff_V + (i  )*cs_V, rs_V );
            }
        }
    }

    return n_iter_prev;
}



FLA_Error FLA_Bsvd_v_opd_var1( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d,
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               dcomplex* buff_H, integer rs_H, integer cs_H,
                               double*   buff_U, integer rs_U, integer cs_U,
                               double*   buff_V, integer rs_V, integer cs_V,
                               integer       b_alg )
{
    dcomplex  one        = bl1_z1();
    double    rzero      = bl1_d0();

    integer       maxitr     = 6;

    double    eps;
    double    tolmul;
    double    tol;
    double    thresh;

    dcomplex* G;
    dcomplex* H;
    double*   d1;
    double*   e1;
    integer       r_val;
    integer       done;
    integer       m_GH_sweep_max;
    integer       ij_begin;
    integer       ijTL, ijBR;
    integer       m_A11;
    integer       n_iter_perf;
    integer       n_UV_apply;
    integer       total_deflations;
    integer       n_deflations;
    integer       n_iter_prev;
    integer       n_iter_perf_sweep_max;

    // Compute some convergence constants.
    eps    = FLA_Mach_params_opd( FLA_MACH_EPS );
    tolmul = max( 10.0, min( 100.0, pow( eps, -0.125 ) ) );
    FLA_Bsvd_compute_tol_thresh_opd( min_m_n,
                                     tolmul,
                                     maxitr,
                                     buff_d, inc_d,
                                     buff_e, inc_e,
                                     &tol,
                                     &thresh );
#ifdef PRINTF
    printf( "FLA_Bsvd_v_opd_var1: tolmul = %12.6e\n", tolmul );
    printf( "FLA_Bsvd_v_opd_var1: tol    = %12.6e\n", tol );
    printf( "FLA_Bsvd_v_opd_var1: thresh = %12.6e\n", thresh );
#endif

    // Initialize our completion flag.
    done = FALSE;

    // Initialize a counter that holds the maximum number of rows of G
    // and H that we would need to initialize for the next sweep.
    m_GH_sweep_max = min_m_n - 1;

    // Initialize a counter for the total number of iterations performed.
    n_iter_prev = 0;

    // Iterate until the matrix has completely deflated.
    for ( total_deflations = 0; done != TRUE; )
    {

        // Initialize G and H to contain only identity rotations.
        bl1_zsetm( m_GH_sweep_max,
                   n_GH,
                   &one,
                   buff_G, rs_G, cs_G );
        bl1_zsetm( m_GH_sweep_max,
                   n_GH,
                   &one,
                   buff_H, rs_H, cs_H );

        // Keep track of the maximum number of iterations performed in the
        // current sweep. This is used when applying the sweep's Givens
        // rotations.
        n_iter_perf_sweep_max = 0;

        // Perform a sweep: Move through the matrix and perform a bidiagonal
        // SVD on each non-zero submatrix that is encountered. During the
        // first time through, ijTL will be 0 and ijBR will be min_m_n - 1.
        for ( ij_begin = 0; ij_begin < min_m_n;  )
        {

#ifdef PRINTF
            if ( ij_begin == 0 )
                printf( "FLA_Bsvd_v_opd_var1: beginning new sweep (ij_begin = %d)\n", ij_begin );
#endif

            // Search for the first submatrix along the diagonal that is
            // bounded by zeroes (or endpoints of the matrix). If no
            // submatrix is found (ie: if the entire superdiagonal is zero
            // then FLA_FAILURE is returned. This function also inspects
            // superdiagonal elements for proximity to zero. If a given
            // element is close enough to zero, then it is deemed
            // converged and manually set to zero.
            r_val = FLA_Bsvd_find_submatrix_opd( min_m_n,
                                                 ij_begin,
                                                 buff_d, inc_d,
                                                 buff_e, inc_e,
                                                 &ijTL,
                                                 &ijBR );

            // Verify that a submatrix was found. If one was not found,
            // then we are done with the current sweep. Furthermore, if
            // a submatrix was not found AND we began our search at the
            // beginning of the matrix (ie: ij_begin == 0), then the
            // matrix has completely deflated and so we are done with
            // Francis step iteration.
            if ( r_val == FLA_FAILURE )
            {
                if ( ij_begin == 0 )
                {
#ifdef PRINTF
                    printf( "FLA_Bsvd_v_opd_var1: superdiagonal is completely zero.\n" );
                    printf( "FLA_Bsvd_v_opd_var1: Francis iteration is done!\n" );
#endif
                    done = TRUE;
                }

                // Break out of the current sweep so we can apply the last
                // remaining Givens rotations.
                break;
            }

            // If we got this far, then:
            //   (a) ijTL refers to the index of the first non-zero
            //       superdiagonal along the diagonal, and
            //   (b) ijBR refers to either:
            //       - the first zero element that occurs after ijTL, or
            //       - the the last diagonal element.
            // Note that ijTL and ijBR also correspond to the first and
            // last diagonal elements of the submatrix of interest. Thus,
            // we may compute the dimension of this submatrix as:
            m_A11 = ijBR - ijTL + 1;

#ifdef PRINTF
            printf( "FLA_Bsvd_v_opd_var1: ij_begin = %d\n", ij_begin );
            printf( "FLA_Bsvd_v_opd_var1: ijTL     = %d\n", ijTL );
            printf( "FLA_Bsvd_v_opd_var1: ijBR     = %d\n", ijBR );
            printf( "FLA_Bsvd_v_opd_var1: m_A11    = %d\n", m_A11 );
#endif

            // Adjust ij_begin, which gets us ready for the next submatrix
            // search in the current sweep.
            ij_begin = ijBR + 1;

            // Index to the submatrices upon which we will operate.
            d1 = buff_d + ijTL * inc_d;
            e1 = buff_e + ijTL * inc_e;
            G  = buff_G + ijTL * rs_G;
            H  = buff_H + ijTL * rs_H;

            // Search for a batch of singular values, recursing on deflated
            // subproblems whenever a split occurs. Iteration continues as
            // long as
            //   (a) there is still matrix left to operate on, and
            //   (b) the number of iterations performed in this batch is
            //       less than n_G.
            // If/when either of the two above conditions fails to hold,
            // the function returns.
            n_deflations = FLA_Bsvd_iteracc_v_opd_var1( m_A11,
                           n_GH,
                           ijTL,
                           tol,
                           thresh,
                           d1, inc_d,
                           e1, inc_e,
                           G,  rs_G, cs_G,
                           H,  rs_H, cs_H,
                           &n_iter_perf );

            // Record the number of deflations that were observed.
            total_deflations += n_deflations;

            // Update the maximum number of iterations performed in the
            // current sweep.
            n_iter_perf_sweep_max = max( n_iter_perf_sweep_max, n_iter_perf );

#ifdef PRINTF
            printf( "FLA_Bsvd_v_opd_var1: deflations observed       = %d\n", n_deflations );
            printf( "FLA_Bsvd_v_opd_var1: total deflations observed = %d\n", total_deflations );
            printf( "FLA_Bsvd_v_opd_var1: num iterations performed  = %d\n", n_iter_perf );
#endif

            // Store the most recent value of ijBR in m_G_sweep_max.
            // When the sweep is done, this value will contain the minimum
            // number of rows of G and H we can apply and safely include all
            // non-identity rotations that were computed during the
            // singular value searches.
            m_GH_sweep_max = ijBR;

            // Make sure we haven't exceeded our maximum iteration count.
            if ( n_iter_prev >= n_iter_max * min_m_n )
            {
#ifdef PRINTF
                printf( "FLA_Bsvd_v_opd_var1: reached maximum total number of iterations: %d\n", n_iter_prev );
#endif
                FLA_Abort();
                //return FLA_FAILURE;
            }
        }

        // The sweep is complete. Now we must apply the Givens rotations
        // that were accumulated during the sweep.

        // Recall that the number of columns of U and V to which we apply
        // rotations is one more than the number of rotations.
        n_UV_apply = m_GH_sweep_max + 1;

#ifdef PRINTF
        printf( "FLA_Bsvd_v_opd_var1: applying %d sets of Givens rotations\n", n_iter_perf_sweep_max );
        printf( "FLA_Bsvd_v_opd_var1: m_U = %d\n", m_U );
        printf( "FLA_Bsvd_v_opd_var1: m_V = %d\n", m_V );
        printf( "FLA_Bsvd_v_opd_var1: napp= %d\n", n_UV_apply );
#endif
        // Apply the Givens rotations. Note that we only apply k sets of
        // rotations, where k = n_iter_perf_sweep_max. Also note that we only
        // apply to n_UV_apply columns of U and V since this is the most we
        // need to touch given the most recent value stored to m_GH_sweep_max.
        //FLA_Apply_G_rf_bld_var5( n_iter_perf_sweep_max,
        FLA_Apply_G_rf_bld_var3( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_bld_var9( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_bld_var6( n_iter_perf_sweep_max,
                                 m_U,
                                 n_UV_apply,
                                 buff_G, rs_G, cs_G,
                                 buff_U, rs_U, cs_U,
                                 b_alg );
        //FLA_Apply_G_rf_bld_var5( n_iter_perf_sweep_max,
        FLA_Apply_G_rf_bld_var3( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_bld_var9( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_bld_var6( n_iter_perf_sweep_max,
                                 m_V,
                                 n_UV_apply,
                                 buff_H, rs_H, cs_H,
                                 buff_V, rs_V, cs_V,
                                 b_alg );

        // Increment the total number of iterations previously performed.
        n_iter_prev += n_iter_perf_sweep_max;

#ifdef PRINTF
        printf( "FLA_Bsvd_v_opd_var1: total number of iterations performed: %d\n", n_iter_prev );
#endif
    }

    // Make all the singular values positive.
    {
        integer    i;
        double minus_one = bl1_dm1();

        for ( i = 0; i < min_m_n; ++i )
        {
            if ( buff_d[ (i  )*inc_d ] < rzero )
            {
                buff_d[ (i  )*inc_d ] = -buff_d[ (i  )*inc_d ];

                // Scale the right singular vectors.
                bl1_dscalv( BLIS1_NO_CONJUGATE,
                            m_V,
                            &minus_one,
                            buff_V + (i  )*cs_V, rs_V );
            }
        }
    }

    return n_iter_prev;
}

FLA_Error FLA_Bsvd_v_opc_var1( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               float*    buff_d, integer inc_d,
                               float*    buff_e, integer inc_e,
                               scomplex* buff_G, integer rs_G, integer cs_G,
                               scomplex* buff_H, integer rs_H, integer cs_H,
                               scomplex* buff_U, integer rs_U, integer cs_U,
                               scomplex* buff_V, integer rs_V, integer cs_V,
                               integer       b_alg )
{
    scomplex  one        = bl1_c1();
    float     rzero      = bl1_s0();

    integer       maxitr     = 6;

    float     eps;
    float     tolmul;
    float     tol;
    float     thresh;

    scomplex* G;
    scomplex* H;
    float*    d1;
    float*    e1;
    integer       r_val;
    integer       done;
    integer       m_GH_sweep_max;
    integer       ij_begin;
    integer       ijTL, ijBR;
    integer       m_A11;
    integer       n_iter_perf;
    integer       n_UV_apply;
    integer       total_deflations;
    integer       n_deflations;
    integer       n_iter_prev;
    integer       n_iter_perf_sweep_max;

    // Compute some convergence constants.
    eps    = FLA_Mach_params_ops( FLA_MACH_EPS );
    tolmul = max( 10.0F, min( 100.0F, powf( eps, -0.125F ) ) );
    FLA_Bsvd_compute_tol_thresh_ops( min_m_n,
                                     tolmul,
                                     maxitr,
                                     buff_d, inc_d,
                                     buff_e, inc_e,
                                     &tol,
                                     &thresh );

    // Initialize our completion flag.
    done = FALSE;

    // Initialize a counter that holds the maximum number of rows of G
    // and H that we would need to initialize for the next sweep.
    m_GH_sweep_max = min_m_n - 1;

    // Initialize a counter for the total number of iterations performed.
    n_iter_prev = 0;

    // Iterate until the matrix has completely deflated.
    for ( total_deflations = 0; done != TRUE; )
    {

        // Initialize G and H to contain only identity rotations.
        bl1_csetm( m_GH_sweep_max,
                   n_GH,
                   &one,
                   buff_G, rs_G, cs_G );
        bl1_csetm( m_GH_sweep_max,
                   n_GH,
                   &one,
                   buff_H, rs_H, cs_H );

        // Keep track of the maximum number of iterations performed in the
        // current sweep. This is used when applying the sweep's Givens
        // rotations.
        n_iter_perf_sweep_max = 0;

        // Perform a sweep: Move through the matrix and perform a bidiagonal
        // SVD on each non-zero submatrix that is encountered. During the
        // first time through, ijTL will be 0 and ijBR will be min_m_n - 1.
        for ( ij_begin = 0; ij_begin < min_m_n;  )
        {
            // Search for the first submatrix along the diagonal that is
            // bounded by zeroes (or endpoints of the matrix). If no
            // submatrix is found (ie: if the entire superdiagonal is zero
            // then FLA_FAILURE is returned. This function also inspects
            // superdiagonal elements for proximity to zero. If a given
            // element is close enough to zero, then it is deemed
            // converged and manually set to zero.
            r_val = FLA_Bsvd_find_submatrix_ops( min_m_n,
                                                 ij_begin,
                                                 buff_d, inc_d,
                                                 buff_e, inc_e,
                                                 &ijTL,
                                                 &ijBR );

            // Verify that a submatrix was found. If one was not found,
            // then we are done with the current sweep. Furthermore, if
            // a submatrix was not found AND we began our search at the
            // beginning of the matrix (ie: ij_begin == 0), then the
            // matrix has completely deflated and so we are done with
            // Francis step iteration.
            if ( r_val == FLA_FAILURE )
            {
                if ( ij_begin == 0 )
                {
                    done = TRUE;
                }

                // Break out of the current sweep so we can apply the last
                // remaining Givens rotations.
                break;
            }

            // If we got this far, then:
            //   (a) ijTL refers to the index of the first non-zero
            //       superdiagonal along the diagonal, and
            //   (b) ijBR refers to either:
            //       - the first zero element that occurs after ijTL, or
            //       - the the last diagonal element.
            // Note that ijTL and ijBR also correspond to the first and
            // last diagonal elements of the submatrix of interest. Thus,
            // we may compute the dimension of this submatrix as:
            m_A11 = ijBR - ijTL + 1;

            // Adjust ij_begin, which gets us ready for the next submatrix
            // search in the current sweep.
            ij_begin = ijBR + 1;

            // Index to the submatrices upon which we will operate.
            d1 = buff_d + ijTL * inc_d;
            e1 = buff_e + ijTL * inc_e;
            G  = buff_G + ijTL * rs_G;
            H  = buff_H + ijTL * rs_H;

            // Search for a batch of singular values, recursing on deflated
            // subproblems whenever a split occurs. Iteration continues as
            // long as
            //   (a) there is still matrix left to operate on, and
            //   (b) the number of iterations performed in this batch is
            //       less than n_G.
            // If/when either of the two above conditions fails to hold,
            // the function returns.
            n_deflations = FLA_Bsvd_iteracc_v_ops_var1( m_A11,
                                                        n_GH,
                                                        ijTL,
                                                        tol,
                                                        thresh,
                                                        d1, inc_d,
                                                        e1, inc_e,
                                                        G,  rs_G, cs_G,
                                                        H,  rs_H, cs_H,
                                                        &n_iter_perf );

            // Record the number of deflations that were observed.
            total_deflations += n_deflations;

            // Update the maximum number of iterations performed in the
            // current sweep.
            n_iter_perf_sweep_max = max( n_iter_perf_sweep_max, n_iter_perf );

            // Store the most recent value of ijBR in m_G_sweep_max.
            // When the sweep is done, this value will contain the minimum
            // number of rows of G and H we can apply and safely include all
            // non-identity rotations that were computed during the
            // singular value searches.
            m_GH_sweep_max = ijBR;

            // Make sure we haven't exceeded our maximum iteration count.
            if ( n_iter_prev >= n_iter_max * min_m_n )
            {
                FLA_Abort();
                //return FLA_FAILURE;
            }
        }

        // The sweep is complete. Now we must apply the Givens rotations
        // that were accumulated during the sweep.

        // Recall that the number of columns of U and V to which we apply
        // rotations is one more than the number of rotations.
        n_UV_apply = m_GH_sweep_max + 1;

        // Apply the Givens rotations. Note that we only apply k sets of
        // rotations, where k = n_iter_perf_sweep_max. Also note that we only
        // apply to n_UV_apply columns of U and V since this is the most we
        // need to touch given the most recent value stored to m_GH_sweep_max.
        //FLA_Apply_G_rf_blc_var5( n_iter_perf_sweep_max,
        FLA_Apply_G_rf_blc_var3( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_blz_var9( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_blz_var6( n_iter_perf_sweep_max,
                                 m_U,
                                 n_UV_apply,
                                 buff_G, rs_G, cs_G,
                                 buff_U, rs_U, cs_U,
                                 b_alg );
        //FLA_Apply_G_rf_blc_var5( n_iter_perf_sweep_max,
        FLA_Apply_G_rf_blc_var3( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_blz_var9( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_blz_var6( n_iter_perf_sweep_max,
                                 m_V,
                                 n_UV_apply,
                                 buff_H, rs_H, cs_H,
                                 buff_V, rs_V, cs_V,
                                 b_alg );

        // Increment the total number of iterations previously performed.
        n_iter_prev += n_iter_perf_sweep_max;
    }

    // Make all the singular values positive.
    {
        integer    i;
        float  minus_one = bl1_sm1();

        for ( i = 0; i < min_m_n; ++i )
        {
            if ( buff_d[ (i  )*inc_d ] < rzero )
            {
                buff_d[ (i  )*inc_d ] = -buff_d[ (i  )*inc_d ];

                // Scale the right singular vectors.
                bl1_csscalv( BLIS1_NO_CONJUGATE,
                             m_V,
                             &minus_one,
                             buff_V + (i  )*cs_V, rs_V );
            }
        }
    }

    return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Bsvd_v_opz_var1( integer       min_m_n,
                               integer       m_U,
                               integer       m_V,
                               integer       n_GH,
                               integer       n_iter_max,
                               double*   buff_d, integer inc_d,
                               double*   buff_e, integer inc_e,
                               dcomplex* buff_G, integer rs_G, integer cs_G,
                               dcomplex* buff_H, integer rs_H, integer cs_H,
                               dcomplex* buff_U, integer rs_U, integer cs_U,
                               dcomplex* buff_V, integer rs_V, integer cs_V,
                               integer       b_alg )
{
    dcomplex  one        = bl1_z1();
    double    rzero      = bl1_d0();

    integer       maxitr     = 6;

    double    eps;
    double    tolmul;
    double    tol;
    double    thresh;

    dcomplex* G;
    dcomplex* H;
    double*   d1;
    double*   e1;
    integer       r_val;
    integer       done;
    integer       m_GH_sweep_max;
    integer       ij_begin;
    integer       ijTL, ijBR;
    integer       m_A11;
    integer       n_iter_perf;
    integer       n_UV_apply;
    integer       total_deflations;
    integer       n_deflations;
    integer       n_iter_prev;
    integer       n_iter_perf_sweep_max;

    // Compute some convergence constants.
    eps    = FLA_Mach_params_opd( FLA_MACH_EPS );
    tolmul = max( 10.0, min( 100.0, pow( eps, -0.125 ) ) );
    FLA_Bsvd_compute_tol_thresh_opd( min_m_n,
                                     tolmul,
                                     maxitr,
                                     buff_d, inc_d,
                                     buff_e, inc_e,
                                     &tol,
                                     &thresh );
#ifdef PRINTF
    printf( "FLA_Bsvd_v_opz_var1: tolmul = %12.6e\n", tolmul );
    printf( "FLA_Bsvd_v_opz_var1: tol    = %12.6e\n", tol );
    printf( "FLA_Bsvd_v_opz_var1: thresh = %12.6e\n", thresh );
#endif

    // Initialize our completion flag.
    done = FALSE;

    // Initialize a counter that holds the maximum number of rows of G
    // and H that we would need to initialize for the next sweep.
    m_GH_sweep_max = min_m_n - 1;

    // Initialize a counter for the total number of iterations performed.
    n_iter_prev = 0;

    // Iterate until the matrix has completely deflated.
    for ( total_deflations = 0; done != TRUE; )
    {

        // Initialize G and H to contain only identity rotations.
        bl1_zsetm( m_GH_sweep_max,
                   n_GH,
                   &one,
                   buff_G, rs_G, cs_G );
        bl1_zsetm( m_GH_sweep_max,
                   n_GH,
                   &one,
                   buff_H, rs_H, cs_H );

        // Keep track of the maximum number of iterations performed in the
        // current sweep. This is used when applying the sweep's Givens
        // rotations.
        n_iter_perf_sweep_max = 0;

        // Perform a sweep: Move through the matrix and perform a bidiagonal
        // SVD on each non-zero submatrix that is encountered. During the
        // first time through, ijTL will be 0 and ijBR will be min_m_n - 1.
        for ( ij_begin = 0; ij_begin < min_m_n;  )
        {

#ifdef PRINTF
            if ( ij_begin == 0 )
                printf( "FLA_Bsvd_v_opz_var1: beginning new sweep (ij_begin = %d)\n", ij_begin );
#endif

            // Search for the first submatrix along the diagonal that is
            // bounded by zeroes (or endpoints of the matrix). If no
            // submatrix is found (ie: if the entire superdiagonal is zero
            // then FLA_FAILURE is returned. This function also inspects
            // superdiagonal elements for proximity to zero. If a given
            // element is close enough to zero, then it is deemed
            // converged and manually set to zero.
            r_val = FLA_Bsvd_find_submatrix_opd( min_m_n,
                                                 ij_begin,
                                                 buff_d, inc_d,
                                                 buff_e, inc_e,
                                                 &ijTL,
                                                 &ijBR );

            // Verify that a submatrix was found. If one was not found,
            // then we are done with the current sweep. Furthermore, if
            // a submatrix was not found AND we began our search at the
            // beginning of the matrix (ie: ij_begin == 0), then the
            // matrix has completely deflated and so we are done with
            // Francis step iteration.
            if ( r_val == FLA_FAILURE )
            {
                if ( ij_begin == 0 )
                {
#ifdef PRINTF
                    printf( "FLA_Bsvd_v_opz_var1: superdiagonal is completely zero.\n" );
                    printf( "FLA_Bsvd_v_opz_var1: Francis iteration is done!\n" );
#endif
                    done = TRUE;
                }

                // Break out of the current sweep so we can apply the last
                // remaining Givens rotations.
                break;
            }

            // If we got this far, then:
            //   (a) ijTL refers to the index of the first non-zero
            //       superdiagonal along the diagonal, and
            //   (b) ijBR refers to either:
            //       - the first zero element that occurs after ijTL, or
            //       - the the last diagonal element.
            // Note that ijTL and ijBR also correspond to the first and
            // last diagonal elements of the submatrix of interest. Thus,
            // we may compute the dimension of this submatrix as:
            m_A11 = ijBR - ijTL + 1;

#ifdef PRINTF
            printf( "FLA_Bsvd_v_opz_var1: ij_begin = %d\n", ij_begin );
            printf( "FLA_Bsvd_v_opz_var1: ijTL     = %d\n", ijTL );
            printf( "FLA_Bsvd_v_opz_var1: ijBR     = %d\n", ijBR );
            printf( "FLA_Bsvd_v_opz_var1: m_A11    = %d\n", m_A11 );
#endif

            // Adjust ij_begin, which gets us ready for the next submatrix
            // search in the current sweep.
            ij_begin = ijBR + 1;

            // Index to the submatrices upon which we will operate.
            d1 = buff_d + ijTL * inc_d;
            e1 = buff_e + ijTL * inc_e;
            G  = buff_G + ijTL * rs_G;
            H  = buff_H + ijTL * rs_H;

            // Search for a batch of singular values, recursing on deflated
            // subproblems whenever a split occurs. Iteration continues as
            // long as
            //   (a) there is still matrix left to operate on, and
            //   (b) the number of iterations performed in this batch is
            //       less than n_G.
            // If/when either of the two above conditions fails to hold,
            // the function returns.
            n_deflations = FLA_Bsvd_iteracc_v_opd_var1( m_A11,
                           n_GH,
                           ijTL,
                           tol,
                           thresh,
                           d1, inc_d,
                           e1, inc_e,
                           G,  rs_G, cs_G,
                           H,  rs_H, cs_H,
                           &n_iter_perf );

            // Record the number of deflations that were observed.
            total_deflations += n_deflations;

            // Update the maximum number of iterations performed in the
            // current sweep.
            n_iter_perf_sweep_max = max( n_iter_perf_sweep_max, n_iter_perf );

#ifdef PRINTF
            printf( "FLA_Bsvd_v_opz_var1: deflations observed       = %d\n", n_deflations );
            printf( "FLA_Bsvd_v_opz_var1: total deflations observed = %d\n", total_deflations );
            printf( "FLA_Bsvd_v_opz_var1: num iterations performed  = %d\n", n_iter_perf );
#endif

            // Store the most recent value of ijBR in m_G_sweep_max.
            // When the sweep is done, this value will contain the minimum
            // number of rows of G and H we can apply and safely include all
            // non-identity rotations that were computed during the
            // singular value searches.
            m_GH_sweep_max = ijBR;

            // Make sure we haven't exceeded our maximum iteration count.
            if ( n_iter_prev >= n_iter_max * min_m_n )
            {
#ifdef PRINTF
                printf( "FLA_Bsvd_v_opz_var1: reached maximum total number of iterations: %d\n", n_iter_prev );
#endif
                FLA_Abort();
                //return FLA_FAILURE;
            }
        }

        // The sweep is complete. Now we must apply the Givens rotations
        // that were accumulated during the sweep.

        // Recall that the number of columns of U and V to which we apply
        // rotations is one more than the number of rotations.
        n_UV_apply = m_GH_sweep_max + 1;

#ifdef PRINTF
        printf( "FLA_Bsvd_v_opz_var1: applying %d sets of Givens rotations\n", n_iter_perf_sweep_max );
        printf( "FLA_Bsvd_v_opz_var1: m_U = %d\n", m_U );
        printf( "FLA_Bsvd_v_opz_var1: m_V = %d\n", m_V );
        printf( "FLA_Bsvd_v_opz_var1: napp= %d\n", n_UV_apply );
#endif
        // Apply the Givens rotations. Note that we only apply k sets of
        // rotations, where k = n_iter_perf_sweep_max. Also note that we only
        // apply to n_UV_apply columns of U and V since this is the most we
        // need to touch given the most recent value stored to m_GH_sweep_max.
        //FLA_Apply_G_rf_blz_var5( n_iter_perf_sweep_max,
        FLA_Apply_G_rf_blz_var3( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_blz_var9( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_blz_var6( n_iter_perf_sweep_max,
                                 m_U,
                                 n_UV_apply,
                                 buff_G, rs_G, cs_G,
                                 buff_U, rs_U, cs_U,
                                 b_alg );
        //FLA_Apply_G_rf_blz_var5( n_iter_perf_sweep_max,
        FLA_Apply_G_rf_blz_var3( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_blz_var9( n_iter_perf_sweep_max,
                                 //FLA_Apply_G_rf_blz_var6( n_iter_perf_sweep_max,
                                 m_V,
                                 n_UV_apply,
                                 buff_H, rs_H, cs_H,
                                 buff_V, rs_V, cs_V,
                                 b_alg );

        // Increment the total number of iterations previously performed.
        n_iter_prev += n_iter_perf_sweep_max;

#ifdef PRINTF
        printf( "FLA_Bsvd_v_opz_var1: total number of iterations performed: %d\n", n_iter_prev );
#endif
    }

    // Make all the singular values positive.
    {
        integer    i;
        double minus_one = bl1_dm1();

        for ( i = 0; i < min_m_n; ++i )
        {
            if ( buff_d[ (i  )*inc_d ] < rzero )
            {
                buff_d[ (i  )*inc_d ] = -buff_d[ (i  )*inc_d ];

                // Scale the right singular vectors.
                bl1_zdscalv( BLIS1_NO_CONJUGATE,
                             m_V,
                             &minus_one,
                             buff_V + (i  )*cs_V, rs_V );
            }
        }
    }

    return n_iter_prev;
}

