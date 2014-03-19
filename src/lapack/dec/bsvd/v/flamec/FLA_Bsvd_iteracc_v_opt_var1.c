/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bsvd_iteracc_v_ops_var1( int       m_A,
                                       int       n_GH,
                                       int       ijTL,
                                       float     tol,
                                       float     thresh,
                                       float*    buff_d, int inc_d,
                                       float*    buff_e, int inc_e,
                                       scomplex* buff_G, int rs_G, int cs_G,
                                       scomplex* buff_H, int rs_H, int cs_H,
                                       int*      n_iter_perf )
{
    FLA_Error r_val;
    int       i, k;
    int       k_iter       = 0;
    int       n_deflations = 0;

    // Iterate from back to front until all that is left is a 2x2.
    for ( i = m_A - 1; i > 1; --i )
    {
        scomplex* G1     = buff_G + (k_iter)*cs_G;
        scomplex* H1     = buff_H + (k_iter)*cs_H;
        int       m_ATL  = i + 1;
        int       k_left = n_GH - k_iter;

        /*------------------------------------------------------------*/

        // Find a singular value of ATL submatrix.
        r_val = FLA_Bsvd_sinval_v_ops_var1( m_ATL,
                                            n_GH,
                                            k_left,
                                            tol,
                                            thresh,
                                            G1,     rs_G, cs_G,
                                            H1,     rs_H, cs_H,
                                            buff_d, inc_d,
                                            buff_e, inc_e,
                                            &k );

        // Update local counters according to the results of the singular
        // value search.
        k_iter       += k;
        n_deflations += 1;

        if ( r_val == FLA_FAILURE )
        {
            *n_iter_perf = k_iter;
            return n_deflations;
        }

        // If the most recent singular value search put us at our
        // limit for accumulated Givens rotation sets, return.
        if ( k_iter == n_GH )
        {
            *n_iter_perf = k_iter;
            return n_deflations;
        }

        // If r_val != i, then a split occurred somewhere within ATL.
        // Therefore, we must recurse with subproblems.
        if ( r_val != i )
        {
            int       m_TLr = r_val + 1;
            int       m_BRr = m_ATL - m_TLr;
            int       ijTLr = 0;
            int       ijBRr = m_TLr;
            int       n_GHr = n_GH - k_iter;
            float*    dTL   = buff_d + (0    )*inc_d;
            float*    eTL   = buff_e + (0    )*inc_e;
            scomplex* GT    = buff_G + (0    )*rs_G + (k_iter)*cs_G;
            scomplex* HT    = buff_H + (0    )*rs_H + (k_iter)*cs_H;
            float*    dBR   = buff_d + (ijBRr)*inc_d;
            float*    eBR   = buff_e + (ijBRr)*inc_e;
            scomplex* GB    = buff_G + (ijBRr)*rs_G + (k_iter)*cs_G;
            scomplex* HB    = buff_H + (ijBRr)*rs_H + (k_iter)*cs_H;

            int       n_deflationsTL;
            int       n_deflationsBR;
            int       n_iter_perfTL;
            int       n_iter_perfBR;

            n_deflationsTL = FLA_Bsvd_iteracc_v_ops_var1( m_TLr,
                                                          n_GHr,
                                                          ijTL + ijTLr,
                                                          tol,
                                                          thresh,
                                                          dTL, inc_d,
                                                          eTL, inc_e,
                                                          GT,  rs_G, cs_G,
                                                          HT,  rs_H, cs_H,
                                                          &n_iter_perfTL );
            n_deflationsBR = FLA_Bsvd_iteracc_v_ops_var1( m_BRr,
                                                          n_GHr,
                                                          ijTL + ijBRr,
                                                          tol,
                                                          thresh,
                                                          dBR, inc_d,
                                                          eBR, inc_e,
                                                          GB,  rs_G, cs_G,
                                                          HB,  rs_H, cs_H,
                                                          &n_iter_perfBR );

            *n_iter_perf = k_iter + max( n_iter_perfTL, n_iter_perfBR );

            return n_deflations + n_deflationsTL + n_deflationsBR;
        }

        /*------------------------------------------------------------*/
    }

    // Skip 1x1 matrices (and submatrices) entirely.
    if ( m_A > 1 )
    {
        scomplex* g1 = buff_G + (k_iter)*cs_G;
        scomplex* h1 = buff_H + (k_iter)*cs_H;

        float*    alpha11 = buff_d + (0  )*inc_d;
        float*    alpha12 = buff_e + (0  )*inc_e;
        float*    alpha22 = buff_d + (1  )*inc_d;

        float     smin;
        float     smax;

        float     gammaL;
        float     sigmaL;
        float     gammaR;
        float     sigmaR;

        // Find the singular value decomposition of the remaining (or only)
        // 2x2 submatrix.
        FLA_Svv_2x2_ops( alpha11,
                         alpha12,
                         alpha22,
                         &smin,
                         &smax,
                         &gammaL,
                         &sigmaL,
                         &gammaR,
                         &sigmaR );

        *alpha11 = smax;
        *alpha22 = smin;

        // Zero out the remaining diagonal.
        *alpha12 = 0.0F;

        // Store the rotations.
        g1[0].real = gammaL;
        g1[0].imag = sigmaL;
        h1[0].real = gammaR;
        h1[0].imag = sigmaR;

        // Update the local counters.
        k_iter       += 1;
        n_deflations += 1;

    }

    *n_iter_perf = k_iter;
    return n_deflations;
}

//#define PRINTF

FLA_Error FLA_Bsvd_iteracc_v_opd_var1( int       m_A,
                                       int       n_GH,
                                       int       ijTL,
                                       double    tol,
                                       double    thresh,
                                       double*   buff_d, int inc_d,
                                       double*   buff_e, int inc_e,
                                       dcomplex* buff_G, int rs_G, int cs_G,
                                       dcomplex* buff_H, int rs_H, int cs_H,
                                       int*      n_iter_perf )
{
    FLA_Error r_val;
    int       i, k;
    int       k_iter       = 0;
    int       n_deflations = 0;

    // Iterate from back to front until all that is left is a 2x2.
    for ( i = m_A - 1; i > 1; --i )
    {
        dcomplex* G1     = buff_G + (k_iter)*cs_G;
        dcomplex* H1     = buff_H + (k_iter)*cs_H;
        int       m_ATL  = i + 1;
        int       k_left = n_GH - k_iter;

        /*------------------------------------------------------------*/

        // Find a singular value of ATL submatrix.
        r_val = FLA_Bsvd_sinval_v_opd_var1( m_ATL,
                                            n_GH,
                                            k_left,
                                            tol,
                                            thresh,
                                            G1,     rs_G, cs_G,
                                            H1,     rs_H, cs_H,
                                            buff_d, inc_d,
                                            buff_e, inc_e,
                                            &k );

        // Update local counters according to the results of the singular
        // value search.
        k_iter       += k;
        n_deflations += 1;

        if ( r_val == FLA_FAILURE )
        {
#ifdef PRINTF
            printf( "FLA_Bsvd_iteracc_v_opd_var1: failed to converge (m_A11 = %d) after %2d iters k_total=%d/%d\n", i, k, k_iter, n_G );
#endif
            *n_iter_perf = k_iter;
            return n_deflations;
        }

#ifdef PRINTF
        if ( r_val == i )
            printf( "FLA_Bsvd_iteracc_v_opd_var1: found sv %22.15e in col %3d (n=%d) after %2d it  k_tot=%d/%d\n", buff_d[ r_val*inc_d ], ijTL+r_val, m_ATL, k, k_iter, n_GH );
        else
            printf( "FLA_Bsvd_iteracc_v_opd_var1: split occurred in col %3d. (n=%d) after %2d it  k_tot=%d/%d\n", r_val, m_ATL, k, k_iter, n_GH );
#endif

        // If the most recent singular value search put us at our
        // limit for accumulated Givens rotation sets, return.
        if ( k_iter == n_GH )
        {
            *n_iter_perf = k_iter;
            return n_deflations;
        }

        // If r_val != i, then a split occurred somewhere within ATL.
        // Therefore, we must recurse with subproblems.
        if ( r_val != i )
        {
            int       m_TLr = r_val + 1;
            int       m_BRr = m_ATL - m_TLr;
            int       ijTLr = 0;
            int       ijBRr = m_TLr;
            int       n_GHr = n_GH - k_iter;
            double*   dTL   = buff_d + (0    )*inc_d;
            double*   eTL   = buff_e + (0    )*inc_e;
            dcomplex* GT    = buff_G + (0    )*rs_G + (k_iter)*cs_G;
            dcomplex* HT    = buff_H + (0    )*rs_H + (k_iter)*cs_H;
            double*   dBR   = buff_d + (ijBRr)*inc_d;
            double*   eBR   = buff_e + (ijBRr)*inc_e;
            dcomplex* GB    = buff_G + (ijBRr)*rs_G + (k_iter)*cs_G;
            dcomplex* HB    = buff_H + (ijBRr)*rs_H + (k_iter)*cs_H;

            int       n_deflationsTL;
            int       n_deflationsBR;
            int       n_iter_perfTL;
            int       n_iter_perfBR;

#ifdef PRINTF
            printf( "FLA_Bsvd_iteracc_v_opd_var1: Deflation occurred in col %d\n", r_val );
            printf( "FLA_Bsvd_iteracc_v_opd_var1: alpha11 alpha12 = %22.15e %22.15e\n", buff_d[r_val*inc_d], buff_e[(r_val)*inc_e] );
            printf( "FLA_Bsvd_iteracc_v_opd_var1:         alpha22 =         %37.15e\n", buff_d[(r_val+1)*inc_d] );

            printf( "FLA_Bsvd_iteracc_v_opd_var1: recursing: ijTLr m_TLr: %d %d\n", ijTLr, m_TLr );
            printf( "FLA_Bsvd_iteracc_v_opd_var1:            GB(0,0) i,j: %d %d\n", ijTL + m_TLr+1, k_iter );
#endif
            n_deflationsTL = FLA_Bsvd_iteracc_v_opd_var1( m_TLr,
                             n_GHr,
                             ijTL + ijTLr,
                             tol,
                             thresh,
                             dTL, inc_d,
                             eTL, inc_e,
                             GT,  rs_G, cs_G,
                             HT,  rs_H, cs_H,
                             &n_iter_perfTL );
#ifdef PRINTF
            printf( "FLA_Bsvd_iteracc_v_opd_var1: returning: ijTLr m_TLr: %d %d\n", ijTLr, m_TLr );
            printf( "FLA_Bsvd_iteracc_v_opd_var1: recursing: ijBRr m_BRr: %d %d\n", ijBRr, m_BRr );
            printf( "FLA_Bsvd_iteracc_v_opd_var1:            GB(0,0) i,j: %d %d\n", ijTL + m_TLr+1, k_iter );
#endif
            n_deflationsBR = FLA_Bsvd_iteracc_v_opd_var1( m_BRr,
                             n_GHr,
                             ijTL + ijBRr,
                             tol,
                             thresh,
                             dBR, inc_d,
                             eBR, inc_e,
                             GB,  rs_G, cs_G,
                             HB,  rs_H, cs_H,
                             &n_iter_perfBR );
#ifdef PRINTF
            printf( "FLA_Bsvd_iteracc_v_opd_var1: returning: ijBRr m_BRr: %d %d\n", ijBRr, m_BRr );
#endif

            *n_iter_perf = k_iter + max( n_iter_perfTL, n_iter_perfBR );

            return n_deflations + n_deflationsTL + n_deflationsBR;
        }

        /*------------------------------------------------------------*/
    }

    // Skip 1x1 matrices (and submatrices) entirely.
    if ( m_A > 1 )
    {
        dcomplex* g1 = buff_G + (k_iter)*cs_G;
        dcomplex* h1 = buff_H + (k_iter)*cs_H;

        double*   alpha11 = buff_d + (0  )*inc_d;
        double*   alpha12 = buff_e + (0  )*inc_e;
        double*   alpha22 = buff_d + (1  )*inc_d;

        double    smin;
        double    smax;

        double    gammaL;
        double    sigmaL;
        double    gammaR;
        double    sigmaR;

        // Find the singular value decomposition of the remaining (or only)
        // 2x2 submatrix.
        FLA_Svv_2x2_opd( alpha11,
                         alpha12,
                         alpha22,
                         &smin,
                         &smax,
                         &gammaL,
                         &sigmaL,
                         &gammaR,
                         &sigmaR );

        *alpha11 = smax;
        *alpha22 = smin;

        // Zero out the remaining diagonal.
        *alpha12 = 0.0;

        // Store the rotations.
        g1[0].real = gammaL;
        g1[0].imag = sigmaL;
        h1[0].real = gammaR;
        h1[0].imag = sigmaR;

        // Update the local counters.
        k_iter       += 1;
        n_deflations += 1;

#ifdef PRINTF
        printf( "FLA_Bsvd_iteracc_v_opd_var1: Svv sval %22.15e in col %3d (n=%d) after %2d it  k_tot=%d/%d\n", buff_d[ 1*inc_d ], ijTL+1, 2, 1, k_iter, n_GH );
        printf( "FLA_Bsvd_iteracc_v_opd_var1: Svv sval %22.15e in col %3d (n=%d) after %2d it  k_tot=%d/%d\n", buff_d[ 0*inc_d ], ijTL+0, 2, 0, k_iter, n_GH );
#endif
    }

    *n_iter_perf = k_iter;
    return n_deflations;
}
