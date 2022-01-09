/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Svv_2x2( FLA_Obj alpha11, FLA_Obj alpha12, FLA_Obj alpha22,
                       FLA_Obj sigma1, FLA_Obj sigma2,
                       FLA_Obj gammaL, FLA_Obj sigmaL,
                       FLA_Obj gammaR, FLA_Obj sigmaR )
/*
  Compute the singular value decomposition of a 2x2 triangular matrix A
  such that

    / alpha11 alpha12 \
    \    0    alpha22 /

  is equal to

    / gammaL -sigmaL \ / sigma1    0    \ / gammaR -sigmaR \'
    \ sigmaL  gammaL / \    0    sigma2 / \ sigmaR  gammaR /

  Upon completion, sigma1 and sigma2 are overwritten with the
  singular values of smaller and larger absolute values, respectively,
  while gammaL, sigmaL, gammaR, and sigmaR determine the corresponding
  left and right singular vector elements.

  This routine is a nearly-verbatim translation of slasv2() and dlasv2()
  from the netlib distribution of LAPACK.

  -FGVZ
*/
{
    FLA_Datatype datatype;

    datatype = FLA_Obj_datatype( alpha11 );

    switch ( datatype )
    {
    case FLA_FLOAT:
    {
        float*  buff_alpha11 = FLA_FLOAT_PTR( alpha11 );
        float*  buff_alpha12 = FLA_FLOAT_PTR( alpha12 );
        float*  buff_alpha22 = FLA_FLOAT_PTR( alpha22 );
        float*  buff_sigma1  = FLA_FLOAT_PTR( sigma1 );
        float*  buff_sigma2  = FLA_FLOAT_PTR( sigma2 );
        float*  buff_gammaL  = FLA_FLOAT_PTR( gammaL );
        float*  buff_sigmaL  = FLA_FLOAT_PTR( sigmaL );
        float*  buff_gammaR  = FLA_FLOAT_PTR( gammaR );
        float*  buff_sigmaR  = FLA_FLOAT_PTR( sigmaR );

        FLA_Svv_2x2_ops( buff_alpha11,
                         buff_alpha12,
                         buff_alpha22,
                         buff_sigma1,
                         buff_sigma2,
                         buff_gammaL,
                         buff_sigmaL,
                         buff_gammaR,
                         buff_sigmaR );

        break;
    }

    case FLA_DOUBLE:
    {
        double* buff_alpha11 = FLA_DOUBLE_PTR( alpha11 );
        double* buff_alpha12 = FLA_DOUBLE_PTR( alpha12 );
        double* buff_alpha22 = FLA_DOUBLE_PTR( alpha22 );
        double* buff_sigma1  = FLA_DOUBLE_PTR( sigma1 );
        double* buff_sigma2  = FLA_DOUBLE_PTR( sigma2 );
        double* buff_gammaL  = FLA_DOUBLE_PTR( gammaL );
        double* buff_sigmaL  = FLA_DOUBLE_PTR( sigmaL );
        double* buff_gammaR  = FLA_DOUBLE_PTR( gammaR );
        double* buff_sigmaR  = FLA_DOUBLE_PTR( sigmaR );

        FLA_Svv_2x2_opd( buff_alpha11,
                         buff_alpha12,
                         buff_alpha22,
                         buff_sigma1,
                         buff_sigma2,
                         buff_gammaL,
                         buff_sigmaL,
                         buff_gammaR,
                         buff_sigmaR );

        break;
    }
    }

    return FLA_SUCCESS;
}



FLA_Error FLA_Svv_2x2_ops( float*    alpha11,
                           float*    alpha12,
                           float*    alpha22,
                           float*    sigma1,
                           float*    sigma2,
                           float*    gammaL,
                           float*    sigmaL,
                           float*    gammaR,
                           float*    sigmaR )
{
    float  zero = 0.0F;
    float  half = 0.5F;
    float  one  = 1.0F;
    float  two  = 2.0F;
    float  four = 4.0F;

    float  eps;

    float  f, g, h;
    float  clt, crt, slt, srt;
    float  a, d, fa, ft, ga, gt, ha, ht, l;
    float  m, mm, r, s, t, temp, tsign, tt;
    float  ssmin, ssmax;
    float  csl, snl;
    float  csr, snr;

    integer    gasmal, swap;
    integer    pmax;

    f = *alpha11;
    g = *alpha12;
    h = *alpha22;

    eps = FLA_Mach_params_ops( FLA_MACH_EPS );

    ft = f;
    fa = fabsf( f );
    ht = h;
    ha = fabsf( h );

    // pmax points to the maximum absolute element of matrix.
    //   pmax = 1 if f largest in absolute values.
    //   pmax = 2 if g largest in absolute values.
    //   pmax = 3 if h largest in absolute values.

    pmax = 1;

    swap = ( ha > fa );
    if ( swap )
    {
        pmax = 3;

        temp = ft;
        ft = ht;
        ht = temp;

        temp = fa;
        fa = ha;
        ha = temp;
    }

    gt = g;
    ga = fabsf( g );

    if ( ga == zero )
    {
        // Diagonal matrix case.

        ssmin = ha;
        ssmax = fa;
        clt   = one;
        slt   = zero;
        crt   = one;
        srt   = zero;
    }
    else
    {
        gasmal = TRUE;

        if ( ga > fa )
        {
            pmax = 2;

            if ( ( fa / ga ) < eps )
            {
                // Case of very large ga.

                gasmal = FALSE;

                ssmax  = ga;

                if ( ha > one ) ssmin = fa / ( ga / ha );
                else            ssmin = ( fa / ga ) * ha;

                clt = one;
                slt = ht / gt;
                crt = ft / gt;
                srt = one;
            }
        }

        if ( gasmal )
        {
            // Normal case.

            d = fa - ha;

            if ( d == fa ) l = one;
            else           l = d / fa;

            m = gt / ft;

            t = two - l;

            mm = m * m;
            tt = t * t;
            s = sqrtf( tt + mm );

            if ( l == zero ) r = fabsf( m );
            else             r = sqrtf( l * l + mm );

            a = half * ( s + r );

            ssmin = ha / a;
            ssmax = fa * a;

            if ( mm == zero )
            {
                // Here, m is tiny.

                if ( l == zero ) t = signof( two, ft ) * signof( one, gt );
                else             t = gt / signof( d, ft ) + m / t;
            }
            else
            {
                t = ( m / ( s + t ) + m / ( r + l ) ) * ( one + a );
            }

            l = sqrtf( t*t + four );
            crt = two / l;
            srt = t / l;
            clt = ( crt + srt * m ) / a;
            slt = ( ht / ft ) * srt / a;
        }
    }

    if ( swap )
    {
        csl = srt;
        snl = crt;
        csr = slt;
        snr = clt;
    }
    else
    {
        csl = clt;
        snl = slt;
        csr = crt;
        snr = srt;
    }


    // Correct the signs of ssmax and ssmin.

    if      ( pmax == 1 )
        tsign = signof( one, csr ) * signof( one, csl ) * signof( one, f );
    else if ( pmax == 2 )
        tsign = signof( one, snr ) * signof( one, csl ) * signof( one, g );
    else // if ( pmax == 3 )
        tsign = signof( one, snr ) * signof( one, snl ) * signof( one, h );

    ssmax = signof( ssmax, tsign );
    ssmin = signof( ssmin, tsign * signof( one, f ) * signof( one, h ) );

    // Save the output values.

    *sigma1 = ssmin;
    *sigma2 = ssmax;
    *gammaL = csl;
    *sigmaL = snl;
    *gammaR = csr;
    *sigmaR = snr;

    return FLA_SUCCESS;
}



FLA_Error FLA_Svv_2x2_opd( double*   alpha11,
                           double*   alpha12,
                           double*   alpha22,
                           double*   sigma1,
                           double*   sigma2,
                           double*   gammaL,
                           double*   sigmaL,
                           double*   gammaR,
                           double*   sigmaR )
{
    double zero = 0.0;
    double half = 0.5;
    double one  = 1.0;
    double two  = 2.0;
    double four = 4.0;

    double eps;

    double f, g, h;
    double clt, crt, slt, srt;
    double a, d, fa, ft, ga, gt, ha, ht, l;
    double m, mm, r, s, t, temp, tsign, tt;
    double ssmin, ssmax;
    double csl, snl;
    double csr, snr;

    integer    gasmal, swap;
    integer    pmax;

    f = *alpha11;
    g = *alpha12;
    h = *alpha22;

    eps = FLA_Mach_params_opd( FLA_MACH_EPS );

    ft = f;
    fa = fabs( f );
    ht = h;
    ha = fabs( h );

    // pmax points to the maximum absolute element of matrix.
    //   pmax = 1 if f largest in absolute values.
    //   pmax = 2 if g largest in absolute values.
    //   pmax = 3 if h largest in absolute values.

    pmax = 1;

    swap = ( ha > fa );
    if ( swap )
    {
        pmax = 3;

        temp = ft;
        ft = ht;
        ht = temp;

        temp = fa;
        fa = ha;
        ha = temp;
    }

    gt = g;
    ga = fabs( g );

    if ( ga == zero )
    {
        // Diagonal matrix case.

        ssmin = ha;
        ssmax = fa;
        clt   = one;
        slt   = zero;
        crt   = one;
        srt   = zero;
    }
    else
    {
        gasmal = TRUE;

        if ( ga > fa )
        {
            pmax = 2;

            if ( ( fa / ga ) < eps )
            {
                // Case of very large ga.

                gasmal = FALSE;

                ssmax  = ga;

                if ( ha > one ) ssmin = fa / ( ga / ha );
                else            ssmin = ( fa / ga ) * ha;

                clt = one;
                slt = ht / gt;
                crt = ft / gt;
                srt = one;
            }
        }

        if ( gasmal )
        {
            // Normal case.

            d = fa - ha;

            if ( d == fa ) l = one;
            else           l = d / fa;

            m = gt / ft;

            t = two - l;

            mm = m * m;
            tt = t * t;
            s = sqrt( tt + mm );

            if ( l == zero ) r = fabs( m );
            else             r = sqrt( l * l + mm );

            a = half * ( s + r );

            ssmin = ha / a;
            ssmax = fa * a;

            if ( mm == zero )
            {
                // Here, m is tiny.

                if ( l == zero ) t = signof( two, ft ) * signof( one, gt );
                else             t = gt / signof( d, ft ) + m / t;
            }
            else
            {
                t = ( m / ( s + t ) + m / ( r + l ) ) * ( one + a );
            }

            l = sqrt( t*t + four );
            crt = two / l;
            srt = t / l;
            clt = ( crt + srt * m ) / a;
            slt = ( ht / ft ) * srt / a;
        }
    }

    if ( swap )
    {
        csl = srt;
        snl = crt;
        csr = slt;
        snr = clt;
    }
    else
    {
        csl = clt;
        snl = slt;
        csr = crt;
        snr = srt;
    }


    // Correct the signs of ssmax and ssmin.

    if      ( pmax == 1 )
        tsign = signof( one, csr ) * signof( one, csl ) * signof( one, f );
    else if ( pmax == 2 )
        tsign = signof( one, snr ) * signof( one, csl ) * signof( one, g );
    else // if ( pmax == 3 )
        tsign = signof( one, snr ) * signof( one, snl ) * signof( one, h );

    ssmax = signof( ssmax, tsign );
    ssmin = signof( ssmin, tsign * signof( one, f ) * signof( one, h ) );

    // Save the output values.

    *sigma1 = ssmin;
    *sigma2 = ssmax;
    *gammaL = csl;
    *sigmaL = snl;
    *gammaR = csr;
    *sigmaR = snr;

    return FLA_SUCCESS;
}

