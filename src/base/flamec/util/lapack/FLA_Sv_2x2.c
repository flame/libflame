/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Sv_2x2( FLA_Obj alpha11, FLA_Obj alpha12, FLA_Obj alpha22,
                      FLA_Obj sigma1, FLA_Obj sigma2 )
/*
  Compute the singular value decomposition of a 2x2 triangular matrix A
  such that

    / alpha11 alpha12 \
    \    0    alpha22 /

  is equal to

    / gammaL -sigmaL \ / sigma1    0    \ / gammaR -sigmaR \'
    \ sigmaL  gammaL / \    0    sigma2 / \ sigmaR  gammaR /

  Upon completion, sigma1 and sigma2 are overwritten with the
  singular values of smaller and larger absolute values, respectively.
  gammaL, sigmaL, gammaR, and sigmaR are not computed.

  This routine is a nearly-verbatim translation of slas2() and dlas2()
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

        FLA_Sv_2x2_ops( buff_alpha11,
                        buff_alpha12,
                        buff_alpha22,
                        buff_sigma1,
                        buff_sigma2 );

        break;
    }

    case FLA_DOUBLE:
    {
        double* buff_alpha11 = FLA_DOUBLE_PTR( alpha11 );
        double* buff_alpha12 = FLA_DOUBLE_PTR( alpha12 );
        double* buff_alpha22 = FLA_DOUBLE_PTR( alpha22 );
        double* buff_sigma1  = FLA_DOUBLE_PTR( sigma1 );
        double* buff_sigma2  = FLA_DOUBLE_PTR( sigma2 );

        FLA_Sv_2x2_opd( buff_alpha11,
                        buff_alpha12,
                        buff_alpha22,
                        buff_sigma1,
                        buff_sigma2 );

        break;
    }
    }

    return FLA_SUCCESS;
}



FLA_Error FLA_Sv_2x2_ops( float*    alpha11,
                          float*    alpha12,
                          float*    alpha22,
                          float*    sigma1,
                          float*    sigma2 )
{
    float  zero = 0.0F;
    float  one  = 1.0F;
    float  two  = 2.0F;

    float  f, g, h;
    float  as, at, au, c, fa, fhmin, fhmax, ga, ha;
    float  ssmin, ssmax;
    float  temp, temp2;

    f = *alpha11;
    g = *alpha12;
    h = *alpha22;

    fa = fabsf( f );
    ga = fabsf( g );
    ha = fabsf( h );

    fhmin = min( fa, ha );
    fhmax = max( fa, ha );

    if ( fhmin == zero )
    {
        ssmin = zero;

        if ( fhmax == zero )
            ssmax = ga;
        else
        {
            temp = min( fhmax, ga ) / max( fhmax, ga );
            ssmax = max( fhmax, ga ) * sqrtf( one + temp * temp );
        }
    }
    else
    {
        if ( ga < fhmax )
        {
            as = one + fhmin / fhmax;
            at = ( fhmax - fhmin ) / fhmax;
            au = ( ga / fhmax ) * ( ga / fhmax );
            c  = two / ( sqrtf( as * as + au ) + sqrtf( at * at + au ) );
            ssmin = fhmin * c;
            ssmax = fhmax / c;
        }
        else
        {
            au = fhmax / ga;

            if ( au == zero )
            {
                ssmin = ( fhmin * fhmax ) / ga;
                ssmax = ga;
            }
            else
            {
                as = one + fhmin / fhmax;
                at = ( fhmax - fhmin ) / fhmax;
                temp  = as * au;
                temp2 = at * au;
                c  = one / ( sqrtf( one + temp * temp ) +
                             sqrtf( one + temp2 * temp2 ) );
                ssmin = ( fhmin * c ) * au;
                ssmin = ssmin + ssmin;
                ssmax = ga / ( c + c );
            }
        }
    }

    // Save the output values.

    *sigma1 = ssmin;
    *sigma2 = ssmax;

    return FLA_SUCCESS;
}



FLA_Error FLA_Sv_2x2_opd( double*   alpha11,
                          double*   alpha12,
                          double*   alpha22,
                          double*   sigma1,
                          double*   sigma2 )
{
    double zero = 0.0;
    double one  = 1.0;
    double two  = 2.0;

    double f, g, h;
    double as, at, au, c, fa, fhmin, fhmax, ga, ha;
    double ssmin, ssmax;
    double temp, temp2;

    f = *alpha11;
    g = *alpha12;
    h = *alpha22;

    fa = fabs( f );
    ga = fabs( g );
    ha = fabs( h );

    fhmin = min( fa, ha );
    fhmax = max( fa, ha );

    if ( fhmin == zero )
    {
        ssmin = zero;

        if ( fhmax == zero )
            ssmax = ga;
        else
        {
            temp = min( fhmax, ga ) / max( fhmax, ga );
            ssmax = max( fhmax, ga ) * sqrt( one + temp * temp );
        }
    }
    else
    {
        if ( ga < fhmax )
        {
            as = one + fhmin / fhmax;
            at = ( fhmax - fhmin ) / fhmax;
            au = ( ga / fhmax ) * ( ga / fhmax );
            c  = two / ( sqrt( as * as + au ) + sqrt( at * at + au ) );
            ssmin = fhmin * c;
            ssmax = fhmax / c;
        }
        else
        {
            au = fhmax / ga;

            if ( au == zero )
            {
                ssmin = ( fhmin * fhmax ) / ga;
                ssmax = ga;
            }
            else
            {
                as = one + fhmin / fhmax;
                at = ( fhmax - fhmin ) / fhmax;
                temp  = as * au;
                temp2 = at * au;
                c  = one / ( sqrt( one + temp * temp ) +
                             sqrt( one + temp2 * temp2 ) );
                ssmin = ( fhmin * c ) * au;
                ssmin = ssmin + ssmin;
                ssmax = ga / ( c + c );
            }
        }
    }

    // Save the output values.

    *sigma1 = ssmin;
    *sigma2 = ssmax;

    return FLA_SUCCESS;
}

