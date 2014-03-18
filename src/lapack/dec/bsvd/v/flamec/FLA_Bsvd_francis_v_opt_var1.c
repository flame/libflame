
#include "FLAME.h"

FLA_Error FLA_Bsvd_francis_v_opt_var1( FLA_Obj shift, FLA_Obj g, FLA_Obj h, FLA_Obj d, FLA_Obj e )
{
    FLA_Datatype datatype;
    int          m_A;
    int          inc_g;
    int          inc_h;
    int          inc_d;
    int          inc_e;

    datatype = FLA_Obj_datatype( d );

    m_A      = FLA_Obj_vector_dim( d );

    inc_g    = FLA_Obj_vector_inc( g );
    inc_h    = FLA_Obj_vector_inc( h );
    inc_d    = FLA_Obj_vector_inc( d );
    inc_e    = FLA_Obj_vector_inc( e );


    switch ( datatype )
    {
    case FLA_FLOAT:
    {
        float*    buff_shift  = FLA_FLOAT_PTR( shift );
        scomplex* buff_g      = FLA_COMPLEX_PTR( g );
        scomplex* buff_h      = FLA_COMPLEX_PTR( h );
        float*    buff_d      = FLA_FLOAT_PTR( d );
        float*    buff_e      = FLA_FLOAT_PTR( e );

        FLA_Bsvd_francis_v_ops_var1( m_A,
                                     *buff_shift,
                                     buff_g, inc_g,
                                     buff_h, inc_h,
                                     buff_d, inc_d,
                                     buff_e, inc_e );

        break;
    }

    case FLA_DOUBLE:
    {
        double*   buff_shift  = FLA_DOUBLE_PTR( shift );
        dcomplex* buff_g      = FLA_DOUBLE_COMPLEX_PTR( g );
        dcomplex* buff_h      = FLA_DOUBLE_COMPLEX_PTR( h );
        double*   buff_d      = FLA_DOUBLE_PTR( d );
        double*   buff_e      = FLA_DOUBLE_PTR( e );

        FLA_Bsvd_francis_v_opd_var1( m_A,
                                     *buff_shift,
                                     buff_g, inc_g,
                                     buff_h, inc_h,
                                     buff_d, inc_d,
                                     buff_e, inc_e );

        break;
    }
    }

    return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_francis_v_ops_var1( int       m_A,
                                       float     shift,
                                       scomplex* buff_g, int inc_g,
                                       scomplex* buff_h, int inc_h,
                                       float*    buff_d, int inc_d,
                                       float*    buff_e, int inc_e )
{
    float     one   = bl1_s1();
    float     bulge = 0.0F;
    int       i;

    // If the shift is zero, perform a simplified Francis step.
    if ( shift == 0.0F )
    {
        float cs, oldcs;
        float sn, oldsn = 0.0F;
        float r, temp;
        float a11temp, a22temp;

        float* d_last = buff_d + (m_A-1)*inc_d;
        float* e_last = buff_e + (m_A-2)*inc_e;

        cs    = one;
        oldcs = one;

        for ( i = 0; i < m_A - 1; ++i )
        {
            float*   alpha01  = buff_e + (i-1)*inc_e;
            float*   alpha11  = buff_d + (i  )*inc_d;
            float*   alpha12  = buff_e + (i  )*inc_e;
            float*   alpha22  = buff_d + (i+1)*inc_d;

            float*   gammaL   = &(buff_g[(i  )*inc_g].real);
            float*   sigmaL   = &(buff_g[(i  )*inc_g].imag);
            float*   gammaR   = &(buff_h[(i  )*inc_h].real);
            float*   sigmaR   = &(buff_h[(i  )*inc_h].imag);

            a11temp = *alpha11 * cs;
            MAC_Givens2_ops( &a11temp,
                             alpha12,
                             &cs,
                             &sn,
                             &r );

            if ( i > 0 ) *alpha01 = oldsn * r;

            a11temp = oldcs * r;
            a22temp = *alpha22 * sn;
            MAC_Givens2_ops( &a11temp,
                             &a22temp,
                             &oldcs,
                             &oldsn,
                             alpha11 );

            *gammaR = cs;
            *sigmaR = sn;
            *gammaL = oldcs;
            *sigmaL = oldsn;
        }

        temp    = *d_last * cs;
        *d_last = temp * oldcs;
        *e_last = temp * oldsn;

        return FLA_SUCCESS;
    }


    // Apply Givens rotations in forward order.
    for ( i = 0; i < m_A - 1; ++i )
    {
        float*    alpha01  = buff_e + (i-1)*inc_e;
        float*    alpha11  = buff_d + (i  )*inc_d;
        float*    alpha21  = &bulge;
        float*    alpha02  = &bulge;
        float*    alpha12  = buff_e + (i  )*inc_e;
        float*    alpha22  = buff_d + (i+1)*inc_d;
        float*    alpha13  = &bulge;
        float*    alpha23  = buff_e + (i+1)*inc_e;

        float*    gammaL   = &(buff_g[(i  )*inc_g].real);
        float*    sigmaL   = &(buff_g[(i  )*inc_g].imag);
        float*    gammaR   = &(buff_h[(i  )*inc_h].real);
        float*    sigmaR   = &(buff_h[(i  )*inc_h].imag);

        float     alpha01_new;
        float     alpha11_new;

        int       mn_ahead  = m_A - i - 2;

        /*------------------------------------------------------------*/

        if ( i == 0 )
        {
            float alpha11_temp;
            float alpha12_temp;

            // Induce an iteration that introduces the bulge (from the right).
            //alpha11_temp = *buff_d - shift;
            alpha11_temp = ( fabsf( *buff_d ) - shift ) *
                           ( signof( one, *buff_d ) + ( shift / *buff_d ) );
            alpha12_temp = *buff_e;

            // Compute a Givens rotation that introduces the bulge.
            MAC_Givens2_ops( &alpha11_temp,
                             &alpha12_temp,
                             gammaR,
                             sigmaR,
                             &alpha11_new );

            // Apply the bulge-introducting Givens rotation (from the right)
            // to the top-left 2x2 matrix.
            MAC_Apply_G_2x2_ops( gammaR,
                                 sigmaR,
                                 alpha11,
                                 alpha21,
                                 alpha12,
                                 alpha22 );
        }
        else
        {
            // Compute a new Givens rotation to push the bulge (from the
            // right).
            MAC_Givens2_ops( alpha01,
                             alpha02,
                             gammaR,
                             sigmaR,
                             &alpha01_new );

            // Apply the Givens rotation (from the right) to the 1x2 vector
            // from which it was computed, which annihilates alpha02.
            *alpha01 = alpha01_new;
            *alpha02 = 0.0F;

            // Apply the Givens rotation to the 2x2 matrix below the 1x2
            // vector that was just modified.
            MAC_Apply_G_2x2_ops( gammaR,
                                 sigmaR,
                                 alpha11,
                                 alpha21,
                                 alpha12,
                                 alpha22 );
        }

        // Compute a new Givens rotation to push the bulge (from the left).
        MAC_Givens2_ops( alpha11,
                         alpha21,
                         gammaL,
                         sigmaL,
                         &alpha11_new );

        // Apply the Givens rotation (from the left) to the 2x1 vector
        // from which it was computed, which annihilates alpha11.
        *alpha11 = alpha11_new;
        *alpha21 = 0.0F;

        if ( mn_ahead > 0 )
        {
            // Apply the Givens rotation (from the left) to the 2x2 submatrix
            // at alpha12.
            MAC_Apply_GT_2x2_ops( gammaL,
                                  sigmaL,
                                  alpha12,
                                  alpha22,
                                  alpha13,
                                  alpha23 );
        }
        else
        {
            // Apply the Givens rotation (from the left) to the last 2x1 vector
            // at alpha12.
            MAC_Apply_GT_2x1_ops( gammaL,
                                  sigmaL,
                                  alpha12,
                                  alpha22 );
        }

        /*------------------------------------------------------------*/
    }

    return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_francis_v_opd_var1( int       m_A,
                                       double    shift,
                                       dcomplex* buff_g, int inc_g,
                                       dcomplex* buff_h, int inc_h,
                                       double*   buff_d, int inc_d,
                                       double*   buff_e, int inc_e )
{
    double    one   = bl1_d1();
    double    bulge = 0.0;
    int       i;

    // If the shift is zero, perform a simplified Francis step.
    if ( shift == 0.0 )
    {
        double cs, oldcs;
        double sn, oldsn = 0;
        double r, temp;
        double a11temp, a22temp;

        double* d_last = buff_d + (m_A-1)*inc_d;
        double* e_last = buff_e + (m_A-2)*inc_e;

        cs    = one;
        oldcs = one;

        for ( i = 0; i < m_A - 1; ++i )
        {
            double*   alpha01  = buff_e + (i-1)*inc_e;
            double*   alpha11  = buff_d + (i  )*inc_d;
            double*   alpha12  = buff_e + (i  )*inc_e;
            double*   alpha22  = buff_d + (i+1)*inc_d;

            double*   gammaL   = &(buff_g[(i  )*inc_g].real);
            double*   sigmaL   = &(buff_g[(i  )*inc_g].imag);
            double*   gammaR   = &(buff_h[(i  )*inc_h].real);
            double*   sigmaR   = &(buff_h[(i  )*inc_h].imag);

            a11temp = *alpha11 * cs;
            MAC_Givens2_opd( &a11temp,
                             alpha12,
                             &cs,
                             &sn,
                             &r );

            if ( i > 0 ) *alpha01 = oldsn * r;

            a11temp = oldcs * r;
            a22temp = *alpha22 * sn;
            MAC_Givens2_opd( &a11temp,
                             &a22temp,
                             &oldcs,
                             &oldsn,
                             alpha11 );

            *gammaR = cs;
            *sigmaR = sn;
            *gammaL = oldcs;
            *sigmaL = oldsn;
        }

        temp    = *d_last * cs;
        *d_last = temp * oldcs;
        *e_last = temp * oldsn;

        return FLA_SUCCESS;
    }


    // Apply Givens rotations in forward order.
    for ( i = 0; i < m_A - 1; ++i )
    {
        double*   alpha01  = buff_e + (i-1)*inc_e;
        double*   alpha11  = buff_d + (i  )*inc_d;
        double*   alpha21  = &bulge;
        double*   alpha02  = &bulge;
        double*   alpha12  = buff_e + (i  )*inc_e;
        double*   alpha22  = buff_d + (i+1)*inc_d;
        double*   alpha13  = &bulge;
        double*   alpha23  = buff_e + (i+1)*inc_e;

        double*   gammaL   = &(buff_g[(i  )*inc_g].real);
        double*   sigmaL   = &(buff_g[(i  )*inc_g].imag);
        double*   gammaR   = &(buff_h[(i  )*inc_h].real);
        double*   sigmaR   = &(buff_h[(i  )*inc_h].imag);

        double    alpha01_new;
        double    alpha11_new;

        int       mn_ahead  = m_A - i - 2;

        /*------------------------------------------------------------*/

        if ( i == 0 )
        {
            double alpha11_temp;
            double alpha12_temp;

            // Induce an iteration that introduces the bulge (from the right).
            //alpha11_temp = *buff_d - shift;
            alpha11_temp = ( fabs( *buff_d ) - shift ) *
                           ( signof( one, *buff_d ) + ( shift / *buff_d ) );
            alpha12_temp = *buff_e;

            // Compute a Givens rotation that introduces the bulge.
            MAC_Givens2_opd( &alpha11_temp,
                             &alpha12_temp,
                             gammaR,
                             sigmaR,
                             &alpha11_new );

            // Apply the bulge-introducting Givens rotation (from the right)
            // to the top-left 2x2 matrix.
            MAC_Apply_G_2x2_opd( gammaR,
                                 sigmaR,
                                 alpha11,
                                 alpha21,
                                 alpha12,
                                 alpha22 );
        }
        else
        {
            // Compute a new Givens rotation to push the bulge (from the
            // right).
            MAC_Givens2_opd( alpha01,
                             alpha02,
                             gammaR,
                             sigmaR,
                             &alpha01_new );

            // Apply the Givens rotation (from the right) to the 1x2 vector
            // from which it was computed, which annihilates alpha02.
            *alpha01 = alpha01_new;
            *alpha02 = 0.0;

            // Apply the Givens rotation to the 2x2 matrix below the 1x2
            // vector that was just modified.
            MAC_Apply_G_2x2_opd( gammaR,
                                 sigmaR,
                                 alpha11,
                                 alpha21,
                                 alpha12,
                                 alpha22 );
        }

        // Compute a new Givens rotation to push the bulge (from the left).
        MAC_Givens2_opd( alpha11,
                         alpha21,
                         gammaL,
                         sigmaL,
                         &alpha11_new );

        // Apply the Givens rotation (from the left) to the 2x1 vector
        // from which it was computed, which annihilates alpha11.
        *alpha11 = alpha11_new;
        *alpha21 = 0.0;

        if ( mn_ahead > 0 )
        {
            // Apply the Givens rotation (from the left) to the 2x2 submatrix
            // at alpha12.
            MAC_Apply_GT_2x2_opd( gammaL,
                                  sigmaL,
                                  alpha12,
                                  alpha22,
                                  alpha13,
                                  alpha23 );
        }
        else
        {
            // Apply the Givens rotation (from the left) to the last 2x1 vector
            // at alpha12.
            MAC_Apply_GT_2x1_opd( gammaL,
                                  sigmaL,
                                  alpha12,
                                  alpha22 );
        }

        /*------------------------------------------------------------*/
    }

    return FLA_SUCCESS;
}

