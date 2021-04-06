/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bsvd_compute_shift( FLA_Obj tol, FLA_Obj sminl, FLA_Obj smax, FLA_Obj d, FLA_Obj e, FLA_Obj shift )
{
    FLA_Datatype datatype;
    integer          m_A;
    integer          inc_d;
    integer          inc_e;

    datatype = FLA_Obj_datatype( d );

    m_A      = FLA_Obj_vector_dim( d );

    inc_d    = FLA_Obj_vector_inc( d );
    inc_e    = FLA_Obj_vector_inc( e );

    switch ( datatype )
    {
    case FLA_FLOAT:
    {
        float*    buff_tol       = FLA_FLOAT_PTR( tol );
        float*    buff_sminl     = FLA_FLOAT_PTR( sminl );
        float*    buff_smax      = FLA_FLOAT_PTR( smax );
        float*    buff_d         = FLA_FLOAT_PTR( d );
        float*    buff_e         = FLA_FLOAT_PTR( e );
        float*    buff_shift     = FLA_FLOAT_PTR( shift );

        FLA_Bsvd_compute_shift_ops( m_A,
                                    *buff_tol,
                                    *buff_sminl,
                                    *buff_smax,
                                    buff_d, inc_d,
                                    buff_e, inc_e,
                                    buff_shift );

        break;
    }

    case FLA_DOUBLE:
    {
        double*   buff_tol       = FLA_DOUBLE_PTR( tol );
        double*   buff_sminl     = FLA_DOUBLE_PTR( sminl );
        double*   buff_smax      = FLA_DOUBLE_PTR( smax );
        double*   buff_d         = FLA_DOUBLE_PTR( d );
        double*   buff_e         = FLA_DOUBLE_PTR( e );
        double*   buff_shift     = FLA_DOUBLE_PTR( shift );

        FLA_Bsvd_compute_shift_opd( m_A,
                                    *buff_tol,
                                    *buff_sminl,
                                    *buff_smax,
                                    buff_d, inc_d,
                                    buff_e, inc_e,
                                    buff_shift );

        break;
    }
    }

    return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_compute_shift_ops( integer       m_A,
                                      float     tol,
                                      float     sminl,
                                      float     smax,
                                      float*    buff_d, integer inc_d,
                                      float*    buff_e, integer inc_e,
                                      float*    shift )
{
    float  hndrth = 0.01;
    float  eps;
    float* d_first;
    float* e_last;
    float* d_last_m1;
    float* d_last;
    float  sll, temp;

    eps = FLA_Mach_params_ops( FLA_MACH_EPS );

    d_first   = buff_d + (0    )*inc_d;
    e_last    = buff_e + (m_A-2)*inc_e;
    d_last_m1 = buff_d + (m_A-2)*inc_d;
    d_last    = buff_d + (m_A-1)*inc_d;

    // If the shift would ruin relative accuracy, set it to zero.
    if ( m_A * tol * ( sminl / smax ) <= max( eps, hndrth * tol ) )
    {
        *shift = 0.0;
    }
    else
    {
        // Compute the shift from the last 2x2 matrix.
        FLA_Sv_2x2_ops( d_last_m1,
                        e_last,
                        d_last,
                        shift,
                        &temp );

        sll = fabsf( *d_first );

        // Check if the shift is negligible; if so, set it to zero.
        if ( sll > 0.0F )
        {
            temp = ( *shift / sll );
            if ( temp * temp < eps )
            {
                *shift = 0.0F;
            }
        }
    }

    return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_compute_shift_opd( integer       m_A,
                                      double    tol,
                                      double    sminl,
                                      double    smax,
                                      double*   buff_d, integer inc_d,
                                      double*   buff_e, integer inc_e,
                                      double*   shift )
{
    double  hndrth = 0.01;
    double  eps;
    double* d_first;
    double* e_last;
    double* d_last_m1;
    double* d_last;
    double  sll, temp;

    eps = FLA_Mach_params_opd( FLA_MACH_EPS );

    d_first   = buff_d + (0    )*inc_d;
    e_last    = buff_e + (m_A-2)*inc_e;
    d_last_m1 = buff_d + (m_A-2)*inc_d;
    d_last    = buff_d + (m_A-1)*inc_d;

    // If the shift would ruin relative accuracy, set it to zero.
    if ( m_A * tol * ( sminl / smax ) <= max( eps, hndrth * tol ) )
    {
#ifdef PRINTF
        printf( "FLA_Bsvd_compute_shift_opd: shift would ruin accuracy; setting shift to 0.\n" );
        printf( "                   m_A = %d     \n", m_A );
        printf( "                   tol = %20.15e\n", tol );
        printf( "                 sminl = %20.15e\n", sminl );
        printf( "                  smax = %20.15e\n", smax );
        printf( "                   LHS = %20.15e\n", m_A * tol * ( sminl / smax ) );
        printf( "      max(eps,0.01*tol)= %20.15e\n", max( eps, hndrth * tol ) );
#endif
        *shift = 0.0;
    }
    else
    {
        // Compute the shift from the last 2x2 matrix.
        FLA_Sv_2x2_opd( d_last_m1,
                        e_last,
                        d_last,
                        shift,
                        &temp );

        sll = fabs( *d_first );

        // Check if the shift is negligible; if so, set it to zero.
        if ( sll > 0.0 )
        {
            temp = ( *shift / sll );
            if ( temp * temp < eps )
            {
#ifdef PRINTF
                printf( "FLA_Bsvd_compute_shift_opd: shift is negligible; setting shift to 0.\n" );
#endif
                *shift = 0.0;
            }
        }
    }

    return FLA_SUCCESS;
}

