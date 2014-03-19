/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bsvd_find_max( FLA_Obj d, FLA_Obj e, FLA_Obj smax, FLA_Obj smin )
{
    FLA_Datatype datatype;
    int          m_A;
    int          inc_d;
    int          inc_e;

    datatype = FLA_Obj_datatype( d );

    m_A      = FLA_Obj_vector_dim( d );

    inc_d    = FLA_Obj_vector_inc( d );
    inc_e    = FLA_Obj_vector_inc( e );


    switch ( datatype )
    {
    case FLA_FLOAT:
    {
        float*    buff_d    = FLA_FLOAT_PTR( d );
        float*    buff_e    = FLA_FLOAT_PTR( e );
        float*    buff_smax = FLA_FLOAT_PTR( smax );
        float*    buff_smin = FLA_FLOAT_PTR( smin );

        FLA_Bsvd_find_max_min_ops( m_A,
                                   buff_d, inc_d,
                                   buff_e, inc_e,
                                   buff_smax,
                                   buff_smin );

        break;
    }

    case FLA_DOUBLE:
    {
        double*   buff_d    = FLA_DOUBLE_PTR( d );
        double*   buff_e    = FLA_DOUBLE_PTR( e );
        double*   buff_smax = FLA_DOUBLE_PTR( smax );
        double*   buff_smin = FLA_DOUBLE_PTR( smin );

        FLA_Bsvd_find_max_min_opd( m_A,
                                   buff_d, inc_d,
                                   buff_e, inc_e,
                                   buff_smax,
                                   buff_smin );

        break;
    }
    }

    return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_find_max_min_ops( int       m_A,
                                     float*    buff_d, int inc_d,
                                     float*    buff_e, int inc_e,
                                     float*    smax,
                                     float*    smin )
{
    float  smax_cand;
    float  smin_cand;
    int    i;

    smax_cand = fabsf( buff_d[ (m_A-1)*inc_d ] );
    smin_cand = smax_cand;

    for ( i = 0; i < m_A - 1; ++i )
    {
        float abs_di = fabsf( buff_d[ i*inc_d ] );
        float abs_ei = fabsf( buff_e[ i*inc_e ] );

        // Track the minimum element.
        smin_cand = min( smin_cand, abs_di );

        // Track the maximum element.
        smax_cand = max( smax_cand, abs_di );
        smax_cand = max( smax_cand, abs_ei );
    }

    // Save the results of the search.
    *smax = smax_cand;
    *smin = smin_cand;

    return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_find_max_min_opd( int       m_A,
                                     double*   buff_d, int inc_d,
                                     double*   buff_e, int inc_e,
                                     double*   smax,
                                     double*   smin )
{
    double smax_cand;
    double smin_cand;
    int    i;

    smax_cand = fabs( buff_d[ (m_A-1)*inc_d ] );
    smin_cand = smax_cand;

    for ( i = 0; i < m_A - 1; ++i )
    {
        double abs_di = fabs( buff_d[ i*inc_d ] );
        double abs_ei = fabs( buff_e[ i*inc_e ] );

        // Track the minimum element.
        smin_cand = min( smin_cand, abs_di );

        // Track the maximum element.
        smax_cand = max( smax_cand, abs_di );
        smax_cand = max( smax_cand, abs_ei );
    }

    // Save the results of the search.
    *smax = smax_cand;
    *smin = smin_cand;

    return FLA_SUCCESS;
}

