
#include "FLAME.h"

FLA_Error FLA_Bsvd_find_converged( FLA_Obj tol, FLA_Obj d, FLA_Obj e, FLA_Obj sminl )
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
        float*    buff_tol    = FLA_FLOAT_PTR( tol );
        float*    buff_d      = FLA_FLOAT_PTR( d );
        float*    buff_e      = FLA_FLOAT_PTR( e );
        float*    buff_sminl  = FLA_FLOAT_PTR( sminl );

        FLA_Bsvd_find_converged_ops( m_A,
                                     *buff_tol,
                                     buff_d, inc_d,
                                     buff_e, inc_e,
                                     buff_sminl );

        break;
    }

    case FLA_DOUBLE:
    {
        double*   buff_tol    = FLA_DOUBLE_PTR( tol );
        double*   buff_d      = FLA_DOUBLE_PTR( d );
        double*   buff_e      = FLA_DOUBLE_PTR( e );
        double*   buff_sminl  = FLA_DOUBLE_PTR( sminl );

        FLA_Bsvd_find_converged_opd( m_A,
                                     *buff_tol,
                                     buff_d, inc_d,
                                     buff_e, inc_e,
                                     buff_sminl );

        break;
    }
    }

    return FLA_SUCCESS;
}



FLA_Error FLA_Bsvd_find_converged_ops( int       m_A,
                                       float     tol,
                                       float*    buff_d, int inc_d,
                                       float*    buff_e, int inc_e,
                                       float*    sminl )
{
    float*  epsilon_last;
    float*  delta_last;
    float   mu;
    int     i;

    epsilon_last = buff_e + (m_A-2)*inc_e;
    delta_last   = buff_d + (m_A-1)*inc_d;

    // Check convergence at the bottom of the matrix first.
    if ( MAC_Bsvd_sinval_is_converged_ops( tol, *delta_last, *epsilon_last ) )
    {
        *sminl = 0.0F;
        return m_A - 2;
    }

    // If the last element is not yet converged, check interior elements.
    // Also, accumulate sminl for later use when it comes time to check
    // the shift.

    mu     = fabsf( *buff_d );
    *sminl = mu;

    for ( i = 0; i < m_A - 1; ++i )
    {
        float* epsilon1 = buff_e + (i  )*inc_e;
        float* delta2   = buff_d + (i+1)*inc_d;

        // Check convergence of epsilon1 against the value of mu accumulated
        // so far.
        if ( MAC_Bsvd_sinval_is_converged_ops( tol, mu, *epsilon1 ) )
        {
            return i;
        }

        // Update mu and sminl.
        mu     = fabsf( *delta2 ) * ( mu / ( mu + fabsf( *epsilon1 ) ) );
        *sminl = min( *sminl, mu );
    }

    // Return with no convergence found.
    return FLA_SUCCESS;
}

FLA_Error FLA_Bsvd_find_converged_opd( int       m_A,
                                       double    tol,
                                       double*   buff_d, int inc_d,
                                       double*   buff_e, int inc_e,
                                       double*   sminl )
{
    double* epsilon_last;
    double* delta_last;
    double  mu;
    int     i;

    epsilon_last = buff_e + (m_A-2)*inc_e;
    delta_last   = buff_d + (m_A-1)*inc_d;

    // Check convergence at the bottom of the matrix first.
    if ( MAC_Bsvd_sinval_is_converged_opd( tol, *delta_last, *epsilon_last ) )
    {
        //*epsilon_last = 0.0;
        *sminl = 0.0;
        return m_A - 2;
    }

    // If the last element is not yet converged, check interior elements.
    // Also, accumulate sminl for later use when it comes time to check
    // the shift.

    mu     = fabs( *buff_d );
    *sminl = mu;

    for ( i = 0; i < m_A - 1; ++i )
    {
        double* epsilon1 = buff_e + (i  )*inc_e;
        double* delta2   = buff_d + (i+1)*inc_d;

        // Check convergence of epsilon1 against the value of mu accumulated
        // so far.
        if ( MAC_Bsvd_sinval_is_converged_opd( tol, mu, *epsilon1 ) )
        {
//printf( "FLA_Bsvd_sinval_find_converged: Split occurred in col %d\n", i );
//printf( "FLA_Bsvd_sinval_find_converged: mu      alpha12 = %23.19e %23.19e\n", mu, *epsilon1 );
//printf( "FLA_Bsvd_sinval_find_converged:         alpha22 =         %43.19e\n", *delta2 );
            //*epsilon1 = 0.0;
            //return FLA_FAILURE;
            return i;
        }

        // Update mu and sminl.
        mu     = fabs( *delta2 ) * ( mu / ( mu + fabs( *epsilon1 ) ) );
        *sminl = min( *sminl, mu );
    }

    // Return with no convergence found.
    return FLA_SUCCESS;
}

