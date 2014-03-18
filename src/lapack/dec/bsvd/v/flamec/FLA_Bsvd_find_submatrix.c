
#include "FLAME.h"


FLA_Error FLA_Bsvd_find_submatrix_ops( int       mn_A,
                                       int       ij_begin,
                                       float*    buff_d, int inc_d,
                                       float*    buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR )
{
    float  rzero = bl1_s0();
    int    ij_tl;
    int    ij_br;

    // Search for the first non-zero superdiagonal element starting at
    // the index specified by ij_begin.
    for ( ij_tl = ij_begin; ij_tl < mn_A - 1; ++ij_tl )
    {
        float* e1 = buff_e + (ij_tl  )*inc_e;

        // If we find a non-zero element, record it and break out of this
        // loop.
        if ( *e1 != rzero )
        {
            *ijTL = ij_tl;
            break;
        }
    }

    // If ij_tl was incremented all the way up to mn_A - 1, then we didn't
    // find any non-zeros.
    if ( ij_tl == mn_A - 1 )
    {
        return FLA_FAILURE;
    }

    // If we've gotten this far, then a non-zero superdiagonal element was
    // found. Now we must walk the remaining portion of the superdiagonal
    // to find the first zero element, or if one is not found, we simply
    // use the last element of the superdiagonal.
    for ( ij_br = ij_tl; ij_br < mn_A - 1; ++ij_br )
    {
        float* e1 = buff_e + (ij_br  )*inc_e;

        // If we find a zero element, record it and break out of this
        // loop.
        if ( *e1 == rzero )
        {
            break;
        }
    }

    // If a zero element was found, then ij_br should hold the index of
    // that element. If a zero element was not found, then ij_br should
    // hold mn_A - 1. Either way, we save the value and return success.
    *ijBR = ij_br;

    return FLA_SUCCESS;
}

//#define PRINTF

FLA_Error FLA_Bsvd_find_submatrix_opd( int       mn_A,
                                       int       ij_begin,
                                       double*   buff_d, int inc_d,
                                       double*   buff_e, int inc_e,
                                       int*      ijTL,
                                       int*      ijBR )
{
    double rzero = bl1_d0();
    int    ij_tl;
    int    ij_br;

    // Search for the first non-zero superdiagonal element starting at
    // the index specified by ij_begin.
    for ( ij_tl = ij_begin; ij_tl < mn_A - 1; ++ij_tl )
    {
        double* e1 = buff_e + (ij_tl  )*inc_e;

        // If we find a non-zero element, record it and break out of this
        // loop.
        if ( *e1 != rzero )
        {
#ifdef PRINTF
            printf( "FLA_Bsvd_find_submatrix_opd: found non-zero superdiagonal element\n" );
            printf( "                             e[%3d] = %22.19e\n", ij_tl, *e1 );
#endif
            *ijTL = ij_tl;
            break;
        }
    }

    // If ij_tl was incremented all the way up to mn_A - 1, then we didn't
    // find any non-zeros.
    if ( ij_tl == mn_A - 1 )
    {
#ifdef PRINTF
        printf( "FLA_Bsvd_find_submatrix_opd: no submatrices found.\n" );
#endif
        return FLA_FAILURE;
    }

    // If we've gotten this far, then a non-zero superdiagonal element was
    // found. Now we must walk the remaining portion of the superdiagonal
    // to find the first zero element, or if one is not found, we simply
    // use the last element of the superdiagonal.
    for ( ij_br = ij_tl; ij_br < mn_A - 1; ++ij_br )
    {
        double* e1 = buff_e + (ij_br  )*inc_e;

        // If we find a zero element, record it and break out of this
        // loop.
        if ( *e1 == rzero )
        {
#ifdef PRINTF
            printf( "FLA_Bsvd_find_submatrix_opd: found zero superdiagonal element\n" );
            printf( "                             e[%3d] = %22.19e\n", ij_br, *e1 );
#endif
            break;
        }
    }

    // If a zero element was found, then ij_br should hold the index of
    // that element. If a zero element was not found, then ij_br should
    // hold mn_A - 1. Either way, we save the value and return success.
    *ijBR = ij_br;

    return FLA_SUCCESS;
}

