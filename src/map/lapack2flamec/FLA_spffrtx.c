/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    Oct 24, 2020
*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  SPFFRTX computes partial or incomplete LDL' factorization of
  symmetric matrix in packed storage format.
*/

#define LAPACK_spffrtx(prefix)                                           \
  int F77_ ## prefix ## spffrtx( PREFIX2LAPACK_TYPEDEF(prefix)* buff_AP, \
                                 integer* n,                                 \
                                 integer* ncolm,                             \
                                 PREFIX2LAPACK_TYPEDEF(prefix)* work,    \
                                 PREFIX2LAPACK_TYPEDEF(prefix)* work2 )

#define LAPACK_spffrtx_body(prefix)                                      \
           prefix ## spffrtx_fla( buff_AP,                               \
                                  n, ncolm,                              \
                                  work, work2);


LAPACK_spffrtx(s)
{
    {
        LAPACK_RETURN_CHECK( sspffrtx_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ) )
    }
    {
        LAPACK_spffrtx_body(s)
    }
}
LAPACK_spffrtx(d)
{
    {
        LAPACK_RETURN_CHECK( dspffrtx_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ) )
    }
    {
        LAPACK_spffrtx_body(d)
    }
}
LAPACK_spffrtx(c)
{
    {
        LAPACK_RETURN_CHECK( cspffrtx_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ) )
    }
    {
        LAPACK_spffrtx_body(c)
    }
}
LAPACK_spffrtx(z)
{
    {
        LAPACK_RETURN_CHECK( zspffrtx_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ) )
    }
    {
        LAPACK_spffrtx_body(z)
    }
}
 
#endif
