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
  SPFFRT2 computes partial or incomplete LDL' factorization of
  symmetric matrix in packed storage format.
*/

#define LAPACK_spffrt2(prefix)                                           \
  int F77_ ## prefix ## spffrt2( PREFIX2LAPACK_TYPEDEF(prefix)* buff_AP, \
                                 integer* n,                                 \
                                 integer* ncolm,                             \
                                 PREFIX2LAPACK_TYPEDEF(prefix)* work,    \
                                 PREFIX2LAPACK_TYPEDEF(prefix)* work2 )

#define LAPACK_spffrt2_body(prefix)                                      \
           prefix ## spffrt2_fla( buff_AP,                               \
                                  n, ncolm,                              \
                                  work, work2);


LAPACK_spffrt2(s)
{
    {
        LAPACK_RETURN_CHECK( sspffrt2_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ) )
    }
    {
        LAPACK_spffrt2_body(s)
    }
}
LAPACK_spffrt2(d)
{
    {
        LAPACK_RETURN_CHECK( dspffrt2_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ) )
    }
    {
        LAPACK_spffrt2_body(d)
    }
}
LAPACK_spffrt2(c)
{
    {
        LAPACK_RETURN_CHECK( cspffrt2_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ) )
    }
    {
        LAPACK_spffrt2_body(c)
    }
}
LAPACK_spffrt2(z)
{
    {
        LAPACK_RETURN_CHECK( zspffrt2_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ) )
    }
    {
        LAPACK_spffrt2_body(z)
    }
}

#endif
