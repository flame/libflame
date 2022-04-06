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
                                  work, work2);                          \



LAPACK_spffrt2(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sspffrt2 inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1( sspffrt2_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrt2_body(s)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_spffrt2(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dspffrt2 inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1( dspffrt2_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_spffrt2_body(d)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_spffrt2(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cspffrt2 inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1( cspffrt2_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrt2_body(c)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_spffrt2(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zspffrt2 inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1( zspffrt2_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrt2_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}

#endif
