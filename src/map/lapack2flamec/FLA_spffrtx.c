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

extern void sspffrtx_fla(float *ap, integer *n, integer * ncolm, float *work, float *work2);
extern void dspffrtx_fla(double *ap, integer *n, integer * ncolm, double *work, double *work2);
extern void cspffrtx_fla(scomplex *ap, integer *n, integer * ncolm, scomplex *work, scomplex *work2);
extern void zspffrtx_fla(dcomplex *ap, integer *n, integer * ncolm, dcomplex *work, dcomplex *work2);
extern int sspffrtx_check(float *ap, integer *n, integer * ncolm, float *work, float *work2);
extern int dspffrtx_check(double *ap, integer *n, integer * ncolm, double *work, double *work2);
extern int cspffrtx_check(scomplex *ap, integer *n, integer * ncolm, scomplex *work, scomplex *work2);
extern int zspffrtx_check(dcomplex *ap, integer *n, integer * ncolm, dcomplex *work, dcomplex *work2);

#define LAPACK_spffrtx(prefix)                                           \
  int F77_ ## prefix ## spffrtx( PREFIX2LAPACK_TYPEDEF(prefix)* buff_AP, \
                                 integer* n,                                 \
                                 integer* ncolm,                             \
                                 PREFIX2LAPACK_TYPEDEF(prefix)* work,    \
                                 PREFIX2LAPACK_TYPEDEF(prefix)* work2 )

#define LAPACK_spffrtx_body(prefix)                                      \
           prefix ## spffrtx_fla( buff_AP,                               \
                                  n, ncolm,                              \
                                  work, work2);                          \



LAPACK_spffrtx(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sspffrtx inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1( sspffrtx_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_spffrtx_body(s)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_spffrtx(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dspffrtx inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1( dspffrtx_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrtx_body(d)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_spffrtx(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cspffrtx inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1( cspffrtx_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrtx_body(c)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_spffrtx(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zspffrtx inputs: n %" FLA_IS ", ncolm %" FLA_IS "", *n, *ncolm);
    {
        LAPACK_RETURN_CHECK_VAR1( zspffrtx_check( buff_AP,
                                             n, ncolm,
                                             work, work2 ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_spffrtx_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
 
#endif
