/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
    Copyright (c) 2021-2023 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  GEQRF computes a QR factorization of a M-by-N matrix A: A = Q * R.
*/

extern int dgeqrf_fla(integer *m, integer *n, doublereal *a, integer * lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
extern int sgeqrf_fla(integer *m, integer *n, real *a, integer *lda, real *tau, real *work, integer *lwork, integer *info);
extern int sgeqrfp_fla(integer *m, integer *n, real *a, integer *lda, real *tau, real *work, integer *lwork, integer *info);
extern int dgeqrfp_fla(integer *m, integer *n, doublereal *a, integer * lda, doublereal *tau, doublereal *work, integer *lwork, integer *info);
extern int sgeqr2p_fla(integer *m, integer *n, real *a, integer *lda, real *tau, real *work, integer *info);
extern int dgeqr2p_fla(integer *m, integer *n, doublereal *a, integer * lda, doublereal *tau, doublereal *work, integer *info);

extern void DTL_Trace(
		    uint8 ui8LogLevel,
		    uint8 ui8LogType,
		    const int8 *pi8FileName,
		    const int8 *pi8FunctionName,
		    uint32 ui32LineNumber,
		    const int8 *pi8Message);

#define FLA_ENABLE_ALT_PATHS  0

// GEQRF and GEQR2
#define LAPACK_geqrf(prefix)                                                          \
  int F77_ ## prefix ## geqrf(integer* m,                                             \
                              integer* n,                                             \
                              PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                              PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,                  \
                              PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, integer* lwork,  \
                              integer* info )

#define LAPACK_geqrf_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Obj      A, t, T;                                         \
  integer      min_m_n  = fla_min( *m, *n );                        \
  FLA_Error    init_result;                                     \
                                                                \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &t );            \
  FLA_Obj_attach_buffer( buff_t, 1, min_m_n, &t );                      \
                                                                        \
  FLA_Set( FLA_ZERO, t );                                               \
                                                                        \
  FLA_QR_UT_create_T( A, &T );                                          \
  FLA_QR_UT( A, T );                                                    \
  FLA_QR_UT_recover_tau( T, t );                                        \
  PREFIX2FLAME_INVERT_TAU(prefix,t);                                    \
                                                                        \
  FLA_Obj_free_without_buffer( &A );                                    \
  FLA_Obj_free_without_buffer( &t );                                    \
  FLA_Obj_free( &T );                                                   \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \


LAPACK_geqrf(s)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgeqrf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", lwork %" FLA_IS "", *m, *n, *ldim_A, *lwork);
#if !FLA_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( sgeqrf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(s)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
      sgeqrf_fla(m, n, buff_A, ldim_A, buff_t, buff_w, lwork, info);
      AOCL_DTL_TRACE_LOG_EXIT
      return 0;
    }
#endif
}
LAPACK_geqrf(d)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgeqrf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", lwork %" FLA_IS "", *m, *n, *ldim_A, *lwork);
#if !FLA_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(dgeqrf_check(m, n,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w, lwork,
                                              info),fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(d)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
      /* Initialize global context data */
      aocl_fla_init();

      dgeqrf_fla(m, n, buff_A, ldim_A, buff_t, buff_w, lwork, info);
      AOCL_DTL_TRACE_LOG_EXIT
      return 0;
    }
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqrf(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgeqrf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", lwork %" FLA_IS "", *m, *n, *ldim_A, *lwork);
    {
        LAPACK_RETURN_CHECK_VAR1( cgeqrf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(c)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_geqrf(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgeqrf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", lwork %" FLA_IS "", *m, *n, *ldim_A, *lwork);
    {
        LAPACK_RETURN_CHECK_VAR1(zgeqrf_check(m, n,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w, lwork,
                                              info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(z)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#define LAPACK_geqr2(prefix)                                                           \
  int F77_ ## prefix ## geqr2( integer* m,                                             \
                               integer* n,                                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,                  \
                               integer* info )

LAPACK_geqr2(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgeqr2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( sgeqr2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w,
                                           info ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(s)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_geqr2(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgeqr2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dgeqr2_check(m, n,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w,
                                              info),
                                 fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(d)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqr2(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgeqr2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(cgeqr2_check(m, n,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w,
                                              info),
                                 fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(c)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_geqr2(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgeqr2 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zgeqr2_check(m, n,
                                              buff_A, ldim_A,
                                              buff_t,
                                              buff_w,
                                              info),fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(z)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

// GEQRFP and GEQR2P
#define LAPACK_geqrfp(prefix)                                                          \
  int F77_ ## prefix ## geqrfp(integer* m,                                             \
                               integer* n,                                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, integer* lwork,  \
                               integer* info )
LAPACK_geqrfp(s)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgeqrfp inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
#if !FLA_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( sgeqrfp_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w, lwork,
                                            info ), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(s)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
      sgeqrfp_fla( m, n,
                   buff_A, ldim_A,
                   buff_t,
                   buff_w, lwork,
                   info );
      AOCL_DTL_TRACE_LOG_EXIT
      return 0;
    }
#endif
}
LAPACK_geqrfp(d)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgeqrfp inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
#if !FLA_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(dgeqrfp_check(m, n,
                                               buff_A, ldim_A,
                                               buff_t,
                                               buff_w, lwork,
                                               info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(d)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
      dgeqrfp_fla( m, n,
                   buff_A, ldim_A,
                   buff_t,
                   buff_w, lwork,
                   info );
      AOCL_DTL_TRACE_LOG_EXIT
      return 0;
    }
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqrfp(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgeqrfp inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( cgeqrfp_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w, lwork,
                                            info ), fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(c)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_geqrfp(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgeqrfp inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zgeqrfp_check(m, n,
                                               buff_A, ldim_A,
                                               buff_t,
                                               buff_w, lwork,
                                               info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(z)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#define LAPACK_geqr2p(prefix)                                                           \
  int F77_ ## prefix ## geqr2p( integer* m,                                             \
                                integer* n,                                             \
                                PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                                PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,                  \
                                PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,                  \
                                integer* info )

LAPACK_geqr2p(s)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgeqr2p inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
#if !FLA_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( sgeqr2p_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w,
                                            info ), fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(s)
      /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        sgeqr2p_fla( m, n,
                     buff_A, ldim_A,
                     buff_t,
                     buff_w,
                     info );
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}
LAPACK_geqr2p(d)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgeqr2p inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
#if !FLA_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(dgeqr2p_check(m, n,
                                               buff_A, ldim_A,
                                               buff_t,
                                               buff_w,
                                               info),
                                 fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(d)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        dgeqr2p_fla(m, n,
                    buff_A, ldim_A,
                    buff_t,
                    buff_w,
                    info);
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqr2p(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgeqr2p inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( cgeqr2p_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w,
                                            info ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(c)
       /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_geqr2p(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgeqr2p inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( zgeqr2p_check( m, n,
                                            buff_A, ldim_A,
                                            buff_t,
                                            buff_w,
                                            info ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_geqrf_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#endif
