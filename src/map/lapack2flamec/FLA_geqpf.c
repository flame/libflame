/******************************************************************************
* Copyright (C) 2023, Advanced Micro Devices, Inc. All rights reserved.
*******************************************************************************/
/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  - GEQPF computes a QR factorization with column pivoting of a
  real M-by-N matrix A: A*P = Q*R.
  This routine is deprecated and has been replaced by routine GEQP3.

  GEQPF is used for generalized eigenvalue decomposition;
  JPVT is not only output but also input. Tricky to handle this
  in flame pivoting. So, use the original implementation. It still
  use our QR.

  - GEQP3 computes a QR factorization with column pivoting of a
  matrix A: A*P = Q*R using Level 3 BLAS.

  FLA_QR_UT_piv uses level 3 BLAS.
*/

extern int sgeqpf_fla(integer *m, integer *n, real *a, integer *lda, integer *jpvt, real *tau, real *work, integer *info);
extern int dgeqpf_fla(integer *m, integer *n, doublereal *a, integer * lda, integer *jpvt, doublereal *tau, doublereal *work, integer *info);

// GEQPF
#define LAPACK_geqpf(prefix)                                            \
  int F77_ ## prefix ## geqpf( integer* m,                                   \
                               integer* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               integer* buff_p,                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,   \
                               integer* info )

// Notation for LAPACK column pvioting is not consistent to pivoting in LU.
// This does not perform pre-ordering when jpiv include non-zero pivots.
//
#define LAPACK_geqpf_body(prefix)                                       \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                \
  FLA_Obj      A, t, T, w, p, jpiv;                                     \
  dim_t        min_m_n  = fla_min( *m, *n );                                \
  FLA_Error    init_result;                                             \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, min_m_n, 1, &t );            \
  FLA_Obj_attach_buffer( buff_t, 1, min_m_n, &t );                      \
  FLA_Set( FLA_ZERO, t );                                               \
                                                                        \
  FLA_Obj_create_without_buffer( PREFIX2FLAME_REALTYPE(prefix), *n, 1, &w ); \
  FLA_Obj_attach_buffer( buff_w, 1, *n, &w );                           \
                                                                        \
  /* LAPACK pivot storage */                                            \
  FLA_Obj_create_without_buffer( FLA_INT, *n, 1, &jpiv );               \
  FLA_Obj_attach_buffer( buff_p, 1, *n, &jpiv );                        \
                                                                        \
  /* FLAME pivot storage */                                             \
  FLA_Obj_create( FLA_INT, *n, 1, 0, 0, &p );                           \
  FLA_Set( FLA_ZERO, p );                                               \
                                                                        \
  /* QR_UT_piv */                                                       \
  FLA_QR_UT_create_T( A, &T );                                          \
  FLA_Set( FLA_ZERO, T );                                               \
                                                                        \
  FLA_QR_UT_piv( A, T, w, p );                                          \
  FLA_QR_UT_recover_tau( T, t );                                        \
  PREFIX2FLAME_INVERT_TAU(prefix,t);                                    \
                                                                        \
  /* Transform FLAME column pivots to LAPACK pivots */                  \
  FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, jpiv );              \
                                                                        \
  /* Cleaning */                                                        \
  FLA_Obj_free_without_buffer( &A );                                    \
  FLA_Obj_free_without_buffer( &t );                                    \
  FLA_Obj_free_without_buffer( &w );                                    \
  FLA_Obj_free_without_buffer( &jpiv );                                 \
  FLA_Obj_free( &p );                                                   \
  FLA_Obj_free( &T );                                                   \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \

/* 
    LAPACK path is enabled for both {S,D}GEQPF when FLA_ENABLE_AMD_OPT
    is set to to fix the incorrect results and NANs observed while
    testing xGEQPF.
 */

LAPACK_geqpf(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgeqpf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
#if !FLA_ENABLE_AMD_OPT
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK_VAR1( sgeqpf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqpf_body(s)
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        sgeqpf_fla( m, n,
            buff_A, ldim_A,
            buff_p,
            buff_t,
            buff_w,
            info );
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#endif
}

LAPACK_geqpf(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgeqpf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);    
#if !FLA_ENABLE_AMD_OPT
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK_VAR1( dgeqpf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w,
                                           info ), fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqpf_body(d)
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        dgeqpf_fla( m, n,
            buff_A, ldim_A,
            buff_p,
            buff_t,
            buff_w,
            info );
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#endif  
}

#define LAPACK_geqpf_complex(prefix)                                    \
  int F77_ ## prefix ## geqpf( integer* m,                                  \
                               integer* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               integer* buff_p,                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,   \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_r,   \
                               integer* info )

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqpf_complex(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgeqpf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK_VAR1( cgeqpf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w,
                                           buff_r,
                                           info ), fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqpf_body(c)
            fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_geqpf_complex(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgeqpf inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK_VAR1( zgeqpf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w,
                                           buff_r,
                                           info ), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqpf_body(z)
            fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

// GEQP3
#define LAPACK_geqp3(prefix)                                            \
  int F77_ ## prefix ## geqp3( integer* m,                                  \
                               integer* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               integer* buff_p,                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, integer* lwork, \
                               integer* info )
LAPACK_geqp3(s)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("sgeqp3 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    extern int sgeqp3_fla(integer *m, integer *n, real *a, integer *lda, integer *jpvt, real *tau, real *work, integer *lwork, integer *info);

#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1(sgeqp3_check(m, n,
                                             buff_A, ldim_A,
                                             buff_p,
                                             buff_t,
                                             buff_w, lwork,
                                             info), fla_error)
    }
    {
        if( *lwork == -1 )
        {
            AOCL_DTL_TRACE_LOG_EXIT
            return 0;
        }
        for (int i = 0; i < *n; ++i) buff_p[i] = (i + 1);
        if( *m == 0 || *n == 0 )
        {
            AOCL_DTL_TRACE_LOG_EXIT
            return 0;
        }
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_geqpf_body(s)
            fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        sgeqp3_fla( m, n,
                    buff_A, ldim_A,
                    buff_p,
                    buff_t,
                    buff_w, lwork,
                    info ) ;
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}

LAPACK_geqp3(d)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dgeqp3 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *ldim_A);
    extern int dgeqp3_fla(integer *m, integer *n, doublereal *a, integer * lda, integer *jpvt, doublereal *tau, doublereal *work, integer *lwork, integer *info);

#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( dgeqp3_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w, lwork,
                                           info ),fla_error )
    }
    {
        if( *lwork == -1 )
        {
            AOCL_DTL_TRACE_LOG_EXIT
            return 0;
        }
        for (int i = 0; i < *n; ++i) buff_p[i] = (i + 1);
        if( *m == 0 || *n == 0 )
        {
            AOCL_DTL_TRACE_LOG_EXIT
            return 0;
        }
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_geqpf_body(d)
            fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        dgeqp3_fla( m, n,
                    buff_A, ldim_A,
                    buff_p,
                    buff_t,
                    buff_w, lwork,
                    info ) ;
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}


#define LAPACK_geqp3_complex(prefix)                                            \
  int F77_ ## prefix ## geqp3( integer* m,                                  \
                               integer* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               integer* buff_p,                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_t,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, integer* lwork, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_r,   \
                               integer* info )


#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_geqp3_complex(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cgeqp3 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *lda);
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK_VAR1( cgeqp3_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           buff_t,
                                           buff_w, lwork,
                                           buff_r,
                                           info ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_geqpf_body(c)
            fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_geqp3_complex(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zgeqp3 inputs: m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS "", *m, *n, *lda);
    {
        for ( int i=0; i<*n; ++i) buff_p[i] = (i+1);
    }
    {
        LAPACK_RETURN_CHECK_VAR1(zgeqp3_check(m, n,
                                              buff_A, ldim_A,
                                              buff_p,
                                              buff_t,
                                              buff_w, lwork,
                                              buff_r,
                                              info),fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_geqpf_body(z)
            fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#endif
