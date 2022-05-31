/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
    Copyright (c) 2021 Advanced Micro Devices, Inc. All rights reserved.
    Aug 4, 2021
*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME 

#define FLA_AMD_SVD2_OPT 0

#include "FLA_lapack2flame_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  GESVD computes the singular value decomposition (SVD) of a M-by-N
  matrix A, optionally computing the left and/or right singular vectors.
  The SVD is written
  A = U * S * transpose(V)
  where S is an M-by-N matrix which is zero except for its min(m,n)
  diagonal elements, U is an M-by-M orthogonal matrix, and V is an N-by-N
  orthogonal matrix.  The diagonal elements of S are the singular values
  of A; they are real and non-negative, and are returned in descending order.
  The first min(m,n) columns of U and V are the left and right singular
  vectors of A.

  Note that the routine returns V**T, not V.
*/

#define LAPACK_gesvd_real(prefix)                                       \
  int F77_ ## prefix ## gesvd( char* jobu,                              \
                               char* jobv,                              \
                               integer*  m,                                 \
                               integer*  n,                                 \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A,  integer* ldim_A, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_s,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_U,  integer* ldim_U, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_Vh, integer* ldim_Vh, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,  integer* lwork, \
                               integer* info )

#define LAPACK_gesvd_complex(prefix)                                    \
  int F77_ ## prefix ## gesvd( char* jobu,                              \
                               char* jobv,                              \
                               integer*  m,                                 \
                               integer*  n,                                 \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A,  integer* ldim_A, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_s,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_U,  integer* ldim_U, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_Vh, integer* ldim_Vh, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,  integer* lwork, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_r,   \
                               integer* info )


#define LAPACK_gesvd_body_d(prefix)                                                                                \
  if(FLA_AMD_SVD2_OPT && *m==2 && *n==2)  {                     \
  dgesvd2x2( jobu, jobv, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U,buff_Vh , ldim_Vh, buff_w, lwork,info );} \
  else {                                                                                \
  lapack_dgesvd ( jobu, jobv, m, n, buff_A, ldim_A, buff_s, buff_U, ldim_U,buff_Vh , ldim_Vh, buff_w, lwork,info );}   \
                                                                                                     

#define LAPACK_gesvd_body(prefix)                                       \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                \
  FLA_Datatype dtype_re = PREFIX2FLAME_REALTYPE(prefix);                \
  dim_t        min_m_n  = min( *m, *n );                                \
  FLA_Svd_type jobu_fla;                                                \
  FLA_Svd_type jobv_fla;                                                \
  FLA_Bool     create_U;                                                \
  FLA_Bool     create_V;                                                \
  dim_t        m_U, n_U;                                                \
  dim_t        m_V, n_V;                                                \
  FLA_Obj      A, s, U, V;                                              \
  FLA_Error    e_val, init_result;                                      \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  /* Parameters */                                                      \
  FLA_Param_map_netlib_to_flame_svd_type( jobu, &jobu_fla );            \
  FLA_Param_map_netlib_to_flame_svd_type( jobv, &jobv_fla );            \
                                                                        \
  m_U = *m; n_U = ( jobu_fla == FLA_SVD_VECTORS_ALL ? *m : min_m_n );   \
  n_V = *n; m_V = ( jobv_fla == FLA_SVD_VECTORS_ALL ? *n : min_m_n );   \
                                                                        \
  create_U = ( jobu_fla == FLA_SVD_VECTORS_ALL         ||               \
               jobu_fla == FLA_SVD_VECTORS_MIN_COPY );                  \
  create_V = ( jobv_fla == FLA_SVD_VECTORS_ALL         ||               \
               jobv_fla == FLA_SVD_VECTORS_MIN_COPY );                  \
                                                                        \
  /* Given A */                                                         \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  /* Singular values are stored in s */                                 \
  FLA_Obj_create_without_buffer( dtype_re, min_m_n, 1, &s );            \
  FLA_Obj_attach_buffer( buff_s, 1, min_m_n, &s );                      \
                                                                        \
  /* U */                                                               \
  if ( create_U ) {                                                     \
    FLA_Obj_create_without_buffer( datatype, m_U, n_U, &U );            \
    FLA_Obj_attach_buffer( buff_U, 1, *ldim_U, &U );                    \
  } else {                                                              \
    FLA_Obj_nullify( &U );                                              \
  }                                                                     \
  /* V^H */                                                             \
  if ( create_V ) {                                                     \
    FLA_Obj_create_without_buffer( datatype, m_V, n_V, &V );            \
    FLA_Obj_attach_buffer( buff_Vh, 1, *ldim_Vh, &V );                  \
  } else {                                                              \
    FLA_Obj_nullify( &V );                                              \
  }                                                                     \
  /* Compute SVD */                                                     \
  e_val = FLA_Svd_ext( jobu_fla, FLA_NO_TRANSPOSE,                      \
                       jobv_fla, FLA_CONJ_TRANSPOSE,                    \
                       A, s, U, V );                                    \
                                                                        \
  /* Clean up */                                                        \
  if ( create_U ) FLA_Obj_free_without_buffer( &U );                    \
  if ( create_V ) FLA_Obj_free_without_buffer( &V );                    \
                                                                        \
  FLA_Obj_free_without_buffer( &A );                                    \
  FLA_Obj_free_without_buffer( &s );                                    \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
  *info = 0;                                                            \
                                                                        \



LAPACK_gesvd_real(s)
{
  int fla_error = LAPACK_SUCCESS;
  AOCL_DTL_TRACE_LOG_INIT
  AOCL_DTL_SNPRINTF("sgesvd inputs: jobu %c, jobvt %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldu %" FLA_IS ", ldvt %" FLA_IS "", *jobu, *jobv, *m, *n, *ldim_A, *ldim_U, *ldim_Vh);
  {
      LAPACK_RETURN_CHECK_VAR1(sgesvd_check(jobu, jobv,
                                       m, n,
                                       buff_A, ldim_A,
                                       buff_s,
                                       buff_U, ldim_U,
                                       buff_Vh, ldim_Vh,
                                       buff_w, lwork,
                                       info),fla_error)
  }
  if(fla_error == LAPACK_SUCCESS)
  {
    LAPACK_gesvd_body(s)
        /** fla_error set to e_val on LAPACK_SUCCESS */
        fla_error = e_val;
  }
  AOCL_DTL_TRACE_LOG_EXIT
  return fla_error;
}

LAPACK_gesvd_real(d)
{
  int fla_error = LAPACK_SUCCESS;
  AOCL_DTL_TRACE_LOG_INIT
  AOCL_DTL_SNPRINTF("dgesvd inputs: jobu %c, jobvt %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldu %" FLA_IS ", ldvt %" FLA_IS "", *jobu, *jobv, *m, *n, *ldim_A, *ldim_U, *ldim_Vh);
#if FLA_AMD_OPT
    {
      LAPACK_gesvd_body_d(d)
          /** fla_error set to *info on LAPACK_SUCCESS */
          fla_error = *info;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
       LAPACK_RETURN_CHECK_VAR1( dgesvd_check (  jobu, jobv,
                                            m, n,
                                            buff_A, ldim_A,
                                            buff_s,
                                            buff_U, ldim_U,
                                            buff_Vh, ldim_Vh,
                                            buff_w, lwork,
                                            info ),fla_error )

    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gesvd_body(d)
        /** fla_error set to e_val on LAPACK_SUCCESS */
      fla_error = e_val;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#endif
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gesvd_complex(c)
{
  int fla_error = LAPACK_SUCCESS;
  AOCL_DTL_TRACE_LOG_INIT
  AOCL_DTL_SNPRINTF("cgesvd inputs: jobu %c, jobvt %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldu %" FLA_IS ", ldvt %" FLA_IS "", *jobu, *jobv, *m, *n, *ldim_A, *ldim_U, *ldim_Vh);
    {
        LAPACK_RETURN_CHECK_VAR1( cgesvd_check( jobu, jobv,
                                           m, n,
                                           buff_A, ldim_A,
                                           buff_s,
                                           buff_U, ldim_U,
                                           buff_Vh, ldim_Vh,
                                           buff_w, lwork,
                                           buff_r,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gesvd_body(c)
        /** fla_error set to e_val on LAPACK_SUCCESS */
      fla_error = e_val;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_gesvd_complex(z)
{
  int fla_error = LAPACK_SUCCESS;
  AOCL_DTL_TRACE_LOG_INIT
  AOCL_DTL_SNPRINTF("zgesvd inputs: jobu %c, jobvt %c, m %" FLA_IS ", n %" FLA_IS ", lda %" FLA_IS ", ldu %" FLA_IS ", ldvt %" FLA_IS "", *jobu, *jobv, *m, *n, *ldim_A, *ldim_U, *ldim_Vh);
    {
        LAPACK_RETURN_CHECK_VAR1( zgesvd_check( jobu, jobv,
                                           m, n,
                                           buff_A, ldim_A,
                                           buff_s,
                                           buff_U, ldim_U,
                                           buff_Vh, ldim_Vh,
                                           buff_w, lwork,
                                           buff_r,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_gesvd_body(z)
       /** fla_error set to e_val on LAPACK_SUCCESS */
      fla_error = e_val;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif

#endif
