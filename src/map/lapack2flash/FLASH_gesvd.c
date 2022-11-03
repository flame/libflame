/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLASH

#include "FLASH_lapack2flash_util_defs.h"
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
                               int*  m,                                 \
                               int*  n,                                 \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A,  int* ldim_A, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_s,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_U,  int* ldim_U, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_Vh, int* ldim_Vh, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,  int* lwork, \
                               int* info )

#define LAPACK_gesvd_complex(prefix)                                    \
  int F77_ ## prefix ## gesvd( char* jobu,                              \
                               char* jobv,                              \
                               int*  m,                                 \
                               int*  n,                                 \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A,  int* ldim_A, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_s,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_U,  int* ldim_U, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_Vh, int* ldim_Vh, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,  int* lwork, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_r,   \
                               int* info )

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
  return e_val;


LAPACK_gesvd_real(s)
{
    {
        LAPACK_RETURN_CHECK( sgesvd_check( jobu, jobv,
                                           m, n,
                                           buff_A, ldim_A,
                                           buff_s,
                                           buff_U, ldim_U,
                                           buff_Vh, ldim_Vh,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gesvd_body(s)
    }
}

LAPACK_gesvd_real(d)
{
    {
        LAPACK_RETURN_CHECK( dgesvd_check( jobu, jobv,
                                           m, n,
                                           buff_A, ldim_A,
                                           buff_s,
                                           buff_U, ldim_U,
                                           buff_Vh, ldim_Vh,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gesvd_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gesvd_complex(c)
{
    {
        LAPACK_RETURN_CHECK( cgesvd_check( jobu, jobv,
                                           m, n,
                                           buff_A, ldim_A,
                                           buff_s,
                                           buff_U, ldim_U,
                                           buff_Vh, ldim_Vh,
                                           buff_w, lwork,
                                           buff_r,
                                           info ) )
    }
    {
        LAPACK_gesvd_body(c)
    }
}
LAPACK_gesvd_complex(z)
{
    {
        LAPACK_RETURN_CHECK( zgesvd_check( jobu, jobv,
                                           m, n,
                                           buff_A, ldim_A,
                                           buff_s,
                                           buff_U, ldim_U,
                                           buff_Vh, ldim_Vh,
                                           buff_w, lwork,
                                           buff_r,
                                           info ) )
    }
    {
        LAPACK_gesvd_body(z)
    }
}
#endif

#endif
