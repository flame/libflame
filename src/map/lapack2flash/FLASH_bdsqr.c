/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_YET_LAPACK2FLAME

#include "FLASH_lapack2flash_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  BDSQR computes the singular values and, optionally, the right and/or
  left singular vectors from the singular value decomposition (SVD) of
  a N-by-N (upper or lower) bidiagonal matrix B using the implicit
  zero-shift QR algorithm. The SVD of B has the form
       B = Q * S * P**T
  where S is the diagonal matrix of singular values, Q is an orthogonal
  matrix of left singular vectors, and P is an orthogonal matrix of
  right singular vectors. If left singular vectors are requested, this
  subroutine actually returns U*Q instead of Q, and, if right singular
  vectors are requested, this subroutine returns P**T*VT instead of
  P**T, for given real input matrices U and VT. When U and VT are the
  orthogonal matrices that reduce a general matrix A to bidiagonal
  form:  A = U*B*VT, as computed by SGEBRD, then
      A = (U*Q) * S * (P**T*VT)
  is the SVD of A. Optionally, the subroutine may also compute Q**T*C
  for a given real input matrix C.

  Detailed algorithm is described in the following reference:

  Restructuring the Tridiagonal and Bidiagonal QR Algorithms for
  Performance

  FIELD G. VAN ZEE, The University of Texas at Austin
  ROBERT A. VAN DE GEIJN, The University of Texas at Austin
  GREGORIO QUINTANA-ORT´ı, Universitat Jaume I
*/

#define LAPACK_bdsqr(prefix)                                            \
  int F77_ ## prefix ## bdsqr( char* uplo,                              \
                               int*  m_d,                               \
                               int*  n_Vt,                              \
                               int*  m_U,                               \
                               int*  n_C,                               \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_d,   \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_e,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_Vt, int* ldim_Vt, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_U,  int* ldim_U, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_C,  int* ldim_C, \
                               PREFIX2LAPACK_REALDEF(prefix)* rwork,    \
                               int*  info )

#define LAPACK_bdsqr_body(prefix)                                       \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                \
  FLA_Datatype dtype_re = PREFIX2FLAME_REALTYPE(prefix);                \
  FLA_Obj      d, e, U, Vt, G, H, C;                                    \
  FLA_Svd_type jobu, jobv;                                              \
  FLA_Uplo     uplo_fla;                                                \
  FLA_Error    init_result;                                             \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  /* Param mapping */                                                   \
  FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );                \
                                                                        \
  /* Create diagonals */                                                \
  FLA_Obj_create_without_buffer( dtype_re, *m_d, 1, &d );               \
  FLA_Obj_attach_buffer( buff_d, 1, *m_d, &d );                         \
  if ( *m_d > 1 ) {                                                     \
    FLA_Obj_create_without_buffer( dtype_re, *m_d - 1, 1, &e );         \
    FLA_Obj_attach_buffer( buff_d, 1, *m_d - 1, &e );                   \
  }                                                                     \
  /* Create U and Vt */                                                 \
  if ( *m_U > 0 ) {                                                     \
    FLA_Obj_create_without_buffer( datatype, *m_U, *m_d, &U  );         \
    FLA_Obj_attach_buffer( buff_U, 1, *ldim_U, &U );                    \
    jobu = FLA_SVD_VECTORS_ALL;                                         \
  } else {                                                              \
    jobu = FLA_SVD_VECTORS_NONE;                                        \
  }                                                                     \
  if ( *n_Vt > 0 ) {                                                    \
    FLA_Obj_create_without_buffer( datatype, *m_d, *n_Vt, &Vt );        \
    FLA_Obj_attach_buffer( buff_Vt, 1, *ldim_Vt, &Vt );                 \
    jobv = FLA_SVD_VECTORS_ALL;                                         \
    /* V should be flipped */                                           \
    FLA_Obj_flip_base( &Vt );                                           \
    FLA_Obj_flip_view( &Vt );                                           \
  } else {                                                              \
    jobv = FLA_SVD_VECTORS_NONE;                                        \
  }                                                                     \
  if ( *n_C > 0 ) {                                                     \
    FLA_Obj_create_without_buffer( datatype, *m_d, *n_C, &C );          \
    FLA_Obj_attach_buffer( buff_C, 1, *ldim_C, &C );                    \
  }                                                                     \
                                                                        \
  /* Create workspace */                                                \
  FLA_Bsvd_create_workspace( d, &G, &H );                               \
                                                                        \
  /* Perform Bidiagonal SVD */                                          \
  FLA_Bsvd_ext( uplo_fla, d, e, G, H,                                   \
                jobu, U,                                                \
                jobv, Vt,                                               \
                *n_C, C );                                              \
                                                                        \
  /* Free matrices */                                                   \
  if ( jobv != FLA_SVD_VECTORS_NONE ) {                                 \
    /* if V is flipped abd complex, then apply conjugate */             \
    if ( FLA_Obj_is_complex( Vt ) == TRUE )                             \
      FLA_Conjugate( Vt );                                              \
    FLA_Obj_free_without_buffer( &Vt );                                 \
  }                                                                     \
  if ( jobu != FLA_SVD_VECTORS_NONE ) FLA_Obj_free_without_buffer( &U ); \
  if ( *n_C > 0 )                     FLA_Obj_free_without_buffer( &C ); \
  FLA_Obj_free( &H );                                                   \
  FLA_Obj_free( &G );                                                   \
  FLA_Obj_free_without_buffer( &d );                                    \
  if ( *m_d > 1 )                                                       \
    FLA_Obj_free_without_buffer( &e );                                  \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;


LAPACK_bdsqr(s)
{
    {
        LAPACK_RETURN_CHECK( sbdsqr_check( uplo, m_d, n_Vt, m_U, n_C,
                                           buff_d, buff_e,
                                           buff_Vt, ldim_Vt,
                                           buff_U,  ldim_U,
                                           buff_C,  ldim_C,
                                           rwork,
                                           info ) )
    }
    {
        LAPACK_bdsqr_body(s)
    }
}
LAPACK_bdsqr(d)
{
    {
        LAPACK_RETURN_CHECK( dbdsqr_check( uplo, m_d, n_Vt, m_U, n_C,
                                           buff_d, buff_e,
                                           buff_Vt, ldim_Vt,
                                           buff_U,  ldim_U,
                                           buff_C,  ldim_C,
                                           rwork,
                                           info ) )
    }
    {
        LAPACK_bdsqr_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_bdsqr(c)
{
    {
        LAPACK_RETURN_CHECK( cbdsqr_check( uplo, m_d, n_Vt, m_U, n_C,
                                           buff_d, buff_e,
                                           buff_Vt, ldim_Vt,
                                           buff_U,  ldim_U,
                                           buff_C,  ldim_C,
                                           rwork,
                                           info ) )
    }
    {
        LAPACK_bdsqr_body(c)
    }
}
LAPACK_bdsqr(z)
{
    {
        LAPACK_RETURN_CHECK( zbdsqr_check( uplo, m_d, n_Vt, m_U, n_C,
                                           buff_d, buff_e,
                                           buff_Vt, ldim_Vt,
                                           buff_U,  ldim_U,
                                           buff_C,  ldim_C,
                                           rwork,
                                           info ) )
    }
    {
        LAPACK_bdsqr_body(z)
    }
}
#endif

#endif
