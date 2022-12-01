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
  ORGBR generates one of the real orthogonal matrices Q or P**T
  determined by GEBRD when reducing a real matrix A to bidiagonal
  form: A = Q * B * P**T.  Q and P**T are defined as products of
  elementary reflectors H(i) or G(i) respectively.

  If VECT = 'Q', A is assumed to have been an M-by-K matrix, and
  Q is of order M: if m >= k, Q = H(1) H(2) . . . H(k) and ORGBR
  returns the first n columns of Q, where m >= n >= k; if m < k,
  Q = H(1) H(2) . . . H(m-1) and ORGBR returns Q as an
  M-by-M matrix.

  If VECT = 'P', A is assumed to have been a K-by-N matrix, and
  P**T is of order N: if k < n, P**T = G(k) . . . G(2) G(1) and
  ORGBR returns the first m rows of P**T, where n >= m >= k;
  if k >= n, P**T = G(n-1) . . . G(2) G(1) and ORGBR returns
  P**T as an N-by-N matrix.
*/

#define LAPACK_orgbr(prefix, name)                                      \
  int F77_ ## prefix ## name ## br( char* vect,                         \
                                    int*  m,                            \
                                    int*  n,                            \
                                    int*  k,                            \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, \
                                    int* lwork,                         \
                                    int* info )

// buff_t shoud not include any zero. if it has one, that is the right dimension to go.
#define LAPACK_orgbr_body(prefix)                                       \
  FLA_Datatype datatype   = PREFIX2FLAME_DATATYPE(prefix);              \
  FLA_Obj      A, ATL, ATR, ABL, ABR, A1, A2, Ah, T, TL, TR, t;         \
  FLA_Error    init_result;                                             \
  FLA_Uplo     uplo;                                                    \
  dim_t        m_A, n_A, m_t;                                           \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  m_A = *m; n_A = *n;                                                   \
                                                                        \
  /* A is assumed to m x k for 'Q' and k x n for 'Q'. */                \
  if        ( *vect == 'Q' ) {                                          \
    uplo = ( *m >= *k ? FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR );  \
    m_t = min( *m, *k );                                                \
  } else {/*( *vect == 'P' ) */                                         \
    uplo = ( *k >= *n ? FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR );  \
    m_t = min( *k, *n );                                                \
  }                                                                     \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, m_A, n_A, &A );              \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );                \
  if ( m_t > 0 ) {                                                      \
    FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );                        \
    PREFIX2FLAME_INVERT_TAU(prefix,t);                                  \
  }                                                                     \
                                                                        \
  FLA_Part_2x2( A, &ATL, &ATR,                                          \
                   &ABL, &ABR, m_t, m_t, FLA_TL );                      \
                                                                        \
  /* Accumulate house-holder vectors */                                 \
  if        ( *vect == 'Q' )    {  /* U */                              \
    FLA_Part_2x1( ATL, &A1,                                             \
                       &A2, ( uplo == FLA_LOWER_TRIANGULAR ), FLA_TOP );\
    FLA_Merge_2x1( A2,                                                  \
                   ABL, &Ah );                                          \
    FLA_Merge_2x1( ATL,                                                 \
                   ABL, &A2 );                                          \
    FLA_Part_2x1( t, &t,                                                \
                     &T, FLA_Obj_min_dim( Ah ), FLA_TOP );              \
                                                                        \
    FLA_Bidiag_UT_create_T( A2, &T, NULL );                             \
    FLA_Set( FLA_ZERO, T );                                             \
    FLA_Part_1x2( T, &TL, &TR, FLA_Obj_length( t ), FLA_LEFT );         \
    FLA_Accum_T_UT( FLA_FORWARD, FLA_COLUMNWISE, Ah, t, TL );           \
  } else {/*( *vect == 'P' )          V */                              \
    FLA_Part_1x2( ATL, &A1, &A2, ( uplo == FLA_UPPER_TRIANGULAR ), FLA_LEFT ); \
    FLA_Merge_1x2( A2, ATR, &Ah );                                      \
    FLA_Merge_1x2( ATL, ATR, &A2 );                                     \
    FLA_Part_2x1( t, &t,                                                \
                     &T, FLA_Obj_min_dim( Ah ), FLA_TOP );              \
                                                                        \
    FLA_Bidiag_UT_create_T( A2, NULL, &T );                             \
    FLA_Set( FLA_ZERO, T );                                             \
    FLA_Part_1x2( T, &TL, &TR, FLA_Obj_length( t ), FLA_LEFT );         \
    FLA_Accum_T_UT( FLA_FORWARD, FLA_ROWWISE, Ah, t, TL );              \
  }                                                                     \
  if ( m_t > 0 ) {                                                      \
    PREFIX2FLAME_INVERT_TAU(prefix,t);                                  \
  }                                                                     \
                                                                        \
  /* Form Q or P^T ( U or V ) */                                        \
  if ( FLA_Obj_is_complex( A ) == TRUE && m_t > 0 ) {                   \
    FLA_Obj d2, e2, rL, rR;                                             \
                                                                        \
    /* Temporary vectors to store diagonal and subdiagonal */           \
    FLA_Obj_create( datatype, m_t, 1, 0, 0, &d2 );                      \
    if ( m_t > 1 ) FLA_Obj_create( datatype, m_t - 1, 1, 0, 0, &e2 );   \
                                                                        \
    /* Temporary vectors to store realifying transformation */          \
    FLA_Obj_create( datatype, m_t, 1, 0, 0, &rL );                      \
    FLA_Obj_create( datatype, m_t, 1, 0, 0, &rR );                      \
                                                                        \
    /* Extract diagonals (complex) and realify them. */                 \
    /* This is tricky as the shape of A is explicitly */                \
    /* assumed by m, n, k. */                                           \
    if ( uplo == FLA_UPPER_TRIANGULAR )                                 \
      FLA_Bidiag_UT_u_extract_diagonals( A2, d2, e2 );                  \
    else                                                                \
      FLA_Bidiag_UT_l_extract_diagonals( A2, d2, e2 );                  \
    FLA_Bidiag_UT_realify_diagonals( uplo, d2, e2, rL, rR );            \
                                                                        \
    if        ( *vect == 'Q' ) {                                        \
      /* Overwrite A to compute Q */                                    \
      FLA_Bidiag_UT_form_U_ext( uplo, A, T, FLA_NO_TRANSPOSE, A );      \
                                                                        \
      /* Applying rL */                                                 \
      FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE, rL, A2 );        \
    } else {/*( *vect == 'P' ) */                                       \
      /* Overwrite A to compute P */                                    \
      FLA_Bidiag_UT_form_V_ext( uplo, A, T, FLA_CONJ_TRANSPOSE, A );    \
                                                                        \
      /* Applying rR */                                                 \
      FLA_Apply_diag_matrix( FLA_LEFT, FLA_CONJUGATE, rR, A2 );         \
    }                                                                   \
                                                                        \
    /* Clean up */                                                      \
    FLA_Obj_free( &rR );                                                \
    FLA_Obj_free( &rL );                                                \
    if ( m_t > 1 ) FLA_Obj_free( &e2 );                                 \
    FLA_Obj_free( &d2 );                                                \
  } else {                                                              \
    if        ( *vect == 'Q' ) {                                        \
      FLA_Bidiag_UT_form_U_ext( uplo, A, T, FLA_NO_TRANSPOSE, A );      \
    } else {/*( *vect == 'P' ) */                                       \
      FLA_Bidiag_UT_form_V_ext( uplo, A, T, FLA_CONJ_TRANSPOSE, A );    \
    }                                                                   \
  }                                                                     \
                                                                        \
  FLA_Obj_free( &T );                                                   \
  FLA_Obj_free_without_buffer( &t );                                    \
  FLA_Obj_free_without_buffer( &A );                                    \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;

LAPACK_orgbr(s, org)
{
    {
        LAPACK_RETURN_CHECK( sorgbr_check( vect,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orgbr_body(s)
    }
}
LAPACK_orgbr(d, org)
{
    {
        LAPACK_RETURN_CHECK( dorgbr_check( vect,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orgbr_body(d)
    }
}
#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_orgbr(c, ung)
{
    {
        LAPACK_RETURN_CHECK( cungbr_check( vect,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orgbr_body(c)
    }
}
LAPACK_orgbr(z, ung)
{
    {
        LAPACK_RETURN_CHECK( zungbr_check( vect,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_orgbr_body(z)
    }
}
#endif

#endif
