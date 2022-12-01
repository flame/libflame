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
  If VECT = 'Q', SORMBR overwrites the general real M-by-N matrix C
  with              SIDE = 'L'     SIDE = 'R'
  TRANS = 'N':      Q * C          C * Q
  TRANS = 'T':      Q**T * C       C * Q**T

  If VECT = 'P', SORMBR overwrites the general real M-by-N matrix C
  with              SIDE = 'L'     SIDE = 'R'
  TRANS = 'N':      P * C          C * P
  TRANS = 'T':      P**T * C       C * P**T

  Here Q and P**T are the orthogonal matrices determined by SGEBRD when
  reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
  P**T are defined as products of elementary reflectors H(i) and G(i)
  respectively.

  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
  order of the orthogonal matrix Q or P**T that is applied.

  If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
  if nq >= k, Q = H(1) H(2) . . . H(k);
  if nq < k, Q = H(1) H(2) . . . H(nq-1).

  If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
  if k < nq, P = G(1) G(2) . . . G(k);
  if k >= nq, P = G(1) G(2) . . . G(nq-1).

  Here dimenions m and n are defined w.r.t C.
*/

#define LAPACK_ormbr(prefix, name)                                      \
  int F77_ ## prefix ## name ## br( char* vect,                         \
                                    char* side,                         \
                                    char* trans,                        \
                                    int*  m,                            \
                                    int*  n,                            \
                                    int*  k,                            \
                                    PREFIX2LAPACK_TYPEDEF(prefix) *buff_A, int *ldim_A, \
                                    PREFIX2LAPACK_TYPEDEF(prefix) *buff_t, \
                                    PREFIX2LAPACK_TYPEDEF(prefix) *buff_C, int *ldim_C, \
                                    PREFIX2LAPACK_TYPEDEF(prefix) *buff_w, int *lwork, \
                                    int *info )

#define LAPACK_ormbr_body(prefix)                                       \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                \
  FLA_Side     side_fla;                                                \
  FLA_Trans    trans_fla;                                               \
  dim_t        nq, /* nw, */ m_t, mm, nn;                               \
  FLA_Obj      A, C, T, W, t;                                           \
  FLA_Obj      d2, e2, rL, rR;                                          \
  FLA_Uplo     uplo;                                                    \
  FLA_Error    init_result;                                             \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Param_map_netlib_to_flame_side( side, &side_fla );                \
  FLA_Param_map_netlib_to_flame_trans( trans, &trans_fla );             \
                                                                        \
  if      ( side_fla == FLA_LEFT )     {  nq = *m; /* nw = *n; */ }     \
  else /* ( side_fla == FLA_RIGHT) */  {  nq = *n; /* nw = *m; */ }     \
                                                                        \
  m_t = min( nq, *k );                                                  \
                                                                        \
  if      ( *vect == 'Q' )       { mm = nq; nn = *k; }                  \
  else /* ( *vect == 'P' ) */    { mm = *k; nn = nq; }                  \
                                                                        \
  /* Bidiag is assumed to be applied to */                              \
  /* A w.r.t. the following dimensions. */                              \
  FLA_Obj_create_without_buffer( datatype, mm, nn, &A );                \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  uplo = ( mm >= nn ? FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR );    \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, m_t, 1, &t );                \
  FLA_Obj_attach_buffer( buff_t, 1, m_t, &t );                          \
  PREFIX2FLAME_INVERT_TAU(prefix,t);                                    \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &C );                \
  FLA_Obj_attach_buffer( buff_C, 1, *ldim_C, &C );                      \
                                                                        \
  if ( FLA_Obj_is_complex( A ) == TRUE ) {                              \
    /* Temporary vectors to store diagonal and subdiagonal */           \
    FLA_Obj_create( datatype, m_t, 1, 0, 0, &d2 );                      \
    if ( m_t > 1 ) FLA_Obj_create( datatype, m_t - 1, 1, 0, 0, &e2 );   \
                                                                        \
    /* Temporary vectors to store realifying transformation */          \
    FLA_Obj_create( datatype, m_t, 1, 0, 0, &rL );                      \
    FLA_Obj_create( datatype, m_t, 1, 0, 0, &rR );                      \
                                                                        \
    /* Extract diagonals (complex) and realify them. */                 \
    FLA_Bidiag_UT_extract_diagonals( A, d2, e2 );                       \
    FLA_Bidiag_UT_realify_diagonals( uplo, d2, e2, rL, rR );            \
  }                                                                     \
                                                                        \
  if          ( *vect == 'Q'  ) {                                       \
    if ( mm < nn ) {                                                    \
      FLA_Part_2x1( A, &W, &A, 1, FLA_TOP );                            \
      if ( side_fla == FLA_LEFT )                                       \
        FLA_Part_2x1( C, &W, &C, 1, FLA_TOP );                          \
      else                                                              \
        FLA_Part_1x2( C, &W, &C, 1, FLA_LEFT );                         \
    }                                                                   \
    if ( FLA_Obj_min_dim( A ) > 0 ) {                                   \
      FLA_Part_1x2( A, &A, &W, FLA_Obj_min_dim( A ), FLA_LEFT );        \
      FLA_Part_2x1( t, &t,                                              \
                       &W, FLA_Obj_min_dim( A ), FLA_TOP );             \
      FLA_QR_UT_create_T( A, &T ); FLA_Set( FLA_ZERO, T );              \
      FLA_Apply_Q_UT_create_workspace_side( side_fla, T, C, &W);        \
      FLA_Accum_T_UT( FLA_FORWARD, FLA_COLUMNWISE, A, t, T );           \
                                                                        \
      if ( FLA_Obj_is_complex( A ) == TRUE ) {                          \
                                                                        \
        if      ( side_fla  == FLA_LEFT &&                              \
                  trans_fla == FLA_NO_TRANSPOSE )                       \
          FLA_Apply_diag_matrix( FLA_LEFT, FLA_CONJUGATE, rL, C );      \
        else if ( side_fla  == FLA_RIGHT &&                             \
                  trans_fla == FLA_CONJ_TRANSPOSE )                     \
          FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, rL, C );  \
                                                                        \
        FLA_Apply_Q_UT( side_fla, trans_fla, FLA_FORWARD, FLA_COLUMNWISE, \
                        A, T, W, C );                                   \
                                                                        \
        if      ( side_fla  == FLA_LEFT &&                              \
                  trans_fla == FLA_CONJ_TRANSPOSE )                     \
          FLA_Apply_diag_matrix( FLA_LEFT, FLA_NO_CONJUGATE, rL, C );   \
        else if ( side_fla  == FLA_RIGHT &&                             \
                  trans_fla == FLA_NO_TRANSPOSE )                       \
          FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE, rL, C );     \
                                                                        \
      } else {                                                          \
        FLA_Apply_Q_UT( side_fla, trans_fla, FLA_FORWARD, FLA_COLUMNWISE, \
                        A, T, W, C );                                   \
      }                                                                 \
                                                                        \
      FLA_Obj_free( &T );                                               \
      FLA_Obj_free( &W );                                               \
    }                                                                   \
  } else { /* ( *vect == 'P'  ) */                                      \
    if ( mm >= nn ) {                                                   \
      FLA_Part_1x2( A, &W, &A, 1, FLA_LEFT );                           \
      if ( side_fla == FLA_LEFT )                                       \
        FLA_Part_2x1( C, &W, &C, 1, FLA_TOP );                          \
      else                                                              \
        FLA_Part_1x2( C, &W, &C, 1, FLA_LEFT );                         \
    }                                                                   \
    if ( FLA_Obj_min_dim( A ) > 0 ) {                                   \
      FLA_Part_2x1( A, &A,                                              \
                       &W, FLA_Obj_min_dim( A ), FLA_TOP );             \
      FLA_Part_2x1( t, &t,                                              \
                       &W, FLA_Obj_min_dim( A ), FLA_TOP );             \
      FLA_LQ_UT_create_T( A, &T ); FLA_Set( FLA_ZERO, T );              \
      FLA_Apply_Q_UT_create_workspace_side( side_fla, T, C, &W);        \
      FLA_Accum_T_UT( FLA_FORWARD, FLA_ROWWISE, A, t, T );              \
                                                                        \
      if ( FLA_Obj_is_complex( A ) == TRUE ) {                          \
                                                                        \
        if      ( side_fla  == FLA_LEFT &&                              \
                  trans_fla == FLA_NO_TRANSPOSE )                       \
          FLA_Apply_diag_matrix( FLA_LEFT, FLA_CONJUGATE, rR, C );      \
        else if ( side_fla  == FLA_RIGHT &&                             \
                  trans_fla == FLA_CONJ_TRANSPOSE )                     \
          FLA_Apply_diag_matrix( FLA_RIGHT, FLA_NO_CONJUGATE, rR, C );  \
                                                                        \
        FLA_Apply_Q_UT( side_fla, trans_fla, FLA_BACKWARD, FLA_ROWWISE, \
                        A, T, W, C );                                   \
                                                                        \
        if      ( side_fla  == FLA_LEFT &&                              \
                  trans_fla == FLA_CONJ_TRANSPOSE )                     \
          FLA_Apply_diag_matrix( FLA_LEFT, FLA_NO_CONJUGATE, rR, C );   \
        else if ( side_fla  == FLA_RIGHT &&                             \
                  trans_fla == FLA_NO_TRANSPOSE )                       \
          FLA_Apply_diag_matrix( FLA_RIGHT, FLA_CONJUGATE, rR, C );     \
                                                                        \
      } else {                                                          \
        FLA_Apply_Q_UT( side_fla, trans_fla, FLA_BACKWARD, FLA_ROWWISE, \
                        A, T, W, C );                                   \
      }                                                                 \
                                                                        \
      FLA_Obj_free( &T );                                               \
      FLA_Obj_free( &W );                                               \
    }                                                                   \
  }                                                                     \
                                                                        \
  if ( FLA_Obj_is_complex( A ) == TRUE ) {                              \
    /* Clean up */                                                      \
    FLA_Obj_free( &rR );                                                \
    FLA_Obj_free( &rL );                                                \
    if ( m_t > 1 ) FLA_Obj_free( &e2 );                                 \
    FLA_Obj_free( &d2 );                                                \
  }                                                                     \
                                                                        \
  PREFIX2FLAME_INVERT_TAU(prefix,t);                                    \
  FLA_Obj_free_without_buffer( &t );                                    \
                                                                        \
  FLA_Obj_free_without_buffer( &A );                                    \
  FLA_Obj_free_without_buffer( &C );                                    \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;

LAPACK_ormbr(s, orm)
{
    {
        LAPACK_RETURN_CHECK( sormbr_check( vect, side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_C, ldim_C,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_ormbr_body(s)
    }
}
LAPACK_ormbr(d, orm)
{
    {
        LAPACK_RETURN_CHECK( dormbr_check( vect, side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_C, ldim_C,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_ormbr_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_ormbr(c, unm)
{
    {
        LAPACK_RETURN_CHECK( cunmbr_check( vect, side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_C, ldim_C,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_ormbr_body(c)
    }
}
LAPACK_ormbr(z, unm)
{
    {
        LAPACK_RETURN_CHECK( zunmbr_check( vect, side, trans,
                                           m, n, k,
                                           buff_A, ldim_A,
                                           buff_t,
                                           buff_C, ldim_C,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_ormbr_body(z)
    }
}
#endif

#endif
