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
  GEBRD reduces a general complex M-by-N matrix A to upper or lower
  bidiagonal form B by a unitary transformation: Q**H * A * P = B.
  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.

  The interface of this function is different from LAPACK as
  FLA_Bidiag_UT does not produce real (sub)diagonals. LAPACK should
  be fixed to use complex datatypes for those diagonals.
*/

extern fla_bidiagut_t* fla_bidiagut_cntl_plain;

#define LAPACK_gebrd(prefix)                                            \
  int F77_ ## prefix ## gebrd( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_d,   \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_e,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_tu,  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_tv,  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                               int* info)

#define LAPACK_gebrd_body(prefix)                                       \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                \
  FLA_Datatype dtype_re = PREFIX2FLAME_REALTYPE(prefix);                \
  dim_t        min_m_n  = min( *m, *n );                                \
  dim_t        m_d      = min_m_n;                                      \
  dim_t        m_e      = min_m_n - 1;                                  \
  dim_t        m_t      = min_m_n;                                      \
  FLA_Obj      A, d, e, tu, tv, TU, TV, alpha;                          \
  FLA_Error    init_result;                                             \
  FLA_Uplo     uplo;                                                    \
  int          apply_scale;                                             \
                                                                        \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  uplo = ( *m >= *n ? FLA_UPPER_TRIANGULAR : FLA_LOWER_TRIANGULAR );    \
                                                                        \
  FLA_Obj_create_without_buffer( dtype_re, m_d, 1, &d );                \
  FLA_Obj_attach_buffer( buff_d, 1, m_d, &d );                          \
                                                                        \
  FLA_Obj_create_without_buffer( dtype_re, m_e, 1, &e );                \
  if ( m_e > 0 ) FLA_Obj_attach_buffer( buff_e, 1, m_e, &e );           \
                                                                        \
  /* m_t is assumed to be same although it is different */              \
  FLA_Obj_create_without_buffer( datatype, m_t, 1, &tu );               \
  FLA_Obj_attach_buffer( buff_tu, 1, m_t, &tu );                        \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, m_t, 1, &tv );               \
  FLA_Obj_attach_buffer( buff_tv, 1, m_t, &tv );                        \
                                                                        \
  FLA_Obj_create( dtype_re, 1, 1, 0, 0, &alpha );                       \
  FLA_Max_abs_value( A, alpha );                                        \
                                                                        \
  apply_scale =                                                         \
    ( FLA_Obj_gt( alpha, FLA_OVERFLOW_SQUARE_THRES  ) == TRUE ) -       \
    ( FLA_Obj_lt( alpha, FLA_UNDERFLOW_SQUARE_THRES ) == TRUE );        \
                                                                        \
  if ( apply_scale )                                                    \
    FLA_Scal( apply_scale > 0 ? FLA_SAFE_MIN : FLA_SAFE_INV_MIN, A );   \
                                                                        \
  FLA_Bidiag_UT_create_T( A, &TU, &TV );                                \
  FLA_Set( FLA_ZERO, TU );FLA_Set( FLA_ZERO, TV );                      \
                                                                        \
  FLA_Bidiag_UT_internal( A, TU, TV, fla_bidiagut_cntl_plain );         \
                                                                        \
  if ( apply_scale )                                                    \
    FLA_Bidiag_UT_scale_diagonals( apply_scale < 0 ? FLA_SAFE_MIN : FLA_SAFE_INV_MIN, A ); \
                                                                        \
  if ( FLA_Obj_is_complex( A ) == TRUE ) {                              \
    FLA_Obj d2, e2, rL, rR;                                             \
                                                                        \
    /* Temporary vectors to store diagonal and subdiagonal */           \
    FLA_Obj_create( datatype, m_d, 1, 0, 0, &d2 );                      \
    if ( m_e > 0 ) FLA_Obj_create( datatype, m_e, 1, 0, 0, &e2 );       \
                                                                        \
    /* Temporary vectors to store realifying transformation */          \
    FLA_Obj_create( datatype, m_d, 1, 0, 0, &rL );                      \
    FLA_Obj_create( datatype, m_d, 1, 0, 0, &rR );                      \
                                                                        \
    /* Do not touch factors in A */                                     \
    FLA_Bidiag_UT_extract_diagonals( A, d2, e2 );                       \
    FLA_Bidiag_UT_realify_diagonals( uplo, d2, e2, rL, rR );            \
                                                                        \
    FLA_Obj_extract_real_part( d2, d );                                 \
    if ( m_e > 0 ) FLA_Obj_extract_real_part( e2, e );                  \
                                                                        \
    /* Clean up */                                                      \
    FLA_Obj_free( &rL );                                                \
    FLA_Obj_free( &rR );                                                \
    FLA_Obj_free( &d2 );                                                \
    if ( m_e > 0 ) FLA_Obj_free( &e2 );                                 \
  } else {                                                              \
    FLA_Bidiag_UT_extract_real_diagonals( A, d, e );                    \
  }                                                                     \
  FLA_Bidiag_UT_recover_tau( TU, TV, tu, tv );                          \
                                                                        \
  PREFIX2FLAME_INVERT_TAU(prefix,tu);                                   \
  PREFIX2FLAME_INVERT_TAU(prefix,tv);                                   \
                                                                        \
  FLA_Obj_free( &alpha );                                               \
  FLA_Obj_free( &TU );                                                  \
  FLA_Obj_free( &TV );                                                  \
                                                                        \
  FLA_Obj_free_without_buffer( &A );                                    \
  FLA_Obj_free_without_buffer( &d );                                    \
  FLA_Obj_free_without_buffer( &e );                                    \
  FLA_Obj_free_without_buffer( &tu );                                   \
  FLA_Obj_free_without_buffer( &tv );                                   \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;


LAPACK_gebrd(s)
{
    {
        LAPACK_RETURN_CHECK( sgebrd_check( m, n,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_tu, buff_tv,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gebrd_body(s)
    }
}
LAPACK_gebrd(d)
{
    {
        LAPACK_RETURN_CHECK( dgebrd_check( m, n,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_tu, buff_tv,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gebrd_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gebrd(c)
{
    {
        LAPACK_RETURN_CHECK( cgebrd_check( m, n,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_tu, buff_tv,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gebrd_body(c)
    }
}
LAPACK_gebrd(z)
{
    {
        LAPACK_RETURN_CHECK( zgebrd_check( m, n,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_tu, buff_tv,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_gebrd_body(z)
    }
}
#endif

#define LAPACK_gebd2(prefix)                                            \
  int F77_ ## prefix ## gebd2( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_d,   \
                               PREFIX2LAPACK_REALDEF(prefix)* buff_e,   \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_tu,  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_tv,  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_w,   \
                               int* info )

LAPACK_gebd2(s)
{
    {
        LAPACK_RETURN_CHECK( sgebd2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_tu, buff_tv,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gebrd_body(s)
    }
}
LAPACK_gebd2(d)
{
    {
        LAPACK_RETURN_CHECK( dgebd2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_tu, buff_tv,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gebrd_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_gebd2(c)
{
    {
        LAPACK_RETURN_CHECK( cgebd2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_tu, buff_tv,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gebrd_body(c)
    }
}
LAPACK_gebd2(z)
{
    {
        LAPACK_RETURN_CHECK( zgebd2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_tu, buff_tv,
                                           buff_w,
                                           info ) )
    }
    {
        LAPACK_gebrd_body(z)
    }
}
#endif

#endif
