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
  A is array, dimension (LDA,N). On entry,
  the hermitian matrix A.  If UPLO = 'U', the
  leading N-by-N upper triangular part of A
  contains the upper triangular part of the matrix A,
  and the strictly lower triangular part of A
  is not referenced. If UPLO = 'L', the leading N-by-N
  lower triangular part of A contains the lower
  triangular part of the matrix A, and the strictly
  upper triangular part of A is not referenced.
  On exit, if UPLO = 'U', the diagonal and first
  superdiagonal of A are overwritten by the corresponding
  elements of the tridiagonal matrix T, and the elements
  above the first superdiagonal, with the array TAU,
  represent the orthogonal matrix Q as a product of elementary
  reflectors; if UPLO = 'L', the diagonal and first subdiagonal
  of A are over-written by the corresponding elements of the
  tridiagonal matrix T, and the elements below the first
  subdiagonal, with the array TAU, represent the orthogonal
  matrix Q as a product of elementary reflectors.

                                                                \
  TODO:: To interface upper triangular, QL (or storing house
  holder vectors backward) is required.
*/

#define LAPACK_hetrd(prefix, name)                                      \
  int F77_ ## prefix ## name ## trd( char* uplo,                        \
                                     int*  m,                           \
                                     PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                                     PREFIX2LAPACK_REALDEF(prefix)* buff_d, \
                                     PREFIX2LAPACK_REALDEF(prefix)* buff_e, \
                                     PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                     PREFIX2LAPACK_TYPEDEF(prefix)* buff_w, int* lwork, \
                                     int*  info )

#define LAPACK_hetrd_body(prefix)                                     \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);              \
  FLA_Datatype dtype_re = PREFIX2FLAME_REALTYPE(prefix);              \
  dim_t        m_d      = *m;                                         \
  dim_t        m_e      = m_d - 1;                                    \
  FLA_Uplo     uplo_fla;                                              \
  FLA_Obj      A, d, e, t, T;                                         \
  FLA_Error    init_result;                                           \
                                                                      \
  FLA_Init_safe( &init_result );                                      \
                                                                      \
  FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );              \
                                                                      \
  FLA_Obj_create_without_buffer( datatype, *m, *m, &A );              \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                    \
                                                                      \
  FLA_Obj_create_without_buffer( dtype_re, m_d, 1, &d );              \
  FLA_Obj_attach_buffer( buff_d, 1, m_d, &d );                        \
                                                                      \
  if ( m_e > 0 ) {                                                      \
    FLA_Obj_create_without_buffer( dtype_re, m_e, 1, &e );              \
    FLA_Obj_attach_buffer( buff_e, 1, m_e, &e );                        \
                                                                        \
    FLA_Obj_create_without_buffer( datatype, m_e, 1, &t );              \
    FLA_Obj_attach_buffer( buff_t, 1, m_e, &t );                        \
  }                                                                     \
                                                                        \
  FLA_Tridiag_UT_create_T( A, &T );                                     \
  FLA_Set( FLA_ZERO, T );                                               \
  FLA_Tridiag_UT( uplo_fla, A, T );                                     \
                                                                        \
  if ( FLA_Obj_is_complex( A ) == TRUE && m_e > 0 ) {                   \
    FLA_Obj d2, e2, r;                                                  \
                                                                        \
    /* Temporary vectors to store the subidagonal */                    \
    FLA_Obj_create( datatype, m_d, 1, 0, 0, &d2 );                      \
    FLA_Obj_create( datatype, m_e, 1, 0, 0, &e2 );                      \
                                                                        \
    /* Temporary vectors to store realifying transformation */          \
    FLA_Obj_create( datatype, m_d, 1, 0, 0, &r );                       \
                                                                        \
    /* Do not touch factors in A */                                     \
    FLA_Tridiag_UT_extract_diagonals( uplo_fla, A, d2, e2 );            \
    FLA_Tridiag_UT_realify_subdiagonal( e2, r );                        \
                                                                        \
    FLA_Obj_extract_real_part( d2, d );                                 \
    FLA_Obj_extract_real_part( e2, e );                                 \
                                                                        \
    /* Clean up */                                                      \
    FLA_Obj_free( &r  );                                                \
    FLA_Obj_free( &e2 );                                                \
    FLA_Obj_free( &d2 );                                                \
  } else {                                                              \
    FLA_Tridiag_UT_extract_real_diagonals( uplo_fla, A, d, e );         \
  }                                                                     \
                                                                        \
  if ( m_e > 0 ) {                                                      \
    FLA_Tridiag_UT_recover_tau( T, t );                                 \
    PREFIX2FLAME_INVERT_TAU(prefix,t);                                  \
  }                                                                     \
  FLA_Obj_free( &T );                                                   \
                                                                        \
  if ( m_e > 0 ) {                                                      \
    FLA_Obj_free_without_buffer( &e );                                  \
    FLA_Obj_free_without_buffer( &t );                                  \
  }                                                                     \
  FLA_Obj_free_without_buffer( &d );                                    \
  FLA_Obj_free_without_buffer( &A );                                    \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
                                                                        \
  return 0;


// Original lapack implementation for upper triangular versions.
// Upper triangular versions are not yet implemented in libflame.
// Thus, those routines should be isolated from others.
extern int chetd2_fla(char *uplo, integer *n, complex       *a, integer *lda, real       *d__, real       *e, complex       *tau, integer *info);
extern int dsytd2_fla(char *uplo, integer *n, doublereal    *a, integer *lda, doublereal *d__, doublereal *e, doublereal    *tau, integer *info);
extern int ssytd2_fla(char *uplo, integer *n, real          *a, integer *lda, real       *d__, real       *e, real          *tau, integer *info);
extern int zhetd2_fla(char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau, integer *info);

extern int chetrd_fla(char *uplo, integer *n, complex       *a, integer *lda, real       *d__, real       *e, complex       *tau, complex       *work, integer *lwork, integer *info);
extern int dsytrd_fla(char *uplo, integer *n, doublereal    *a, integer *lda, doublereal *d__, doublereal *e, doublereal    *tau, doublereal    *work, integer *lwork, integer *info);
extern int ssytrd_fla(char *uplo, integer *n, real          *a, integer *lda, real       *d__, real       *e, real          *tau, real          *work, integer *lwork, integer *info);
extern int zhetrd_fla(char *uplo, integer *n, doublecomplex *a, integer *lda, doublereal *d__, doublereal *e, doublecomplex *tau, doublecomplex *work, integer *lwork, integer *info);

LAPACK_hetrd(s,sy)
{
    {
        if ( *uplo == 'U' )
        {
            ssytrd_fla( uplo, m,
                        buff_A, ldim_A,
                        buff_d, buff_e,
                        buff_t,
                        buff_w, lwork,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( ssytrd_check( uplo, m,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_hetrd_body(s)
    }
}
LAPACK_hetrd(d,sy)
{
    {
        if ( *uplo == 'U' )
        {
            dsytrd_fla( uplo, m,
                        buff_A, ldim_A,
                        buff_d, buff_e,
                        buff_t,
                        buff_w, lwork,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( dsytrd_check( uplo, m,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_hetrd_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_hetrd(c,he)
{
    {
        if ( *uplo == 'U' )
        {
            chetrd_fla( uplo, m,
                        (complex*)buff_A, ldim_A,
                        (real*)buff_d, (real*)buff_e,
                        (complex*)buff_t,
                        (complex*)buff_w, lwork,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( chetrd_check( uplo, m,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_hetrd_body(c)
    }
}
LAPACK_hetrd(z,he)
{
    {
        if ( *uplo == 'U' )
        {
            zhetrd_fla( uplo, m,
                        (doublecomplex*)buff_A, ldim_A,
                        (doublereal*)buff_d, (doublereal*)buff_e,
                        (doublecomplex*)buff_t,
                        (doublecomplex*)buff_w, lwork,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( zhetrd_check( uplo, m,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_t,
                                           buff_w, lwork,
                                           info ) )
    }
    {
        LAPACK_hetrd_body(z)
    }
}
#endif

#define LAPACK_hetd2(prefix, name)                                      \
  int F77_ ## prefix ## name ## td2( char* uplo,                        \
                                     int*  m,                           \
                                     PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                                     PREFIX2LAPACK_REALDEF(prefix)* buff_d, \
                                     PREFIX2LAPACK_REALDEF(prefix)* buff_e, \
                                     PREFIX2LAPACK_TYPEDEF(prefix)* buff_t, \
                                     int*  info )

LAPACK_hetd2(s,sy)
{
    {
        if ( *uplo == 'U' )
        {
            ssytd2_fla( uplo, m,
                        buff_A, ldim_A,
                        buff_d, buff_e,
                        buff_t,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( ssytd2_check( uplo, m,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_t,
                                           info ) )
    }
    {
        LAPACK_hetrd_body(s)
    }
}
LAPACK_hetd2(d,sy)
{
    {
        if ( *uplo == 'U' )
        {
            dsytd2_fla( uplo, m,
                        buff_A, ldim_A,
                        buff_d, buff_e,
                        buff_t,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( dsytd2_check( uplo, m,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_t,
                                           info ) )
    }
    {
        LAPACK_hetrd_body(d)
    }
}

#ifdef FLA_LAPACK2FLAME_SUPPORT_COMPLEX
LAPACK_hetd2(c,he)
{
    {
        if ( *uplo == 'U' )
        {
            chetd2_fla( uplo, m,
                        (complex*)buff_A, ldim_A,
                        (real*)buff_d, (real*)buff_e,
                        (complex*)buff_t,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( chetd2_check( uplo, m,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_t,
                                           info ) )
    }
    {
        LAPACK_hetrd_body(c)
    }
}
LAPACK_hetd2(z,he)
{
    {
        if ( *uplo == 'U' )
        {
            zhetd2_fla( uplo, m,
                        (doublecomplex*)buff_A, ldim_A,
                        (doublereal*)buff_d, (doublereal*)buff_e,
                        (doublecomplex*)buff_t,
                        info );
            return 0;
        }
    }
    {
        LAPACK_RETURN_CHECK( zhetd2_check( uplo, m,
                                           buff_A, ldim_A,
                                           buff_d, buff_e,
                                           buff_t,
                                           info ) )
    }
    {
        LAPACK_hetrd_body(z)
    }
}
#endif


#endif
