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
  GETRF computes an LU factorization of a general M-by-N matrix A
  using partial pivoting with row interchanges.

  INFO
  = 0: successful exit
  < 0: if INFO = -i, the i-th argument had an illegal value - LAPACK_getrf_op_check
  > 0: if INFO = i, U(i,i) is exactly zero. The factorization - FLA_LU_piv
  has been completed, but the factor U is exactly
  singular, and division by zero will occur if it is used
  to solve a system of equations.
*/

#define LAPACK_getrf(prefix)                                            \
  int F77_ ## prefix ## getrf( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               int* buff_p,                             \
                               int* info )

// Note that p should be set zero.
#define LAPACK_getrf_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Obj      A, p;                                            \
  int          min_m_n    = min( *m, *n );                      \
  FLA_Error    e_val;                                           \
  FLA_Error    init_result;                                     \
                                                                \
  FLA_Init_safe( &init_result );                                \
                                                                \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );        \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );              \
                                                                \
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &p );     \
  FLA_Obj_attach_buffer( buff_p, 1, min_m_n, &p );              \
  FLA_Set( FLA_ZERO, p );                                       \
                                                                \
  e_val = FLA_LU_piv( A, p );                                   \
  FLA_Shift_pivots_to( FLA_LAPACK_PIVOTS, p );                  \
                                                                \
  FLA_Obj_free_without_buffer( &A );                            \
  FLA_Obj_free_without_buffer( &p );                            \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \
  if ( e_val != FLA_SUCCESS ) *info = e_val + 1;                \
  else                        *info = 0;                        \
                                                                \
  return 0;

LAPACK_getrf(s)
{
    {
        LAPACK_RETURN_CHECK( sgetrf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           info ) )
    }
    {
        LAPACK_getrf_body(s)
    }
}
LAPACK_getrf(d)
{
    {
        LAPACK_RETURN_CHECK( dgetrf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           info ) )
    }
    {
        LAPACK_getrf_body(d)
    }
}
LAPACK_getrf(c)
{
    {
        LAPACK_RETURN_CHECK( cgetrf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           info ) )
    }
    {
        LAPACK_getrf_body(c)
    }

}
LAPACK_getrf(z)
{
    {
        LAPACK_RETURN_CHECK( zgetrf_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           info ) )
    }
    {
        LAPACK_getrf_body(z)
    }
}


#define LAPACK_getf2(prefix)                                            \
  int F77_ ## prefix ## getf2( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               int* buff_p,                             \
                               int* info )

LAPACK_getf2(s)
{
    {
        LAPACK_RETURN_CHECK( sgetf2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           info ) )
    }
    {
        LAPACK_getrf_body(s)
    }
}
LAPACK_getf2(d)
{
    {
        LAPACK_RETURN_CHECK( dgetf2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           info ) )
    }
    {
        LAPACK_getrf_body(d)
    }
}
LAPACK_getf2(c)
{
    {
        LAPACK_RETURN_CHECK( cgetf2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           info ) )
    }
    {
        LAPACK_getrf_body(c)
    }
}
LAPACK_getf2(z)
{
    {
        LAPACK_RETURN_CHECK( zgetf2_check( m, n,
                                           buff_A, ldim_A,
                                           buff_p,
                                           info ) )
    }
    {
        LAPACK_getrf_body(z)
    }
}

#endif
