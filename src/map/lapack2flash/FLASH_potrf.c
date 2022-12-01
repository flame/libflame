/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022, Advanced Micro Devices, Inc.

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
  POTRF computes the Cholesky factorization of a symmetric (hermitian)
  positive definite matrix A.

  INFO is INTEGER
  = 0: successful exit
  < 0: if INFO = -i, the i-th argument had an illegal value -  LAPACK_potrf_op_check
  > 0: if INFO = i, the leading minor of order i is not - FLA_Chol
  positive definite, and the factorization could not be
  completed.
*/

#define LAPACK_potrf(prefix)                                            \
  int F77_ ## prefix ## potrf( char* uplo,                              \
                               int*  n,                                 \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A,   \
                               int*  ldim_A,                            \
                               int*  info )

#define LAPACK_potrf_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Uplo     uplo_fla;                                        \
  FLA_Obj      A;                                               \
  FLA_Error    e_val;                                           \
  FLA_Error    init_result;                                     \
  dim_t        blocksize = min( FLASH_get_preferred_blocksize(),\
                                *ldim_A );                      \
                                                                \
  FLA_Init_safe( &init_result );                                \
  FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );        \
                                                                \
  FLASH_Obj_create_without_buffer( datatype,                    \
                                   *n,                          \
                                   *n,                          \
                                   FLASH_get_depth(),           \
                                   &blocksize,                  \
                                   &A );                        \
  FLASH_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );            \
                                                                \
  e_val = FLASH_Chol( uplo_fla, A );                            \
                                                                \
  FLASH_Obj_free_without_buffer( &A );                          \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \
  if ( e_val != FLA_SUCCESS ) *info = e_val + 1;                \
  else                        *info = 0;                        \
                                                                \
  return 0;

LAPACK_potrf(s)
{
    {
        LAPACK_RETURN_CHECK( spotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(s)
    }
}
LAPACK_potrf(d)
{
    {
        LAPACK_RETURN_CHECK( dpotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(d)
    }
}
LAPACK_potrf(c)
{
    {
        LAPACK_RETURN_CHECK( cpotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(c)
    }
}
LAPACK_potrf(z)
{
    {
        LAPACK_RETURN_CHECK( zpotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(z)
    }
}

#define LAPACK_potf2(prefix)                                    \
  int F77_ ## prefix ## potf2( char* uplo,                      \
                               int*  n,                         \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                               int*  ldim_A,                    \
                               int*  info )

LAPACK_potf2(s)
{
    {
        LAPACK_RETURN_CHECK( spotf2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(s)
    }
}
LAPACK_potf2(d)
{
    {
        LAPACK_RETURN_CHECK( dpotf2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(d)
    }
}
LAPACK_potf2(c)
{
    {
        LAPACK_RETURN_CHECK( cpotf2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(c)
    }
}
LAPACK_potf2(z)
{
    {
        LAPACK_RETURN_CHECK( zpotf2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potrf_body(z)
    }
}

#endif
