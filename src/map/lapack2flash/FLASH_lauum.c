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
  LAUUM computes the product U * U**H or L**H * L, where the triangular
  factor U or L is stored in the upper or lower triangular part of
  the array A.

  INFO is INTEGER
  = 0: successful exit
  < 0: if INFO = -k, the k-th argument had an illegal value - LAPACK_lauum_op_check
*/

#define LAPACK_lauum(prefix)                                                        \
  int F77_ ## prefix ## lauum( char* uplo,                                          \
                               int* n,                                              \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A,  \
                               int* info )

#define LAPACK_lauum_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Uplo     uplo_fla;                                        \
  FLA_Obj      A;                                               \
  FLA_Error    init_result;                                     \
  dim_t        blocksize = min( FLASH_get_preferred_blocksize(),\
                                *ldim_A );                      \
                                                                \
  FLA_Init_safe( &init_result );                                \
                                                                \
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
  FLASH_Ttmm( uplo_fla, A );                                    \
                                                                \
  FLASH_Obj_free_without_buffer( &A );                          \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \
  *info = 0;                                                    \
                                                                \
  return 0;

LAPACK_lauum(s)
{
    {
        LAPACK_RETURN_CHECK( slauum_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_lauum_body(s)
    }
}
LAPACK_lauum(d)
{
    {
        LAPACK_RETURN_CHECK( dlauum_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_lauum_body(d)
    }
}
LAPACK_lauum(c)
{
    {
        LAPACK_RETURN_CHECK( clauum_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_lauum_body(c)
    }
}
LAPACK_lauum(z)
{
    {
        LAPACK_RETURN_CHECK( zlauum_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_lauum_body(z)
    }
}

#define LAPACK_lauu2(prefix)                                            \
  int F77_ ## prefix ## lauu2( char* uplo,                              \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               int* info )


LAPACK_lauu2(s)
{
    {
        LAPACK_RETURN_CHECK( slauu2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_lauum_body(s)
    }
}
LAPACK_lauu2(d)
{
    {
        LAPACK_RETURN_CHECK( dlauu2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_lauum_body(d)
    }
}
LAPACK_lauu2(c)
{
    {
        LAPACK_RETURN_CHECK( clauu2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_lauum_body(c)
    }
}
LAPACK_lauu2(z)
{
    {
        LAPACK_RETURN_CHECK( zlauu2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_lauum_body(z)
    }
}


#endif
