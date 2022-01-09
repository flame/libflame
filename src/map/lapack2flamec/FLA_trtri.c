/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME


#include "FLA_lapack2flame_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
   TRTRI computes the inverse of a upper or lower triangular
   matrix A.

   This is the Level 3 BLAS version of the algorithm.

   INFO
   = 0: successful exit
   < 0: if INFO = -i, the i-th argument had an illegal value -> LAPACK_trtri_op_check
   > 0: if INFO = i, A(i,i) is exactly zero. The triangular -> LAPACK_check_diag_zero
   matrix is singular and its inverse can not be computed.

   Hence, if the routine passes above error checking, FLA_Trinv should not produce any error.
*/

#define LAPACK_trtri(prefix)                                    \
  int F77_ ## prefix ## trtri( char* uplo,                      \
                               char* diag,                      \
                               integer* n,                          \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                               integer* ldim_A,                     \
                               integer* info )

#define LAPACK_trtri_body(prefix)                               \
  AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);                 \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Uplo     uplo_fla;                                        \
  FLA_Diag     diag_fla;                                        \
  FLA_Obj      A;                                               \
  FLA_Error    init_result;                                     \
                                                                \
  FLA_Init_safe( &init_result );                                        \
                                                                        \
  FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );                \
  FLA_Param_map_netlib_to_flame_diag( diag, &diag_fla );                \
                                                                        \
  FLA_Obj_create_without_buffer( datatype, *n, *n, &A );                \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                      \
                                                                        \
  FLA_Trinv( uplo_fla, diag_fla, A );                                   \
                                                                        \
  FLA_Obj_free_without_buffer( &A );                                    \
                                                                        \
  FLA_Finalize_safe( init_result );                                     \
                                                                        \
  *info = 0;                                                            \
  AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);                  	      \
                                                                        \
  return 0;

LAPACK_trtri(s)
{
    {
        LAPACK_RETURN_CHECK( strtri_check( uplo, diag, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_trtri_body(s)
    }
}
LAPACK_trtri(d)
{
    {
        LAPACK_RETURN_CHECK( dtrtri_check( uplo, diag, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_trtri_body(d)
    }
}
LAPACK_trtri(c)
{
    {
        LAPACK_RETURN_CHECK( ctrtri_check( uplo, diag, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_trtri_body(c)
    }
}
LAPACK_trtri(z)
{
    {
        LAPACK_RETURN_CHECK( ztrtri_check( uplo, diag, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_trtri_body(z)
    }
}


#define LAPACK_trti2(prefix)                                    \
  int F77_ ## prefix ## trti2( char* uplo,                      \
                               char* diag,                      \
                               integer* n,                          \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                               integer* ldim_A,                     \
                               integer* info )

LAPACK_trti2(s)
{
    {
        LAPACK_RETURN_CHECK( strti2_check( uplo, diag, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_trtri_body(s)
    }
}
LAPACK_trti2(d)
{
    {
        LAPACK_RETURN_CHECK( dtrti2_check( uplo, diag, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_trtri_body(d)
    }
}
LAPACK_trti2(c)
{
    {
        LAPACK_RETURN_CHECK( ctrti2_check( uplo, diag, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_trtri_body(c)
    }
}
LAPACK_trti2(z)
{
    {
        LAPACK_RETURN_CHECK( ztrti2_check( uplo, diag, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_trtri_body(z)
    }
}


#endif
