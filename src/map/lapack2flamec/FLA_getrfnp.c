/*

    Copyright (c) 2020 Advanced Micro Devices, Inc.

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
  GETRFNP computes an LU factorization of a general M-by-N matrix A
  without pivoting.

  INFO
  = 0: successful exit
  < 0: if INFO = -i, the i-th argument had an illegal value - LAPACK_getrfnp_op_check
  > 0: if INFO = i, U(i,i) is exactly zero. The factorization - FLA_LU_nopiv
  has been completed, but the factor U is exactly
  singular, and division by zero will occur if it is used
  to solve a system of equations.
*/


extern void DTL_Trace(
		    uint8 ui8LogLevel,
		    uint8 ui8LogType,
		    const int8 *pi8FileName,
		    const int8 *pi8FunctionName,
		    uint32 ui32LineNumber,
		    const int8 *pi8Message);

#define LAPACK_getrfnp(prefix)                                            \
  int F77_ ## prefix ## getrfnp( int* m,                                  \
                               int* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               int* info )

#define LAPACK_getrfnp_body(prefix)                               \
  AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5); 		\
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Obj      A;                                            \
  FLA_Error    e_val;                                           \
  FLA_Error    init_result;                                     \
  extern TLS_CLASS_SPEC fla_lu_t*  fla_lu_nopiv_cntl2; 		\
  								\
  FLA_Init_safe( &init_result );                                \
                                                                \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );        \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );              \
  								\
  e_val = FLA_LU_nopiv_internal( A,  fla_lu_nopiv_cntl2);       \
                                                                \
  FLA_Obj_free_without_buffer( &A );                            \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \
  if ( e_val != FLA_SUCCESS ) *info = e_val + 1;                \
  else                        *info = 0;                        \
                                                                \
  AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);                  \
  return 0;

LAPACK_getrfnp(s)
{
    {
        LAPACK_RETURN_CHECK( sgetrfnp_check( m, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_getrfnp_body(s)
    }
}
LAPACK_getrfnp(d)
{
    {
        LAPACK_RETURN_CHECK( dgetrfnp_check( m, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_getrfnp_body(d)
    }
}
LAPACK_getrfnp(c)
{
    {
        LAPACK_RETURN_CHECK( cgetrfnp_check( m, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_getrfnp_body(c)
    }
}
LAPACK_getrfnp(z)
{
    {
        LAPACK_RETURN_CHECK( zgetrfnp_check( m, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_getrfnp_body(z)
    }
}


#endif
