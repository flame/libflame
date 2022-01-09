/*
 *  Copyright (c) 2020-21 Advanced Micro Devices, Inc. All rights reserved.
 * */

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  GETRFNPI computes an LU factorization(complete/incomplete) based on given parameter
  nfact , of a general M-by-N matrix A without pivoting.

  INFO
  = 0: successful exit
  < 0: if INFO = -i, the i-th argument had an illegal value - LAPACK_getrfnp_op_check
  > 0: if INFO = i, U(i,i) is exactly zero. The factorization - FLA_LU_nopiv
  has been completed, but the factor U is exactly
  singular, and division by zero will occur if it is used
  to solve a system of equations.

  m
  INTEGER. The number of rows in matrix A; m≥ 0.
  n
  INTEGER. The number of columns in matrix A; n≥ 0.
  nfact
  INTEGER. The number of rows and columns to factor; 0 ≤nfact≤ min(m, n). Note that if nfact < min(m, n), incomplete factorization is performed.

   A
   REAL for sgetrfnpi
   DOUBLE PRECISION for dgetrfnpi
   COMPLEX for cgetrfnpi
   DOUBLE COMPLEX for zgetrfnpi
   Array of size (lda,*). Contains the matrix A. The second dimension of a must be at least max(1, n).
   lda
   INTEGER. The leading dimension of array a. lda≥ max(1, m).

*/

extern void DTL_Trace(
		uint8 ui8LogLevel,
		uint8 ui8LogType,
		const int8 *pi8FileName,
		const int8 *pi8FunctionName,
		uint32 ui32LineNumber,
		const int8 *pi8Message);

#define LAPACK_getrfnpi(prefix)                                                       \
  int F77_ ## prefix ## getrfnpi( integer* m,                                             \
                                  integer* n,                                             \
                                  integer* nfact,                                         \
                                  PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                                  integer* info )

#define LAPACK_getrfnpi_body(prefix)                                                           \
  AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);                                                \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                                       \
  FLA_Error e_val ;                                                                            \
                                                                                               \
  if ( (datatype == FLA_DOUBLE) && ((*m) * (*n) <=  FLA_MN_SIZE) && ((*nfact) <= (FLA_NFACT_PERCENT * (*m))) && ((*m) == (*n)) )                 \
  {                                                                                                                                              \
    if( (*n)*(*nfact)-(*nfact-1)/4 <= FLA_FULL_DGER_CONSTANT  )                                                                                  \
    e_val = FLA_LU_nopiv_id_unblk_var2( *m, *n, buff_A, *nfact, 1, *ldim_A);                                                                     \
    else                                                                                                                                         \
    e_val =  FLA_LU_nopiv_id_unblk_var1( *m, *n, buff_A,  *nfact, 1, *ldim_A);                                                                   \
  }                                                                                                                                              \
  else                                                                                                                                           \
  {                                                                                                                                              \
     FLA_Obj   A;                                                                                                                                \
     FLA_Error init_result;                                                                                                                      \
     FLA_Init_safe( &init_result );                                                                                                              \
     FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                                                                                      \
     FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                                                                                            \
     switch( datatype )                                                                                                                          \
     {                                                                                                                                           \
        case FLA_FLOAT:                                                                                                                          \
        { e_val = FLA_LU_nopiv_is_blk_var1( *m, *n, A, buff_A, *nfact, 1, *ldim_A); break; }                                                     \
        case FLA_DOUBLE:                                                                                                                         \
        { e_val = FLA_LU_nopiv_id_blk_var1( *m, *n, A, buff_A, *nfact, 1, *ldim_A); break; }                                                     \
        case FLA_COMPLEX:                                                                                                                        \
        { e_val = FLA_LU_nopiv_ic_blk_var1( *m, *n, A, buff_A, *nfact, 1, *ldim_A); break; }                                                     \
        case FLA_DOUBLE_COMPLEX:                                                                                                                 \
        { e_val = FLA_LU_nopiv_iz_blk_var1( *m, *n, A, buff_A, *nfact, 1, *ldim_A); break; }                                                     \
     }                                                                                                                                           \
     FLA_Obj_free_without_buffer( &A );                                                                                                          \
     FLA_Finalize_safe( init_result );                                                                                                           \
  }                                                                                                                                              \
  if ( e_val != FLA_SUCCESS ) *info = e_val + 1;                                                                                                 \
  else                        *info = 0;                                                                                                         \
  AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);                                                                                                   \
                                                                                                                                                 \
  return 0;




LAPACK_getrfnpi(s)
{
    {
        LAPACK_RETURN_CHECK( sgetrfnpi_check( m, n, nfact,
                                              buff_A, ldim_A,
                                              info ) )
    }
    {
        LAPACK_getrfnpi_body(s)
    }
}
LAPACK_getrfnpi(d)
{
    {
        LAPACK_RETURN_CHECK( dgetrfnpi_check( m, n, nfact,
                                              buff_A, ldim_A,
                                              info ) )
    }
    {
        LAPACK_getrfnpi_body(d)
    }
}
LAPACK_getrfnpi(c)
{
    {
        LAPACK_RETURN_CHECK( cgetrfnpi_check( m, n, nfact,
                                              buff_A, ldim_A,
                                              info ) )
    }
    {
        LAPACK_getrfnpi_body(c)
    }
}
LAPACK_getrfnpi(z)
{
    {
        LAPACK_RETURN_CHECK( zgetrfnpi_check( m, n, nfact,
                                              buff_A, ldim_A,
                                              info ) )
    }
    {
        LAPACK_getrfnpi_body(z)
    }
}


#endif
