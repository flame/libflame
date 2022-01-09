/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
    Copyright (c) 2021 Advanced Micro Devices, Inc.  All rights reserved.
    Mar 16, 2021
*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_util_defs.h"
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

extern void DTL_Trace(
		    uint8 ui8LogLevel,
		    uint8 ui8LogType,
		    const int8 *pi8FileName,
		    const int8 *pi8FunctionName,
		    uint32 ui32LineNumber,
		    const int8 *pi8Message);

#define FLA_ENABLE_ALT_PATH 0

#define LAPACK_getrf(prefix)                                                           \
  int F77_ ## prefix ## getrf( integer* m,                                             \
                               integer* n,                                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               integer* buff_p,                                        \
                               integer* info )

#ifndef FLA_ENABLE_MULTITHREADING

#if FLA_AMD_OPT /* FLA_AMD_OPT */

/* FLA_AMD_OPT enables the code which selects algorithm variants based on size */
#define LAPACK_getrf_body_d(prefix)                                                    \
  AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);                                        \
  if( *m <= FLA_DGETRF_SMALL_THRESH0 && *n <= FLA_DGETRF_SMALL_THRESH0 )               \
  {                                                                                    \
    FLA_LU_piv_small_d_var0( m, n, buff_A, ldim_A, buff_p, info );                     \
  }                                                                                    \
  else if( *m < FLA_DGETRF_SMALL_THRESH1 && *n < FLA_DGETRF_SMALL_THRESH1 )            \
  {                                                                                    \
    FLA_LU_piv_small_d_var1( m, n, buff_A, ldim_A, buff_p, info );                     \
  }                                                                                    \
  else                                                                                 \
  {                                                                                    \
    dgetrf2_( m, n, buff_A, ldim_A, buff_p, info);                                     \
  }                                                                                    \
  AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);                                         \
  return 0;

#else /* FLA_AMD_OPT */

/* Original FLA path */
#define LAPACK_getrf_body_d(prefix)                                                    \
  AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);                                        \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                               \
  FLA_Obj      A, p;                                                                   \
  integer      min_m_n    = min( *m, *n );                                             \
  FLA_Error    e_val = FLA_SUCCESS;                                                    \
  FLA_Error    init_result;                                                            \
  FLA_Bool skip = FALSE;                                                               \
                                                                                       \
  FLA_Init_safe( &init_result );                                                       \
                                                                                       \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                               \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                                     \
                                                                                       \
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &p );                            \
  FLA_Obj_attach_buffer( buff_p, 1, min_m_n, &p );                                     \
                                                                                       \
  e_val = FLA_LU_piv( A, p );                                                          \
  FLA_Shift_pivots_to( FLA_LAPACK_PIVOTS, p );                                         \
                                                                                       \
  FLA_Obj_free_without_buffer( &A );                                                   \
  FLA_Obj_free_without_buffer( &p );                                                   \
                                                                                       \
  FLA_Finalize_safe( init_result );                                                    \
                                                                                       \
  if ( e_val != FLA_SUCCESS ) *info = e_val + 1;                                       \
  AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);                                         \
  return 0;

#endif /* FLA_AMD_OPT */

// Note that p should be set zero.
#define LAPACK_getrf_body(prefix)                               \
  AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);                 \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Obj      A, p;                                            \
  integer      min_m_n    = min( *m, *n );                      \
  FLA_Error    e_val = FLA_SUCCESS;                             \
  FLA_Error    init_result;                                     \
  FLA_Bool skip = FALSE;                                        \
                                                                \
                                                                \
  if( *m < FLA_GETRF_SMALL  && *n < FLA_GETRF_SMALL && !FLA_ENABLE_ALT_PATH )  /* Small sizes- lapack path */       \
  {                                                                                                                 \
    switch(datatype)                                                                                                \
    {                                                                                                               \
       case FLA_FLOAT:                                                                                              \
       { lapack_sgetrf( m, n, buff_A, ldim_A, buff_p, info); break; }                                               \
       case FLA_COMPLEX:                                                                                            \
       { lapack_cgetrf( m, n, buff_A, ldim_A, buff_p, info); break; }                                               \
       case FLA_DOUBLE_COMPLEX:                                                                                     \
       { lapack_zgetrf( m, n, buff_A, ldim_A, buff_p, info); break; }                                               \
    }  if ( *info != 0 ) skip  = TRUE;                                                                              \
                                                                                                                    \
  }                                                                                                                 \
  else if( ( datatype == FLA_FLOAT && *m < FLA_GETRF_FLOAT && * n < FLA_GETRF_FLOAT  )||                            \
           ( datatype == FLA_COMPLEX && *m < FLA_GETRF_COMPLEX  && *n < FLA_GETRF_COMPLEX  ) ||                     \
           ( datatype == FLA_DOUBLE_COMPLEX && *m < FLA_GETRF_DOUBLE_COMPLEX  && *n < FLA_GETRF_DOUBLE_COMPLEX  ) ) \
  {                                                                                    \
    FLA_Init_safe( &init_result );                                                     \
                                                                                       \
    FLA_Obj_create_without_buffer( datatype, *m, *n, &A );                             \
    FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                                   \
                                                                                       \
    FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &p );                          \
    FLA_Obj_attach_buffer( buff_p, 1, min_m_n, &p );                                   \
                                                                                       \
    e_val = FLA_LU_piv( A, p );                                                        \
    FLA_Shift_pivots_to( FLA_LAPACK_PIVOTS, p );                                       \
                                                                                       \
    FLA_Obj_free_without_buffer( &A );                                                 \
    FLA_Obj_free_without_buffer( &p );                                                 \
                                                                                       \
    FLA_Finalize_safe( init_result );                                                  \
 }                                                                                     \
  else                                                                                 \
  {                                                                                    \
    switch(datatype)                                                                   \
    {                                                                                  \
       case FLA_FLOAT:                                                                 \
       { sgetrf2_( m, n, buff_A, ldim_A, buff_p, info); break; }                       \
       case FLA_COMPLEX:                                                               \
       { cgetrf2_( m, n, buff_A, ldim_A, buff_p, info); break; }                       \
       case FLA_DOUBLE_COMPLEX:                                                        \
       { zgetrf2_( m, n, buff_A, ldim_A, buff_p, info); break; }                       \
    }  if ( *info != 0 ) skip  = TRUE;                                                 \
                                                                                       \
  }                                                                                    \
                                                                                       \
  if ( e_val != FLA_SUCCESS ) *info = e_val + 1;                                       \
  else if( skip != TRUE )       *info = 0;                                             \
  AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);                                         \
  return 0;

#else /* FLA_ENABLE_MULTITHREADING */

#define LAPACK_getrf_body_d LAPACK_getrf_body

// Note that p should be set zero.
#define LAPACK_getrf_body(prefix)                               \
  AOCL_DTL_TRACE_ENTRY(AOCL_DTL_LEVEL_TRACE_5);                 \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Obj      A, p, AH, ph;                                    \
  integer      min_m_n    = min( *m, *n );                      \
  dim_t        nth, b_flash;                                    \
  FLA_Error    e_val;                                           \
  FLA_Error    init_result;                                     \
                                                                \
  nth = FLASH_get_num_threads( 1 );                             \
  b_flash = FLA_EXT_HIER_BLOCKSIZE;                             \
                                                                \
  FLA_Init_safe( &init_result );                                \
  FLASH_Queue_set_num_threads( nth );                           \
                                                                \
  FLA_Obj_create_without_buffer( datatype, *m, *n, &A );        \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );              \
                                                                \
  FLA_Obj_create_without_buffer( FLA_INT, min_m_n, 1, &p );     \
  FLA_Obj_attach_buffer( buff_p, 1, min_m_n, &p );              \
                                                                \
  FLA_Set( FLA_ZERO, p );                                       \
                                                                \
  FLASH_Obj_create_hier_copy_of_flat( A, 1, &b_flash, &AH );    \
  FLASH_Obj_create_hier_copy_of_flat( p, 1, &b_flash, &ph );    \
                                                                \
  e_val = FLASH_LU_piv( AH, ph );                               \
                                                                \
  FLASH_Obj_flatten( AH, A );                                   \
  FLASH_Obj_flatten( ph, p );                                   \
  FLA_Shift_pivots_to( FLA_LAPACK_PIVOTS, p );                  \
                                                                \
  FLA_Obj_free( &AH );                                          \
  FLA_Obj_free( &ph );                                          \
                                                                \
  FLA_Obj_free_without_buffer( &A );                            \
  FLA_Obj_free_without_buffer( &p );                            \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \
  if ( e_val != FLA_SUCCESS ) *info = e_val + 1;                \
  else                        *info = 0;                        \
                                                                \
  AOCL_DTL_TRACE_EXIT(AOCL_DTL_LEVEL_TRACE_5);                  \
  return 0;

#endif /* FLA_ENABLE_MULTITHREADING */

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
        LAPACK_getrf_body_d(d)
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
  int F77_ ## prefix ## getf2( integer* m,                                  \
                               integer* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               integer* buff_p,                             \
                               integer* info )

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
        LAPACK_getrf_body_d(d)
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
