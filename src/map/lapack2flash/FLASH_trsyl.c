/*

    Copyright (C) 2014, The University of Texas at Austin
    Copyright (C) 2022, Advanced Micro Devices, Inc.

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_YET_LAPACK2FLAME

#include "FLASH_lapack2flash_util_defs.h"
#include "check/FLA_trsyl.in"

/*
   TRSYL solves the Sylvester matrix equation:

   op(A)*X + X*op(B) = scale*C or
   op(A)*X - X*op(B) = scale*C,

   where op(A) = A or A**H, and A and B are both upper triangular. A is
   M-by-M and B is N-by-N; the right hand side C and the solution X are
   M-by-N; and scale is an output scale factor, set <= 1 to avoid
   overflow in X.

   INFO is INTEGER
   = 0: successful exit
   < 0: if INFO = -i, the i-th argument had an illegal value - LAPACK_trsyl_op_check
   = 1: A and B have common or very close eigenvalues; perturbed
   values were used to solve the equation (but the matrices
   A and B are unchanged).

   FLAME does not check what LAPACK does. So, let's exclude at this moment.
*/

#define LAPACK_trsyl(prefix)                                                       \
  int F77_ ## prefix ## trsyl( char* transa,                                       \
                               char* transb,                                       \
                               int* sgn,                                           \
                               int* m,                                             \
                               int* n,                                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, int* ldim_B, \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_C, int* ldim_C, \
                               PREFIX2LAPACK_REALDEF(prefix)* scale,               \
                               int* info )

#define LAPACK_trsyl_body(prefix, srname)                       \
  FLA_Datatype datatype       =  PREFIX2FLAME_DATATYPE(prefix); \
  FLA_Datatype datatype_scale =  PREFIX2FLAME_REALTYPE(prefix); \
  FLA_Trans    transa_fla;                                      \
  FLA_Trans    transb_fla;                                      \
  FLA_Obj      sgn_fla;                                         \
  FLA_Obj      A, B, C;                                         \
  FLA_Obj      scale_fla;                                       \
  FLA_Error    e_val;                                           \
  FLA_Error    init_result;                                     \
  dim_t        blocksize;                                       \
                                                                \
  LAPACK_trsyl_op_check(prefix,srname)                          \
                                                                \
  FLA_Init_safe( &init_result );                                \
                                                                \
  FLA_Param_map_netlib_to_flame_trans( transa, &transa_fla );   \
  FLA_Param_map_netlib_to_flame_trans( transb, &transb_fla );   \
                                                                \
  if ( *sgn == 1 ) sgn_fla = FLA_ONE;                           \
  else             sgn_fla = FLA_MINUS_ONE;                     \
                                                                \
  blocksize = min( FLASH_get_preferred_blocksize(), *ldim_A );  \
                                                                \
  FLA_Bool toggle = FLASH_Check_offload_to_gpu( blocksize, *m,  \
                          *n, FLASH_get_tile_offload() );       \
                                                                \
  FLASH_Obj_create_without_buffer( datatype,                    \
                                   *m,                          \
                                   *m,                          \
                                   FLASH_get_depth(),           \
                                   &blocksize,                  \
                                   &A );                        \
  FLASH_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );            \
                                                                \
  blocksize = min( FLASH_get_preferred_blocksize(), *ldim_B );  \
  FLASH_Obj_create_without_buffer( datatype,                    \
                                   *n,                          \
                                   *n,                          \
                                   FLASH_get_depth(),           \
                                   &blocksize,                  \
                                   &B );                        \
  FLASH_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );            \
                                                                \
  blocksize = min( FLASH_get_preferred_blocksize(), *ldim_C );  \
  FLASH_Obj_create_without_buffer( datatype,                    \
                                   *m,                          \
                                   *n,                          \
                                   FLASH_get_depth(),           \
                                   &blocksize,                  \
                                   &C );                        \
  FLASH_Obj_attach_buffer( buff_C, 1, *ldim_C, &C );            \
                                                                \
  blocksize = min( FLASH_get_preferred_blocksize(), 1 );        \
  FLASH_Obj_create_without_buffer( datatype_scale,              \
                                   1,                           \
                                   1,                           \
                                   FLASH_get_depth(),           \
                                   &blocksize,                  \
                                   &scale_fla );                \
  FLASH_Obj_attach_buffer( &scale, 1, 1, &scale_fla );          \
                                                                \
  e_val = FLASH_Sylv( transa_fla,                               \
                      transb_fla,                               \
                      sgn_fla,                                  \
                      A,                                        \
                      B,                                        \
                      C,                                        \
                      scale_fla );                              \
                                                                \
  FLASH_Obj_free_without_buffer( &A );                          \
  FLASH_Obj_free_without_buffer( &B );                          \
  FLASH_Obj_free_without_buffer( &C );                          \
  FLASH_Obj_free_without_buffer( &scale_fla );                  \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \
  FLASH_Toggle_gpu_offload( toggle );                           \
                                                                \
  if ( e_val != FLA_SUCCESS ) *info = 1;                        \
  else                        *info = 0;                        \
                                                                \
  return 0;


LAPACK_trsyl(s)
{
    LAPACK_trsyl_body(s, STRSYL)
}
LAPACK_trsyl(d)
{
    LAPACK_trsyl_body(d, DTRSYL)
}
LAPACK_trsyl(c)
{
    LAPACK_trsyl_body(c, CTRSYL)
}
LAPACK_trsyl(z)
{
    LAPACK_trsyl_body(z, ZTRSYL)
}

#endif
