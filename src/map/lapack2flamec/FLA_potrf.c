/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
/*
    Copyright (c) 2021 Advanced Micro Devices, Inc.  All rights reserved.
    Aug 5, 2021
*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLAME

#include "FLA_lapack2flame_util_defs.h"
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

extern void DTL_Trace(
		    uint8 ui8LogLevel,
		    uint8 ui8LogType,
		    const int8 *pi8FileName,
		    const int8 *pi8FunctionName,
		    uint32 ui32LineNumber,
		    const int8 *pi8Message);


#define LAPACK_potrf(prefix)                                          \
  int F77_ ## prefix ## potrf( char* uplo,                            \
                               integer*  n,                           \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                               integer*  ldim_A,                      \
                               integer*  info )

#if FLA_AMD_OPT
#define LAPACK_potrf_body_s(prefix)                                                          \
  if( *n < FLA_POTRF_FLOAT_SMALL )                                                           \
  lapack_spotf2( uplo, n, buff_A, ldim_A,  info );                                           \
  else                                                                                       \
  lapack_spotrf( uplo, n, buff_A, ldim_A,  info );                                           \

#define LAPACK_potrf_body_d(prefix)                                                          \
  if( *n < FLA_POTRF_DOUBLE_SMALL )                                                          \
  lapack_dpotf2( uplo, n, buff_A, ldim_A,  info );                                           \
  else                                                                                       \
  lapack_dpotrf( uplo, n, buff_A, ldim_A,  info );                                           \

#endif

#define LAPACK_potrf_body(prefix)                                                            \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);                                     \
  FLA_Uplo     uplo_fla;                                                                     \
  FLA_Obj      A;                                                                            \
  FLA_Error    e_val = FLA_SUCCESS;                                                          \
  FLA_Error    init_result;                                                                  \
  FLA_Init_safe( &init_result );                                                             \
  FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );                                     \
                                                                                             \
  FLA_Obj_create_without_buffer( datatype, *n, *n, &A );                                     \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );                                           \
                                                                                             \
  e_val = FLA_Chol( uplo_fla, A );                                                           \
                                                                                             \
  FLA_Obj_free_without_buffer( &A );                                                         \
                                                                                             \
  FLA_Finalize_safe( init_result );                                                          \
                                                                                             \
  if ( e_val != FLA_SUCCESS ) *info = e_val + 1;                                             \
  else                        *info = 0;                                                     \
                                                                                             \


LAPACK_potrf(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spotrf inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);

    {
        LAPACK_RETURN_CHECK_VAR1( spotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ), fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
#if FLA_AMD_OPT
        {   
            LAPACK_potrf_body_s(s);
        }
#else
        {
            LAPACK_potrf_body(s)
        }
#endif
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}

LAPACK_potrf(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dpotrf inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( dpotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
#if FLA_AMD_OPT
        {
        LAPACK_potrf_body_d(d)
        }
#else
        {
            LAPACK_potrf_body(d)
        }
#endif
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_potrf(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cpotrf inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( cpotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potrf_body(c)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_potrf(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zpotrf inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( zpotrf_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ), fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potrf_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}


#define LAPACK_potf2(prefix)                                    \
  int F77_ ## prefix ## potf2( char* uplo,                      \
                               integer*  n,                         \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                               integer*  ldim_A,                    \
                               integer*  info )

LAPACK_potf2(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spotf2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( spotf2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
#if FLA_AMD_OPT
        {
            LAPACK_potrf_body_s(s)
        }
#else
        {
            LAPACK_potrf_body(s)
        }
#endif
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_potf2(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dpotf2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dpotf2_check(uplo, n,
                                              buff_A, ldim_A,
                                              info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
#if FLA_AMD_OPT
        {
            LAPACK_potrf_body_d(d)
        }
#else
        {
            LAPACK_potrf_body(d)
        }
#endif
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_potf2(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cpotf2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( cpotf2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ), fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potrf_body(c)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_potf2(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zpotf2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zpotf2_check(uplo, n,
                                              buff_A, ldim_A,
                                              info), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potrf_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
#endif
