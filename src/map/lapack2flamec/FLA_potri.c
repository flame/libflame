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
   POTRI computes the inverse of a symmetric (hermitian) positive definite
   matrix A using the Cholesky factorization A = U**H*U or A = L*L**H
   computed by ZPOTRF.
*/

#define LAPACK_potri(prefix)                                    \
  int F77_ ## prefix ## potri( char* uplo,                      \
                               integer*  n,                         \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                               integer*  ldim_A,                         \
                               integer*  info )

#define LAPACK_potri_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Uplo     uplo_fla;                                        \
  FLA_Obj      A;                                               \
  FLA_Error    e_val;                                           \
  FLA_Error    init_result;                                     \
                                                                \
  FLA_Init_safe( &init_result );                                \
  FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );        \
                                                                \
  FLA_Obj_create_without_buffer( datatype, *n, *n, &A );        \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );              \
                                                                \
  e_val = FLA_Trinv( uplo_fla, FLA_NONUNIT_DIAG, A );           \
                                                                \
  if ( e_val != FLA_SUCCESS )            *info = e_val + 1;     \
  else {                                                        \
    e_val = FLA_Ttmm( uplo_fla, A );                            \
    if ( e_val != FLA_SUCCESS )          *info = e_val + 1;     \
  }                                                             \
  FLA_Obj_free_without_buffer( &A );                            \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \


LAPACK_potri(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("spotri inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( spotri_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_potri_body(s)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_potri(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dpotri inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( dpotri_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potri_body(d)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_potri(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("cpotri inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( cpotri_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ), fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potri_body(c)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_potri(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zpotri inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( zpotri_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ), fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_potri_body(z)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}

#endif
