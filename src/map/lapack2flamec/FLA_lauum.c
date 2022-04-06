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
  LAUUM computes the product U * U**H or L**H * L, where the triangular
  factor U or L is stored in the upper or lower triangular part of
  the array A.

  INFO is INTEGER
  = 0: successful exit
  < 0: if INFO = -k, the k-th argument had an illegal value - LAPACK_lauum_op_check
*/

#define LAPACK_lauum(prefix)                                            \
  int F77_ ## prefix ## lauum( char* uplo,                              \
                               integer* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               integer* info )

#define LAPACK_lauum_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Uplo     uplo_fla;                                        \
  FLA_Obj      A;                                               \
  FLA_Error    init_result;                                     \
                                                                \
  FLA_Init_safe( &init_result );                                \
                                                                \
  FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );        \
                                                                \
  FLA_Obj_create_without_buffer( datatype, *n, *n, &A );        \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );              \
                                                                \
  FLA_Ttmm( uplo_fla, A );                                      \
                                                                \
  FLA_Obj_free_without_buffer( &A );                            \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \
  *info = 0;                                                    \
                                                                \


LAPACK_lauum(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slauum inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);

    {
        LAPACK_RETURN_CHECK_VAR1( slauum_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ),fla_error )
    }if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_lauum_body(s)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_lauum(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlauum inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);

    {
        LAPACK_RETURN_CHECK_VAR1(dlauum_check(uplo, n,
                                              buff_A, ldim_A,
                                              info),
                                 fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_lauum_body(d)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_lauum(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("clauum inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);

    {
        LAPACK_RETURN_CHECK_VAR1(clauum_check(uplo, n,
                                              buff_A, ldim_A,
                                              info),
                                 fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_lauum_body(c)
      /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_lauum(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlauum inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);

    {
        LAPACK_RETURN_CHECK_VAR1(zlauum_check(uplo, n,
                                              buff_A, ldim_A,
                                              info),
                                 fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_lauum_body(z)
    /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}

#define LAPACK_lauu2(prefix)                                            \
  int F77_ ## prefix ## lauu2( char* uplo,                              \
                               integer* n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                               integer* info )


LAPACK_lauu2(s)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("slauu2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1( slauu2_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_lauum_body(s)
    /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_lauu2(d)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("dlauu2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(dlauu2_check(uplo, n,
                                              buff_A, ldim_A,
                                              info),
                                 fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_lauum_body(d)
   /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_lauu2(c)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("clauu2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(clauu2_check(uplo, n,
                                              buff_A, ldim_A,
                                              info),
                                 fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_lauum_body(c)
    /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_lauu2(z)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zlauu2 inputs: uplo %c, n %" FLA_IS ", lda %" FLA_IS "", *uplo, *n, *ldim_A);
    {
        LAPACK_RETURN_CHECK_VAR1(zlauu2_check(uplo, n,
                                              buff_A, ldim_A,
                                              info),
                                 fla_error)
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_lauum_body(z)
    /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}


#endif
