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

extern int zhegst_fla(integer *itype, char *uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, integer *info);
extern int chegst_fla(integer *itype, char *uplo, integer *n, complex * a, integer *lda, complex *b, integer *ldb, integer *info);
extern int chegs2_fla(integer *itype, char *uplo, integer *n, complex * a, integer *lda, complex *b, integer *ldb, integer *info);
extern int zhegs2_fla(integer *itype, char *uplo, integer *n, doublecomplex *a, integer *lda, doublecomplex *b, integer *ldb, integer *info);

/*
  ZHEGST reduces a complex Hermitian-definite generalized
  eigenproblem to standard form.
*/

#define LAPACK_hegst(prefix, name)                                      \
  int F77_ ## prefix ## name ## gst( integer*  itype,                       \
                                     char* uplo,                        \
                                     integer*  m,                           \
                                     PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                                     PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, integer* ldim_B, \
                                     integer*  info )

#define LAPACK_hegst_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Inv      inv_fla;                                         \
  FLA_Uplo     uplo_fla;                                        \
  FLA_Obj      A, B;                                            \
  FLA_Error    init_result;                                     \
                                                                \
  FLA_Init_safe( &init_result );                                \
                                                                \
  FLA_Param_map_netlib_to_flame_inv( (int *) itype, &inv_fla );         \
  FLA_Param_map_netlib_to_flame_uplo( uplo, &uplo_fla );        \
                                                                \
  FLA_Obj_create_without_buffer( datatype, *m, *m, &A );        \
  FLA_Obj_attach_buffer( buff_A, 1, *ldim_A, &A );              \
                                                                \
  FLA_Obj_create_without_buffer( datatype, *m, *m, &B );        \
  FLA_Obj_attach_buffer( buff_B, 1, *ldim_B, &B );              \
                                                                \
  FLA_Eig_gest( inv_fla, uplo_fla, A, B );                      \
                                                                \
  FLA_Obj_free_without_buffer( &A );                            \
  FLA_Obj_free_without_buffer( &B );                            \
                                                                \
  FLA_Finalize_safe( init_result );                             \
                                                                \
  *info = 0;                                                    \
                                                                \



LAPACK_hegst(s,sy)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hegst-ssygst inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "", *itype, *uplo, *m, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1( ssygst_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ),fla_error )
    }
    if(fla_error==LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(s)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_hegst(d,sy)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hegst-dsygst inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "", *itype, *uplo, *m, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1( dsygst_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(d)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_hegst(c,he)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("chegst inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "", *itype, *uplo, *m, *ldim_A, *ldim_B);
#if !FLA_ENABLE_AMD_OPT 
    int fla_error = LAPACK_SUCCESS;   
    {
        LAPACK_RETURN_CHECK_VAR1( chegst_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(c)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        chegst_fla( itype, uplo,
                    m,
                    (complex *) buff_A, ldim_A,
                    (complex *) buff_B, ldim_B,
                    info );
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}
LAPACK_hegst(z,he)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhegst inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "", *itype, *uplo, *m, *ldim_A, *ldim_B);
#if !FLA_ENABLE_AMD_OPT  
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( zhegst_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ), fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(z)
        /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        zhegst_fla( itype, uplo,
                    m,
                    (doublecomplex *) buff_A, ldim_A,
                    (doublecomplex *) buff_B, ldim_B,
                    info );
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}


#define LAPACK_hegs2(prefix, name)                                      \
  int F77_ ## prefix ## name ## gs2(integer*  itype,                        \
                                    char* uplo,                         \
                                    integer*  m,                            \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, integer* ldim_A, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, integer* ldim_B, \
                                    integer*  info )

LAPACK_hegs2(s,sy)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hegs2-ssygs2 inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "", *itype, *uplo, *m, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1( ssygs2_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ), fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(s)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_hegs2(d,sy)
{
    int fla_error = LAPACK_SUCCESS;
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("hegs2-dsygs2 inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "", *itype, *uplo, *m, *ldim_A, *ldim_B);
    {
        LAPACK_RETURN_CHECK_VAR1( dsygs2_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(d)
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
}
LAPACK_hegs2(c,he)
{
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("chegs2 inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "", *itype, *uplo, *m, *ldim_A, *ldim_B);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS; 
    {
        LAPACK_RETURN_CHECK_VAR1( chegs2_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(c)        
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        chegs2_fla( itype, uplo,
                    m,
                    (complex *) buff_A, ldim_A,
                    (complex *) buff_B, ldim_B,
                    info );
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}
LAPACK_hegs2(z,he)
{ 
    AOCL_DTL_TRACE_LOG_INIT
    AOCL_DTL_SNPRINTF("zhegs2 inputs: itype %" FLA_IS ", uplo %c, n %" FLA_IS ", lda %" FLA_IS ", ldb %" FLA_IS "", *itype, *uplo, *m, *ldim_A, *ldim_B);
#if !FLA_ENABLE_AMD_OPT
    int fla_error = LAPACK_SUCCESS;
    {
        LAPACK_RETURN_CHECK_VAR1( zhegs2_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ),fla_error )
    }
    if (fla_error == LAPACK_SUCCESS)
    {
        LAPACK_hegst_body(z)        
         /** fla_error set to 0 on LAPACK_SUCCESS */
        fla_error = 0;
    }
    AOCL_DTL_TRACE_LOG_EXIT
    return fla_error;
#else
    {
        zhegs2_fla( itype, uplo,
                    m,
                    (doublecomplex *) buff_A, ldim_A,
                    (doublecomplex *) buff_B, ldim_B,
                    info );
        AOCL_DTL_TRACE_LOG_EXIT
        return 0;
    }
#endif
}



#endif
