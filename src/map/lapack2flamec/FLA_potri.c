
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
                               int*  n,                         \
                               PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, \
                               int*  ldim_A,                         \
                               int*  info )

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
  return 0;

LAPACK_potri(s)
{
    {
        LAPACK_RETURN_CHECK( spotri_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potri_body(s)
    }
}
LAPACK_potri(d)
{
    {
        LAPACK_RETURN_CHECK( dpotri_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potri_body(d)
    }
}
LAPACK_potri(c)
{
    {
        LAPACK_RETURN_CHECK( cpotri_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potri_body(c)
    }
}
LAPACK_potri(z)
{
    {
        LAPACK_RETURN_CHECK( zpotri_check( uplo, n,
                                           buff_A, ldim_A,
                                           info ) )
    }
    {
        LAPACK_potri_body(z)
    }
}

#endif
