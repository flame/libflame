/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_LAPACK2FLASH

#include "FLASH_lapack2flash_util_defs.h"
#include "FLA_lapack2flame_return_defs.h"
#include "FLA_lapack2flame_prototypes.h"

/*
  ZHEGST reduces a complex Hermitian-definite generalized
  eigenproblem to standard form.
*/

#define LAPACK_hegst(prefix, name)                                      \
  int F77_ ## prefix ## name ## gst( int*  itype,                       \
                                     char* uplo,                        \
                                     int*  m,                           \
                                     PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                                     PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, int* ldim_B, \
                                     int*  info )

#define LAPACK_hegst_body(prefix)                               \
  FLA_Datatype datatype = PREFIX2FLAME_DATATYPE(prefix);        \
  FLA_Inv      inv_fla;                                         \
  FLA_Uplo     uplo_fla;                                        \
  FLA_Obj      A, B;                                            \
  FLA_Error    init_result;                                     \
                                                                \
  FLA_Init_safe( &init_result );                                \
                                                                \
  FLA_Param_map_netlib_to_flame_inv( itype, &inv_fla );         \
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
  return 0;


LAPACK_hegst(s,sy)
{
    {
        LAPACK_RETURN_CHECK( ssygst_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ) )
    }
    {
        LAPACK_hegst_body(s)
    }
}
LAPACK_hegst(d,sy)
{
    {
        LAPACK_RETURN_CHECK( dsygst_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ) )
    }
    {
        LAPACK_hegst_body(d)
    }
}
LAPACK_hegst(c,he)
{
    {
        LAPACK_RETURN_CHECK( chegst_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ) )
    }
    {
        LAPACK_hegst_body(c)
    }
}
LAPACK_hegst(z,he)
{
    {
        LAPACK_RETURN_CHECK( zhegst_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ) )
    }
    {
        LAPACK_hegst_body(z)
    }
}


#define LAPACK_hegs2(prefix, name)                                      \
  int F77_ ## prefix ## name ## gs2(int*  itype,                        \
                                    char* uplo,                         \
                                    int*  m,                            \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_A, int* ldim_A, \
                                    PREFIX2LAPACK_TYPEDEF(prefix)* buff_B, int* ldim_B, \
                                    int*  info )

LAPACK_hegs2(s,sy)
{
    {
        LAPACK_RETURN_CHECK( ssygs2_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ) )
    }
    {
        LAPACK_hegst_body(s)
    }
}
LAPACK_hegs2(d,sy)
{
    {
        LAPACK_RETURN_CHECK( dsygs2_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ) )
    }
    {
        LAPACK_hegst_body(d)
    }
}
LAPACK_hegs2(c,he)
{
    {
        LAPACK_RETURN_CHECK( chegs2_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ) )
    }
    {
        LAPACK_hegst_body(c)
    }
}
LAPACK_hegs2(z,he)
{
    {
        LAPACK_RETURN_CHECK( zhegs2_check( itype, uplo,
                                           m,
                                           buff_A, ldim_A,
                                           buff_B, ldim_B,
                                           info ) )
    }
    {
        LAPACK_hegst_body(z)
    }
}



#endif
