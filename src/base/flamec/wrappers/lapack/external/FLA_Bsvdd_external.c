/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bsvdd_external( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, FLA_Obj U, FLA_Obj V )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  integer          cs_U;
  integer          cs_V;
  integer          min_m_n;
  integer          lwork, liwork;
  FLA_Obj      work, iwork;
  char         blas_uplo;
  char         blas_compq = 'I';

  //if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
  //  FLA_Bsvd_check( uplo, d, e, U, V );

  if ( FLA_Obj_has_zero_dim( d ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( U );
  dt_real  = FLA_Obj_datatype_proj_to_real( U );

  cs_U     = FLA_Obj_col_stride( U );

  cs_V     = FLA_Obj_col_stride( V );

  min_m_n  = FLA_Obj_vector_dim( d );

  lwork   = fla_max( 1, 3*min_m_n*min_m_n + 4*min_m_n );
  liwork  = 8*min_m_n;

  FLA_Obj_create( dt_real, lwork,  1, 0, 0, &work );
  FLA_Obj_create( FLA_INT, liwork, 1, 0, 0, &iwork );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );

    switch( datatype ) {

    case FLA_FLOAT:
    {
      float*    buff_d     = ( float * ) FLA_FLOAT_PTR( d );
      float*    buff_e     = ( float * ) FLA_FLOAT_PTR( e );
      float*    buff_U     = ( float * ) FLA_FLOAT_PTR( U );
      float*    buff_V     = ( float * ) FLA_FLOAT_PTR( V );
      float*    buff_Q     = ( float * ) NULL;
      float*    buff_IQ    = ( float * ) NULL;
      float*    buff_work  = ( float * ) FLA_FLOAT_PTR( work );
      integer*      buff_iwork = ( integer   * ) FLA_INT_PTR( iwork );
  
      F77_sbdsdc( &blas_uplo,
                  &blas_compq,
                  &min_m_n,
                  buff_d,
                  buff_e,
                  buff_U, &cs_U,
                  buff_V, &cs_V,
                  buff_Q,
                  buff_IQ,
                  buff_work,
                  buff_iwork,
                  &info );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_d     = ( double * ) FLA_DOUBLE_PTR( d );
      double*   buff_e     = ( double * ) FLA_DOUBLE_PTR( e );
      double*   buff_U     = ( double * ) FLA_DOUBLE_PTR( U );
      double*   buff_V     = ( double * ) FLA_DOUBLE_PTR( V );
      double*   buff_Q     = ( double * ) NULL;
      double*   buff_IQ    = ( double * ) NULL;
      double*   buff_work  = ( double * ) FLA_DOUBLE_PTR( work );
      integer*      buff_iwork = ( integer    * ) FLA_INT_PTR( iwork );
  
      F77_dbdsdc( &blas_uplo,
                  &blas_compq,
                  &min_m_n,
                  buff_d,
                  buff_e,
                  buff_U, &cs_U,
                  buff_V, &cs_V,
                  buff_Q,
                  buff_IQ,
                  buff_work,
                  buff_iwork,
                  &info );

      break;
    } 
  
    }

  FLA_Obj_free( &work );
  FLA_Obj_free( &iwork );

#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

