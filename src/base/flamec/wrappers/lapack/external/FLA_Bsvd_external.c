/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bsvd_external( FLA_Uplo uplo, FLA_Obj d, FLA_Obj e, FLA_Obj U, FLA_Obj V )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  integer          m_U, cs_U;
  integer          n_V, cs_V;
  integer          n_C, cs_C;
  integer          min_m_n;
  integer          lrwork;
  FLA_Obj      rwork;
  char         blas_uplo;

  //if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
  //  FLA_Hevd_check( jobz, uplo, A, e );

  if ( FLA_Obj_has_zero_dim( d ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( U );
  dt_real  = FLA_Obj_datatype_proj_to_real( U );

  m_U      = FLA_Obj_length( U );
  cs_U     = FLA_Obj_col_stride( U );

  n_V      = FLA_Obj_length( V );
  cs_V     = FLA_Obj_col_stride( V );

  n_C      = 0;
  cs_C     = 1;

  min_m_n  = FLA_Obj_vector_dim( d );

  lrwork   = fla_max( 1, 4 * min_m_n - 4 );
  FLA_Obj_create( dt_real, lrwork, 1, 0, 0, &rwork );

  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );

    switch( datatype ) {

    case FLA_FLOAT:
    {
      float*    buff_d     = ( float * ) FLA_FLOAT_PTR( d );
      float*    buff_e     = ( float * ) FLA_FLOAT_PTR( e );
      float*    buff_U     = ( float * ) FLA_FLOAT_PTR( U );
      float*    buff_V     = ( float * ) FLA_FLOAT_PTR( V );
      float*    buff_C     = ( float * ) NULL;
      float*    buff_rwork = ( float * ) FLA_FLOAT_PTR( rwork );
  
      F77_sbdsqr( &blas_uplo,
                  &min_m_n,
                  &n_V,
                  &m_U,
                  &n_C,
                  buff_d,
                  buff_e,
                  buff_V, &cs_V,
                  buff_U, &cs_U,
                  buff_C, &cs_C,
                  buff_rwork,
                  &info );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_d     = ( double * ) FLA_DOUBLE_PTR( d );
      double*   buff_e     = ( double * ) FLA_DOUBLE_PTR( e );
      double*   buff_U     = ( double * ) FLA_DOUBLE_PTR( U );
      double*   buff_V     = ( double * ) FLA_DOUBLE_PTR( V );
      double*   buff_C     = ( double * ) NULL;
      double*   buff_rwork = ( double * ) FLA_DOUBLE_PTR( rwork );
  
      F77_dbdsqr( &blas_uplo,
                  &min_m_n,
                  &n_V,
                  &m_U,
                  &n_C,
                  buff_d,
                  buff_e,
                  buff_V, &cs_V,
                  buff_U, &cs_U,
                  buff_C, &cs_C,
                  buff_rwork,
                  &info );

      break;
    } 
  
    case FLA_COMPLEX:
    {
      float*    buff_d     = ( float    * ) FLA_FLOAT_PTR( d );
      float*    buff_e     = ( float    * ) FLA_FLOAT_PTR( e );
      scomplex* buff_U     = ( scomplex * ) FLA_COMPLEX_PTR( U );
      scomplex* buff_V     = ( scomplex * ) FLA_COMPLEX_PTR( V );
      scomplex* buff_C     = ( scomplex * ) NULL;
      float*    buff_rwork = ( float    * ) FLA_FLOAT_PTR( rwork );
  
      F77_cbdsqr( &blas_uplo,
                  &min_m_n,
                  &n_V,
                  &m_U,
                  &n_C,
                  buff_d,
                  buff_e,
                  buff_V, &cs_V,
                  buff_U, &cs_U,
                  buff_C, &cs_C,
                  buff_rwork,
                  &info );

      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      double*   buff_d     = ( double   * ) FLA_DOUBLE_PTR( d );
      double*   buff_e     = ( double   * ) FLA_DOUBLE_PTR( e );
      dcomplex* buff_U     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( U );
      dcomplex* buff_V     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( V );
      dcomplex* buff_C     = ( dcomplex * ) NULL;
      double*   buff_rwork = ( double   * ) FLA_DOUBLE_PTR( rwork );
  
      F77_zbdsqr( &blas_uplo,
                  &min_m_n,
                  &n_V,
                  &m_U,
                  &n_C,
                  buff_d,
                  buff_e,
                  buff_V, &cs_V,
                  buff_U, &cs_U,
                  buff_C, &cs_C,
                  buff_rwork,
                  &info );

      break;
    } 

    }

  FLA_Obj_free( &rwork );

#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

