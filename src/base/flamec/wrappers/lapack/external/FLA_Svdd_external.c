/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Svdd_external( FLA_Svd_type jobz, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  integer      info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  FLA_Datatype dt_int;
  integer          m_A, n_A, cs_A;
  integer          cs_U;
  integer          cs_V;
  integer          min_m_n;
  integer          lwork, lrwork, liwork;
  FLA_Obj      work, rwork, iwork;
  char         blas_jobz;
  integer          i;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Svdd_check( jobz, A, s, U, V );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );
  dt_real  = FLA_Obj_datatype_proj_to_real( A );
  dt_int   = FLA_INT;

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  cs_U     = FLA_Obj_col_stride( U );

  cs_V     = FLA_Obj_col_stride( V );

  min_m_n  = fla_min( m_A, n_A );

  // Allocate the rwork and iwork arrays up front.
  if ( jobz == FLA_SVD_VECTORS_NONE ) lrwork   = 5 * min_m_n;
  else                                lrwork   = 5 * min_m_n * min_m_n +
                                                 7 * min_m_n;
  liwork = 8 * min_m_n;

  FLA_Obj_create( dt_int,  liwork, 1, 0, 0, &iwork );
  if ( FLA_Obj_is_complex( A ) )
    FLA_Obj_create( dt_real, lrwork, 1, 0, 0, &rwork );

  FLA_Param_map_flame_to_netlib_svd_type( jobz, &blas_jobz );

  // Make a workspace query the first time through. This will provide us with
  // and ideal workspace size based on an internal block size.
  lwork = -1;
  FLA_Obj_create( datatype, 1, 1, 0, 0, &work );

  for ( i = 0; i < 2; ++i )
  {
    if ( i == 1 )
    {
      // Grab the queried ideal workspace size from the work array, free the
      // work object, and then re-allocate the workspace with the ideal size.
      if      ( datatype == FLA_FLOAT || datatype == FLA_COMPLEX )
        lwork = ( integer ) *FLA_FLOAT_PTR( work );
      else if ( datatype == FLA_DOUBLE || datatype == FLA_DOUBLE_COMPLEX )
        lwork = ( integer ) *FLA_DOUBLE_PTR( work );

      FLA_Obj_free( &work );
      FLA_Obj_create( datatype, lwork, 1, 0, 0, &work );
    }

    switch( datatype ) {

    case FLA_FLOAT:
    {
      float*    buff_A     = ( float*    ) FLA_FLOAT_PTR( A );
      float*    buff_s     = ( float*    ) FLA_FLOAT_PTR( s );
      float*    buff_U     = ( float*    ) FLA_FLOAT_PTR( U );
      float*    buff_V     = ( float*    ) FLA_FLOAT_PTR( V );
      float*    buff_work  = ( float*    ) FLA_FLOAT_PTR( work );
      integer*      buff_iwork = ( integer*      ) FLA_INT_PTR( iwork );
  
      F77_sgesdd( &blas_jobz,
                  &m_A,
                  &n_A,
                  buff_A,    &cs_A,
                  buff_s,
                  buff_U,    &cs_U,
                  buff_V,    &cs_V,
                  buff_work, &lwork,
                  buff_iwork,
                  &info );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A     = ( double*   ) FLA_DOUBLE_PTR( A );
      double*   buff_s     = ( double*   ) FLA_DOUBLE_PTR( s );
      double*   buff_U     = ( double*   ) FLA_DOUBLE_PTR( U );
      double*   buff_V     = ( double*   ) FLA_DOUBLE_PTR( V );
      double*   buff_work  = ( double*   ) FLA_DOUBLE_PTR( work );
      integer*      buff_iwork = ( integer*      ) FLA_INT_PTR( iwork );
  
      F77_dgesdd( &blas_jobz,
                  &m_A,
                  &n_A,
                  buff_A,    &cs_A,
                  buff_s,
                  buff_U,    &cs_U,
                  buff_V,    &cs_V,
                  buff_work, &lwork,
                  buff_iwork,
                  &info );
  
      break;
    } 
  
    case FLA_COMPLEX:
    {
      scomplex* buff_A     = ( scomplex* ) FLA_COMPLEX_PTR( A );
      float*    buff_s     = ( float*    ) FLA_FLOAT_PTR( s );
      scomplex* buff_U     = ( scomplex* ) FLA_COMPLEX_PTR( U );
      scomplex* buff_V     = ( scomplex* ) FLA_COMPLEX_PTR( V );
      scomplex* buff_work  = ( scomplex* ) FLA_COMPLEX_PTR( work );
      float*    buff_rwork = ( float*    ) FLA_FLOAT_PTR( rwork );
      integer*      buff_iwork = ( integer*      ) FLA_INT_PTR( iwork );
  
      F77_cgesdd( &blas_jobz,
                  &m_A,
                  &n_A,
                  buff_A,    &cs_A,
                  buff_s,
                  buff_U,    &cs_U,
                  buff_V,    &cs_V,
                  buff_work, &lwork,
                  buff_rwork,
                  buff_iwork,
                  &info );
  
      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A     = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_s     = ( double*   ) FLA_DOUBLE_PTR( s );
      dcomplex* buff_U     = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( U );
      dcomplex* buff_V     = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( V );
      dcomplex* buff_work  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( work );
      double*   buff_rwork = ( double*   ) FLA_DOUBLE_PTR( rwork );
      integer*      buff_iwork = ( integer*      ) FLA_INT_PTR( iwork );
  
      F77_zgesdd( &blas_jobz,
                  &m_A,
                  &n_A,
                  buff_A,    &cs_A,
                  buff_s,
                  buff_U,    &cs_U,
                  buff_V,    &cs_V,
                  buff_work, &lwork,
                  buff_rwork,
                  buff_iwork,
                  &info );
  
      break;
    } 

    }
  }

  FLA_Obj_free( &work );
  FLA_Obj_free( &iwork );
  if ( FLA_Obj_is_complex( A ) )
    FLA_Obj_free( &rwork );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

