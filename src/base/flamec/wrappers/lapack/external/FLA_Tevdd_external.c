/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tevdd_external( FLA_Evd_type jobz, FLA_Obj d, FLA_Obj e, FLA_Obj A )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  int          n_A, cs_A;
  int          lwork, lrwork, liwork;
  FLA_Obj      work, rwork, iwork;
  char         blas_jobz;
  int          i;

  //if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
  //  FLA_Tevdd_check( jobz, d, e, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );
  dt_real  = FLA_Obj_datatype_proj_to_real( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Param_map_flame_to_netlib_evd_type( jobz, &blas_jobz );

  // Make a workspace query the first time through. This will provide us with
  // and ideal workspace size.
  lwork = -1;
  lrwork = -1;
  liwork = -1;
  FLA_Obj_create( datatype, 1, 1, 0, 0, &work );
  FLA_Obj_create( datatype, 1, 1, 0, 0, &rwork );
  FLA_Obj_create( FLA_INT,  1, 1, 0, 0, &iwork );

  for ( i = 0; i < 2; ++i )
  {
    if ( i == 1 )
    {
      // Grab the queried ideal workspace size from the work arrays, free the
      // work object, and then re-allocate the workspace with the ideal size.
      if      ( datatype == FLA_FLOAT || datatype == FLA_COMPLEX )
      {
        lwork  = ( int ) *FLA_FLOAT_PTR( work );
        lrwork = ( int ) *FLA_FLOAT_PTR( rwork );
        liwork = ( int ) *FLA_INT_PTR( iwork );
      }
      else if ( datatype == FLA_DOUBLE || datatype == FLA_DOUBLE_COMPLEX )
      {
        lwork  = ( int ) *FLA_DOUBLE_PTR( work );
        lrwork = ( int ) *FLA_DOUBLE_PTR( rwork );
        liwork = ( int ) *FLA_INT_PTR( iwork );
      }

//printf( "ideal workspace for n = %d\n", n_A );
//printf( "                lwork = %d\n", lwork );
//printf( "               lrwork = %d\n", lrwork );
//printf( "               liwork = %d\n", liwork );

      FLA_Obj_free( &work );
      FLA_Obj_free( &iwork );
      FLA_Obj_free( &rwork );
      FLA_Obj_create( datatype, lwork,  1, 0, 0, &work );
      FLA_Obj_create( datatype, liwork, 1, 0, 0, &iwork );
      if ( FLA_Obj_is_complex( A ) )
        FLA_Obj_create( datatype, lrwork, 1, 0, 0, &rwork );
    }

    switch( datatype ) {

    case FLA_FLOAT:
    {
      float* buff_d     = ( float* ) FLA_FLOAT_PTR( d );
      float* buff_e     = ( float* ) FLA_FLOAT_PTR( e );
      float* buff_A     = ( float* ) FLA_FLOAT_PTR( A );
      float* buff_work  = ( float* ) FLA_FLOAT_PTR( work );
      int*   buff_iwork = ( int*   ) FLA_INT_PTR( iwork );

      F77_sstedc( &blas_jobz,
                  &n_A,
                  buff_d,
                  buff_e,
                  buff_A,     &cs_A,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_d     = ( double* ) FLA_DOUBLE_PTR( d );
      double* buff_e     = ( double* ) FLA_DOUBLE_PTR( e );
      double* buff_A     = ( double* ) FLA_DOUBLE_PTR( A );
      double* buff_work  = ( double* ) FLA_DOUBLE_PTR( work );
      int*    buff_iwork = ( int*    ) FLA_INT_PTR( iwork );
  
      F77_dstedc( &blas_jobz,
                  &n_A,
                  buff_d,
                  buff_e,
                  buff_A,     &cs_A,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );
  
      break;
    } 
  
    case FLA_COMPLEX:
    {
      float*    buff_d     = ( float*    ) FLA_FLOAT_PTR( d );
      float*    buff_e     = ( float*    ) FLA_FLOAT_PTR( e );
      scomplex* buff_A     = ( scomplex* ) FLA_COMPLEX_PTR( A );
      scomplex* buff_work  = ( scomplex* ) FLA_COMPLEX_PTR( work );
      float*    buff_rwork = ( float*    ) FLA_FLOAT_PTR( rwork );
      int*      buff_iwork = ( int*      ) FLA_INT_PTR( iwork );
  
      F77_cstedc( &blas_jobz,
                  &n_A,
                  buff_d,
                  buff_e,
                  buff_A,     &cs_A,
                  buff_work,  &lwork,
                  buff_rwork, &lrwork,
                  buff_iwork, &liwork,
                  &info );
  
      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      double*   buff_d     = ( double*   ) FLA_DOUBLE_PTR( d );
      double*   buff_e     = ( double*   ) FLA_DOUBLE_PTR( e );
      dcomplex* buff_A     = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_work  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( work );
      double*   buff_rwork = ( double*   ) FLA_DOUBLE_PTR( rwork );
      int*      buff_iwork = ( int*      ) FLA_INT_PTR( iwork );
  
      F77_zstedc( &blas_jobz,
                  &n_A,
                  buff_d,
                  buff_e,
                  buff_A,     &cs_A,
                  buff_work,  &lwork,
                  buff_rwork, &lrwork,
                  buff_iwork, &liwork,
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

