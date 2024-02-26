/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tevdr_external( FLA_Evd_type jobz, FLA_Obj d, FLA_Obj e, FLA_Obj l, FLA_Obj A )
{
  integer      info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  integer          n_A, cs_A;
  integer          lisuppz, lwork, liwork;
  FLA_Obj      isuppz, work, iwork;
  char         blas_jobz;
  char         blas_range;
  integer          i;
  integer          vl, vu;
  integer          il, iu;
  integer          nzc;
  integer          try_rac;
  integer          n_eig_found;

  //if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
  //  FLA_Tevdd_check( jobz, d, e, A );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );
  dt_real  = FLA_Obj_datatype_proj_to_real( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Param_map_flame_to_netlib_evd_type( jobz, &blas_jobz );

  vl = 0;
  vu = 0;

  // Hard-code some parameters.
  blas_range = 'A';
  nzc        = n_A;
  try_rac    = TRUE;

  // Allocate space for the isuppz array.
  lisuppz = 2 * n_A;
  FLA_Obj_create( FLA_INT, lisuppz, 1, 0, 0, &isuppz );

  // Make a workspace query the first time through. This will provide us with
  // and ideal workspace size.
  lwork = -1;
  liwork = -1;
  FLA_Obj_create( dt_real, 1, 1, 0, 0, &work );
  FLA_Obj_create( FLA_INT, 1, 1, 0, 0, &iwork );

  for ( i = 0; i < 2; ++i )
  {
    if ( i == 1 )
    {
      // Grab the queried ideal workspace size from the work arrays, free the
      // work object, and then re-allocate the workspace with the ideal size.
      if      ( datatype == FLA_FLOAT || datatype == FLA_COMPLEX )
      {
        lwork  = ( integer ) *FLA_FLOAT_PTR( work );
        liwork = ( integer ) *FLA_INT_PTR( iwork );
      }
      else if ( datatype == FLA_DOUBLE || datatype == FLA_DOUBLE_COMPLEX )
      {
        lwork  = ( integer ) *FLA_DOUBLE_PTR( work );
        liwork = ( integer ) *FLA_INT_PTR( iwork );
      }
//printf( "ideal workspace for n = %d\n", n_A );
//printf( "                lwork = %d\n", lwork );
//printf( "               liwork = %d\n", liwork );
      FLA_Obj_free( &work );
      FLA_Obj_free( &iwork );
      FLA_Obj_create( dt_real, lwork,  1, 0, 0, &work );
      FLA_Obj_create( FLA_INT, liwork, 1, 0, 0, &iwork );
    }

    switch( datatype ) {

    case FLA_FLOAT:
    {
      float*    buff_d      = ( float*    ) FLA_FLOAT_PTR( d );
      float*    buff_e      = ( float*    ) FLA_FLOAT_PTR( e );
      float*    buff_l      = ( float*    ) FLA_FLOAT_PTR( l );
      float*    buff_A      = ( float*    ) FLA_FLOAT_PTR( A );
      integer*      buff_isuppz = ( integer*      ) FLA_INT_PTR( isuppz );
      float*    buff_work   = ( float*    ) FLA_FLOAT_PTR( work );
      integer*      buff_iwork  = ( integer*      ) FLA_INT_PTR( iwork );
      float vlf = (float) vl;
      float vuf = (float) vu;

      F77_sstemr( &blas_jobz,
                  &blas_range,
                  &n_A,
                  buff_d,
                  buff_e,
                  &vlf, &vuf,
                  &il, &iu,
                  &n_eig_found,
                  buff_l,
                  buff_A,     &cs_A,
                  &nzc,
                  buff_isuppz,
                  &try_rac,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_d      = ( double*   ) FLA_DOUBLE_PTR( d );
      double*   buff_e      = ( double*   ) FLA_DOUBLE_PTR( e );
      double*   buff_l      = ( double*   ) FLA_DOUBLE_PTR( l );
      double*   buff_A      = ( double*   ) FLA_DOUBLE_PTR( A );
      integer*      buff_isuppz = ( integer*      ) FLA_INT_PTR( isuppz );
      double*   buff_work   = ( double*   ) FLA_DOUBLE_PTR( work );
      integer*      buff_iwork  = ( integer*      ) FLA_INT_PTR( iwork );
      double vlf = (double) vl;
      double vuf = (double) vu;

      F77_dstemr( &blas_jobz,
                  &blas_range,
                  &n_A,
                  buff_d,
                  buff_e,
                  &vlf, &vuf,
                  &il, &iu,
                  &n_eig_found,
                  buff_l,
                  buff_A,     &cs_A,
                  &nzc,
                  buff_isuppz,
                  &try_rac,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );
  
      break;
    } 
  
    case FLA_COMPLEX:
    {
      float*    buff_d      = ( float*    ) FLA_FLOAT_PTR( d );
      float*    buff_e      = ( float*    ) FLA_FLOAT_PTR( e );
      float*    buff_l      = ( float*    ) FLA_FLOAT_PTR( l );
      scomplex* buff_A      = ( scomplex* ) FLA_COMPLEX_PTR( A );
      integer*      buff_isuppz = ( integer*      ) FLA_INT_PTR( isuppz );
      float*    buff_work   = ( float*    ) FLA_FLOAT_PTR( work );
      integer*      buff_iwork  = ( integer*      ) FLA_INT_PTR( iwork );
      float vlf = (float) vl;
      float vuf = (float) vu;

      F77_cstemr( &blas_jobz,
                  &blas_range,
                  &n_A,
                  buff_d,
                  buff_e,
                  &vlf, &vuf,
                  &il, &iu,
                  &n_eig_found,
                  buff_l,
                  buff_A,     &cs_A,
                  &nzc,
                  buff_isuppz,
                  &try_rac,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );
  
      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      double*   buff_d      = ( double*   ) FLA_DOUBLE_PTR( d );
      double*   buff_e      = ( double*   ) FLA_DOUBLE_PTR( e );
      double*   buff_l      = ( double*   ) FLA_DOUBLE_PTR( l );
      dcomplex* buff_A      = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );
      integer*      buff_isuppz = ( integer*      ) FLA_INT_PTR( isuppz );
      double*   buff_work   = ( double*   ) FLA_DOUBLE_PTR( work );
      integer*      buff_iwork  = ( integer*      ) FLA_INT_PTR( iwork );
      double vlf = (double) vl;
      double vuf = (double) vu;

      F77_zstemr( &blas_jobz,
                  &blas_range,
                  &n_A,
                  buff_d,
                  buff_e,
                  &vlf, &vuf,
                  &il, &iu,
                  &n_eig_found,
                  buff_l,
                  buff_A,     &cs_A,
                  &nzc,
                  buff_isuppz,
                  &try_rac,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );
  
      break;
    } 

    }
  }

  FLA_Obj_free( &isuppz );
  FLA_Obj_free( &work );
  FLA_Obj_free( &iwork );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

