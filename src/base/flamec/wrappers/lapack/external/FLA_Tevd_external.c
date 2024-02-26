/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Tevd_external( FLA_Evd_type jobz, FLA_Obj d, FLA_Obj e, FLA_Obj A )
{
  integer      info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  integer          n_A, cs_A;
  integer          lwork;
  FLA_Obj      work, d_use, e_use;
  char         blas_jobz;

  //if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
  //  FLA_Hevd_check( jobz, uplo, A, e );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );
  dt_real  = FLA_Obj_datatype_proj_to_real( A );

  n_A      = FLA_Obj_vector_dim( d );
  cs_A     = FLA_Obj_col_stride( A );

  if ( FLA_Obj_vector_inc( d ) != 1 )
  {
    FLA_Obj_create( dt_real, n_A, 1, 0, 0, &d_use );
    FLA_Copy( d, d_use );
  }
  else { d_use = d; }

  if ( FLA_Obj_vector_inc( e ) != 1 )
  {
    FLA_Obj_create( dt_real, n_A-1, 1, 0, 0, &e_use );
    FLA_Copy( e, e_use );
  }
  else { e_use = e; }

  // Allocate thw work array up front.
  lwork   = fla_max( 1.0, 2.0 * n_A - 2 );
  FLA_Obj_create( dt_real, lwork, 1, 0, 0, &work );

  FLA_Param_map_flame_to_netlib_evd_type( jobz, &blas_jobz );

    switch( datatype ) {

    case FLA_FLOAT:
    {
      float* buff_A     = ( float * ) FLA_FLOAT_PTR( A );
      float* buff_d     = ( float * ) FLA_FLOAT_PTR( d_use );
      float* buff_e     = ( float * ) FLA_FLOAT_PTR( e_use );
      float* buff_work  = ( float * ) FLA_FLOAT_PTR( work );

      F77_ssteqr( &blas_jobz,
                  &n_A,
                  buff_d,
                  buff_e,
                  buff_A, &cs_A,
                  buff_work,
                  &info );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
      double* buff_d     = ( double * ) FLA_DOUBLE_PTR( d_use );
      double* buff_e     = ( double * ) FLA_DOUBLE_PTR( e_use );
      double* buff_work  = ( double * ) FLA_DOUBLE_PTR( work );
  
      F77_dsteqr( &blas_jobz,
                  &n_A,
                  buff_d,
                  buff_e,
                  buff_A, &cs_A,
                  buff_work,
                  &info );
  
      break;
    } 
  
    case FLA_COMPLEX:
    {
      scomplex* buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
      float*    buff_d     = ( float    * ) FLA_FLOAT_PTR( d_use );
      float*    buff_e     = ( float    * ) FLA_FLOAT_PTR( e_use );
      float*    buff_work  = ( float    * ) FLA_FLOAT_PTR( work );
  
      F77_csteqr( &blas_jobz,
                  &n_A,
                  buff_d,
                  buff_e,
                  buff_A, &cs_A,
                  buff_work,
                  &info );
  
      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A     = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_d     = ( double*   ) FLA_DOUBLE_PTR( d_use );
      double*   buff_e     = ( double*   ) FLA_DOUBLE_PTR( e_use );
      double*   buff_work  = ( double*   ) FLA_DOUBLE_PTR( work );
  
      F77_zsteqr( &blas_jobz,
                  &n_A,
                  buff_d,
                  buff_e,
                  buff_A, &cs_A,
                  buff_work,
                  &info );
  
      break;
    } 

    }

  if ( FLA_Obj_vector_inc( d ) != 1 )
  {
    FLA_Copy( d_use, d );
    FLA_Obj_free( &d_use );
  }
  if ( FLA_Obj_vector_inc( e ) != 1 )
  {
    FLA_Copy( e_use, e );
    FLA_Obj_free( &e_use );
  }

  FLA_Obj_free( &work );

#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

