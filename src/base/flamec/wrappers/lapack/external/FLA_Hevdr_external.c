
#include "FLAME.h"

FLA_Error FLA_Hevdr_external( FLA_Evd_type jobz, FLA_Uplo uplo, FLA_Obj A, FLA_Obj l, FLA_Obj Z )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  int          n_A, cs_A;
  int               cs_Z;
  int          lwork, lrwork, liwork, lisuppz;
  FLA_Obj      work, rwork, iwork, isuppz, abstol;
  char         blas_jobz;
  char         blas_uplo;
  int          i;

  char         blas_range = 'A';
  int          il, iu;
  int          eigs_found;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Hevdr_check( jobz, uplo, A, l, Z );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );
  dt_real  = FLA_Obj_datatype_proj_to_real( A );

  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );
  cs_Z     = FLA_Obj_col_stride( Z );

  FLA_Param_map_flame_to_netlib_evd_type( jobz, &blas_jobz );
  FLA_Param_map_flame_to_netlib_uplo( uplo, &blas_uplo );

  lisuppz = 2 * n_A;
  FLA_Obj_create( FLA_INT,  lisuppz, 1, 0, 0, &isuppz );
  FLA_Obj_create( dt_real,  1,       1, 0, 0, &abstol );

  // Query the safe minimum to use as the abstol parameter.
  FLA_Mach_params( FLA_MACH_SFMIN, abstol );

  // Make a workspace query the first time through. This will provide us with
  // and ideal workspace size.
  lwork = -1;
  lrwork = -1;
  liwork = -1;
  FLA_Obj_create( datatype, 1, 1, 0, 0, &work );
  FLA_Obj_create( dt_real,  1, 1, 0, 0, &rwork );
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

      FLA_Obj_free( &work );
      FLA_Obj_free( &iwork );
      FLA_Obj_free( &rwork );
      FLA_Obj_create( datatype, lwork,  1, 0, 0, &work );
      FLA_Obj_create( FLA_INT,  liwork, 1, 0, 0, &iwork );
      if ( FLA_Obj_is_complex( A ) )
        FLA_Obj_create( dt_real,  lrwork, 1, 0, 0, &rwork );
    }

    switch( datatype ) {

    case FLA_FLOAT:
    {
      float*    buff_A      = ( float*    ) FLA_FLOAT_PTR( A );
      float*    buff_l      = ( float*    ) FLA_FLOAT_PTR( l );
      float*    buff_Z      = ( float*    ) FLA_FLOAT_PTR( Z );
      float*    buff_work   = ( float*    ) FLA_FLOAT_PTR( work );
      float*    buff_abstol = ( float*    ) FLA_FLOAT_PTR( abstol );
      int*      buff_iwork  = ( int*      ) FLA_INT_PTR( iwork );
      int*      buff_isuppz = ( int*      ) FLA_INT_PTR( isuppz );
      float     vl, vu;

      F77_ssyevr( &blas_jobz,
                  &blas_range,
                  &blas_uplo,
                  &n_A,
                  buff_A,     &cs_A,
                  &vl, &vu,
                  &il, &iu,
                  buff_abstol,
                  &eigs_found,
                  buff_l,
                  buff_Z,     &cs_Z,
                  buff_isuppz,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );
      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_A      = ( double*   ) FLA_DOUBLE_PTR( A );
      double*   buff_l      = ( double*   ) FLA_DOUBLE_PTR( l );
      double*   buff_Z      = ( double*   ) FLA_DOUBLE_PTR( Z );
      double*   buff_work   = ( double*   ) FLA_DOUBLE_PTR( work );
      double*   buff_abstol = ( double*   ) FLA_DOUBLE_PTR( abstol );
      int*      buff_iwork  = ( int*      ) FLA_INT_PTR( iwork );
      int*      buff_isuppz = ( int*      ) FLA_INT_PTR( isuppz );
      double    vl, vu;
  
      F77_dsyevr( &blas_jobz,
                  &blas_range,
                  &blas_uplo,
                  &n_A,
                  buff_A,     &cs_A,
                  &vl, &vu,
                  &il, &iu,
                  buff_abstol,
                  &eigs_found,
                  buff_l,
                  buff_Z,     &cs_Z,
                  buff_isuppz,
                  buff_work,  &lwork,
                  buff_iwork, &liwork,
                  &info );
      break;
    } 
  
    case FLA_COMPLEX:
    {
      scomplex* buff_A      = ( scomplex* ) FLA_COMPLEX_PTR( A );
      float*    buff_l      = ( float*    ) FLA_FLOAT_PTR( l );
      scomplex* buff_Z      = ( scomplex* ) FLA_COMPLEX_PTR( Z );
      scomplex* buff_work   = ( scomplex* ) FLA_COMPLEX_PTR( work );
      float*    buff_rwork  = ( float*    ) FLA_FLOAT_PTR( rwork );
      float*    buff_abstol = ( float*    ) FLA_FLOAT_PTR( abstol );
      int*      buff_iwork  = ( int*      ) FLA_INT_PTR( iwork );
      int*      buff_isuppz = ( int*      ) FLA_INT_PTR( isuppz );
      float     vl, vu;
  
      F77_cheevr( &blas_jobz,
                  &blas_range,
                  &blas_uplo,
                  &n_A,
                  buff_A,     &cs_A,
                  &vl, &vu,
                  &il, &iu,
                  buff_abstol,
                  &eigs_found,
                  buff_l,
                  buff_Z,     &cs_Z,
                  buff_isuppz,
                  buff_work,  &lwork,
                  buff_rwork, &lrwork,
                  buff_iwork, &liwork,
                  &info );
      break;
    } 
  
    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A      = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );
      double*   buff_l      = ( double*   ) FLA_DOUBLE_PTR( l );
      dcomplex* buff_Z      = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( Z );
      dcomplex* buff_work   = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( work );
      double*   buff_rwork  = ( double*   ) FLA_DOUBLE_PTR( rwork );
      double*   buff_abstol = ( double*   ) FLA_DOUBLE_PTR( abstol );
      int*      buff_iwork  = ( int*      ) FLA_INT_PTR( iwork );
      int*      buff_isuppz = ( int*      ) FLA_INT_PTR( isuppz );
      double    vl, vu;
  
      F77_zheevr( &blas_jobz,
                  &blas_range,
                  &blas_uplo,
                  &n_A,
                  buff_A,     &cs_A,
                  &vl, &vu,
                  &il, &iu,
                  buff_abstol,
                  &eigs_found,
                  buff_l,
                  buff_Z,     &cs_Z,
                  buff_isuppz,
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
  FLA_Obj_free( &isuppz );
  FLA_Obj_free( &abstol );
  if ( FLA_Obj_is_complex( A ) )
    FLA_Obj_free( &rwork );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

