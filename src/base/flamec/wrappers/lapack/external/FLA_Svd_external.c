
#include "FLAME.h"

FLA_Error FLA_Svd_external( FLA_Svd_type jobu, FLA_Svd_type jobv, FLA_Obj A, FLA_Obj s, FLA_Obj U, FLA_Obj V )
{
  int          info = 0;
#ifdef FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES
  FLA_Datatype datatype;
  FLA_Datatype dt_real;
  int          m_A, n_A, cs_A;
  int          cs_U;
  int          cs_V;
  int          min_m_n;
  int          lwork, lrwork;
  FLA_Obj      work, rwork;
  char         blas_jobu;
  char         blas_jobv;
  int          i;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Svd_check( jobu, jobv, A, s, U, V );

  if ( FLA_Obj_has_zero_dim( A ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );
  dt_real  = FLA_Obj_datatype_proj_to_real( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  cs_A     = FLA_Obj_col_stride( A );

  cs_U     = FLA_Obj_col_stride( U );

  cs_V     = FLA_Obj_col_stride( V );

  min_m_n  = min( m_A, n_A );

  // Allocate the rwork array up front since its size is not dependent on
  // internal block sizes.
  lrwork   = 5 * min_m_n;
  if ( FLA_Obj_is_complex( A ) )
    FLA_Obj_create( dt_real, lrwork, 1, 0, 0, &rwork );

  FLA_Param_map_flame_to_netlib_svd_type( jobu, &blas_jobu );
  FLA_Param_map_flame_to_netlib_svd_type( jobv, &blas_jobv );

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
        lwork = ( int ) *FLA_FLOAT_PTR( work );
      else if ( datatype == FLA_DOUBLE || datatype == FLA_DOUBLE_COMPLEX )
        lwork = ( int ) *FLA_DOUBLE_PTR( work );

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
  
      F77_sgesvd( &blas_jobu,
                  &blas_jobv,
                  &m_A,
                  &n_A,
                  buff_A,    &cs_A,
                  buff_s,
                  buff_U,    &cs_U,
                  buff_V,    &cs_V,
                  buff_work, &lwork,
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
  
      F77_dgesvd( &blas_jobu,
                  &blas_jobv,
                  &m_A,
                  &n_A,
                  buff_A,    &cs_A,
                  buff_s,
                  buff_U,    &cs_U,
                  buff_V,    &cs_V,
                  buff_work, &lwork,
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
  
      F77_cgesvd( &blas_jobu,
                  &blas_jobv,
                  &m_A,
                  &n_A,
                  buff_A,    &cs_A,
                  buff_s,
                  buff_U,    &cs_U,
                  buff_V,    &cs_V,
                  buff_work, &lwork,
                  buff_rwork,
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
  
      F77_zgesvd( &blas_jobu,
                  &blas_jobv,
                  &m_A,
                  &n_A,
                  buff_A,    &cs_A,
                  buff_s,
                  buff_U,    &cs_U,
                  buff_V,    &cs_V,
                  buff_work, &lwork,
                  buff_rwork,
                  &info );
  
      break;
    } 

    }
  }

  FLA_Obj_free( &work );
  if ( FLA_Obj_is_complex( A ) )
    FLA_Obj_free( &rwork );
#else
  FLA_Check_error_code( FLA_EXTERNAL_LAPACK_NOT_IMPLEMENTED );
#endif

  return info;
}

