
#include "FLAME.h"

FLA_Error FLA_Gemvc_external( FLA_Trans transa, FLA_Conj conjx, FLA_Obj alpha, FLA_Obj A, FLA_Obj x, FLA_Obj beta, FLA_Obj y )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_x;
  int          inc_y;
  trans1_t      blis_transa;
  conj1_t       blis_conjx;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Gemvc_check( transa, conjx, alpha, A, x, beta, y );

  if ( FLA_Obj_has_zero_dim( A ) )
  {
    FLA_Scal_external( beta, y );
    return FLA_SUCCESS;
  }

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_x    = FLA_Obj_vector_inc( x );
  inc_y    = FLA_Obj_vector_inc( y );

  FLA_Param_map_flame_to_blis_trans( transa, &blis_transa );
  FLA_Param_map_flame_to_blis_conj( conjx, &blis_conjx );


  switch( datatype ){
  
  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_x     = ( float * ) FLA_FLOAT_PTR( x );
    float *buff_y     = ( float * ) FLA_FLOAT_PTR( y );
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    bl1_sgemv( blis_transa,
               blis_conjx,
               m_A,
               n_A, 
               buff_alpha,  
               buff_A, rs_A, cs_A, 
               buff_x, inc_x,
               buff_beta,  
               buff_y, inc_y );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_x     = ( double * ) FLA_DOUBLE_PTR( x );
    double *buff_y     = ( double * ) FLA_DOUBLE_PTR( y );
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );

    bl1_dgemv( blis_transa,
               blis_conjx,
               m_A,
               n_A, 
               buff_alpha,  
               buff_A, rs_A, cs_A, 
               buff_x, inc_x,
               buff_beta,  
               buff_y, inc_y );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_x     = ( scomplex * ) FLA_COMPLEX_PTR( x );
    scomplex *buff_y     = ( scomplex * ) FLA_COMPLEX_PTR( y );
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );
    scomplex *buff_beta  = ( scomplex * ) FLA_COMPLEX_PTR( beta );

    bl1_cgemv( blis_transa,
               blis_conjx,
               m_A,
               n_A, 
               buff_alpha,  
               buff_A, rs_A, cs_A, 
               buff_x, inc_x,
               buff_beta,  
               buff_y, inc_y );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_x     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( x );
    dcomplex *buff_y     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( y );
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    dcomplex *buff_beta  = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    bl1_zgemv( blis_transa,
               blis_conjx,
               m_A,
               n_A, 
               buff_alpha,  
               buff_A, rs_A, cs_A, 
               buff_x, inc_x,
               buff_beta,  
               buff_y, inc_y );

    break;
  }

  }
  
  return FLA_SUCCESS;
}

