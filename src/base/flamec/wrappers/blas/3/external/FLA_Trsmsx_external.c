
#include "FLAME.h"

FLA_Error FLA_Trsmsx_external( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B, FLA_Obj beta, FLA_Obj C )
{
  FLA_Datatype datatype;
  int          m_B, n_B;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  int          rs_C, cs_C;
  side1_t       blis_side; 
  uplo1_t       blis_uplo;
  trans1_t      blis_trans;
  diag1_t       blis_diag;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsmsx_check( side, uplo, trans, diag, alpha, A, B, beta, C );

  if ( FLA_Obj_has_zero_dim( B ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );

  rs_C     = FLA_Obj_row_stride( C );
  cs_C     = FLA_Obj_col_stride( C );

  FLA_Param_map_flame_to_blis_side( side, &blis_side );
  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );
  FLA_Param_map_flame_to_blis_trans( trans, &blis_trans );
  FLA_Param_map_flame_to_blis_diag( diag, &blis_diag );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_B     = ( float * ) FLA_FLOAT_PTR( B );
    float *buff_C     = ( float * ) FLA_FLOAT_PTR( C );
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );
    float *buff_beta  = ( float * ) FLA_FLOAT_PTR( beta );

    bl1_strsmsx( blis_side,
                 blis_uplo, 
                 blis_trans,
                 blis_diag,
                 m_B,
                 n_B,
                 buff_alpha,
                 buff_A, rs_A, cs_A, 
                 buff_B, rs_B, cs_B, 
                 buff_beta,
                 buff_C, rs_C, cs_C );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_B     = ( double * ) FLA_DOUBLE_PTR( B );
    double *buff_C     = ( double * ) FLA_DOUBLE_PTR( C );
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );
    double *buff_beta  = ( double * ) FLA_DOUBLE_PTR( beta );

    bl1_dtrsmsx( blis_side,
                 blis_uplo, 
                 blis_trans,
                 blis_diag,
                 m_B,
                 n_B,
                 buff_alpha,
                 buff_A, rs_A, cs_A, 
                 buff_B, rs_B, cs_B, 
                 buff_beta,
                 buff_C, rs_C, cs_C );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_B     = ( scomplex * ) FLA_COMPLEX_PTR( B );
    scomplex *buff_C     = ( scomplex * ) FLA_COMPLEX_PTR( C );
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );
    scomplex *buff_beta  = ( scomplex * ) FLA_COMPLEX_PTR( beta );

    bl1_ctrsmsx( blis_side,
                 blis_uplo, 
                 blis_trans,
                 blis_diag,
                 m_B,
                 n_B,
                 buff_alpha,
                 buff_A, rs_A, cs_A, 
                 buff_B, rs_B, cs_B, 
                 buff_beta,
                 buff_C, rs_C, cs_C );

    break;
  }


  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_B     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
    dcomplex *buff_C     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( C );
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );
    dcomplex *buff_beta  = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( beta );

    bl1_ztrsmsx( blis_side,
                 blis_uplo, 
                 blis_trans,
                 blis_diag,
                 m_B,
                 n_B,
                 buff_alpha,
                 buff_A, rs_A, cs_A, 
                 buff_B, rs_B, cs_B, 
                 buff_beta,
                 buff_C, rs_C, cs_C );

    break;
  }

  }

  return FLA_SUCCESS;
}

