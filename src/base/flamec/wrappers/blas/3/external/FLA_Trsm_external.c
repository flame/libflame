
#include "FLAME.h"

FLA_Error FLA_Trsm_external( FLA_Side side, FLA_Uplo uplo, FLA_Trans trans, FLA_Diag diag, FLA_Obj alpha, FLA_Obj A, FLA_Obj B )
{
  FLA_Datatype datatype;
  int          m_B, n_B;
  int          rs_A, cs_A;
  int          rs_B, cs_B;
  side1_t       blis_side; 
  uplo1_t       blis_uplo;
  trans1_t      blis_trans;
  diag1_t       blis_diag;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING ) 
    FLA_Trsm_check( side, uplo, trans, diag, alpha, A, B );

  if ( FLA_Obj_has_zero_dim( B ) ) return FLA_SUCCESS;

  datatype = FLA_Obj_datatype( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  m_B      = FLA_Obj_length( B );
  n_B      = FLA_Obj_width( B );
  rs_B     = FLA_Obj_row_stride( B );
  cs_B     = FLA_Obj_col_stride( B );

  FLA_Param_map_flame_to_blis_side( side, &blis_side );
  FLA_Param_map_flame_to_blis_uplo( uplo, &blis_uplo );
  FLA_Param_map_flame_to_blis_trans( trans, &blis_trans );
  FLA_Param_map_flame_to_blis_diag( diag, &blis_diag );


  switch( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A     = ( float * ) FLA_FLOAT_PTR( A );
    float *buff_B     = ( float * ) FLA_FLOAT_PTR( B );
    float *buff_alpha = ( float * ) FLA_FLOAT_PTR( alpha );

    bl1_strsm( blis_side,
               blis_uplo, 
               blis_trans,
               blis_diag,
               m_B,
               n_B,
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_B, rs_B, cs_B );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A     = ( double * ) FLA_DOUBLE_PTR( A );
    double *buff_B     = ( double * ) FLA_DOUBLE_PTR( B );
    double *buff_alpha = ( double * ) FLA_DOUBLE_PTR( alpha );

    bl1_dtrsm( blis_side,
               blis_uplo, 
               blis_trans,
               blis_diag,
               m_B,
               n_B,
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_B, rs_B, cs_B );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A     = ( scomplex * ) FLA_COMPLEX_PTR( A );
    scomplex *buff_B     = ( scomplex * ) FLA_COMPLEX_PTR( B );
    scomplex *buff_alpha = ( scomplex * ) FLA_COMPLEX_PTR( alpha );

    bl1_ctrsm( blis_side,
               blis_uplo, 
               blis_trans,
               blis_diag,
               m_B,
               n_B,
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_B, rs_B, cs_B );

    break;
  }


  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );
    dcomplex *buff_B     = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( B );
    dcomplex *buff_alpha = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( alpha );

    bl1_ztrsm( blis_side,
               blis_uplo, 
               blis_trans,
               blis_diag,
               m_B,
               n_B,
               buff_alpha,
               buff_A, rs_A, cs_A, 
               buff_B, rs_B, cs_B );

    break;
  }

  }

  return FLA_SUCCESS;
}

