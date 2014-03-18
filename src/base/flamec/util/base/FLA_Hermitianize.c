
#include "FLAME.h"

FLA_Error FLA_Hermitianize( FLA_Uplo uplo, FLA_Obj A )
{
  FLA_Datatype datatype;
  dim_t        n_A;
  dim_t        rs_A, cs_A;
  conj1_t       blis_conj;
  uplo1_t       blis_uplo;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Hermitianize_check( uplo, A );

  datatype = FLA_Obj_datatype( A );

  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  FLA_Param_map_flame_to_blis_conj( FLA_CONJUGATE, &blis_conj );
  FLA_Param_map_flame_to_blis_uplo( uplo,          &blis_uplo );


  switch ( datatype ){

  case FLA_FLOAT:
  {
    float *buff_A = ( float * ) FLA_FLOAT_PTR( A );

    bl1_ssymmize( blis_conj,
                  blis_uplo,
                  n_A,
                  buff_A, rs_A, cs_A );

    break;
  }

  case FLA_DOUBLE:
  {
    double *buff_A = ( double * ) FLA_DOUBLE_PTR( A );

    bl1_dsymmize( blis_conj,
                  blis_uplo,
                  n_A,
                  buff_A, rs_A, cs_A );

    break;
  }

  case FLA_COMPLEX:
  {
    scomplex *buff_A = ( scomplex * ) FLA_COMPLEX_PTR( A );

    bl1_csymmize( blis_conj,
                  blis_uplo,
                  n_A,
                  buff_A, rs_A, cs_A );

    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex *buff_A = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( A );

    bl1_zsymmize( blis_conj,
                  blis_uplo,
                  n_A,
                  buff_A, rs_A, cs_A );

    break;
  }

  }

  return FLA_SUCCESS;
}

