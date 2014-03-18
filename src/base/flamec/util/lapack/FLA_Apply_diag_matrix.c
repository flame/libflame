
#include "FLAME.h"

FLA_Error FLA_Apply_diag_matrix( FLA_Side side, FLA_Conj conj, FLA_Obj x, FLA_Obj A )
{
  FLA_Datatype dt_x, dt_A;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_x;
  side1_t       blis_side; 
  conj1_t       blis_conj;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING ) 
    FLA_Apply_diag_matrix_check( side, conj, x, A );

  dt_x     = FLA_Obj_datatype( x );
  dt_A     = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_x    = FLA_Obj_vector_inc( x );

  FLA_Param_map_flame_to_blis_side( side, &blis_side );
  FLA_Param_map_flame_to_blis_conj( conj, &blis_conj );


  switch ( dt_A )
  {
    case FLA_FLOAT:
    {
      float* buff_x  = ( float* ) FLA_FLOAT_PTR( x );
      float* buff_A  = ( float* ) FLA_FLOAT_PTR( A );

      bl1_sapdiagmv( blis_side,
                     blis_conj,
                     m_A,
                     n_A,
                     buff_x, inc_x,
                     buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double*   buff_x  = ( double* ) FLA_DOUBLE_PTR( x );
      double*   buff_A  = ( double* ) FLA_DOUBLE_PTR( A );

      bl1_dapdiagmv( blis_side,
                     blis_conj,
                     m_A,
                     n_A,
                     buff_x, inc_x,
                     buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      if ( dt_x == FLA_FLOAT )
      {
        float*    buff_x  = ( float*    ) FLA_FLOAT_PTR( x );
        scomplex* buff_A  = ( scomplex* ) FLA_COMPLEX_PTR( A );

        bl1_csapdiagmv( blis_side,
                        blis_conj,
                        m_A,
                        n_A,
                        buff_x, inc_x,
                        buff_A, rs_A, cs_A );
      }
      else if ( dt_x == FLA_COMPLEX )
      {
        scomplex* buff_x  = ( scomplex* ) FLA_COMPLEX_PTR( x );
        scomplex* buff_A  = ( scomplex* ) FLA_COMPLEX_PTR( A );

        bl1_capdiagmv( blis_side,
                       blis_conj,
                       m_A,
                       n_A,
                       buff_x, inc_x,
                       buff_A, rs_A, cs_A );
      }
      
      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      if ( dt_x == FLA_DOUBLE )
      {
        double*   buff_x  = ( double*   ) FLA_DOUBLE_PTR( x );
        dcomplex* buff_A  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );

        bl1_zdapdiagmv( blis_side,
                        blis_conj,
                        m_A,
                        n_A,
                        buff_x, inc_x,
                        buff_A, rs_A, cs_A );
      }
      else if ( dt_x == FLA_DOUBLE_COMPLEX )
      {
        dcomplex* buff_x  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( x );
        dcomplex* buff_A  = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( A );

        bl1_zapdiagmv( blis_side,
                       blis_conj,
                       m_A,
                       n_A,
                       buff_x, inc_x,
                       buff_A, rs_A, cs_A );
      }

      break;
    }
  }

  return FLA_SUCCESS;
}

