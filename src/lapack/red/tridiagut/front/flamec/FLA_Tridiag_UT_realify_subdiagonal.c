
#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_realify_subdiagonal( FLA_Obj b, FLA_Obj d )
{
  FLA_Error r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Tridiag_UT_realify_subdiagonal_check( b, d );

  r_val = FLA_Tridiag_UT_realify_subdiagonal_opt( b, d );

  return r_val;
}

FLA_Error FLA_Tridiag_UT_realify_subdiagonal_opt( FLA_Obj b, FLA_Obj d )
{
  FLA_Datatype datatype;
  int          m, inc_b, inc_d;
  int          i;

  datatype = FLA_Obj_datatype( d );

  m        = FLA_Obj_vector_dim( d );
  inc_d    = FLA_Obj_vector_inc( d );

  inc_b    = ( m > 1 ? FLA_Obj_vector_inc( b ) : 0 );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_d = FLA_FLOAT_PTR( d );
      float* buff_1 = FLA_FLOAT_PTR( FLA_ONE );

      bl1_ssetv( m,
                 buff_1,
                 buff_d, inc_d );

      break;
    }
    case FLA_DOUBLE:
    {
      double* buff_d = FLA_DOUBLE_PTR( d );
      double* buff_1 = FLA_DOUBLE_PTR( FLA_ONE );

      bl1_dsetv( m,
                 buff_1,
                 buff_d, inc_d );

      break;
    }
    case FLA_COMPLEX:
    {
      scomplex* buff_b = ( m > 1 ? FLA_COMPLEX_PTR( b ) : NULL );
      scomplex* buff_d = FLA_COMPLEX_PTR( d );
      scomplex* buff_1 = FLA_COMPLEX_PTR( FLA_ONE );

      bl1_csetv( 1,
                 buff_1,
                 buff_d, inc_d );

      for ( i = 1; i < m; ++i )
      {
        scomplex* beta1  = buff_b + (i-1)*inc_b;
        scomplex* delta1 = buff_d + (i  )*inc_d;
        scomplex  absv;
        scomplex  conj_delta1;

        if ( beta1->imag == 0.0F )
          *delta1 = *buff_1;
        else
        {
          bl1_ccopys( BLIS1_CONJUGATE, beta1, delta1 );
          bl1_cabsval2( beta1, &absv );
          bl1_cinvscals( &absv, delta1 );
          
          *beta1 = absv;
        }
        if ( i < ( m - 1 ) )
        {
          scomplex* beta2 = buff_b + (i  )*inc_b;
          bl1_ccopyconj( delta1, &conj_delta1 );
          bl1_cscals( &conj_delta1, beta2 );
        }
      }
      break;
    }
    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_b = ( m > 1 ? FLA_DOUBLE_COMPLEX_PTR( b ) : NULL );
      dcomplex* buff_d = FLA_DOUBLE_COMPLEX_PTR( d );
      dcomplex* buff_1 = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );

      bl1_zsetv( 1,
                 buff_1,
                 buff_d, inc_d );

      for ( i = 1; i < m; ++i )
      {
        dcomplex* beta1  = buff_b + (i-1)*inc_b;
        dcomplex* delta1 = buff_d + (i  )*inc_d;
        dcomplex  absv;
        dcomplex  conj_delta1;

        if ( beta1->imag == 0.0 )
          *delta1 = *buff_1;
        else
        {
          bl1_zcopys( BLIS1_CONJUGATE, beta1, delta1 );
          bl1_zabsval2( beta1, &absv );
          bl1_zinvscals( &absv, delta1 );

          *beta1 = absv;
        }
        if ( i < ( m - 1 ) )
        {
          dcomplex* beta2 = buff_b + (i  )*inc_b;
          bl1_zcopyconj( delta1, &conj_delta1 );
          bl1_zscals( &conj_delta1, beta2 );
        }
      }
      break;
    }
  }

  return FLA_SUCCESS;
}




