/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_realify_diagonals( FLA_Uplo uplo, FLA_Obj a, FLA_Obj b, FLA_Obj d, FLA_Obj e )
{
  FLA_Error r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Bidiag_UT_realify_diagonals_check( uplo, a, b, d, e );

  if ( uplo == FLA_LOWER_TRIANGULAR )
    r_val = FLA_Bidiag_UT_realify_diagonals_opt( a, b, d, e );
  else
    r_val = FLA_Bidiag_UT_realify_diagonals_opt( a, b, e, d );

  return r_val;
}

FLA_Error FLA_Bidiag_UT_realify_diagonals_opt( FLA_Obj a, FLA_Obj b, FLA_Obj d, FLA_Obj e ) 
{
  FLA_Datatype datatype;
  integer          i, m, inc_a, inc_b, inc_d, inc_e;

  datatype = FLA_Obj_datatype( a );

  m        = FLA_Obj_vector_dim( a );  

  inc_a    = FLA_Obj_vector_inc( a );
  inc_b    = ( m > 1 ? FLA_Obj_vector_inc( b ) : 0 );

  inc_d    = FLA_Obj_vector_inc( d );
  inc_e    = FLA_Obj_vector_inc( e );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_d = FLA_FLOAT_PTR( d );
      float* buff_e = FLA_FLOAT_PTR( e );
      float* buff_1 = FLA_FLOAT_PTR( FLA_ONE );
      
      bl1_ssetv( m, 
                 buff_1,
                 buff_d, inc_d );

      bl1_ssetv( m,
                 buff_1,
                 buff_e, inc_e );

      break;
    }
    case FLA_DOUBLE:
    {
      double* buff_d = FLA_DOUBLE_PTR( d );
      double* buff_e = FLA_DOUBLE_PTR( e );
      double* buff_1 = FLA_DOUBLE_PTR( FLA_ONE );

      bl1_dsetv( m,
                 buff_1,
                 buff_d, inc_d );

      bl1_dsetv( m,
                 buff_1,
                 buff_e, inc_e );

      break;
    }
    case FLA_COMPLEX:
    {
      scomplex* buff_a = FLA_COMPLEX_PTR( a );    
      scomplex* buff_b = ( m > 1 ? FLA_COMPLEX_PTR( b ) : NULL );    
      scomplex* buff_d = FLA_COMPLEX_PTR( d );    
      scomplex* buff_e = FLA_COMPLEX_PTR( e );    
      scomplex* buff_1 = FLA_COMPLEX_PTR( FLA_ONE );
      float*    buff_0 = FLA_FLOAT_PTR( FLA_ZERO ); 

      for ( i = 0; i < m; ++i )
      {
        scomplex* alpha1   = buff_a + (i  )*inc_a;
        scomplex* delta1   = buff_d + (i  )*inc_d;
        scomplex* epsilon1 = buff_e + (i  )*inc_e;

        scomplex  absv;

        if ( i == 0 )
        {
          *delta1 = *buff_1;
        }
        else
        {
          scomplex* beta1 = buff_b + (i-1)*inc_b;
          if ( beta1->imag == 0.0F )
            *delta1 = *buff_1;
          else
          {
            bl1_ccopys( BLIS1_CONJUGATE, beta1, delta1 );
            bl1_cabsval2( beta1, &absv );
            bl1_cinvscals( &absv, delta1 );

            bl1_cscals( delta1, beta1 );
            beta1->imag = *buff_0;

            bl1_cscals( delta1, alpha1 );
          }
        }

        if ( alpha1->imag == 0.0F )
          *epsilon1 = *buff_1;          
        else
        {
          bl1_ccopys( BLIS1_CONJUGATE, alpha1, epsilon1 );
          bl1_cabsval2( alpha1, &absv );
          bl1_cinvscals( &absv, epsilon1 );
          
          bl1_cscals( epsilon1, alpha1 );
          alpha1->imag = *buff_0;
        }

        if ( i < ( m - 1 ) )
        {
          scomplex* beta2 = buff_b + (i )*inc_b;
          bl1_cscals( epsilon1, beta2 );
        }
      }
      break;
    }
    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_a = FLA_DOUBLE_COMPLEX_PTR( a );
      dcomplex* buff_b = ( m > 1 ? FLA_DOUBLE_COMPLEX_PTR( b ) : NULL );
      dcomplex* buff_d = FLA_DOUBLE_COMPLEX_PTR( d );
      dcomplex* buff_e = FLA_DOUBLE_COMPLEX_PTR( e );
      dcomplex* buff_1 = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
      double*   buff_0 = FLA_DOUBLE_PTR( FLA_ZERO );

      for ( i = 0; i < m; ++i )
      {
        dcomplex* alpha1   = buff_a + (i  )*inc_a;
        dcomplex* delta1   = buff_d + (i  )*inc_d;
        dcomplex* epsilon1 = buff_e + (i  )*inc_e;

        dcomplex  absv;

        if ( i == 0 )
        {
          *delta1 = *buff_1;
        }
        else
        {
          dcomplex* beta1    = buff_b + (i-1)*inc_b;
          bl1_zcopys( BLIS1_CONJUGATE, beta1, delta1 );
          bl1_zabsval2( beta1, &absv );
          bl1_zinvscals( &absv, delta1 );

          bl1_zscals( delta1, beta1 );
          beta1->imag = *buff_0;

          bl1_zscals( delta1, alpha1 );
        }

        bl1_zcopys( BLIS1_CONJUGATE, alpha1, epsilon1 );
        bl1_zabsval2( alpha1, &absv );
        bl1_zinvscals( &absv, epsilon1 );

        bl1_zscals( epsilon1, alpha1 );
        alpha1->imag = *buff_0;

        if ( i < ( m - 1 ) )
        {
          dcomplex* beta2 = buff_b + (i  )*inc_b;
          bl1_zscals( epsilon1, beta2 );
        }
      }
      break;
    }
  }
  return FLA_SUCCESS;
}

