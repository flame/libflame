/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define ssign( x ) ( (x) < 0.0F ? -1.0F : 1.0F )
#define dsign( x ) ( (x) < 0.0  ? -1.0  : 1.0  )

FLA_Error FLA_Househ3UD_UT( FLA_Obj chi_0, FLA_Obj x1, FLA_Obj y2, FLA_Obj tau )
/*
  Compute an up-and-downdating UT Householder transformation

          / / 1 0 0 \            / 1 0  0 \  / 1  \ ( 1  u1' v2' ) \
    H  =  | | 0 I 0 | - inv(tau) | 0 I  0 |  | u1 |                |
          \ \ 0 0 I /            \ 0 0 -I /  \ v2 /                /

  by computing tau, u1, and v2 such that the following is satisfied:

      / chi_0 \   / alpha \
    H |  x1   | = |   0   |
      \  y2   /   \   0   /

  where

    alpha  = - lambda * chi_0 / | chi_0 |

    lambda =  sqrt( conj(chi0) chi0 + x1' x1 - y2' y2 )

              / chi_0 \
    x      =  |  x1   |
              \  y2   /

    tau    =  ( 1 + u1' u1 - v2' v2 ) / 2

    u1     =  x1 / ( chi_0 - alpha )

    v2     =  -y2 / ( chi_0 - alpha )

  Upon completion, alpha, u1, and v2 have overwritten objects chi_0, x1,
  and y2, respectively.

  -FGVZ
*/
{
  FLA_Datatype datatype;
  integer          m_x1;
  integer          m_y2;
  integer          inc_x1;
  integer          inc_y2;

  datatype = FLA_Obj_datatype( x1 );

  m_x1     = FLA_Obj_vector_dim( x1 );
  m_y2     = FLA_Obj_vector_dim( y2 );
  inc_x1   = FLA_Obj_vector_inc( x1 );
  inc_y2   = FLA_Obj_vector_inc( y2 );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Househ3UD_UT_check( chi_0, x1, y2, tau );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* chi_0_p = ( float* ) FLA_FLOAT_PTR( chi_0 );
      float* x1_p    = ( float* ) FLA_FLOAT_PTR( x1 );
      float* y2_p    = ( float* ) FLA_FLOAT_PTR( y2 );
      float* tau_p   = ( float* ) FLA_FLOAT_PTR( tau );

      FLA_Househ3UD_UT_ops( m_x1,
                            m_y2,
                            chi_0_p,
                            x1_p, inc_x1,
                            y2_p, inc_y2,
                            tau_p );
      break;
    }

    case FLA_DOUBLE:
    {
      double* chi_0_p = ( double* ) FLA_DOUBLE_PTR( chi_0 );
      double* x1_p    = ( double* ) FLA_DOUBLE_PTR( x1 );
      double* y2_p    = ( double* ) FLA_DOUBLE_PTR( y2 );
      double* tau_p   = ( double* ) FLA_DOUBLE_PTR( tau );

      FLA_Househ3UD_UT_opd( m_x1,
                            m_y2,
                            chi_0_p,
                            x1_p, inc_x1,
                            y2_p, inc_y2,
                            tau_p );
      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* chi_0_p = ( scomplex* ) FLA_COMPLEX_PTR( chi_0 );
      scomplex* x1_p    = ( scomplex* ) FLA_COMPLEX_PTR( x1 );
      scomplex* y2_p    = ( scomplex* ) FLA_COMPLEX_PTR( y2 );
      scomplex* tau_p   = ( scomplex* ) FLA_COMPLEX_PTR( tau );

      FLA_Househ3UD_UT_opc( m_x1,
                            m_y2,
                            chi_0_p,
                            x1_p, inc_x1,
                            y2_p, inc_y2,
                            tau_p );
      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* chi_0_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( chi_0 );
      dcomplex* x1_p    = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( x1 );
      dcomplex* y2_p    = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( y2 );
      dcomplex* tau_p   = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( tau );

      FLA_Househ3UD_UT_opz( m_x1,
                            m_y2,
                            chi_0_p,
                            x1_p, inc_x1,
                            y2_p, inc_y2,
                            tau_p );
      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Househ3UD_UT_ops( integer       m_x1,
                                integer       m_y2,
                                float*    chi_0,
                                float*    x1, integer inc_x1,
                                float*    y2, integer inc_y2,
                                float*    tau )
{
  float    one_half = *FLA_FLOAT_PTR( FLA_ONE_HALF );
  float    alpha;
  float    chi_0_minus_alpha;
  float    neg_chi_0_minus_alpha;
  float    abs_chi_0;
  float    norm_x_1;
  float    norm_y_2;
  float    lambda;
  float    abs_sq_chi_0_minus_alpha;
  integer      i_one = 1;

  //
  // Compute the 2-norms of x_1 and y_2:
  //
  //   norm_x_1 := || x_1 ||_2
  //   norm_y_2 := || y_2 ||_2
  //

  bl1_snrm2( m_x1,
             x1, inc_x1,
             &norm_x_1 );

  bl1_snrm2( m_y2,
             y2, inc_y2,
             &norm_y_2 );

  //
  // If 2-norms of x_1, y_2 are zero, then return with trivial tau, chi_0 values.
  //

  if ( norm_x_1 == 0.0F && 
       norm_y_2 == 0.0F )
  {
    *chi_0 = -(*chi_0);
    *tau   = one_half;

    return FLA_SUCCESS;
  }

  //
  // Compute the absolute value (magnitude) of chi_0, which equals the 2-norm
  // of chi_0:
  //
  //   abs_chi_0 :=  | chi_0 |  =  || chi_0 ||_2
  //

  bl1_snrm2( i_one,
             chi_0, i_one,
             &abs_chi_0 );

  //
  // Compute lambda:
  //
  //   lambda := sqrt( conj(chi0) chi0 + x1' x1 - y2' y2 )
  //

  lambda = ( float ) sqrt( abs_chi_0 * abs_chi_0 + 
                           norm_x_1  * norm_x_1  -
                           norm_y_2  * norm_y_2 );

  // Compute alpha:
  //
  //   alpha := - lambda * chi_0 / | chi_0 |
  //          = -sign( chi_0 ) * lambda
  //

  alpha = -ssign( *chi_0 ) * lambda;


  //
  // Overwrite x_1 and y_2 with u_1 and v_2, respectively:
  //
  //   x_1 := x_1 / ( chi_0 - alpha )
  //   y_2 := y_2 / -( chi_0 - alpha )
  //

  chi_0_minus_alpha = (*chi_0) - alpha;

  bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                 m_x1,
                 &chi_0_minus_alpha,
                 x1, inc_x1 );

  neg_chi_0_minus_alpha = -chi_0_minus_alpha;

  bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                 m_y2,
                 &neg_chi_0_minus_alpha,
                 y2, inc_y2 );

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_1' * u_1 - v_2' * v_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_1' * x_1 - y_2' * y_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = ( | chi_1 - alpha |^2 + || x_2 ||_2^2 - || y_2 ||_2^2 ) /
  //          ( 2 * | chi_1 - alpha |^2 )
  //

  abs_sq_chi_0_minus_alpha = chi_0_minus_alpha * chi_0_minus_alpha;

  *tau = ( abs_sq_chi_0_minus_alpha +
           norm_x_1 * norm_x_1 -
           norm_y_2 * norm_y_2 ) /
         ( 2.0F * abs_sq_chi_0_minus_alpha );

  //
  // Overwrite chi_0 with alpha:
  //
  //   chi_0 := alpha
  //

  *chi_0 = alpha;

  return FLA_SUCCESS;
}



FLA_Error FLA_Househ3UD_UT_opd( integer       m_x1,
                                integer       m_y2,
                                double*   chi_0,
                                double*   x1, integer inc_x1,
                                double*   y2, integer inc_y2,
                                double*   tau )
{
  double   one_half = *FLA_DOUBLE_PTR( FLA_ONE_HALF );
  double   alpha;
  double   chi_0_minus_alpha;
  double   neg_chi_0_minus_alpha;
  double   abs_chi_0;
  double   norm_x_1;
  double   norm_y_2;
  double   lambda;
  double   abs_sq_chi_0_minus_alpha;
  integer      i_one = 1;

  //
  // Compute the 2-norms of x_1 and y_2:
  //
  //   norm_x_1 := || x_1 ||_2
  //   norm_y_2 := || y_2 ||_2
  //

  bl1_dnrm2( m_x1,
             x1, inc_x1,
             &norm_x_1 );

  bl1_dnrm2( m_y2,
             y2, inc_y2,
             &norm_y_2 );

  //
  // If 2-norms of x_1, y_2 are zero, then return with trivial tau, chi_0 values.
  //

  if ( norm_x_1 == 0.0 && 
       norm_y_2 == 0.0 )
  {
    *chi_0 = -(*chi_0);
    *tau   = one_half;

    return FLA_SUCCESS;
  }

  //
  // Compute the absolute value (magnitude) of chi_0, which equals the 2-norm
  // of chi_0:
  //
  //   abs_chi_0 :=  | chi_0 |  =  || chi_0 ||_2
  //

  bl1_dnrm2( i_one,
             chi_0, i_one,
             &abs_chi_0 );

  //
  // Compute lambda:
  //
  //   lambda := sqrt( conj(chi0) chi0 + x1' x1 - y2' y2 )
  //

  lambda = sqrt( abs_chi_0 * abs_chi_0 + 
                 norm_x_1  * norm_x_1  -
                 norm_y_2  * norm_y_2 );

  // Compute alpha:
  //
  //   alpha := - lambda * chi_0 / | chi_0 |
  //          = -sign( chi_0 ) * lambda
  //

  alpha = -dsign( *chi_0 ) * lambda;

  //
  // Overwrite x_1 and y_2 with u_1 and v_2, respectively:
  //
  //   x_1 := x_1 / ( chi_0 - alpha )
  //   y_2 := y_2 / -( chi_0 - alpha )
  //

  chi_0_minus_alpha = (*chi_0) - alpha;

  bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                 m_x1,
                 &chi_0_minus_alpha,
                 x1, inc_x1 );

  neg_chi_0_minus_alpha = -chi_0_minus_alpha;

  bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                 m_y2,
                 &neg_chi_0_minus_alpha,
                 y2, inc_y2 );

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_1' * u_1 - v_2' * v_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_1' * x_1 - y_2' * y_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = ( | chi_1 - alpha |^2 + || x_2 ||_2^2 - || y_2 ||_2^2 ) /
  //          ( 2 * | chi_1 - alpha |^2 )
  //

  abs_sq_chi_0_minus_alpha = chi_0_minus_alpha * chi_0_minus_alpha;

  *tau = ( abs_sq_chi_0_minus_alpha +
           norm_x_1 * norm_x_1 -
           norm_y_2 * norm_y_2 ) /
         ( 2.0 * abs_sq_chi_0_minus_alpha );

  //
  // Overwrite chi_0 with alpha:
  //
  //   chi_0 := alpha
  //

  *chi_0 = alpha;

  return FLA_SUCCESS;
}



FLA_Error FLA_Househ3UD_UT_opc( integer       m_x1,
                                integer       m_y2,
                                scomplex* chi_0,
                                scomplex* x1, integer inc_x1,
                                scomplex* y2, integer inc_y2,
                                scomplex* tau )
{
  scomplex one_half = *FLA_COMPLEX_PTR( FLA_ONE_HALF );
  scomplex alpha;
  scomplex chi_0_minus_alpha;
  scomplex neg_chi_0_minus_alpha;
  float    abs_chi_0;
  float    norm_x_1;
  float    norm_y_2;
  float    lambda;
  float    abs_sq_chi_0_minus_alpha;
  integer      i_one = 1;

  //
  // Compute the 2-norms of x_1 and y_2:
  //
  //   norm_x_1 := || x_1 ||_2
  //   norm_y_2 := || y_2 ||_2
  //

  bl1_cnrm2( m_x1,
             x1, inc_x1,
             &norm_x_1 );

  bl1_cnrm2( m_y2,
             y2, inc_y2,
             &norm_y_2 );

  //
  // If 2-norms of x_1, y_2 are zero, then return with trivial tau, chi_0 values.
  //

  if ( norm_x_1 == 0.0F && 
       norm_y_2 == 0.0F )
  {
    chi_0->real = -(chi_0->real);
    chi_0->imag = -(chi_0->imag);
    tau->real   = one_half.real;
    tau->imag   = one_half.imag;

    return FLA_SUCCESS;
  }

  //
  // Compute the absolute value (magnitude) of chi_0, which equals the 2-norm
  // of chi_0:
  //
  //   abs_chi_0 :=  | chi_0 |  =  || chi_0 ||_2
  //

  bl1_cnrm2( i_one,
             chi_0, i_one,
             &abs_chi_0 );

  //
  // Compute lambda:
  //
  //   lambda := sqrt( conj(chi0) chi0 + x1' x1 - y2' y2 )
  //

  lambda = ( float ) sqrt( abs_chi_0 * abs_chi_0 + 
                           norm_x_1  * norm_x_1  -
                           norm_y_2  * norm_y_2 );
  
  //
  // Compute alpha:
  //
  //   alpha := - lambda * chi_0 / | chi_0 |
  //

  alpha.real = -chi_0->real * lambda / abs_chi_0;
  alpha.imag = -chi_0->imag * lambda / abs_chi_0;

  //
  // Overwrite x_1 and y_2 with u_1 and v_2, respectively:
  //
  //   x_1 := x_1 / ( chi_0 - alpha )
  //   y_2 := y_2 / -( chi_0 - alpha )
  //

  chi_0_minus_alpha.real = chi_0->real - alpha.real;
  chi_0_minus_alpha.imag = chi_0->imag - alpha.imag;

  bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                 m_x1,
                 &chi_0_minus_alpha,
                 x1, inc_x1 );

  neg_chi_0_minus_alpha.real = -chi_0_minus_alpha.real;
  neg_chi_0_minus_alpha.imag = -chi_0_minus_alpha.imag;

  bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                 m_y2,
                 &neg_chi_0_minus_alpha,
                 y2, inc_y2 );

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_1' * u_1 - v_2' * v_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_1' * x_1 - y_2' * y_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = ( | chi_1 - alpha |^2 + || x_2 ||_2^2 - || y_2 ||_2^2 ) /
  //          ( 2 * | chi_1 - alpha |^2 )
  //

  abs_sq_chi_0_minus_alpha = chi_0_minus_alpha.real * chi_0_minus_alpha.real +
                             chi_0_minus_alpha.imag * chi_0_minus_alpha.imag;

  tau->real = ( abs_sq_chi_0_minus_alpha +
                norm_x_1 * norm_x_1 -
                norm_y_2 * norm_y_2 ) /
              ( 2.0F * abs_sq_chi_0_minus_alpha );
  tau->imag = 0.0F;

  //
  // Overwrite chi_0 with alpha:
  //
  //   chi_0 := alpha
  //

  chi_0->real = alpha.real;
  chi_0->imag = alpha.imag;

  return FLA_SUCCESS;
}



FLA_Error FLA_Househ3UD_UT_opz( integer       m_x1,
                                integer       m_y2,
                                dcomplex* chi_0,
                                dcomplex* x1, integer inc_x1,
                                dcomplex* y2, integer inc_y2,
                                dcomplex* tau )
{
  dcomplex one_half = *FLA_DOUBLE_COMPLEX_PTR( FLA_ONE_HALF );
  dcomplex alpha;
  dcomplex chi_0_minus_alpha;
  dcomplex neg_chi_0_minus_alpha;
  double   abs_chi_0;
  double   norm_x_1;
  double   norm_y_2;
  double   lambda;
  double   abs_sq_chi_0_minus_alpha;
  integer      i_one = 1;

  //
  // Compute the 2-norms of x_1 and y_2:
  //
  //   norm_x_1 := || x_1 ||_2
  //   norm_y_2 := || y_2 ||_2
  //

  bl1_znrm2( m_x1,
             x1, inc_x1,
             &norm_x_1 );

  bl1_znrm2( m_y2,
             y2, inc_y2,
             &norm_y_2 );

  //
  // If 2-norms of x_1, y_2 are zero, then return with trivial tau, chi_0 values.
  //

  if ( norm_x_1 == 0.0 && 
       norm_y_2 == 0.0 )
  {
    chi_0->real = -(chi_0->real);
    chi_0->imag = -(chi_0->imag);
    tau->real   = one_half.real;
    tau->imag   = one_half.imag;

    return FLA_SUCCESS;
  }

  //
  // Compute the absolute value (magnitude) of chi_0, which equals the 2-norm
  // of chi_0:
  //
  //   abs_chi_0 :=  | chi_0 |  =  || chi_0 ||_2
  //

  bl1_znrm2( i_one,
             chi_0, i_one,
             &abs_chi_0 );

  //
  // Compute lambda:
  //
  //   lambda := sqrt( conj(chi0) chi0 + x1' x1 - y2' y2 )
  //

  lambda = sqrt( abs_chi_0 * abs_chi_0 + 
                 norm_x_1  * norm_x_1  -
                 norm_y_2  * norm_y_2 );
  
  //
  // Compute alpha:
  //
  //   alpha := - lambda * chi_0 / | chi_0 |
  //

  alpha.real = -chi_0->real * lambda / abs_chi_0;
  alpha.imag = -chi_0->imag * lambda / abs_chi_0;

  //
  // Overwrite x_1 and y_2 with u_1 and v_2, respectively:
  //
  //   x_1 := x_1 / ( chi_0 - alpha )
  //   y_2 := y_2 / -( chi_0 - alpha )
  //

  chi_0_minus_alpha.real = chi_0->real - alpha.real;
  chi_0_minus_alpha.imag = chi_0->imag - alpha.imag;

  bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                 m_x1,
                 &chi_0_minus_alpha,
                 x1, inc_x1 );

  neg_chi_0_minus_alpha.real = -chi_0_minus_alpha.real;
  neg_chi_0_minus_alpha.imag = -chi_0_minus_alpha.imag;

  bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                 m_y2,
                 &neg_chi_0_minus_alpha,
                 y2, inc_y2 );

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_1' * u_1 - v_2' * v_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_1' * x_1 - y_2' * y_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = ( | chi_1 - alpha |^2 + || x_2 ||_2^2 - || y_2 ||_2^2 ) /
  //          ( 2 * | chi_1 - alpha |^2 )
  //

  abs_sq_chi_0_minus_alpha = chi_0_minus_alpha.real * chi_0_minus_alpha.real +
                             chi_0_minus_alpha.imag * chi_0_minus_alpha.imag;

  tau->real = ( abs_sq_chi_0_minus_alpha +
                norm_x_1 * norm_x_1 -
                norm_y_2 * norm_y_2 ) /
              ( 2.0 * abs_sq_chi_0_minus_alpha );
  tau->imag = 0.0;

  //
  // Overwrite chi_0 with alpha:
  //
  //   chi_0 := alpha
  //

  chi_0->real = alpha.real;
  chi_0->imag = alpha.imag;

  return FLA_SUCCESS;
}

