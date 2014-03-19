/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#define ssign( x ) ( (x) < 0.0F ? -1.0F : 1.0F )
#define dsign( x ) ( (x) < 0.0  ? -1.0  : 1.0  )

FLA_Error FLA_Househ2s_UT( FLA_Side side, FLA_Obj chi_1, FLA_Obj x2, FLA_Obj alpha, FLA_Obj chi_1_minus_alpha, FLA_Obj tau )
{
  FLA_Datatype datatype;
  int          m_x2;
  int          inc_x2;

  datatype = FLA_Obj_datatype( x2 );

  m_x2     = FLA_Obj_vector_dim( x2 );
  inc_x2   = FLA_Obj_vector_inc( x2 );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Househ2s_UT_check( side, chi_1, x2, alpha, chi_1_minus_alpha, tau );

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* chi_1_p             = ( float* ) FLA_FLOAT_PTR( chi_1 );
      float* x2_p                = ( float* ) FLA_FLOAT_PTR( x2 );
      float* alpha_p             = ( float* ) FLA_FLOAT_PTR( alpha );
      float* chi_1_minus_alpha_p = ( float* ) FLA_FLOAT_PTR( chi_1_minus_alpha );
      float* tau_p               = ( float* ) FLA_FLOAT_PTR( tau );

      if ( side == FLA_LEFT )
        FLA_Househ2s_UT_l_ops( m_x2,
                               chi_1_p,
                               x2_p, inc_x2,
                               alpha_p,
                               chi_1_minus_alpha_p,
                               tau_p );
      else // if ( side == FLA_RIGHT )
        FLA_Househ2s_UT_r_ops( m_x2,
                               chi_1_p,
                               x2_p, inc_x2,
                               alpha_p,
                               chi_1_minus_alpha_p,
                               tau_p );

      break;
    }

    case FLA_DOUBLE:
    {
      double* chi_1_p             = ( double* ) FLA_DOUBLE_PTR( chi_1 );
      double* x2_p                = ( double* ) FLA_DOUBLE_PTR( x2 );
      double* alpha_p             = ( double* ) FLA_DOUBLE_PTR( alpha );
      double* chi_1_minus_alpha_p = ( double* ) FLA_DOUBLE_PTR( chi_1_minus_alpha );
      double* tau_p               = ( double* ) FLA_DOUBLE_PTR( tau );

      if ( side == FLA_LEFT )
        FLA_Househ2s_UT_l_opd( m_x2,
                               chi_1_p,
                               x2_p, inc_x2,
                               alpha_p,
                               chi_1_minus_alpha_p,
                               tau_p );
      else // if ( side == FLA_RIGHT )
        FLA_Househ2s_UT_r_opd( m_x2,
                               chi_1_p,
                               x2_p, inc_x2,
                               alpha_p,
                               chi_1_minus_alpha_p,
                               tau_p );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* chi_1_p             = ( scomplex* ) FLA_COMPLEX_PTR( chi_1 );
      scomplex* x2_p                = ( scomplex* ) FLA_COMPLEX_PTR( x2 );
      scomplex* alpha_p             = ( scomplex* ) FLA_COMPLEX_PTR( alpha );
      scomplex* chi_1_minus_alpha_p = ( scomplex* ) FLA_COMPLEX_PTR( chi_1_minus_alpha );
      scomplex* tau_p               = ( scomplex* ) FLA_COMPLEX_PTR( tau );

      if ( side == FLA_LEFT )
        FLA_Househ2s_UT_l_opc( m_x2,
                               chi_1_p,
                               x2_p, inc_x2,
                               alpha_p,
                               chi_1_minus_alpha_p,
                               tau_p );
      else // if ( side == FLA_RIGHT )
        FLA_Househ2s_UT_r_opc( m_x2,
                               chi_1_p,
                               x2_p, inc_x2,
                               alpha_p,
                               chi_1_minus_alpha_p,
                               tau_p );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* chi_1_p             = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( chi_1 );
      dcomplex* x2_p                = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( x2 );
      dcomplex* alpha_p             = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( alpha );
      dcomplex* chi_1_minus_alpha_p = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( chi_1_minus_alpha );
      dcomplex* tau_p               = ( dcomplex* ) FLA_DOUBLE_COMPLEX_PTR( tau );

      if ( side == FLA_LEFT )
        FLA_Househ2s_UT_l_opz( m_x2,
                               chi_1_p,
                               x2_p, inc_x2,
                               alpha_p,
                               chi_1_minus_alpha_p,
                               tau_p );
      else // if ( side == FLA_RIGHT )
        FLA_Househ2s_UT_r_opz( m_x2,
                               chi_1_p,
                               x2_p, inc_x2,
                               alpha_p,
                               chi_1_minus_alpha_p,
                               tau_p );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Househ2s_UT_l_ops( int    m_x2,
                                 float* chi_1,
                                 float* x2, int inc_x2,
                                 float* alpha,
                                 float* chi_1_minus_alpha,
                                 float* tau )
{
  float    one_half = *FLA_FLOAT_PTR( FLA_ONE_HALF );
  float    y[2];
  float    abs_chi_1;
  float    norm_x_2;
  float    norm_x;
  float    abs_chi_1_minus_alpha;
  float    norm_x_2_div_abs_chi_1_minus_alpha;
  int      i_one = 1;
  int      i_two = 2;

  //
  // Compute the 2-norm of x_2:
  //
  //   norm_x_2 := || x_2 ||_2
  //

  bl1_snrm2( m_x2,
             x2, inc_x2,
             &norm_x_2 );

  //
  // If 2-norm of x_2 is zero, then return with trivial values.
  //

  if ( norm_x_2 == 0.0F )
  {
    *alpha             = -(*chi_1);
    *chi_1_minus_alpha = 2.0F * (*chi_1);
    *tau               = one_half;

    return FLA_SUCCESS;
  }

  //
  // Compute the absolute value (magnitude) of chi_1, which equals the 2-norm
  // of chi_1:
  //
  //   abs_chi_1 :=  | chi_1 |  =  || chi_1 ||_2
  //

  bl1_snrm2( i_one,
             chi_1, i_one,
             &abs_chi_1 );

  //
  // Compute the 2-norm of x via the two norms previously computed above:
  //
  //   norm_x :=  || x ||_2  =  || / chi_1 \ ||   =  || / || chi_1 ||_2 \ ||
  //                            || \  x_2  / ||_2    || \  || x_2 ||_2  / ||_2
  //

  y[0] = abs_chi_1;
  y[1] = norm_x_2;

  bl1_snrm2( i_two,
             y, i_one,
             &norm_x );

  //
  // Compute alpha:
  //
  //   alpha := - || x ||_2 * chi_1 / | chi_1 |
  //          = -sign( chi_1 ) * || x ||_2
  //

  *alpha = -ssign( *chi_1 ) * norm_x;

  *chi_1_minus_alpha = (*chi_1) - (*alpha);

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_2' * u_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_2' * x_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = 1/2 + ( || x ||_2 / | chi_1 - alpha | )^2
  //                                                                                                         
  bl1_sabsval2( chi_1_minus_alpha, &abs_chi_1_minus_alpha );

  norm_x_2_div_abs_chi_1_minus_alpha = norm_x_2 / abs_chi_1_minus_alpha;
  *tau = one_half + one_half*( norm_x_2_div_abs_chi_1_minus_alpha *
                               norm_x_2_div_abs_chi_1_minus_alpha );

  return FLA_SUCCESS;
}



FLA_Error FLA_Househ2s_UT_l_opd( int     m_x2,
                                 double* chi_1,
                                 double* x2, int inc_x2,
                                 double* alpha,
                                 double* chi_1_minus_alpha,
                                 double* tau )
{
  double   one_half = *FLA_DOUBLE_PTR( FLA_ONE_HALF );
  double   y[2];
  double   abs_chi_1;
  double   norm_x_2;
  double   norm_x;
  double   abs_chi_1_minus_alpha;
  double   norm_x_2_div_abs_chi_1_minus_alpha;
  int      i_one = 1;
  int      i_two = 2;

  //
  // Compute the 2-norm of x_2:
  //
  //   norm_x_2 := || x_2 ||_2
  //

  bl1_dnrm2( m_x2,
             x2, inc_x2,
             &norm_x_2 );

  //
  // If 2-norm of x_2 is zero, then return with trivial values.
  //

  if ( norm_x_2 == 0.0 )
  {
    *alpha             = -(*chi_1);
    *chi_1_minus_alpha = 2.0 * (*chi_1);
    *tau               = one_half;

    return FLA_SUCCESS;
  }

  //
  // Compute the absolute value (magnitude) of chi_1, which equals the 2-norm
  // of chi_1:
  //
  //   abs_chi_1 :=  | chi_1 |  =  || chi_1 ||_2
  //

  bl1_dnrm2( i_one,
             chi_1, i_one,
             &abs_chi_1 );

  //
  // Compute the 2-norm of x via the two norms previously computed above:
  //
  //   norm_x :=  || x ||_2  =  || / chi_1 \ ||   =  || / || chi_1 ||_2 \ ||
  //                            || \  x_2  / ||_2    || \  || x_2 ||_2  / ||_2
  //

  y[0] = abs_chi_1;
  y[1] = norm_x_2;

  bl1_dnrm2( i_two,
             y, i_one,
             &norm_x );

  //
  // Compute alpha:
  //
  //   alpha := - || x ||_2 * chi_1 / | chi_1 |
  //          = -sign( chi_1 ) * || x ||_2
  //

  *alpha = -dsign( *chi_1 ) * norm_x;

  *chi_1_minus_alpha = (*chi_1) - (*alpha);

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_2' * u_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_2' * x_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = 1/2 + ( || x ||_2 / | chi_1 - alpha | )^2
  //                                                                                                         
  bl1_dabsval2( chi_1_minus_alpha, &abs_chi_1_minus_alpha );

  norm_x_2_div_abs_chi_1_minus_alpha = norm_x_2 / abs_chi_1_minus_alpha;
  *tau = one_half + one_half*( norm_x_2_div_abs_chi_1_minus_alpha *
                               norm_x_2_div_abs_chi_1_minus_alpha );

  return FLA_SUCCESS;
}



FLA_Error FLA_Househ2s_UT_l_opc( int       m_x2,
                                 scomplex* chi_1,
                                 scomplex* x2, int inc_x2,
                                 scomplex* alpha,
                                 scomplex* chi_1_minus_alpha,
                                 scomplex* tau )
{
  scomplex one_half = *FLA_COMPLEX_PTR( FLA_ONE_HALF );
  scomplex y[2];
  float    abs_chi_1;
  float    norm_x_2;
  float    norm_x;
  float    abs_chi_1_minus_alpha;
  float    norm_x_2_div_abs_chi_1_minus_alpha;
  int      i_one = 1;
  int      i_two = 2;

  //
  // Compute the 2-norm of x_2:
  //
  //   norm_x_2 := || x_2 ||_2
  //

  bl1_cnrm2( m_x2,
             x2, inc_x2,
             &norm_x_2 );

  //
  // If 2-norm of x_2 is zero, then return with trivial values.
  //

  if ( norm_x_2 == 0.0F )
  {
    alpha->real             = -(chi_1->real);
    alpha->imag             = -(chi_1->imag);
    chi_1_minus_alpha->real = 2.0F * chi_1->real;
    chi_1_minus_alpha->imag = 2.0F * chi_1->imag;
    tau->real               = one_half.real;
    tau->imag               = one_half.imag;

    return FLA_SUCCESS;
  }

  //
  // Compute the absolute value (magnitude) of chi_1, which equals the 2-norm
  // of chi_1:
  //
  //   abs_chi_1 :=  | chi_1 |  =  || chi_1 ||_2
  //

  bl1_cnrm2( i_one,
             chi_1, i_one,
             &abs_chi_1 );

  //
  // Compute the 2-norm of x via the two norms previously computed above:
  //
  //   norm_x :=  || x ||_2  =  || / chi_1 \ ||   =  || / || chi_1 ||_2 \ ||
  //                            || \  x_2  / ||_2    || \  || x_2 ||_2  / ||_2
  //

  y[0].real = abs_chi_1;
  y[0].imag = 0.0;

  y[1].real = norm_x_2;
  y[1].imag = 0.0F;

  bl1_cnrm2( i_two,
             y, i_one,
             &norm_x );

  //
  // Compute alpha:
  //
  //   alpha := - || x ||_2 * chi_1 / | chi_1 |
  //

  if ( abs_chi_1 == 0.0F )
  {
    alpha->real = norm_x * ( -1.0F );
    alpha->imag = norm_x * ( -1.0F );
  }
  else
  {
    alpha->real = norm_x * ( -chi_1->real / abs_chi_1 );
    alpha->imag = norm_x * ( -chi_1->imag / abs_chi_1 );
  }

  chi_1_minus_alpha->real = chi_1->real - alpha->real;
  chi_1_minus_alpha->imag = chi_1->imag - alpha->imag;

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_2' * u_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_2' * x_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = 1/2 + ( || x ||_2 / | chi_1 - alpha | )^2
  //                                                                                                         
  bl1_csabsval2( chi_1_minus_alpha, &abs_chi_1_minus_alpha );

  norm_x_2_div_abs_chi_1_minus_alpha = norm_x_2 / abs_chi_1_minus_alpha;
  tau->real = one_half.real + one_half.real*( norm_x_2_div_abs_chi_1_minus_alpha *
                                              norm_x_2_div_abs_chi_1_minus_alpha );
  tau->imag = 0.0F;

  return FLA_SUCCESS;
}



FLA_Error FLA_Househ2s_UT_l_opz( int       m_x2,
                                 dcomplex* chi_1,
                                 dcomplex* x2, int inc_x2,
                                 dcomplex* alpha,
                                 dcomplex* chi_1_minus_alpha,
                                 dcomplex* tau )
{
  dcomplex one_half = *FLA_DOUBLE_COMPLEX_PTR( FLA_ONE_HALF );
  dcomplex y[2];
  double   abs_chi_1;
  double   norm_x_2;
  double   norm_x;
  double   abs_chi_1_minus_alpha;
  double   norm_x_2_div_abs_chi_1_minus_alpha;
  int      i_one = 1;
  int      i_two = 2;

  //
  // Compute the 2-norm of x_2:
  //
  //   norm_x_2 := || x_2 ||_2
  //

  bl1_znrm2( m_x2,
             x2, inc_x2,
             &norm_x_2 );

  //
  // If 2-norm of x_2 is zero, then return with trivial values.
  //

  if ( norm_x_2 == 0.0 )
  {
    alpha->real             = -(chi_1->real);
    alpha->imag             = -(chi_1->imag);
    chi_1_minus_alpha->real = 2.0 * chi_1->real;
    chi_1_minus_alpha->imag = 2.0 * chi_1->imag;
    tau->real               = one_half.real;
    tau->imag               = one_half.imag;

    return FLA_SUCCESS;
  }

  //
  // Compute the absolute value (magnitude) of chi_1, which equals the 2-norm
  // of chi_1:
  //
  //   abs_chi_1 :=  | chi_1 |  =  || chi_1 ||_2
  //

  bl1_znrm2( i_one,
             chi_1, i_one,
             &abs_chi_1 );

  //
  // Compute the 2-norm of x via the two norms previously computed above:
  //
  //   norm_x :=  || x ||_2  =  || / chi_1 \ ||   =  || / || chi_1 ||_2 \ ||
  //                            || \  x_2  / ||_2    || \  || x_2 ||_2  / ||_2
  //

  y[0].real = abs_chi_1;
  y[0].imag = 0.0;

  y[1].real = norm_x_2;
  y[1].imag = 0.0;

  bl1_znrm2( i_two,
             y, i_one,
             &norm_x );

  //
  // Compute alpha:
  //
  //   alpha := - || x ||_2 * chi_1 / | chi_1 |
  //

  if ( abs_chi_1 == 0.0 ) 
  {
    alpha->real = norm_x * ( -1.0 );
    alpha->imag = norm_x * ( -1.0 );
  }
  else
  {
    alpha->real = norm_x * ( -chi_1->real / abs_chi_1 );
    alpha->imag = norm_x * ( -chi_1->imag / abs_chi_1 );
  }

  chi_1_minus_alpha->real = chi_1->real - alpha->real;
  chi_1_minus_alpha->imag = chi_1->imag - alpha->imag;

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_2' * u_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_2' * x_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = 1/2 + ( || x ||_2 / | chi_1 - alpha | )^2
  //                                                                                                         
  bl1_zdabsval2( chi_1_minus_alpha, &abs_chi_1_minus_alpha );

  norm_x_2_div_abs_chi_1_minus_alpha = norm_x_2 / abs_chi_1_minus_alpha;
  tau->real = one_half.real + one_half.real*( norm_x_2_div_abs_chi_1_minus_alpha *
                                              norm_x_2_div_abs_chi_1_minus_alpha );
  tau->imag = 0.0;

  return FLA_SUCCESS;
}




FLA_Error FLA_Househ2s_UT_r_ops( int    m_x2,
                                 float* chi_1,
                                 float* x2, int inc_x2,
                                 float* alpha,
                                 float* chi_1_minus_alpha,
                                 float* tau )
{
  FLA_Househ2s_UT_l_ops( m_x2,
                         chi_1,
                         x2, inc_x2,
                         alpha,
                         chi_1_minus_alpha,
                         tau );

  return FLA_SUCCESS;
}

FLA_Error FLA_Househ2s_UT_r_opd( int     m_x2,
                                 double* chi_1,
                                 double* x2, int inc_x2,
                                 double* alpha,
                                 double* chi_1_minus_alpha,
                                 double* tau )
{
  FLA_Househ2s_UT_l_opd( m_x2,
                         chi_1,
                         x2, inc_x2,
                         alpha,
                         chi_1_minus_alpha,
                         tau );

  return FLA_SUCCESS;
}

FLA_Error FLA_Househ2s_UT_r_opc( int       m_x2,
                                 scomplex* chi_1,
                                 scomplex* x2, int inc_x2,
                                 scomplex* alpha,
                                 scomplex* chi_1_minus_alpha,
                                 scomplex* tau )
{
  FLA_Househ2s_UT_l_opc( m_x2,
                         chi_1,
                         x2, inc_x2,
                         alpha,
                         chi_1_minus_alpha,
                         tau );

  //chi_1_minus_alpha->real = chi_1->real -  alpha->real;
  //chi_1_minus_alpha->imag = chi_1->imag - -alpha->imag;

  return FLA_SUCCESS;
}


FLA_Error FLA_Househ2s_UT_r_opz( int       m_x2,
                                 dcomplex* chi_1,
                                 dcomplex* x2, int inc_x2,
                                 dcomplex* alpha,
                                 dcomplex* chi_1_minus_alpha,
                                 dcomplex* tau )
{
  FLA_Househ2s_UT_l_opz( m_x2,
                         chi_1,
                         x2, inc_x2,
                         alpha,
                         chi_1_minus_alpha,
                         tau );

  //chi_1_minus_alpha->real = chi_1->real -  alpha->real;
  //chi_1_minus_alpha->imag = chi_1->imag - -alpha->imag;

  return FLA_SUCCESS;
}

