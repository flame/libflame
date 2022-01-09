/*
    Copyright (c) 2021 Advanced Micro Devices, Inc.  All rights reserved.
    May 09, 2021
*/

#include "FLAME.h"

#define ssign( x ) ( (x) < 0.0F ? -1.0F : 1.0F )
#define dsign( x ) ( (x) < 0.0  ? -1.0  : 1.0  )

int slarf_(char *, integer *, integer *, float *, integer *, float *, float *, integer *, float *);
int dlarf_(char *, integer *, integer *, double *, integer *, double *, double *, integer *, double *);

FLA_Error FLA_EXT_Househ2_l_ops( integer  m_x2,
                                 float*   chi_1,
                                 float*   x2, integer inc_x2,
                                 float*   tau )
{
  float   one_half = 1.0F/2.0F;
  float   y[2];
  float   alpha;
  float   chi_1_minus_alpha, inv_chi_1_minus_alpha;
  float   norm_x_2;
  float   norm_x;
  float   safmin, rsafmn, lchi1;
  integer i_one = 1;
  integer i_two = 2;
  integer kn;

  //
  // Compute the 2-norm of x_2:
  //
  //   norm_x_2 := || x_2 ||_2
  //

  norm_x_2 = snrm2_( &m_x2,
                     x2, &inc_x2 );

  //
  // If 2-norm of x_2 is zero, then return with trivial values.
  //

  if ( norm_x_2 == 0.0 )
  {
    *chi_1 = -(*chi_1);
    *tau   = one_half;

    return FLA_SUCCESS;
  }

  //
  // Compute the 2-norm of x via the two norms previously computed above:
  //
  //   norm_x :=  || x ||_2  =  || / chi_1 \ ||   =  || / || chi_1 ||_2 \ ||
  //                            || \  x_2  / ||_2    || \  || x_2 ||_2  / ||_2
  //

  lchi1 = *chi_1;
  y[0] = lchi1;
  y[1] = norm_x_2;

  norm_x = snrm2_( &i_two,
                   y, &i_one );

  //
  // Compute alpha:
  //
  //   alpha := - || x ||_2 * chi_1 / | chi_1 |
  //          = -sign( chi_1 ) * || x ||_2
  //

  alpha = -ssign( lchi1 ) * norm_x;

  //
  // Overwrite x_2 with u_2:
  //
  //   x_2 := x_2 / ( chi_1 - alpha )
  //

  chi_1_minus_alpha = lchi1 - alpha;

  //
  // Scale X2, chi_1_minus_alpha, alpha
  // if norm factor is very less
  //
  safmin = fla_slamch("S", 1) / fla_slamch("E", 1);
  kn = 0;
  if( fabs( chi_1_minus_alpha ) < safmin )
  {
    rsafmn = 1. / safmin;

    for( kn = 1; kn < 20; kn++ )
    {
      sscal_( &m_x2,
              &rsafmn,
              x2, &inc_x2 );
      alpha = alpha * rsafmn;
      chi_1_minus_alpha = chi_1_minus_alpha * rsafmn;
      
      if( fabs( chi_1_minus_alpha ) > safmin )
        break;
    }
  }

  inv_chi_1_minus_alpha = 1.0F / chi_1_minus_alpha;

  sscal_( &m_x2,
          &inv_chi_1_minus_alpha,
          x2, &inc_x2 );

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_2' * u_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_2' * x_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = 1/2 + ( || x ||_2 / | chi_1 - alpha | )^2
  //        = alpha / ( alpha - chi_1 )
  //

  *tau = -1.0F * alpha * inv_chi_1_minus_alpha;

  //
  // Scale back alpha
  //
  for( ; kn > 0; kn-- )
  {
    alpha = alpha * safmin;
  }

  //
  // Overwrite chi_1 with alpha:
  //
  //   chi_1 := alpha
  //

  *chi_1 = alpha;

  return FLA_SUCCESS;
}

FLA_Error FLA_EXT_sgeqrf( integer  m_A, integer n_A,
                          float*   buff_A, integer cs_A,
                          float*   buff_t,
                          float*   buff_w,
                          integer* lwork,
                          integer* info )
{
  integer min_m_n = min( m_A, n_A );
  integer i, rs_A = 1;

  for ( i = 0; i < min_m_n; ++i )
  {
    float* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    float* tau1     = buff_t + i;
    float  alphat = *alpha11;

    integer m_curr   = m_A - i;

    integer m_ahead  = m_A - i - 1;
    integer n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // Compute Householder transformation for current column
    FLA_EXT_Househ2_l_ops( m_ahead,
                           &alphat,
                           a21, rs_A,
                           tau1 );

    *alpha11 = 1.0F;
    if( *tau1 != 0.0F )
        *tau1 = 1.0F / *tau1;

    // Apply the computed Householder transformation on the matrix
    slarf_( "Left",
            &m_curr, &n_ahead,
            alpha11, &rs_A,
            tau1,
            a12t, &cs_A,
            buff_w );

    *alpha11 = alphat;

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

FLA_Error FLA_EXT_Househ2_l_opd( integer   m_x2,
                                 double*   chi_1,
                                 double*   x2, integer inc_x2,
                                 double*   tau )
{
  double   one_half = 1.0/2.0;
  double   y[2];
  double   alpha;
  double   chi_1_minus_alpha, inv_chi_1_minus_alpha;
  double   norm_x_2;
  double   norm_x;
  double   safmin, rsafmn, lchi1;
  integer  i_one = 1;
  integer  i_two = 2;
  integer  kn;

  //
  // Compute the 2-norm of x_2:
  //
  //   norm_x_2 := || x_2 ||_2
  //

  norm_x_2 = dnrm2_( &m_x2,
                     x2, &inc_x2 );

  //
  // If 2-norm of x_2 is zero, then return with trivial values.
  //

  if ( norm_x_2 == 0.0 )
  {
    *chi_1 = -(*chi_1);
    *tau   = one_half;

    return FLA_SUCCESS;
  }

  //
  // Compute the 2-norm of x via the two norms previously computed above:
  //
  //   norm_x :=  || x ||_2  =  || / chi_1 \ ||   =  || / || chi_1 ||_2 \ ||
  //                            || \  x_2  / ||_2    || \  || x_2 ||_2  / ||_2
  //

  lchi1 = *chi_1;
  y[0] = lchi1;
  y[1] = norm_x_2;

  norm_x = dnrm2_( &i_two,
                   y, &i_one );

  //
  // Compute alpha:
  //
  //   alpha := - || x ||_2 * chi_1 / | chi_1 |
  //          = -sign( chi_1 ) * || x ||_2
  //

  alpha = -dsign( lchi1 ) * norm_x;

  //
  // Overwrite x_2 with u_2:
  //
  //   x_2 := x_2 / ( chi_1 - alpha )
  //

  chi_1_minus_alpha = lchi1 - alpha;

  //
  // Scale X2, chi_1_minus_alpha, alpha
  // if norm factor is very less
  //
  safmin = fla_dlamch("S", 1) / fla_dlamch("E", 1);
  kn = 0;
  if( fabs( chi_1_minus_alpha ) < safmin )
  {
    rsafmn = 1. / safmin;

    for( kn = 1; kn < 20; kn++ )
    {
      dscal_( &m_x2,
               &rsafmn,
               x2, &inc_x2 );
      alpha = alpha * rsafmn;
      chi_1_minus_alpha = chi_1_minus_alpha * rsafmn;
      
      if( fabs( chi_1_minus_alpha ) > safmin )
        break;
    }
  }

  inv_chi_1_minus_alpha = 1.0 / chi_1_minus_alpha;

  dscal_( &m_x2,
          &inv_chi_1_minus_alpha,
          x2, &inc_x2 );

  //
  // Compute tau:
  //
  //   tau := ( 1 + u_2' * u_2 ) / 2
  //        = ( ( chi_1 - alpha ) * conj( chi_1 - alpha ) + x_2' * x_2 ) /
  //          ( 2 * ( chi_1 - alpha ) * conj( chi_1 - alpha ) )
  //        = 1/2 + ( || x ||_2 / | chi_1 - alpha | )^2
  //        = alpha / ( alpha - chi_1 )
  //

  *tau = -1.0 * alpha * inv_chi_1_minus_alpha;

  //
  // Scale back alpha
  //
  for( ; kn > 0; kn-- )
  {
    alpha = alpha * safmin;
  }

  //
  // Overwrite chi_1 with alpha:
  //
  //   chi_1 := alpha
  //

  *chi_1 = alpha;

  return FLA_SUCCESS;
}

FLA_Error FLA_EXT_dgeqrf( integer  m_A, integer n_A,
                          double*  buff_A, integer cs_A,
                          double*  buff_t,
                          double*  buff_w,
                          integer* lwork,
                          integer* info )
{
  integer min_m_n = min( m_A, n_A );
  integer i, rs_A = 1;

  for ( i = 0; i < min_m_n; ++i )
  {
    double* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double* a12t     = buff_A + (i+1)*cs_A + (i  )*rs_A;

    double* tau1     = buff_t + i;
    double  alphat   = *alpha11;

    integer m_curr   = m_A - i;

    integer m_ahead  = m_A - i - 1;
    integer n_ahead  = n_A - i - 1;

    /*------------------------------------------------------------*/

    // Compute Householder transformation for current column
    FLA_EXT_Househ2_l_opd( m_ahead,
                           &alphat,
                           a21, rs_A,
                           tau1 );

    *alpha11 = 1.0;
    if( *tau1 != 0 )
        *tau1 = 1.0 / *tau1;

    // Apply the computed Householder transformation on the matrix
    dlarf_( "L",
            &m_curr, &n_ahead,
            alpha11, &rs_A,
            tau1,
            a12t, &cs_A,
            buff_w );

    *alpha11 = alphat;

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}
