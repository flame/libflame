
#include "FLAME.h"

FLA_Error FLA_Fused_Her2_Ax_l_opt_var1( FLA_Obj beta, FLA_Obj u, FLA_Obj z, FLA_Obj A, FLA_Obj x, FLA_Obj w )
{
/*
   Effective computation:
   A = A + beta * ( u * z' + z * u' );
   w = A * x;
*/
  FLA_Datatype datatype;
  int          m_A;
  int          rs_A, cs_A;
  int          inc_u, inc_z, inc_x, inc_w;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_u    = FLA_Obj_vector_inc( u );
  inc_z    = FLA_Obj_vector_inc( z );
  inc_x    = FLA_Obj_vector_inc( x );
  inc_w    = FLA_Obj_vector_inc( w );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_u = FLA_FLOAT_PTR( u );
      float* buff_z = FLA_FLOAT_PTR( z );
      float* buff_x = FLA_FLOAT_PTR( x );
      float* buff_w = FLA_FLOAT_PTR( w );
      float* buff_beta = FLA_FLOAT_PTR( beta );

      FLA_Fused_Her2_Ax_l_ops_var1( m_A,
                                    buff_beta,
                                    buff_u, inc_u,
                                    buff_z, inc_z,
                                    buff_A, rs_A, cs_A,
                                    buff_x, inc_x,
                                    buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_u = FLA_DOUBLE_PTR( u );
      double* buff_z = FLA_DOUBLE_PTR( z );
      double* buff_x = FLA_DOUBLE_PTR( x );
      double* buff_w = FLA_DOUBLE_PTR( w );
      double* buff_beta = FLA_DOUBLE_PTR( beta );

      FLA_Fused_Her2_Ax_l_opd_var1( m_A,
                                    buff_beta,
                                    buff_u, inc_u,
                                    buff_z, inc_z,
                                    buff_A, rs_A, cs_A,
                                    buff_x, inc_x,
                                    buff_w, inc_w );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_u = FLA_COMPLEX_PTR( u );
      scomplex* buff_z = FLA_COMPLEX_PTR( z );
      scomplex* buff_x = FLA_COMPLEX_PTR( x );
      scomplex* buff_w = FLA_COMPLEX_PTR( w );
      scomplex* buff_beta = FLA_COMPLEX_PTR( beta );

      FLA_Fused_Her2_Ax_l_opc_var1( m_A,
                                    buff_beta,
                                    buff_u, inc_u,
                                    buff_z, inc_z,
                                    buff_A, rs_A, cs_A,
                                    buff_x, inc_x,
                                    buff_w, inc_w );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_u = FLA_DOUBLE_COMPLEX_PTR( u );
      dcomplex* buff_z = FLA_DOUBLE_COMPLEX_PTR( z );
      dcomplex* buff_x = FLA_DOUBLE_COMPLEX_PTR( x );
      dcomplex* buff_w = FLA_DOUBLE_COMPLEX_PTR( w );
      dcomplex* buff_beta = FLA_DOUBLE_COMPLEX_PTR( beta );

      FLA_Fused_Her2_Ax_l_opz_var1( m_A,
                                    buff_beta,
                                    buff_u, inc_u,
                                    buff_z, inc_z,
                                    buff_A, rs_A, cs_A,
                                    buff_x, inc_x,
                                    buff_w, inc_w );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Her2_Ax_l_ops_var1( int m_A,
                                        float* buff_beta, 
                                        float* buff_u, int inc_u, 
                                        float* buff_z, int inc_z, 
                                        float* buff_A, int rs_A, int cs_A, 
                                        float* buff_x, int inc_x, 
                                        float* buff_w, int inc_w )
{
  float*    buff_0  = FLA_FLOAT_PTR( FLA_ZERO );
  int       i;

  bl1_ssetv( m_A,
             buff_0,
             buff_w, inc_w );

  for ( i = 0; i < m_A; ++i )
  {
    float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    float*    upsilon1 = buff_u + (i  )*inc_u;
    float*    u2       = buff_u + (i+1)*inc_u;

    float*    zeta1    = buff_z + (i  )*inc_z;
    float*    z2       = buff_z + (i+1)*inc_z;

    float*    chi1     = buff_x + (i  )*inc_x;
    float*    x2       = buff_x + (i+1)*inc_x;

    float*    omega1   = buff_w + (i  )*inc_w;
    float*    w2       = buff_w + (i+1)*inc_w;

    // float*    beta     = buff_beta;

    float    minus_conj_upsilon1;
    float    minus_conj_zeta1;
    float    temp;

    int      m_ahead   = m_A - i - 1;

    /*------------------------------------------------------------*/

    // bl1_scopyconj( zeta1, &conj_zeta1 );
    // bl1_smult3( beta, &conj_zeta1, &minus_conj_zeta1 );
    // bl1_smult3( &minus_conj_zeta1, upsilon1, &temp );
    // bl1_sadd3( &temp, alpha11, alpha11 );

    //bl1_scopyconj( upsilon1, &conj_upsilon1 );
    //bl1_smult3( beta, &conj_upsilon1, &minus_conj_upsilon1 );
    //bl1_smult3( &minus_conj_upsilon1, zeta1, &temp );
    //bl1_sadd3( &temp, alpha11, alpha11 );
    minus_conj_zeta1    = - *zeta1;
    minus_conj_upsilon1 = - *upsilon1;

    *alpha11 -= 2.0F * *zeta1 * *upsilon1;

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &minus_conj_zeta1,
                u2,  inc_u,
                a21, rs_A );
/*
    F77_saxpy( &m_ahead,
               &minus_conj_zeta1,
               u2,  &inc_u,
               a21, &rs_A );
*/


    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &minus_conj_upsilon1,
                z2,  inc_z,
                a21, rs_A );
/*
    F77_saxpy( &m_ahead,
               &minus_conj_upsilon1,
               z2,  &inc_z,
               a21, &rs_A );
*/

    // bl1_smult3( alpha11, chi1, &temp );
    // bl1_sadd3( &temp, omega1, omega1 );
    *omega1 += *alpha11 * *chi1;

    bl1_sdot( BLIS1_CONJUGATE,
              m_ahead,
              a21, rs_A,
              x2,  inc_x,
              &temp );
/*
    temp = F77_sdot( &m_ahead,
                     a21, &rs_A,
                     x2,  &inc_x );
*/

    // bl1_sadd3( &temp, omega1, omega1 );
    *omega1 += temp;

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                chi1,
                a21, rs_A,
                w2,  inc_w );
/*
    F77_saxpy( &m_ahead,
               chi1,
               a21, &rs_A,
               w2,  &inc_w );
*/

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Her2_Ax_l_opd_var1( int m_A,
                                        double* buff_beta, 
                                        double* buff_u, int inc_u, 
                                        double* buff_z, int inc_z, 
                                        double* buff_A, int rs_A, int cs_A, 
                                        double* buff_x, int inc_x, 
                                        double* buff_w, int inc_w )
{
  double*   buff_0  = FLA_DOUBLE_PTR( FLA_ZERO );
  int       i;

  bl1_dsetv( m_A,
             buff_0,
             buff_w, inc_w );

  for ( i = 0; i < m_A; ++i )
  {
    double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    double*   upsilon1 = buff_u + (i  )*inc_u;
    double*   u2       = buff_u + (i+1)*inc_u;

    double*   zeta1    = buff_z + (i  )*inc_z;
    double*   z2       = buff_z + (i+1)*inc_z;

    double*   chi1     = buff_x + (i  )*inc_x;
    double*   x2       = buff_x + (i+1)*inc_x;

    double*   omega1   = buff_w + (i  )*inc_w;
    double*   w2       = buff_w + (i+1)*inc_w;

    // double*   beta     = buff_beta;

    double   minus_conj_upsilon1;
    double   minus_conj_zeta1;
    double   temp;

    int      m_ahead   = m_A - i - 1;

    /*------------------------------------------------------------*/

    // bl1_dcopyconj( zeta1, &conj_zeta1 );
    // bl1_dmult3( beta, &conj_zeta1, &minus_conj_zeta1 );
    // bl1_dmult3( &minus_conj_zeta1, upsilon1, &temp );
    // bl1_dadd3( &temp, alpha11, alpha11 );

    //bl1_dcopyconj( upsilon1, &conj_upsilon1 );
    //bl1_dmult3( beta, &conj_upsilon1, &minus_conj_upsilon1 );
    //bl1_dmult3( &minus_conj_upsilon1, zeta1, &temp );
    //bl1_dadd3( &temp, alpha11, alpha11 );
    minus_conj_zeta1    = - *zeta1;
    minus_conj_upsilon1 = - *upsilon1;

    *alpha11 -= 2.0 * *zeta1 * *upsilon1;

    // bl1_dmult3( alpha11, chi1, &temp );
    // bl1_dadd3( &temp, omega1, omega1 );
    *omega1 += *alpha11 * *chi1;

    bl1_daxpyv2bdotaxpy( m_ahead,
                         &minus_conj_zeta1,
                         u2,  inc_u,
                         &minus_conj_upsilon1,
                         z2,  inc_z,
                         a21, rs_A,
                         x2,  inc_x,
                         chi1,
                         &temp,
                         w2,  inc_w );

    // bl1_dadd3( &temp, omega1, omega1 );
    *omega1 += temp;

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Her2_Ax_l_opc_var1( int m_A,
                                        scomplex* buff_beta, 
                                        scomplex* buff_u, int inc_u, 
                                        scomplex* buff_z, int inc_z, 
                                        scomplex* buff_A, int rs_A, int cs_A, 
                                        scomplex* buff_x, int inc_x, 
                                        scomplex* buff_w, int inc_w )
{
  scomplex* buff_0  = FLA_COMPLEX_PTR( FLA_ZERO );
  int       i;

  bl1_csetv( m_A,
             buff_0,
             buff_w, inc_w );

  for ( i = 0; i < m_A; ++i )
  {
    scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    scomplex* upsilon1 = buff_u + (i  )*inc_u;
    scomplex* u2       = buff_u + (i+1)*inc_u;

    scomplex* zeta1    = buff_z + (i  )*inc_z;
    scomplex* z2       = buff_z + (i+1)*inc_z;

    scomplex* chi1     = buff_x + (i  )*inc_x;
    scomplex* x2       = buff_x + (i+1)*inc_x;

    scomplex* omega1   = buff_w + (i  )*inc_w;
    scomplex* w2       = buff_w + (i+1)*inc_w;

    // scomplex* beta     = buff_beta;

    scomplex minus_conj_upsilon1;
    scomplex minus_conj_zeta1;
    scomplex temp;

    int      m_ahead   = m_A - i - 1;

    /*------------------------------------------------------------*/

    // bl1_ccopyconj( zeta1, &conj_zeta1 );
    // bl1_cmult3( beta, &conj_zeta1, &minus_conj_zeta1 );
    // bl1_cmult3( &minus_conj_zeta1, upsilon1, &temp );
    // bl1_cadd3( &temp, alpha11, alpha11 );

    //bl1_ccopyconj( upsilon1, &conj_upsilon1 );
    //bl1_cmult3( beta, &conj_upsilon1, &minus_conj_upsilon1 );
    //bl1_cmult3( &minus_conj_upsilon1, zeta1, &temp );
    //bl1_cadd3( &temp, alpha11, alpha11 );
    minus_conj_zeta1.real    = -  zeta1->real;
    minus_conj_zeta1.imag    = - -zeta1->imag;
    minus_conj_upsilon1.real = -  upsilon1->real;
    minus_conj_upsilon1.imag = - -upsilon1->imag;

    alpha11->real -=  zeta1->real * upsilon1->real - -zeta1->imag * upsilon1->imag +
                      zeta1->real * upsilon1->real - zeta1->imag * -upsilon1->imag;
    alpha11->imag -= -zeta1->imag * upsilon1->real +  zeta1->real * upsilon1->imag +
                      zeta1->imag * upsilon1->real + zeta1->real * -upsilon1->imag;

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &minus_conj_zeta1,
                u2,  inc_u,
                a21, rs_A );
/*
    F77_caxpy( &m_ahead,
               &minus_conj_zeta1,
               u2,  &inc_u,
               a21, &rs_A );
*/


    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                &minus_conj_upsilon1,
                z2,  inc_z,
                a21, rs_A );
/*
    F77_caxpy( &m_ahead,
               &minus_conj_upsilon1,
               z2,  &inc_z,
               a21, &rs_A );
*/

    // bl1_cmult3( alpha11, chi1, &temp );
    // bl1_cadd3( &temp, omega1, omega1 );
    omega1->real += alpha11->real * chi1->real - alpha11->imag * chi1->imag;
    omega1->imag += alpha11->imag * chi1->real + alpha11->real * chi1->imag;

    bl1_cdot( BLIS1_CONJUGATE,
              m_ahead,
              a21, rs_A,
              x2,  inc_x,
              &temp );
    // bl1_cadd3( &temp, omega1, omega1 );
    omega1->real += temp.real;
    omega1->imag += temp.imag;

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_ahead,
                chi1,
                a21, rs_A,
                w2,  inc_w );
/*
    F77_caxpy( &m_ahead,
               chi1,
               a21, &rs_A,
               w2,  &inc_w );
*/

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Her2_Ax_l_opz_var1( int m_A,
                                        dcomplex* buff_beta, 
                                        dcomplex* buff_u, int inc_u, 
                                        dcomplex* buff_z, int inc_z, 
                                        dcomplex* buff_A, int rs_A, int cs_A, 
                                        dcomplex* buff_x, int inc_x, 
                                        dcomplex* buff_w, int inc_w )
{
  dcomplex  zero = bl1_z0();
  int       i;

  bl1_zsetv( m_A,
             &zero,
             buff_w, inc_w );

  for ( i = 0; i < m_A; ++i )
  {
    dcomplex* restrict alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* restrict a21      = buff_A + (i  )*cs_A + (i+1)*rs_A;

    dcomplex* restrict upsilon1 = buff_u + (i  )*inc_u;
    dcomplex* restrict u2       = buff_u + (i+1)*inc_u;

    dcomplex* restrict zeta1    = buff_z + (i  )*inc_z;
    dcomplex* restrict z2       = buff_z + (i+1)*inc_z;

    dcomplex* restrict chi1     = buff_x + (i  )*inc_x;
    dcomplex* restrict x2       = buff_x + (i+1)*inc_x;

    dcomplex* restrict omega1   = buff_w + (i  )*inc_w;
    dcomplex* restrict w2       = buff_w + (i+1)*inc_w;

    //dcomplex* restrict beta     = buff_beta;

    dcomplex minus_conj_upsilon1;
    dcomplex minus_conj_zeta1;
    dcomplex temp;

    dcomplex ze1;
    dcomplex up1;
    dcomplex a11;
    dcomplex om1;
    dcomplex ch1;
    
    int      m_ahead   = m_A - i - 1;

    /*------------------------------------------------------------*/

    // bl1_zcopyconj( zeta1, &conj_zeta1 );
    // bl1_zmult3( beta, &conj_zeta1, &minus_conj_zeta1 );
    // bl1_zmult3( &minus_conj_zeta1, upsilon1, &temp );
    // bl1_zadd3( &temp, alpha11, alpha11 );

    //bl1_zcopyconj( upsilon1, &conj_upsilon1 );
    //bl1_zmult3( beta, &conj_upsilon1, &minus_conj_upsilon1 );
    //bl1_zmult3( &minus_conj_upsilon1, zeta1, &temp );
    //bl1_zadd3( &temp, alpha11, alpha11 );
    minus_conj_zeta1.real    = -  zeta1->real;
    minus_conj_zeta1.imag    = - -zeta1->imag;
    minus_conj_upsilon1.real = -  upsilon1->real;
    minus_conj_upsilon1.imag = - -upsilon1->imag;

    ze1 = *zeta1;
    up1 = *upsilon1;
    a11 = *alpha11;
    om1 = *omega1;
    ch1 = *chi1;
    
    //alpha11->real -=  zeta1->real * upsilon1->real - -zeta1->imag * upsilon1->imag +
    //                  zeta1->real * upsilon1->real - zeta1->imag * -upsilon1->imag;
    //alpha11->imag -= -zeta1->imag * upsilon1->real +  zeta1->real * upsilon1->imag +
    //                  zeta1->imag * upsilon1->real + zeta1->real * -upsilon1->imag;
    a11.real -= ze1.real * up1.real - -ze1.imag * up1.imag +
                up1.real * ze1.real - -up1.imag * ze1.imag;
    a11.imag -= ze1.real * up1.imag + -ze1.imag * up1.real +
                up1.real * ze1.imag + -up1.imag * ze1.real;

    // bl1_zmult3( alpha11, chi1, &temp );
    // bl1_zadd3( &temp, omega1, omega1 );
    //omega1->real += alpha11->real * chi1->real - alpha11->imag * chi1->imag;
    //omega1->imag += alpha11->imag * chi1->real + alpha11->real * chi1->imag;
    om1.real += a11.real * ch1.real - a11.imag * ch1.imag;
    om1.imag += a11.imag * ch1.real + a11.real * ch1.imag;

    *alpha11 = a11;
    *omega1  = om1;

/*
    bl1_zaxpyv2bdotaxpy( m_ahead,
                         &minus_conj_zeta1,
                         u2,  inc_u,
                         &minus_conj_upsilon1,
                         z2,  inc_z,
                         a21, rs_A,
                         x2,  inc_x,
                         chi1,
                         &temp,
                         w2,  inc_w );
*/

    bl1_zaxpyv2b( m_ahead,
                  &minus_conj_zeta1,
                  &minus_conj_upsilon1,
                  u2,  inc_u,
                  z2,  inc_z,
                  a21, rs_A );

    bl1_zdotaxpy( m_ahead,
                  a21, rs_A,
                  x2,  inc_x,
                  chi1,
                  &temp,
                  w2,  inc_w );


    // bl1_zadd3( &temp, omega1, omega1 );
    omega1->real += temp.real;
    omega1->imag += temp.imag;

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

