
#include "FLAME.h"

FLA_Error FLA_Bidiag_UT_realify( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Error r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Bidiag_UT_realify_check( A, d, e );

  if ( FLA_Obj_is_real( A ) )
  {
    FLA_Set( FLA_ONE, d );
    FLA_Set( FLA_ONE, e );
    return FLA_SUCCESS;
  }

  if ( FLA_Obj_length( A ) < FLA_Obj_width( A ) )
    //r_val = FLA_Bidiag_UT_l_realify_unb( A, d, e );
    r_val = FLA_Bidiag_UT_l_realify_opt( A, d, e );
  else
    //r_val = FLA_Bidiag_UT_u_realify_unb( A, d, e );
    r_val = FLA_Bidiag_UT_u_realify_opt( A, d, e );

  return r_val;
}


FLA_Error FLA_Bidiag_UT_l_realify_unb( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj dT,              d0,
          dB,              delta1,
                           d2;

  FLA_Obj eT,              e0,
          eB,              epsilon1,
                           e2;

  FLA_Obj a10t_l, a10t_r;

  FLA_Obj a21_t,
          a21_b;

  FLA_Obj absv;


  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, &absv );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( d,    &dT,
                      &dB,            0, FLA_TOP );

  FLA_Part_2x1( e,    &eT,
                      &eB,            0, FLA_TOP );

  while ( FLA_Obj_min_dim( ABR ) > 0 )
  {
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( dT,                &d0, 
                        /* ** */            /* ****** */
                                              &delta1, 
                           dB,                &d2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( eT,                &e0, 
                        /* ** */            /* ******** */
                                              &epsilon1, 
                           eB,                &e2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    if ( FLA_Obj_width( a10t ) == 0 )
    {
      // delta1 = 1;
      FLA_Set( FLA_ONE, delta1 );
    }
    else
    {
      FLA_Part_1x2( a10t,   &a10t_l, &a10t_r,    1, FLA_RIGHT );

      // delta1 = conj(a10t_r) / abs(a10t_r); 
      FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a10t_r, delta1 );
      FLA_Copyt( FLA_NO_TRANSPOSE, a10t_r, absv );
      FLA_Absolute_value( absv );
      FLA_Inv_scal( absv, delta1 );

      // a10t_r = delta1 * a10t_r;
      // a10t_r.imag = 0;
      FLA_Scalc( FLA_NO_CONJUGATE, delta1, a10t_r );
      FLA_Obj_set_imag_part( FLA_ZERO, a10t_r );

      // alpha11 = delta1 * alpha11;
      FLA_Scalc( FLA_NO_CONJUGATE, delta1, alpha11 );
    }

    // epsilon1 = conj(alpha11) / abs(alpha11);
    FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha11, epsilon1 );
    FLA_Copyt( FLA_NO_TRANSPOSE, alpha11, absv );
    FLA_Absolute_value( absv );
    FLA_Inv_scal( absv, epsilon1 );

    // alpha11 = epsilon1 * alpha11;
    // alpha11.imag = 0;
    FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, alpha11 );
    FLA_Obj_set_imag_part( FLA_ZERO, alpha11 );

    if ( FLA_Obj_length( a21 ) > 0 )
    {
      FLA_Part_2x1( a21,   &a21_t,
                           &a21_b,    1, FLA_TOP );

      // a21_t = epsilon1 * a21_t;
      FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, a21_t );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &dT,                d0, 
                                                  delta1, 
                            /* ** */           /* ****** */
                              &dB,                d2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &eT,                e0, 
                                                  epsilon1, 
                            /* ** */           /* ******** */
                              &eB,                e2,     FLA_TOP );
  }

  FLA_Obj_free( &absv );

  return FLA_SUCCESS;
}

FLA_Error FLA_Bidiag_UT_l_realify_opt( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          min_m_n;
  int          rs_A, cs_A;
  int          inc_d;
  int          inc_e;
  int          i;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  min_m_n  = FLA_Obj_min_dim( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_d    = FLA_Obj_vector_inc( d );

  inc_e    = FLA_Obj_vector_inc( e );


  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_d = FLA_FLOAT_PTR( d );
      float* buff_e = FLA_FLOAT_PTR( e );
      float* buff_1 = FLA_FLOAT_PTR( FLA_ONE );

      bl1_ssetv( min_m_n,
                 buff_1,
                 buff_d, inc_d );

      bl1_ssetv( min_m_n,
                 buff_1,
                 buff_e, inc_e );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_d = FLA_DOUBLE_PTR( d );
      double* buff_e = FLA_DOUBLE_PTR( e );
      double* buff_1 = FLA_DOUBLE_PTR( FLA_ONE );

      bl1_dsetv( min_m_n,
                 buff_1,
                 buff_d, inc_d );

      bl1_dsetv( min_m_n,
                 buff_1,
                 buff_e, inc_e );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_d = FLA_COMPLEX_PTR( d );
      scomplex* buff_e = FLA_COMPLEX_PTR( e );
      scomplex* buff_1 = FLA_COMPLEX_PTR( FLA_ONE );
      float*    buff_0 = FLA_FLOAT_PTR( FLA_ZERO );

      for ( i = 0; i < min_m_n; ++i )
      {

        scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        scomplex* delta1   = buff_d + (i  )*inc_d;
        scomplex* epsilon1 = buff_e + (i  )*inc_e;
        scomplex  absv;

        int       m_ahead  = m_A - i - 1;
        int       m_behind = i;

        if ( m_behind == 0 )
        {
          // FLA_Set( FLA_ONE, delta1 );
          *delta1 = *buff_1;
        }
        else
        {
          scomplex* a10t_r   = buff_A + (i-1)*cs_A + (i  )*rs_A;
          // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a10t_r, delta1 );
          // FLA_Copyt( FLA_NO_TRANSPOSE, a10t_r, absv );
          // FLA_Absolute_value( absv );
          // FLA_Inv_scal( absv, delta1 );
          bl1_ccopys( BLIS1_CONJUGATE, a10t_r, delta1 );
          bl1_cabsval2( a10t_r, &absv );
          bl1_cinvscals( &absv, delta1 );

          // FLA_Scalc( FLA_NO_CONJUGATE, delta1, a10t_r );
          // FLA_Obj_set_imag_part( FLA_ZERO, a10t_r );
          bl1_cscals( delta1, a10t_r );
          a10t_r->imag = *buff_0;

          // FLA_Scalc( FLA_NO_CONJUGATE, delta1, alpha11 );
          bl1_cscals( delta1, alpha11 );
        }

        // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha11, epsilon1 );
        // FLA_Copyt( FLA_NO_TRANSPOSE, alpha11, absv );
        // FLA_Absolute_value( absv );
        // FLA_Inv_scal( absv, epsilon1 );
        bl1_ccopys( BLIS1_CONJUGATE, alpha11, epsilon1 );
        bl1_cabsval2( alpha11, &absv );
        bl1_cinvscals( &absv, epsilon1 );

        // FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, alpha11 );
        // FLA_Obj_set_imag_part( FLA_ZERO, alpha11 );
        bl1_cscals( epsilon1, alpha11 );
        alpha11->imag = *buff_0;

        if ( m_ahead > 0 )
        {
          scomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
          // FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, a21_t );
          bl1_cscals( epsilon1, a21_t );
        }
      }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_d = FLA_DOUBLE_COMPLEX_PTR( d );
      dcomplex* buff_e = FLA_DOUBLE_COMPLEX_PTR( e );
      dcomplex* buff_1 = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
      double*   buff_0 = FLA_DOUBLE_PTR( FLA_ZERO );

      for ( i = 0; i < min_m_n; ++i )
      {

        dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        dcomplex* delta1   = buff_d + (i  )*inc_d;
        dcomplex* epsilon1 = buff_e + (i  )*inc_e;
        dcomplex  absv;

        int       m_ahead  = m_A - i - 1;
        int       m_behind = i;

        if ( m_behind == 0 )
        {
          // FLA_Set( FLA_ONE, delta1 );
          *delta1 = *buff_1;
        }
        else
        {
          dcomplex* a10t_r   = buff_A + (i-1)*cs_A + (i  )*rs_A;
          // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a10t_r, delta1 );
          // FLA_Copyt( FLA_NO_TRANSPOSE, a10t_r, absv );
          // FLA_Absolute_value( absv );
          // FLA_Inv_scal( absv, delta1 );
          bl1_zcopys( BLIS1_CONJUGATE, a10t_r, delta1 );
          bl1_zabsval2( a10t_r, &absv );
          bl1_zinvscals( &absv, delta1 );

          // FLA_Scalc( FLA_NO_CONJUGATE, delta1, a10t_r );
          // FLA_Obj_set_imag_part( FLA_ZERO, a10t_r );
          bl1_zscals( delta1, a10t_r );
          a10t_r->imag = *buff_0;

          // FLA_Scalc( FLA_NO_CONJUGATE, delta1, alpha11 );
          bl1_zscals( delta1, alpha11 );
        }

        // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha11, epsilon1 );
        // FLA_Copyt( FLA_NO_TRANSPOSE, alpha11, absv );
        // FLA_Absolute_value( absv );
        // FLA_Inv_scal( absv, epsilon1 );
        bl1_zcopys( BLIS1_CONJUGATE, alpha11, epsilon1 );
        bl1_zabsval2( alpha11, &absv );
        bl1_zinvscals( &absv, epsilon1 );

        // FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, alpha11 );
        // FLA_Obj_set_imag_part( FLA_ZERO, alpha11 );
        bl1_zscals( epsilon1, alpha11 );
        alpha11->imag = *buff_0;

        if ( m_ahead > 0 )
        {
          dcomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
          // FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, a21_t );
          bl1_zscals( epsilon1, a21_t );
        }
      }

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Bidiag_UT_u_realify_unb( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj dT,              d0,
          dB,              delta1,
                           d2;

  FLA_Obj eT,              e0,
          eB,              epsilon1,
                           e2;

  FLA_Obj a01_t,
          a01_b;

  FLA_Obj a12t_l, a12t_r;

  FLA_Obj absv;


  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, &absv );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     0, 0, FLA_TL );

  FLA_Part_2x1( d,    &dT,
                      &dB,            0, FLA_TOP );

  FLA_Part_2x1( e,    &eT,
                      &eB,            0, FLA_TOP );

  while ( FLA_Obj_min_dim( ABR ) > 0 )
  {
    FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00,  /**/ &a01,     &A02,
                        /* ************* */   /* ************************** */
                                                &a10t, /**/ &alpha11, &a12t,
                           ABL, /**/ ABR,       &A20,  /**/ &a21,     &A22,
                           1, 1, FLA_BR );

    FLA_Repart_2x1_to_3x1( dT,                &d0, 
                        /* ** */            /* ****** */
                                              &delta1, 
                           dB,                &d2,        1, FLA_BOTTOM );

    FLA_Repart_2x1_to_3x1( eT,                &e0, 
                        /* ** */            /* ******** */
                                              &epsilon1, 
                           eB,                &e2,        1, FLA_BOTTOM );

    /*------------------------------------------------------------*/

    if ( FLA_Obj_length( a01 ) == 0 )
    {
      // epsilon1 = 1;
      FLA_Set( FLA_ONE, epsilon1 );
    }
    else
    {
      FLA_Part_2x1( a01,   &a01_t,
                           &a01_b,    1, FLA_BOTTOM );

      // epsilon1 = conj(a01_b) / abs(a01_b); 
      FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a01_b, epsilon1 );
      FLA_Copyt( FLA_NO_TRANSPOSE, a01_b, absv );
      FLA_Absolute_value( absv );
      FLA_Inv_scal( absv, epsilon1 );

      // a01_b = epsilon1 * a01_b;
      // a01_b.imag = 0;
      FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, a01_b );
      FLA_Obj_set_imag_part( FLA_ZERO, a01_b );

      // alpha11 = epsilon1 * alpha11;
      FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, alpha11 );
    }

    // delta1 = conj(alpha11) / abs(alpha11);
    FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha11, delta1 );
    FLA_Copyt( FLA_NO_TRANSPOSE, alpha11, absv );
    FLA_Absolute_value( absv );
    FLA_Inv_scal( absv, delta1 );

    // alpha11 = delta1 * alpha11;
    // alpha11.imag = 0;
    FLA_Scalc( FLA_NO_CONJUGATE, delta1, alpha11 );
    FLA_Obj_set_imag_part( FLA_ZERO, alpha11 );

    if ( FLA_Obj_width( a12t ) > 0 )
    {
      FLA_Part_1x2( a12t,   &a12t_l, &a12t_r,    1, FLA_LEFT );

      // a12t_l = delta1 * a12t_l;
      FLA_Scalc( FLA_NO_CONJUGATE, delta1, a12t_l );
    }

    /*------------------------------------------------------------*/

    FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00,  a01,     /**/ A02,
                                                     a10t, alpha11, /**/ a12t,
                            /* ************** */  /* ************************ */
                              &ABL, /**/ &ABR,       A20,  a21,     /**/ A22,
                              FLA_TL );

    FLA_Cont_with_3x1_to_2x1( &dT,                d0, 
                                                  delta1, 
                            /* ** */           /* ****** */
                              &dB,                d2,     FLA_TOP );

    FLA_Cont_with_3x1_to_2x1( &eT,                e0, 
                                                  epsilon1, 
                            /* ** */           /* ******** */
                              &eB,                e2,     FLA_TOP );
  }

  FLA_Obj_free( &absv );

  return FLA_SUCCESS;
}

FLA_Error FLA_Bidiag_UT_u_realify_opt( FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          min_m_n;
  int          rs_A, cs_A;
  int          inc_d;
  int          inc_e;
  int          i;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  min_m_n  = FLA_Obj_min_dim( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_d    = FLA_Obj_vector_inc( d );

  inc_e    = FLA_Obj_vector_inc( e );


  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_d = FLA_FLOAT_PTR( d );
      float* buff_e = FLA_FLOAT_PTR( e );
      float* buff_1 = FLA_FLOAT_PTR( FLA_ONE );

      bl1_ssetv( min_m_n,
                 buff_1,
                 buff_d, inc_d );

      bl1_ssetv( min_m_n,
                 buff_1,
                 buff_e, inc_e );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_d = FLA_DOUBLE_PTR( d );
      double* buff_e = FLA_DOUBLE_PTR( e );
      double* buff_1 = FLA_DOUBLE_PTR( FLA_ONE );

      bl1_dsetv( min_m_n,
                 buff_1,
                 buff_d, inc_d );

      bl1_dsetv( min_m_n,
                 buff_1,
                 buff_e, inc_e );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_d = FLA_COMPLEX_PTR( d );
      scomplex* buff_e = FLA_COMPLEX_PTR( e );
      scomplex* buff_1 = FLA_COMPLEX_PTR( FLA_ONE );
      float*    buff_0 = FLA_FLOAT_PTR( FLA_ZERO );

      for ( i = 0; i < min_m_n; ++i )
      {
        scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        scomplex* delta1   = buff_d + (i  )*inc_d;
        scomplex* epsilon1 = buff_e + (i  )*inc_e;
        scomplex  absv;

        int       n_ahead  = n_A - i - 1;
        int       n_behind = i;

        if ( n_behind == 0 )
        {
          // FLA_Set( FLA_ONE, epsilon1 );
          *epsilon1 = *buff_1;
        }
        else
        {
          scomplex* a01_b    = buff_A + (i  )*cs_A + (i-1)*rs_A;
          // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a01_b, epsilon1 );
          // FLA_Copyt( FLA_NO_TRANSPOSE, a01_b, absv );
          // FLA_Absolute_value( absv );
          // FLA_Inv_scal( absv, epsilon1 );
          bl1_ccopys( BLIS1_CONJUGATE, a01_b, epsilon1 );
          bl1_cabsval2( a01_b, &absv );
          bl1_cinvscals( &absv, epsilon1 );

          // FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, a01_b );
          // FLA_Obj_set_imag_part( FLA_ZERO, a01_b );
          bl1_cscals( epsilon1, a01_b );
          a01_b->imag = *buff_0;

          // FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, alpha11 );
          bl1_cscals( epsilon1, alpha11 );
        }

        // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha11, delta1 );
        // FLA_Copyt( FLA_NO_TRANSPOSE, alpha11, absv );
        // FLA_Absolute_value( absv );
        // FLA_Inv_scal( absv, delta1 );
        bl1_ccopys( BLIS1_CONJUGATE, alpha11, delta1 );
        bl1_cabsval2( alpha11, &absv );
        bl1_cinvscals( &absv, delta1 );

        // FLA_Scalc( FLA_NO_CONJUGATE, delta1, alpha11 );
        // FLA_Obj_set_imag_part( FLA_ZERO, alpha11 );
        bl1_cscals( delta1, alpha11 );
        alpha11->imag = *buff_0;

        if ( n_ahead > 0 )
        {
          scomplex* a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
          // FLA_Scalc( FLA_NO_CONJUGATE, delta1, a12t_l );
          bl1_cscals( delta1, a12t_l );
        }
      }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_d = FLA_DOUBLE_COMPLEX_PTR( d );
      dcomplex* buff_e = FLA_DOUBLE_COMPLEX_PTR( e );
      dcomplex* buff_1 = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
      double*   buff_0 = FLA_DOUBLE_PTR( FLA_ZERO );

      for ( i = 0; i < min_m_n; ++i )
      {
        dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
        dcomplex* delta1   = buff_d + (i  )*inc_d;
        dcomplex* epsilon1 = buff_e + (i  )*inc_e;
        dcomplex  absv;

        int       n_ahead  = n_A - i - 1;
        int       n_behind = i;

        if ( n_behind == 0 )
        {
          // FLA_Set( FLA_ONE, epsilon1 );
          *epsilon1 = *buff_1;
        }
        else
        {
          dcomplex* a01_b    = buff_A + (i  )*cs_A + (i-1)*rs_A;
          // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a01_b, epsilon1 );
          // FLA_Copyt( FLA_NO_TRANSPOSE, a01_b, absv );
          // FLA_Absolute_value( absv );
          // FLA_Inv_scal( absv, epsilon1 );
          bl1_zcopys( BLIS1_CONJUGATE, a01_b, epsilon1 );
          bl1_zabsval2( a01_b, &absv );
          bl1_zinvscals( &absv, epsilon1 );

          // FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, a01_b );
          // FLA_Obj_set_imag_part( FLA_ZERO, a01_b );
          bl1_zscals( epsilon1, a01_b );
          a01_b->imag = *buff_0;

          // FLA_Scalc( FLA_NO_CONJUGATE, epsilon1, alpha11 );
          bl1_zscals( epsilon1, alpha11 );
        }

        // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, alpha11, delta1 );
        // FLA_Copyt( FLA_NO_TRANSPOSE, alpha11, absv );
        // FLA_Absolute_value( absv );
        // FLA_Inv_scal( absv, delta1 );
        bl1_zcopys( BLIS1_CONJUGATE, alpha11, delta1 );
        bl1_zabsval2( alpha11, &absv );
        bl1_zinvscals( &absv, delta1 );

        // FLA_Scalc( FLA_NO_CONJUGATE, delta1, alpha11 );
        // FLA_Obj_set_imag_part( FLA_ZERO, alpha11 );
        bl1_zscals( delta1, alpha11 );
        alpha11->imag = *buff_0;

        if ( n_ahead > 0 )
        {
          dcomplex* a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
          // FLA_Scalc( FLA_NO_CONJUGATE, delta1, a12t_l );
          bl1_zscals( delta1, a12t_l );
        }
      }

      break;
    }
  }

  return FLA_SUCCESS;
}

