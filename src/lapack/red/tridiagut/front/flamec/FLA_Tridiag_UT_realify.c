
#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_realify( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d )
{
  FLA_Error r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Tridiag_UT_realify_check( uplo, A, d );

  if ( FLA_Obj_is_real( A ) )
  {
    FLA_Set( FLA_ONE, d );
    return FLA_SUCCESS;
  }

  if ( uplo == FLA_LOWER_TRIANGULAR )
    //r_val = FLA_Tridiag_UT_l_realify_unb( A, d );
    r_val = FLA_Tridiag_UT_l_realify_opt( A, d );
  else
    //r_val = FLA_Tridiag_UT_u_realify_unb( A, d );
    r_val = FLA_Tridiag_UT_u_realify_opt( A, d );

  return r_val;
}


FLA_Error FLA_Tridiag_UT_l_realify_unb( FLA_Obj A, FLA_Obj d )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj dT,              d0,
          dB,              delta1,
                           d2;

  FLA_Obj a10t_l, a10t_r;

  FLA_Obj a21_t,
          a21_b;

  FLA_Obj absv;


  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, &absv );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     1, 1, FLA_TL );

  FLA_Part_2x1( d,    &dT,
                      &dB,            1, FLA_TOP );

  // Set first element of vector d to one.
  FLA_Set( FLA_ONE, dT );

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

    /*------------------------------------------------------------*/

    FLA_Part_1x2( a10t,   &a10t_l, &a10t_r,    1, FLA_RIGHT );

    FLA_Part_2x1( a21,    &a21_t,
                          &a21_b,              1, FLA_TOP );

    // delta1 = conj(a10t_r) / abs(a10t_r);
    FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a10t_r, delta1 );
    FLA_Copyt( FLA_NO_TRANSPOSE, a10t_r, absv );
    FLA_Absolute_value( absv );
    FLA_Inv_scal( absv, delta1 );

    // a10t_r  = delta1 * a10t_r;
    //         = abs(a10t_r);
    // alpha11 = delta1 * alpha11 * conj(delta1);
    //         = alpha11;
    // a21_t   = a21_t * conj(delta1);
    FLA_Copyt( FLA_NO_TRANSPOSE, absv, a10t_r );
    FLA_Scalc( FLA_CONJUGATE, delta1, a21_t );

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
  }

  FLA_Obj_free( &absv );

  return FLA_SUCCESS;
}

FLA_Error FLA_Tridiag_UT_l_realify_opt( FLA_Obj A, FLA_Obj d )
{
  FLA_Datatype datatype;
  int          m_A;
  int          rs_A, cs_A;
  int          inc_d;
  int          i;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_d    = FLA_Obj_vector_inc( d );


  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_d = FLA_FLOAT_PTR( d );
      float* buff_1 = FLA_FLOAT_PTR( FLA_ONE );

      bl1_ssetv( m_A,
                 buff_1,
                 buff_d, inc_d );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_d = FLA_DOUBLE_PTR( d );
      double* buff_1 = FLA_DOUBLE_PTR( FLA_ONE );

      bl1_dsetv( m_A,
                 buff_1,
                 buff_d, inc_d );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_d = FLA_COMPLEX_PTR( d );
      scomplex* buff_1 = FLA_COMPLEX_PTR( FLA_ONE );

      bl1_csetv( 1,
                 buff_1,
                 buff_d, inc_d );

      for ( i = 1; i < m_A; ++i )
      {
        scomplex* a10t_r   = buff_A + (i-1)*cs_A + (i  )*rs_A;
        scomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        scomplex* delta1   = buff_d + (i  )*inc_d;
        scomplex  absv;
        scomplex  conj_delta1;

        int       m_ahead  = m_A - i - 1;

        // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a10t_r, delta1 );
        // FLA_Copyt( FLA_NO_TRANSPOSE, a10t_r, absv );
        // FLA_Absolute_value( absv );
        // FLA_Inv_scal( absv, delta1 );
        bl1_ccopys( BLIS1_CONJUGATE, a10t_r, delta1 );
        bl1_cabsval2( a10t_r, &absv );
        bl1_cinvscals( &absv, delta1 );

        // FLA_Copyt( FLA_NO_TRANSPOSE, absv, a10t_r );
        // FLA_Scalc( FLA_CONJUGATE, delta1, a21_t );
        *a10t_r = absv;
        if ( m_ahead > 0 )
        {
          bl1_ccopyconj( delta1, &conj_delta1 );
          bl1_cscals( &conj_delta1, a21_t );
        }
      }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_d = FLA_DOUBLE_COMPLEX_PTR( d );
      dcomplex* buff_1 = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );

      bl1_zsetv( 1,
                 buff_1,
                 buff_d, inc_d );

      for ( i = 1; i < m_A; ++i )
      {
        dcomplex* a10t_r   = buff_A + (i-1)*cs_A + (i  )*rs_A;
        dcomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
        dcomplex* delta1   = buff_d + (i  )*inc_d;
        dcomplex  absv;
        dcomplex  conj_delta1;

        int       m_ahead  = m_A - i - 1;

        // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a10t_r, delta1 );
        // FLA_Copyt( FLA_NO_TRANSPOSE, a10t_r, absv );
        // FLA_Absolute_value( absv );
        // FLA_Inv_scal( absv, delta1 );
        bl1_zcopys( BLIS1_CONJUGATE, a10t_r, delta1 );
        bl1_zabsval2( a10t_r, &absv );
        bl1_zinvscals( &absv, delta1 );

        // FLA_Copyt( FLA_NO_TRANSPOSE, absv, a10t_r );
        // FLA_Scalc( FLA_CONJUGATE, delta1, a21_t );
        *a10t_r = absv;
        if ( m_ahead > 0 )
        {
          bl1_zcopyconj( delta1, &conj_delta1 );
          bl1_zscals( &conj_delta1, a21_t );
        }
      }

      break;
    }
  }

  return FLA_SUCCESS;
}




FLA_Error FLA_Tridiag_UT_u_realify_unb( FLA_Obj A, FLA_Obj d )
{
  FLA_Obj ATL,   ATR,      A00,  a01,     A02, 
          ABL,   ABR,      a10t, alpha11, a12t,
                           A20,  a21,     A22;

  FLA_Obj dT,              d0,
          dB,              delta1,
                           d2;

  FLA_Obj a01_t,
          a01_b;

  FLA_Obj a12t_l, a12t_r;


  FLA_Obj absv;


  FLA_Obj_create( FLA_Obj_datatype( A ), 1, 1, 0, 0, &absv );

  FLA_Part_2x2( A,    &ATL, &ATR,
                      &ABL, &ABR,     1, 1, FLA_TL );

  FLA_Part_2x1( d,    &dT,
                      &dB,            1, FLA_TOP );

  // Set first element of vector d to one.
  FLA_Set( FLA_ONE, dT );

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

    /*------------------------------------------------------------*/

    FLA_Part_2x1( a01,    &a01_t,
                          &a01_b,              1, FLA_BOTTOM );

    FLA_Part_1x2( a12t,   &a12t_l, &a12t_r,    1, FLA_LEFT );

    // delta1 = conj(a01_b) / abs(a01_b);
    FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a01_b, delta1 );
    FLA_Copyt( FLA_NO_TRANSPOSE, a01_b, absv );
    FLA_Absolute_value( absv );
    FLA_Inv_scal( absv, delta1 );

    // a01_b   = delta1 * a01_b;
    //         = abs(a01_b);
    // alpha11 = delta1 * alpha11 * conj(delta1);
    //         = alpha11;
    // a12t_l  = a12t_l * conj(delta1);
    FLA_Copyt( FLA_NO_TRANSPOSE, absv, a01_b );
    FLA_Scalc( FLA_CONJUGATE, delta1, a12t_l );

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
  }

  FLA_Obj_free( &absv );

  return FLA_SUCCESS;
}

FLA_Error FLA_Tridiag_UT_u_realify_opt( FLA_Obj A, FLA_Obj d )
{
  FLA_Datatype datatype;
  int          m_A;
  int          rs_A, cs_A;
  int          inc_d;
  int          i;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_d    = FLA_Obj_vector_inc( d );


  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_d = FLA_FLOAT_PTR( d );
      float* buff_1 = FLA_FLOAT_PTR( FLA_ONE );

      bl1_ssetv( m_A,
                 buff_1,
                 buff_d, inc_d );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_d = FLA_DOUBLE_PTR( d );
      double* buff_1 = FLA_DOUBLE_PTR( FLA_ONE );

      bl1_dsetv( m_A,
                 buff_1,
                 buff_d, inc_d );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_d = FLA_COMPLEX_PTR( d );
      scomplex* buff_1 = FLA_COMPLEX_PTR( FLA_ONE );

      bl1_csetv( 1,
                 buff_1,
                 buff_d, inc_d );

      for ( i = 1; i < m_A; ++i )
      {
        scomplex* a01_b    = buff_A + (i  )*cs_A + (i-1)*rs_A;
        scomplex* a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        scomplex* delta1   = buff_d + (i  )*inc_d;
        scomplex  absv;
        scomplex  conj_delta1;

        int       m_ahead  = m_A - i - 1;

        // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a01_b, delta1 );
        // FLA_Copyt( FLA_NO_TRANSPOSE, a01_b, absv );
        // FLA_Absolute_value( absv );
        // FLA_Inv_scal( absv, delta1 );
        bl1_ccopys( BLIS1_CONJUGATE, a01_b, delta1 );
        bl1_cabsval2( a01_b, &absv );
        bl1_cinvscals( &absv, delta1 );

        // FLA_Copyt( FLA_NO_TRANSPOSE, absv, a01_b );
        // FLA_Scalc( FLA_CONJUGATE, delta1, a12t_l );
        *a01_b = absv;
        if ( m_ahead > 0 )
        {
          bl1_ccopyconj( delta1, &conj_delta1 );
          bl1_cscals( &conj_delta1, a12t_l );
        }
      }

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_d = FLA_DOUBLE_COMPLEX_PTR( d );
      dcomplex* buff_1 = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );

      bl1_zsetv( 1,
                 buff_1,
                 buff_d, inc_d );

      for ( i = 1; i < m_A; ++i )
      {
        dcomplex* a01_b    = buff_A + (i  )*cs_A + (i-1)*rs_A;
        dcomplex* a12t_l   = buff_A + (i+1)*cs_A + (i  )*rs_A;
        dcomplex* delta1   = buff_d + (i  )*inc_d;
        dcomplex  absv;
        dcomplex  conj_delta1;

        int       m_ahead  = m_A - i - 1;

        // FLA_Copyt( FLA_CONJ_NO_TRANSPOSE, a01_b, delta1 );
        // FLA_Copyt( FLA_NO_TRANSPOSE, a01_b, absv );
        // FLA_Absolute_value( absv );
        // FLA_Inv_scal( absv, delta1 );
        bl1_zcopys( BLIS1_CONJUGATE, a01_b, delta1 );
        bl1_zabsval2( a01_b, &absv );
        bl1_zinvscals( &absv, delta1 );

        // FLA_Copyt( FLA_NO_TRANSPOSE, absv, a01_b );
        // FLA_Scalc( FLA_CONJUGATE, delta1, a12t_l );
        *a01_b = absv;
        if ( m_ahead > 0 )
        {
          bl1_zcopyconj( delta1, &conj_delta1 );
          bl1_zscals( &conj_delta1, a12t_l );
        }
      }

      break;
    }
  }

  return FLA_SUCCESS;
}




