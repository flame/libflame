
#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_LU_piv_opt_var3( FLA_Obj A, FLA_Obj p )
{
  FLA_Error    r_val = FLA_SUCCESS;
  FLA_Datatype datatype;
  int          m_A, n_A;
  int          rs_A, cs_A;
  int          inc_p;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_p    = FLA_Obj_vector_inc( p );


  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      int*   buff_p = FLA_INT_PTR( p );

      r_val = FLA_LU_piv_ops_var3( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_p, inc_p );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      int*    buff_p = FLA_INT_PTR( p );

      r_val = FLA_LU_piv_opd_var3( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_p, inc_p );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      int*      buff_p = FLA_INT_PTR( p );
      
      r_val = FLA_LU_piv_opc_var3( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_p, inc_p );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      int*      buff_p = FLA_INT_PTR( p );

      r_val = FLA_LU_piv_opz_var3( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_p, inc_p );

      break;
    }
  }

  return r_val;
}



FLA_Error FLA_LU_piv_ops_var3( int m_A,
                               int n_A,
                               float*    buff_A, int rs_A, int cs_A,
                               int*      buff_p, int inc_p )
{
  FLA_Error r_val   = FLA_SUCCESS;
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    float     pivot_val = fzero;
    float*    A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    float*    a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    float*    a1        = buff_A + (i  )*cs_A + (0  )*rs_A;

    int*      p0        = buff_p;
    int*      pi1       = buff_p + i*inc_p;

    int       m_ahead   = m_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p0, a1 );
    FLA_Apply_pivots_ln_ops_var1( 1,
                                  a1, rs_A, cs_A,
                                  0,
                                  mn_behind - 1,
                                  p0, inc_p );

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A00, a01 );
    bl1_strsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bl1_sdots( BLIS1_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, a01, FLA_ONE, a21 );
    bl1_sgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               mn_behind,
               buff_m1,
               A20, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a21, rs_A );

    // FLA_Merge_2x1( alpha11,
    //                    a21,      &aB1 );

    // FLA_Amax_external( aB1, pi1 );
    bl1_samax( m_ahead + 1,
               alpha11, rs_A,
               pi1 );

    // If a null pivot is encountered, return the index. 
    pivot_val = *(alpha11 + *pi1);
    if ( pivot_val == fzero ) r_val = ( r_val == FLA_SUCCESS ? i : r_val );
    else 
    {
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
      FLA_Apply_pivots_ln_ops_var1( 1,
                                    alpha11, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Inv_scal_external( alpha11, a21 );
      bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     alpha11,
                     a21, rs_A );
      
      // FLA_Merge_2x1( a10t,
      //                A20,      &AB0 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
      FLA_Apply_pivots_ln_ops_var1( mn_behind,
                                    a10t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
    }
    /*------------------------------------------------------------*/

  }

  if ( m_A < n_A )
  {
    float*    ATL = buff_A;
    float*    ATR = buff_A + m_A*cs_A;

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, ATR );
    FLA_Apply_pivots_ln_ops_var1( n_A - m_A,
                                  ATR, rs_A, cs_A,
                                  0,
                                  m_A - 1,
                                  buff_p, inc_p );

    // FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
    //                    FLA_ONE, ATL, ATR );
    bl1_strsm( BLIS1_LEFT,
               BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               m_A,
               n_A - m_A,
               buff_1,
               ATL, rs_A, cs_A,
               ATR, rs_A, cs_A );
  }

  return r_val;
}



FLA_Error FLA_LU_piv_opd_var3( int m_A,
                               int n_A,
                               double*   buff_A, int rs_A, int cs_A,
                               int*      buff_p, int inc_p )
{
  FLA_Error r_val   = FLA_SUCCESS;
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    double    pivot_val = dzero;
    double*   A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    double*   a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    double*   a1        = buff_A + (i  )*cs_A + (0  )*rs_A;

    int*      p0        = buff_p;
    int*      pi1       = buff_p + i*inc_p;

    int       m_ahead   = m_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p0, a1 );
    FLA_Apply_pivots_ln_opd_var1( 1,
                                  a1, rs_A, cs_A,
                                  0,
                                  mn_behind - 1,
                                  p0, inc_p );

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A00, a01 );
    bl1_dtrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bl1_ddots( BLIS1_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, a01, FLA_ONE, a21 );
    bl1_dgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               mn_behind,
               buff_m1,
               A20, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a21, rs_A );

    // FLA_Merge_2x1( alpha11,
    //                    a21,      &aB1 );

    // FLA_Amax_external( aB1, pi1 );
    bl1_damax( m_ahead + 1,
               alpha11, rs_A,
               pi1 );

    // If a null pivot is encountered, return the index. 
    pivot_val =*(alpha11 + *pi1);
    if ( pivot_val == dzero )  r_val = ( r_val == FLA_SUCCESS ? i : r_val );
    else
    {
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
      FLA_Apply_pivots_ln_opd_var1( 1,
                                    alpha11, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Inv_scal_external( alpha11, a21 );
      bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     alpha11,
                     a21, rs_A );
      
      // FLA_Merge_2x1( a10t,
      //                A20,      &AB0 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
      FLA_Apply_pivots_ln_opd_var1( mn_behind,
                                    a10t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
    }
    /*------------------------------------------------------------*/

  }

  if ( m_A < n_A )
  {
    double*   ATL = buff_A;
    double*   ATR = buff_A + m_A*cs_A;

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, ATR );
    FLA_Apply_pivots_ln_opd_var1( n_A - m_A,
                                  ATR, rs_A, cs_A,
                                  0,
                                  m_A - 1,
                                  buff_p, inc_p );

    // FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
    //                    FLA_ONE, ATL, ATR );
    bl1_dtrsm( BLIS1_LEFT,
               BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               m_A,
               n_A - m_A,
               buff_1,
               ATL, rs_A, cs_A,
               ATR, rs_A, cs_A );
  }

  return r_val;
}



FLA_Error FLA_LU_piv_opc_var3( int m_A,
                               int n_A,
                               scomplex* buff_A, int rs_A, int cs_A,
                               int*      buff_p, int inc_p )
{
  FLA_Error r_val   = FLA_SUCCESS;
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    scomplex  pivot_val = czero;
    scomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    scomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    scomplex* a1        = buff_A + (i  )*cs_A + (0  )*rs_A;

    int*      p0        = buff_p;
    int*      pi1       = buff_p + i*inc_p;

    int       m_ahead   = m_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p0, a1 );
    FLA_Apply_pivots_ln_opc_var1( 1,
                                  a1, rs_A, cs_A,
                                  0,
                                  mn_behind - 1,
                                  p0, inc_p );

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A00, a01 );
    bl1_ctrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bl1_cdots( BLIS1_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, a01, FLA_ONE, a21 );
    bl1_cgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               mn_behind,
               buff_m1,
               A20, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a21, rs_A );

    // FLA_Merge_2x1( alpha11,
    //                    a21,      &aB1 );

    // FLA_Amax_external( aB1, pi1 );
    bl1_camax( m_ahead + 1,
               alpha11, rs_A,
               pi1 );

    // If a null pivot is encountered, return the index. 
    pivot_val =*(alpha11 + *pi1);
    if ( pivot_val.real == czero.real && 
         pivot_val.imag == czero.imag )  r_val = ( r_val == FLA_SUCCESS ? i : r_val );
    else
    {
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
      FLA_Apply_pivots_ln_opc_var1( 1,
                                    alpha11, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Inv_scal_external( alpha11, a21 );
      bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     alpha11,
                     a21, rs_A );
      
      // FLA_Merge_2x1( a10t,
      //                A20,      &AB0 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
      FLA_Apply_pivots_ln_opc_var1( mn_behind,
                                    a10t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
    }
    /*------------------------------------------------------------*/

  }

  if ( m_A < n_A )
  {
    scomplex* ATL = buff_A;
    scomplex* ATR = buff_A + m_A*cs_A;

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, ATR );
    FLA_Apply_pivots_ln_opc_var1( n_A - m_A,
                                  ATR, rs_A, cs_A,
                                  0,
                                  m_A - 1,
                                  buff_p, inc_p );

    // FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
    //                    FLA_ONE, ATL, ATR );
    bl1_ctrsm( BLIS1_LEFT,
               BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               m_A,
               n_A - m_A,
               buff_1,
               ATL, rs_A, cs_A,
               ATR, rs_A, cs_A );
  }

  return r_val;
}



FLA_Error FLA_LU_piv_opz_var3( int m_A,
                               int n_A,
                               dcomplex* buff_A, int rs_A, int cs_A,
                               int*      buff_p, int inc_p )
{
  FLA_Error r_val   = FLA_SUCCESS;
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  int       min_m_n = min( m_A, n_A );
  int       i;

  for ( i = 0; i < min_m_n; ++i )
  {
    dcomplex  pivot_val = zzero;
    dcomplex* A00       = buff_A + (0  )*cs_A + (0  )*rs_A;
    dcomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;

    dcomplex* a1        = buff_A + (i  )*cs_A + (0  )*rs_A;

    int*      p0        = buff_p;
    int*      pi1       = buff_p + i*inc_p;

    int       m_ahead   = m_A - i - 1;
    int       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p0, a1 );
    FLA_Apply_pivots_ln_opz_var1( 1,
                                  a1, rs_A, cs_A,
                                  0,
                                  mn_behind - 1,
                                  p0, inc_p );

    // FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG, A00, a01 );
    bl1_ztrsv( BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               mn_behind,
               A00, rs_A, cs_A,
               a01, rs_A );

    // FLA_Dots_external( FLA_MINUS_ONE, a10t, a01, FLA_ONE, alpha11 );
    bl1_zdots( BLIS1_NO_CONJUGATE,
               mn_behind,
               buff_m1,
               a10t, cs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemv_external( FLA_NO_TRANSPOSE, FLA_MINUS_ONE, A20, a01, FLA_ONE, a21 );
    bl1_zgemv( BLIS1_NO_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               m_ahead,
               mn_behind,
               buff_m1,
               A20, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a21, rs_A );

    // FLA_Merge_2x1( alpha11,
    //                    a21,      &aB1 );

    // FLA_Amax_external( aB1, pi1 );
    bl1_zamax( m_ahead + 1,
               alpha11, rs_A,
               pi1 );

    // If a null pivot is encountered, return the index. 
    pivot_val =*(alpha11 + *pi1);
    if ( pivot_val.real == zzero.real &&
         pivot_val.imag == zzero.imag )  r_val = ( r_val == FLA_SUCCESS ? i : r_val );
    else
    {
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
      FLA_Apply_pivots_ln_opz_var1( 1,
                                    alpha11, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Inv_scal_external( alpha11, a21 );
      bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     alpha11,
                     a21, rs_A );
      
      // FLA_Merge_2x1( a10t,
      //                A20,      &AB0 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
      FLA_Apply_pivots_ln_opz_var1( mn_behind,
                                    a10t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
    }
    /*------------------------------------------------------------*/

  }

  if ( m_A < n_A )
  {
    dcomplex* ATL = buff_A;
    dcomplex* ATR = buff_A + m_A*cs_A;

    // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, p, ATR );
    FLA_Apply_pivots_ln_opz_var1( n_A - m_A,
                                  ATR, rs_A, cs_A,
                                  0,
                                  m_A - 1,
                                  buff_p, inc_p );

    // FLA_Trsm_external( FLA_LEFT, FLA_LOWER_TRIANGULAR,
    //                    FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
    //                    FLA_ONE, ATL, ATR );
    bl1_ztrsm( BLIS1_LEFT,
               BLIS1_LOWER_TRIANGULAR,
               BLIS1_NO_TRANSPOSE,
               BLIS1_UNIT_DIAG,
               m_A,
               n_A - m_A,
               buff_1,
               ATL, rs_A, cs_A,
               ATR, rs_A, cs_A );
  }

  return r_val;
}

#endif
