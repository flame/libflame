/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_NON_CRITICAL_CODE

FLA_Error FLA_LU_piv_opt_var4( FLA_Obj A, FLA_Obj p )
{
  FLA_Error    r_val = FLA_SUCCESS;
  FLA_Datatype datatype;
  integer          m_A, n_A;
  integer          rs_A, cs_A;
  integer          inc_p;

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
      integer*   buff_p = FLA_INT_PTR( p );

      r_val = FLA_LU_piv_ops_var4( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_p, inc_p );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      integer*    buff_p = FLA_INT_PTR( p );

      r_val = FLA_LU_piv_opd_var4( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_p, inc_p );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      integer*      buff_p = FLA_INT_PTR( p );

      r_val = FLA_LU_piv_opc_var4( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_p, inc_p );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      integer*      buff_p = FLA_INT_PTR( p );

      r_val = FLA_LU_piv_opz_var4( m_A,
                                   n_A,
                                   buff_A, rs_A, cs_A,
                                   buff_p, inc_p );

      break;
    }
  }

  return r_val;
}



FLA_Error FLA_LU_piv_ops_var4( integer m_A,
                               integer n_A,
                               float*    buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p )
{
  FLA_Error r_val = FLA_SUCCESS;
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  integer       min_m_n = min( m_A, n_A );
  integer       i, is_null_pivot;


  for ( i = 0; i < min_m_n; ++i )
  {
    float     pivot_val = fzero;
    float*    a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    float*    A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    float*    a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    float*    A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    integer*      pi1       = buff_p + i*inc_p;

    integer       m_ahead   = m_A - i - 1;
    integer       n_ahead   = n_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

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

    is_null_pivot = (pivot_val == fzero);
    if ( is_null_pivot ) 
    {
        r_val = ( r_val == FLA_SUCCESS ? i : r_val );
    }
    else
    {
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
      FLA_Apply_pivots_ln_ops_var1( 1,
                                    alpha11, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );

      // FLA_Merge_2x1( a10t,
      //                A20,      &AB0 );

      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
      FLA_Apply_pivots_ln_ops_var1( mn_behind,
                                    a10t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Merge_2x1( a12t,
      //                A22,      &AB2 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB2 );
      FLA_Apply_pivots_ln_ops_var1( n_ahead,
                                    a12t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
    }

    // FLA_Gemv_external( FLA_TRANSPOSE, FLA_MINUS_ONE, A02, a10t, FLA_ONE, a12t );
    bl1_sgemv( BLIS1_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               mn_behind,
               n_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );

    if ( ! is_null_pivot )
    {
      // FLA_Inv_scal_external( alpha11, a21 );
      bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     alpha11,
                     a21, rs_A );
    }
    /*------------------------------------------------------------*/

  }

  return r_val;
}



FLA_Error FLA_LU_piv_opd_var4( integer m_A,
                               integer n_A,
                               double*   buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p )
{
  FLA_Error r_val   = FLA_SUCCESS;
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  integer       min_m_n = min( m_A, n_A );
  integer       i, is_null_pivot;

  for ( i = 0; i < min_m_n; ++i )
  {
    double    pivot_val = dzero;
    double*   a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    double*   A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    double*   a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    double*   A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    integer*      pi1       = buff_p + i*inc_p;

    integer       m_ahead   = m_A - i - 1;
    integer       n_ahead   = n_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

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

    is_null_pivot = (pivot_val == dzero);
    if ( is_null_pivot ) 
    {
      r_val = ( r_val == FLA_SUCCESS ? i : r_val );
    }
    else
    {
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
      FLA_Apply_pivots_ln_opd_var1( 1,
                                    alpha11, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Merge_2x1( a10t,
      //                A20,      &AB0 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
      FLA_Apply_pivots_ln_opd_var1( mn_behind,
                                    a10t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Merge_2x1( a12t,
      //                A22,      &AB2 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB2 );
      FLA_Apply_pivots_ln_opd_var1( n_ahead,
                                    a12t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
    }

    // FLA_Gemv_external( FLA_TRANSPOSE, FLA_MINUS_ONE, A02, a10t, FLA_ONE, a12t );
    bl1_dgemv( BLIS1_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               mn_behind,
               n_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );

    if ( ! is_null_pivot )
    {
      // FLA_Inv_scal_external( alpha11, a21 );
      bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     alpha11,
                     a21, rs_A );
    }
    /*------------------------------------------------------------*/

  }

  return r_val;
}



FLA_Error FLA_LU_piv_opc_var4( integer m_A,
                               integer n_A,
                               scomplex* buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p )
{
  FLA_Error r_val   = FLA_SUCCESS;
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  integer       min_m_n = min( m_A, n_A );
  integer       i, is_null_pivot;

  for ( i = 0; i < min_m_n; ++i )
  {
    scomplex  pivot_val = czero;
    scomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    scomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    scomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    scomplex* A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    integer*      pi1       = buff_p + i*inc_p;

    integer       m_ahead   = m_A - i - 1;
    integer       n_ahead   = n_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

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

    is_null_pivot = (pivot_val.real == czero.real && pivot_val.imag == czero.imag);
    if ( is_null_pivot ) 
    {
      r_val = ( r_val == FLA_SUCCESS ? i : r_val );
    }
    else
    {
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
      FLA_Apply_pivots_ln_opc_var1( 1,
                                    alpha11, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Merge_2x1( a10t,
      //                A20,      &AB0 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
      FLA_Apply_pivots_ln_opc_var1( mn_behind,
                                    a10t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Merge_2x1( a12t,
      //                A22,      &AB2 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB2 );
      FLA_Apply_pivots_ln_opc_var1( n_ahead,
                                    a12t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
    }

    // FLA_Gemv_external( FLA_TRANSPOSE, FLA_MINUS_ONE, A02, a10t, FLA_ONE, a12t );
    bl1_cgemv( BLIS1_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               mn_behind,
               n_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );

    if ( ! is_null_pivot )
    {
      // FLA_Inv_scal_external( alpha11, a21 );
      bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     alpha11,
                     a21, rs_A );
    }
    /*------------------------------------------------------------*/

  }

  return r_val;
}



FLA_Error FLA_LU_piv_opz_var4( integer m_A,
                               integer n_A,
                               dcomplex* buff_A, integer rs_A, integer cs_A,
                               integer*      buff_p, integer inc_p )
{
  FLA_Error r_val   = FLA_SUCCESS;
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  integer       min_m_n = min( m_A, n_A );
  integer       i, is_null_pivot;

  for ( i = 0; i < min_m_n; ++i )
  {
    dcomplex  pivot_val = zzero;
    dcomplex* a10t      = buff_A + (0  )*cs_A + (i  )*rs_A;
    dcomplex* A20       = buff_A + (0  )*cs_A + (i+1)*rs_A;
    dcomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* a21       = buff_A + (i  )*cs_A + (i+1)*rs_A;
    dcomplex* A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    integer*      pi1       = buff_p + i*inc_p;

    integer       m_ahead   = m_A - i - 1;
    integer       n_ahead   = n_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

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

    is_null_pivot = (pivot_val.real == zzero.real && pivot_val.imag == zzero.imag);
    if ( is_null_pivot ) 
    {
      r_val = ( r_val == FLA_SUCCESS ? i : r_val );
    }
    else
    {
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, aB1 );
      FLA_Apply_pivots_ln_opz_var1( 1,
                                    alpha11, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Merge_2x1( a10t,
      //                A20,      &AB0 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB0 );
      FLA_Apply_pivots_ln_opz_var1( mn_behind,
                                    a10t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
      
      // FLA_Merge_2x1( a12t,
      //                A22,      &AB2 );
      
      // FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE, pi1, AB2 );
      FLA_Apply_pivots_ln_opz_var1( n_ahead,
                                    a12t, rs_A, cs_A,
                                    0,
                                    0,
                                    pi1, inc_p );
    }

    // FLA_Gemv_external( FLA_TRANSPOSE, FLA_MINUS_ONE, A02, a10t, FLA_ONE, a12t );
    bl1_zgemv( BLIS1_TRANSPOSE,
               BLIS1_NO_CONJUGATE,
               mn_behind,
               n_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a10t, cs_A,
               buff_1,
               a12t, cs_A );

    if ( ! is_null_pivot ) 
    {
      // FLA_Inv_scal_external( alpha11, a21 );
      bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                     m_ahead,
                     alpha11,
                     a21, rs_A );
    }
    /*------------------------------------------------------------*/

  }

  return r_val;
}

#endif
