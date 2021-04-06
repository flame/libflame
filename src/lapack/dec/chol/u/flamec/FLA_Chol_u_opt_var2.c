/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Chol_u_opt_var2( FLA_Obj A )
{
  FLA_Error    r_val = FLA_SUCCESS;
  FLA_Datatype datatype;
  integer          mn_A;
  integer          rs_A, cs_A;

  datatype = FLA_Obj_datatype( A );

  mn_A     = FLA_Obj_length( A );
  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );

      r_val = FLA_Chol_u_ops_var2( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );

      r_val = FLA_Chol_u_opd_var2( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );

      r_val = FLA_Chol_u_opc_var2( mn_A,
                                   buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );

      r_val = FLA_Chol_u_opz_var2( mn_A,
                                   buff_A, rs_A, cs_A );
      
      break;
    }
  }

  return r_val;
}



FLA_Error FLA_Chol_u_ops_var2( integer mn_A,
                               float* buff_A, integer rs_A, integer cs_A )
{
  float*    buff_1  = FLA_FLOAT_PTR( FLA_ONE );
  float*    buff_m1 = FLA_FLOAT_PTR( FLA_MINUS_ONE );
  integer       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    float*    a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    float*    A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    float*    a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    integer       mn_ahead  = mn_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, a01, FLA_ONE, alpha11 );
    bl1_sdots( BLIS1_CONJUGATE,
               mn_behind,
               buff_m1,
               a01, rs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, a01, FLA_ONE, a12t );
    bl1_sgemv( BLIS1_TRANSPOSE,
               BLIS1_CONJUGATE,
               mn_behind,
               mn_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a12t, cs_A );

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_ssqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a12t );
    bl1_sinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_u_opd_var2( integer mn_A,
                               double* buff_A, integer rs_A, integer cs_A )
{
  double*   buff_1  = FLA_DOUBLE_PTR( FLA_ONE );
  double*   buff_m1 = FLA_DOUBLE_PTR( FLA_MINUS_ONE );
  integer       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    double*   a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    double*   A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    double*   a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    integer       mn_ahead  = mn_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, a01, FLA_ONE, alpha11 );
    bl1_ddots( BLIS1_CONJUGATE,
               mn_behind,
               buff_m1,
               a01, rs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, a01, FLA_ONE, a12t );
    bl1_dgemv( BLIS1_TRANSPOSE,
               BLIS1_CONJUGATE,
               mn_behind,
               mn_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a12t, cs_A );

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_dsqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a12t );
    bl1_dinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_u_opc_var2( integer mn_A,
                               scomplex* buff_A, integer rs_A, integer cs_A )
{
  scomplex* buff_1  = FLA_COMPLEX_PTR( FLA_ONE );
  scomplex* buff_m1 = FLA_COMPLEX_PTR( FLA_MINUS_ONE );
  integer       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    scomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    scomplex* A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    scomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    integer       mn_ahead  = mn_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, a01, FLA_ONE, alpha11 );
    bl1_cdots( BLIS1_CONJUGATE,
               mn_behind,
               buff_m1,
               a01, rs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, a01, FLA_ONE, a12t );
    bl1_cgemv( BLIS1_TRANSPOSE,
               BLIS1_CONJUGATE,
               mn_behind,
               mn_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a12t, cs_A );

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_csqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a12t );
    bl1_cinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Chol_u_opz_var2( integer mn_A,
                               dcomplex* buff_A, integer rs_A, integer cs_A )
{
  dcomplex* buff_1  = FLA_DOUBLE_COMPLEX_PTR( FLA_ONE );
  dcomplex* buff_m1 = FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
  integer       i;
  FLA_Error e_val;

  for ( i = 0; i < mn_A; ++i )
  {
    dcomplex* a01       = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* alpha11   = buff_A + (i  )*cs_A + (i  )*rs_A;
    dcomplex* A02       = buff_A + (i+1)*cs_A + (0  )*rs_A;
    dcomplex* a12t      = buff_A + (i+1)*cs_A + (i  )*rs_A;

    integer       mn_ahead  = mn_A - i - 1;
    integer       mn_behind = i;

    /*------------------------------------------------------------*/

    // FLA_Dotcs_external( FLA_CONJUGATE, FLA_MINUS_ONE, a01, a01, FLA_ONE, alpha11 );
    bl1_zdots( BLIS1_CONJUGATE,
               mn_behind,
               buff_m1,
               a01, rs_A,
               a01, rs_A,
               buff_1,
               alpha11 );

    // FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_MINUS_ONE, A02, a01, FLA_ONE, a12t );
    bl1_zgemv( BLIS1_TRANSPOSE,
               BLIS1_CONJUGATE,
               mn_behind,
               mn_ahead,
               buff_m1,
               A02, rs_A, cs_A,
               a01, rs_A,
               buff_1,
               a12t, cs_A );

    // r_val = FLA_Sqrt( alpha11 );
    // if ( r_val != FLA_SUCCESS )
    //   return ( FLA_Obj_length( A00 ) + 1 );
    bl1_zsqrte( alpha11, &e_val );
    if ( e_val != FLA_SUCCESS ) return mn_behind;

    // FLA_Inv_scal_external( alpha11, a12t );
    bl1_zinvscalv( BLIS1_NO_CONJUGATE,
                   mn_ahead,
                   alpha11,
                   a12t, cs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

