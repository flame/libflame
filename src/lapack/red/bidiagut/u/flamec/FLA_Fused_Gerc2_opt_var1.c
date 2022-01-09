/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Fused_Gerc2_opt_var1( FLA_Obj alpha, FLA_Obj u, FLA_Obj y, FLA_Obj z, FLA_Obj v, FLA_Obj A )
{
/*
   Effective computation:
   A = A + alpha * ( u * y' + z * v' );
*/
  FLA_Datatype datatype;
  integer          m_A, n_A;
  integer          rs_A, cs_A;
  integer          inc_u, inc_y, inc_z, inc_v;

  datatype = FLA_Obj_datatype( A );

  m_A      = FLA_Obj_length( A );
  n_A      = FLA_Obj_width( A );

  rs_A     = FLA_Obj_row_stride( A );
  cs_A     = FLA_Obj_col_stride( A );

  inc_u    = FLA_Obj_vector_inc( u );
  inc_y    = FLA_Obj_vector_inc( y );
  inc_z    = FLA_Obj_vector_inc( z );
  inc_v    = FLA_Obj_vector_inc( v );
  

  switch ( datatype )
  {
    case FLA_FLOAT:
    {
      float* buff_A = FLA_FLOAT_PTR( A );
      float* buff_u = FLA_FLOAT_PTR( u );
      float* buff_y = FLA_FLOAT_PTR( y );
      float* buff_z = FLA_FLOAT_PTR( z );
      float* buff_v = FLA_FLOAT_PTR( v );
      float* buff_alpha = FLA_FLOAT_PTR( alpha );

      FLA_Fused_Gerc2_ops_var1( m_A,
                                n_A,
                                buff_alpha,
                                buff_u, inc_u,
                                buff_y, inc_y,
                                buff_z, inc_z,
                                buff_v, inc_v,
                                buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE:
    {
      double* buff_A = FLA_DOUBLE_PTR( A );
      double* buff_u = FLA_DOUBLE_PTR( u );
      double* buff_y = FLA_DOUBLE_PTR( y );
      double* buff_z = FLA_DOUBLE_PTR( z );
      double* buff_v = FLA_DOUBLE_PTR( v );
      double* buff_alpha = FLA_DOUBLE_PTR( alpha );

      FLA_Fused_Gerc2_opd_var1( m_A,
                                n_A,
                                buff_alpha,
                                buff_u, inc_u,
                                buff_y, inc_y,
                                buff_z, inc_z,
                                buff_v, inc_v,
                                buff_A, rs_A, cs_A );

      break;
    }

    case FLA_COMPLEX:
    {
      scomplex* buff_A = FLA_COMPLEX_PTR( A );
      scomplex* buff_u = FLA_COMPLEX_PTR( u );
      scomplex* buff_y = FLA_COMPLEX_PTR( y );
      scomplex* buff_z = FLA_COMPLEX_PTR( z );
      scomplex* buff_v = FLA_COMPLEX_PTR( v );
      scomplex* buff_alpha = FLA_COMPLEX_PTR( alpha );

      FLA_Fused_Gerc2_opc_var1( m_A,
                                n_A,
                                buff_alpha,
                                buff_u, inc_u,
                                buff_y, inc_y,
                                buff_z, inc_z,
                                buff_v, inc_v,
                                buff_A, rs_A, cs_A );

      break;
    }

    case FLA_DOUBLE_COMPLEX:
    {
      dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
      dcomplex* buff_u = FLA_DOUBLE_COMPLEX_PTR( u );
      dcomplex* buff_y = FLA_DOUBLE_COMPLEX_PTR( y );
      dcomplex* buff_z = FLA_DOUBLE_COMPLEX_PTR( z );
      dcomplex* buff_v = FLA_DOUBLE_COMPLEX_PTR( v );
      dcomplex* buff_alpha = FLA_DOUBLE_COMPLEX_PTR( alpha );

      FLA_Fused_Gerc2_opz_var1( m_A,
                                n_A,
                                buff_alpha,
                                buff_u, inc_u,
                                buff_y, inc_y,
                                buff_z, inc_z,
                                buff_v, inc_v,
                                buff_A, rs_A, cs_A );

      break;
    }
  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Gerc2_ops_var1( integer m_A,
                                    integer n_A,
                                    float* buff_alpha, 
                                    float* buff_u, integer inc_u, 
                                    float* buff_y, integer inc_y, 
                                    float* buff_z, integer inc_z, 
                                    float* buff_v, integer inc_v, 
                                    float* buff_A, integer rs_A, integer cs_A )
{
  integer       i;

  for ( i = 0; i < n_A; ++i )
  {
    float*    a1       = buff_A + (i  )*cs_A + (0  )*rs_A;
    float*    u        = buff_u;
    float*    psi1     = buff_y + (i  )*inc_y;
    float*    z        = buff_z;
    float*    nu1      = buff_v + (i  )*inc_v;
    float*    alpha    = buff_alpha;
    float     temp1;
    float     temp2;

    /*------------------------------------------------------------*/

    // bl1_smult3( alpha, psi1, &temp1 );
    temp1 = *alpha * *psi1;

    // bl1_smult3( alpha, nu1, &temp2 );
    temp2 = *alpha * *nu1;

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp1,
                u,  inc_u,
                a1, rs_A );
/*
    F77_saxpy( &m_A,
               &temp1,
               u,  &inc_u,
               a1, &rs_A );
*/

    bl1_saxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp2,
                z,  inc_z,
                a1, rs_A );
/*
    F77_saxpy( &m_A,
               &temp2,
               z,  &inc_z,
               a1, &rs_A );
*/

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Gerc2_opd_var1( integer m_A,
                                    integer n_A,
                                    double* buff_alpha, 
                                    double* buff_u, integer inc_u, 
                                    double* buff_y, integer inc_y, 
                                    double* buff_z, integer inc_z, 
                                    double* buff_v, integer inc_v, 
                                    double* buff_A, integer rs_A, integer cs_A ) 
{
  integer       i;

  for ( i = 0; i < n_A; ++i )
  {
/*
   Effective computation:
   A = A + alpha * ( u * y' + z * v' );
*/
    double*   restrict a1       = buff_A + (i  )*cs_A + (0  )*rs_A;
    double*   restrict u        = buff_u;
    double*   restrict psi1     = buff_y + (i  )*inc_y;
    double*   restrict z        = buff_z;
    double*   restrict nu1      = buff_v + (i  )*inc_v;
    double*   restrict alpha    = buff_alpha;
    double    alpha_conj_psi1;
    double    alpha_conj_nu1;

    /*------------------------------------------------------------*/

    bl1_dmult3( alpha, psi1, &alpha_conj_psi1 );

    bl1_dmult3( alpha, nu1, &alpha_conj_nu1 );

    bl1_daxpyv2b( m_A,
	              &alpha_conj_psi1,
	              &alpha_conj_nu1,
	              u,  inc_u,
	              z,  inc_z,
	              a1, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Gerc2_opc_var1( integer m_A,
                                    integer n_A,
                                    scomplex* buff_alpha, 
                                    scomplex* buff_u, integer inc_u, 
                                    scomplex* buff_y, integer inc_y, 
                                    scomplex* buff_z, integer inc_z, 
                                    scomplex* buff_v, integer inc_v, 
                                    scomplex* buff_A, integer rs_A, integer cs_A ) 
{
  integer       i;

  for ( i = 0; i < n_A; ++i )
  {
    scomplex* a1       = buff_A + (i  )*cs_A + (0  )*rs_A;
    scomplex* u        = buff_u;
    scomplex* psi1     = buff_y + (i  )*inc_y;
    scomplex* z        = buff_z;
    scomplex* nu1      = buff_v + (i  )*inc_v;
    scomplex* alpha    = buff_alpha;
    scomplex  psi1_conj;
    scomplex  nu1_conj;
    scomplex  temp1;
    scomplex  temp2;

    /*------------------------------------------------------------*/

    bl1_ccopyconj( psi1, &psi1_conj );
    bl1_cmult3( alpha, &psi1_conj, &temp1 );

    bl1_ccopyconj( nu1, &nu1_conj );
    bl1_cmult3( alpha, &nu1_conj, &temp2 );

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp1,
                u,  inc_u,
                a1, rs_A );
/*
    F77_caxpy( &m_A,
               &temp1,
               u,  &inc_u,
               a1, &rs_A );
*/

    bl1_caxpyv( BLIS1_NO_CONJUGATE,
                m_A,
                &temp2,
                z,  inc_z,
                a1, rs_A );
/*
    F77_caxpy( &m_A,
               &temp2,
               z,  &inc_z,
               a1, &rs_A );
*/

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}



FLA_Error FLA_Fused_Gerc2_opz_var1( integer m_A,
                                    integer n_A,
                                    dcomplex* buff_alpha, 
                                    dcomplex* buff_u, integer inc_u, 
                                    dcomplex* buff_y, integer inc_y, 
                                    dcomplex* buff_z, integer inc_z, 
                                    dcomplex* buff_v, integer inc_v, 
                                    dcomplex* buff_A, integer rs_A, integer cs_A ) 
{
  integer i;

  for ( i = 0; i < n_A; ++i )
  {
    dcomplex* restrict a1       = buff_A + (i  )*cs_A + (0  )*rs_A;
    dcomplex* restrict u        = buff_u;
    dcomplex* restrict psi1     = buff_y + (i  )*inc_y;
    dcomplex* restrict z        = buff_z;
    dcomplex* restrict nu1      = buff_v + (i  )*inc_v;
    dcomplex* restrict alpha    = buff_alpha;
    dcomplex  conj_psi1;
    dcomplex  conj_nu1;
    dcomplex  alpha_conj_psi1;
    dcomplex  alpha_conj_nu1;

    /*------------------------------------------------------------*/

    bl1_zcopyconj( psi1, &conj_psi1 );
    bl1_zmult3( alpha, &conj_psi1, &alpha_conj_psi1 );

    bl1_zcopyconj( nu1, &conj_nu1 );
    bl1_zmult3( alpha, &conj_nu1, &alpha_conj_nu1 );

    bl1_zaxpyv2b( m_A,
	              &alpha_conj_psi1,
	              &alpha_conj_nu1,
	              u,  inc_u,
	              z,  inc_z,
	              a1, rs_A );

    /*------------------------------------------------------------*/

  }

  return FLA_SUCCESS;
}

