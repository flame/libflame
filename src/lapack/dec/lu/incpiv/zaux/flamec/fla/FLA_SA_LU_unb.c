/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_SA_LU_unb( FLA_Obj U, FLA_Obj D, FLA_Obj p, FLA_Obj L )
{
  FLA_Datatype datatype;
  integer          m_U, cs_U;
  integer          m_D, cs_D;
  integer               cs_L;
  // integer               rs_U;
  integer               rs_D;
  // integer               rs_L;
  integer          m_U_min_j, m_U_min_j_min_1; 
  integer          j, ipiv;
  integer*         buff_p;

  if ( FLA_Obj_has_zero_dim( U ) ) return FLA_SUCCESS;
  
  datatype = FLA_Obj_datatype( U );

  m_U      = FLA_Obj_length( U );
  // rs_U     = FLA_Obj_row_stride( U );
  cs_U     = FLA_Obj_col_stride( U );

  m_D      = FLA_Obj_length( D );
  rs_D     = FLA_Obj_row_stride( D );
  cs_D     = FLA_Obj_col_stride( D );
  
  // rs_L     = FLA_Obj_row_stride( L );
  cs_L     = FLA_Obj_col_stride( L );

  FLA_Copy_external( U, L );
  FLA_Triangularize( FLA_UPPER_TRIANGULAR, FLA_NONUNIT_DIAG, L );

  buff_p     = ( integer * ) FLA_INT_PTR( p );

  switch ( datatype ){

  case FLA_FLOAT:
  {
    float* buff_U      = ( float * ) FLA_FLOAT_PTR( U );
    float* buff_D      = ( float * ) FLA_FLOAT_PTR( D );
    float* buff_L      = ( float * ) FLA_FLOAT_PTR( L );
    float* buff_minus1 = ( float * ) FLA_FLOAT_PTR( FLA_MINUS_ONE );
    float  L_tmp;
    float  D_tmp;
    float  d_inv_Ljj;

    for ( j = 0; j < m_U; ++j )
    {
      bl1_samax( m_D, 
                 buff_D + j*cs_D + 0*rs_D,
                 rs_D,
                 &ipiv );

      L_tmp = buff_L[ j*cs_L + j    ];
      D_tmp = buff_D[ j*cs_D + ipiv ];

      if ( fabsf( L_tmp ) < fabsf( D_tmp ) )
      {
        bl1_sswap( m_U,
                   buff_L + 0*cs_L + j,    cs_L,
                   buff_D + 0*cs_D + ipiv, cs_D ); 

        buff_p[ j ] = ipiv + m_U - j;
      }        
      else
      {
        buff_p[ j ] = 0;
      }

      d_inv_Ljj = 1.0F / buff_L[ j*cs_L + j ];

      bl1_sscal( m_D,
                 &d_inv_Ljj,
                 buff_D + j*cs_D + 0, rs_D ); 

      m_U_min_j_min_1 = m_U - j - 1;

      if ( m_U_min_j_min_1 > 0  )
      {
        bl1_sger( BLIS1_NO_CONJUGATE,
                  BLIS1_NO_CONJUGATE,
                  m_D,
                  m_U_min_j_min_1,
                  buff_minus1, 
                  buff_D + (j+0)*cs_D + 0, rs_D,
                  buff_L + (j+1)*cs_L + j, cs_L,
                  buff_D + (j+1)*cs_D + 0, rs_D, cs_D );
      }

      m_U_min_j = m_U - j;

      if ( m_U_min_j > 0 ) 
      {
        bl1_scopy( m_U_min_j,
                   buff_L + j*cs_L + j, cs_L,
                   buff_U + j*cs_U + j, cs_U );
      }
    }                 
    break;
  }

  case FLA_DOUBLE:
  {
    double* buff_U      = ( double * ) FLA_DOUBLE_PTR( U );
    double* buff_D      = ( double * ) FLA_DOUBLE_PTR( D );
    double* buff_L      = ( double * ) FLA_DOUBLE_PTR( L );
    double* buff_minus1 = ( double * ) FLA_DOUBLE_PTR( FLA_MINUS_ONE );
    double  L_tmp;
    double  D_tmp;
    double  d_inv_Ljj;

    for ( j = 0; j < m_U; ++j )
    {
      bl1_damax( m_D, 
                 buff_D + j*cs_D + 0*rs_D,
                 rs_D,
                 &ipiv );

      L_tmp = buff_L[ j*cs_L + j    ];
      D_tmp = buff_D[ j*cs_D + ipiv ];

      if ( fabs( L_tmp ) < fabs( D_tmp ) )
      {
        bl1_dswap( m_U,
                   buff_L + 0*cs_L + j,    cs_L,
                   buff_D + 0*cs_D + ipiv, cs_D ); 

        buff_p[ j ] = ipiv + m_U - j;
      }        
      else
      {
        buff_p[ j ] = 0;
      }

      d_inv_Ljj = 1.0 / buff_L[ j*cs_L + j ];

      bl1_dscal( m_D,
                 &d_inv_Ljj,
                 buff_D + j*cs_D + 0, rs_D ); 

      m_U_min_j_min_1 = m_U - j - 1;

      if ( m_U_min_j_min_1 > 0  )
      {
        bl1_dger( BLIS1_NO_CONJUGATE,
                  BLIS1_NO_CONJUGATE,
                  m_D,
                  m_U_min_j_min_1,
                  buff_minus1, 
                  buff_D + (j+0)*cs_D + 0, rs_D,
                  buff_L + (j+1)*cs_L + j, cs_L,
                  buff_D + (j+1)*cs_D + 0, rs_D, cs_D );
      }

      m_U_min_j = m_U - j;

      if ( m_U_min_j > 0 ) 
      {
        bl1_dcopy( m_U_min_j,
                   buff_L + j*cs_L + j, cs_L,
                   buff_U + j*cs_U + j, cs_U );
      }
    }                 
    break;
  }

  case FLA_COMPLEX:
  {
    scomplex* buff_U      = ( scomplex * ) FLA_COMPLEX_PTR( U );
    scomplex* buff_D      = ( scomplex * ) FLA_COMPLEX_PTR( D );
    scomplex* buff_L      = ( scomplex * ) FLA_COMPLEX_PTR( L );
    scomplex* buff_minus1 = ( scomplex * ) FLA_COMPLEX_PTR( FLA_MINUS_ONE );
    scomplex  L_tmp;
    scomplex  D_tmp;
    scomplex  d_inv_Ljj;
    scomplex  Ljj;
    float     temp;

    for ( j = 0; j < m_U; ++j )
    {
      bl1_camax( m_D, 
                 buff_D + j*cs_D + 0*rs_D,
                 rs_D,
                 &ipiv );

      L_tmp = buff_L[ j*cs_L + j    ];
      D_tmp = buff_D[ j*cs_D + ipiv ];

      if ( fabsf( L_tmp.real + L_tmp.imag ) < fabsf( D_tmp.real + D_tmp.imag ) )
      {
        bl1_cswap( m_U,
                   buff_L + 0*cs_L + j,    cs_L,
                   buff_D + 0*cs_D + ipiv, cs_D ); 

        buff_p[ j ] = ipiv + m_U - j;
      }        
      else
      {
        buff_p[ j ] = 0;
      }

      Ljj = buff_L[ j*cs_L + j ];

      // d_inv_Ljj = 1.0 / Ljj
      temp = 1.0F / ( Ljj.real * Ljj.real +
                      Ljj.imag * Ljj.imag );
      d_inv_Ljj.real = Ljj.real *  temp;
      d_inv_Ljj.imag = Ljj.imag * -temp;

      bl1_cscal( m_D,
                 &d_inv_Ljj,
                 buff_D + j*cs_D + 0, rs_D ); 

      m_U_min_j_min_1 = m_U - j - 1;

      if ( m_U_min_j_min_1 > 0  )
      {
        bl1_cger( BLIS1_NO_CONJUGATE,
                  BLIS1_NO_CONJUGATE,
                  m_D,
                  m_U_min_j_min_1,
                  buff_minus1, 
                  buff_D + (j+0)*cs_D + 0, rs_D,
                  buff_L + (j+1)*cs_L + j, cs_L,
                  buff_D + (j+1)*cs_D + 0, rs_D, cs_D );
      }

      m_U_min_j = m_U - j;

      if ( m_U_min_j > 0 ) 
      {
        bl1_ccopy( m_U_min_j,
                   buff_L + j*cs_L + j, cs_L,
                   buff_U + j*cs_U + j, cs_U );
      }
    }                 
    break;
  }

  case FLA_DOUBLE_COMPLEX:
  {
    dcomplex* buff_U      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( U );
    dcomplex* buff_D      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( D );
    dcomplex* buff_L      = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( L );
    dcomplex* buff_minus1 = ( dcomplex * ) FLA_DOUBLE_COMPLEX_PTR( FLA_MINUS_ONE );
    dcomplex  L_tmp;
    dcomplex  D_tmp;
    dcomplex  d_inv_Ljj;
    dcomplex  Ljj;
    double    temp;

    for ( j = 0; j < m_U; ++j )
    {
      bl1_zamax( m_D, 
                 buff_D + j*cs_D + 0*rs_D,
                 rs_D,
                 &ipiv );

      L_tmp = buff_L[ j*cs_L + j    ];
      D_tmp = buff_D[ j*cs_D + ipiv ];

      if ( fabs( L_tmp.real + L_tmp.imag ) < fabs( D_tmp.real + D_tmp.imag ) )
      {
        bl1_zswap( m_U,
                   buff_L + 0*cs_L + j,    cs_L,
                   buff_D + 0*cs_D + ipiv, cs_D ); 

        buff_p[ j ] = ipiv + m_U - j;
      }        
      else
      {
        buff_p[ j ] = 0;
      }

      Ljj = buff_L[ j*cs_L + j ];

      // d_inv_Ljj = 1.0 / Ljj
      temp = 1.0  / ( Ljj.real * Ljj.real +
                      Ljj.imag * Ljj.imag );
      d_inv_Ljj.real = Ljj.real *  temp;
      d_inv_Ljj.imag = Ljj.imag * -temp;

      bl1_zscal( m_D,
                 &d_inv_Ljj,
                 buff_D + j*cs_D + 0, rs_D ); 

      m_U_min_j_min_1 = m_U - j - 1;

      if ( m_U_min_j_min_1 > 0  )
      {
        bl1_zger( BLIS1_NO_CONJUGATE,
                  BLIS1_NO_CONJUGATE,
                  m_D,
                  m_U_min_j_min_1,
                  buff_minus1, 
                  buff_D + (j+0)*cs_D + 0, rs_D,
                  buff_L + (j+1)*cs_L + j, cs_L,
                  buff_D + (j+1)*cs_D + 0, rs_D, cs_D );
      }

      m_U_min_j = m_U - j;

      if ( m_U_min_j > 0 ) 
      {
        bl1_zcopy( m_U_min_j,
                   buff_L + j*cs_L + j, cs_L,
                   buff_U + j*cs_U + j, cs_U );
      }
    }                 
    break;
  }

  }

  return FLA_SUCCESS;
}
