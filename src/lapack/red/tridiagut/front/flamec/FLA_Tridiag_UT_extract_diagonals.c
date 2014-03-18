
#include "FLAME.h"

FLA_Error FLA_Tridiag_UT_extract_diagonals( FLA_Uplo uplo, FLA_Obj A, FLA_Obj d, FLA_Obj e )
{
  FLA_Error r_val = FLA_SUCCESS;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Tridiag_UT_extract_diagonals_check( uplo, A, d, e );

  if ( uplo == FLA_LOWER_TRIANGULAR )
    r_val = FLA_Bidiag_UT_l_extract_diagonals( A, d, e );
  else
    r_val = FLA_Bidiag_UT_u_extract_diagonals( A, d, e );

  return r_val;
}
//// 
//// FLA_Error FLA_Tridiag_UT_l_extract_diagonals( FLA_Obj A, FLA_Obj d, FLA_Obj e )
//// {
////   FLA_Datatype datatype;
////   int          m_A;
////   int          rs_A, cs_A;
////   int          inc_d;
////   int          inc_e;
////   int          i;
//// 
////   datatype = FLA_Obj_datatype( A );
//// 
////   m_A      = FLA_Obj_length( A );
//// 
////   rs_A     = FLA_Obj_row_stride( A );
////   cs_A     = FLA_Obj_col_stride( A );
//// 
////   inc_d    = FLA_Obj_vector_inc( d );
//// 
////   inc_e    = FLA_Obj_vector_inc( e );
//// 
////   switch ( datatype )
////   {
////     case FLA_FLOAT:
////     {
////       float*    buff_A = FLA_FLOAT_PTR( A );
////       float*    buff_d = FLA_FLOAT_PTR( d );
////       float*    buff_e = FLA_FLOAT_PTR( e );
//// 
////       for ( i = 0; i < m_A; ++i )
////       {
////         float*    alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
////         float*    a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
////         float*    delta1   = buff_d + (i  )*inc_d;
////         float*    epsilon1 = buff_e + (i  )*inc_e;
//// 
////         int       m_ahead  = m_A - i - 1;
//// 
////         // delta1 = alpha11;
////         *delta1 = *alpha11;
//// 
////         // epsilon1 = a21_t;
////         if ( m_ahead > 0 )
////           *epsilon1 = *a21_t;
////       }
//// 
////       break;
////     }
//// 
////     case FLA_DOUBLE:
////     {
////       double*   buff_A = FLA_DOUBLE_PTR( A );
////       double*   buff_d = FLA_DOUBLE_PTR( d );
////       double*   buff_e = FLA_DOUBLE_PTR( e );
//// 
////       for ( i = 0; i < m_A; ++i )
////       {
////         double*   alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
////         double*   a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
////         double*   delta1   = buff_d + (i  )*inc_d;
////         double*   epsilon1 = buff_e + (i  )*inc_e;
//// 
////         int       m_ahead  = m_A - i - 1;
//// 
////         // delta1 = alpha11;
////         *delta1 = *alpha11;
//// 
////         // epsilon1 = a21_t;
////         if ( m_ahead > 0 )
////           *epsilon1 = *a21_t;
////       }
//// 
////       break;
////     }
//// 
////     case FLA_COMPLEX:
////     {
////       scomplex* buff_A = FLA_COMPLEX_PTR( A );
////       float*    buff_d = FLA_FLOAT_PTR( d );
////       float*    buff_e = FLA_FLOAT_PTR( e );
//// 
////       for ( i = 0; i < m_A; ++i )
////       {
////         scomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
////         scomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
////         float*    delta1   = buff_d + (i  )*inc_d;
////         float*    epsilon1 = buff_e + (i  )*inc_e;
//// 
////         int       m_ahead  = m_A - i - 1;
//// 
////         // delta1 = alpha11;
////         *delta1 = alpha11->real;
//// 
////         // epsilon1 = a21_t;
////         if ( m_ahead > 0 )
////           *epsilon1 = a21_t->real;
////       }
//// 
////       break;
////     }
//// 
////     case FLA_DOUBLE_COMPLEX:
////     {
////       dcomplex* buff_A = FLA_DOUBLE_COMPLEX_PTR( A );
////       double*   buff_d = FLA_DOUBLE_PTR( d );
////       double*   buff_e = FLA_DOUBLE_PTR( e );
//// 
////       for ( i = 0; i < m_A; ++i )
////       {
////         dcomplex* alpha11  = buff_A + (i  )*cs_A + (i  )*rs_A;
////         dcomplex* a21_t    = buff_A + (i  )*cs_A + (i+1)*rs_A;
////         double*   delta1   = buff_d + (i  )*inc_d;
////         double*   epsilon1 = buff_e + (i  )*inc_e;
//// 
////         int       m_ahead  = m_A - i - 1;
//// 
////         // delta1 = alpha11;
////         *delta1 = alpha11->real;
//// 
////         // epsilon1 = a21_t;
////         if ( m_ahead > 0 )
////           *epsilon1 = a21_t->real;
////       }
//// 
////       break;
////     }
////   }
//// 
////   return FLA_SUCCESS;
//// }
//// 
