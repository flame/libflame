/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


static TLS_CLASS_SPEC unsigned int fla_error_checking_level = FLA_INTERNAL_ERROR_CHECKING_LEVEL;



unsigned int FLA_Check_error_level()
{
  return fla_error_checking_level;
}

unsigned int FLA_Check_error_level_set( unsigned int new_level )
{
  FLA_Error    e_val;
  unsigned int old_level;

  e_val = FLA_Check_valid_error_level( new_level );
  FLA_Check_error_code( e_val );

  old_level = fla_error_checking_level;

  fla_error_checking_level = new_level;

  return old_level;
}

FLA_Error FLA_Check_error_code_helper( int code, char* file, int line )
{
  if ( code == FLA_SUCCESS )
    return code;

  //if ( /* fatal error checking enabled */ )
  if ( TRUE )
  {
    if ( FLA_ERROR_CODE_MAX <= code && code <= FLA_ERROR_CODE_MIN )
    {
      FLA_Print_message( FLA_Error_string_for_code( code ),
                         file, line );
      FLA_Abort();
    }
    else
    {
      FLA_Print_message( FLA_Error_string_for_code( FLA_UNDEFINED_ERROR_CODE ),
                         file, line );
      FLA_Abort();
    }
  }

  return code;
}

FLA_Error FLA_Check_valid_side( FLA_Side side )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( side != FLA_LEFT && 
       side != FLA_RIGHT && 
       side != FLA_TOP &&
       side != FLA_BOTTOM )
    e_val = FLA_INVALID_SIDE;

  return e_val;
}

FLA_Error FLA_Check_valid_uplo( FLA_Uplo uplo )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( uplo != FLA_LOWER_TRIANGULAR && 
       uplo != FLA_UPPER_TRIANGULAR )
    e_val = FLA_INVALID_UPLO;

  return e_val;
}

FLA_Error FLA_Check_valid_trans( FLA_Trans trans )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( trans != FLA_NO_TRANSPOSE &&
       trans != FLA_TRANSPOSE &&
       trans != FLA_CONJ_TRANSPOSE &&
       trans != FLA_CONJ_NO_TRANSPOSE )
    e_val = FLA_INVALID_TRANS;

  return e_val;
}

FLA_Error FLA_Check_valid_diag( FLA_Diag diag )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( diag != FLA_NONUNIT_DIAG && 
       diag != FLA_UNIT_DIAG &&
       diag != FLA_ZERO_DIAG )
    e_val = FLA_INVALID_DIAG;

  return e_val;
}

FLA_Error FLA_Check_valid_conj( FLA_Conj conj )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( conj != FLA_NO_CONJUGATE &&
       conj != FLA_CONJUGATE )
    e_val = FLA_INVALID_CONJ;

  return e_val;
}

FLA_Error FLA_Check_valid_direct( FLA_Direct direct )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( direct != FLA_FORWARD &&
       direct != FLA_BACKWARD )
    e_val = FLA_INVALID_DIRECT;

  return e_val;
}

FLA_Error FLA_Check_valid_storev( FLA_Store storev )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( storev != FLA_COLUMNWISE &&
       storev != FLA_ROWWISE )
    e_val = FLA_INVALID_STOREV;

  return e_val;
}

FLA_Error FLA_Check_valid_inverse( FLA_Inv inv )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( inv != FLA_NO_INVERSE &&
       inv != FLA_INVERSE )
    e_val = FLA_INVALID_INVERSE;

  return e_val;
}

FLA_Error FLA_Check_valid_datatype( FLA_Datatype datatype )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( datatype != FLA_INT &&
       datatype != FLA_FLOAT && 
       datatype != FLA_DOUBLE && 
       datatype != FLA_COMPLEX && 
       datatype != FLA_DOUBLE_COMPLEX && 
       datatype != FLA_CONSTANT )
    e_val = FLA_INVALID_DATATYPE;

  return e_val;
}

FLA_Error FLA_Check_valid_object_datatype( FLA_Obj A )
{
  FLA_Error    e_val;
  FLA_Datatype datatype;

  datatype = FLA_Obj_datatype( A );

  e_val = FLA_Check_valid_datatype( datatype );

  return e_val;
}

FLA_Error FLA_Check_floating_datatype( FLA_Datatype datatype )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( datatype != FLA_CONSTANT &&
       datatype != FLA_FLOAT && 
       datatype != FLA_DOUBLE && 
       datatype != FLA_COMPLEX && 
       datatype != FLA_DOUBLE_COMPLEX )
    e_val = FLA_INVALID_FLOATING_DATATYPE;
  
  return e_val;
}

FLA_Error FLA_Check_int_datatype( FLA_Datatype datatype )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( datatype != FLA_CONSTANT &&
       datatype != FLA_INT )
    e_val = FLA_INVALID_INTEGER_DATATYPE;

  return e_val;
}

FLA_Error FLA_Check_real_datatype( FLA_Datatype datatype )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( datatype != FLA_CONSTANT &&
       datatype != FLA_FLOAT &&
       datatype != FLA_DOUBLE )
    e_val = FLA_INVALID_REAL_DATATYPE;

  return e_val;
}

FLA_Error FLA_Check_complex_datatype( FLA_Datatype datatype )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( datatype != FLA_CONSTANT &&
       datatype != FLA_COMPLEX &&
       datatype != FLA_DOUBLE_COMPLEX )
    e_val = FLA_INVALID_COMPLEX_DATATYPE;

  return e_val;
}

FLA_Error FLA_Check_floating_object( FLA_Obj A )
{
  FLA_Error    e_val = FLA_SUCCESS;
  FLA_Datatype datatype;

  datatype = FLA_Obj_datatype( A );

  if ( FLA_Check_floating_datatype( datatype ) != FLA_SUCCESS )
    e_val = FLA_OBJECT_NOT_FLOATING_POINT;

  return e_val;
}

FLA_Error FLA_Check_int_object( FLA_Obj A )
{
  FLA_Error    e_val = FLA_SUCCESS;
  FLA_Datatype datatype;

  datatype = FLA_Obj_datatype( A );

  if ( FLA_Check_int_datatype( datatype ) != FLA_SUCCESS )
    e_val = FLA_OBJECT_NOT_INTEGER;

  return e_val;
}

FLA_Error FLA_Check_real_object( FLA_Obj A )
{
  FLA_Error    e_val = FLA_SUCCESS;
  FLA_Datatype datatype;

  datatype = FLA_Obj_datatype( A );

  if ( FLA_Check_real_datatype( datatype ) != FLA_SUCCESS )
    e_val = FLA_OBJECT_NOT_REAL;

  return e_val;
}

FLA_Error FLA_Check_comparable_object( FLA_Obj A )
{
  FLA_Error    e_val = FLA_SUCCESS;
  FLA_Datatype datatype;

  datatype = FLA_Obj_datatype( A );

  if ( FLA_Check_int_datatype( datatype ) != FLA_SUCCESS &&
       FLA_Check_real_datatype( datatype ) != FLA_SUCCESS )
    e_val = FLA_OBJECT_NOT_COMPARABLE;
  
  return e_val;
}

FLA_Error FLA_Check_complex_object( FLA_Obj A )
{
  FLA_Error    e_val = FLA_SUCCESS;
  FLA_Datatype datatype;

  datatype = FLA_Obj_datatype( A );

  if ( FLA_Check_complex_datatype( datatype ) != FLA_SUCCESS )
    e_val = FLA_OBJECT_NOT_COMPLEX;

  return e_val;
}

FLA_Error FLA_Check_identical_object_precision( FLA_Obj A, FLA_Obj B )
{
  FLA_Error    e_val = FLA_SUCCESS;
  FLA_Datatype datatype_A;
  FLA_Datatype datatype_B;
  dim_t        precision_A;
  dim_t        precision_B;

  datatype_A = FLA_Obj_datatype( A );
  datatype_B = FLA_Obj_datatype( B );

  if ( datatype_A == FLA_CONSTANT ||
       datatype_B == FLA_CONSTANT )
  {
    return FLA_SUCCESS;
  }

  if ( FLA_Check_floating_object( A ) != FLA_SUCCESS ||
       FLA_Check_floating_object( B ) != FLA_SUCCESS )
  {
    return FLA_OBJECT_NOT_FLOATING_POINT;
  }

  datatype_A = FLA_Obj_datatype( A );
  datatype_B = FLA_Obj_datatype( B );

  precision_A = FLA_Obj_datatype_size( datatype_A );
  precision_B = FLA_Obj_datatype_size( datatype_B );

  if ( FLA_Obj_is_complex( A ) )
    precision_A = precision_A / 2;

  if ( FLA_Obj_is_complex( B ) )
    precision_B = precision_B / 2;

  if ( precision_A != precision_B )
    e_val = FLA_INCONSISTENT_OBJECT_PRECISION;

  return e_val;
}

FLA_Error FLA_Check_consistent_object_datatype( FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_datatype( A ) != FLA_CONSTANT &&
       FLA_Obj_datatype( B ) != FLA_CONSTANT )
    if ( FLA_Obj_datatype( A ) != FLA_Obj_datatype( B ) )
      e_val = FLA_INCONSISTENT_DATATYPES;

  return e_val;
}

FLA_Error FLA_Check_consistent_datatype( FLA_Datatype datatype, FLA_Obj A )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_datatype( A ) != FLA_CONSTANT &&
                    datatype != FLA_CONSTANT )
    if ( FLA_Obj_datatype( A ) != datatype )
      e_val = FLA_INCONSISTENT_DATATYPES;

  return e_val;
}

FLA_Error FLA_Check_square( FLA_Obj A )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_length( A ) != FLA_Obj_width( A ) )
    e_val = FLA_OBJECT_NOT_SQUARE;

  return e_val;
}

FLA_Error FLA_Check_if_scalar( FLA_Obj A )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_length( A ) != 1 || FLA_Obj_width( A ) != 1 )
    e_val = FLA_OBJECT_NOT_SCALAR;

  return e_val;
}

FLA_Error FLA_Check_if_vector( FLA_Obj A )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_length( A ) != 1 && FLA_Obj_width( A ) != 1 )
    e_val = FLA_OBJECT_NOT_VECTOR;

  return e_val;
}

FLA_Error FLA_Check_conformal_dims( FLA_Trans trans, FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( trans == FLA_NO_TRANSPOSE || trans == FLA_CONJ_NO_TRANSPOSE )
  {
    if ( FLA_Obj_length( A ) != FLA_Obj_length( B ) )
      e_val = FLA_NONCONFORMAL_DIMENSIONS;

    if ( FLA_Obj_width( A ) != FLA_Obj_width( B ) )
      e_val = FLA_NONCONFORMAL_DIMENSIONS;
  }
  else
  {
    if ( FLA_Obj_width( A ) != FLA_Obj_length( B ) )
      e_val = FLA_NONCONFORMAL_DIMENSIONS;

    if ( FLA_Obj_length( A ) != FLA_Obj_width( B ) )
      e_val = FLA_NONCONFORMAL_DIMENSIONS;
  }

  return e_val;
}

FLA_Error FLA_Check_matrix_matrix_dims( FLA_Trans transa, FLA_Trans transb, FLA_Obj A, FLA_Obj B, FLA_Obj C )
{
  FLA_Error e_val = FLA_SUCCESS;
  dim_t     k_A, k_B;
  dim_t     m_A, m_C;
  dim_t     n_B, n_C;

  m_A = ( transa == FLA_NO_TRANSPOSE || 
          transa == FLA_CONJ_NO_TRANSPOSE ? FLA_Obj_length( A ) :
                                            FLA_Obj_width( A )  );
  k_A = ( transa == FLA_NO_TRANSPOSE ||
          transa == FLA_CONJ_NO_TRANSPOSE ? FLA_Obj_width( A )  :
                                            FLA_Obj_length( A ) );

  k_B = ( transb == FLA_NO_TRANSPOSE ||
          transb == FLA_CONJ_NO_TRANSPOSE ? FLA_Obj_length( B ) :
                                            FLA_Obj_width( B )  );
  n_B = ( transb == FLA_NO_TRANSPOSE ||
          transb == FLA_CONJ_NO_TRANSPOSE ? FLA_Obj_width( B )  :
                                            FLA_Obj_length( B ) );

  m_C = FLA_Obj_length( C );
  n_C = FLA_Obj_width( C );

  if ( m_A != m_C )
    e_val = FLA_NONCONFORMAL_DIMENSIONS;

  if ( k_A != k_B )
    e_val = FLA_NONCONFORMAL_DIMENSIONS;

  if ( n_B != n_C )
    e_val = FLA_NONCONFORMAL_DIMENSIONS;

  return e_val;
}

FLA_Error FLA_Check_matrix_vector_dims( FLA_Trans trans, FLA_Obj A, FLA_Obj x, FLA_Obj y )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( trans == FLA_NO_TRANSPOSE || trans == FLA_CONJ_NO_TRANSPOSE )
  {
    if ( FLA_Obj_width( A ) != FLA_Obj_vector_dim( x ) )
      e_val = FLA_NONCONFORMAL_DIMENSIONS;

    if ( FLA_Obj_length( A ) != FLA_Obj_vector_dim( y ) )
      e_val = FLA_NONCONFORMAL_DIMENSIONS;
  }
  else
  {
    if ( FLA_Obj_length( A ) != FLA_Obj_vector_dim( x ) )
      e_val = FLA_NONCONFORMAL_DIMENSIONS;

    if ( FLA_Obj_width( A ) != FLA_Obj_vector_dim( y ) )
      e_val = FLA_NONCONFORMAL_DIMENSIONS;
  }

  return e_val;
}

FLA_Error FLA_Check_equal_vector_dims( FLA_Obj x, FLA_Obj y )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_vector_dim( x ) != FLA_Obj_vector_dim( y ) )
    e_val = FLA_UNEQUAL_VECTOR_DIMS;

  return e_val;
}

FLA_Error FLA_Check_conj1_trans_and_datatype( FLA_Trans trans, FLA_Obj A )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( trans == FLA_CONJ_TRANSPOSE || trans == FLA_CONJ_NO_TRANSPOSE ){
    if ( FLA_Obj_is_complex( A ) == FALSE )
      e_val = FLA_INVALID_TRANS_GIVEN_DATATYPE;
  }

  return e_val;
}

FLA_Error FLA_Check_hess_indices( FLA_Obj A, int ilo, int ihi )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_width( A ) == 0 && ilo != 0 && ihi != -1 )
    e_val = FLA_INVALID_HESSENBERG_INDICES;

  if ( ilo < 0 || FLA_Obj_width( A ) - 1 < ilo )
    e_val = FLA_INVALID_HESSENBERG_INDICES;

  if ( ihi < 0 || FLA_Obj_width( A ) - 1 < ihi )
    e_val = FLA_INVALID_HESSENBERG_INDICES;

  if ( ihi < ilo )
    e_val = FLA_INVALID_HESSENBERG_INDICES;

  return e_val;
}

FLA_Error FLA_Check_null_pointer( void* ptr )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( ptr == NULL )
    e_val = FLA_NULL_POINTER;

  return e_val;
}

FLA_Error FLA_Check_object_dims( FLA_Trans trans, dim_t m, dim_t n, FLA_Obj A )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( trans == FLA_NO_TRANSPOSE || trans == FLA_CONJ_NO_TRANSPOSE )
  {
    if ( FLA_Obj_length( A ) != m )
      e_val = FLA_SPECIFIED_OBJ_DIM_MISMATCH;

    if ( FLA_Obj_width( A ) != n )
      e_val = FLA_SPECIFIED_OBJ_DIM_MISMATCH;
  }
  else
  {
    if ( FLA_Obj_length( A ) != n )
      e_val = FLA_SPECIFIED_OBJ_DIM_MISMATCH;

    if ( FLA_Obj_width( A ) != m )
      e_val = FLA_SPECIFIED_OBJ_DIM_MISMATCH;
  }

  return e_val;
}

FLA_Error FLA_Check_valid_pivot_type( FLA_Pivot_type ptype )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( ptype != FLA_NATIVE_PIVOTS && ptype != FLA_LAPACK_PIVOTS )
    e_val = FLA_INVALID_PIVOT_TYPE;

  return e_val;
}

FLA_Error FLA_Check_malloc_pointer( void* ptr )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( ptr == NULL )
    e_val = FLA_MALLOC_RETURNED_NULL_POINTER;

  return e_val;
}

FLA_Error FLA_Check_base_buffer_mismatch( FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( A.base->buffer != B.base->buffer )
    e_val = FLA_OBJECT_BASE_BUFFER_MISMATCH;

  return e_val;
}

FLA_Error FLA_Check_adjacent_objects_2x2( FLA_Obj ATL, FLA_Obj ATR,
                                          FLA_Obj ABL, FLA_Obj ABR )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_length( ATL ) != FLA_Obj_length( ATR ) || 
       FLA_Obj_length( ABL ) != FLA_Obj_length( ABR ) ||
       FLA_Obj_width( ATL )  != FLA_Obj_width( ABL )  ||
       FLA_Obj_width( ATR )  != FLA_Obj_width( ABR ) )
    e_val = FLA_ADJACENT_OBJECT_DIM_MISMATCH;

  if ( ATL.offm != ABL.offm + FLA_Obj_length( ABL ) ||
       ATR.offm != ABR.offm + FLA_Obj_length( ABL ) )
    e_val = FLA_OBJECTS_NOT_VERTICALLY_ADJ;

  if ( ATL.offn != ABL.offn ||
       ATR.offn != ABR.offn )
    e_val = FLA_OBJECTS_NOT_VERTICALLY_ALIGNED;

  if ( ATL.offn != ATR.offn + FLA_Obj_width( ATR ) ||
       ABL.offn != ABR.offn + FLA_Obj_width( ATR ) )
    e_val = FLA_OBJECTS_NOT_HORIZONTALLY_ADJ;

  if ( ATL.offm != ATR.offm ||
       ABL.offm != ABR.offm )
    e_val = FLA_OBJECTS_NOT_HORIZONTALLY_ALIGNED;

  return e_val;
}

FLA_Error FLA_Check_adjacent_objects_2x1( FLA_Obj AT,
                                          FLA_Obj AB )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_width( AT ) != FLA_Obj_width( AB ) )
    e_val = FLA_ADJACENT_OBJECT_DIM_MISMATCH;

  if ( AB.offm != AT.offm + FLA_Obj_length( AT ) )
    e_val = FLA_OBJECTS_NOT_VERTICALLY_ADJ;

  if ( AB.offn != AT.offn )
    e_val = FLA_OBJECTS_NOT_VERTICALLY_ALIGNED;

  return e_val;
}

FLA_Error FLA_Check_adjacent_objects_1x2( FLA_Obj AL, FLA_Obj AR )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_length( AL ) != FLA_Obj_length( AR ) )
    e_val = FLA_ADJACENT_OBJECT_DIM_MISMATCH;

  if ( AR.offn != AL.offn + FLA_Obj_width( AL ) )
    e_val = FLA_OBJECTS_NOT_HORIZONTALLY_ADJ;

  if ( AL.offm != AR.offm )
    e_val = FLA_OBJECTS_NOT_HORIZONTALLY_ALIGNED;

  return e_val;
}

FLA_Error FLA_Check_blocksize_value( dim_t b )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( b <= 0 )
    e_val = FLA_INVALID_BLOCKSIZE_VALUE;

  return e_val;
}

FLA_Error FLA_Check_blocksize_object( FLA_Datatype datatype, fla_blocksize_t* bp )
{
  FLA_Error e_val = FLA_SUCCESS;
  dim_t     b;

  b = FLA_Blocksize_extract( datatype, bp );
  if ( b <= 0 )
    e_val = FLA_INVALID_BLOCKSIZE_OBJ;

  return e_val;
}

FLA_Error FLA_Check_file_descriptor( int fd )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( fd == -1 )
    e_val = FLA_OPEN_RETURNED_ERROR;

  return e_val;
}

FLA_Error FLA_Check_lseek_result( int requested_offset, int lseek_r_val )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( lseek_r_val != requested_offset )
    e_val = FLA_LSEEK_RETURNED_ERROR;

  return e_val;
}

FLA_Error FLA_Check_close_result( int close_r_val )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( close_r_val == -1 )
    e_val = FLA_CLOSE_RETURNED_ERROR;

  return e_val;
}

FLA_Error FLA_Check_unlink_result( int unlink_r_val )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( unlink_r_val == -1 )
    e_val = FLA_UNLINK_RETURNED_ERROR;

  return e_val;
}

FLA_Error FLA_Check_read_result( int requested_size, int read_r_val )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( read_r_val == -1 )
    e_val = FLA_READ_RETURNED_ERROR;

  return e_val;
}

FLA_Error FLA_Check_write_result( int requested_size, int write_r_val )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( write_r_val != requested_size )
    e_val = FLA_WRITE_RETURNED_ERROR;

  return e_val;
}

FLA_Error FLA_Check_valid_quadrant( FLA_Quadrant quad )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( quad != FLA_TL && 
       quad != FLA_TR && 
       quad != FLA_BL &&
       quad != FLA_BR )
    e_val = FLA_INVALID_QUADRANT;

  return e_val;
}

FLA_Error FLA_Check_vector_dim_min( FLA_Obj x, dim_t min_dim )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_vector_dim( x ) < min_dim )
    e_val = FLA_VECTOR_DIM_BELOW_MIN;

  return e_val;
}

FLA_Error FLA_Check_pthread_create_result( int pthread_create_r_val )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( pthread_create_r_val != 0 )
    e_val = FLA_PTHREAD_CREATE_RETURNED_ERROR;

  return e_val;
}

FLA_Error FLA_Check_pthread_join_result( int pthread_join_r_val )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( pthread_join_r_val != 0 )
    e_val = FLA_PTHREAD_JOIN_RETURNED_ERROR;

  return e_val;
}

FLA_Error FLA_Check_valid_isgn_value( FLA_Obj isgn )
{
  FLA_Error e_val = FLA_SUCCESS;
  
  if ( !FLA_Obj_is( isgn, FLA_ONE ) && 
       !FLA_Obj_is( isgn, FLA_MINUS_ONE ) )
    e_val = FLA_INVALID_ISGN_VALUE;

  return e_val;
}

FLA_Error FLA_Check_sylv_matrix_dims( FLA_Obj A, FLA_Obj B, FLA_Obj C )
{
  FLA_Error e_val = FLA_SUCCESS;
  dim_t     m_A, m_C;
  dim_t     n_B, n_C;

  m_A = FLA_Obj_length( A );

  n_B = FLA_Obj_width( B );

  m_C = FLA_Obj_length( C );
  n_C = FLA_Obj_width( C );

  if ( m_A != m_C )
    e_val = FLA_NONCONFORMAL_DIMENSIONS;

  if ( n_B != n_C )
    e_val = FLA_NONCONFORMAL_DIMENSIONS;

  return e_val;
}

FLA_Error FLA_Check_chol_failure( FLA_Error r_val )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( r_val > 0 )
    e_val = FLA_CHOL_FAILED_MATRIX_NOT_SPD;

  return e_val;
}

FLA_Error FLA_Check_valid_elemtype( FLA_Elemtype elemtype )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( elemtype != FLA_SCALAR &&
       elemtype != FLA_MATRIX )
    e_val = FLA_INVALID_ELEMTYPE;

  return e_val;
}

FLA_Error FLA_Check_posix_memalign_failure( int r_val )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( r_val != 0 )
    e_val = FLA_POSIX_MEMALIGN_FAILED;

  return e_val;
}

FLA_Error FLA_Check_submatrix_dims_and_offset( dim_t m, dim_t n, dim_t i, dim_t j, FLA_Obj A )
{
  FLA_Error e_val = FLA_SUCCESS;
  dim_t     m_A, n_A;

  if ( FLA_Obj_elemtype( A ) == FLA_MATRIX )
  {
    m_A = FLASH_Obj_scalar_length( A );
    n_A = FLASH_Obj_scalar_width( A );
  }
  else
  {
    m_A = FLA_Obj_length( A );
    n_A = FLA_Obj_width( A );
  }

  if      ( i     > m_A || j     > n_A )
    e_val = FLA_INVALID_SUBMATRIX_OFFSET;
  else if ( i + m > m_A || j + n > n_A )
    e_val = FLA_INVALID_SUBMATRIX_DIMS;

  return e_val;
}

FLA_Error FLA_Check_object_scalar_elemtype( FLA_Obj A )
{
  FLA_Error     e_val = FLA_SUCCESS;
  FLA_Elemtype  elemtype;

  elemtype = FLA_Obj_elemtype( A );

  if ( elemtype != FLA_SCALAR )
    e_val = FLA_OBJECT_NOT_SCALAR_ELEMTYPE;

  return e_val;
}

FLA_Error FLA_Check_object_matrix_elemtype( FLA_Obj A )
{
  FLA_Error     e_val = FLA_SUCCESS;
  FLA_Elemtype  elemtype;

  elemtype = FLA_Obj_elemtype( A );

  if ( elemtype != FLA_MATRIX )
    e_val = FLA_OBJECT_NOT_MATRIX_ELEMTYPE;

  return e_val;
}

FLA_Error FLA_Check_num_threads( unsigned int n_threads )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( n_threads < 1 )
    e_val = FLA_ENCOUNTERED_NON_POSITIVE_NTHREADS;

  return e_val;
}

FLA_Error FLA_Check_conj_and_datatype( FLA_Conj conj, FLA_Obj A )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( conj == FLA_CONJUGATE ){
    if ( FLA_Obj_is_complex( A ) == FALSE )
      e_val = FLA_INVALID_CONJ_GIVEN_DATATYPE;
  }

  return e_val;
}

FLA_Error FLA_Check_valid_complex_trans( FLA_Trans trans )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( trans != FLA_NO_TRANSPOSE &&
       trans != FLA_CONJ_TRANSPOSE )
    e_val = FLA_INVALID_COMPLEX_TRANS;

  return e_val;
}

FLA_Error FLA_Check_valid_real_trans( FLA_Trans trans )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( trans != FLA_NO_TRANSPOSE &&
       trans != FLA_TRANSPOSE )
    e_val = FLA_INVALID_REAL_TRANS;

  return e_val;
}

FLA_Error FLA_Check_valid_blas_trans( FLA_Trans trans )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( trans != FLA_NO_TRANSPOSE &&
       trans != FLA_TRANSPOSE &&
       trans != FLA_CONJ_TRANSPOSE )
    e_val = FLA_INVALID_BLAS_TRANS;

  return e_val;
}

FLA_Error FLA_Check_nonconstant_datatype( FLA_Datatype datatype )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( datatype != FLA_INT &&
       datatype != FLA_FLOAT && 
       datatype != FLA_DOUBLE && 
       datatype != FLA_COMPLEX && 
       datatype != FLA_DOUBLE_COMPLEX )
    e_val = FLA_INVALID_NONCONSTANT_DATATYPE;
  
  return e_val;
}

FLA_Error FLA_Check_nonconstant_object( FLA_Obj A )
{
  FLA_Error    e_val = FLA_SUCCESS;
  FLA_Datatype datatype;

  datatype = FLA_Obj_datatype( A );

  if ( FLA_Check_nonconstant_datatype( datatype ) != FLA_SUCCESS )
    e_val = FLA_OBJECT_NOT_NONCONSTANT;

  return e_val;
}

FLA_Error FLA_Check_identical_object_datatype( FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_datatype( A ) != FLA_Obj_datatype( B ) )
    e_val = FLA_OBJECT_DATATYPES_NOT_EQUAL;

  return e_val;
}

FLA_Error FLA_Check_divide_by_zero( FLA_Obj alpha )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_equals( alpha, FLA_ZERO ) )
    e_val = FLA_DIVIDE_BY_ZERO;

  return e_val;
}

FLA_Error FLA_Check_identical_object_elemtype( FLA_Obj A, FLA_Obj B )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_elemtype( A ) != FLA_Obj_elemtype( B ) )
    e_val = FLA_OBJECT_ELEMTYPES_NOT_EQUAL;

  return e_val;
}

FLA_Error FLA_Check_pivot_index_range( FLA_Obj p, dim_t k1, dim_t k2 )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_has_zero_dim( p ) )
    return e_val;

  // FGVZ: Note that we don't need to test if k1 < 0 since the type is an
  // unsigned integer. Same goes for k2 < 0. If we ever wanted to allow
  // the caller to specify the index range k1:k2 == -1:0, which would
  // result in no action. If this "no-op" behavior is needed, then k1 and
  // k2 should be type re-declared as ints.
  //if ( k1 < 0 || FLA_Obj_length( p ) - 1 < k1 )
  if ( FLA_Obj_length( p ) - 1 < k1 )
    e_val = FLA_INVALID_PIVOT_INDEX_RANGE;

  //if ( k2 < 0 || FLA_Obj_length( p ) - 1 < k2 )
  if ( FLA_Obj_length( p ) - 1 < k2 )
    e_val = FLA_INVALID_PIVOT_INDEX_RANGE;

  if ( k2 < k1 )
    e_val = FLA_INVALID_PIVOT_INDEX_RANGE;

  return e_val;
}

FLA_Error FLA_Check_householder_panel_dims( FLA_Obj A, FLA_Obj T )
{
  FLA_Error e_val = FLA_SUCCESS;
  dim_t     nb_alg;

  nb_alg = FLA_Query_blocksize( FLA_Obj_datatype( A ), FLA_DIMENSION_MIN );

  if ( FLA_Obj_length( T ) < nb_alg )
    e_val = FLA_HOUSEH_PANEL_MATRIX_TOO_SMALL;

  if ( FLA_Obj_width( T ) < FLA_Obj_min_dim( A ) )
    e_val = FLA_HOUSEH_PANEL_MATRIX_TOO_SMALL;

  return e_val;
}

FLA_Error FLA_Check_object_length_equals( FLA_Obj A, dim_t m )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_length( A ) != m )
    e_val = FLA_INVALID_OBJECT_LENGTH;

  return e_val;
}

FLA_Error FLA_Check_object_width_equals( FLA_Obj A, dim_t n )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_width( A ) != n )
    e_val = FLA_INVALID_OBJECT_WIDTH;

  return e_val;
}

FLA_Error FLA_Check_object_length_min( FLA_Obj A, dim_t m )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_length( A ) < m )
    e_val = FLA_INVALID_OBJECT_LENGTH;

  return e_val;
}

FLA_Error FLA_Check_object_width_min( FLA_Obj A, dim_t n )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_width( A ) < n )
    e_val = FLA_INVALID_OBJECT_WIDTH;

  return e_val;
}

FLA_Error FLA_Check_valid_error_level( unsigned int level )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( level != FLA_NO_ERROR_CHECKING &&
       level != FLA_MIN_ERROR_CHECKING &&
       level != FLA_FULL_ERROR_CHECKING )
    e_val = FLA_INVALID_ERROR_CHECKING_LEVEL;

  return e_val;
}

FLA_Error FLA_Check_attempted_repart_2x2( FLA_Obj A_quad, dim_t b_m, dim_t b_n )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( b_m > FLA_Obj_length( A_quad ) )
    e_val = FLA_ATTEMPTED_OVER_REPART_2X2;

  if ( b_n > FLA_Obj_width( A_quad ) )
    e_val = FLA_ATTEMPTED_OVER_REPART_2X2;

  return e_val;
}

FLA_Error FLA_Check_attempted_repart_2x1( FLA_Obj A_side, dim_t b_m )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( b_m > FLA_Obj_length( A_side ) )
    e_val = FLA_ATTEMPTED_OVER_REPART_2X1;

  return e_val;
}

FLA_Error FLA_Check_attempted_repart_1x2( FLA_Obj A_side, dim_t b_n )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( b_n > FLA_Obj_width( A_side ) )
    e_val = FLA_ATTEMPTED_OVER_REPART_1X2;

  return e_val;
}

FLA_Error FLA_Check_valid_leftright_side( FLA_Side side )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( side != FLA_LEFT && 
       side != FLA_RIGHT )
    e_val = FLA_INVALID_SIDE;

  return e_val;
}

FLA_Error FLA_Check_valid_topbottom_side( FLA_Side side )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( side != FLA_TOP && 
       side != FLA_BOTTOM )
    e_val = FLA_INVALID_SIDE;

  return e_val;
}

FLA_Error FLA_Check_matrix_strides( dim_t m, dim_t n, dim_t rs, dim_t cs )
{
  FLA_Error e_val = FLA_SUCCESS;

  // Note: The default case (whereby we interpret rs == cs == 0 as a request
  // for column-major order) is handled prior to calling this function, so we
  // never see zero strides here.

  // Disallow either of the strides to be zero.
  if ( ( rs == 0 || cs == 0 ) )
    return FLA_INVALID_STRIDE_COMBINATION;

  // Check stride consistency in cases of general stride.
  if ( rs != 1 && cs != 1 )
  {
    // We apply different tests depending on which way the strides "tilt".
    if      ( rs == cs )
    {
      // If rs == cs, then we must be dealing with an m-by-1 or a 1-by-n matrix
      // and thus at least one of the dimensions, m or n, must be unit (even if
      // the other is zero).
      if ( m != 1 && n != 1 )
        return FLA_INVALID_STRIDE_COMBINATION;
    }
    else if ( rs < cs )
    {
      // For column-major tilt, cs must be equal or larger than m * rs.
      if ( m * rs > cs )
        return FLA_INVALID_STRIDE_COMBINATION;
    }
    else if ( cs < rs )
    {
      // For row-major tilt, rs must be equal or larger than n * cs.
      if ( n * cs > rs )
        return FLA_INVALID_STRIDE_COMBINATION;
    }
  }

  if ( rs == 1 && cs == 1 )
  {
    // Only allow rs == cs == 1 for scalars and "empty" objects.
    if ( !( m == 1 && n == 1 ) &&
         !( m == 0 )           && 
         !( n == 0 )           )
      return FLA_INVALID_STRIDE_COMBINATION;
  }
  else // perform additional stride/dimension checks on non-scalars.
  {
    if ( rs == 1 )
    {
      // For column-major storage, don't allow the column stride to be less than
      // the m dimension.
      if ( cs < m )
        e_val = FLA_INVALID_COL_STRIDE;
    }
    else if ( cs == 1 )
    {
      // For row-major storage, don't allow the row stride to be less than
      // the n dimension.
      if ( rs < n )
        e_val = FLA_INVALID_ROW_STRIDE;
    }
  }

  return e_val;
}

FLA_Error FLA_Check_vector_dim( FLA_Obj x, dim_t expected_length )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_vector_dim( x ) != expected_length )
    e_val = FLA_INVALID_VECTOR_DIM;

  return e_val;
}

FLA_Error FLA_Check_row_vector( FLA_Obj x )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_length( x ) != 1 )
    e_val = FLA_EXPECTED_ROW_VECTOR;

  return e_val;
}

FLA_Error FLA_Check_col_vector( FLA_Obj x )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_width( x ) != 1 )
    e_val = FLA_EXPECTED_COL_VECTOR;

  return e_val;
}

FLA_Error FLA_Check_valid_evd_type( FLA_Evd_type evd_type )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( evd_type != FLA_EVD_WITHOUT_VECTORS && 
       evd_type != FLA_EVD_WITH_VECTORS )
    e_val = FLA_INVALID_EVD_TYPE;

  return e_val;
}

FLA_Error FLA_Check_valid_svd_type( FLA_Svd_type svd_type )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( svd_type != FLA_SVD_VECTORS_ALL && 
       svd_type != FLA_SVD_VECTORS_MIN_COPY &&
       svd_type != FLA_SVD_VECTORS_MIN_OVERWRITE &&
       svd_type != FLA_SVD_VECTORS_NONE )
    e_val = FLA_INVALID_SVD_TYPE;

  return e_val;
}

FLA_Error FLA_Check_valid_svd_type_combination( FLA_Svd_type svd_type_u, FLA_Svd_type svd_type_v )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( svd_type_u == FLA_SVD_VECTORS_MIN_OVERWRITE &&
       svd_type_v == FLA_SVD_VECTORS_MIN_OVERWRITE )
    e_val = FLA_INVALID_SVD_TYPE_COMBINATION;
  
  return e_val;
}

FLA_Error FLA_Check_valid_svd_type_and_trans_combination( FLA_Svd_type svd_type_u, FLA_Trans transu,
                                                          FLA_Svd_type svd_type_v, FLA_Trans transv )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( svd_type_u == FLA_SVD_VECTORS_MIN_OVERWRITE )
    if ( transu == FLA_TRANSPOSE || 
         transu == FLA_CONJ_TRANSPOSE )
      e_val = FLA_INVALID_SVD_TYPE_AND_TRANS_COMBINATION;
  if ( svd_type_v == FLA_SVD_VECTORS_MIN_OVERWRITE )
    if ( transv == FLA_NO_TRANSPOSE || 
         transv == FLA_CONJ_NO_TRANSPOSE )
      e_val = FLA_INVALID_SVD_TYPE_AND_TRANS_COMBINATION;
  
  return e_val;
}

FLA_Error FLA_Check_valid_machval( FLA_Machval val )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( val != FLA_MACH_EPS      &&
       val != FLA_MACH_SFMIN    &&
       val != FLA_MACH_BASE     &&
       val != FLA_MACH_PREC     &&
       val != FLA_MACH_NDIGMANT &&
       val != FLA_MACH_RND      &&
       val != FLA_MACH_EMIN     &&
       val != FLA_MACH_RMIN     &&
       val != FLA_MACH_EMAX     &&
       val != FLA_MACH_RMAX     &&
       val != FLA_MACH_EPS2 )
    e_val = FLA_INVALID_MACHVAL;

  return e_val;
}

FLA_Error FLA_Check_valid_diag_offset( FLA_Obj A, FLA_Diag_off offset )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_min_dim( A ) <= f2c_abs( offset ) )
    e_val = FLA_INVALID_DIAG_OFFSET;

  return e_val;
}

FLA_Error FLA_Check_col_storage( FLA_Obj A )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_row_stride( A ) != 1 )
    e_val = FLA_EXPECTED_COL_STORAGE;

  return e_val;
}

FLA_Error FLA_Check_row_storage( FLA_Obj A )
{
  FLA_Error e_val = FLA_SUCCESS;

  if ( FLA_Obj_col_stride( A ) != 1 )
    e_val = FLA_EXPECTED_ROW_STORAGE;

  return e_val;
}

