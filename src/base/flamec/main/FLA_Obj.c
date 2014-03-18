
#include "FLAME.h"


#ifdef FLA_ENABLE_SCC
typedef volatile unsigned char* t_vcharp;
t_vcharp RCCE_shmalloc(size_t);
void     RCCE_shfree(t_vcharp);
int      RCCE_ue(void);


void* FLA_shmalloc( size_t size )
{
  return ( void * ) RCCE_shmalloc( size );
}


void FLA_shfree( void* ptr )
{
  RCCE_shfree( ( t_vcharp ) ptr );
}


FLA_Bool FLA_is_owner( void )
{
  if ( RCCE_ue() == 0 )
    return TRUE;
  return FALSE;
}
#endif

FLA_Error FLA_Obj_nullify( FLA_Obj *obj )
{
  // Nullify the fields in the view object.
  obj->m       = 0;
  obj->n       = 0;
  obj->offm    = 0;
  obj->offn    = 0;
  obj->m_inner = 0;
  obj->n_inner = 0;
  obj->base    = NULL;

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_create( FLA_Datatype datatype, dim_t m, dim_t n, dim_t rs, dim_t cs, FLA_Obj *obj )
{
  FLA_Obj_create_ext( datatype, FLA_SCALAR, m, n, m, n, rs, cs, obj );

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_create_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t m, dim_t n, dim_t m_inner, dim_t n_inner, dim_t rs, dim_t cs, FLA_Obj *obj )
{
  size_t buffer_size;
  size_t n_elem;

  // Adjust the strides, if necessary.
  FLA_adjust_strides( m, n, &rs, &cs );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_ext_check( datatype, elemtype, m, n, m_inner, n_inner, rs, cs, obj );

  // Populate the fields in the view object.
  obj->m                = m;
  obj->n                = n;
  obj->offm             = 0;
  obj->offn             = 0;
  obj->m_inner          = m_inner;
  obj->n_inner          = n_inner;

  // Allocate the base object field.
  obj->base             = ( FLA_Base_obj * ) FLA_malloc( sizeof( FLA_Base_obj ) );

  // Populate the fields in the base object.
  obj->base->datatype   = datatype;
  obj->base->elemtype   = elemtype;
  obj->base->m          = m;
  obj->base->n          = n;
  obj->base->m_inner    = m_inner;
  obj->base->n_inner    = n_inner;
  obj->base->id         = ( unsigned long ) obj->base;
  obj->base->m_index    = 0;
  obj->base->n_index    = 0;

  // Compute the number of elements needed for the buffer, adjusting
  // the strides for alignment if needed.
  n_elem = FLA_compute_num_elem( FLA_Obj_elem_size( *obj ),
                                 m, n, &rs, &cs );

  // Compute the buffer size in bytes.
  buffer_size = ( size_t ) n_elem *
                ( size_t ) FLA_Obj_elem_size( *obj );

  // Allocate the base object's element buffer.
#ifdef FLA_ENABLE_SCC
  obj->base->buffer = ( FLA_Obj_elemtype( *obj ) == FLA_MATRIX ? FLA_malloc( buffer_size ) : FLA_shmalloc( buffer_size ) );
#else
  obj->base->buffer = FLA_malloc( buffer_size );
#endif
  obj->base->buffer_info = 0;

  // Just in case this is a FLASH object, save the number of elements
  // allocated so that we can more easily free the elements later on.
  obj->base->n_elem_alloc = n_elem;

  // Save the row and column strides used in the memory allocation.
  obj->base->rs     = rs;
  obj->base->cs     = cs;

#ifdef FLA_ENABLE_SUPERMATRIX
  // Initialize SuperMatrix fields.
  obj->base->n_read_tasks   = 0;
  obj->base->read_task_head = NULL;
  obj->base->read_task_tail = NULL;
  obj->base->write_task     = NULL;
#endif

  return FLA_SUCCESS;
}


dim_t FLA_compute_num_elem( dim_t elem_size, dim_t m, dim_t n, dim_t* rs, dim_t* cs )
{
  dim_t n_elem;

  // Determine the amount of space we need to allocate based on the values of
  // the row and column strides.
  if ( m == 0 || n == 0 )
  {
    // For empty objects, set the length of the buffer to 0. Row and column
    // strides should remain unchanged (because alignment is not needed).
    n_elem = 0;
  }
  else if ( *rs == 1 )
  {
    // For column-major storage, use cs for computing the length of the buffer
    // to allocate.

    // Align the leading dimension to some user-defined address multiple,
    // if requested at configure-time.
    *cs = FLA_align_ldim( *cs, elem_size );

    // Compute the length of the buffer needed for the object we're creating.
    n_elem = ( size_t ) *cs *
             ( size_t ) n;
  }
  else if ( *cs == 1 )
  {
    // For row-major storage, use rs for computing the length of the buffer
    // to allocate.

    // Align the leading dimension to some user-defined address multiple,
    // if requested at configure-time.
    *rs = FLA_align_ldim( *rs, elem_size );

    // Compute the length of the buffer needed for the object we're creating.
    n_elem = ( size_t ) m *
             ( size_t ) *rs;
  }
  else
  {
    // For general storage, use rs and cs to compute the length of the buffer
    // to allocate.

    // Compute the size of the buffer needed for the object we're creating.
    if ( *rs < *cs )
    {
      *cs = FLA_align_ldim( *cs, elem_size );

      n_elem = ( size_t ) *cs *
               ( size_t ) n;
    }
    else if ( *rs > *cs )
    {
      *rs = FLA_align_ldim( *rs, elem_size );

      n_elem = ( size_t ) m *
               ( size_t ) *rs;
    }
    else // if ( rs == cs )
    {
      //rs = FLA_align_ldim( rs, FLA_Obj_elem_size( *obj ) );
      *cs = FLA_align_ldim( *cs, elem_size );

      // Note that if rs == cs, then we must be creating either a 1-by-n matrix
      // or a m-by-1 matrix. This constraint is enforced in
      // FLA_Check_matrix_strides(). Thus, we can compute the buffer length:
      // m * n * (rs|cs).
      n_elem = ( size_t ) m *
               ( size_t ) n *
               ( size_t ) *cs;
    }
  }

  return n_elem;
}


dim_t FLA_align_ldim( dim_t ldim, dim_t elem_size )
{
#ifdef FLA_ENABLE_MEMORY_ALIGNMENT
  #ifdef FLA_ENABLE_LDIM_ALIGNMENT
    // Increase ldim so that ( ldim * elem_size ) is a multiple of the desired
    // alignment.
    ldim = ( ( ldim * elem_size + FLA_MEMORY_ALIGNMENT_BOUNDARY - 1 ) / 
             FLA_MEMORY_ALIGNMENT_BOUNDARY ) *
           FLA_MEMORY_ALIGNMENT_BOUNDARY /
           elem_size;
  #endif
#endif

  return ldim;
}


void FLA_adjust_strides( dim_t m, dim_t n, dim_t* rs, dim_t* cs )
{
  // Check the strides, and modify them if needed.
  if ( *rs == 0 && *cs == 0 )
  {
    // Default values induce column-major storage, except when m == 1,
    // because we dont want both strides to be unit.
    if ( m == 1 && n > 1 )
    {
      *rs = n;
      *cs = 1;
    }
    else
    {
      *rs = 1;
      *cs = m;
    }
  }
  else if ( *rs == 1 && *cs == 1 )
  {
    // If both strides are unit, this is probably a "lazy" request for a
    // single vector (but could also be a request for a 1xn matrix in column-
    // major order or an mx1 matrix in row-major order). In libflame, we have
    // decided to "reserve" the case where rs == cs == 1 for scalars only, as
    // having unit strides can upset the BLAS error checking when attempting
    // to induce a row-major operation. Also, there is another special case
    // where rs == cs == 1 and one or both of m and n equal zero. This last
    // case is supported to allow creating "empty" objects.

    if ( m == 0 || n == 0 )
    {
      // Nothing needs to be done for the "empty" case where m and/or n
      // equal zero.
    }
    else if ( m == 1 && n == 1 )
    {
      // Nothing needs to be done for the scalar case where m == n == 1.
    }
    else if ( m > 1 && n == 1 )
    {
      // Set the column stride to indicate that this is a column vector stored
      // in column-major order. This is necessary because, in some cases, we
      // have to satisify the error checking in the underlying BLAS library,
      // which expects the leading dimension to be set to at least m, even if
      // it will never be used for indexing since it is a vector and thus only
      // has one column of data.
      *cs = m;
    }
    else if ( m == 1 && n > 1 )
    {
      // Set the row stride to indicate that this is a row vector stored
      // in row-major order.
      *rs = n;
    }
  }
}


FLA_Error FLA_Obj_create_conf_to( FLA_Trans trans, FLA_Obj obj_cur, FLA_Obj *obj_new )
{
  FLA_Datatype datatype;
  FLA_Elemtype elemtype;
  dim_t        m, n;
  dim_t        rs, cs;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_conf_to_check( trans, obj_cur, obj_new );

  datatype = FLA_Obj_datatype( obj_cur );
  elemtype = FLA_Obj_elemtype( obj_cur );

  // Query the dimensions of the existing object.
  if ( trans == FLA_NO_TRANSPOSE || trans == FLA_CONJ_NO_TRANSPOSE )
  {
    m = FLA_Obj_length( obj_cur );
    n = FLA_Obj_width( obj_cur );
  }
  else // if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
  {
    m = FLA_Obj_width( obj_cur );
    n = FLA_Obj_length( obj_cur );
  }

  // Query the row and column strides of the existing object. We don't care
  // about the actual leading dimension of the existing object, only whether
  // it is in row- or column-major format.
  rs = FLA_Obj_row_stride( obj_cur );
  cs = FLA_Obj_col_stride( obj_cur );

  if ( ( rs == 1 && cs == 1 ) )
  {
    // Do nothing. This special case will be handled by FLA_adjust_strides().
    ;
  }
  else if ( rs == 1 )
  {
    // For column-major storage, use the m dimension as the column stride.
    // Row stride is already set to 1.
    cs = m;
  }
  else if ( cs == 1 )
  {
    // For row-major storage, use the n dimension as the row stride.
    // Column stride is already set to 1.
    rs = n;
  }

  // Handle empty views.
  if ( m == 0 ) cs = 1;
  if ( n == 0 ) rs = 1;

  FLA_Obj_create_ext( datatype, elemtype, m, n, m, n, rs, cs, obj_new );

  return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_copy_of( FLA_Trans trans, FLA_Obj obj_cur, FLA_Obj *obj_new )
{
  // Create a new object conformal to the current object.
  FLA_Obj_create_conf_to( trans, obj_cur, obj_new );

#ifdef FLA_ENABLE_SCC
  if ( !FLA_is_owner() )
    return FLA_SUCCESS;
#endif

  // Copy the contents of the current object to the new object.
  FLA_Copyt_external( trans, obj_cur, *obj_new );

  return FLA_SUCCESS;
}


FLA_Error FLA_Obj_create_without_buffer( FLA_Datatype datatype, dim_t m, dim_t n, FLA_Obj *obj )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_without_buffer_check( datatype, m, n, obj );

  // Populate the fields in the view object.
  obj->m                = m;
  obj->n                = n;
  obj->offm             = 0;
  obj->offn             = 0;
  obj->m_inner          = m;
  obj->n_inner          = n;

  // Allocate the base object field.
  obj->base             = ( FLA_Base_obj * ) FLA_malloc( sizeof( FLA_Base_obj ) );

  // Populate the fields in the base object.
  obj->base->datatype   = datatype;
  obj->base->elemtype   = FLA_SCALAR;
  obj->base->m          = m;
  obj->base->n          = n;
  obj->base->m_inner    = m;
  obj->base->n_inner    = n;
  obj->base->id         = ( unsigned long ) obj->base;
  obj->base->m_index    = 0;
  obj->base->n_index    = 0;

  // Set the row and column strides to invalid values.
  obj->base->rs         = 0;
  obj->base->cs         = 0;

  // Initialize the base object's element buffer to NULL.
  obj->base->buffer       = NULL;
  obj->base->buffer_info  = 0;
  obj->base->n_elem_alloc = 0;

#ifdef FLA_ENABLE_SUPERMATRIX
  // Initialize SuperMatrix fields.
  obj->base->n_read_tasks   = 0;
  obj->base->read_task_head = NULL;
  obj->base->read_task_tail = NULL;
  obj->base->write_task     = NULL;
#endif

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_create_constant( double const_real, FLA_Obj *obj )
{
  int*      temp_i;
  float*    temp_s;
  double*   temp_d;
  scomplex* temp_c;
  dcomplex* temp_z;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_constant_check( const_real, obj );

  FLA_Obj_create( FLA_CONSTANT, 1, 1, 0, 0, obj );

#ifdef FLA_ENABLE_SCC
  if ( !FLA_is_owner() )
    return FLA_SUCCESS;
#endif

  temp_i       = FLA_INT_PTR( *obj );
  temp_s       = FLA_FLOAT_PTR( *obj );
  temp_d       = FLA_DOUBLE_PTR( *obj );
  temp_c       = FLA_COMPLEX_PTR( *obj );
  temp_z       = FLA_DOUBLE_COMPLEX_PTR( *obj );

  *temp_i      = ( int   ) const_real;
  *temp_s      = ( float ) const_real;
  *temp_d      =           const_real;
  temp_c->real = ( float ) const_real;
  temp_c->imag = ( float ) 0.0;
  temp_z->real =           const_real;
  temp_z->imag =           0.0;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_create_constant_ext( float const_s, double const_d, FLA_Obj *obj )
{
  int*      temp_i;
  float*    temp_s;
  double*   temp_d;
  scomplex* temp_c;
  dcomplex* temp_z;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_constant_ext_check( const_s, const_d, obj );

  FLA_Obj_create( FLA_CONSTANT, 1, 1, 0, 0, obj );

#ifdef FLA_ENABLE_SCC
  if ( !FLA_is_owner() )
    return FLA_SUCCESS;
#endif

  temp_i       = FLA_INT_PTR( *obj );
  temp_s       = FLA_FLOAT_PTR( *obj );
  temp_d       = FLA_DOUBLE_PTR( *obj );
  temp_c       = FLA_COMPLEX_PTR( *obj );
  temp_z       = FLA_DOUBLE_COMPLEX_PTR( *obj );

  *temp_i      = ( int   ) const_s;
  *temp_s      =           const_s;
  *temp_d      =           const_d;
  temp_c->real =           const_s;
  temp_c->imag =           0.0F;
  temp_z->real =           const_d;
  temp_z->imag =           0.0;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_create_complex_constant( double const_real, double const_imag, FLA_Obj *obj )
{
  int*      temp_i;
  float*    temp_s;
  double*   temp_d;
  scomplex* temp_c;
  dcomplex* temp_z;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_complex_constant_check( const_real, const_imag, obj );

  FLA_Obj_create( FLA_CONSTANT, 1, 1, 0, 0, obj );

#ifdef FLA_ENABLE_SCC
  if ( !FLA_is_owner() )
    return FLA_SUCCESS;
#endif

  temp_i       = FLA_INT_PTR( *obj );
  temp_s       = FLA_FLOAT_PTR( *obj );
  temp_d       = FLA_DOUBLE_PTR( *obj );
  temp_c       = FLA_COMPLEX_PTR( *obj );
  temp_z       = FLA_DOUBLE_COMPLEX_PTR( *obj );

  *temp_i      = ( int   ) const_real;
  *temp_s      = ( float ) const_real;
  *temp_d      =           const_real;
  temp_c->real = ( float ) const_real;
  temp_c->imag = ( float ) const_imag;
  temp_z->real =           const_real;
  temp_z->imag =           const_imag;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_attach_buffer( void *buffer, dim_t rs, dim_t cs, FLA_Obj *obj )
{
  dim_t m, n;

  m = FLA_Obj_length( *obj );
  n = FLA_Obj_width( *obj );

  // Adjust the strides, if necessary.
  FLA_adjust_strides( m, n, &rs, &cs );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_attach_buffer_check( buffer, rs, cs, obj );

  obj->base->buffer      = buffer;
  obj->base->rs          = rs;
  obj->base->cs          = cs;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_create_buffer( dim_t rs, dim_t cs, FLA_Obj *obj )
{
  size_t buffer_size;
  size_t n_elem;
  dim_t  m, n;

  m = FLA_Obj_length( *obj );
  n = FLA_Obj_width( *obj );

  // Adjust the strides, if necessary.
  FLA_adjust_strides( m, n, &rs, &cs );

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_create_buffer_check( rs, cs, obj );

  // Compute the number of elements needed for the buffer, adjusting
  // the strides for alignment if needed.
  n_elem = FLA_compute_num_elem( FLA_Obj_elem_size( *obj ),
                                 m, n, &rs, &cs );

  // Compute the buffer size in bytes.
  buffer_size = ( size_t ) n_elem *
                ( size_t ) FLA_Obj_elem_size( *obj );

  // Allocate the base object's element buffer.
#ifdef FLA_ENABLE_SCC
  obj->base->buffer = ( FLA_Obj_elemtype( *obj ) == FLA_MATRIX ? FLA_malloc( buffer_size ) : FLA_shmalloc( buffer_size ) );
#else
  obj->base->buffer = FLA_malloc( buffer_size );
#endif
  obj->base->buffer_info = 0;

  // Save the number of elements allocated (for use with FLASH).
  obj->base->n_elem_alloc = n_elem;

  // Save the row and column strides used in the memory allocation.
  obj->base->rs     = rs;
  obj->base->cs     = cs;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_free( FLA_Obj *obj )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_free_check( obj );

  if ( obj->base != NULL ) 
  {
#ifdef FLA_ENABLE_SCC
    ( FLA_Obj_elemtype( *obj ) == FLA_MATRIX ? FLA_free( obj->base->buffer ) : FLA_shfree( obj->base->buffer ) );
#else
    //printf( "freeing buff %p\n", obj->base->buffer ); fflush( stdout );
    FLA_free( obj->base->buffer );
#endif
    //printf( "freeing base %p\n", obj->base ); fflush( stdout );
    FLA_free( ( void * ) obj->base );
  }

  obj->offm = 0;
  obj->offn = 0;
  obj->m    = 0;
  obj->n    = 0;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_free_without_buffer( FLA_Obj *obj )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_free_without_buffer_check( obj );

  FLA_free( ( void * ) obj->base );

  obj->offm = 0;
  obj->offn = 0;
  obj->m    = 0;
  obj->n    = 0;

  return FLA_SUCCESS;
}



FLA_Error FLA_Obj_free_buffer( FLA_Obj *obj )
{
  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Obj_free_buffer_check( obj );

#ifdef FLA_ENABLE_SCC
  ( FLA_Obj_elemtype( *obj ) == FLA_MATRIX ? FLA_free( obj->base->buffer ) : FLA_shfree( obj->base->buffer ) );
#else
  FLA_free( obj->base->buffer );
#endif
  obj->base->buffer = NULL;

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_flip_base( FLA_Obj *obj )
{
  FLA_Error e_val;
  dim_t temp;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
  {
    e_val = FLA_Check_null_pointer( obj );
    FLA_Check_error_code( e_val );
    
    e_val = FLA_Check_null_pointer( obj->base );
    FLA_Check_error_code( e_val );
  }

  exchange( obj->base->m,       obj->base->n,       temp );
  exchange( obj->base->cs,      obj->base->rs,      temp );
  exchange( obj->base->m_inner, obj->base->n_inner, temp );
  exchange( obj->base->m_index, obj->base->n_index, temp );

  return FLA_SUCCESS;
}

FLA_Error FLA_Obj_flip_view( FLA_Obj *obj )
{
  FLA_Error e_val;
  dim_t temp;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
  {
    e_val = FLA_Check_null_pointer( obj );
    FLA_Check_error_code( e_val );
  }

  exchange( obj->offm,    obj->offn,    temp );
  exchange( obj->m,       obj->n,       temp );
  exchange( obj->m_inner, obj->n_inner, temp );

  return FLA_SUCCESS;
}
