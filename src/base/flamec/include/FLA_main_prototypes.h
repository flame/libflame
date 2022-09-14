/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifdef FLA_ENABLE_HIP
#include <rocblas/rocblas.h>
#include <rocsolver/rocsolver.h>
#endif

// -----------------------------------------------------------------------------

fla_blocksize_t* FLA_Blocksize_create( dim_t b_s, dim_t b_d, dim_t b_c, dim_t b_z );
fla_blocksize_t* FLA_Blocksize_create_copy( fla_blocksize_t* bp );
void             FLA_Blocksize_set( fla_blocksize_t* bp, dim_t b_s, dim_t b_d, dim_t b_c, dim_t b_z );
void             FLA_Blocksize_scale( fla_blocksize_t* bp, double factor );
void             FLA_Blocksize_free( fla_blocksize_t* bp );
dim_t            FLA_Blocksize_extract( FLA_Datatype dt, fla_blocksize_t* bp );

fla_blocksize_t* FLA_Query_blocksizes( FLA_Dimension dim );
dim_t            FLA_Query_blocksize( FLA_Datatype dt, FLA_Dimension dim );

dim_t            FLA_Determine_blocksize( FLA_Obj A_unproc, FLA_Quadrant to_dir, fla_blocksize_t* cntl_blocksizes );
dim_t            FLA_determine_matrix_size( FLA_Obj A_unproc, FLA_Quadrant to_dir );



// -----------------------------------------------------------------------------

unsigned int  FLA_Check_error_level( void );
unsigned int  FLA_Check_error_level_set( unsigned int level );
FLA_Error     FLA_Check_error_code_helper( int code, char* file, int line );
FLA_Error     FLA_Check_valid_side( FLA_Side side );
FLA_Error     FLA_Check_valid_uplo( FLA_Uplo uplo );
FLA_Error     FLA_Check_valid_trans( FLA_Trans trans );
FLA_Error     FLA_Check_valid_diag( FLA_Diag diag );
FLA_Error     FLA_Check_valid_conj( FLA_Conj conj );
FLA_Error     FLA_Check_valid_direct( FLA_Conj direct );
FLA_Error     FLA_Check_valid_storev( FLA_Conj storev );
FLA_Error     FLA_Check_valid_inverse( FLA_Inv inv );
FLA_Error     FLA_Check_valid_datatype( FLA_Datatype datatype );
FLA_Error     FLA_Check_valid_object_datatype( FLA_Obj A );
FLA_Error     FLA_Check_valid_evd_type( FLA_Evd_type evd_type );
FLA_Error     FLA_Check_valid_svd_type( FLA_Svd_type svd_type );
FLA_Error     FLA_Check_valid_svd_type_combination( FLA_Svd_type svd_type_u, FLA_Svd_type svd_type_v );
FLA_Error     FLA_Check_valid_svd_type_and_trans_combination( FLA_Svd_type svd_type_u, FLA_Trans transu,
                                                              FLA_Svd_type svd_type_v, FLA_Trans transv );
FLA_Error     FLA_Check_floating_datatype( FLA_Datatype datatype );
FLA_Error     FLA_Check_int_datatype( FLA_Datatype datatype );
FLA_Error     FLA_Check_real_datatype( FLA_Datatype datatype );
FLA_Error     FLA_Check_complex_datatype( FLA_Datatype datatype );
FLA_Error     FLA_Check_floating_object( FLA_Obj A );
FLA_Error     FLA_Check_int_object( FLA_Obj A );
FLA_Error     FLA_Check_real_object( FLA_Obj A );
FLA_Error     FLA_Check_comparable_object( FLA_Obj A );
FLA_Error     FLA_Check_complex_object( FLA_Obj A );
FLA_Error     FLA_Check_consistent_datatype( FLA_Datatype datatype, FLA_Obj A );
FLA_Error     FLA_Check_consistent_object_datatype( FLA_Obj A, FLA_Obj B );
FLA_Error     FLA_Check_identical_object_precision( FLA_Obj A, FLA_Obj B );
FLA_Error     FLA_Check_square( FLA_Obj A );
FLA_Error     FLA_Check_if_scalar( FLA_Obj A );
FLA_Error     FLA_Check_if_vector( FLA_Obj A );
FLA_Error     FLA_Check_conformal_dims( FLA_Trans trans, FLA_Obj A, FLA_Obj B );
FLA_Error     FLA_Check_matrix_matrix_dims( FLA_Trans transa, FLA_Trans transb, FLA_Obj A, FLA_Obj B, FLA_Obj C );
FLA_Error     FLA_Check_matrix_vector_dims( FLA_Trans trans, FLA_Obj A, FLA_Obj x, FLA_Obj y );
FLA_Error     FLA_Check_equal_vector_dims( FLA_Obj x, FLA_Obj y );
FLA_Error     FLA_Check_conj1_trans_and_datatype( FLA_Trans trans, FLA_Obj A );
FLA_Error     FLA_Check_hess_indices( FLA_Obj A, int ilo, int ihi );
FLA_Error     FLA_Check_null_pointer( void* ptr );
FLA_Error     FLA_Check_object_dims( FLA_Trans trans, dim_t m, dim_t n, FLA_Obj A );
FLA_Error     FLA_Check_valid_pivot_type( FLA_Pivot_type ptype );
FLA_Error     FLA_Check_malloc_pointer( void* ptr );
FLA_Error     FLA_Check_base_buffer_mismatch( FLA_Obj A, FLA_Obj B );
FLA_Error     FLA_Check_adjacent_objects_2x2( FLA_Obj A11, FLA_Obj A12,
                                              FLA_Obj A21, FLA_Obj A22 );
FLA_Error     FLA_Check_adjacent_objects_2x1( FLA_Obj AT,
                                              FLA_Obj AB );
FLA_Error     FLA_Check_adjacent_objects_1x2( FLA_Obj AL, FLA_Obj AR );
FLA_Error     FLA_Check_blocksize_value( dim_t b );
FLA_Error     FLA_Check_blocksize_object( FLA_Datatype datatype, fla_blocksize_t* bp );
FLA_Error     FLA_Check_file_descriptor( int fd );
FLA_Error     FLA_Check_lseek_result( int requested_offset, int lseek_r_val );
FLA_Error     FLA_Check_close_result( int close_r_val );
FLA_Error     FLA_Check_unlink_result( int unlink_r_val );
FLA_Error     FLA_Check_read_result( int requested_size, int read_r_val );
FLA_Error     FLA_Check_write_result( int requested_size, int write_r_val );
FLA_Error     FLA_Check_valid_quadrant( FLA_Quadrant quad );
FLA_Error     FLA_Check_vector_dim_min( FLA_Obj x, dim_t min_dim );
FLA_Error     FLA_Check_pthread_create_result( int pthread_create_r_val );
FLA_Error     FLA_Check_pthread_join_result( int pthread_join_r_val );
FLA_Error     FLA_Check_valid_isgn_value( FLA_Obj isgn );
FLA_Error     FLA_Check_sylv_matrix_dims( FLA_Obj A, FLA_Obj B, FLA_Obj C );
FLA_Error     FLA_Check_chol_failure( FLA_Error r_val );
FLA_Error     FLA_Check_valid_elemtype( FLA_Elemtype elemtype );
FLA_Error     FLA_Check_posix_memalign_failure( int r_val );
FLA_Error     FLA_Check_submatrix_dims_and_offset( dim_t m, dim_t n, dim_t i, dim_t j, FLA_Obj A );
FLA_Error     FLA_Check_object_scalar_elemtype( FLA_Obj A );
FLA_Error     FLA_Check_object_matrix_elemtype( FLA_Obj A );
FLA_Error     FLA_Check_num_threads( unsigned int n_threads );
FLA_Error     FLA_Check_conj_and_datatype( FLA_Conj conj, FLA_Obj A );
FLA_Error     FLA_Check_valid_complex_trans( FLA_Trans trans );
FLA_Error     FLA_Check_valid_real_trans( FLA_Trans trans );
FLA_Error     FLA_Check_valid_blas_trans( FLA_Trans trans );
FLA_Error     FLA_Check_nonconstant_datatype( FLA_Datatype datatype );
FLA_Error     FLA_Check_nonconstant_object( FLA_Obj A );
FLA_Error     FLA_Check_identical_object_datatype( FLA_Obj A, FLA_Obj B );
FLA_Error     FLA_Check_divide_by_zero( FLA_Obj alpha );
FLA_Error     FLA_Check_identical_object_elemtype( FLA_Obj A, FLA_Obj B );
FLA_Error     FLA_Check_pivot_index_range( FLA_Obj p, dim_t k1, dim_t k2 );
FLA_Error     FLA_Check_householder_panel_dims( FLA_Obj A, FLA_Obj T );
FLA_Error     FLA_Check_object_length_equals( FLA_Obj A, dim_t m );
FLA_Error     FLA_Check_object_width_equals( FLA_Obj A, dim_t n );
FLA_Error     FLA_Check_object_length_min( FLA_Obj A, dim_t m );
FLA_Error     FLA_Check_object_width_min( FLA_Obj A, dim_t n );
FLA_Error     FLA_Check_valid_error_level( unsigned int level );
FLA_Error     FLA_Check_attempted_repart_2x2( FLA_Obj A_quad, dim_t b_m, dim_t b_n );
FLA_Error     FLA_Check_attempted_repart_2x1( FLA_Obj A_side, dim_t b_m );
FLA_Error     FLA_Check_attempted_repart_1x2( FLA_Obj A_side, dim_t b_n );
FLA_Error     FLA_Check_valid_leftright_side( FLA_Side side );
FLA_Error     FLA_Check_valid_topbottom_side( FLA_Side side );
FLA_Error     FLA_Check_matrix_strides( dim_t m, dim_t n, dim_t rs, dim_t cs );
FLA_Error     FLA_Check_vector_dim( FLA_Obj x, dim_t expected_length );
FLA_Error     FLA_Check_row_vector( FLA_Obj x );
FLA_Error     FLA_Check_col_vector( FLA_Obj x );
FLA_Error     FLA_Check_valid_machval( FLA_Machval val );
FLA_Error     FLA_Check_valid_evd_type( FLA_Evd_type evd_type );
FLA_Error     FLA_Check_valid_svd_type( FLA_Svd_type svd_type );
FLA_Error     FLA_Check_valid_diag_offset( FLA_Obj A, FLA_Diag_off offset );
FLA_Error     FLA_Check_col_storage( FLA_Obj A );
FLA_Error     FLA_Check_row_storage( FLA_Obj A );




// -----------------------------------------------------------------------------

char*         FLA_Error_string_for_code( int code );
void          FLA_Error_messages_init( void );
void          FLA_Print_message( char *str, char *file, int line );
void          FLA_Abort( void );



// -----------------------------------------------------------------------------

void          FLA_Init( void );
void          FLA_Finalize( void );
FLA_Bool      FLA_Initialized( void );

void          FLA_Init_safe( FLA_Error* init_result );
void          FLA_Finalize_safe( FLA_Error init_result );

void          FLA_Init_constants( void );
void          FLA_Finalize_constants( void );

void          FLA_Init_numerical_constants( void );
void          FLA_Finalize_numerical_constants( void );



//------------------------------------------------------------------------------

void          FLA_Lock_init( FLA_Lock* fla_lock_ptr );
void          FLA_Lock_destroy( FLA_Lock* fla_lock_ptr );
void          FLA_Lock_acquire( FLA_Lock* fla_lock_ptr );
void          FLA_Lock_release( FLA_Lock* fla_lock_ptr );



// -----------------------------------------------------------------------------

void          FLA_Memory_leak_counter_init( void );
void          FLA_Memory_leak_counter_finalize( void );
FLA_Bool      FLA_Memory_leak_counter_status( void );
FLA_Bool      FLA_Memory_leak_counter_set( FLA_Bool new_status );

void*         FLA_malloc( size_t size );
void*         FLA_realloc( void* old_ptr, size_t size );
void*         FLA_buff_malloc( size_t size );
void          FLA_free( void *ptr );
void          FLA_buff_free( void *ptr );
 


// -----------------------------------------------------------------------------

FLA_Error     FLA_Obj_copy_view( FLA_Obj A, FLA_Obj* B );
void          FLA_Obj_extract_real_scalar( FLA_Obj alpha, double* alpha_value );
void          FLA_Obj_extract_complex_scalar( FLA_Obj alpha, dcomplex* alpha_value );
void          FLA_Obj_extract_real_part( FLA_Obj alpha, FLA_Obj beta );
void          FLA_Obj_extract_imag_part( FLA_Obj alpha, FLA_Obj beta );
void          FLA_Obj_set_real_part( FLA_Obj alpha, FLA_Obj beta );
void          FLA_Obj_set_imag_part( FLA_Obj alpha, FLA_Obj beta );
FLA_Error     FLA_Obj_show( char *s1, FLA_Obj A, char *format, char *s2 );
FLA_Error     FLA_Obj_fshow( FILE* file, char *s1, FLA_Obj A, char *format, char *s2 );

FLA_Error     FLA_Obj_copy_view_check( FLA_Obj A, FLA_Obj* B );
FLA_Error     FLA_Obj_extract_real_scalar_check( FLA_Obj alpha, double* alpha_value );
FLA_Error     FLA_Obj_extract_complex_scalar_check( FLA_Obj alpha, dcomplex* alpha_value );
FLA_Error     FLA_Obj_extract_real_part_check( FLA_Obj alpha, FLA_Obj beta );
FLA_Error     FLA_Obj_extract_imag_part_check( FLA_Obj alpha, FLA_Obj beta );
FLA_Error     FLA_Obj_set_real_part_check( FLA_Obj alpha, FLA_Obj beta );
FLA_Error     FLA_Obj_set_imag_part_check( FLA_Obj alpha, FLA_Obj beta );
FLA_Error     FLA_Obj_show_check( char* s1, FLA_Obj obj, char* format, char* s2 );
FLA_Error     FLA_Obj_fshow_check( FILE* file, char* s1, FLA_Obj obj, char* format, char* s2 );


// -----------------------------------------------------------------------------

FLA_Error     FLA_Copy_buffer_to_object( FLA_Trans trans, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj obj );
FLA_Error     FLA_Copy_object_to_buffer( FLA_Trans trans, dim_t i, dim_t j, FLA_Obj obj, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs );
FLA_Error     FLA_Copy_buffer_to_object_check( FLA_Trans trans, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj obj );
FLA_Error     FLA_Copy_object_to_buffer_check( FLA_Trans trans, dim_t i, dim_t j, FLA_Obj obj, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs );



// -----------------------------------------------------------------------------

FLA_Error     FLA_Axpy_buffer_to_object( FLA_Trans trans, FLA_Obj alpha, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj C );
FLA_Error     FLA_Axpy_object_to_buffer( FLA_Trans trans, FLA_Obj alpha, dim_t i, dim_t j, FLA_Obj C, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs );

FLA_Error     FLA_Axpy_buffer_to_object_check( FLA_Trans trans, FLA_Obj alpha, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj C );
FLA_Error     FLA_Axpy_object_to_buffer_check( FLA_Trans trans, FLA_Obj alpha, dim_t i, dim_t j, FLA_Obj C, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs );



// -----------------------------------------------------------------------------

#ifdef FLA_ENABLE_SCC
void*         FLA_shmalloc( size_t size );
void          FLA_shfree( void* ptr );
FLA_Bool      FLA_is_owner( void );
#endif
FLA_Error     FLA_Obj_nullify( FLA_Obj *obj );
FLA_Error     FLA_Obj_create( FLA_Datatype datatype, dim_t m, dim_t n, dim_t rs, dim_t cs, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_ext( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t m, dim_t n, dim_t m_inner, dim_t n_inner, dim_t rs, dim_t cs, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_conf_to( FLA_Trans trans, FLA_Obj old, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_copy_of( FLA_Trans trans, FLA_Obj old, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_without_buffer( FLA_Datatype datatype, dim_t m, dim_t n, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_constant( double const_real, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_constant_ext( float const_s, double const_d, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_complex_constant( double const_real, double const_imag, FLA_Obj *obj );
FLA_Error     FLA_Obj_attach_buffer( void *buffer, dim_t rs, dim_t cs, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_buffer( dim_t rs, dim_t cs, FLA_Obj *obj );
FLA_Error     FLA_Obj_free( FLA_Obj *obj );
FLA_Error     FLA_Obj_free_without_buffer( FLA_Obj *obj );
FLA_Error     FLA_Obj_free_buffer( FLA_Obj *obj );
dim_t         FLA_align_ldim( dim_t ldim, dim_t elem_size );
dim_t         FLA_compute_num_elem( dim_t elem_size, dim_t m, dim_t n, dim_t* rs, dim_t* cs );
void          FLA_adjust_strides( dim_t m, dim_t n, dim_t* rs, dim_t* cs );

FLA_Error     FLA_Obj_flip_base( FLA_Obj *obj );
FLA_Error     FLA_Obj_flip_view( FLA_Obj *obj );

FLA_Error     FLA_Obj_create_ext_check( FLA_Datatype datatype, FLA_Elemtype elemtype, dim_t m, dim_t n, dim_t m_inner, dim_t n_inner, dim_t rs, dim_t cs, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_conf_to_check( FLA_Trans trans, FLA_Obj obj_old, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_without_buffer_check( FLA_Datatype datatype, dim_t m, dim_t n, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_constant_check( double const_real, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_constant_ext_check( float const_s, double const_d, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_complex_constant_check( double const_real, double const_imag, FLA_Obj *obj );
FLA_Error     FLA_Obj_attach_buffer_check( void *buffer, dim_t rs, dim_t cs, FLA_Obj *obj );
FLA_Error     FLA_Obj_create_buffer_check( dim_t rs, dim_t cs, FLA_Obj *obj );
FLA_Error     FLA_Obj_free_check( FLA_Obj *obj );
FLA_Error     FLA_Obj_free_without_buffer_check( FLA_Obj *obj );
FLA_Error     FLA_Obj_free_buffer_check( FLA_Obj *obj );

FLA_Error     FLA_Obj_create_buffer_task( dim_t rs, dim_t cs, FLA_Obj obj, void* cntl );
FLA_Error     FLA_Obj_free_buffer_task( FLA_Obj obj, void* cntl );


// -----------------------------------------------------------------------------

FLA_Datatype  FLA_Obj_datatype( FLA_Obj obj );
FLA_Datatype  FLA_Obj_datatype_proj_to_real( FLA_Obj A );
FLA_Datatype  FLA_Obj_datatype_proj_to_complex( FLA_Obj A );
FLA_Elemtype  FLA_Obj_elemtype( FLA_Obj obj );
dim_t         FLA_Obj_datatype_size( FLA_Datatype datatype );
dim_t         FLA_Obj_elem_size( FLA_Obj obj );
dim_t         FLA_Obj_length( FLA_Obj obj );
dim_t         FLA_Obj_width( FLA_Obj obj );
FLA_Uplo      FLA_Obj_structure( FLA_Obj obj );
dim_t         FLA_Obj_vector_dim( FLA_Obj obj );
dim_t         FLA_Obj_vector_inc( FLA_Obj obj );
dim_t         FLA_Obj_min_dim( FLA_Obj obj );
dim_t         FLA_Obj_max_dim( FLA_Obj obj );
dim_t         FLA_Obj_row_stride( FLA_Obj obj );
dim_t         FLA_Obj_col_stride( FLA_Obj obj );
dim_t         FLA_Obj_row_offset( FLA_Obj obj );
dim_t         FLA_Obj_col_offset( FLA_Obj obj );
dim_t         FLA_Obj_base_length( FLA_Obj obj );
dim_t         FLA_Obj_base_width( FLA_Obj obj );
dim_t         FLA_Obj_num_elem_alloc( FLA_Obj obj );
void*         FLA_Obj_base_buffer( FLA_Obj obj );
void*         FLA_Obj_buffer_at_view( FLA_Obj obj );
FLA_Bool      FLA_Obj_buffer_is_null( FLA_Obj obj );
FLA_Bool      FLA_Obj_is_int( FLA_Obj A );
FLA_Bool      FLA_Obj_is_floating_point( FLA_Obj A );
FLA_Bool      FLA_Obj_is_constant( FLA_Obj A );
FLA_Bool      FLA_Obj_is_real( FLA_Obj A );
FLA_Bool      FLA_Obj_is_complex( FLA_Obj A );
FLA_Bool      FLA_Obj_is_single_precision( FLA_Obj A );
FLA_Bool      FLA_Obj_is_double_precision( FLA_Obj A );
FLA_Bool      FLA_Obj_is_scalar( FLA_Obj A );
FLA_Bool      FLA_Obj_is_vector( FLA_Obj A );
FLA_Bool      FLA_Obj_has_zero_dim( FLA_Obj A );
FLA_Bool      FLA_Obj_is_row_major( FLA_Obj A );
FLA_Bool      FLA_Obj_is_col_major( FLA_Obj A );
FLA_Bool      FLA_Obj_is_conformal_to( FLA_Trans trans, FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_is( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_is_identical( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_is_overlapped( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_equals( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_gt( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_ge( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_lt( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_le( FLA_Obj A, FLA_Obj B );
void*         FLA_Submatrix_at( FLA_Datatype datatype, void* buffer, dim_t i, dim_t j, dim_t rs, dim_t cs );
FLA_Bool      FLA_Obj_has_nan( FLA_Obj A );

FLA_Error     FLA_Obj_datatype_check( FLA_Obj obj );
FLA_Error     FLA_Obj_datatype_proj_to_real_check( FLA_Obj obj );
FLA_Error     FLA_Obj_elemtype_check( FLA_Obj obj );
FLA_Error     FLA_Obj_datatype_size_check( FLA_Datatype datatype );
FLA_Error     FLA_Obj_elem_size_check( FLA_Obj obj );
FLA_Error     FLA_Obj_buffer_at_view_check( FLA_Obj obj );
FLA_Error     FLA_Obj_equals_check( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_gt_check( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_ge_check( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_lt_check( FLA_Obj A, FLA_Obj B );
FLA_Bool      FLA_Obj_le_check( FLA_Obj A, FLA_Obj B );
FLA_Error     FLA_Submatrix_at_check( FLA_Datatype datatype, void* buffer, dim_t i, dim_t j, dim_t rs, dim_t cs );
FLA_Error     FLA_Obj_has_nan_check( FLA_Obj A );


// -----------------------------------------------------------------------------

void          FLA_Param_map_flame_to_netlib_trans( FLA_Trans trans, void* blas_trans );
void          FLA_Param_map_flame_to_netlib_uplo( FLA_Uplo uplo, void* blas_uplo );
void          FLA_Param_map_flame_to_netlib_side( FLA_Uplo side, void* blas_side );
void          FLA_Param_map_flame_to_netlib_diag( FLA_Diag diag, void* blas_diag );
void          FLA_Param_map_flame_to_netlib_direct( FLA_Direct direct, void* lapack_direct );
void          FLA_Param_map_flame_to_netlib_storev( FLA_Store storev, void* lapack_storev );
void          FLA_Param_map_flame_to_netlib_evd_type( FLA_Evd_type evd_type, void* lapack_evd_type );
void          FLA_Param_map_flame_to_netlib_svd_type( FLA_Svd_type svd_type, void* lapack_svd_type );
void          FLA_Param_map_flame_to_netlib_machval( FLA_Machval machval, void* blas_machval );

#ifdef FLA_ENABLE_HIP
rocblas_operation FLA_Param_map_flame_to_rocblas_trans( FLA_Trans trans, FLA_Bool is_real );
rocblas_fill      FLA_Param_map_flame_to_rocblas_uplo( FLA_Uplo uplo );
rocblas_side      FLA_Param_map_flame_to_rocblas_side( FLA_Side side );
rocblas_diagonal  FLA_Param_map_flame_to_rocblas_diag( FLA_Diag diag );
rocblas_evect     FLA_Param_map_flame_to_rocblas_evd_type( FLA_Evd_type evd_type );
rocblas_svect     FLA_Param_map_flame_to_rocblas_svd_type( FLA_Svd_type svd_type );
#endif

void          FLA_Param_map_flame_to_blis_trans( FLA_Trans trans, trans1_t* blis_trans );
void          FLA_Param_map_flame_to_blis_conj( FLA_Conj conj, conj1_t* blis_conj );
void          FLA_Param_map_flame_to_blis_uplo( FLA_Uplo uplo, uplo1_t* blis_uplo );
void          FLA_Param_map_flame_to_blis_side( FLA_Uplo side, side1_t* blis_side );
void          FLA_Param_map_flame_to_blis_diag( FLA_Diag diag, diag1_t* blis_diag );
#if 0
void          FLA_Param_map_flame_to_blis2_trans( FLA_Trans trans, trans_t* blis_trans );
void          FLA_Param_map_flame_to_blis2_conj( FLA_Conj conj, conj_t* blis_conj );
void          FLA_Param_map_flame_to_blis2_uplo( FLA_Uplo uplo, uplo_t* blis_uplo );
void          FLA_Param_map_flame_to_blis2_side( FLA_Uplo side, side_t* blis_side );
void          FLA_Param_map_flame_to_blis2_diag( FLA_Diag diag, diag_t* blis_diag );
#endif

void          FLA_Param_map_blis_to_flame_trans( trans1_t trans, FLA_Trans* flame_trans );
void          FLA_Param_map_blis_to_flame_uplo( uplo1_t uplo, FLA_Uplo* flame_uplo );
void          FLA_Param_map_blis_to_flame_side( side1_t side, FLA_Side* flame_side );
void          FLA_Param_map_blis_to_flame_diag( diag1_t diag, FLA_Diag* flame_diag );

void          FLA_Param_map_char_to_flame_trans( char* trans, FLA_Trans* flame_trans );
void          FLA_Param_map_char_to_flame_uplo( char* uplo, FLA_Uplo* flame_uplo );
void          FLA_Param_map_char_to_flame_side( char* side, FLA_Side* flame_side );
void          FLA_Param_map_char_to_flame_diag( char* diag, FLA_Diag* flame_diag );
void          FLA_Param_map_char_to_flame_storev( char* storev, FLA_Direct* flame_storev );
void          FLA_Param_map_char_to_flame_direct( char* direct, FLA_Direct* flame_direct );
void          FLA_Param_map_char_to_flame_inv( char* inv, FLA_Inv* flame_inv );

/*
void          FLA_Param_map_blis_to_netlib_trans( char blis_trans, void* blas_trans );
void          FLA_Param_map_blis_to_netlib_uplo( char blis_uplo, void* blas_uplo );
void          FLA_Param_map_blis_to_netlib_side( char blis_side, void* blas_side );
void          FLA_Param_map_blis_to_netlib_diag( char blis_diag, void* blas_diag );
*/

void          FLA_Param_map_netlib_to_flame_trans( char* trans, FLA_Trans* flame_trans );
void          FLA_Param_map_netlib_to_flame_uplo( char* uplo, FLA_Uplo* flame_uplo );
void          FLA_Param_map_netlib_to_flame_side( char* side, FLA_Side* flame_side );
void          FLA_Param_map_netlib_to_flame_diag( char* diag, FLA_Diag* flame_diag );
void          FLA_Param_map_netlib_to_flame_inv( int* itype, FLA_Inv* flame_inv );
void          FLA_Param_map_netlib_to_flame_svd_type( char* svd, FLA_Svd_type* flame_svd );




// -----------------------------------------------------------------------------

FLA_Error     FLA_Part_2x2( FLA_Obj A,  FLA_Obj *A11, FLA_Obj *A12,
                                        FLA_Obj *A21, FLA_Obj *A22,
                            dim_t  mb,  dim_t     nb, FLA_Quadrant quadrant );

FLA_Error     FLA_Part_2x1 ( FLA_Obj A,  FLA_Obj *A1,
                                         FLA_Obj *A2,
                             dim_t  mb,  FLA_Side side );

FLA_Error     FLA_Part_1x2( FLA_Obj A,  FLA_Obj *A1, FLA_Obj *A2,
                                        dim_t    nb, FLA_Side side );
 
FLA_Error     FLA_Merge_2x2( FLA_Obj A11, FLA_Obj A12,
                             FLA_Obj A21, FLA_Obj A22,  FLA_Obj *A );
 
FLA_Error     FLA_Merge_2x1( FLA_Obj AT,
                             FLA_Obj AB,  FLA_Obj *A );

FLA_Error     FLA_Merge_1x2( FLA_Obj AL, FLA_Obj AR,  FLA_Obj *A );

FLA_Error     FLA_Repart_2x2_to_3x3( FLA_Obj ATL, FLA_Obj ATR,  FLA_Obj *A00, FLA_Obj *A01, FLA_Obj *A02,
                                                                FLA_Obj *A10, FLA_Obj *A11, FLA_Obj *A12,
                                     FLA_Obj ABL, FLA_Obj ABR,  FLA_Obj *A20, FLA_Obj *A21, FLA_Obj *A22,
                                     dim_t   mb,  dim_t    nb,  FLA_Quadrant quadrant );

FLA_Error     FLA_Repart_2x1_to_3x1( FLA_Obj AT,  FLA_Obj *A0,
                                                  FLA_Obj *A1,
                                     FLA_Obj AB,  FLA_Obj *A2,
                                     dim_t   mb,  FLA_Side side );

FLA_Error     FLA_Repart_1x2_to_1x3( FLA_Obj  AL,              FLA_Obj  AR,
                                     FLA_Obj *A0, FLA_Obj *A1, FLA_Obj *A2,
                                                  dim_t    nb, FLA_Side side );

FLA_Error     FLA_Cont_with_3x3_to_2x2( FLA_Obj *ATL, FLA_Obj *ATR,  FLA_Obj A00, FLA_Obj A01, FLA_Obj A02,
                                                                     FLA_Obj A10, FLA_Obj A11, FLA_Obj A12,
                                        FLA_Obj *ABL, FLA_Obj *ABR,  FLA_Obj A20, FLA_Obj A21, FLA_Obj A22,
                                                                     FLA_Quadrant quadrant );

FLA_Error     FLA_Cont_with_3x1_to_2x1( FLA_Obj *AT,  FLA_Obj A0,
                                                      FLA_Obj A1,
                                        FLA_Obj *AB,  FLA_Obj A2,
                                                      FLA_Side side );

FLA_Error     FLA_Cont_with_1x3_to_1x2( FLA_Obj *AL,              FLA_Obj *AR,
                                        FLA_Obj  A0, FLA_Obj  A1, FLA_Obj  A2,
                                                                  FLA_Side side );

FLA_Error     FLA_Repart_3x3_to_5x5( FLA_Obj ATL, FLA_Obj ATM, FLA_Obj ATR,
                                     FLA_Obj AML, FLA_Obj AMM, FLA_Obj AMR,
                                     FLA_Obj ABL, FLA_Obj ABM, FLA_Obj ABR,
                                     FLA_Obj *A00, FLA_Obj *A01, FLA_Obj *A02, FLA_Obj *A03, FLA_Obj *A04,
                                     FLA_Obj *A10, FLA_Obj *A11, FLA_Obj *A12, FLA_Obj *A13, FLA_Obj *A14,
                                     FLA_Obj *A20, FLA_Obj *A21, FLA_Obj *A22, FLA_Obj *A23, FLA_Obj *A24,
                                     FLA_Obj *A30, FLA_Obj *A31, FLA_Obj *A32, FLA_Obj *A33, FLA_Obj *A34,
                                     FLA_Obj *A40, FLA_Obj *A41, FLA_Obj *A42, FLA_Obj *A43, FLA_Obj *A44,
                                     dim_t b, FLA_Quadrant quadrant );

FLA_Error     FLA_Cont_with_5x5_to_3x3( FLA_Obj *ATL, FLA_Obj *ATM, FLA_Obj *ATR,
                                        FLA_Obj *AML, FLA_Obj *AMM, FLA_Obj *AMR,
                                        FLA_Obj *ABL, FLA_Obj *ABM, FLA_Obj *ABR,
                                        FLA_Obj A00, FLA_Obj A01, FLA_Obj A02, FLA_Obj A03, FLA_Obj A04,
                                        FLA_Obj A10, FLA_Obj A11, FLA_Obj A12, FLA_Obj A13, FLA_Obj A14,
                                        FLA_Obj A20, FLA_Obj A21, FLA_Obj A22, FLA_Obj A23, FLA_Obj A24,
                                        FLA_Obj A30, FLA_Obj A31, FLA_Obj A32, FLA_Obj A33, FLA_Obj A34,
                                        FLA_Obj A40, FLA_Obj A41, FLA_Obj A42, FLA_Obj A43, FLA_Obj A44,
                                        FLA_Quadrant quadrant );



FLA_Error     FLA_Part_2x2_check( FLA_Obj A,  FLA_Obj *A11, FLA_Obj *A12,
                                              FLA_Obj *A21, FLA_Obj *A22,
                                  dim_t  mb,  dim_t     nb, FLA_Quadrant quadrant );

FLA_Error     FLA_Part_2x1_check( FLA_Obj A,  FLA_Obj *A1,
                                               FLA_Obj *A2,
                                   dim_t  mb,  FLA_Side side );

FLA_Error     FLA_Part_1x2_check( FLA_Obj A,  FLA_Obj *A1, FLA_Obj *A2,
                                              dim_t    nb, FLA_Side side );
 
FLA_Error     FLA_Merge_2x2_check( FLA_Obj A11, FLA_Obj A12,
                                   FLA_Obj A21, FLA_Obj A22,  FLA_Obj *A );
 
FLA_Error     FLA_Merge_2x1_check( FLA_Obj AT,
                                   FLA_Obj AB,  FLA_Obj *A );

FLA_Error     FLA_Merge_1x2_check( FLA_Obj AL, FLA_Obj AR,  FLA_Obj *A );

FLA_Error     FLA_Repart_2x2_to_3x3_check( FLA_Obj ATL, FLA_Obj ATR,  FLA_Obj *A00, FLA_Obj *A01, FLA_Obj *A02,
                                                                      FLA_Obj *A10, FLA_Obj *A11, FLA_Obj *A12,
                                           FLA_Obj ABL, FLA_Obj ABR,  FLA_Obj *A20, FLA_Obj *A21, FLA_Obj *A22,
                                           dim_t   mb,  dim_t    nb,  FLA_Quadrant quadrant );

FLA_Error     FLA_Repart_2x1_to_3x1_check( FLA_Obj AT,  FLA_Obj *A0,
                                                        FLA_Obj *A1,
                                           FLA_Obj AB,  FLA_Obj *A2,
                                           dim_t   mb,  FLA_Side side );

FLA_Error     FLA_Repart_1x2_to_1x3_check( FLA_Obj  AL,              FLA_Obj  AR,
                                           FLA_Obj *A0, FLA_Obj *A1, FLA_Obj *A2,
                                                        dim_t    nb, FLA_Side side );

FLA_Error     FLA_Cont_with_3x3_to_2x2_check( FLA_Obj *ATL, FLA_Obj *ATR,  FLA_Obj A00, FLA_Obj A01, FLA_Obj A02,
                                                                           FLA_Obj A10, FLA_Obj A11, FLA_Obj A12,
                                              FLA_Obj *ABL, FLA_Obj *ABR,  FLA_Obj A20, FLA_Obj A21, FLA_Obj A22,
                                                                           FLA_Quadrant quadrant );

FLA_Error     FLA_Cont_with_3x1_to_2x1_check( FLA_Obj *AT,  FLA_Obj A0,
                                                            FLA_Obj A1,
                                              FLA_Obj *AB,  FLA_Obj A2,
                                                            FLA_Side side );

FLA_Error     FLA_Cont_with_1x3_to_1x2_check( FLA_Obj *AL,              FLA_Obj *AR,
                                              FLA_Obj  A0, FLA_Obj  A1, FLA_Obj  A2,
                                                                        FLA_Side side );
