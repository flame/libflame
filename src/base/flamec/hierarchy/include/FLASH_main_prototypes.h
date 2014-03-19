/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// -----------------------------------------------------------------------------

FLA_Error    FLASH_Obj_blocksizes_check( FLA_Obj H, dim_t* b_m, dim_t* b_n );

FLA_Error    FLASH_Obj_create_helper_check( FLA_Bool without_buffer, FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_hierarchy_check( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* elem_sizes_m, dim_t* elem_sizes_n, FLA_Obj flat_matrix, FLA_Obj* H, unsigned long id, dim_t depth_overall, dim_t* depth_sizes_m, dim_t* depth_sizes_n, dim_t* m_offsets, dim_t* n_offsets );

FLA_Error    FLASH_Obj_create_conf_to_check( FLA_Trans trans, FLA_Obj H_cur, FLA_Obj* H_new );

FLA_Error    FLASH_Obj_create_hier_conf_to_flat_check( FLA_Trans trans, FLA_Obj F, dim_t depth, dim_t* b_mn, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_hier_conf_to_flat_ext_check( FLA_Trans trans, FLA_Obj F, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_flat_conf_to_hier_check( FLA_Trans trans, FLA_Obj H, FLA_Obj* F );
FLA_Error    FLASH_Obj_create_hier_copy_of_flat_check( FLA_Obj F, dim_t depth, dim_t* b_mn, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_hier_copy_of_flat_ext_check( FLA_Obj F, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_flat_copy_of_hier_check( FLA_Obj H, FLA_Obj* F );

FLA_Error    FLASH_Obj_free_check( FLA_Obj* H );
FLA_Error    FLASH_Obj_free_without_buffer_check( FLA_Obj* H );
FLA_Error    FLASH_Obj_free_hierarchy_check( FLA_Obj* H );

FLA_Error    FLASH_Obj_attach_buffer_check( void *buffer, dim_t rs, dim_t cs, FLA_Obj* H );
FLA_Error    FLASH_Obj_attach_buffer_hierarchy_check( FLA_Obj F, FLA_Obj* H );

// -----------------------------------------------------------------------------

FLA_Error FLASH_Part_create_2x1( FLA_Obj A,    FLA_Obj* AT,
                                               FLA_Obj* AB,
                                 dim_t n_rows, FLA_Side side );
FLA_Error FLASH_Part_create_1x2( FLA_Obj A,    FLA_Obj* AL, FLA_Obj* AR,
                                 dim_t n_cols, FLA_Side side );
FLA_Error FLASH_Part_create_2x2( FLA_Obj A,    FLA_Obj* ATL, FLA_Obj* ATR,
                                               FLA_Obj* ABL, FLA_Obj* ABR,
                                 dim_t n_rows, dim_t n_cols, FLA_Side side );

FLA_Error FLASH_Part_free_2x1( FLA_Obj* AT,
                               FLA_Obj* AB );
FLA_Error FLASH_Part_free_1x2( FLA_Obj* AL, FLA_Obj* AR );
FLA_Error FLASH_Part_free_2x2( FLA_Obj* ATL, FLA_Obj* ATR,
                               FLA_Obj* ABL, FLA_Obj* ABR );

FLA_Error FLASH_Obj_adjust_views( FLA_Bool attach_buffer, dim_t offm, dim_t offn, dim_t m, dim_t n, FLA_Obj A, FLA_Obj* S );
FLA_Error FLASH_Obj_adjust_views_hierarchy( FLA_Bool attach_buffer, dim_t offm, dim_t offn, dim_t m, dim_t n, FLA_Obj A, FLA_Obj* S );

dim_t FLASH_Obj_scalar_length( FLA_Obj H );
dim_t FLASH_Obj_scalar_width( FLA_Obj H );
dim_t FLASH_Obj_scalar_min_dim( FLA_Obj H );
dim_t FLASH_Obj_scalar_max_dim( FLA_Obj H );
dim_t FLASH_Obj_scalar_vector_dim( FLA_Obj H );
dim_t FLASH_Obj_scalar_row_offset( FLA_Obj H );
dim_t FLASH_Obj_scalar_col_offset( FLA_Obj H );
dim_t FLASH_Obj_scalar_length_tl( FLA_Obj H );
dim_t FLASH_Obj_scalar_width_tl( FLA_Obj H );
dim_t FLASH_Obj_base_scalar_length( FLA_Obj H );
dim_t FLASH_Obj_base_scalar_width( FLA_Obj H );

FLA_Error FLASH_Obj_show( char* header, FLA_Obj H, char* elem_format, char* footer );
FLA_Error FLASH_Obj_show_hierarchy( FLA_Obj H, dim_t i, char* elem_format );

// -----------------------------------------------------------------------------

FLA_Error    FLASH_Axpy_buffer_to_hier( FLA_Obj alpha, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj H );
FLA_Error    FLASH_Axpy_hier_to_buffer( FLA_Obj alpha, dim_t i, dim_t j, FLA_Obj H, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs );
FLA_Error    FLASH_Axpy_flat_to_hier( FLA_Obj alpha, FLA_Obj F, dim_t i, dim_t j, FLA_Obj H );
FLA_Error    FLASH_Axpy_hier_to_flat( FLA_Obj alpha, dim_t i, dim_t j, FLA_Obj H, FLA_Obj F );

FLA_Error    FLASH_Axpy_hierarchy( int direction, FLA_Obj alpha, FLA_Obj F, FLA_Obj* H );

// -----------------------------------------------------------------------------

FLA_Error    FLASH_Copy_buffer_to_hier( dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj H );
FLA_Error    FLASH_Copy_hier_to_buffer( dim_t i, dim_t j, FLA_Obj H, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs );
FLA_Error    FLASH_Copy_flat_to_hier( FLA_Obj F, dim_t i, dim_t j, FLA_Obj H );
FLA_Error    FLASH_Copy_hier_to_flat( dim_t i, dim_t j, FLA_Obj H, FLA_Obj F );

FLA_Error    FLASH_Copy_hierarchy( int direction, FLA_Obj F, FLA_Obj* H );

// -----------------------------------------------------------------------------

FLA_Datatype FLASH_Obj_datatype( FLA_Obj H );
dim_t        FLASH_Obj_depth( FLA_Obj H );
dim_t        FLASH_Obj_blocksizes( FLA_Obj H, dim_t* b_m, dim_t* b_n );

FLA_Error    FLASH_Obj_create( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_mn, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_ext( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_without_buffer( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_mn, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_without_buffer_ext( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H );

FLA_Error    FLASH_Obj_create_helper( FLA_Bool without_buffer, FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_hierarchy( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* elem_sizes_m, dim_t* elem_sizes_n, FLA_Obj flat_matrix, FLA_Obj* H, unsigned long id, dim_t depth_overall, dim_t* depth_sizes_m, dim_t* depth_sizes_n, dim_t* m_offsets, dim_t* n_offsets );

FLA_Error    FLASH_Obj_create_conf_to( FLA_Trans trans, FLA_Obj H_cur, FLA_Obj* H_new );
FLA_Error    FLASH_Obj_create_hier_conf_to_flat( FLA_Trans trans, FLA_Obj F, dim_t depth, dim_t* b_mn, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_hier_conf_to_flat_ext( FLA_Trans trans, FLA_Obj F, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_flat_conf_to_hier( FLA_Trans trans, FLA_Obj H, FLA_Obj* F );
FLA_Error    FLASH_Obj_create_copy_of( FLA_Trans trans, FLA_Obj H_cur, FLA_Obj* H_new );
FLA_Error    FLASH_Obj_create_hier_copy_of_flat( FLA_Obj F, dim_t depth, dim_t* b_mn, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_hier_copy_of_flat_ext( FLA_Obj F, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H );
FLA_Error    FLASH_Obj_create_flat_copy_of_hier( FLA_Obj H, FLA_Obj* F );

void         FLASH_Obj_free( FLA_Obj* H );
void         FLASH_Obj_free_hierarchy( FLA_Obj* H );
void         FLASH_Obj_free_without_buffer( FLA_Obj* H );

FLA_Error    FLASH_Obj_attach_buffer( void* buffer, dim_t rs, dim_t cs, FLA_Obj* H );
FLA_Error    FLASH_Obj_attach_buffer_hierarchy( FLA_Obj F, FLA_Obj* H );

FLA_Error    FLASH_Obj_flatten( FLA_Obj H, FLA_Obj F );
FLA_Error    FLASH_Obj_hierarchify( FLA_Obj F, FLA_Obj H );

void*        FLASH_Obj_extract_buffer( FLA_Obj H );

FLA_Error    FLASH_Obj_show( char* header, FLA_Obj H, char* elem_format, char* footer );

void         FLASH_print_struct( FLA_Obj H );
void         FLASH_print_struct_helper( FLA_Obj H, int indent );
