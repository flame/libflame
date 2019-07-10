/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

float     FLA_random_float( void );
double    FLA_random_double( void );
scomplex  FLA_random_scomplex( void );
dcomplex  FLA_random_dcomplex( void );

FLA_Error FLA_Absolute_square( FLA_Obj alpha );
FLA_Error FLA_Absolute_value( FLA_Obj alpha );
double    FLA_Clock( void );
FLA_Error FLA_Conjugate( FLA_Obj A );
FLA_Error FLA_Conjugate_r( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Fill_with_linear_dist( FLA_Obj shift, FLA_Obj delta, FLA_Obj x );
FLA_Error FLA_Fill_with_inverse_dist( FLA_Obj alpha, FLA_Obj x );
FLA_Error FLA_Fill_with_geometric_dist( FLA_Obj alpha, FLA_Obj x );
FLA_Error FLA_Fill_with_random_dist( FLA_Obj shift, FLA_Obj max, FLA_Obj x );
FLA_Error FLA_Fill_with_logarithmic_dist( FLA_Obj max, FLA_Obj x );
FLA_Error FLA_Fill_with_cluster_dist( FLA_Obj n_clusters, FLA_Obj cluster_width, FLA_Obj x );
FLA_Error FLA_Hermitianize( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Invert( FLA_Conj conj, FLA_Obj x );
FLA_Error FLA_Inv_scal_elemwise( FLA_Trans trans, FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Max_abs_value( FLA_Obj A, FLA_Obj amax );
FLA_Error FLA_Max_abs_value_herm( FLA_Uplo uplo, FLA_Obj A, FLA_Obj maxabs );
double    FLA_Max_elemwise_diff( FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Mult_add( FLA_Obj alpha, FLA_Obj beta, FLA_Obj gamma );
FLA_Error FLA_Negate( FLA_Obj x );
FLA_Error FLA_Norm1( FLA_Obj A, FLA_Obj norm );
FLA_Error FLA_Norm_inf( FLA_Obj A, FLA_Obj norm );
FLA_Error FLA_Norm_frob( FLA_Obj A, FLA_Obj norm );
FLA_Error FLA_Pow( FLA_Obj base, FLA_Obj exp, FLA_Obj btoe );
FLA_Error FLA_Random_matrix( FLA_Obj A );
FLA_Error FLA_Random_herm_matrix( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Random_symm_matrix( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Random_spd_matrix( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Random_tri_matrix( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLA_Random_unitary_matrix( FLA_Obj A );
FLA_Error FLA_Scal_elemwise( FLA_Trans trans, FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Setr( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Shift_pivots_to_check( FLA_Pivot_type ptype, FLA_Obj p );
FLA_Error FLA_Sqrt( FLA_Obj alpha );
FLA_Error FLA_Symmetrize( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Triangularize( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLA_Transpose( FLA_Obj A );

FLA_Error FLA_Set( FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Set_diag( FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Set_offdiag( int offset, FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Set_to_identity( FLA_Obj A );
FLA_Error FLA_Add_to_diag( void *diag_value, FLA_Obj A );
FLA_Error FLA_Shift_diag( FLA_Conj conj, FLA_Obj sigma, FLA_Obj A );
FLA_Error FLA_Scale_diag( FLA_Conj conj, FLA_Obj alpha, FLA_Obj A );

FLA_Error FLA_Set_diagonal_vector( FLA_Obj A, FLA_Obj d );
FLA_Error FLA_Set_diagonal_matrix( FLA_Obj d, FLA_Obj A );

// -----------------------------------------------------------------------------

FLA_Error FLA_Absolute_square_check( FLA_Obj alpha );
FLA_Error FLA_Absolute_value_check( FLA_Obj alpha );
FLA_Error FLA_Conjugate_check( FLA_Obj A );
FLA_Error FLA_Conjugate_r_check( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Fill_with_linear_dist_check( FLA_Obj shift, FLA_Obj delta, FLA_Obj x );
FLA_Error FLA_Fill_with_inverse_dist_check( FLA_Obj alpha, FLA_Obj x );
FLA_Error FLA_Fill_with_geometric_dist_check( FLA_Obj alpha, FLA_Obj x );
FLA_Error FLA_Fill_with_random_dist_check( FLA_Obj shift, FLA_Obj max, FLA_Obj x );
FLA_Error FLA_Fill_with_logarithmic_dist_check( FLA_Obj alpha, FLA_Obj x );
FLA_Error FLA_Fill_with_cluster_dist_check( FLA_Obj n_clusters, FLA_Obj cluster_width, FLA_Obj x );
FLA_Error FLA_Hermitianize_check( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Invert_check( FLA_Conj conj, FLA_Obj x );
FLA_Error FLA_Inv_scal_elemwise_check( FLA_Trans trans, FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Max_abs_value_check( FLA_Obj A, FLA_Obj amax );
FLA_Error FLA_Max_abs_value_herm_check( FLA_Uplo uplo, FLA_Obj A, FLA_Obj maxabs );
FLA_Error FLA_Max_elemwise_diff_check( FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Mult_add_check( FLA_Obj alpha, FLA_Obj beta, FLA_Obj gamma );
FLA_Error FLA_Negate_check( FLA_Obj x );
FLA_Error FLA_Norm1_check( FLA_Obj A, FLA_Obj norm );
FLA_Error FLA_Norm_inf_check( FLA_Obj A, FLA_Obj norm );
FLA_Error FLA_Norm_frob_check( FLA_Obj A, FLA_Obj norm );
FLA_Error FLA_Pow_check( FLA_Obj base, FLA_Obj exp, FLA_Obj btoe );
FLA_Error FLA_Random_matrix_check( FLA_Obj A );
FLA_Error FLA_Random_herm_matrix_check( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Random_symm_matrix_check( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Random_spd_matrix_check( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Random_tri_matrix_check( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLA_Random_unitary_matrix_check( FLA_Obj A );
FLA_Error FLA_Scal_elemwise_check( FLA_Trans trans, FLA_Obj A, FLA_Obj B );
FLA_Error FLA_Setr_check( FLA_Uplo uplo, FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Sort_check( FLA_Direct direct, FLA_Obj x );
FLA_Error FLA_Sqrt_check( FLA_Obj alpha );
FLA_Error FLA_Symmetrize_check( FLA_Uplo uplo, FLA_Obj A );
FLA_Error FLA_Triangularize_check( FLA_Uplo uplo, FLA_Diag diag, FLA_Obj A );
FLA_Error FLA_Transpose_check( FLA_Obj A );

FLA_Error FLA_Set_check( FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Set_diag_check( FLA_Obj alpha, FLA_Obj A );
FLA_Error FLA_Set_to_identity_check( FLA_Obj A );
FLA_Error FLA_Add_to_diag_check( void *diag_value, FLA_Obj A );
FLA_Error FLA_Shift_diag_check( FLA_Conj conj, FLA_Obj sigma, FLA_Obj A );
FLA_Error FLA_Scale_diag_check( FLA_Conj conj, FLA_Obj alpha, FLA_Obj A );

// -----------------------------------------------------------------------------

FLA_Error FLA_Transpose_blk_var1( FLA_Obj A, fla_tpose_t* cntl );
FLA_Error FLA_Transpose_blk_var2( FLA_Obj A, fla_tpose_t* cntl );
FLA_Error FLA_Transpose_unb_var1( FLA_Obj A );
FLA_Error FLA_Transpose_unb_var2( FLA_Obj A );
FLA_Error FLA_Swap_t_blk_var1( FLA_Obj A, FLA_Obj B, fla_swap_t* cntl );
FLA_Error FLA_Swap_t_blk_var2( FLA_Obj A, FLA_Obj B, fla_swap_t* cntl );

FLA_Error FLA_Sort( FLA_Direct direct, FLA_Obj x );
FLA_Error FLA_Sort_f_ops( int     m_x,
                          float*  x, int inc_x );
FLA_Error FLA_Sort_b_ops( int     m_x,
                          float*  x, int inc_x );
FLA_Error FLA_Sort_f_opd( int     m_x,
                          double* x, int inc_x );
FLA_Error FLA_Sort_b_opd( int     m_x,
                          double* x, int inc_x );

