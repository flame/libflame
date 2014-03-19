/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

// --- Query routine prototypes ------------------------------------------------

// --- trans ---

int bl1_does_trans( trans1_t trans );
int bl1_does_notrans( trans1_t trans );
int bl1_does_conj( trans1_t trans );

int bl1_is_notrans( trans1_t trans );
int bl1_is_trans( trans1_t trans );
int bl1_is_conjnotrans( trans1_t trans );
int bl1_is_conjtrans( trans1_t trans );

// --- conj ---

int bl1_is_noconj( conj1_t conj );
int bl1_is_conj( conj1_t conj );

// --- uplo ---

int bl1_is_lower( uplo1_t uplo );
int bl1_is_upper( uplo1_t uplo );

// --- side ---

int bl1_is_left( side1_t side );
int bl1_is_right( side1_t side );

// --- diag ---

int bl1_is_nonunit_diag( diag1_t diag );
int bl1_is_unit_diag( diag1_t diag );
int bl1_is_zero_diag( diag1_t diag );

// --- mapping-related ---

conj1_t bl1_proj_trans1_to_conj( trans1_t trans );

// --- storage-related ---

void bl1_check_storage_3m( int a_rs, int a_cs, int b_rs, int b_cs, int c_rs, int c_cs );
void bl1_check_storage_2m( int a_rs, int a_cs, int b_rs, int b_cs );
int bl1_is_row_or_col_storage( int rs, int cs );
int bl1_is_row_storage( int rs, int cs );
int bl1_is_col_storage( int rs, int cs );
int bl1_is_gen_storage( int rs, int cs );
int bl1_is_vector( int m, int n );

// --- vector-related ---

int bl1_vector_dim( int m, int n );
int bl1_vector_inc( trans1_t trans, int m, int n, int rs, int cs );

// --- dimension-related ---

int bl1_zero_dim1( int m );
int bl1_zero_dim2( int m, int n );
int bl1_zero_dim3( int m, int k, int n );

