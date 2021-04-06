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

void bl1_check_storage_3m( integer a_rs, integer a_cs, integer b_rs, integer b_cs, integer c_rs, integer c_cs );
void bl1_check_storage_2m( integer a_rs, integer a_cs, integer b_rs, integer b_cs );
int bl1_is_row_or_col_storage( integer rs, integer cs );
int bl1_is_row_storage( integer rs, integer cs );
int bl1_is_col_storage( integer rs, integer cs );
int bl1_is_gen_storage( integer rs, integer cs );
int bl1_is_vector( integer m, integer n );

// --- vector-related ---

int bl1_vector_dim( integer m, integer n );
int bl1_vector_inc( trans1_t trans, integer m, integer n, integer rs, integer cs );

// --- dimension-related ---

int bl1_zero_dim1( integer m );
int bl1_zero_dim2( integer m, integer n );
int bl1_zero_dim3( integer m, integer k, integer n );

