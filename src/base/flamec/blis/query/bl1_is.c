/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

// --- trans -------------------------------------------------------------------

int bl1_is_notrans( trans1_t trans )
{
	return ( trans == BLIS1_NO_TRANSPOSE ); 
}

int bl1_is_trans( trans1_t trans )
{
	return ( trans == BLIS1_TRANSPOSE ); 
}

int bl1_is_conjnotrans( trans1_t trans )
{
	return ( trans == BLIS1_CONJ_NO_TRANSPOSE ); 
}

int bl1_is_conjtrans( trans1_t trans )
{
	return ( trans == BLIS1_CONJ_TRANSPOSE ); 
}

// --- conj --------------------------------------------------------------------

int bl1_is_noconj( conj1_t conj )
{
	return ( conj == BLIS1_NO_CONJUGATE ); 
}

int bl1_is_conj( conj1_t conj )
{
	return ( conj == BLIS1_CONJUGATE ); 
}

// --- uplo --------------------------------------------------------------------

int bl1_is_lower( uplo1_t uplo )
{
	return ( uplo == BLIS1_LOWER_TRIANGULAR ); 
}

int bl1_is_upper( uplo1_t uplo )
{
	return ( uplo == BLIS1_UPPER_TRIANGULAR );
}

// --- side --------------------------------------------------------------------

int bl1_is_left( side1_t side )
{
	return ( side == BLIS1_LEFT ); 
}

int bl1_is_right( side1_t side )
{
	return ( side == BLIS1_RIGHT ); 
}

// --- diag --------------------------------------------------------------------

int bl1_is_nonunit_diag( diag1_t diag )
{
	return ( diag == BLIS1_NONUNIT_DIAG ); 
}

int bl1_is_unit_diag( diag1_t diag )
{
	return ( diag == BLIS1_UNIT_DIAG ); 
}

int bl1_is_zero_diag( diag1_t diag )
{
	return ( diag == BLIS1_ZERO_DIAG ); 
}

// --- storage-related ---------------------------------------------------------

int bl1_is_col_storage( int rs, int cs )
{
	return ( rs == 1 ); 
}

int bl1_is_row_storage( int rs, int cs )
{
	return ( cs == 1 ); 
}

int bl1_is_gen_storage( int rs, int cs )
{
	return ( !bl1_is_col_storage( rs, cs ) && 
             !bl1_is_row_storage( rs, cs ) ); 
}

int bl1_is_vector( int m, int n )
{
	return ( m == 1 || n == 1 ); 
}

// --- dimension-related -------------------------------------------------------

int bl1_zero_dim1( int m )
{
	return ( m == 0 ); 
}

int bl1_zero_dim2( int m, int n )
{
	return ( m == 0 || n == 0 ); 
}

int bl1_zero_dim3( int m, int k, int n )
{
	return ( m == 0 || k == 0 || n == 0 ); 
}

