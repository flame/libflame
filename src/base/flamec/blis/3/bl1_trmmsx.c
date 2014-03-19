/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_strmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, float* alpha, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs, float* beta, float* c, int c_rs, int c_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	float*    a_save    = a;
	float*    b_save    = b;
	float*    c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       b_rs_save = b_rs;
	int       b_cs_save = b_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	float     one = bl1_s1();
	float*    b_copy;
	int       dim_a;
	int       b_copy_rs, b_copy_cs;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_set_dim_with_side( side, m, n, &dim_a );
	bl1_screate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_screate_contigm( m,
	                     n,
	                     b_save, b_rs_save, b_cs_save,
	                     &b,     &b_rs,     &b_cs );

	bl1_screate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// Create a copy of B to use in the computation so the original matrix is
	// left untouched.
	b_copy = bl1_sallocm( m, n );

	// Match the strides of B_copy to that of B.
	if ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		b_copy_rs = 1;
		b_copy_cs = m;
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		b_copy_rs = n;
		b_copy_cs = 1;
	}

	// Copy the contents of B to B_copy.
	bl1_scopymt( BLIS1_NO_TRANSPOSE,
	             m,
	             n,
	             b,      b_rs,      b_cs,
	             b_copy, b_copy_rs, b_copy_cs );
	
	// Perform the operation on B_copy.
	bl1_strmm( side,
	           uplo,
	           trans,
	           diag,
	           m,
	           n,
		       alpha,
	           a,      a_rs,      a_cs,
	           b_copy, b_copy_rs, b_copy_cs );

	// Scale C by beta.
	bl1_sscalm( BLIS1_NO_CONJUGATE,
	            m,
	            n,
	            beta,
	            c, c_rs, c_cs );

	// Add B_copy into C.
	bl1_saxpymt( BLIS1_NO_TRANSPOSE,
	             m,
	             n,
	             &one,
	             b_copy, b_copy_rs, b_copy_cs,
	             c,      c_rs,      c_cs );

	// Free the copy of B.
	bl1_sfree( b_copy );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_sfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_sfree_contigm( b_save, b_rs_save, b_cs_save,
	                   &b,     &b_rs,     &b_cs );

	bl1_sfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

void bl1_dtrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, double* alpha, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs, double* beta, double* c, int c_rs, int c_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	double*   a_save    = a;
	double*   b_save    = b;
	double*   c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       b_rs_save = b_rs;
	int       b_cs_save = b_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	double    one = bl1_d1();
	double*   b_copy;
	int       dim_a;
	int       b_copy_rs, b_copy_cs;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_set_dim_with_side( side, m, n, &dim_a );
	bl1_dcreate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_dcreate_contigm( m,
	                     n,
	                     b_save, b_rs_save, b_cs_save,
	                     &b,     &b_rs,     &b_cs );

	bl1_dcreate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// Create a copy of B to use in the computation so the original matrix is
	// left untouched.
	b_copy = bl1_dallocm( m, n );

	// Match the strides of B_copy to that of B.
	if ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		b_copy_rs = 1;
		b_copy_cs = m;
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		b_copy_rs = n;
		b_copy_cs = 1;
	}

	// Copy the contents of B to B_copy.
	bl1_dcopymt( BLIS1_NO_TRANSPOSE,
	             m,
	             n,
	             b,      b_rs,      b_cs,
	             b_copy, b_copy_rs, b_copy_cs );
	
	// Perform the operation on B_copy.
	bl1_dtrmm( side,
	           uplo,
	           trans,
	           diag,
	           m,
	           n,
		       alpha,
	           a,      a_rs,      a_cs,
	           b_copy, b_copy_rs, b_copy_cs );

	// Scale C by beta.
	bl1_dscalm( BLIS1_NO_CONJUGATE,
	            m,
	            n,
	            beta,
	            c, c_rs, c_cs );

	// Add B_copy into C.
	bl1_daxpymt( BLIS1_NO_TRANSPOSE,
	             m,
	             n,
	             &one,
	             b_copy, b_copy_rs, b_copy_cs,
	             c,      c_rs,      c_cs );

	// Free the copy of B.
	bl1_dfree( b_copy );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_dfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_dfree_contigm( b_save, b_rs_save, b_cs_save,
	                   &b,     &b_rs,     &b_cs );

	bl1_dfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

void bl1_ctrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	scomplex* a_save    = a;
	scomplex* b_save    = b;
	scomplex* c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       b_rs_save = b_rs;
	int       b_cs_save = b_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	scomplex  one = bl1_c1();
	scomplex* b_copy;
	int       dim_a;
	int       b_copy_rs, b_copy_cs;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_set_dim_with_side( side, m, n, &dim_a );
	bl1_ccreate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_ccreate_contigm( m,
	                     n,
	                     b_save, b_rs_save, b_cs_save,
	                     &b,     &b_rs,     &b_cs );

	bl1_ccreate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// Create a copy of B to use in the computation so the original matrix is
	// left untouched.
	b_copy = bl1_callocm( m, n );

	// Match the strides of B_copy to that of B.
	if ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		b_copy_rs = 1;
		b_copy_cs = m;
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		b_copy_rs = n;
		b_copy_cs = 1;
	}

	// Copy the contents of B to B_copy.
	bl1_ccopymt( BLIS1_NO_TRANSPOSE,
	             m,
	             n,
	             b,      b_rs,      b_cs,
	             b_copy, b_copy_rs, b_copy_cs );
	
	// Perform the operation on B_copy.
	bl1_ctrmm( side,
	           uplo,
	           trans,
	           diag,
	           m,
	           n,
		       alpha,
	           a,      a_rs,      a_cs,
	           b_copy, b_copy_rs, b_copy_cs );

	// Scale C by beta.
	bl1_cscalm( BLIS1_NO_CONJUGATE,
	            m,
	            n,
	            beta,
	            c, c_rs, c_cs );

	// Add B_copy into C.
	bl1_caxpymt( BLIS1_NO_TRANSPOSE,
	             m,
	             n,
	             &one,
	             b_copy, b_copy_rs, b_copy_cs,
	             c,      c_rs,      c_cs );

	// Free the copy of B.
	bl1_cfree( b_copy );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_cfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_cfree_contigm( b_save, b_rs_save, b_cs_save,
	                   &b,     &b_rs,     &b_cs );

	bl1_cfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

void bl1_ztrmmsx( side1_t side, uplo1_t uplo, trans1_t trans, diag1_t diag, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs )
{
	int       m_save    = m;
	int       n_save    = n;
	dcomplex* a_save    = a;
	dcomplex* b_save    = b;
	dcomplex* c_save    = c;
	int       a_rs_save = a_rs;
	int       a_cs_save = a_cs;
	int       b_rs_save = b_rs;
	int       b_cs_save = b_cs;
	int       c_rs_save = c_rs;
	int       c_cs_save = c_cs;
	dcomplex  one = bl1_z1();
	dcomplex* b_copy;
	int       dim_a;
	int       b_copy_rs, b_copy_cs;

	// Return early if possible.
	if ( bl1_zero_dim2( m, n ) ) return;

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_set_dim_with_side( side, m, n, &dim_a );
	bl1_zcreate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_zcreate_contigm( m,
	                     n,
	                     b_save, b_rs_save, b_cs_save,
	                     &b,     &b_rs,     &b_cs );

	bl1_zcreate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// Create a copy of B to use in the computation so the original matrix is
	// left untouched.
	b_copy = bl1_zallocm( m, n );

	// Match the strides of B_copy to that of B.
	if ( bl1_is_col_storage( b_rs, b_cs ) )
	{
		b_copy_rs = 1;
		b_copy_cs = m;
	}
	else // if ( bl1_is_row_storage( b_rs, b_cs ) )
	{
		b_copy_rs = n;
		b_copy_cs = 1;
	}

	// Copy the contents of B to B_copy.
	bl1_zcopymt( BLIS1_NO_TRANSPOSE,
	             m,
	             n,
	             b,      b_rs,      b_cs,
	             b_copy, b_copy_rs, b_copy_cs );
	
	// Perform the operation on B_copy.
	bl1_ztrmm( side,
	           uplo,
	           trans,
	           diag,
	           m,
	           n,
		       alpha,
	           a,      a_rs,      a_cs,
	           b_copy, b_copy_rs, b_copy_cs );

	// Scale C by beta.
	bl1_zscalm( BLIS1_NO_CONJUGATE,
	            m,
	            n,
	            beta,
	            c, c_rs, c_cs );

	// Add B_copy into C.
	bl1_zaxpymt( BLIS1_NO_TRANSPOSE,
	             m,
	             n,
	             &one,
	             b_copy, b_copy_rs, b_copy_cs,
	             c,      c_rs,      c_cs );

	// Free the copy of B.
	bl1_zfree( b_copy );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_zfree_contigm( a_save, a_rs_save, a_cs_save,
	                   &a,     &a_rs,     &a_cs );

	bl1_zfree_contigm( b_save, b_rs_save, b_cs_save,
	                   &b,     &b_rs,     &b_cs );

	bl1_zfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

