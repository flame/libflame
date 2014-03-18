
#include "blis1.h"

void bl1_ssymm( side1_t side, uplo1_t uplo, int m, int n, float* alpha, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs, float* beta, float* c, int c_rs, int c_cs )
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
	float     zero = bl1_s0();
	float     one  = bl1_s1();
	float*    b_copy;
	float*    c_trans;
	int       dim_a;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       ldb_copy, incb_copy;
	int       ldc_trans, incc_trans;
	int       symm_needs_copyb  = FALSE;
	int       symm_needs_transb = FALSE;
	int       symm_needs_axpyt  = FALSE;

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

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;
	
	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_c ) * B_c
				// effective operation: C_c += uplo( A_c ) * B_c
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_c ) * B_r
				// effective operation: C_c += uplo( A_c ) * B_c
				symm_needs_copyb = TRUE;
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=  uplo( A_r ) * B_c
				// effective operation: C_c += ~uplo( conj( A_c ) ) * B_c
				bl1_swap_ints( lda, inca );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_r ) * B_r
				// effective operation: C_c += ( B_c * ~uplo( conj( A_c ) ) )^T
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_side( side );
				bl1_toggle_uplo( uplo );

				symm_needs_axpyt = TRUE;
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_c ) * B_c
				// effective operation: C_c += ( uplo( A_c ) * B_c )^T
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );

				symm_needs_axpyt = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_c ) * B_r
				// effective operation: C_c += B_c * ~uplo( conj( A_c ) )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( ldb, incb );

				bl1_swap_ints( m, n );

				bl1_toggle_side( side );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_r ) * B_c
				// effective operation: C_c += B_c^T * ~uplo( A_c )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );

				bl1_swap_ints( m, n );

				bl1_toggle_side( side );
				bl1_toggle_uplo( uplo );

				symm_needs_copyb  = TRUE;
				symm_needs_transb = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_r ) * B_r
				// effective operation: C_c += B_c * conj( ~uplo( A_c ) )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_swap_ints( m, n );

				bl1_toggle_uplo( uplo );
				bl1_toggle_side( side );
			}
		}
	}

	// We need a temporary matrix for the cases where B needs to be copied.
	b_copy    = b;
	ldb_copy  = ldb;
	incb_copy = incb;
	
	// There are two cases where we need to make a copy of B: one where the
	// copy's dimensions are transposed from the original B, and one where
	// the dimensions are not swapped.
	if ( symm_needs_copyb )
	{
		trans1_t transb;

		// Set transb, which determines whether or not we need to copy from B
		// as if it needs a transposition. If a transposition is needed, then
		// m and n and have already been swapped. So in either case m
		// represents the leading dimension of the copy.
		if ( symm_needs_transb ) transb = BLIS1_TRANSPOSE;
		else                     transb = BLIS1_NO_TRANSPOSE;
		
		b_copy    = bl1_sallocm( m, n );
		ldb_copy  = m;
		incb_copy = 1;

		bl1_scopymt( transb,
		             m,
		             n,
		             b,      incb,      ldb,
		             b_copy, incb_copy, ldb_copy );
	}

	// There are two cases where we need to perform the symm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( symm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, and thus C_trans is n-by-m
		// (interpreting both as column-major matrices). So the leading
		// dimension of the temporary matrix holding C^T is n.
		c_trans    = bl1_sallocm( n, m );
		ldc_trans  = n;
		incc_trans = 1;

		// Compute A * B (or B * A) and store the result in C_trans.
		// Note that there is no overlap between the axpyt cases and
		// the conja/copyb cases, hence the use of a, b, lda, and ldb.
		bl1_ssymm_blas( side,
		                uplo,
		                n,
		                m,
		                alpha,
		                a,       lda,
		                b,       ldb,
		                &zero,
		                c_trans, ldc_trans );

		// Scale C by beta.
		bl1_sscalm( BLIS1_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, incc, ldc );
		
		// And finally, accumulate the matrix product in C_trans into C
		// with a transpose.
		bl1_saxpymt( BLIS1_TRANSPOSE,
		             m,
		             n,
		             &one,
		             c_trans, incc_trans, ldc_trans,
		             c,       incc,       ldc );

		// Free the temporary matrix for C.
		bl1_sfree( c_trans );
	}
	else // no extra axpyt step needed
	{
		bl1_ssymm_blas( side,
		                uplo,
		                m,
		                n,
		                alpha,
		                a,      lda,
		                b_copy, ldb_copy,
		                beta,
		                c,      ldc );
	}

	if ( symm_needs_copyb )
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

void bl1_dsymm( side1_t side, uplo1_t uplo, int m, int n, double* alpha, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs, double* beta, double* c, int c_rs, int c_cs )
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
	double    zero = bl1_d0();
	double    one  = bl1_d1();
	double*   b_copy;
	double*   c_trans;
	int       dim_a;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       ldb_copy, incb_copy;
	int       ldc_trans, incc_trans;
	int       symm_needs_copyb  = FALSE;
	int       symm_needs_transb = FALSE;
	int       symm_needs_axpyt  = FALSE;

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

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;
	
	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_c ) * B_c
				// effective operation: C_c += uplo( A_c ) * B_c
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_c ) * B_r
				// effective operation: C_c += uplo( A_c ) * B_c
				symm_needs_copyb = TRUE;
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=  uplo( A_r ) * B_c
				// effective operation: C_c += ~uplo( conj( A_c ) ) * B_c
				bl1_swap_ints( lda, inca );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_r ) * B_r
				// effective operation: C_c += ( B_c * ~uplo( conj( A_c ) ) )^T
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_side( side );
				bl1_toggle_uplo( uplo );

				symm_needs_axpyt = TRUE;
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_c ) * B_c
				// effective operation: C_c += ( uplo( A_c ) * B_c )^T
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );

				symm_needs_axpyt = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_c ) * B_r
				// effective operation: C_c += B_c * ~uplo( conj( A_c ) )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( ldb, incb );

				bl1_swap_ints( m, n );

				bl1_toggle_side( side );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_r ) * B_c
				// effective operation: C_c += B_c^T * ~uplo( A_c )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );

				bl1_swap_ints( m, n );

				bl1_toggle_side( side );
				bl1_toggle_uplo( uplo );

				symm_needs_copyb  = TRUE;
				symm_needs_transb = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_r ) * B_r
				// effective operation: C_c += B_c * conj( ~uplo( A_c ) )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_swap_ints( m, n );

				bl1_toggle_uplo( uplo );
				bl1_toggle_side( side );
			}
		}
	}

	// We need a temporary matrix for the cases where B needs to be copied.
	b_copy    = b;
	ldb_copy  = ldb;
	incb_copy = incb;
	
	// There are two cases where we need to make a copy of B: one where the
	// copy's dimensions are transposed from the original B, and one where
	// the dimensions are not swapped.
	if ( symm_needs_copyb )
	{
		trans1_t transb;

		// Set transb, which determines whether or not we need to copy from B
		// as if it needs a transposition. If a transposition is needed, then
		// m and n and have already been swapped. So in either case m
		// represents the leading dimension of the copy.
		if ( symm_needs_transb ) transb = BLIS1_TRANSPOSE;
		else                     transb = BLIS1_NO_TRANSPOSE;
		
		b_copy    = bl1_dallocm( m, n );
		ldb_copy  = m;
		incb_copy = 1;

		bl1_dcopymt( transb,
		             m,
		             n,
		             b,      incb,      ldb,
		             b_copy, incb_copy, ldb_copy );
	}

	// There are two cases where we need to perform the symm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( symm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, and thus C_trans is n-by-m
		// (interpreting both as column-major matrices). So the leading
		// dimension of the temporary matrix holding C^T is n.
		c_trans    = bl1_dallocm( n, m );
		ldc_trans  = n;
		incc_trans = 1;

		// Compute A * B (or B * A) and store the result in C_trans.
		// Note that there is no overlap between the axpyt cases and
		// the conja/copyb cases, hence the use of a, b, lda, and ldb.
		bl1_dsymm_blas( side,
		                uplo,
		                n,
		                m,
		                alpha,
		                a,       lda,
		                b,       ldb,
		                &zero,
		                c_trans, ldc_trans );

		// Scale C by beta.
		bl1_dscalm( BLIS1_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, incc, ldc );
		
		// And finally, accumulate the matrix product in C_trans into C
		// with a transpose.
		bl1_daxpymt( BLIS1_TRANSPOSE,
		             m,
		             n,
		             &one,
		             c_trans, incc_trans, ldc_trans,
		             c,       incc,       ldc );

		// Free the temporary matrix for C.
		bl1_dfree( c_trans );
	}
	else // no extra axpyt step needed
	{
		bl1_dsymm_blas( side,
		                uplo,
		                m,
		                n,
		                alpha,
		                a,      lda,
		                b_copy, ldb_copy,
		                beta,
		                c,      ldc );
	}

	if ( symm_needs_copyb )
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

void bl1_csymm( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs )
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
	scomplex  zero = bl1_c0();
	scomplex  one  = bl1_c1();
	scomplex* b_copy;
	scomplex* c_trans;
	int       dim_a;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       ldb_copy, incb_copy;
	int       ldc_trans, incc_trans;
	int       symm_needs_copyb  = FALSE;
	int       symm_needs_transb = FALSE;
	int       symm_needs_axpyt  = FALSE;

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

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;
	
	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_c ) * B_c
				// effective operation: C_c += uplo( A_c ) * B_c
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_c ) * B_r
				// effective operation: C_c += uplo( A_c ) * B_c
				symm_needs_copyb = TRUE;
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=  uplo( A_r ) * B_c
				// effective operation: C_c += ~uplo( conj( A_c ) ) * B_c
				bl1_swap_ints( lda, inca );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_r ) * B_r
				// effective operation: C_c += ( B_c * ~uplo( conj( A_c ) ) )^T
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_side( side );
				bl1_toggle_uplo( uplo );

				symm_needs_axpyt = TRUE;
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_c ) * B_c
				// effective operation: C_c += ( uplo( A_c ) * B_c )^T
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );

				symm_needs_axpyt = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_c ) * B_r
				// effective operation: C_c += B_c * ~uplo( conj( A_c ) )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( ldb, incb );

				bl1_swap_ints( m, n );

				bl1_toggle_side( side );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_r ) * B_c
				// effective operation: C_c += B_c^T * ~uplo( A_c )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );

				bl1_swap_ints( m, n );

				bl1_toggle_side( side );
				bl1_toggle_uplo( uplo );

				symm_needs_copyb  = TRUE;
				symm_needs_transb = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_r ) * B_r
				// effective operation: C_c += B_c * conj( ~uplo( A_c ) )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_swap_ints( m, n );

				bl1_toggle_uplo( uplo );
				bl1_toggle_side( side );
			}
		}
	}

	// We need a temporary matrix for the cases where B needs to be copied.
	b_copy    = b;
	ldb_copy  = ldb;
	incb_copy = incb;
	
	// There are two cases where we need to make a copy of B: one where the
	// copy's dimensions are transposed from the original B, and one where
	// the dimensions are not swapped.
	if ( symm_needs_copyb )
	{
		trans1_t transb;

		// Set transb, which determines whether or not we need to copy from B
		// as if it needs a transposition. If a transposition is needed, then
		// m and n and have already been swapped. So in either case m
		// represents the leading dimension of the copy.
		if ( symm_needs_transb ) transb = BLIS1_TRANSPOSE;
		else                     transb = BLIS1_NO_TRANSPOSE;
		
		b_copy    = bl1_callocm( m, n );
		ldb_copy  = m;
		incb_copy = 1;

		bl1_ccopymt( transb,
		             m,
		             n,
		             b,      incb,      ldb,
		             b_copy, incb_copy, ldb_copy );
	}

	// There are two cases where we need to perform the symm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( symm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, and thus C_trans is n-by-m
		// (interpreting both as column-major matrices). So the leading
		// dimension of the temporary matrix holding C^T is n.
		c_trans    = bl1_callocm( n, m );
		ldc_trans  = n;
		incc_trans = 1;

		// Compute A * B (or B * A) and store the result in C_trans.
		// Note that there is no overlap between the axpyt cases and
		// the conja/copyb cases, hence the use of a, b, lda, and ldb.
		bl1_csymm_blas( side,
		                uplo,
		                n,
		                m,
		                alpha,
		                a,       lda,
		                b,       ldb,
		                &zero,
		                c_trans, ldc_trans );

		// Scale C by beta.
		bl1_cscalm( BLIS1_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, incc, ldc );
		
		// And finally, accumulate the matrix product in C_trans into C
		// with a transpose.
		bl1_caxpymt( BLIS1_TRANSPOSE,
		             m,
		             n,
		             &one,
		             c_trans, incc_trans, ldc_trans,
		             c,       incc,       ldc );

		// Free the temporary matrix for C.
		bl1_cfree( c_trans );
	}
	else // no extra axpyt step needed
	{
		bl1_csymm_blas( side,
		                uplo,
		                m,
		                n,
		                alpha,
		                a,      lda,
		                b_copy, ldb_copy,
		                beta,
		                c,      ldc );
	}

	if ( symm_needs_copyb )
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

void bl1_zsymm( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs )
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
	dcomplex  zero = bl1_z0();
	dcomplex  one  = bl1_z1();
	dcomplex* b_copy;
	dcomplex* c_trans;
	int       dim_a;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       ldb_copy, incb_copy;
	int       ldc_trans, incc_trans;
	int       symm_needs_copyb  = FALSE;
	int       symm_needs_transb = FALSE;
	int       symm_needs_axpyt  = FALSE;

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

	// Initialize with values assuming column-major storage.
	lda  = a_cs;
	inca = a_rs;
	ldb  = b_cs;
	incb = b_rs;
	ldc  = c_cs;
	incc = c_rs;
	
	// Adjust the parameters based on the storage of each matrix.
	if ( bl1_is_col_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_c ) * B_c
				// effective operation: C_c += uplo( A_c ) * B_c
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_c ) * B_r
				// effective operation: C_c += uplo( A_c ) * B_c
				symm_needs_copyb = TRUE;
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=  uplo( A_r ) * B_c
				// effective operation: C_c += ~uplo( conj( A_c ) ) * B_c
				bl1_swap_ints( lda, inca );

				bl1_toggle_uplo( uplo );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_r ) * B_r
				// effective operation: C_c += ( B_c * ~uplo( conj( A_c ) ) )^T
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_side( side );
				bl1_toggle_uplo( uplo );

				symm_needs_axpyt = TRUE;
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_c ) * B_c
				// effective operation: C_c += ( uplo( A_c ) * B_c )^T
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );

				symm_needs_axpyt = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_c ) * B_r
				// effective operation: C_c += B_c * ~uplo( conj( A_c ) )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( ldb, incb );

				bl1_swap_ints( m, n );

				bl1_toggle_side( side );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_r ) * B_c
				// effective operation: C_c += B_c^T * ~uplo( A_c )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );

				bl1_swap_ints( m, n );

				bl1_toggle_side( side );
				bl1_toggle_uplo( uplo );

				symm_needs_copyb  = TRUE;
				symm_needs_transb = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_r ) * B_r
				// effective operation: C_c += B_c * conj( ~uplo( A_c ) )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_swap_ints( m, n );

				bl1_toggle_uplo( uplo );
				bl1_toggle_side( side );
			}
		}
	}

	// We need a temporary matrix for the cases where B needs to be copied.
	b_copy    = b;
	ldb_copy  = ldb;
	incb_copy = incb;
	
	// There are two cases where we need to make a copy of B: one where the
	// copy's dimensions are transposed from the original B, and one where
	// the dimensions are not swapped.
	if ( symm_needs_copyb )
	{
		trans1_t transb;

		// Set transb, which determines whether or not we need to copy from B
		// as if it needs a transposition. If a transposition is needed, then
		// m and n and have already been swapped. So in either case m
		// represents the leading dimension of the copy.
		if ( symm_needs_transb ) transb = BLIS1_TRANSPOSE;
		else                     transb = BLIS1_NO_TRANSPOSE;
		
		b_copy    = bl1_zallocm( m, n );
		ldb_copy  = m;
		incb_copy = 1;

		bl1_zcopymt( transb,
		             m,
		             n,
		             b,      incb,      ldb,
		             b_copy, incb_copy, ldb_copy );
	}

	// There are two cases where we need to perform the symm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( symm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, and thus C_trans is n-by-m
		// (interpreting both as column-major matrices). So the leading
		// dimension of the temporary matrix holding C^T is n.
		c_trans    = bl1_zallocm( n, m );
		ldc_trans  = n;
		incc_trans = 1;

		// Compute A * B (or B * A) and store the result in C_trans.
		// Note that there is no overlap between the axpyt cases and
		// the conja/copyb cases, hence the use of a, b, lda, and ldb.
		bl1_zsymm_blas( side,
		                uplo,
		                n,
		                m,
		                alpha,
		                a,       lda,
		                b,       ldb,
		                &zero,
		                c_trans, ldc_trans );

		// Scale C by beta.
		bl1_zscalm( BLIS1_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, incc, ldc );
		
		// And finally, accumulate the matrix product in C_trans into C
		// with a transpose.
		bl1_zaxpymt( BLIS1_TRANSPOSE,
		             m,
		             n,
		             &one,
		             c_trans, incc_trans, ldc_trans,
		             c,       incc,       ldc );

		// Free the temporary matrix for C.
		bl1_zfree( c_trans );
	}
	else // no extra axpyt step needed
	{
		bl1_zsymm_blas( side,
		                uplo,
		                m,
		                n,
		                alpha,
		                a,      lda,
		                b_copy, ldb_copy,
		                beta,
		                c,      ldc );
	}

	if ( symm_needs_copyb )
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

// --- Classic routine wrappers ---

void bl1_ssymm_blas( side1_t side, uplo1_t uplo, int m, int n, float* alpha, float* a, int lda, float* b, int ldb, float* beta, float* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_side( side, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_ssymm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             m,
	             n,
	             *alpha,
	             a, lda,
	             b, ldb,
	             *beta,
	             c, ldc );
#else
	char blas_side;
	char blas_uplo;

	bl1_param_map_to_netlib_side( side, &blas_side );
	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_ssymm( &blas_side,
	           &blas_uplo,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

void bl1_dsymm_blas( side1_t side, uplo1_t uplo, int m, int n, double* alpha, double* a, int lda, double* b, int ldb, double* beta, double* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_side( side, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_dsymm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             m,
	             n,
	             *alpha,
	             a, lda,
	             b, ldb,
	             *beta,
	             c, ldc );
#else
	char blas_side;
	char blas_uplo;

	bl1_param_map_to_netlib_side( side, &blas_side );
	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_dsymm( &blas_side,
	           &blas_uplo,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

void bl1_csymm_blas( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_side( side, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_csymm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             m,
	             n,
	             alpha,
	             a, lda,
	             b, ldb,
	             beta,
	             c, ldc );
#else
	char blas_side;
	char blas_uplo;

	bl1_param_map_to_netlib_side( side, &blas_side );
	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_csymm( &blas_side,
	           &blas_uplo,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

void bl1_zsymm_blas( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_side( side, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_zsymm( cblas_order,
	             cblas_side,
	             cblas_uplo,
	             m,
	             n,
	             alpha,
	             a, lda,
	             b, ldb,
	             beta,
	             c, ldc );
#else
	char blas_side;
	char blas_uplo;

	bl1_param_map_to_netlib_side( side, &blas_side );
	bl1_param_map_to_netlib_uplo( uplo, &blas_uplo );

	F77_zsymm( &blas_side,
	           &blas_uplo,
	           &m,
	           &n,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

