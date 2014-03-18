
#include "blis1.h"

void bl1_sgemm( trans1_t transa, trans1_t transb, int m, int k, int n, float* alpha, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs, float* beta, float* c, int c_rs, int c_cs )
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
	float*    a_unswap;
	float*    b_unswap;
	float*    c_trans;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       ldc_trans, incc_trans;
	int       m_gemm, n_gemm;
	int       gemm_needs_axpyt = FALSE;

	// Return early if possible.
	if ( bl1_zero_dim3( m, k, n ) )
	{
		bl1_sscalm( BLIS1_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, c_rs, c_cs );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_screate_contigmt( transa,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_screate_contigmt( transb,
	                      k,
	                      n,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bl1_screate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// These are used to track the original values of a and b prior to any
	// operand swapping that might take place. This is necessary for proper
	// freeing of memory when one is a temporary contiguous matrix.
	a_unswap = a;
	b_unswap = b;

	// These are used to track the dimensions of the product of the
	// A and B operands to the BLAS invocation of gemm. These differ
	// from m and n when the operands need to be swapped.
	m_gemm = m;
	n_gemm = n;

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
				// requested operation: C_c += tr( A_c ) * tr( B_c )
				// effective operation: C_c += tr( A_c ) * tr( B_c )
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				
				// requested operation: C_c += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( A_c ) * tr( B_c )^T
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( transb );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( A_r )^T * tr( B_c )
				bl1_swap_ints( lda, inca );

				bl1_toggle_trans( transa );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=   tr( A_r ) * tr( B_r )
				// effective operation: C_c += ( tr( B_c ) * tr( A_c ) )^T
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_sswap_pointers( a, b );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );

				gemm_needs_axpyt = TRUE;
				bl1_swap_ints( m_gemm, n_gemm );
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r +=   tr( A_c ) * tr( B_c )
				// effective operation: C_c += ( tr( A_c ) * tr( B_c ) )^T
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );

				gemm_needs_axpyt = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )^T
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( transa );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_sswap_pointers( a, b );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( B_c )^T * tr( A_c )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );

				bl1_toggle_trans( transb );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_sswap_pointers( a, b );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_sswap_pointers( a, b );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
		}
	}

	// There are two cases where we need to perform the gemm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( gemm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, while m_gemm and n_gemm are the
		// dimensions of the actual product op(A)*op(B), which may be n-by-m
		// since the operands may have been swapped.
		c_trans    = bl1_sallocm( m_gemm, n_gemm );
		ldc_trans  = m_gemm;
		incc_trans = 1;

		// Compute tr( A ) * tr( B ), where A and B may have been swapped
		// to reference the other, and store the result in C_trans.
		bl1_sgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
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
		bl1_sgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a, lda,
		                b, ldb,
		                beta,
		                c, ldc );
	}

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_sfree_contigm( a_save,    a_rs_save, a_cs_save,
	                   &a_unswap, &a_rs,     &a_cs );

	bl1_sfree_contigm( b_save,    b_rs_save, b_cs_save,
	                   &b_unswap, &b_rs,     &b_cs );

	bl1_sfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

void bl1_dgemm( trans1_t transa, trans1_t transb, int m, int k, int n, double* alpha, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs, double* beta, double* c, int c_rs, int c_cs )
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
	double*   a_unswap;
	double*   b_unswap;
	double*   c_trans;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       ldc_trans, incc_trans;
	int       m_gemm, n_gemm;
	int       gemm_needs_axpyt = FALSE;

	// Return early if possible.
	if ( bl1_zero_dim3( m, k, n ) )
	{
		bl1_dscalm( BLIS1_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, c_rs, c_cs );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_dcreate_contigmt( transa,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_dcreate_contigmt( transb,
	                      k,
	                      n,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bl1_dcreate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// These are used to track the original values of a and b prior to any
	// operand swapping that might take place. This is necessary for proper
	// freeing of memory when one is a temporary contiguous matrix.
	a_unswap = a;
	b_unswap = b;

	// These are used to track the dimensions of the product of the
	// A and B operands to the BLAS invocation of gemm. These differ
	// from m and n when the operands need to be swapped.
	m_gemm = m;
	n_gemm = n;

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
				// requested operation: C_c += tr( A_c ) * tr( B_c )
				// effective operation: C_c += tr( A_c ) * tr( B_c )
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				
				// requested operation: C_c += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( A_c ) * tr( B_c )^T
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( transb );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( A_r )^T * tr( B_c )
				bl1_swap_ints( lda, inca );

				bl1_toggle_trans( transa );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=   tr( A_r ) * tr( B_r )
				// effective operation: C_c += ( tr( B_c ) * tr( A_c ) )^T
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_dswap_pointers( a, b );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );

				gemm_needs_axpyt = TRUE;
				bl1_swap_ints( m_gemm, n_gemm );
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r +=   tr( A_c ) * tr( B_c )
				// effective operation: C_c += ( tr( A_c ) * tr( B_c ) )^T
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );

				gemm_needs_axpyt = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )^T
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( transa );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_dswap_pointers( a, b );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( B_c )^T * tr( A_c )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );

				bl1_toggle_trans( transb );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_dswap_pointers( a, b );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_dswap_pointers( a, b );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
		}
	}

	// There are two cases where we need to perform the gemm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( gemm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, while m_gemm and n_gemm are the
		// dimensions of the actual product op(A)*op(B), which may be n-by-m
		// since the operands may have been swapped.
		c_trans    = bl1_dallocm( m_gemm, n_gemm );
		ldc_trans  = m_gemm;
		incc_trans = 1;

		// Compute tr( A ) * tr( B ), where A and B may have been swapped
		// to reference the other, and store the result in C_trans.
		bl1_dgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
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
		bl1_dgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a, lda,
		                b, ldb,
		                beta,
		                c, ldc );
	}

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_dfree_contigm( a_save,    a_rs_save, a_cs_save,
	                   &a_unswap, &a_rs,     &a_cs );

	bl1_dfree_contigm( b_save,    b_rs_save, b_cs_save,
	                   &b_unswap, &b_rs,     &b_cs );

	bl1_dfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

void bl1_cgemm( trans1_t transa, trans1_t transb, int m, int k, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs )
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
	scomplex* a_unswap;
	scomplex* b_unswap;
	scomplex* a_conj;
	scomplex* b_conj;
	scomplex* c_trans;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       lda_conj, inca_conj;
	int       ldb_conj, incb_conj;
	int       ldc_trans, incc_trans;
	int       m_gemm, n_gemm;
	int       gemm_needs_axpyt = FALSE;
	int       a_was_copied;
	int       b_was_copied;

	// Return early if possible.
	if ( bl1_zero_dim3( m, k, n ) )
	{
		bl1_cscalm( BLIS1_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, c_rs, c_cs );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_ccreate_contigmt( transa,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_ccreate_contigmt( transb,
	                      k,
	                      n,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bl1_ccreate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// Figure out whether A and/or B was copied to contiguous memory. This
	// is used later to prevent redundant copying.
	a_was_copied = ( a != a_save );
	b_was_copied = ( b != b_save );

	// These are used to track the original values of a and b prior to any
	// operand swapping that might take place. This is necessary for proper
	// freeing of memory when one is a temporary contiguous matrix.
	a_unswap = a;
	b_unswap = b;

	// These are used to track the dimensions of the product of the
	// A and B operands to the BLAS invocation of gemm. These differ
	// from m and n when the operands need to be swapped.
	m_gemm = m;
	n_gemm = n;

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
				// requested operation: C_c += tr( A_c ) * tr( B_c )
				// effective operation: C_c += tr( A_c ) * tr( B_c )
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				
				// requested operation: C_c += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( A_c ) * tr( B_c )^T
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( transb );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( A_r )^T * tr( B_c )
				bl1_swap_ints( lda, inca );

				bl1_toggle_trans( transa );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=   tr( A_r ) * tr( B_r )
				// effective operation: C_c += ( tr( B_c ) * tr( A_c ) )^T
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_cswap_pointers( a, b );
				bl1_swap_ints( a_was_copied, b_was_copied );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );

				gemm_needs_axpyt = TRUE;
				bl1_swap_ints( m_gemm, n_gemm );
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r +=   tr( A_c ) * tr( B_c )
				// effective operation: C_c += ( tr( A_c ) * tr( B_c ) )^T
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );

				gemm_needs_axpyt = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )^T
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( transa );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_cswap_pointers( a, b );
				bl1_swap_ints( a_was_copied, b_was_copied );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( B_c )^T * tr( A_c )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );

				bl1_toggle_trans( transb );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_cswap_pointers( a, b );
				bl1_swap_ints( a_was_copied, b_was_copied );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_cswap_pointers( a, b );
				bl1_swap_ints( a_was_copied, b_was_copied );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
		}
	}

	// We need a temporary matrix for the case where A is conjugated.
	a_conj    = a;
	lda_conj  = lda;
	inca_conj = inca;

	// If transa indicates conjugate-no-transpose and A was not already
	// copied, then copy and conjugate it to a temporary matrix. Otherwise,
	// if transa indicates conjugate-no-transpose and A was already copied,
	// just conjugate it.
	if ( bl1_is_conjnotrans( transa ) && !a_was_copied )
	{
		a_conj    = bl1_callocm( m_gemm, k );
		lda_conj  = m_gemm;
		inca_conj = 1;

		bl1_ccopymt( BLIS1_CONJ_NO_TRANSPOSE,
		             m_gemm,
		             k,
		             a,      inca,      lda,
		             a_conj, inca_conj, lda_conj );
	}
	else if ( bl1_is_conjnotrans( transa ) && a_was_copied )
	{
		bl1_cconjm( m_gemm,
		            k,
		            a_conj, inca_conj, lda_conj );
	}

	// We need a temporary matrix for the case where B is conjugated.
	b_conj    = b;
	ldb_conj  = ldb;
	incb_conj = incb;

	// If transb indicates conjugate-no-transpose and B was not already
	// copied, then copy and conjugate it to a temporary matrix. Otherwise,
	// if transb indicates conjugate-no-transpose and B was already copied,
	// just conjugate it.
	if ( bl1_is_conjnotrans( transb ) && !b_was_copied )
	{
		b_conj    = bl1_callocm( k, n_gemm );
		ldb_conj  = k;
		incb_conj = 1;

		bl1_ccopymt( BLIS1_CONJ_NO_TRANSPOSE,
		             k,
		             n_gemm,
		             b,      incb,      ldb,
		             b_conj, incb_conj, ldb_conj );
	}
	else if ( bl1_is_conjnotrans( transb ) && b_was_copied )
	{
		bl1_cconjm( k,
		            n_gemm,
		            b_conj, incb_conj, ldb_conj );
	}

	// There are two cases where we need to perform the gemm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( gemm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, while m_gemm and n_gemm are the
		// dimensions of the actual product op(A)*op(B), which may be n-by-m
		// since the operands may have been swapped.
		c_trans    = bl1_callocm( m_gemm, n_gemm );
		ldc_trans  = m_gemm;
		incc_trans = 1;

		// Compute tr( A ) * tr( B ), where A and B may have been swapped
		// to reference the other, and store the result in C_trans.
		bl1_cgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a_conj,  lda_conj,
		                b_conj,  ldb_conj,
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
		bl1_cgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a_conj, lda_conj,
		                b_conj, ldb_conj,
		                beta,
		                c,      ldc );
	}

	if ( bl1_is_conjnotrans( transa ) && !a_was_copied )
		bl1_cfree( a_conj );

	if ( bl1_is_conjnotrans( transb ) && !b_was_copied )
		bl1_cfree( b_conj );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_cfree_contigm( a_save,    a_rs_save, a_cs_save,
	                   &a_unswap, &a_rs,     &a_cs );

	bl1_cfree_contigm( b_save,    b_rs_save, b_cs_save,
	                   &b_unswap, &b_rs,     &b_cs );

	bl1_cfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

void bl1_zgemm( trans1_t transa, trans1_t transb, int m, int k, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs )
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
	dcomplex* a_unswap;
	dcomplex* b_unswap;
	dcomplex* a_conj;
	dcomplex* b_conj;
	dcomplex* c_trans;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       lda_conj, inca_conj;
	int       ldb_conj, incb_conj;
	int       ldc_trans, incc_trans;
	int       m_gemm, n_gemm;
	int       gemm_needs_axpyt = FALSE;
	int       a_was_copied;
	int       b_was_copied;

	// Return early if possible.
	if ( bl1_zero_dim3( m, k, n ) )
	{
		bl1_zscalm( BLIS1_NO_CONJUGATE,
		            m,
		            n,
		            beta,
		            c, c_rs, c_cs );
		return;
	}

	// If necessary, allocate, initialize, and use a temporary contiguous
	// copy of each matrix rather than the original matrices.
	bl1_zcreate_contigmt( transa,
	                      m,
	                      k,
	                      a_save, a_rs_save, a_cs_save,
	                      &a,     &a_rs,     &a_cs );

	bl1_zcreate_contigmt( transb,
	                      k,
	                      n,
	                      b_save, b_rs_save, b_cs_save,
	                      &b,     &b_rs,     &b_cs );

	bl1_zcreate_contigm( m,
	                     n,
	                     c_save, c_rs_save, c_cs_save,
	                     &c,     &c_rs,     &c_cs );

	// Figure out whether A and/or B was copied to contiguous memory. This
	// is used later to prevent redundant copying.
	a_was_copied = ( a != a_save );
	b_was_copied = ( b != b_save );

	// These are used to track the original values of a and b prior to any
	// operand swapping that might take place. This is necessary for proper
	// freeing of memory when one is a temporary contiguous matrix.
	a_unswap = a;
	b_unswap = b;

	// These are used to track the dimensions of the product of the
	// A and B operands to the BLAS invocation of gemm. These differ
	// from m and n when the operands need to be swapped.
	m_gemm = m;
	n_gemm = n;

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
				// requested operation: C_c += tr( A_c ) * tr( B_c )
				// effective operation: C_c += tr( A_c ) * tr( B_c )
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				
				// requested operation: C_c += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( A_c ) * tr( B_c )^T
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( transb );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( A_r )^T * tr( B_c )
				bl1_swap_ints( lda, inca );

				bl1_toggle_trans( transa );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c +=   tr( A_r ) * tr( B_r )
				// effective operation: C_c += ( tr( B_c ) * tr( A_c ) )^T
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_zswap_pointers( a, b );
				bl1_swap_ints( a_was_copied, b_was_copied );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );

				gemm_needs_axpyt = TRUE;
				bl1_swap_ints( m_gemm, n_gemm );
			}
		}
	}
	else // if ( bl1_is_row_storage( c_rs, c_cs ) )
	{
		if ( bl1_is_col_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r +=   tr( A_c ) * tr( B_c )
				// effective operation: C_c += ( tr( A_c ) * tr( B_c ) )^T
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );

				gemm_needs_axpyt = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_c ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )^T
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_trans( transa );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_zswap_pointers( a, b );
				bl1_swap_ints( a_was_copied, b_was_copied );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
		}
		else // if ( bl1_is_row_storage( a_rs, a_cs ) )
		{
			if ( bl1_is_col_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r )   * tr( B_c )
				// effective operation: C_c += tr( B_c )^T * tr( A_c )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( lda, inca );

				bl1_toggle_trans( transb );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_zswap_pointers( a, b );
				bl1_swap_ints( a_was_copied, b_was_copied );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += tr( A_r ) * tr( B_r )
				// effective operation: C_c += tr( B_c ) * tr( A_c )
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );
				bl1_swap_ints( ldc, incc );

				bl1_swap_ints( m, n );
				bl1_swap_ints( m_gemm, n_gemm );
				bl1_zswap_pointers( a, b );
				bl1_swap_ints( a_was_copied, b_was_copied );
				bl1_swap_ints( lda, ldb );
				bl1_swap_ints( inca, incb );
				bl1_swap_trans( transa, transb );
			}
		}
	}

	// We need a temporary matrix for the case where A is conjugated.
	a_conj    = a;
	lda_conj  = lda;
	inca_conj = inca;

	// If transa indicates conjugate-no-transpose and A was not already
	// copied, then copy and conjugate it to a temporary matrix. Otherwise,
	// if transa indicates conjugate-no-transpose and A was already copied,
	// just conjugate it.
	if ( bl1_is_conjnotrans( transa ) && !a_was_copied )
	{
		a_conj    = bl1_zallocm( m_gemm, k );
		lda_conj  = m_gemm;
		inca_conj = 1;

		bl1_zcopymt( BLIS1_CONJ_NO_TRANSPOSE,
		             m_gemm,
		             k,
		             a,      inca,      lda,
		             a_conj, inca_conj, lda_conj );
	}
	else if ( bl1_is_conjnotrans( transa ) && a_was_copied )
	{
		bl1_zconjm( m_gemm,
		            k,
		            a_conj, inca_conj, lda_conj );
	}

	// We need a temporary matrix for the case where B is conjugated.
	b_conj    = b;
	ldb_conj  = ldb;
	incb_conj = incb;

	// If transb indicates conjugate-no-transpose and B was not already
	// copied, then copy and conjugate it to a temporary matrix. Otherwise,
	// if transb indicates conjugate-no-transpose and B was already copied,
	// just conjugate it.
	if ( bl1_is_conjnotrans( transb ) && !b_was_copied )
	{
		b_conj    = bl1_zallocm( k, n_gemm );
		ldb_conj  = k;
		incb_conj = 1;

		bl1_zcopymt( BLIS1_CONJ_NO_TRANSPOSE,
		             k,
		             n_gemm,
		             b,      incb,      ldb,
		             b_conj, incb_conj, ldb_conj );
	}
	else if ( bl1_is_conjnotrans( transb ) && b_was_copied )
	{
		bl1_zconjm( k,
		            n_gemm,
		            b_conj, incb_conj, ldb_conj );
	}

	// There are two cases where we need to perform the gemm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( gemm_needs_axpyt )
	{
		// We need a temporary matrix for holding C^T. Notice that m and n
		// represent the dimensions of C, while m_gemm and n_gemm are the
		// dimensions of the actual product op(A)*op(B), which may be n-by-m
		// since the operands may have been swapped.
		c_trans    = bl1_zallocm( m_gemm, n_gemm );
		ldc_trans  = m_gemm;
		incc_trans = 1;

		// Compute tr( A ) * tr( B ), where A and B may have been swapped
		// to reference the other, and store the result in C_trans.
		bl1_zgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a_conj,  lda_conj,
		                b_conj,  ldb_conj,
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
		bl1_zgemm_blas( transa,
		                transb,
		                m_gemm,
		                n_gemm,
		                k,
		                alpha,
		                a_conj, lda_conj,
		                b_conj, ldb_conj,
		                beta,
		                c,      ldc );
	}

	if ( bl1_is_conjnotrans( transa ) && !a_was_copied )
		bl1_zfree( a_conj );

	if ( bl1_is_conjnotrans( transb ) && !b_was_copied )
		bl1_zfree( b_conj );

	// Free any temporary contiguous matrices, copying the result back to
	// the original matrix.
	bl1_zfree_contigm( a_save,    a_rs_save, a_cs_save,
	                   &a_unswap, &a_rs,     &a_cs );

	bl1_zfree_contigm( b_save,    b_rs_save, b_cs_save,
	                   &b_unswap, &b_rs,     &b_cs );

	bl1_zfree_saved_contigm( m_save,
	                         n_save,
	                         c_save, c_rs_save, c_cs_save,
	                         &c,     &c_rs,     &c_cs );
}

// --- Classic routine wrappers ---

void bl1_sgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, float* alpha, float* a, int lda, float* b, int ldb, float* beta, float* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;
	enum CBLAS_TRANSPOSE cblas_transb;

	bl1_param_map_to_netlib_trans( transa, &cblas_transa );
	bl1_param_map_to_netlib_trans( transb, &cblas_transb );

	cblas_sgemm( cblas_order,
	             cblas_transa,
	             cblas_transb,
	             m,
	             n,
	             k,
	             *alpha,
	             a, lda,
	             b, ldb,
	             *beta,
	             c, ldc );
#else
	char blas_transa;
	char blas_transb;

	bl1_param_map_to_netlib_trans( transa, &blas_transa );
	bl1_param_map_to_netlib_trans( transb, &blas_transb );

	F77_sgemm( &blas_transa,
	           &blas_transb,
	           &m,
	           &n,
	           &k,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

void bl1_dgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, double* alpha, double* a, int lda, double* b, int ldb, double* beta, double* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;
	enum CBLAS_TRANSPOSE cblas_transb;

	bl1_param_map_to_netlib_trans( transa, &cblas_transa );
	bl1_param_map_to_netlib_trans( transb, &cblas_transb );

	cblas_dgemm( cblas_order,
	             cblas_transa,
	             cblas_transb,
	             m,
	             n,
	             k,
	             *alpha,
	             a, lda,
	             b, ldb,
	             *beta,
	             c, ldc );
#else
	char blas_transa;
	char blas_transb;

	bl1_param_map_to_netlib_trans( transa, &blas_transa );
	bl1_param_map_to_netlib_trans( transb, &blas_transb );

	F77_dgemm( &blas_transa,
	           &blas_transb,
	           &m,
	           &n,
	           &k,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

void bl1_cgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;
	enum CBLAS_TRANSPOSE cblas_transb;

	bl1_param_map_to_netlib_trans( transa, &cblas_transa );
	bl1_param_map_to_netlib_trans( transb, &cblas_transb );

	cblas_cgemm( cblas_order,
	             cblas_transa,
	             cblas_transb,
	             m,
	             n,
	             k,
	             alpha,
	             a, lda,
	             b, ldb,
	             beta,
	             c, ldc );
#else
	char blas_transa;
	char blas_transb;

	bl1_param_map_to_netlib_trans( transa, &blas_transa );
	bl1_param_map_to_netlib_trans( transb, &blas_transb );

	F77_cgemm( &blas_transa,
	           &blas_transb,
	           &m,
	           &n,
	           &k,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

void bl1_zgemm_blas( trans1_t transa, trans1_t transb, int m, int n, int k, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER     cblas_order = CblasColMajor;
	enum CBLAS_TRANSPOSE cblas_transa;
	enum CBLAS_TRANSPOSE cblas_transb;

	bl1_param_map_to_netlib_trans( transa, &cblas_transa );
	bl1_param_map_to_netlib_trans( transb, &cblas_transb );

	cblas_zgemm( cblas_order,
	             cblas_transa,
	             cblas_transb,
	             m,
	             n,
	             k,
	             alpha,
	             a, lda,
	             b, ldb,
	             beta,
	             c, ldc );
#else
	char blas_transa;
	char blas_transb;

	bl1_param_map_to_netlib_trans( transa, &blas_transa );
	bl1_param_map_to_netlib_trans( transb, &blas_transb );

	F77_zgemm( &blas_transa,
	           &blas_transb,
	           &m,
	           &n,
	           &k,
	           alpha,
	           a, &lda,
	           b, &ldb,
	           beta,
	           c, &ldc );
#endif
}

