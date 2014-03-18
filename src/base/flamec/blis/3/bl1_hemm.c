
#include "blis1.h"

void bl1_shemm( side1_t side, uplo1_t uplo, int m, int n, float* alpha, float* a, int a_rs, int a_cs, float* b, int b_rs, int b_cs, float* beta, float* c, int c_rs, int c_cs )
{
	bl1_ssymm( side,
	           uplo,
	           m,
	           n,
	           alpha,
	           a, a_rs, a_cs,
	           b, b_rs, b_cs,
	           beta,
	           c, c_rs, c_cs );
}
void bl1_dhemm( side1_t side, uplo1_t uplo, int m, int n, double* alpha, double* a, int a_rs, int a_cs, double* b, int b_rs, int b_cs, double* beta, double* c, int c_rs, int c_cs )
{
	bl1_dsymm( side,
	           uplo,
	           m,
	           n,
	           alpha,
	           a, a_rs, a_cs,
	           b, b_rs, b_cs,
	           beta,
	           c, c_rs, c_cs );
}

void bl1_chemm( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int a_rs, int a_cs, scomplex* b, int b_rs, int b_cs, scomplex* beta, scomplex* c, int c_rs, int c_cs )
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
	scomplex* a_conj;
	scomplex* b_copy;
	scomplex* c_trans;
	int       dim_a;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       lda_conj, inca_conj;
	int       ldb_copy, incb_copy;
	int       ldc_trans, incc_trans;
	int       hemm_needs_conja  = FALSE;
	int       hemm_needs_copyb  = FALSE;
	int       hemm_needs_transb = FALSE;
	int       hemm_needs_axpyt  = FALSE;
	int       a_was_copied;

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

	// Figure out whether A was copied to contiguous memory. This is used to
	// prevent redundant copying.
	a_was_copied = ( a != a_save );

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
				hemm_needs_copyb = TRUE;
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

				hemm_needs_conja = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_r ) * B_r
				// effective operation: C_c += ( B_c * ~uplo( conj( A_c ) ) )^T
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_side( side );
				bl1_toggle_uplo( uplo );

				hemm_needs_axpyt = TRUE;
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

				hemm_needs_axpyt = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_c ) * B_r
				// effective operation: C_c += B_c * ~uplo( conj( A_c ) )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( ldb, incb );

				bl1_swap_ints( m, n );

				bl1_toggle_side( side );

				hemm_needs_conja = TRUE;
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

				hemm_needs_copyb  = TRUE;
				hemm_needs_transb = TRUE;
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

	// We need a temporary matrix for the cases where A is conjugated.
	a_conj    = a;
	lda_conj  = lda;
	inca_conj = inca;
	
	if ( hemm_needs_conja && !a_was_copied )
	{
		int dim_a;

		bl1_set_dim_with_side( side, m, n, &dim_a );
		
		a_conj    = bl1_callocm( dim_a, dim_a );
		lda_conj  = dim_a;
		inca_conj = 1;

		bl1_ccopymrt( uplo,
		              BLIS1_CONJ_NO_TRANSPOSE,
		              dim_a,
		              dim_a,
		              a,      inca,      lda,
		              a_conj, inca_conj, lda_conj );
	}
	else if ( hemm_needs_conja && a_was_copied )
	{
		int dim_a;

		bl1_set_dim_with_side( side, m, n, &dim_a );

		bl1_cconjmr( uplo,
		             dim_a,
		             dim_a,
		             a_conj, inca_conj, lda_conj );
	}
	
	// We need a temporary matrix for the cases where B needs to be copied.
	b_copy    = b;
	ldb_copy  = ldb;
	incb_copy = incb;
	
	// There are two cases where we need to make a copy of B: one where the
	// copy's dimensions are transposed from the original B, and one where
	// the dimensions are not swapped.
	if ( hemm_needs_copyb )
	{
		trans1_t transb;

		// Set transb, which determines whether or not we need to copy from B
		// as if it needs a transposition. If a transposition is needed, then
		// m and n and have already been swapped. So in either case m
		// represents the leading dimension of the copy.
		if ( hemm_needs_transb ) transb = BLIS1_TRANSPOSE;
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

	// There are two cases where we need to perform the hemm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( hemm_needs_axpyt )
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
		bl1_chemm_blas( side,
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
		bl1_chemm_blas( side,
		                uplo,
		                m,
		                n,
		                alpha,
		                a_conj, lda_conj,
		                b_copy, ldb_copy,
		                beta,
		                c,      ldc );
	}

	if ( hemm_needs_conja && !a_was_copied )
		bl1_cfree( a_conj );

	if ( hemm_needs_copyb )
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

void bl1_zhemm( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int a_rs, int a_cs, dcomplex* b, int b_rs, int b_cs, dcomplex* beta, dcomplex* c, int c_rs, int c_cs )
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
	dcomplex* a_conj;
	dcomplex* b_copy;
	dcomplex* c_trans;
	int       dim_a;
	int       lda, inca;
	int       ldb, incb;
	int       ldc, incc;
	int       lda_conj, inca_conj;
	int       ldb_copy, incb_copy;
	int       ldc_trans, incc_trans;
	int       hemm_needs_conja  = FALSE;
	int       hemm_needs_copyb  = FALSE;
	int       hemm_needs_transb = FALSE;
	int       hemm_needs_axpyt  = FALSE;
	int       a_was_copied;

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

	// Figure out whether A was copied to contiguous memory. This is used to
	// prevent redundant copying.
	a_was_copied = ( a != a_save );

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
				hemm_needs_copyb = TRUE;
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

				hemm_needs_conja = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_c += uplo( A_r ) * B_r
				// effective operation: C_c += ( B_c * ~uplo( conj( A_c ) ) )^T
				bl1_swap_ints( lda, inca );
				bl1_swap_ints( ldb, incb );

				bl1_toggle_side( side );
				bl1_toggle_uplo( uplo );

				hemm_needs_axpyt = TRUE;
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

				hemm_needs_axpyt = TRUE;
			}
			else // if ( bl1_is_row_storage( b_rs, b_cs ) )
			{
				// requested operation: C_r += uplo( A_c ) * B_r
				// effective operation: C_c += B_c * ~uplo( conj( A_c ) )
				bl1_swap_ints( ldc, incc );
				bl1_swap_ints( ldb, incb );

				bl1_swap_ints( m, n );

				bl1_toggle_side( side );

				hemm_needs_conja = TRUE;
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

				hemm_needs_copyb  = TRUE;
				hemm_needs_transb = TRUE;
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

	// We need a temporary matrix for the cases where A is conjugated.
	a_conj    = a;
	lda_conj  = lda;
	inca_conj = inca;
	
	if ( hemm_needs_conja && !a_was_copied )
	{
		int dim_a;

		bl1_set_dim_with_side( side, m, n, &dim_a );
		
		a_conj    = bl1_zallocm( dim_a, dim_a );
		lda_conj  = dim_a;
		inca_conj = 1;

		bl1_zcopymrt( uplo,
		              BLIS1_CONJ_NO_TRANSPOSE,
		              dim_a,
		              dim_a,
		              a,      inca,      lda,
		              a_conj, inca_conj, lda_conj );
	}
	else if ( hemm_needs_conja && a_was_copied )
	{
		int dim_a;

		bl1_set_dim_with_side( side, m, n, &dim_a );
		
		bl1_zconjmr( uplo,
		             dim_a,
		             dim_a,
		             a_conj, inca_conj, lda_conj );
	}
	
	// We need a temporary matrix for the cases where B needs to be copied.
	b_copy    = b;
	ldb_copy  = ldb;
	incb_copy = incb;
	
	// There are two cases where we need to make a copy of B: one where the
	// copy's dimensions are transposed from the original B, and one where
	// the dimensions are not swapped.
	if ( hemm_needs_copyb )
	{
		trans1_t transb;

		// Set transb, which determines whether or not we need to copy from B
		// as if it needs a transposition. If a transposition is needed, then
		// m and n and have already been swapped. So in either case m
		// represents the leading dimension of the copy.
		if ( hemm_needs_transb ) transb = BLIS1_TRANSPOSE;
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

	// There are two cases where we need to perform the hemm and then axpy
	// the result into C with a transposition. We handle those cases here.
	if ( hemm_needs_axpyt )
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
		bl1_zhemm_blas( side,
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
		bl1_zhemm_blas( side,
		                uplo,
		                m,
		                n,
		                alpha,
		                a_conj, lda_conj,
		                b_copy, ldb_copy,
		                beta,
		                c,      ldc );
	}

	if ( hemm_needs_conja && !a_was_copied )
		bl1_zfree( a_conj );

	if ( hemm_needs_copyb )
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

void bl1_chemm_blas( side1_t side, uplo1_t uplo, int m, int n, scomplex* alpha, scomplex* a, int lda, scomplex* b, int ldb, scomplex* beta, scomplex* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_side( side, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_chemm( cblas_order,
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

	F77_chemm( &blas_side,
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

void bl1_zhemm_blas( side1_t side, uplo1_t uplo, int m, int n, dcomplex* alpha, dcomplex* a, int lda, dcomplex* b, int ldb, dcomplex* beta, dcomplex* c, int ldc )
{
#ifdef BLIS1_ENABLE_CBLAS_INTERFACES
	enum CBLAS_ORDER cblas_order = CblasColMajor;
	enum CBLAS_SIDE  cblas_side;
	enum CBLAS_UPLO  cblas_uplo;

	bl1_param_map_to_netlib_side( side, &cblas_side );
	bl1_param_map_to_netlib_uplo( uplo, &cblas_uplo );

	cblas_zhemm( cblas_order,
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

	F77_zhemm( &blas_side,
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

