
#include "FLAME.h"

FLA_Error FLASH_Part_create_2x1( FLA_Obj A,    FLA_Obj* AT,
                                               FLA_Obj* AB,
                                 dim_t n_rows, FLA_Side side )
{
	FLA_Datatype dt_A;
	dim_t        m_A,  n_A;
	dim_t        m_A_base, n_A_base;
	dim_t        m_AT, n_AT;
	dim_t        m_AB, n_AB;
	dim_t        depth;
	dim_t*       b_m;
	dim_t*       b_n;
	dim_t        offm_A,  offn_A;
	dim_t        offm_AT, offn_AT;
	dim_t        offm_AB, offn_AB;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Part_2x1_check( A,    AT,
		                          AB,     n_rows, side );

	// Safeguard: if n_rows > m, reduce n_rows to m.
	if ( n_rows > FLASH_Obj_scalar_length( A ) )
		n_rows = FLASH_Obj_scalar_length( A );

	// Acquire various properties of the hierarchical matrix view.
	dt_A     = FLASH_Obj_datatype( A );
	m_A      = FLASH_Obj_scalar_length( A );
	n_A      = FLASH_Obj_scalar_width( A );
	offm_A   = FLASH_Obj_scalar_row_offset( A );
	offn_A   = FLASH_Obj_scalar_col_offset( A );
	m_A_base = FLASH_Obj_base_scalar_length( A );
	n_A_base = FLASH_Obj_base_scalar_width( A );
	depth    = FLASH_Obj_depth( A );

	// Allocate a pair of temporary arrays for the blocksizes, whose lengths
	// are equal to the object's hierarchical depth.
	b_m = ( dim_t* ) FLA_malloc( depth * sizeof( dim_t ) );
	b_n = ( dim_t* ) FLA_malloc( depth * sizeof( dim_t ) );

	// Accumulate the blocksizes into the blocksize buffers.
	FLASH_Obj_blocksizes( A, b_m, b_n );

	// Adjust n_rows to be (m - n_rows) if the side specified is on the
	// bottom so that the right values get assigned below.
	if ( side == FLA_BOTTOM ) n_rows = m_A - n_rows;

	// Set the dimensions of the partitions.
	m_AT = n_rows;
	n_AT = n_A;
	m_AB = m_A - n_rows;
	n_AB = n_A;

	// Set the offsets.
	offm_AT = offm_A + 0;
	offn_AT = offn_A + 0;
	offm_AB = offm_A + m_AT;
	offn_AB = offn_A + 0;
	
	// Create bufferless hierarhical objects that have the desired dimensions
	// for the views.
	FLASH_Obj_create_without_buffer_ext( dt_A, m_A_base, n_A_base, depth, b_m, b_n, AT );
	FLASH_Obj_create_without_buffer_ext( dt_A, m_A_base, n_A_base, depth, b_m, b_n, AB );

/*
printf( "depth      %d\n", depth );
printf( "b_m/n[0]   %d %d\n", b_m[0], b_n[0] );
printf( "b_m/n_tl   %d %d\n", FLASH_Obj_scalar_length_tl( A ), FLASH_Obj_scalar_width_tl( A ) );
printf( "m/n_A_base %d %d\n", m_A_base, n_A_base );
printf( "offm/n_AT: %d %d\n", offm_AT, offn_AT );
printf( "m/n_AT:    %d %d\n", m_AT, n_AT );
printf( "offm/n_AB: %d %d\n", offm_AB, offn_AB );
printf( "m/n_AB:    %d %d\n", m_AB, n_AB );
printf( "A is       %d %d\n", FLA_Obj_length( A ), FLA_Obj_width( A ) );
printf( "AT is      %d %d\n", FLA_Obj_length( *AT ), FLA_Obj_width( *AT ) );
printf( "AB is      %d %d\n", FLA_Obj_length( *AB ), FLA_Obj_width( *AB ) );
*/

	// Recursively walk the hierarchy and adjust the views so that they
	// collectively refer to the absolute offsets given, and attach the
	// leaf-level numerical buffers of A to the new views.
	FLASH_Obj_adjust_views( TRUE, offm_AT, offn_AT, m_AT, n_AT, A, AT );
	FLASH_Obj_adjust_views( TRUE, offm_AB, offn_AB, m_AB, n_AB, A, AB );

	// Free the temporary blocksize buffers.
	FLA_free( b_m );
	FLA_free( b_n );

	return FLA_SUCCESS;
}

FLA_Error FLASH_Part_create_1x2( FLA_Obj A,    FLA_Obj* AL, FLA_Obj* AR,
                                 dim_t n_cols, FLA_Side side )
{
	FLA_Datatype dt_A;
	dim_t        m_A,  n_A;
	dim_t        m_A_base, n_A_base;
	dim_t        m_AL, n_AL;
	dim_t        m_AR, n_AR;
	dim_t        depth;
	dim_t*       b_m;
	dim_t*       b_n;
	dim_t        offm_A,  offn_A;
	dim_t        offm_AL, offn_AL;
	dim_t        offm_AR, offn_AR;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Part_1x2_check( A,    AL, AR,     n_cols, side );

	// Safeguard: if n_cols > n, reduce n_cols to n.
	if ( n_cols > FLASH_Obj_scalar_width( A ) )
		n_cols = FLASH_Obj_scalar_width( A );

	// Acquire various properties of the hierarchical matrix object.
	dt_A     = FLASH_Obj_datatype( A );
	m_A      = FLASH_Obj_scalar_length( A );
	n_A      = FLASH_Obj_scalar_width( A );
	offm_A   = FLASH_Obj_scalar_row_offset( A );
	offn_A   = FLASH_Obj_scalar_col_offset( A );
	m_A_base = FLASH_Obj_base_scalar_length( A );
	n_A_base = FLASH_Obj_base_scalar_width( A );
	depth    = FLASH_Obj_depth( A );

	// Allocate a pair of temporary arrays for the blocksizes, whose lengths
	// are equal to the object's hierarchical depth.
	b_m = ( dim_t* ) FLA_malloc( depth * sizeof( dim_t ) );
	b_n = ( dim_t* ) FLA_malloc( depth * sizeof( dim_t ) );

	// Accumulate the blocksizes into the blocksize buffers.
	FLASH_Obj_blocksizes( A, b_m, b_n );

	// Adjust n_cols to be (n - n_cols) if the side specified is on the
	// right so that the right values get assigned below.
	if ( side == FLA_RIGHT ) n_cols = n_A - n_cols;

	// Set the dimensions of the partitions.
	m_AL = m_A;
	n_AL = n_cols;
	m_AR = m_A;
	n_AR = n_A - n_cols;

	// Set the offsets.
	offm_AL = offm_A + 0;
	offn_AL = offn_A + 0;
	offm_AR = offm_A + 0;
	offn_AR = offn_A + n_AL;
	
	// Create bufferless hierarhical objects that have the desired dimensions
	// for the views.
	FLASH_Obj_create_without_buffer_ext( dt_A, m_A_base, n_A_base, depth, b_m, b_n, AL );
	FLASH_Obj_create_without_buffer_ext( dt_A, m_A_base, n_A_base, depth, b_m, b_n, AR );

	// Recursively walk the hierarchy and adjust the views so that they
	// collectively refer to the absolute offsets given, and attach the
	// leaf-level numerical buffers of A to the new views.
	FLASH_Obj_adjust_views( TRUE, offm_AL, offn_AL, m_AL, n_AL, A, AL );
	FLASH_Obj_adjust_views( TRUE, offm_AR, offn_AR, m_AR, n_AR, A, AR );

	// Free the temporary blocksize buffers.
	FLA_free( b_m );
	FLA_free( b_n );

	return FLA_SUCCESS;
}

FLA_Error FLASH_Part_create_2x2( FLA_Obj A,    FLA_Obj* ATL, FLA_Obj* ATR,
                                               FLA_Obj* ABL, FLA_Obj* ABR,
                                 dim_t n_rows, dim_t n_cols, FLA_Side side )
{
	FLA_Datatype dt_A;
	dim_t        m_A_base, n_A_base;
	dim_t        m_A,  n_A;
	dim_t        m_ATL, n_ATL;
	dim_t        m_ABL, n_ABL;
	dim_t        m_ATR, n_ATR;
	dim_t        m_ABR, n_ABR;
	dim_t        depth;
	dim_t*       b_m;
	dim_t*       b_n;
	dim_t        offm_A,  offn_A;
	dim_t        offm_ATL, offn_ATL;
	dim_t        offm_ABL, offn_ABL;
	dim_t        offm_ATR, offn_ATR;
	dim_t        offm_ABR, offn_ABR;

	if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
		FLA_Part_2x2_check( A,   ATL, ATR,
		                         ABL, ABR,   n_rows, n_cols, side );

	// Safeguard: if n_rows > m, reduce n_rows to m.
	if ( n_rows > FLASH_Obj_scalar_length( A ) )
		n_rows = FLASH_Obj_scalar_length( A );

	// Safeguard: if n_cols > n, reduce n_cols to n.
	if ( n_cols > FLASH_Obj_scalar_width( A ) )
		n_cols = FLASH_Obj_scalar_width( A );

	// Acquire various properties of the hierarchical matrix object.
	dt_A     = FLASH_Obj_datatype( A );
	m_A      = FLASH_Obj_scalar_length( A );
	n_A      = FLASH_Obj_scalar_width( A );
	offm_A   = FLASH_Obj_scalar_row_offset( A );
	offn_A   = FLASH_Obj_scalar_col_offset( A );
	m_A_base = FLASH_Obj_base_scalar_length( A );
	n_A_base = FLASH_Obj_base_scalar_width( A );
	depth    = FLASH_Obj_depth( A );

	// Allocate a pair of temporary arrays for the blocksizes, whose lengths
	// are equal to the object's hierarchical depth.
	b_m = ( dim_t* ) FLA_malloc( depth * sizeof( dim_t ) );
	b_n = ( dim_t* ) FLA_malloc( depth * sizeof( dim_t ) );

	// Accumulate the blocksizes into the blocksize buffers.
	FLASH_Obj_blocksizes( A, b_m, b_n );

	// Adjust n_rows to be (m - n_rows) if the quadrant specified is on
	// the bottom so that the right values get assigned below. Do the same
	// for n_cols.
	if ( side == FLA_BL || side == FLA_BR ) n_rows = m_A - n_rows;
	if ( side == FLA_TR || side == FLA_BR ) n_cols = n_A - n_cols;

	// Set the dimensions of the partitions.
	m_ATL = n_rows;
	n_ATL = n_cols;
	m_ABL = m_A - n_rows;
	n_ABL = n_cols;
	m_ATR = n_rows;
	n_ATR = n_A - n_cols;
	m_ABR = m_A - n_rows;
	n_ABR = n_A - n_cols;

	// Set the offsets.
	offm_ATL = offm_A + 0;
	offn_ATL = offn_A + 0;
	offm_ABL = offm_A + m_ATL;
	offn_ABL = offn_A + 0;
	offm_ATR = offm_A + 0;
	offn_ATR = offn_A + n_ATL;
	offm_ABR = offm_A + m_ATL;
	offn_ABR = offn_A + n_ATL;
	
	// Create bufferless hierarhical objects that have the desired dimensions
	// for the views.
	FLASH_Obj_create_without_buffer_ext( dt_A, m_A_base, n_A_base, depth, b_m, b_n, ATL );
	FLASH_Obj_create_without_buffer_ext( dt_A, m_A_base, n_A_base, depth, b_m, b_n, ABL );
	FLASH_Obj_create_without_buffer_ext( dt_A, m_A_base, n_A_base, depth, b_m, b_n, ATR );
	FLASH_Obj_create_without_buffer_ext( dt_A, m_A_base, n_A_base, depth, b_m, b_n, ABR );

	// Recursively walk the hierarchy and adjust the views so that they
	// collectively refer to the absolute offsets given, and attach the
	// leaf-level numerical buffers of A to the new views.
	FLASH_Obj_adjust_views( TRUE, offm_ATL, offn_ATL, m_ATL, n_ATL, A, ATL );
	FLASH_Obj_adjust_views( TRUE, offm_ABL, offn_ABL, m_ABL, n_ABL, A, ABL );
	FLASH_Obj_adjust_views( TRUE, offm_ATR, offn_ATR, m_ATR, n_ATR, A, ATR );
	FLASH_Obj_adjust_views( TRUE, offm_ABR, offn_ABR, m_ABR, n_ABR, A, ABR );

	// Free the temporary blocksize buffers.
	FLA_free( b_m );
	FLA_free( b_n );

	return FLA_SUCCESS;
}

FLA_Error FLASH_Obj_adjust_views( FLA_Bool attach_buffer, dim_t offm, dim_t offn, dim_t m, dim_t n, FLA_Obj A, FLA_Obj* S )
{
	
	FLASH_Obj_adjust_views_hierarchy( attach_buffer, offm, offn, m, n, A, S );

	return FLA_SUCCESS;
}

FLA_Error FLASH_Obj_adjust_views_hierarchy( FLA_Bool attach_buffer, dim_t offm, dim_t offn, dim_t m, dim_t n, FLA_Obj A, FLA_Obj* S )
{
	FLA_Obj ATL, ATR,
	        ABL, ABR;

	FLA_Obj STL, STR,
	        SBL, SBR;

	// Base case.
	if ( FLA_Obj_elemtype( A ) == FLA_SCALAR )
	{
		// Repartition to exclude elements above and to the left of our
		// submatrix of interest.
		FLA_Part_2x2( A,    &ATL, &ATR,
		                    &ABL, &ABR,    offm, offn, FLA_TL );
		FLA_Part_2x2( *S,   &STL, &STR,
		                    &SBL, &SBR,    offm, offn, FLA_TL );

		// Overwrite the existing views with ones that have updated offsets.
		A  = ABR;
		*S = SBR;

		// Repartition to exclude elements below and to the right of our
		// submatrix of interest.
		FLA_Part_2x2( A,    &ATL, &ATR,
		                    &ABL, &ABR,    m, n, FLA_TL );
		FLA_Part_2x2( *S,   &STL, &STR,
		                    &SBL, &SBR,    m, n, FLA_TL );

		// Overwrite the existing view of S with the view of A so that S
		// Refers to the correct base object.
		A  = ATL;
		*S = STL;

		// Adjust the _inner fields in the view to reflect the number of
		// elements we have in each dimension.
		S->m_inner = m;
		S->n_inner = n;

		// Copy over buffer, stride, and object ID information if requested.
		if ( attach_buffer )
		{
			// Copy over the address of the numerical data buffer and its
			// corresponding row and column strides. This is obviously
			// necessary since we are creating a hierarchial view into an
			// existing hierarhical matrix, not a separate/new matrix
			// altogether.
			S->base->buffer = A.base->buffer;
			S->base->rs     = A.base->rs;
			S->base->cs     = A.base->cs;

			// Copy over the id field of the original matrix. This is used
			// by SuperMatrix to distinguish between distinct hierarchical
			// matrices. Since again, we are not creating a new matrix, we
			// will use the original object's id value.
			S->base->id = A.base->id;
		}
	}
	else // if ( FLA_Obj_elemtype( A ) == FLA_MATRIX )
	{
		FLA_Obj AL,  AR,       A0,  A1,  A2;

		FLA_Obj SL,  SR,       S0,  S1,  S2;

		FLA_Obj A1T,           A01,
		        A1B,           A11,
		                       A21;

		FLA_Obj S1T,           S01,
		        S1B,           S11,
		                       S21;
		dim_t b_m_full;
		dim_t b_n_full;
		dim_t offm_relA;
		dim_t offn_relA;
		dim_t offm_abs;
		dim_t offn_abs;
		dim_t offm_cur;
		dim_t offn_cur;
		dim_t offm_rem;
		dim_t offn_rem;
		dim_t offm_next;
		dim_t offn_next;
		dim_t m_next;
		dim_t n_next;
		dim_t m_ahead;
		dim_t n_ahead;
		dim_t m_behind;
		dim_t n_behind;

		// Acquire the scalar length and width of the top-left (full) block
		// at the current hierarchical level.
		b_m_full = FLASH_Obj_scalar_length_tl( A );
		b_n_full = FLASH_Obj_scalar_width_tl( A );
/*
printf( "-----------------\n" );
printf( "b_m/n_full:    %d %d\n", b_m_full, b_n_full );
printf( "offm/n:        %d %d\n", offm, offn );
printf( "r/c offsets:   %d %d\n", FLA_Obj_row_offset( A ), FLA_Obj_col_offset( A ) );
*/		
		// Compute the offsets for the top-left corner of the submatrix of
		// interest relative to the view at the current level of the
		// hierarchy of A.
		offm_relA = offm / b_m_full - FLA_Obj_row_offset( A );
		offn_relA = offn / b_n_full - FLA_Obj_col_offset( A );

		// Compute the offsets for the top-left corner of the submatrix of
		// interest in absolute units, from the top-left edge of the 
		// overall allocated matrix. This will be used to partition into S
		// Since its view has (presumably) not yet been changed since it
		// was created.
		offm_abs  = offm / b_m_full;
		offn_abs  = offn / b_n_full;
/*
printf( "offm/n_relA:   %d %d\n", offm_relA, offn_relA );
printf( "offm/n_abs:    %d %d\n", offm_abs, offn_abs );
*/
		// Repartition to exclude blocks above and to the left of our
		// submatrix of interest.
		FLA_Part_2x2( A,    &ATL, &ATR,
		                    &ABL, &ABR,    offm_relA, offn_relA, FLA_TL );
		FLA_Part_2x2( *S,   &STL, &STR,
		                    &SBL, &SBR,    offm_abs, offn_abs, FLA_TL );
/*
printf( "ABR.offm/n     %d %d\n", FLA_Obj_row_offset( ABR ), FLA_Obj_col_offset( ABR ) );
printf( "ABR is         %d %d\n", FLA_Obj_length( ABR ), FLA_Obj_width( ABR ) );
printf( "SBR.offm/n     %d %d\n", FLA_Obj_row_offset( SBR ), FLA_Obj_col_offset( SBR ) );
printf( "SBR is         %d %d\n", FLA_Obj_length( SBR ), FLA_Obj_width( SBR ) );
*/
		// Overwrite the existing views with ones that have updated offsets
		// (for this level in the hierarchy).
		A = ABR;
		*S = SBR;

		// Compute the new offsets within SBR, which is the remaining
		// distance after you subtract out the distance spanned by the
		// partitioning we just did.
		offm_rem = offm - offm_abs * b_m_full;
		offn_rem = offn - offn_abs * b_n_full;

//printf( "offm/n_rem:    %d %d\n", offm_rem, offn_rem );

		// Compute a new set of offsets corresponding to the bottom-right
		// edge of the desired submatrix. We'll use this to partition away
		// the remaining (bottom and right) parts of the FLASH matrix at
		// this level.
		offm_cur = ( offm_rem + m ) / b_m_full;
		offn_cur = ( offn_rem + n ) / b_n_full;
		offm_cur += ( (offm_rem + m) % b_m_full ? 1 : 0 );
		offn_cur += ( (offn_rem + n) % b_n_full ? 1 : 0 );

//printf( "offm/n_cur:    %d %d\n", offm_cur, offn_cur );

		// Repartition to exclude blocks below and to the right of our
		// submatrix of interest.
		FLA_Part_2x2( A,    &ATL, &ATR,
		                    &ABL, &ABR,    offm_cur, offn_cur, FLA_TL );
		FLA_Part_2x2( *S,   &STL, &STR,
		                    &SBL, &SBR,    offm_cur, offn_cur, FLA_TL );
/*
printf( "ATL.offm/n     %d %d\n", FLA_Obj_row_offset( ATL ), FLA_Obj_col_offset( ATL ) );
printf( "ATL is         %d %d\n", FLA_Obj_length( ATL ), FLA_Obj_width( ATL ) );
printf( "STL.offm/n     %d %d\n", FLA_Obj_row_offset( STL ), FLA_Obj_col_offset( STL ) );
printf( "STL is         %d %d\n", FLA_Obj_length( STL ), FLA_Obj_width( STL ) );
*/

		// Overwrite the existing views with ones that have updated offsets
		// (for this level in the hierarchy).
		A = ATL;
		*S = STL;

		// Adjust the _inner fields in the view to reflect the number of
		// elements we will eventually have in each dimension.
		S->m_inner = m;
		S->n_inner = n;

		// Initialize a counter that keeps track of the n offset relative to
		// the top-left most edge of the submatrix of interest.
		n_behind = 0;

		FLA_Part_1x2(  A,    &AL,  &AR,      0, FLA_LEFT );
		FLA_Part_1x2( *S,    &SL,  &SR,      0, FLA_LEFT );

		while ( FLA_Obj_width( AL ) < FLA_Obj_width( A ) )
		{
			FLA_Repart_1x2_to_1x3( AL,  /**/ AR,        &A0, /**/ &A1, &A2,
			                       1, FLA_RIGHT );
			FLA_Repart_1x2_to_1x3( SL,  /**/ SR,        &S0, /**/ &S1, &S2,
			                       1, FLA_RIGHT );

			// -------------------------------------------------------------

			// Set the n offset for the next levels of recursion based
			// on which panel of A we are in.
			if ( FLA_Obj_width( AL ) == 0 ) offn_next = offn_rem;
			else                            offn_next = 0;

			// Compute the number of columns left to be visited in the
			// submatrix of interset.
			n_ahead = n - n_behind;

			// Set the n dimensions for the next level of recursion
			// depending on whether the submatrix continues beyond the
			// current block.
			if ( offn_next + n_ahead > b_n_full ) n_next = b_n_full - offn_next;
			else                                  n_next = n_ahead;

			// Initialize a counter that keeps track of the m offset relative
			// to the top-left most edge of the submatrix of interest.
			m_behind = 0;

			FLA_Part_2x1( A1,    &A1T,
			                     &A1B,       0, FLA_TOP );
			FLA_Part_2x1( S1,    &S1T,
			                     &S1B,       0, FLA_TOP );

			while ( FLA_Obj_length( A1T ) < FLA_Obj_length( A1 ) )
			{
				FLA_Repart_2x1_to_3x1( A1T,               &A01,
				                    /* ** */            /* ** */
				                                          &A11,
				                       A1B,               &A21,        1, FLA_BOTTOM );
				FLA_Repart_2x1_to_3x1( S1T,               &S01,
				                    /* ** */            /* ** */
				                                          &S11,
				                       S1B,               &S21,        1, FLA_BOTTOM );

				// -------------------------------------------------------------

				// Set the m offset for the next levels of recursion based
				// on which block of A1 we are in.
				if ( FLA_Obj_length( A1T ) == 0 ) offm_next = offm_rem;
				else                              offm_next = 0;

				// Compute the number of rows left to be visited in the
				// submatrix of interset.
				m_ahead = m - m_behind;

				// Set the m dimensions for the next level of recursion
				// depending on whether the submatrix continues beyond the
				// current block.
				if ( offm_next + m_ahead > b_m_full ) m_next = b_m_full - offm_next;
				else                                  m_next = m_ahead;

//printf( "offm/n_next m/n_next:    %d %d %d %d\n", offm_next, offn_next, m_next, n_next );
				// Recursively call ourselves with new, smaller offsets
				// and the submatrix corresponding to FLASH blocks captured by ABR.
				FLASH_Obj_adjust_views_hierarchy( attach_buffer,
				                                  offm_next,
				                                  offn_next,
				                                  m_next,
				                                  n_next,
				                                  *FLASH_OBJ_PTR_AT( A11 ),
				                                  FLASH_OBJ_PTR_AT( S11 ) );

				// Increment m_behind to keep track of our absolute m offset.
				m_behind += m_next;

				// -------------------------------------------------------------

				FLA_Cont_with_3x1_to_2x1( &A1T,               A01,
				                                              A11,
				                        /* ** */           /* ** */
				                          &A1B,               A21,     FLA_TOP );
				FLA_Cont_with_3x1_to_2x1( &S1T,               S01,
				                                              S11,
				                        /* ** */           /* ** */
				                          &S1B,               S21,     FLA_TOP );
			}

			// Increment n_behind to keep track of our absolute n offset.
			n_behind += n_next;

			// -------------------------------------------------------------

			FLA_Cont_with_1x3_to_1x2( &AL,  /**/ &AR,        A0, A1, /**/ A2,
			                          FLA_LEFT );
			FLA_Cont_with_1x3_to_1x2( &SL,  /**/ &SR,        S0, S1, /**/ S2,
			                          FLA_LEFT );
		}
	}

	return FLA_SUCCESS;
}





FLA_Error FLASH_Part_free_2x1( FLA_Obj* AT,
                               FLA_Obj* AB )
{
	FLASH_Obj_free_without_buffer( AT );
	FLASH_Obj_free_without_buffer( AB );

	return FLA_SUCCESS;
}

FLA_Error FLASH_Part_free_1x2( FLA_Obj* AL, FLA_Obj* AR )
{
	FLASH_Obj_free_without_buffer( AL );
	FLASH_Obj_free_without_buffer( AR );

	return FLA_SUCCESS;
}

FLA_Error FLASH_Part_free_2x2( FLA_Obj* ATL, FLA_Obj* ATR,
                               FLA_Obj* ABL, FLA_Obj* ABR )
{
	FLASH_Obj_free_without_buffer( ATL );
	FLASH_Obj_free_without_buffer( ATR );
	FLASH_Obj_free_without_buffer( ABL );
	FLASH_Obj_free_without_buffer( ABR );

	return FLA_SUCCESS;
}

dim_t FLASH_Obj_scalar_length( FLA_Obj H )
{
	FLA_Obj  HT,              H0,
	         HB,              H1,
	                          H2;
	FLA_Obj* H1p;

	dim_t b = 0;

	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
		return FLA_Obj_length( H );

	if ( FLA_Obj_length( H ) == 0 )
		return 0;

	FLA_Part_2x1( H,    &HT,
	                    &HB,            0, FLA_TOP );

	while ( FLA_Obj_length( HT ) < FLA_Obj_length( H ) )
	{
		FLA_Repart_2x1_to_3x1( HT,                &H0,
		                    /* ** */            /* ** */
		                                          &H1,
		                       HB,                &H2,        1, FLA_BOTTOM );

		/*------------------------------------------------------------*/

		H1p = FLASH_OBJ_PTR_AT( H1 );
		b += H1p->m_inner;

		/*------------------------------------------------------------*/

		FLA_Cont_with_3x1_to_2x1( &HT,                H0,
		                                              H1,
		                        /* ** */           /* ** */   
		                          &HB,                H2,     FLA_TOP );
	}
  
	return b;
}

dim_t FLASH_Obj_scalar_width( FLA_Obj H )
{
	FLA_Obj  HL,    HR,       H0,  H1,  H2;
	FLA_Obj* H1p;

	dim_t b = 0;

	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
		return FLA_Obj_width( H );

	if ( FLA_Obj_width( H ) == 0 )
		return 0;

	FLA_Part_1x2( H,    &HL,  &HR,      0, FLA_LEFT );

	while ( FLA_Obj_width( HL ) < FLA_Obj_width( H ) )
	{
		FLA_Repart_1x2_to_1x3( HL,  /**/ HR,        &H0, /**/ &H1, &H2,
		                       1, FLA_RIGHT );

		/*------------------------------------------------------------*/

		H1p = FLASH_OBJ_PTR_AT( H1 );
		b += H1p->n_inner;

		/*------------------------------------------------------------*/

		FLA_Cont_with_1x3_to_1x2( &HL,  /**/ &HR,        H0, H1, /**/ H2,
		                          FLA_LEFT );
	}

	return b;
}

dim_t FLASH_Obj_scalar_min_dim( FLA_Obj H )
{
	return min( FLASH_Obj_scalar_length( H ),
	            FLASH_Obj_scalar_width( H ) );
}

dim_t FLASH_Obj_scalar_max_dim( FLA_Obj H )
{
	return max( FLASH_Obj_scalar_length( H ),
	            FLASH_Obj_scalar_width( H ) );
}

dim_t FLASH_Obj_scalar_vector_dim( FLA_Obj H )
{
	return ( FLASH_Obj_scalar_length( H ) == 1 ? FLASH_Obj_scalar_width( H )
	                                           : FLASH_Obj_scalar_length( H ) );
}

dim_t FLASH_Obj_scalar_row_offset( FLA_Obj H )
{
	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
	{
		return FLA_Obj_row_offset( H );
	}
	else
	{
		dim_t b_m = FLASH_Obj_scalar_length_tl( H );

		return FLA_Obj_row_offset( H ) * b_m + 
		       FLASH_Obj_scalar_row_offset( *FLASH_OBJ_PTR_AT( H ) );
	}
}

dim_t FLASH_Obj_scalar_col_offset( FLA_Obj H )
{
	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
	{
		return FLA_Obj_col_offset( H );
	}
	else
	{
		dim_t b_n = FLASH_Obj_scalar_width_tl( H );

		return FLA_Obj_col_offset( H ) * b_n + 
		       FLASH_Obj_scalar_col_offset( *FLASH_OBJ_PTR_AT( H ) );
	}
}

dim_t FLASH_Obj_scalar_length_tl( FLA_Obj H )
{
	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
	{
		return FLA_Obj_base_length( H );
	}
	else
	{
		FLA_Obj* H00 = FLA_Obj_base_buffer( H );

		return FLASH_Obj_base_scalar_length( *H00 );
	}
}

dim_t FLASH_Obj_scalar_width_tl( FLA_Obj H )
{
	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
	{
		return FLA_Obj_base_width( H );
	}
	else
	{
		FLA_Obj* H00 = FLA_Obj_base_buffer( H );

		return FLASH_Obj_base_scalar_width( *H00 );
	}
}

FLA_Error FLASH_Obj_show( char* header, FLA_Obj H, char* elem_format, char* footer )
{
	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
	{
		// Display the flat object.
		FLA_Obj_show( header, H, elem_format, footer );
	}
	else
	{
		dim_t m_scalar;
		dim_t i_view;
		dim_t i_abs;
		dim_t offm_scalar;

		// We want to print all m rows in the FLASH view.
		m_scalar = FLASH_Obj_scalar_length( H );

		// Get the scalar offset of the overall FLASH view relative to the
		// top-left corner of the overall object to which the view belongs.
		offm_scalar = FLASH_Obj_scalar_row_offset( H );

//printf( "flash_view_show: %d\n", m_scalar );
//printf( "flash_view_show: %d\n", offm_scalar );
	
		printf( "%s\n", header );

		for ( i_view = 0; i_view < m_scalar; ++i_view )
		{
			// Convert the relative view index to an absolute index.
			i_abs = offm_scalar + i_view; 
		
			// Print the ith row of the FLASH object H.
			FLASH_Obj_show_hierarchy( H, i_abs, elem_format );
			printf( "\n" );
		}

		printf( "%s\n", footer );
	}

	return FLA_SUCCESS;
}

FLA_Error FLASH_Obj_show_hierarchy( FLA_Obj H, dim_t i, char* elem_format )
{
	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
	{
		FLA_Datatype datatype = FLA_Obj_datatype( H );
		dim_t        m        = FLA_Obj_width( H );
		dim_t        rs       = FLA_Obj_row_stride( H );
		dim_t        cs       = FLA_Obj_col_stride( H );
		dim_t        j;

		// At this point, i is an absolute row index. We subtract out the
		// row offset of the view so that the index is relative to the view.
		i = i - FLA_Obj_row_offset( H );

		if ( datatype == FLA_INT )
		{
			int* buffer = FLA_Obj_buffer_at_view( H );

			for ( j = 0; j < m; ++j )
			{
				printf( elem_format, buffer[ j*cs + i*rs ] );
				printf( " " );
			}
		}
		else if ( datatype == FLA_FLOAT )
		{
			float* buffer = FLA_Obj_buffer_at_view( H );

			for ( j = 0; j < m; ++j )
			{
				printf( elem_format, buffer[ j*cs + i*rs ] );
				printf( " " );
			}
		}
		else if ( datatype == FLA_DOUBLE )
		{
			double* buffer = FLA_Obj_buffer_at_view( H );

			for ( j = 0; j < m; ++j )
			{
				printf( elem_format, buffer[ j*cs + i*rs ] );
				printf( " " );
			}
		}
		else if ( datatype == FLA_COMPLEX )
		{
			scomplex* buffer = FLA_Obj_buffer_at_view( H );

			for ( j = 0; j < m; ++j )
			{
				printf( elem_format, buffer[ j*cs + i*rs ].real,
				                     buffer[ j*cs + i*rs ].imag );
				printf( " " );
			}
		}
		else if ( datatype == FLA_DOUBLE_COMPLEX )
		{
			dcomplex* buffer = FLA_Obj_buffer_at_view( H );

			for ( j = 0; j < m; ++j )
			{
				printf( elem_format, buffer[ j*cs + i*rs ].real,
				                     buffer[ j*cs + i*rs ].imag );
				printf( " " );
			}
		}
		else
		{
			FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
		}
	}
	else
	{
		FLA_Obj  HT,
		         HB;
		FLA_Obj  HBL,  HBR,     H10, H11, H12;
		dim_t b_m_scalar;
		dim_t offm_local;
		dim_t i_next;

		// Get the scalar length of the top-left block.
		b_m_scalar = FLASH_Obj_scalar_length_tl( H );

#if 0
printf( "\n------------------------\n" );
printf( "b_m_scalar      %d\n", b_m_scalar );
printf( "i               %d\n", i );
#endif

		// Compute the offset of the matrix block, relative to the current
		// view, that contains the ith row of the matrix.
		offm_local      = ( i ) / b_m_scalar - FLA_Obj_row_offset( H );
		i_next          = ( i ) % b_m_scalar;

#if 0
printf( "row offset        %d\n", FLA_Obj_row_offset( H ) );
printf( "offm_local        %d\n", offm_local );
printf( "i_next            %d\n", i_next );
#endif

		FLA_Part_2x1( H,    &HT,
		                    &HB,       offm_local, FLA_TOP );

		FLA_Part_1x2(  HB,    &HBL,  &HBR,      0, FLA_LEFT );

		while ( FLA_Obj_width( HBL ) < FLA_Obj_width( HB ) )
		{
			FLA_Repart_1x2_to_1x3( HBL,  /**/ HBR,        &H10, /**/ &H11, &H12,
			                       1, FLA_RIGHT );

			// ------------------------------------------------------

			FLASH_Obj_show_hierarchy( *FLASH_OBJ_PTR_AT( H11 ),
			                          i_next, elem_format );

			// ------------------------------------------------------

			FLA_Cont_with_1x3_to_1x2( &HBL,  /**/ &HBR,        H10, H11, /**/ H12,
			                          FLA_LEFT );
		}
	}

	return FLA_SUCCESS;
}

