/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


FLA_Datatype FLASH_Obj_datatype( FLA_Obj H )
{
	return FLA_Obj_datatype( H );
}


dim_t FLASH_Obj_depth( FLA_Obj H )
{
	FLA_Elemtype elemtype;
	FLA_Obj*     buffer_H;
	dim_t        depth = 0;

	// Recurse through the hierarchy to the first leaf node. We initialize
	// the recursion here:
	elemtype = FLA_Obj_elemtype( H );
	buffer_H = ( FLA_Obj* ) FLA_Obj_base_buffer( H );

	while ( elemtype == FLA_MATRIX )
	{
		++depth;

		// Get the element type of the top-leftmost underlying object. Also,
		// get a pointer to the first element of the top-leftmost object and
		// assume that it is of type FLA_Obj* in case elemtype is once again
		// FLA_MATRIX.
		elemtype = FLA_Obj_elemtype( buffer_H[0] );
		buffer_H = ( FLA_Obj * ) FLA_Obj_base_buffer( buffer_H[0] );
	}

	// At this point, the value of depth represents the depth of the matrix
	// hierarchy.
	return depth;
}


dim_t FLASH_Obj_blocksizes( FLA_Obj H, dim_t* b_m, dim_t* b_n )
{
	FLA_Elemtype elemtype;
	FLA_Obj*     buffer_H;
	dim_t        depth = 0;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_blocksizes_check( H, b_m, b_n );

	// Recurse through the hierarchy to the first leaf node. We initialize
	// the recursion here:
	elemtype = FLA_Obj_elemtype( H );
	buffer_H = ( FLA_Obj* ) FLA_Obj_base_buffer( H );

	while ( elemtype == FLA_MATRIX )
	{
		b_m[depth] = FLA_Obj_base_length( buffer_H[0] );
		b_n[depth] = FLA_Obj_base_width( buffer_H[0] );
		++depth;

		// Get the element type of the top-leftmost underlying object. Also,
		// get a pointer to the first element of the top-leftmost object and
		// assume that it is of type FLA_Obj* in case elemtype is once again
		// FLA_MATRIX.
		elemtype = FLA_Obj_elemtype( buffer_H[0] );
		buffer_H = ( FLA_Obj * ) FLA_Obj_base_buffer( buffer_H[0] );
	}

	// At this point, the first depth elements of blocksizes have been filled
	// with the blocksizes of H's various hierarchical levels. Return the
	// matrix depth as a confirmation of how many blocksizes were found.
	return depth;
}

dim_t FLASH_Obj_base_scalar_length( FLA_Obj H )
{
	FLA_Obj* buffer;
	dim_t    m;
	dim_t    rs, cs;
	dim_t    i;
	dim_t    m_base = 0;

	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
		return FLA_Obj_base_length( H );

	// Notice we use the base buffer since we are interested in the
	// whole object, not just the part referened by the view.
	buffer = FLA_Obj_base_buffer( H );
	m      = FLA_Obj_base_length( H );
	rs     = FLA_Obj_row_stride( H );
	cs     = FLA_Obj_col_stride( H );

	// Add up the row dimensions of all the base objects in the 0th
	// column of objects.
	for ( i = 0; i < m; ++i )
	{
		FLA_Obj hij = buffer[ i*rs + 0*cs ];

		m_base += (hij.base)->m_inner;
	}

	return m_base;
}

dim_t FLASH_Obj_base_scalar_width( FLA_Obj H )
{
	FLA_Obj* buffer;
	dim_t    n;
	dim_t    rs, cs;
	dim_t    j;
	dim_t    n_base = 0;

	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
		return FLA_Obj_base_width( H );

	// Notice we use the base buffer since we are interested in the
	// whole object, not just the part referened by the view.
	buffer = FLA_Obj_base_buffer( H );
	n      = FLA_Obj_base_width( H );
	rs     = FLA_Obj_row_stride( H );
	cs     = FLA_Obj_col_stride( H );

	// Add up the column dimensions of all the base objects in the 0th
	// row of objects.
	for ( j = 0; j < n; ++j )
	{
		FLA_Obj hij = buffer[ 0*rs + j*cs ];

		n_base += (hij.base)->n_inner;
	}

	return n_base;
}

FLA_Error FLASH_Obj_create( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_mn, FLA_Obj* H )
{
	FLASH_Obj_create_helper( FALSE, datatype, m, n, depth, b_mn, b_mn, H );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_ext( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H )
{
	FLASH_Obj_create_helper( FALSE, datatype, m, n, depth, b_m, b_n, H );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_without_buffer( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_mn, FLA_Obj* H )
{
	FLASH_Obj_create_helper( TRUE, datatype, m, n, depth, b_mn, b_mn, H );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_without_buffer_ext( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H )
{
	FLASH_Obj_create_helper( TRUE, datatype, m, n, depth, b_m, b_n, H );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_helper( FLA_Bool without_buffer, FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H )
{
	dim_t     i;
	FLA_Obj   flat_matrix;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_create_helper_check( without_buffer, datatype, m, n, depth, b_m, b_n, H );

	if ( depth == 0 )
	{
		// Base case: create a single contiguous matrix block. If we are
		// creating an object with a buffer, then we use column-major order.
		if ( without_buffer == FALSE )
			FLA_Obj_create( datatype, m, n, 0, 0, H );
		else
			FLA_Obj_create_without_buffer( datatype, m, n, H );
	}
	else
	{
		// We need temporary arrays the same length as the blocksizes arrays.
		dim_t* elem_sizes_m  = ( dim_t * ) FLA_malloc( depth * sizeof( dim_t ) );
		dim_t* elem_sizes_n  = ( dim_t * ) FLA_malloc( depth * sizeof( dim_t ) );
		dim_t* depth_sizes_m = ( dim_t * ) FLA_malloc( depth * sizeof( dim_t ) );
		dim_t* depth_sizes_n = ( dim_t * ) FLA_malloc( depth * sizeof( dim_t ) );
		dim_t* m_offsets     = ( dim_t * ) FLA_malloc( depth * sizeof( dim_t ) );
		dim_t* n_offsets     = ( dim_t * ) FLA_malloc( depth * sizeof( dim_t ) );
		
		// Fill two sets of arrays: elem_sizes_m/elem_sizes_n and depth_sizes_m/
		// depth_sizes_n.
		//  - elem_sizes_m[i] will contain the number of numerical elements that span
		//    the row dimension of a block at the ith level of the hierarchy. This is
		//    just the product of all row blocksizes "internal" to and including the
		//    current blocking level. (The elem_sizes_n array tracks similar values
		//    in the column dimension.)
		//  - depth_sizes_m[i] is similar to elem_sizes_m[i]. The only difference is
		//    that instead of tracking the number of numerical elements in the row
		//    dimension, it tracks the number of "storage" blocks that span the m
		//    dimension of a block at the ith level, where the m dimension of a
		//    storage block is the block size given in b_m[depth-1], ie:
		//    the inner-most row dimension block size. (The depth_sizes_n array
		//    tracks similar values in the column dimension.)
		elem_sizes_m[depth-1]  = b_m[depth-1];
		elem_sizes_n[depth-1]  = b_n[depth-1];
		depth_sizes_m[depth-1] = 1;
		depth_sizes_n[depth-1] = 1;
		for ( i = depth - 1; i > 0; --i )
		{
			elem_sizes_m[i-1]  = elem_sizes_m[i]  * b_m[i-1];
			elem_sizes_n[i-1]  = elem_sizes_n[i]  * b_n[i-1];
			depth_sizes_m[i-1] = depth_sizes_m[i] * b_m[i-1];
			depth_sizes_n[i-1] = depth_sizes_n[i] * b_n[i-1];
		}
	
		// Initialize the m_offsets and n_offsets arrays to zero.
		for ( i = 0; i < depth; i++ )
		{
			m_offsets[i] = 0;
			n_offsets[i] = 0;
		}

		// Create a "flat" matrix object. All leaf-level child objects will refer
		// to various offsets within this object's buffer. Whether we create the
		// object with row- or column-major storage is moot, since either way it
		// will be a 1-by-mn length matrix which we will partition through later
		// on in FLASH_Obj_create_hierarchy(). Note that it is IMPORTANT that the
		// matrix be 1-by-mn, and NOT m-by-n, since we want to use the 1x2
		// partitioning routines to walk through it as we attach various parts of
		// the buffer to the matrix hierarchy.
		if ( without_buffer == FALSE )
			FLA_Obj_create( datatype, 1, m*n, 0, 0, &flat_matrix );
		else
			FLA_Obj_create_without_buffer( datatype, m, n, &flat_matrix );
		
		// Recursively create the matrix hierarchy.
		FLASH_Obj_create_hierarchy( datatype, m, n, depth, elem_sizes_m, elem_sizes_n, flat_matrix, H, 0, depth, depth_sizes_m, depth_sizes_n, m_offsets, n_offsets );
		
		// Free the flat_matrix object, but not its buffer. If we created a
		// normal object with a buffer, we don't want to free the buffer because
		// it is being used by the hierarchical objected we just created. If we
		// created a bufferless object, we don't want to free the buffer because
		// there was no buffer allocated in the first place.
		FLA_Obj_free_without_buffer( &flat_matrix );
		
		// Free the local arrays.
		FLA_free( elem_sizes_m );
		FLA_free( elem_sizes_n );
		FLA_free( depth_sizes_m );
		FLA_free( depth_sizes_n );
		FLA_free( m_offsets );
		FLA_free( n_offsets );
	}

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_hierarchy( FLA_Datatype datatype, dim_t m, dim_t n, dim_t depth, dim_t* elem_sizes_m, dim_t* elem_sizes_n, FLA_Obj flat_matrix, FLA_Obj* H, unsigned long id, dim_t depth_overall, dim_t* depth_sizes_m, dim_t* depth_sizes_n, dim_t* m_offsets, dim_t* n_offsets )
{
	dim_t    i, j, b;
	dim_t    next_m, next_n;
	dim_t    num_m, num_n;
	dim_t    m_inner, n_inner;
	dim_t    elem_size_m_cur;
	dim_t    elem_size_n_cur;
	FLA_Obj  FL, FR, F0, F1, F2;
	FLA_Obj* buffer_H;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_create_hierarchy_check( datatype, m, n, depth, elem_sizes_m, elem_sizes_n, flat_matrix, H, id, depth_overall, depth_sizes_m, depth_sizes_n, m_offsets, n_offsets );

	if ( depth == 0 )
	{
		// If we're asked to create a zero-depth matrix, we interpret that as
		// a request to create leaf-level objects using the remaining portion
		// of the segment of the flat_matrix buffer that was passed in.
		FLA_Obj_create_without_buffer( datatype, m, n, H );
		FLA_Obj_attach_buffer( FLA_Obj_buffer_at_view( flat_matrix ), 1, m, H );
#ifdef FLA_ENABLE_SUPERMATRIX
		FLASH_Queue_set_block_size( m * n * FLA_Obj_datatype_size( datatype ) );
#endif
		H->base->id = id;

		// Fill in the m_index and n_index variables, which identify the
		// location of the current leaf node, in units of storage blocks,
		// within the overall matrix.
		for ( i = 0; i < depth_overall; i++ )
		{
			H->base->m_index += m_offsets[i] * depth_sizes_m[i];
			H->base->n_index += n_offsets[i] * depth_sizes_n[i];
		}
	}
	else
	{
		// The "current" level's elem_size value. That is, the number of numerical
		// scalar elements along one side of a full block on the current level,
		// for the row and column dimensions.
		elem_size_m_cur = elem_sizes_m[0];
		elem_size_n_cur = elem_sizes_n[0];

		// Compute the number of rows and columns in the current hierarchical
		// level of blocking.
		num_m = m / elem_size_m_cur + ( (m % elem_size_m_cur) ? 1 : 0 );
		num_n = n / elem_size_n_cur + ( (n % elem_size_n_cur) ? 1 : 0 );

		// The total number of scalar elements contained within/below this level
		// of the hierarchy. (The edge cases are handled by the computation of
		// next_m and next_n below, since they are passed in as the new m and n
		// for the next recursive call.)
		m_inner = m;
		n_inner = n;
		
		// Create a matrix whose elements are FLA_Objs for the current level of
		// blocking.
		FLA_Obj_create_ext( datatype, FLA_MATRIX, num_m, num_n, m_inner, n_inner, 0, 0, H );

		if ( depth == depth_overall )
			id = H->base->id;
		else
			H->base->id = id;

		// Grab the buffer from the new hierarchical object. This is an array of
		// FLA_Objs.
		buffer_H = ( FLA_Obj* ) FLA_Obj_buffer_at_view( *H );
		
		// Prepare to partition through the flat matrix so we can further allocate
		// segments of it to the various hierarchical sub-matrices. (The second
		// case occurs when the current function is called with a flat_matrix
		// argument that was created without a buffer.)
		if ( FLA_Obj_buffer_at_view( flat_matrix ) != NULL )
			FLA_Part_1x2( flat_matrix, &FL, &FR, 0, FLA_LEFT );
		else
			FLA_Obj_create_without_buffer( datatype, 0, 0, &F1 );
		
		for ( j = 0; j < num_n; ++j )
		{
			// Determine the number of elements along the column dimension
			// that will be contained within the submatrix referenced by
			// the (i,j)th FLA_MATRIX element in the current matrix.
			if ( j != num_n-1 || (n % elem_size_n_cur) == 0 )
				next_n = elem_size_n_cur;
			else
				next_n = n % elem_size_n_cur;
			
			n_offsets[depth_overall-depth] = j;
                        
			for ( i = 0; i < num_m; ++i )
			{			
				// Determine the number of elements along the row dimension
				// that will be contained within the submatrix referenced by
				// the (i,j)th FLA_MATRIX element in the current matrix.
				if ( i != num_m-1 || (m % elem_size_m_cur) == 0 )
					next_m = elem_size_m_cur;
				else
					next_m = m % elem_size_m_cur;

				m_offsets[depth_overall-depth] = i;
				
				// Partition the next m*n elements from the flat matrix so we can
				// "attach" them to the hierarchical matrices contained within the
				// (i,j)th FLA_MATRIX object.
				if ( FLA_Obj_buffer_at_view( flat_matrix ) != NULL )
				{
					b = min( FLA_Obj_width( FR ), next_m * next_n );
					FLA_Repart_1x2_to_1x3( FL,  /**/ FR,        &F0, /**/ &F1, &F2,
					                       b, FLA_RIGHT );
				}
				
				// Recursively call ourselves, with the appropriate parameters for
				// the next deeper level in the matrix hierarchy.
				FLASH_Obj_create_hierarchy( datatype, next_m, next_n, depth-1, &elem_sizes_m[1], &elem_sizes_n[1], F1, &buffer_H[j*num_m + i], id, depth_overall, depth_sizes_m, depth_sizes_n, m_offsets, n_offsets );

				// Continue with the repartitioning.
				if ( FLA_Obj_buffer_at_view( flat_matrix ) != NULL )
				{
					FLA_Cont_with_1x3_to_1x2( &FL,  /**/ &FR,        F0, F1, /**/ F2,
					                          FLA_LEFT );
				}
			}
		}

		// Free the temporary flat matrix subpartition object, but only if it was
		// created to begin with. Since it would have been created without a
		// buffer, we must free it in a similar manner.
		if ( FLA_Obj_buffer_at_view( flat_matrix ) == NULL )
			FLA_Obj_free_without_buffer( &F1 );
	}
	
	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_conf_to( FLA_Trans trans, FLA_Obj H, FLA_Obj* H_new )
{
	FLA_Datatype datatype;
	dim_t        m_base, n_base;
	dim_t        m_view, n_view;
	dim_t        offm_scalar, offn_scalar;
	dim_t        depth;
	dim_t*       b_m;
	dim_t*       b_n;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_create_conf_to_check( trans, H, H_new );

	// Acquire some properties of the hierarchical matrix object.
	datatype    = FLASH_Obj_datatype( H );
	m_base      = FLASH_Obj_base_scalar_length( H );
	n_base      = FLASH_Obj_base_scalar_width( H );
	m_view      = FLASH_Obj_scalar_length( H );
	n_view      = FLASH_Obj_scalar_width( H );
	offm_scalar = FLASH_Obj_scalar_row_offset( H );
	offn_scalar = FLASH_Obj_scalar_col_offset( H );
	depth       = FLASH_Obj_depth( H );

	// Allocate a pair of temporary arrays for the blocksizes, whose lengths
	// are equal to the object's hierarchical depth.
	b_m = ( dim_t* ) FLA_malloc( depth * sizeof( dim_t ) );
	b_n = ( dim_t* ) FLA_malloc( depth * sizeof( dim_t ) );

	// Accumulate the blocksizes into the blocksize buffer.
	FLASH_Obj_blocksizes( H, b_m, b_n );

	// Handle the transposition, if requested.
	if ( trans == FLA_TRANSPOSE )
	{
		FLA_Check_error_code( FLA_NOT_YET_IMPLEMENTED );
	}

	// Create the new hierarchical matrix object with the same base dimensions
	// as the original object..
	FLASH_Obj_create_ext( datatype, m_base, n_base, depth, b_m, b_n, H_new ); 

	// Adjust the hierarchical view of the new object to match that of the
	// original object.
	FLASH_Obj_adjust_views( FALSE, offm_scalar, offn_scalar, m_view, n_view, H, H_new );

	// Free the temporary blocksize buffers.
	FLA_free( b_m );
	FLA_free( b_n );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_hier_conf_to_flat( FLA_Trans trans, FLA_Obj F, dim_t depth, dim_t* b_mn, FLA_Obj* H )
{
	FLA_Datatype datatype;
	dim_t        m_H, n_H;
	dim_t        m_F, n_F;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_create_hier_conf_to_flat_check( trans, F, depth, b_mn, H );

	// Acquire the numerical datatype, length, and width of the flat matrix
	// object.
	datatype = FLA_Obj_datatype( F );
	m_F      = FLA_Obj_length( F );
	n_F      = FLA_Obj_width( F );

	// Factor in the transposition, if requested.
	if ( trans == FLA_NO_TRANSPOSE )
	{
		m_H = m_F;
		n_H = n_F;
	}
	else
	{
		m_H = n_F;
		n_H = m_F;
	}

	// Create the hierarchical matrix object.
	FLASH_Obj_create( datatype, m_H, n_H, depth, b_mn, H );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_hier_conf_to_flat_ext( FLA_Trans trans, FLA_Obj F, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H )
{
	FLA_Datatype datatype;
	dim_t        m_H, n_H;
	dim_t        m_F, n_F;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_create_hier_conf_to_flat_ext_check( trans, F, depth, b_m, b_n, H );

	// Acquire the numerical datatype, length, and width of the flat matrix
	// object.
	datatype = FLA_Obj_datatype( F );
	m_F      = FLA_Obj_length( F );
	n_F      = FLA_Obj_width( F );

	// Factor in the transposition, if requested.
	if ( trans == FLA_NO_TRANSPOSE )
	{
		m_H = m_F;
		n_H = n_F;
	}
	else
	{
		m_H = n_F;
		n_H = m_F;
	}

	// Create the hierarchical matrix object.
	FLASH_Obj_create_ext( datatype, m_H, n_H, depth, b_m, b_n, H );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_flat_conf_to_hier( FLA_Trans trans, FLA_Obj H, FLA_Obj* F )
{
	FLA_Datatype datatype;
	dim_t        m_H, n_H;
	dim_t        m_F, n_F;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_create_flat_conf_to_hier_check( trans, H, F );

	// Acquire the numerical datatype, length, and width of the
	// hierarchical matrix object.
	datatype = FLASH_Obj_datatype( H );
	m_H      = FLASH_Obj_scalar_length( H );
	n_H      = FLASH_Obj_scalar_width( H );

	// Factor in the transposition, if requested.
	if ( trans == FLA_NO_TRANSPOSE )
	{
		m_F = m_H;
		n_F = n_H;
	}
	else
	{
		m_F = n_H;
		n_F = m_H;
	}

	// Create the flat matrix object. Default to column-major storage.
	FLA_Obj_create( datatype, m_F, n_F, 0, 0, F );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_copy_of( FLA_Trans trans, FLA_Obj H_cur, FLA_Obj* H_new )
{
	// Create a new object conformal to the current object.
	FLASH_Obj_create_conf_to( trans, H_cur, H_new );

	// This is a workaround until we implement a FLASH_Copyt().
	if ( trans == FLA_NO_TRANSPOSE || trans == FLA_CONJ_NO_TRANSPOSE )
	{
		// Copy the contents of the current object to the new object.
		FLASH_Copy( H_cur, *H_new );

		// NOTE: we don't currently honor requests to conjugate!
		// We could, if we had FLASH_Conj() implemented, but we don't
		// currently.
	}
	else // if ( trans == FLA_TRANSPOSE || trans == FLA_CONJ_TRANSPOSE )
	{
		FLA_Obj F, F_trans;

		FLASH_Obj_create_flat_copy_of_hier( H_cur, &F );
		FLA_Obj_create_copy_of( trans, F, &F_trans );
		FLASH_Obj_hierarchify( F_trans, *H_new );
		FLA_Obj_free( &F );
		FLA_Obj_free( &F_trans );
	}

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_hier_copy_of_flat( FLA_Obj F, dim_t depth, dim_t* b_mn, FLA_Obj* H )
{
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_create_hier_copy_of_flat_check( F, depth, b_mn, H );

	// Create a hierarchical object conformal to the flat object.
	FLASH_Obj_create_hier_conf_to_flat( FLA_NO_TRANSPOSE, F, depth, b_mn, H );

	// Initialize the contents of the hierarchical matrix object with the
	// contents of the flat matrix object.
	FLASH_Copy_flat_to_hier( F, 0, 0, *H );
	
	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_hier_copy_of_flat_ext( FLA_Obj F, dim_t depth, dim_t* b_m, dim_t* b_n, FLA_Obj* H )
{
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_create_hier_copy_of_flat_ext_check( F, depth, b_m, b_n, H );

	// Create a hierarchical object conformal to the flat object.
	FLASH_Obj_create_hier_conf_to_flat_ext( FLA_NO_TRANSPOSE, F, depth, b_m, b_n, H );

	// Initialize the contents of the hierarchical matrix object with the
	// contents of the flat matrix object.
	FLASH_Copy_flat_to_hier( F, 0, 0, *H );
	
	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_create_flat_copy_of_hier( FLA_Obj H, FLA_Obj* F )
{
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_create_flat_copy_of_hier_check( H, F );

	// Create a flat object conformal to the hierarchical object.
	FLASH_Obj_create_flat_conf_to_hier( FLA_NO_TRANSPOSE, H, F );

	// Flatten the hierarchical object's contents into the new flat object.
	FLASH_Copy_hier_to_flat( 0, 0, H, *F );

	return FLA_SUCCESS;
}


void FLASH_Obj_free( FLA_Obj* H )
{
	FLA_Obj* buffer_H;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_free_check( H );

	// Free the object according to whether it contains a hierarchy.
	if ( FLA_Obj_elemtype( *H ) == FLA_MATRIX )
	{
		// Extract a pointer to the data buffer that was parititioned across all
		// leaf-level submatrices.
		buffer_H = ( FLA_Obj * ) FLASH_Obj_extract_buffer( *H );
	
		// Free the data buffer. This works because FLASH_Obj_extract_buffer()
		// returns the starting address of the first element's buffer, which is
		// also the starting address of the entire buffer.
#ifdef FLA_ENABLE_SCC
		FLA_shfree( buffer_H );
#else
		FLA_buff_free( buffer_H );
#endif

		// All that remains now is to free the interior of the matrix hierarchy.
		// This includes non-leaf buffers and their corresponding base objects
		// as well as leaf-level base objects.
		FLASH_Obj_free_hierarchy( H );
	}
	else
	{
		// If the matrix has no hierarchy, treat it like a flat object.
		FLA_Obj_free( H );
	}
}


void FLASH_Obj_free_without_buffer( FLA_Obj* H )
{
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_free_without_buffer_check( H );

	// Free the object according to whether it contains a hierarchy.
	if ( FLA_Obj_elemtype( *H ) == FLA_MATRIX )
	{
		// Skip freeing the numerical data buffer, since the object was
		// presumably created with FLASH_Obj_create_without_buffer().

		// Free the interior of the matrix hierarchy. This includes non-leaf
		// buffers and their corresponding base objects as well as leaf-level
		// base objects.
		FLASH_Obj_free_hierarchy( H );
	}
	else
	{
		// If the matrix has no hierarchy, treat it like a flat object with
		// no internal data buffer.
		FLA_Obj_free_without_buffer( H );
	}
}


void FLASH_Obj_free_hierarchy( FLA_Obj* H )
{
	//dim_t    m_H, n_H, rs, cs, i, j;
	dim_t    i;
	dim_t    n_elem_alloc;
	FLA_Obj* buffer_H;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_free_hierarchy_check( H );

	// If the element type of H is FLA_SCALAR, then it has no children to
	// free, so free the base object. In order to avoid freeing the object's
	// data buffer, which would have already been freed en masse by now if
	// the calling function is FLASH_Obj_free(), we will call
	// FLA_Obj_free_without_buffer().
	if ( FLA_Obj_elemtype( *H ) == FLA_SCALAR )
	{
		FLA_Obj_free_without_buffer( H );
		return;
	}
	else
	{
		// Acquire the number of elements allocated when this node was
		// created.
		n_elem_alloc = FLA_Obj_num_elem_alloc( *H );

		// Acquire the array of objects contained inside of H.
		buffer_H = ( FLA_Obj* ) FLA_Obj_base_buffer( *H );

		// For each allocated submatrix in H...
		for ( i = 0; i < n_elem_alloc; ++i )
		{
			// Recurse with the ith element of the allocated buffer.
			FLASH_Obj_free_hierarchy( &buffer_H[i] );
		}

		// Finally, free the internal array of objects.
		FLA_Obj_free( H );
	}
}


void* FLASH_Obj_extract_buffer( FLA_Obj H )
{
	FLA_Elemtype elemtype;
	FLA_Obj*     buffer_H;

	// Recurse through the hierarchy to the first leaf node to gain access
	// to the address of the actual numerical data buffer (ie: the "flat"
	// matrix object used in FLASH_Obj_create()). We initialize the search
	// here:
	elemtype = FLA_Obj_elemtype( H );
	buffer_H  = ( FLA_Obj* ) FLA_Obj_base_buffer( H );

	while ( elemtype == FLA_MATRIX )
	{
		elemtype = FLA_Obj_elemtype( buffer_H[0] );
		buffer_H = ( FLA_Obj* ) FLA_Obj_base_buffer( buffer_H[0] );
	}

	// At this point, the value in buffer_H is a pointer to the array that
	// holds the numerical data.
	return ( void* ) buffer_H;
}


FLA_Error FLASH_Obj_flatten( FLA_Obj H, FLA_Obj F )
{
	FLASH_Copy_hier_to_flat( 0, 0, H, F );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_hierarchify( FLA_Obj F, FLA_Obj H )
{
	FLASH_Copy_flat_to_hier( F, 0, 0, H );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_attach_buffer( void* buffer, dim_t rs, dim_t cs, FLA_Obj* H )
{
	FLA_Obj      flat_matrix;
	dim_t        m_base, n_base;
	FLA_Datatype datatype;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_attach_buffer_check( buffer, rs, cs, H );

	// Extract the scalar dimensions of the base object(s) and get its
	// numerical datatype. (These fields will be set even if it has a NULL
	// buffer, which it probably does since this function was just invoked.)
	m_base   = FLASH_Obj_base_scalar_length( *H );
	n_base   = FLASH_Obj_base_scalar_width( *H );
	datatype = FLASH_Obj_datatype( *H );

	// Create a temporary conventional object and attach the given buffer.
	// Segments of this buffer will be partitioned out to the various
	// leaf-level matrices of the hierarchical matrix H.
	FLA_Obj_create_without_buffer( datatype, m_base, n_base, &flat_matrix );
	FLA_Obj_attach_buffer( buffer, rs, cs, &flat_matrix );

	// Recurse through the hierarchical matrix, assigning segments of
	// flat_matrix to the various leaf-level matrices, similar to what
	// we would do if we were creating the object outright.
	FLASH_Obj_attach_buffer_hierarchy( flat_matrix, H );

	// Free the object (but don't free the buffer!).
	FLA_Obj_free_without_buffer( &flat_matrix );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Obj_attach_buffer_hierarchy( FLA_Obj F, FLA_Obj* H )
{
	FLA_Obj FL,    FR,       F0,  F1,  F2;

	FLA_Obj HL,    HR,       H0,  H1,  H2;

	FLA_Obj F1T,              F01,
	        F1B,              F11,
	                          F21;

	FLA_Obj H1T,              H01,
	        H1B,              H11,
	                          H21;

	dim_t b_m, b_n;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
		FLASH_Obj_attach_buffer_hierarchy_check( F, H );

	if ( FLA_Obj_elemtype( *H ) == FLA_SCALAR )
	{
		// If we've recursed down to a leaf node, then we can simply attach
		// the matrix buffer to the current leaf-level submatrix.
		// Notice we use FLA_Obj_buffer_at_view() because we want to attach
		// the buffer address referenced by the view of F.
		FLA_Obj_attach_buffer( FLA_Obj_buffer_at_view( F ), 
		                       FLA_Obj_row_stride( F ),
		                       FLA_Obj_col_stride( F ), H );
	}
	else
	{
		FLA_Part_1x2( *H,    &HL,  &HR,      0, FLA_LEFT );

		FLA_Part_1x2(  F,    &FL,  &FR,      0, FLA_LEFT );

		while ( FLA_Obj_width( HL ) < FLA_Obj_width( *H ) )
		{

			FLA_Repart_1x2_to_1x3( HL,  /**/ HR,        &H0, /**/ &H1, &H2,
			                       1, FLA_RIGHT );

			b_n = FLASH_Obj_scalar_width( H1 );

			FLA_Repart_1x2_to_1x3( FL,  /**/ FR,        &F0, /**/ &F1, &F2,
			                       b_n, FLA_RIGHT );

			/*------------------------------------------------------------*/

			FLA_Part_2x1( H1,    &H1T,
			                     &H1B,            0, FLA_TOP );

			FLA_Part_2x1( F1,    &F1T,
			                     &F1B,            0, FLA_TOP );

			while ( FLA_Obj_length( H1T ) < FLA_Obj_length( H1 ) )
			{

				FLA_Repart_2x1_to_3x1( H1T,               &H01,
				                    /* ** */            /* ** */
				                                          &H11,
				                       H1B,               &H21,        1, FLA_BOTTOM );

				b_m = FLASH_Obj_scalar_length( H11 );

				FLA_Repart_2x1_to_3x1( F1T,               &F01,
				                    /* ** */            /* ** */
				                                          &F11,
				                       F1B,               &F21,      b_m, FLA_BOTTOM );

				/*------------------------------------------------------------*/

				FLASH_Obj_attach_buffer_hierarchy( F11,
				                                   FLASH_OBJ_PTR_AT( H11 ) );

				/*------------------------------------------------------------*/

				FLA_Cont_with_3x1_to_2x1( &H1T,               H01,
				                                              H11,
				                        /* ** */           /* ** */
				                          &H1B,               H21,     FLA_TOP );

				FLA_Cont_with_3x1_to_2x1( &F1T,               F01,
				                                              F11,
				                        /* ** */           /* ** */
				                          &F1B,               F21,     FLA_TOP );
			}

			/*------------------------------------------------------------*/

			FLA_Cont_with_1x3_to_1x2( &HL,  /**/ &HR,        H0, H1, /**/ H2,
			                          FLA_LEFT );

			FLA_Cont_with_1x3_to_1x2( &FL,  /**/ &FR,        F0, F1, /**/ F2,
			                          FLA_LEFT );

		}
	}

	return FLA_SUCCESS;
}


void FLASH_print_struct( FLA_Obj H )
{
	dim_t    m_H, n_H, rs, cs, i, j;
	FLA_Obj* buffer_temp;

	m_H = FLA_Obj_length( H );
	n_H = FLA_Obj_width( H );
	rs  = FLA_Obj_row_stride( H );
	cs  = FLA_Obj_col_stride( H );

	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
		FLASH_print_struct_helper( H, 0 );
	else
	{
		for ( j = 0; j < n_H; ++j )
		{
			for ( i = 0; i < m_H; ++i )
			{
				buffer_temp = ( FLA_Obj* ) FLA_Obj_buffer_at_view( H );

				FLASH_print_struct_helper( buffer_temp[ j*cs + i*rs ], 0 );
			}
		}
	}
}


void FLASH_print_struct_helper( FLA_Obj H, int indent )
{
	dim_t    m_H, n_H, rs, cs, i, j, k;
	FLA_Obj* buffer_temp;

	for ( i = 0; i < indent; ++i )
		fprintf( stdout, "  " );

	if ( FLA_Obj_elemtype( H ) == FLA_SCALAR )
	{
		fprintf( stdout, "LEAF (%3d | rs %3lu | cs %3lu | %3lu x %3lu | addr %p)\n",
		                 FLA_Obj_datatype( H ),
		                 FLA_Obj_row_stride( H ), FLA_Obj_col_stride( H ),
		                 FLA_Obj_length( H ), FLA_Obj_width( H ),
		                 FLA_Obj_buffer_at_view( H ) );
		fflush( stdout );
	}
	else
	{
		m_H = FLA_Obj_length( H );
		n_H = FLA_Obj_width( H );
		rs  = FLA_Obj_row_stride( H );
		cs  = FLA_Obj_col_stride( H );
		
		fprintf( stdout, "MATRIX (%lux%lu):%d - %p\n",
		                 m_H, n_H,
		                 FLA_Obj_datatype( H ),
		                 FLA_Obj_buffer_at_view( H ) );
		fflush( stdout );
		
		for ( j = 0; j < n_H; ++j )
		{
			for ( i = 0; i < m_H; ++i )
			{
				for ( k = 0; k < indent; ++k )
					fprintf( stdout, "  " );
				
				buffer_temp = ( FLA_Obj* ) FLA_Obj_buffer_at_view( H );

				FLASH_print_struct_helper( buffer_temp[ j*cs + i*rs ],
				                           indent + 1 );
			}
		}
	}
}


