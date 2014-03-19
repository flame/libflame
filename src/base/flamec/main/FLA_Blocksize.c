/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


// This block of code is enabled by the preprocessor macro FLA_ENABLE_GOTO_INTERFACES.
// It assumes that libgoto provides the BLAS implementation, which includes special
// symbols that are assigned upon calling blas_set_parameter(). These symbols are
// integer variables that hold the ideal blocksizes for each of the four datatypes
// for the host architecture.
#ifdef FLA_ENABLE_GOTO_INTERFACES

static FLA_Bool first_time = TRUE;

extern void blas_set_parameter( void );

extern long sgemm_p, sgemm_q, sgemm_r;
extern long dgemm_p, dgemm_q, dgemm_r;
extern long cgemm_p, cgemm_q, cgemm_r;
extern long zgemm_p, zgemm_q, zgemm_r;

long fla_goto_gemm_blocksize[4][4];

#endif



fla_blocksize_t* FLA_Blocksize_create( dim_t b_s, dim_t b_d, dim_t b_c, dim_t b_z )
{
	fla_blocksize_t* bp;
	
	// Allocate memory for the blocksize structure.
	bp = ( fla_blocksize_t* ) FLA_malloc( sizeof(fla_blocksize_t) );
	
	// Assign the provided blocksize values into the corresponding fields.
	bp->s = b_s;
	bp->d = b_d;
	bp->c = b_c;
	bp->z = b_z;
	
	// Return a pointer to the structure.
	return bp;
}


void FLA_Blocksize_set( fla_blocksize_t* bp, dim_t b_s, dim_t b_d, dim_t b_c, dim_t b_z )
{
	// Assign the provided blocksize values into the corresponding fields.
    if ( b_s != 0 ) bp->s = b_s;
    if ( b_d != 0 ) bp->d = b_d;
    if ( b_c != 0 ) bp->c = b_c;
    if ( b_z != 0 ) bp->z = b_z;
}


void FLA_Blocksize_scale( fla_blocksize_t* bp, double factor )
{
	FLA_Error e_val;
	
	// Verify that the given blocksize pointer is valid.
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_null_pointer( bp );
		FLA_Check_error_code( e_val );
	}
	
	// Assign the provided blocksize values into the corresponding fields.
	bp->s = ( dim_t )( ( double ) bp->s * factor );
	bp->d = ( dim_t )( ( double ) bp->d * factor );
	bp->c = ( dim_t )( ( double ) bp->c * factor );
	bp->z = ( dim_t )( ( double ) bp->z * factor );
}


fla_blocksize_t* FLA_Blocksize_create_copy( fla_blocksize_t* bp )
{
	fla_blocksize_t* bp_copy;
	FLA_Error        e_val;
	
	// Verify that the given blocksize pointer is valid.
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_null_pointer( bp );
		FLA_Check_error_code( e_val );
	}

	// Allocate memory for the blocksize structure.
	bp_copy = ( fla_blocksize_t* ) FLA_malloc( sizeof(fla_blocksize_t) );
	
	// Assign the provided blocksize object's values into the corresponding
	// fields of the new object.
	bp_copy->s = bp->s;
	bp_copy->d = bp->d;
	bp_copy->c = bp->c;
	bp_copy->z = bp->z;
	
	// Return a pointer to the structure.
	return bp_copy;
}


void FLA_Blocksize_free( fla_blocksize_t* bp )
{
	FLA_free( bp );
}


dim_t FLA_Blocksize_extract( FLA_Datatype dt, fla_blocksize_t* bp )
{
	dim_t     b = 0;
	FLA_Error e_val;

	// Verify that the given blocksize pointer is valid.
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_null_pointer( bp );
		FLA_Check_error_code( e_val );
	}

	if      ( dt == FLA_FLOAT )
		b = bp->s;
	else if ( dt == FLA_DOUBLE )
		b = bp->d;
	else if ( dt == FLA_COMPLEX )
		b = bp->c;
	else if ( dt == FLA_DOUBLE_COMPLEX )
		b = bp->z;
	
	// Return the blocksize corresponding with the datatype.
	return b;
}


fla_blocksize_t* FLA_Query_blocksizes( FLA_Dimension dim )
{
	fla_blocksize_t* bp;
	
	// Create an FLA_Blocksize_t object.
	bp = FLA_Blocksize_create( 0, 0, 0, 0 );
	
	// Query the requested blocksize (m, k, or n dimension) for all of the
	// datatypes and package the results in an FLA_Blocksize_t structure.
	bp->s = FLA_Query_blocksize( FLA_FLOAT, dim );
	bp->d = FLA_Query_blocksize( FLA_DOUBLE, dim );
	bp->c = FLA_Query_blocksize( FLA_COMPLEX, dim );
	bp->z = FLA_Query_blocksize( FLA_DOUBLE_COMPLEX, dim );
	
	// Return a pointer to the structure.
	return bp;
}


dim_t FLA_Query_blocksize( FLA_Datatype dt, FLA_Dimension dim )
{
	dim_t b_val = 0;

#ifdef FLA_ENABLE_GOTO_INTERFACES

	int dt_index;
	int dim_index;

	if ( first_time )
	{
		long sgemm_min, dgemm_min, cgemm_min, zgemm_min;
		
		// Find the blocksizes associated with FLA_DIMENSION_MIN.
		sgemm_min = min( sgemm_p, sgemm_q );
		dgemm_min = min( dgemm_p, dgemm_q );
		cgemm_min = min( cgemm_p, cgemm_q );
		zgemm_min = min( zgemm_p, zgemm_q );

		// Set the values for each datatype and dimension constant.
		fla_goto_gemm_blocksize[FLA_S_INDEX][FLA_DIM_M_INDEX]   = sgemm_p;
		fla_goto_gemm_blocksize[FLA_S_INDEX][FLA_DIM_K_INDEX]   = sgemm_q;
		fla_goto_gemm_blocksize[FLA_S_INDEX][FLA_DIM_N_INDEX]   = sgemm_r;
		fla_goto_gemm_blocksize[FLA_S_INDEX][FLA_DIM_MIN_INDEX] = sgemm_min;

		fla_goto_gemm_blocksize[FLA_D_INDEX][FLA_DIM_M_INDEX]   = dgemm_p;
		fla_goto_gemm_blocksize[FLA_D_INDEX][FLA_DIM_K_INDEX]   = dgemm_q;
		fla_goto_gemm_blocksize[FLA_D_INDEX][FLA_DIM_N_INDEX]   = dgemm_r;
		fla_goto_gemm_blocksize[FLA_D_INDEX][FLA_DIM_MIN_INDEX] = dgemm_min;

		fla_goto_gemm_blocksize[FLA_C_INDEX][FLA_DIM_M_INDEX]   = cgemm_p;
		fla_goto_gemm_blocksize[FLA_C_INDEX][FLA_DIM_K_INDEX]   = cgemm_q;
		fla_goto_gemm_blocksize[FLA_C_INDEX][FLA_DIM_N_INDEX]   = cgemm_r;
		fla_goto_gemm_blocksize[FLA_C_INDEX][FLA_DIM_MIN_INDEX] = cgemm_min;

		fla_goto_gemm_blocksize[FLA_Z_INDEX][FLA_DIM_M_INDEX]   = zgemm_p;
		fla_goto_gemm_blocksize[FLA_Z_INDEX][FLA_DIM_K_INDEX]   = zgemm_q;
		fla_goto_gemm_blocksize[FLA_Z_INDEX][FLA_DIM_N_INDEX]   = zgemm_r;
		fla_goto_gemm_blocksize[FLA_Z_INDEX][FLA_DIM_MIN_INDEX] = zgemm_min;

		first_time = FALSE;
	}

	// Compute the index of the requested datatype.
	dt_index  = dt  & FLA_DTYPE_INDEX_MASK;
	dim_index = dim & FLA_DIM_INDEX_MASK;

	// Index into the array and choose the appropriate blocksize.
	b_val = ( dim_t ) fla_goto_gemm_blocksize[dt_index][dim_index];

#else

	// Assign the return value to a default sane blocksize in case
	// we cannot access the libgoto symbols.
	if      ( dim == FLA_DIMENSION_M )
		b_val = FLA_DEFAULT_M_BLOCKSIZE;
	else if ( dim == FLA_DIMENSION_K )
		b_val = FLA_DEFAULT_K_BLOCKSIZE;
	else if ( dim == FLA_DIMENSION_N )
		b_val = FLA_DEFAULT_N_BLOCKSIZE;
	else if ( dim == FLA_DIMENSION_MIN )
    {
		b_val = min( FLA_DEFAULT_M_BLOCKSIZE, FLA_DEFAULT_K_BLOCKSIZE );
		b_val = min( b_val, FLA_DEFAULT_N_BLOCKSIZE );
    }

#endif

	// Return the blocksize.
	return b_val;
}


dim_t FLA_Determine_blocksize( FLA_Obj A_unproc, FLA_Quadrant to_dir, fla_blocksize_t* bp )
{
	FLA_Error    e_val;
	FLA_Datatype datatype;
	dim_t        A_unproc_size;
	dim_t        typed_blocksize;
	dim_t        b;

	// Determine the size of the remaining portion of the matrix.
	A_unproc_size = FLA_determine_matrix_size( A_unproc, to_dir );
	
	// Determine the datatype of the matrix.
	datatype = FLA_Obj_datatype( A_unproc );
	
	// Determine the raw blocksize value from the blocksize structure.
	typed_blocksize = FLA_Blocksize_extract( datatype, bp );

	// Check blocksize for zero value.
	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_blocksize_value( typed_blocksize );
		FLA_Check_error_code( e_val );
	}

	// If the unprocessed partition is smaller than our blocksize allows,
	// we have to use it's length/width instead.
	b = min( A_unproc_size, typed_blocksize );
	
	// Return the computed blocksize.
	return b;
}


dim_t FLA_determine_matrix_size( FLA_Obj A_unproc, FLA_Quadrant to_dir )
{
	dim_t r_val = 0;

	// Determine the size of the matrix dimension along which we are moving.
	switch( to_dir )
	{
		case FLA_TOP:
		case FLA_BOTTOM:
		{
			r_val = FLA_Obj_length( A_unproc );
			break;
		}
		case FLA_LEFT:
		case FLA_RIGHT:
		{
			r_val = FLA_Obj_width( A_unproc );
			break;
		}
		case FLA_TL:
		case FLA_TR:
		case FLA_BL:
		case FLA_BR:
		{
			// We need to use min_dim() here because the matrix might be
			// rectangular.
			r_val = FLA_Obj_min_dim( A_unproc );
			break;
		}
	}
	
	return r_val;
}

