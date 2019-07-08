/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_Copy_buffer_to_hier( dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs, dim_t i, dim_t j, FLA_Obj H )
{
	FLA_Obj      flat_matrix;
	FLA_Datatype datatype;
	FLA_Error    e_val;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_matrix_strides( m, n, rs, cs );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_submatrix_dims_and_offset( m, n, i, j, H );
		FLA_Check_error_code( e_val );
	}

	// Acquire the datatype from the hierarchical matrix object.
	datatype = FLASH_Obj_datatype( H );

	// Create a temporary conventional matrix object of the requested datatype
	// and dimensions and attach the given buffer containing the incoming data.
	FLA_Obj_create_without_buffer( datatype, m, n, &flat_matrix );
	FLA_Obj_attach_buffer( buffer, rs, cs, &flat_matrix );

	// Recurse through H, adding in the corresponding elements of flat_matrix,
	// starting at the (i,j) element offset.
	FLASH_Copy_flat_to_hier( flat_matrix, i, j, H );

	// Free the object (but don't free the buffer!).
	FLA_Obj_free_without_buffer( &flat_matrix );

	return FLA_SUCCESS;
}


FLA_Error FLASH_Copy_hier_to_buffer( dim_t i, dim_t j, FLA_Obj H, dim_t m, dim_t n, void* buffer, dim_t rs, dim_t cs )
{
	FLA_Obj      flat_matrix;
	FLA_Datatype datatype;
	FLA_Error    e_val;

	if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
	{
		e_val = FLA_Check_matrix_strides( m, n, rs, cs );
		FLA_Check_error_code( e_val );

		e_val = FLA_Check_submatrix_dims_and_offset( m, n, i, j, H );
		FLA_Check_error_code( e_val );
	}

	// Acquire the datatype from the hierarchical matrix object.
	datatype = FLASH_Obj_datatype( H );

	// Create a temporary conventional matrix object of the requested datatype
	// and dimensions and attach the given buffer containing the incoming data.
	FLA_Obj_create_without_buffer( datatype, m, n, &flat_matrix );
	FLA_Obj_attach_buffer( buffer, rs, cs, &flat_matrix );

	// Recurse through H, adding in the corresponding elements of flat_matrix,
	// starting at the (i,j) element offset.
	FLASH_Copy_hier_to_flat( i, j, H, flat_matrix );

	// Free the object (but don't free the buffer!).
	FLA_Obj_free_without_buffer( &flat_matrix );

	return FLA_SUCCESS;
}


#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLASH_Copy_flat_to_hier_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj F, dim_t i, dim_t j, FLA_Obj H )
{
	FLA_Obj HTL, HTR,
	        HBL, HBR;
	FLA_Obj HBR_tl, HBR_tr,
	        HBR_bl, HBR_br;
	dim_t   m, n;

	m = FLA_Obj_length( F );
	n = FLA_Obj_width( F );

	FLASH_Part_create_2x2_ts( FLA_cntl_init_i, H,   &HTL, &HTR,
	                            &HBL, &HBR,    i, j, FLA_TL );

	FLASH_Part_create_2x2_ts( FLA_cntl_init_i, HBR,   &HBR_tl, &HBR_tr,
	                              &HBR_bl, &HBR_br,    m, n, FLA_TL );

	FLASH_Copy_hierarchy_ts( FLA_cntl_init_i, FLA_FLAT_TO_HIER, F, &HBR_tl );

	FLASH_Part_free_2x2_ts( FLA_cntl_init_i, &HBR_tl, &HBR_tr,
	                     &HBR_bl, &HBR_br );

	FLASH_Part_free_2x2_ts( FLA_cntl_init_i, &HTL, &HTR,
	                     &HBL, &HBR );

	return FLA_SUCCESS;
}
#endif

FLA_Error FLASH_Copy_flat_to_hier( FLA_Obj F, dim_t i, dim_t j, FLA_Obj H )
{
	FLA_Obj HTL, HTR,
	        HBL, HBR;
	FLA_Obj HBR_tl, HBR_tr,
	        HBR_bl, HBR_br;
	dim_t   m, n;

	m = FLA_Obj_length( F );
	n = FLA_Obj_width( F );

	FLASH_Part_create_2x2( H,   &HTL, &HTR,
	                            &HBL, &HBR,    i, j, FLA_TL );

	FLASH_Part_create_2x2( HBR,   &HBR_tl, &HBR_tr,
	                              &HBR_bl, &HBR_br,    m, n, FLA_TL );

	FLASH_Copy_hierarchy( FLA_FLAT_TO_HIER, F, &HBR_tl );

	FLASH_Part_free_2x2( &HBR_tl, &HBR_tr,
	                     &HBR_bl, &HBR_br );

	FLASH_Part_free_2x2( &HTL, &HTR,
	                     &HBL, &HBR );

	return FLA_SUCCESS;
}


#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLASH_Copy_hier_to_flat_ts( FLA_cntl_init_s *FLA_cntl_init_i, dim_t i, dim_t j, FLA_Obj H, FLA_Obj F )
{
	FLA_Obj HTL, HTR,
	        HBL, HBR;
	FLA_Obj HBR_tl, HBR_tr,
	        HBR_bl, HBR_br;
	dim_t   m, n;

	m = FLA_Obj_length( F );
	n = FLA_Obj_width( F );

	FLASH_Part_create_2x2_ts( FLA_cntl_init_i, H,   &HTL, &HTR,
	                            &HBL, &HBR,    i, j, FLA_TL );

	FLASH_Part_create_2x2_ts( FLA_cntl_init_i, HBR,   &HBR_tl, &HBR_tr,
	                              &HBR_bl, &HBR_br,    m, n, FLA_TL );

	FLASH_Copy_hierarchy_ts( FLA_cntl_init_i, FLA_HIER_TO_FLAT, F, &HBR_tl );

	FLASH_Part_free_2x2_ts( FLA_cntl_init_i, &HBR_tl, &HBR_tr,
	                     &HBR_bl, &HBR_br );

	FLASH_Part_free_2x2_ts( FLA_cntl_init_i, &HTL, &HTR,
	                     &HBL, &HBR );

	return FLA_SUCCESS;
}
#endif

FLA_Error FLASH_Copy_hier_to_flat( dim_t i, dim_t j, FLA_Obj H, FLA_Obj F )
{
	FLA_Obj HTL, HTR,
	        HBL, HBR;
	FLA_Obj HBR_tl, HBR_tr,
	        HBR_bl, HBR_br;
	dim_t   m, n;

	m = FLA_Obj_length( F );
	n = FLA_Obj_width( F );

	FLASH_Part_create_2x2( H,   &HTL, &HTR,
	                            &HBL, &HBR,    i, j, FLA_TL );

	FLASH_Part_create_2x2( HBR,   &HBR_tl, &HBR_tr,
	                              &HBR_bl, &HBR_br,    m, n, FLA_TL );

	FLASH_Copy_hierarchy( FLA_HIER_TO_FLAT, F, &HBR_tl );

	FLASH_Part_free_2x2( &HBR_tl, &HBR_tr,
	                     &HBR_bl, &HBR_br );

	FLASH_Part_free_2x2( &HTL, &HTR,
	                     &HBL, &HBR );

	return FLA_SUCCESS;
}

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES
FLA_Error FLASH_Copy_hierarchy_ts( FLA_cntl_init_s *FLA_cntl_init_i, int direction, FLA_Obj F, FLA_Obj* H )
{
	// Once we get down to a submatrix whose elements are scalars, we are down
	// to our base case.
	if ( FLA_Obj_elemtype( *H ) == FLA_SCALAR )
	{
		// Depending on which top-level function invoked us, we either copy
		// the source data in the flat matrix to the leaf-level submatrix of
		// the hierarchical matrix, or copy the data in the hierarchical
		// submatrix to the flat matrix.
		if      ( direction == FLA_FLAT_TO_HIER )
		{
#ifdef FLA_ENABLE_SCC
			if ( FLA_is_owner() )
#endif
			FLA_Copy_external_ts( FLA_cntl_init_i, F, *H );
		}
		else if ( direction == FLA_HIER_TO_FLAT )
		{
#ifdef FLA_ENABLE_SCC
			if ( FLA_is_owner() )
#endif
			FLA_Copy_external_ts( FLA_cntl_init_i, *H, F );
		}
	}
	else
	{
		FLA_Obj HL,  HR,       H0,  H1,  H2;
		FLA_Obj FL,  FR,       F0,  F1,  F2;

		FLA_Obj H1T,           H01,
		        H1B,           H11,
		                       H21;
		FLA_Obj F1T,           F01,
		        F1B,           F11,
		                       F21;

		dim_t b_m;
		dim_t b_n;

		FLA_Part_1x2_ts( FLA_cntl_init_i, *H,    &HL,  &HR,      0, FLA_LEFT );
		FLA_Part_1x2_ts( FLA_cntl_init_i,  F,    &FL,  &FR,      0, FLA_LEFT );

		while ( FLA_Obj_width( HL ) < FLA_Obj_width( *H ) )
		{
			FLA_Repart_1x2_to_1x3_ts( FLA_cntl_init_i, HL,  /**/ HR,        &H0, /**/ &H1, &H2,
			                       1, FLA_RIGHT );

			// Get the scalar width of H1 and use that to determine the
			// width of F1.
			b_n = FLASH_Obj_scalar_width_ts( FLA_cntl_init_i, H1 );

			FLA_Repart_1x2_to_1x3_ts( FLA_cntl_init_i, FL,  /**/ FR,        &F0, /**/ &F1, &F2,
			                       b_n, FLA_RIGHT );

			// -------------------------------------------------------------

			FLA_Part_2x1_ts( FLA_cntl_init_i, H1,    &H1T,
			                     &H1B,       0, FLA_TOP );
			FLA_Part_2x1_ts( FLA_cntl_init_i, F1,    &F1T,
			                     &F1B,       0, FLA_TOP );

			while ( FLA_Obj_length( H1T ) < FLA_Obj_length( H1 ) )
			{
				FLA_Repart_2x1_to_3x1_ts( FLA_cntl_init_i, H1T,               &H01,
				                    /* ** */            /* *** */
				                                          &H11,
				                       H1B,               &H21,        1, FLA_BOTTOM );

				// Get the scalar length of H11 and use that to determine the
				// length of F11.
				b_m = FLASH_Obj_scalar_length_ts( FLA_cntl_init_i, H11 );

				FLA_Repart_2x1_to_3x1_ts( FLA_cntl_init_i, F1T,               &F01,
				                    /* ** */            /* *** */
				                                          &F11,
				                       F1B,               &F21,        b_m, FLA_BOTTOM );
				// -------------------------------------------------------------

				// Recursively copy between F11 and H11.
				FLASH_Copy_hierarchy_ts( FLA_cntl_init_i, direction, F11,
				                      FLASH_OBJ_PTR_AT_TS( FLA_cntl_init_i, H11 ) );

				// -------------------------------------------------------------

				FLA_Cont_with_3x1_to_2x1_ts( FLA_cntl_init_i, &H1T,               H01,
				                                              H11,
				                        /* ** */           /* *** */
				                          &H1B,               H21,     FLA_TOP );
				FLA_Cont_with_3x1_to_2x1_ts( FLA_cntl_init_i, &F1T,               F01,
				                                              F11,
				                        /* ** */           /* *** */
				                          &F1B,               F21,     FLA_TOP );
			}

			// -------------------------------------------------------------

			FLA_Cont_with_1x3_to_1x2_ts( FLA_cntl_init_i, &HL,  /**/ &HR,        H0, H1, /**/ H2,
			                          FLA_LEFT );
			FLA_Cont_with_1x3_to_1x2_ts( FLA_cntl_init_i, &FL,  /**/ &FR,        F0, F1, /**/ F2,
			                          FLA_LEFT );
		}
	}

	return FLA_SUCCESS;
}
#endif

FLA_Error FLASH_Copy_hierarchy( int direction, FLA_Obj F, FLA_Obj* H )
{
	// Once we get down to a submatrix whose elements are scalars, we are down
	// to our base case.
	if ( FLA_Obj_elemtype( *H ) == FLA_SCALAR )
	{
		// Depending on which top-level function invoked us, we either copy
		// the source data in the flat matrix to the leaf-level submatrix of
		// the hierarchical matrix, or copy the data in the hierarchical
		// submatrix to the flat matrix.
		if      ( direction == FLA_FLAT_TO_HIER )
		{
#ifdef FLA_ENABLE_SCC
			if ( FLA_is_owner() )
#endif
			FLA_Copy_external( F, *H );
		}
		else if ( direction == FLA_HIER_TO_FLAT )
		{
#ifdef FLA_ENABLE_SCC
			if ( FLA_is_owner() )
#endif
			FLA_Copy_external( *H, F );
		}
	}
	else
	{
		FLA_Obj HL,  HR,       H0,  H1,  H2;
		FLA_Obj FL,  FR,       F0,  F1,  F2;

		FLA_Obj H1T,           H01,
		        H1B,           H11,
		                       H21;
		FLA_Obj F1T,           F01,
		        F1B,           F11,
		                       F21;

		dim_t b_m;
		dim_t b_n;

		FLA_Part_1x2( *H,    &HL,  &HR,      0, FLA_LEFT );
		FLA_Part_1x2(  F,    &FL,  &FR,      0, FLA_LEFT );

		while ( FLA_Obj_width( HL ) < FLA_Obj_width( *H ) )
		{
			FLA_Repart_1x2_to_1x3( HL,  /**/ HR,        &H0, /**/ &H1, &H2,
			                       1, FLA_RIGHT );

			// Get the scalar width of H1 and use that to determine the
			// width of F1.
			b_n = FLASH_Obj_scalar_width( H1 );

			FLA_Repart_1x2_to_1x3( FL,  /**/ FR,        &F0, /**/ &F1, &F2,
			                       b_n, FLA_RIGHT );

			// -------------------------------------------------------------

			FLA_Part_2x1( H1,    &H1T,
			                     &H1B,       0, FLA_TOP );
			FLA_Part_2x1( F1,    &F1T,
			                     &F1B,       0, FLA_TOP );

			while ( FLA_Obj_length( H1T ) < FLA_Obj_length( H1 ) )
			{
				FLA_Repart_2x1_to_3x1( H1T,               &H01,
				                    /* ** */            /* *** */
				                                          &H11,
				                       H1B,               &H21,        1, FLA_BOTTOM );

				// Get the scalar length of H11 and use that to determine the
				// length of F11.
				b_m = FLASH_Obj_scalar_length( H11 );

				FLA_Repart_2x1_to_3x1( F1T,               &F01,
				                    /* ** */            /* *** */
				                                          &F11,
				                       F1B,               &F21,        b_m, FLA_BOTTOM );
				// -------------------------------------------------------------

				// Recursively copy between F11 and H11.
				FLASH_Copy_hierarchy( direction, F11,
				                      FLASH_OBJ_PTR_AT( H11 ) );

				// -------------------------------------------------------------

				FLA_Cont_with_3x1_to_2x1( &H1T,               H01,
				                                              H11,
				                        /* ** */           /* *** */
				                          &H1B,               H21,     FLA_TOP );
				FLA_Cont_with_3x1_to_2x1( &F1T,               F01,
				                                              F11,
				                        /* ** */           /* *** */
				                          &F1B,               F21,     FLA_TOP );
			}

			// -------------------------------------------------------------

			FLA_Cont_with_1x3_to_1x2( &HL,  /**/ &HR,        H0, H1, /**/ H2,
			                          FLA_LEFT );
			FLA_Cont_with_1x3_to_1x2( &FL,  /**/ &FR,        F0, F1, /**/ F2,
			                          FLA_LEFT );
		}
	}

	return FLA_SUCCESS;
}

