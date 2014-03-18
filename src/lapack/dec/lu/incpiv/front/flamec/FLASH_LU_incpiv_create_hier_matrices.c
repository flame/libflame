
#include "FLAME.h"

FLA_Error FLASH_LU_incpiv_create_hier_matrices( FLA_Obj A_flat, dim_t depth, dim_t* b_flash, dim_t b_alg, FLA_Obj* A, FLA_Obj* p, FLA_Obj* L )
{
   FLA_Datatype datatype;
   dim_t        m, n;
   dim_t        one = 1;
   
   // *** The current LU_incpiv algorithm implemented assumes that
   // the matrix has a hierarchical depth of 1. We check for that here, because
   // we anticipate that we'll use a more general algorithm in the future, and
   // we don't want to forget to remove the constraint. ***
   if ( depth != 1 )
   {
      FLA_Print_message( "FLASH_LU_incpiv() currently only supports matrices of depth 1",
                         __FILE__, __LINE__ );
      FLA_Abort();
   }

   // Create hierarchical copy of matrix A_flat.
   FLASH_Obj_create_hier_copy_of_flat( A_flat, depth, b_flash, A );

   // Query the datatype of matrix A_flat.
   datatype = FLA_Obj_datatype( A_flat );
   
   // If the user passed in zero for b_alg, then we need to set the algorithmic
   // (inner) blocksize to a reasonable default value.
   if ( b_alg == 0 )
   {
      b_alg = FLASH_LU_incpiv_determine_alg_blocksize( *A );
   }

   // Query the element (not scalar) dimensions of the new hierarchical matrix.
   // This is done so we can create p and L with full blocks for the bottom
   // and right "edge cases" of A.
   m = FLA_Obj_length( *A );
   n = FLA_Obj_width ( *A );

   // Create hierarchical matrices p and L.
   FLASH_Obj_create_ext( FLA_INT,  m * b_flash[0], n, 
                         depth, b_flash, &one, 
                         p );
   
   FLASH_Obj_create_ext( datatype, m * b_flash[0], n * b_alg, 
                         depth, b_flash, &b_alg, 
                         L );
      
   return FLA_SUCCESS;
}


dim_t FLASH_LU_incpiv_determine_alg_blocksize( FLA_Obj A )
{
   dim_t b_alg;
   dim_t b_flash;

   // Acquire the storage blocksize.
   b_flash = FLA_Obj_length( *FLASH_OBJ_PTR_AT( A ) );

   // Scale the storage blocksize by a pre-defined scalar to arrive at a
   // reasonable algorithmic blocksize, but make sure it's at least 1.
   b_alg = ( dim_t ) max( ( double ) b_flash * FLA_LU_INNER_TO_OUTER_B_RATIO, 1 );

   return b_alg;
}
