/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifdef FLA_ENABLE_THREAD_SAFE_INTERFACES

FLA_Error FLA_LU_piv_macro_task_ts( FLA_cntl_init_s *FLA_cntl_init_i, FLA_Obj A, FLA_Obj p, fla_lu_t* cntl )
{
   FLA_Error r_val;
   FLA_Obj   A_flat;

   if ( FLA_Obj_length( A ) > 1 )
   {
      FLASH_Obj_create_flat_copy_of_hier_ts( FLA_cntl_init_i, A, &A_flat );
      
      r_val = FLA_LU_piv_task_ts( FLA_cntl_init_i, A_flat, p, cntl );
      
      FLASH_Copy_flat_to_hier_ts( FLA_cntl_init_i, A_flat, 0, 0, A );
      
      FLA_Obj_free_ts( FLA_cntl_init_i, &A_flat );
   }
   else
   {
      r_val = FLA_LU_piv_task_ts( FLA_cntl_init_i, *FLASH_OBJ_PTR_AT_TS( FLA_cntl_init_i, A ), p, cntl );
   }

   return r_val;
}

#endif

FLA_Error FLA_LU_piv_macro_task( FLA_Obj A, FLA_Obj p, fla_lu_t* cntl )
{
   FLA_Error r_val;
   FLA_Obj   A_flat;

   if ( FLA_Obj_length( A ) > 1 )
   {
      FLASH_Obj_create_flat_copy_of_hier( A, &A_flat );
      
      r_val = FLA_LU_piv_task( A_flat, p, cntl );
      
      FLASH_Copy_flat_to_hier( A_flat, 0, 0, A );
      
      FLA_Obj_free( &A_flat );
   }
   else
   {
      r_val = FLA_LU_piv_task( *FLASH_OBJ_PTR_AT( A ), p, cntl );
   }

   return r_val;
}

