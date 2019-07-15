/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

extern __thread fla_appiv_t* fla_appiv_cntl_leaf;

FLA_Error FLA_Apply_pivots_macro_task( FLA_Side side, FLA_Trans trans, FLA_Obj p, FLA_Obj A, fla_appiv_t* cntl )
{
   FLA_Error r_val;
/*
   FLA_Obj   A_flat;

   FLASH_Obj_create_flat_copy_of_hier( A, &A_flat );

   r_val = FLA_Apply_pivots_unb_external( side, trans, p, A_flat );

   FLASH_Copy_flat_to_hier( A_flat, 0, 0, A );

   FLA_Obj_free( &A_flat );
*/
   if ( FLA_Obj_length( A ) > 1 )
   {
      r_val = FLA_Apply_pivots_macro_external( side, trans, p, A );
   }
   else
   {
      //r_val = FLA_Apply_pivots_unb_external( side, trans, p, 
      //                                       *FLASH_OBJ_PTR_AT( A ) );
      r_val = FLA_Apply_pivots_internal( side, trans, p, 
                                         *FLASH_OBJ_PTR_AT( A ),
                                         fla_appiv_cntl_leaf );
   }
   
   return r_val;
}

