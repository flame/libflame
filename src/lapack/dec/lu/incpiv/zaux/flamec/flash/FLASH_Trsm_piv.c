/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_Trsm_piv( FLA_Obj A, FLA_Obj B, FLA_Obj p, fla_trsm_t* cntl )
{
   FLA_Obj BL,    BR,       B0,  B1,  B2;

   FLA_Part_1x2( B,    &BL,  &BR,      0, FLA_LEFT );

   while ( FLA_Obj_width( BL ) < FLA_Obj_width( B ) )
   {
      FLA_Repart_1x2_to_1x3( BL,  /**/ BR,        &B0, /**/ &B1, &B2,
                             1, FLA_RIGHT );

      /*------------------------------------------------------------*/

      if ( FLASH_Queue_get_enabled( ) )
      {
         // Enqueue
         ENQUEUE_FLASH_Trsm_piv( *FLASH_OBJ_PTR_AT( A ),
                                 *FLASH_OBJ_PTR_AT( B1 ),
                                 *FLASH_OBJ_PTR_AT( p ),
                                 FLA_Cntl_sub_trsm( cntl ) );
      }
      else
      {
         // Execute leaf
         FLA_Trsm_piv_task( *FLASH_OBJ_PTR_AT( A ),
                            *FLASH_OBJ_PTR_AT( B1 ),
                            *FLASH_OBJ_PTR_AT( p ),
                            FLA_Cntl_sub_trsm( cntl ) );
      }

      /*------------------------------------------------------------*/

      FLA_Cont_with_1x3_to_1x2( &BL,  /**/ &BR,        B0, B1, /**/ B2,
                                FLA_LEFT );
   }
   
   return FLA_SUCCESS;
}
