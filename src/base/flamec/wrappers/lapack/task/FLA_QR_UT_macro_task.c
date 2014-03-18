
#include "FLAME.h"

extern fla_qrut_t* fla_qrut_cntl_leaf;

FLA_Error FLA_QR_UT_macro_task( FLA_Obj A, FLA_Obj T, fla_qrut_t* cntl )
{
   FLA_Error r_val;
   FLA_Obj   A_flat;

   if ( FLA_Obj_length( A ) > 1 )
   {
      FLASH_Obj_create_flat_copy_of_hier( A, &A_flat );
  
      r_val = FLA_QR_UT_internal( A_flat, T, fla_qrut_cntl_leaf );
  
      FLASH_Copy_flat_to_hier( A_flat, 0, 0, A );
  
      FLA_Obj_free( &A_flat );
   }
   else
   {
      r_val = FLA_QR_UT_task( *FLASH_OBJ_PTR_AT( A ), T, cntl );
   }

   return r_val;
}

