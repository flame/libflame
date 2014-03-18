
#include "FLAME.h"

FLA_Error FLASH_SA_FS( FLA_Obj L,
                       FLA_Obj D, FLA_Obj p, FLA_Obj C,
                                             FLA_Obj E, dim_t nb_alg, fla_gemm_t* cntl )
{
   FLA_Obj CL,    CR,       C0,  C1,  C2;

   FLA_Obj EL,    ER,       E0,  E1,  E2;

   FLA_Part_1x2( C,    &CL,  &CR,      0, FLA_LEFT );

   FLA_Part_1x2( E,    &EL,  &ER,      0, FLA_LEFT );

   while ( FLA_Obj_width( CL ) < FLA_Obj_width( C ) )
   {
      FLA_Repart_1x2_to_1x3( CL,  /**/ CR,        &C0, /**/ &C1, &C2,
                             1, FLA_RIGHT );

      FLA_Repart_1x2_to_1x3( EL,  /**/ ER,        &E0, /**/ &E1, &E2,
                             1, FLA_RIGHT );

      /*------------------------------------------------------------*/

      if ( FLASH_Queue_get_enabled( ) )
      {
         // Enqueue
         ENQUEUE_FLASH_SA_FS( *FLASH_OBJ_PTR_AT( L ),
                              *FLASH_OBJ_PTR_AT( D ),
                              *FLASH_OBJ_PTR_AT( p ),
                              *FLASH_OBJ_PTR_AT( C1 ),
                              *FLASH_OBJ_PTR_AT( E1 ),
                              nb_alg,
                              FLA_Cntl_sub_gemm( cntl ) );
      }
      else
      {
         // Execute leaf
         FLA_SA_FS_task( *FLASH_OBJ_PTR_AT( L ),
                         *FLASH_OBJ_PTR_AT( D ),
                         *FLASH_OBJ_PTR_AT( p ),
                         *FLASH_OBJ_PTR_AT( C1 ),
                         *FLASH_OBJ_PTR_AT( E1 ),
                         nb_alg,
                         FLA_Cntl_sub_gemm( cntl ) );
      }
      
      /*------------------------------------------------------------*/

      FLA_Cont_with_1x3_to_1x2( &CL,  /**/ &CR,        C0, C1, /**/ C2,
                                FLA_LEFT );

      FLA_Cont_with_1x3_to_1x2( &EL,  /**/ &ER,        E0, E1, /**/ E2,
                                FLA_LEFT );
   }
   
   return FLA_SUCCESS;
}
