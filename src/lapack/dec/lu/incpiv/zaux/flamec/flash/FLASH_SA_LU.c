/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_SA_LU( FLA_Obj B, FLA_Obj C, 
                       FLA_Obj D, FLA_Obj E, FLA_Obj p, FLA_Obj L, dim_t nb_alg, fla_lu_t* cntl )
{
   FLA_Obj DT,              D0,
           DB,              D1,
                            D2;

   FLA_Obj ET,              E0,
           EB,              E1,
                            E2;

   FLA_Obj pT,              p0,
           pB,              p1,
                            p2;

   FLA_Obj LT,              L0,
           LB,              L1,
                            L2;

   FLA_Part_2x1( D,    &DT,
                       &DB,            0, FLA_TOP );

   FLA_Part_2x1( E,    &ET,
                       &EB,            0, FLA_TOP );

   FLA_Part_2x1( p,    &pT,
                       &pB,            0, FLA_TOP );

   FLA_Part_2x1( L,    &LT,
                       &LB,            0, FLA_TOP );

   while ( FLA_Obj_length( DT ) < FLA_Obj_length( D ) )
   {
      FLA_Repart_2x1_to_3x1( DT,                &D0,
                          /* ** */            /* ** */
                                                &D1,
                             DB,                &D2,        1, FLA_BOTTOM );

      FLA_Repart_2x1_to_3x1( ET,                &E0,
                          /* ** */            /* ** */
                                                &E1,
                             EB,                &E2,        1, FLA_BOTTOM );

      FLA_Repart_2x1_to_3x1( pT,                &p0,
                          /* ** */            /* ** */
                                                &p1,
                             pB,                &p2,        1, FLA_BOTTOM );

      FLA_Repart_2x1_to_3x1( LT,                &L0,
                          /* ** */            /* ** */
                                                &L1,
                             LB,                &L2,        1, FLA_BOTTOM );

      /*------------------------------------------------------------*/

      if ( FLASH_Queue_get_enabled( ) )
      {
         // Enqueue
         ENQUEUE_FLASH_SA_LU( *FLASH_OBJ_PTR_AT( B ),
                              *FLASH_OBJ_PTR_AT( D1 ),
                              *FLASH_OBJ_PTR_AT( p1 ),
                              *FLASH_OBJ_PTR_AT( L1 ),
                              nb_alg,
                              FLA_Cntl_sub_lu( cntl ) );
      }
      else
      {
         // Execute leaf
         FLA_SA_LU_task( *FLASH_OBJ_PTR_AT( B ),
                         *FLASH_OBJ_PTR_AT( D1 ),
                         *FLASH_OBJ_PTR_AT( p1 ),
                         *FLASH_OBJ_PTR_AT( L1 ),
                         nb_alg,
                         FLA_Cntl_sub_lu( cntl ) );
      }
      
      FLASH_SA_FS( L1,
                   D1, p1, C,
                           E1, nb_alg, FLA_Cntl_sub_gemm1( cntl ) );

      /*------------------------------------------------------------*/

      FLA_Cont_with_3x1_to_2x1( &DT,                D0,
                                                    D1,
                              /* ** */           /* ** */
                                &DB,                D2,     FLA_TOP );

      FLA_Cont_with_3x1_to_2x1( &ET,                E0,
                                                    E1,
                              /* ** */           /* ** */
                                &EB,                E2,     FLA_TOP );

      FLA_Cont_with_3x1_to_2x1( &pT,                p0,
                                                    p1,
                              /* ** */           /* ** */
                                &pB,                p2,     FLA_TOP );

      FLA_Cont_with_3x1_to_2x1( &LT,                L0,
                                                    L1,
                              /* ** */           /* ** */
                                &LB,                L2,     FLA_TOP );
   }
   
   return FLA_SUCCESS;
}
