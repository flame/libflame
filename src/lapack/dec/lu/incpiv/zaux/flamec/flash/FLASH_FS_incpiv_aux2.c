
#include "FLAME.h"

FLA_Error FLASH_FS_incpiv_aux2( FLA_Obj L,
                                FLA_Obj D, FLA_Obj p, FLA_Obj C,
                                                      FLA_Obj E, dim_t nb_alg )
{
   FLA_Obj LT,              L0,
           LB,              L1,
                            L2;

   FLA_Obj DT,              D0,
           DB,              D1,
                            D2;

   FLA_Obj pT,              p0,
           pB,              p1,
                            p2;

   FLA_Obj ET,              E0,
           EB,              E1,
                            E2;

   FLA_Part_2x1( L,    &LT,
                       &LB,            0, FLA_TOP );

   FLA_Part_2x1( D,    &DT,
                       &DB,            0, FLA_TOP );

   FLA_Part_2x1( p,    &pT,
                       &pB,            0, FLA_TOP );

   FLA_Part_2x1( E,    &ET,
                       &EB,            0, FLA_TOP );

   while ( FLA_Obj_length( DT ) < FLA_Obj_length( D ) )
   {
      FLA_Repart_2x1_to_3x1( LT,                &L0,
                          /* ** */            /* ** */
                                                &L1,
                             LB,                &L2,        1, FLA_BOTTOM );

      FLA_Repart_2x1_to_3x1( DT,                &D0,
                          /* ** */            /* ** */
                                                &D1,
                             DB,                &D2,        1, FLA_BOTTOM );

      FLA_Repart_2x1_to_3x1( pT,                &p0,
                          /* ** */            /* ** */
                                                &p1,
                             pB,                &p2,        1, FLA_BOTTOM );

      FLA_Repart_2x1_to_3x1( ET,                &E0,
                          /* ** */            /* ** */
                                                &E1,
                             EB,                &E2,        1, FLA_BOTTOM );

      /*------------------------------------------------------------*/
      
      FLA_SA_FS_blk( *FLASH_OBJ_PTR_AT( L1 ),
                     *FLASH_OBJ_PTR_AT( D1 ),
                     *FLASH_OBJ_PTR_AT( p1 ),
                     *FLASH_OBJ_PTR_AT( C ),
                     *FLASH_OBJ_PTR_AT( E1 ),
                     nb_alg );
      
      /*------------------------------------------------------------*/

      FLA_Cont_with_3x1_to_2x1( &LT,                L0,
                                                    L1,
                              /* ** */           /* ** */
                                &LB,                L2,     FLA_TOP );

      FLA_Cont_with_3x1_to_2x1( &DT,                D0,
                                                    D1,
                              /* ** */           /* ** */
                                &DB,                D2,     FLA_TOP );

      FLA_Cont_with_3x1_to_2x1( &pT,                p0,
                                                    p1,
                              /* ** */           /* ** */
                                &pB,                p2,     FLA_TOP );

      FLA_Cont_with_3x1_to_2x1( &ET,                E0,
                                                    E1,
                              /* ** */           /* ** */
                                &EB,                E2,     FLA_TOP );
   }
   
   return FLA_SUCCESS;
}
