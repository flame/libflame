/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLASH_FS_incpiv_aux1( FLA_Obj A, FLA_Obj p, FLA_Obj L, FLA_Obj b, dim_t nb_alg )
{
   FLA_Obj ATL,   ATR,      A00, A01, A02,
           ABL,   ABR,      A10, A11, A12,
                            A20, A21, A22;

   FLA_Obj pTL,   pTR,      p00, p01, p02,
           pBL,   pBR,      p10, p11, p12,
                            p20, p21, p22;

   FLA_Obj LTL,   LTR,      L00, L01, L02,
           LBL,   LBR,      L10, L11, L12,
                            L20, L21, L22;

   FLA_Obj bT,              b0,
           bB,              b1,
                            b2;

   FLA_Obj p11_conf,
           p11_rest;

   FLA_Part_2x2( A,    &ATL, &ATR,
                       &ABL, &ABR,     0, 0, FLA_TL );

   FLA_Part_2x2( p,    &pTL, &pTR,
                       &pBL, &pBR,     0, 0, FLA_TL );

   FLA_Part_2x2( L,    &LTL, &LTR,
                       &LBL, &LBR,     0, 0, FLA_TL );

   FLA_Part_2x1( b,    &bT,
                       &bB,            0, FLA_TOP );

   while ( FLA_Obj_length( ATL ) < FLA_Obj_length( A ) &&
           FLA_Obj_width ( ATL ) < FLA_Obj_width ( A ) )
   {
      FLA_Repart_2x2_to_3x3( ATL, /**/ ATR,       &A00, /**/ &A01, &A02,
                          /* ************* */   /* ******************** */
                                                  &A10, /**/ &A11, &A12,
                             ABL, /**/ ABR,       &A20, /**/ &A21, &A22,
                             1, 1, FLA_BR );

      FLA_Repart_2x2_to_3x3( pTL, /**/ pTR,       &p00, /**/ &p01, &p02,
                          /* ************* */   /* ******************** */
                                                  &p10, /**/ &p11, &p12,
                             pBL, /**/ pBR,       &p20, /**/ &p21, &p22,
                             1, 1, FLA_BR );

      FLA_Repart_2x2_to_3x3( LTL, /**/ LTR,       &L00, /**/ &L01, &L02,
                          /* ************* */   /* ******************** */
                                                  &L10, /**/ &L11, &L12,
                             LBL, /**/ LBR,       &L20, /**/ &L21, &L22,
                             1, 1, FLA_BR );

      FLA_Repart_2x1_to_3x1( bT,                  &b0,
                          /* ** */              /* ** */
                                                  &b1,
                             bB,                  &b2,        1, FLA_BOTTOM );

      /*------------------------------------------------------------*/

      FLA_Part_2x1( *FLASH_OBJ_PTR_AT( p11 ),   &p11_conf,
                                                &p11_rest,
                    FLA_Obj_length( *FLASH_OBJ_PTR_AT( b1 ) ), FLA_TOP );


      FLA_Apply_pivots( FLA_LEFT, FLA_NO_TRANSPOSE,
                        p11_conf,
                        *FLASH_OBJ_PTR_AT( b1 ) );

      FLA_Trsv_external( FLA_LOWER_TRIANGULAR, FLA_NO_TRANSPOSE, FLA_UNIT_DIAG,
                         *FLASH_OBJ_PTR_AT( A11 ),
                         *FLASH_OBJ_PTR_AT( b1 ) );

      FLASH_FS_incpiv_aux2( L21,
                            A21, p21, b1,
                                      b2, nb_alg );

      /*------------------------------------------------------------*/

      FLA_Cont_with_3x3_to_2x2( &ATL, /**/ &ATR,       A00, A01, /**/ A02,
                                                       A10, A11, /**/ A12,
                             /* ************** */  /* ****************** */
                                &ABL, /**/ &ABR,       A20, A21, /**/ A22,
                                FLA_TL );

      FLA_Cont_with_3x3_to_2x2( &pTL, /**/ &pTR,       p00, p01, /**/ p02,
                                                       p10, p11, /**/ p12,
                             /* ************** */  /* ****************** */
                                &pBL, /**/ &pBR,       p20, p21, /**/ p22,
                                FLA_TL );

      FLA_Cont_with_3x3_to_2x2( &LTL, /**/ &LTR,       L00, L01, /**/ L02,
                                                       L10, L11, /**/ L12,
                             /* ************** */  /* ****************** */
                                &LBL, /**/ &LBR,       L20, L21, /**/ L22,
                                FLA_TL );

      FLA_Cont_with_3x1_to_2x1( &bT,                   b0,
                                                       b1,
                              /* ** */              /* ** */
                                &bB,                   b2,     FLA_TOP );
   }
   
   return FLA_SUCCESS;
}
