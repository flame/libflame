
#include "FLAME.h"

//
// --- FLA_Part_2x2() ----------------------------------------------------------
//

FLA_Error FLA_Part_2x2( FLA_Obj A,  FLA_Obj *A11, FLA_Obj *A12,
                                    FLA_Obj *A21, FLA_Obj *A22, 
                        dim_t  mb,  dim_t     nb, FLA_Quadrant quadrant )
{
  FLA_Base_obj *base;
  dim_t         m, n, offm, offn;

  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Part_2x2_check( A,    A11, A12,
                              A21, A22,    mb, nb, quadrant );

  // Safeguard: if mb > m, reduce mb to m.
  if ( mb > A.m ) mb = A.m;

  // Safeguard: if nb > n, reduce nb to n.
  if ( nb > A.n ) nb = A.n;

  m        = A.m;
  n        = A.n;
  offm     = A.offm;
  offn     = A.offn;
  base     = A.base;

  // Set mb and nb to be the dimensions of A11.
  if ( quadrant == FLA_BL || quadrant == FLA_BR ) mb = m - mb;
  if ( quadrant == FLA_TR || quadrant == FLA_BR ) nb = n - nb;

  A11->m    = mb;
  A11->n    = nb;
  A11->offm = offm;
  A11->offn = offn;
  A11->base = base;

  A21->m    = m-mb;
  A21->n    = nb;
  A21->offm = offm + mb;
  A21->offn = offn;
  A21->base = base;

  A12->m    = mb;
  A12->n    = n-nb;
  A12->offm = offm;
  A12->offn = offn + nb;
  A12->base = base;

  A22->m    = m-mb;
  A22->n    = n-nb;
  A22->offm = offm + mb;
  A22->offn = offn + nb;
  A22->base = base;

  return FLA_SUCCESS;
}


//
// --- FLA_Part_2x1() ----------------------------------------------------------
//

FLA_Error FLA_Part_2x1( FLA_Obj A,  FLA_Obj *A1, 
                                    FLA_Obj *A2,
                        dim_t  mb,  FLA_Side side )
{ 
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Part_2x1_check( A,    A1,
                              A2,     mb, side );

  // Safeguard: if mb > m, reduce mb to m.
  if ( mb > A.m ) mb = A.m;

  // Set mb to be the dimension of A1.
  if ( side == FLA_BOTTOM ) mb = A.m - mb;

  A1->m    = mb;
  A1->n    = A.n;
  A1->offm = A.offm;
  A1->offn = A.offn;
  A1->base = A.base;

  A2->m    = A.m - mb;
  A2->n    = A.n;
  A2->offm = A.offm + mb;
  A2->offn = A.offn;
  A2->base = A.base;

  return FLA_SUCCESS;
}


//
// --- FLA_Part_1x2() ----------------------------------------------------------
//

FLA_Error FLA_Part_1x2( FLA_Obj A,  FLA_Obj *A1, FLA_Obj *A2,
                                    dim_t    nb, FLA_Side side )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Part_1x2_check( A,    A1,  A2,      nb, side );

  // Safeguard: if nb > n, reduce nb to n.
  if ( nb > A.n ) nb = A.n;

  // Set nb to be the dimension of A1.
  if ( side == FLA_RIGHT ) nb = A.n - nb;

  A1->m    = A.m;
  A1->n    = nb;
  A1->offm = A.offm;
  A1->offn = A.offn;
  A1->base = A.base;

  A2->m    = A.m;
  A2->n    = A.n - nb;
  A2->offm = A.offm;
  A2->offn = A.offn + nb;
  A2->base = A.base;

  return FLA_SUCCESS;
}


//
// --- FLA_Repart_2x2_to_3x3() -------------------------------------------------
//

FLA_Error FLA_Repart_2x2_to_3x3( FLA_Obj ATL, FLA_Obj ATR,  FLA_Obj *A00, FLA_Obj *A01, FLA_Obj *A02,
                                                            FLA_Obj *A10, FLA_Obj *A11, FLA_Obj *A12,
                                 FLA_Obj ABL, FLA_Obj ABR,  FLA_Obj *A20, FLA_Obj *A21, FLA_Obj *A22,
                                 dim_t   mb,  dim_t    nb,  FLA_Quadrant quadrant )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Repart_2x2_to_3x3_check( ATL, ATR,       A00, A01, A02,
                                                 A10, A11, A12,
                                 ABL, ABR,       A20, A21, A22,
                                 mb, nb, quadrant );

  if ( quadrant == FLA_BR )
  { 
    FLA_Part_2x2( ABR,    A11, A12,
                          A21, A22,    mb, nb, FLA_TL );
 
    FLA_Part_1x2( ATR,    A01, A02,    nb, FLA_LEFT );
 
    FLA_Part_2x1( ABL,    A10, 
                          A20,         mb, FLA_TOP );

    A00->m    = ATL.m;
    A00->n    = ATL.n;
    A00->offm = ATL.offm;
    A00->offn = ATL.offn;
    A00->base = ATL.base;
  }
  else if ( quadrant == FLA_BL )
  {
    FLA_Part_2x2( ABL,    A10, A11,
                          A20, A21,    mb, nb, FLA_TR );
 
    FLA_Part_1x2( ATL,    A00, A01,    nb, FLA_RIGHT );
 
    FLA_Part_2x1( ABR,    A12, 
                          A22,         mb, FLA_TOP );

    A02->m    = ATR.m;
    A02->n    = ATR.n;
    A02->offm = ATR.offm;
    A02->offn = ATR.offn;
    A02->base = ATR.base;
  }
  else if ( quadrant == FLA_TL )
  {
    FLA_Part_2x2( ATL,    A00, A01,
                          A10, A11,    mb, nb, FLA_BR );
 
    FLA_Part_1x2( ABL,    A20, A21,    nb, FLA_RIGHT );
 
    FLA_Part_2x1( ATR,    A02, 
                          A12,         mb, FLA_BOTTOM );

    A22->m    = ABR.m;
    A22->n    = ABR.n;
    A22->offm = ABR.offm;
    A22->offn = ABR.offn;
    A22->base = ABR.base;
  }
  else if ( quadrant == FLA_TR )
  {
    FLA_Part_2x2( ATR,    A01, A02,
                          A11, A12,    mb, nb, FLA_BL );

    FLA_Part_1x2( ABR,    A21, A22,    nb, FLA_LEFT );

    FLA_Part_2x1( ATL,    A00,
                          A10,         mb, FLA_BOTTOM );

    A20->m    = ABL.m;
    A20->n    = ABL.n;
    A20->offm = ABL.offm;
    A20->offn = ABL.offn;
    A20->base = ABL.base;
  }

  return FLA_SUCCESS;
}


//
// --- FLA_Repart_2x1_to_3x1() -------------------------------------------------
//

FLA_Error FLA_Repart_2x1_to_3x1( FLA_Obj AT,  FLA_Obj *A0,
                                              FLA_Obj *A1,
                                 FLA_Obj AB,  FLA_Obj *A2,
                                 dim_t   mb,  FLA_Side side )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Repart_2x1_to_3x1_check( AT,     A0, 
                                         A1, 
                                 AB,     A2,     mb, side );

  if ( side == FLA_TOP )
  {
    FLA_Part_2x1 ( AT,    A0, 
                          A1,    mb, FLA_BOTTOM );
 
    A2->m    = AB.m;
    A2->n    = AB.n;
    A2->offm = AB.offm;
    A2->offn = AB.offn;
    A2->base = AB.base;
  }
  else
  {
    A0->m    = AT.m;
    A0->n    = AT.n;
    A0->offm = AT.offm;
    A0->offn = AT.offn;
    A0->base = AT.base;

    FLA_Part_2x1 ( AB,    A1, 
                          A2,    mb, FLA_TOP );
  }

  return FLA_SUCCESS;
}


//
// --- FLA_Repart_1x2_to_1x3() -------------------------------------------------
//

FLA_Error FLA_Repart_1x2_to_1x3( FLA_Obj  AL,              FLA_Obj  AR,
                                 FLA_Obj *A0, FLA_Obj *A1, FLA_Obj *A2,
                                              dim_t    nb, FLA_Side side )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Repart_1x2_to_1x3_check( AL, AR,        A0, A1, A2,
                                 nb, side );

  if ( side == FLA_LEFT )
  {
    FLA_Part_1x2( AL,    A0, A1,    nb, FLA_RIGHT );
 
    A2->m    = AR.m;
    A2->n    = AR.n;
    A2->offm = AR.offm;
    A2->offn = AR.offn;
    A2->base = AR.base;
  }
  else
  {
    A0->m    = AL.m;
    A0->n    = AL.n;
    A0->offm = AL.offm;
    A0->offn = AL.offn;
    A0->base = AL.base;

    FLA_Part_1x2( AR,    A1, A2,    nb, FLA_LEFT );
  }

  return FLA_SUCCESS;
}


//
// --- FLA_Cont_with_3x3_to_2x2() ----------------------------------------------
//

FLA_Error FLA_Cont_with_3x3_to_2x2( FLA_Obj *ATL, FLA_Obj *ATR,  FLA_Obj A00, FLA_Obj A01, FLA_Obj A02,
                                                                 FLA_Obj A10, FLA_Obj A11, FLA_Obj A12,
                                    FLA_Obj *ABL, FLA_Obj *ABR,  FLA_Obj A20, FLA_Obj A21, FLA_Obj A22,
                                                                 FLA_Quadrant quadrant )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Cont_with_3x3_to_2x2_check( ATL, ATR,       A00, A01, A02,
                                                    A10, A11, A12,
                                    ABL, ABR,       A20, A21, A22,
                                    quadrant );

  if ( quadrant == FLA_TL )
  {
    ATL->m    = A00.m + A10.m;
    ATL->n    = A00.n + A01.n;
    ATL->offm = A00.offm;
    ATL->offn = A00.offn;
    ATL->base = A00.base;
 
    ATR->m    = A02.m + A12.m;
    ATR->n    = A02.n;
    ATR->offm = A02.offm;
    ATR->offn = A02.offn;
    ATR->base = A02.base;
 
    ABL->m    = A20.m;
    ABL->n    = A20.n + A21.n;
    ABL->offm = A20.offm;
    ABL->offn = A20.offn;
    ABL->base = A20.base;
 
    ABR->m    = A22.m;
    ABR->n    = A22.n;
    ABR->offm = A22.offm;
    ABR->offn = A22.offn;
    ABR->base = A22.base;
  }
  else if ( quadrant == FLA_TR )
  {
    ATL->m    = A00.m + A10.m;
    ATL->n    = A00.n;
    ATL->offm = A00.offm;
    ATL->offn = A00.offn;
    ATL->base = A00.base;
 
    ATR->m    = A01.m + A11.m;
    ATR->n    = A01.n + A02.n;
    ATR->offm = A01.offm;
    ATR->offn = A01.offn;
    ATR->base = A01.base;
 
    ABL->m    = A20.m;
    ABL->n    = A20.n;
    ABL->offm = A20.offm;
    ABL->offn = A20.offn;
    ABL->base = A20.base;
 
    ABR->m    = A21.m;
    ABR->n    = A21.n + A22.n;
    ABR->offm = A21.offm;
    ABR->offn = A21.offn;
    ABR->base = A21.base;
  }
  else if ( quadrant == FLA_BR )
  {
    ATL->m    = A00.m;
    ATL->n    = A00.n;
    ATL->offm = A00.offm;
    ATL->offn = A00.offn;
    ATL->base = A00.base;
 
    ATR->m    = A01.m;
    ATR->n    = A01.n + A02.n;
    ATR->offm = A01.offm;
    ATR->offn = A01.offn;
    ATR->base = A01.base;
 
    ABL->m    = A10.m + A20.m;
    ABL->n    = A10.n;
    ABL->offm = A10.offm;
    ABL->offn = A10.offn;
    ABL->base = A10.base;
 
    ABR->m    = A11.m + A21.m;
    ABR->n    = A11.n + A12.n;
    ABR->offm = A11.offm;
    ABR->offn = A11.offn;
    ABR->base = A11.base;
  }
  else if ( quadrant == FLA_BL )
  {
    ATL->m    = A00.m;
    ATL->n    = A00.n + A01.n;
    ATL->offm = A00.offm;
    ATL->offn = A00.offn;
    ATL->base = A00.base;
 
    ATR->m    = A02.m;
    ATR->n    = A02.n;
    ATR->offm = A02.offm;
    ATR->offn = A02.offn;
    ATR->base = A02.base;
 
    ABL->m    = A10.m + A20.m;
    ABL->n    = A10.n + A11.n;
    ABL->offm = A10.offm;
    ABL->offn = A10.offn;
    ABL->base = A10.base;
 
    ABR->m    = A12.m + A22.m;
    ABR->n    = A12.n ;
    ABR->offm = A12.offm;
    ABR->offn = A12.offn;
    ABR->base = A12.base;
  }

  return FLA_SUCCESS;
}


//
// --- FLA_Cont_with_3x1_to_2x1() ----------------------------------------------
//

FLA_Error FLA_Cont_with_3x1_to_2x1( FLA_Obj *AT,  FLA_Obj A0,
                                                  FLA_Obj A1,
                                    FLA_Obj *AB,  FLA_Obj A2,
                                                  FLA_Side side )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Cont_with_3x1_to_2x1_check( AT,     A0, 
                                            A1, 
                                    AB,     A2,     side );

  if ( side == FLA_TOP )
  {
    AT->m    = A0.m + A1.m;
    AT->n    = A0.n;
    AT->offm = A0.offm;
    AT->offn = A0.offn;
    AT->base = A0.base;

    AB->m    = A2.m;
    AB->n    = A2.n;
    AB->offm = A2.offm;
    AB->offn = A2.offn;
    AB->base = A2.base;
  }
  else
  {
    AT->m    = A0.m;
    AT->n    = A0.n;
    AT->offm = A0.offm;
    AT->offn = A0.offn;
    AT->base = A0.base;
 
    AB->m    = A1.m + A2.m;
    AB->n    = A1.n;
    AB->offm = A1.offm;
    AB->offn = A1.offn;
    AB->base = A1.base;
  }

  return FLA_SUCCESS;
}


//
// --- FLA_Cont_with_1x3_to_1x2() ----------------------------------------------
//

FLA_Error FLA_Cont_with_1x3_to_1x2( FLA_Obj *AL,              FLA_Obj *AR,
                                    FLA_Obj  A0, FLA_Obj  A1, FLA_Obj  A2,
                                                              FLA_Side side )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Cont_with_1x3_to_1x2_check( AL, AR,        A0, A1, A2,
                                    side );

  if ( side == FLA_LEFT )
  {
    AL->m    = A0.m;
    AL->n    = A0.n + A1.n;
    AL->offm = A0.offm;
    AL->offn = A0.offn;
    AL->base = A0.base;

    AR->m    = A2.m;
    AR->n    = A2.n;
    AR->offm = A2.offm;
    AR->offn = A2.offn;
    AR->base = A2.base;
  }
  else
  {
    AL->m    = A0.m;
    AL->n    = A0.n;
    AL->offm = A0.offm;
    AL->offn = A0.offn;
    AL->base = A0.base;

    AR->m    = A1.m;
    AR->n    = A1.n + A2.n;
    AR->offm = A1.offm;
    AR->offn = A1.offn;
    AR->base = A1.base;
  }

  return FLA_SUCCESS;
}


//
// --- FLA_Merge_2x2() ---------------------------------------------------------
//

FLA_Error FLA_Merge_2x2( FLA_Obj A11, FLA_Obj A12,
                         FLA_Obj A21, FLA_Obj A22,  FLA_Obj *A )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Merge_2x2_check( A11, A12,
                         A21, A22,    A );

  A->m      = A11.m + A21.m;
  A->n      = A11.n + A12.n;
  A->offm   = A11.offm;
  A->offn   = A11.offn;
  A->base   = A11.base;

  return FLA_SUCCESS;
}


//
// --- FLA_Merge_2x1() ---------------------------------------------------------
//

FLA_Error FLA_Merge_2x1( FLA_Obj AT,
                         FLA_Obj AB,  FLA_Obj *A )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Merge_2x1_check( AT,
                         AB,   A );

  A->m    = AT.m + AB.m;
  A->n    = AT.n;
  A->offm = AT.offm;
  A->offn = AT.offn;
  A->base = AT.base;

  return FLA_SUCCESS;
}


//
// --- FLA_Merge_1x2() ---------------------------------------------------------
//

FLA_Error FLA_Merge_1x2( FLA_Obj AL, FLA_Obj AR,   FLA_Obj *A )
{
  if ( FLA_Check_error_level() == FLA_FULL_ERROR_CHECKING )
    FLA_Merge_1x2_check( AL, AR,    A );

  A->m    = AL.m;
  A->n    = AL.n + AR.n;
  A->offm = AL.offm;
  A->offn = AL.offn;
  A->base = AL.base;

  return FLA_SUCCESS;
}

