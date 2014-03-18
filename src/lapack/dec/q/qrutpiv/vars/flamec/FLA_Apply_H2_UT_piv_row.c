
#include "FLAME.h"


FLA_Error FLA_Apply_H2_UT_piv_row( FLA_Obj tau, FLA_Obj a1t, FLA_Obj u1t, FLA_Obj W,
                                   FLA_Obj u2,  FLA_Obj A2,  FLA_Obj U2,  FLA_Obj w1t,
                                   FLA_Obj vt )
/*
  Apply a single Householder transform H' from the left to a row vector a1t and
  a matrix A2T, but the transformation is not applied to A2B (delayed) :

    / a1t \  := H / a1t \                                          (1)
    \ A2  /       \ A2  /


  H is defined as:

    H  =  / I - inv(tau) / 1  \ ( 1  u2') \                            (2)
          \              \ u2 /           /

  a1t is updated as:

      a1t    :=    a1t   -       inv(tau) ( a1t  +  u2' A2 )

  where a1t and A2 on the right hand side are explicitly/implicitly defined as:
 
      a1t    :=    a1t   -   u1t W

      A2     :=    A2    -   U2  W

  Then,

      a1t    :=    a1t   -       inv(tau) ( a1t  +  u2' A2 - u2' U2 W )

  w1t is stored for the next update:

      w1t := inv(tau) ( a1t  +  u2'  A2 )

  -FGVZ
*/
{
  // a1t -= u1t W = 1 a1t -1 W^T u1t;
  FLA_Gemvc_external( FLA_TRANSPOSE, FLA_NO_CONJUGATE, 
                      FLA_MINUS_ONE, W, u1t, FLA_ONE, a1t );

  // w1t := a1t;
  FLA_Copy_external( a1t, w1t );
  
  // w1t += u2' A2 = 1 w1t + 1 A2^T conj(u2);
  FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, FLA_ONE, w1t );

  if ( FLA_Obj_min_dim( U2 ) > 0 ) 
  {
    FLA_Obj vtR;

    // Partition the workspace (a row vector matching the width of a1t)
    FLA_Part_1x2( vt, &vt, &vtR, FLA_Obj_width( U2 ), FLA_LEFT );

    // vt := u2'U2 = 0 vt + 1 U2^T conj(u2);
    FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, 
                        FLA_ONE, U2, u2, FLA_ZERO, vt );

    // w1t -= - vt W = 1 w1t -1 W^T vt;
    FLA_Gemvc_external( FLA_TRANSPOSE, FLA_NO_CONJUGATE, 
                        FLA_MINUS_ONE, W, vt, FLA_ONE, w1t );
  }

  // w1t = w1t / tau;
  FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );

  // a1t  = a1t - w1t;
  FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );

  return FLA_SUCCESS;
}
