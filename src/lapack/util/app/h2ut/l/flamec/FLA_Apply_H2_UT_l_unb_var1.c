
#include "FLAME.h"

FLA_Error FLA_Apply_H2_UT_l_unb_var1( FLA_Obj tau, FLA_Obj u2, FLA_Obj a1t,
                                                               FLA_Obj A2 )
/*
  Apply a single Householder transform H' from the left to a row vector a1t and
  a matrix A2:

    / a1t \  := H / a1t \                                                  (1)
    \ A2  /       \ A2  /

  H is defined as:

    H  =  / I - inv(tau) / 1  \ ( 1  u2' ) \                                (2)
          \              \ u2 /            /

  Typically, a1t and A2 are vertically adjacent views into a larger matrix,
  though this is not always the case.

  Substituting (2) into (1), we have:

    / a1t \  := H / a1t \
    \ A2  /       \ A2  /
              =  / I - inv(tau) / 1  \ ( 1  u2' ) \ / a1t \
                 \              \ u2 /            / \ A2  /
              =  / I - inv(tau) / 1     u2'  \ \ / a1t \
                 \              \ u2  u2 u2' / / \ A2  /
              =  / a1t \ - inv(tau) / 1     u2'  \ / a1t \
                 \ A2  /            \ u2  u2 u2' / \ A2  /
              =  / a1t \ - /   inv(tau)      inv(tau) u2'  \ / a1t \
                 \ A2  /   \ inv(tau) u2   inv(tau) u2 u2' / \ A2  /
              =  / a1t \ - /   inv(tau) a1t  +   inv(tau) u2' A2   \
                 \ A2  /   \ inv(tau) u2 a1t +  inv(tau) u2 u2' A2 /

  Thus, a1t is updated as:

      a1t    :=    a1t  -  inv(tau) a1t  +  inv(tau) u2' A2
              =    a1t  -  inv(tau) ( a1t  +  u2' A2 )

  And A2 is updated as:

      A2     :=    A2   -  inv(tau) u2 a1t  +  inv(tau) u2 u2' A2
              =    A2   -  u2 ( inv(tau) a1t  +  inv(tau) u2' A2 )
              =    A2   -  u2 inv(tau) ( a1t  +  u2' A2 )

  Note that:

    inv(tau) ( a1t  +  u2' A2 )

  is common to both updates, and thus may be computed and stored in
  workspace, and then re-used.
 
  -FGVZ
*/
{
  FLA_Obj w1t;

  if ( FLA_Obj_has_zero_dim( a1t ) || FLA_Obj_equals( tau, FLA_ZERO ) ) return FLA_SUCCESS;

  // w1t = a1t;
  FLA_Obj_create_conf_to( FLA_NO_TRANSPOSE, a1t, &w1t );
  FLA_Copy_external( a1t, w1t );

  // w1t = w1t + u2' * A2;
  // w1t = w1t + A2^T * conj(u2);
  FLA_Gemvc_external( FLA_TRANSPOSE, FLA_CONJUGATE, FLA_ONE, A2, u2, FLA_ONE, w1t );

  // w1t = w1t / tau;
  FLA_Inv_scalc_external( FLA_NO_CONJUGATE, tau, w1t );

  // a1t = a1t - w1t;
  FLA_Axpy_external( FLA_MINUS_ONE, w1t, a1t );

  // A2 = A2 - u2 * w1t;
  FLA_Ger_external( FLA_MINUS_ONE, u2, w1t, A2 );

  FLA_Obj_free( &w1t );

  return FLA_SUCCESS;
}
