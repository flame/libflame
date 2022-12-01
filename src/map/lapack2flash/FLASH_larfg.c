/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"


// Our implementation is fundamentally different from LAPACK implementation
#ifdef FLA_ENABLE_YET_LAPACK2FLAME

#include "FLASH_lapack2flash_util_defs.h"

/*
  LARFG generates a real elementary reflector H of order n,
  such that

  H * ( chi ) = ( alpha ),   H**T * H = I.
      ( x2  )   (   0   )

  where chi and alpha are scalars, and x2 is an (n-1)-element
  vector. H is represented in the form

  H = I - tau * ( 1  ) * ( 1 u2**T ) ,
                ( u2 )
  where tau is a real scalar and u2 is a real (n-1)-element vector.
*/
#define LAPACK_larfg(prefix)                                            \
  int F77_ ## prefix ## larfg( int *n,                                  \
                               PREFIX2LAPACK_TYPEDEF(prefix)* chi,      \
                               PREFIX2LAPACK_TYPEDEF(prefix)* x2,       \
                               int* inc_x2,                             \
                               PREFIX2LAPACK_TYPEDEF(prefix)* tau )

#define LAPACK_larfg_body(prefix)                                       \
  FLA_Error    init_result;                                             \
  FLA_Init_safe( &init_result );                                        \
  if ( *n > 0 ) {                                                       \
    FLA_Househ2_UT_l_op ## prefix ( *n - 1,                             \
                                    chi,                                \
                                    x2,                                 \
                                    *inc_x2,                            \
                                    tau );                              \
      /* TODO: invert tau */
}                                                                     \
FLA_Finalize_safe( init_result );
\
return 0;


       LAPACK_larfg(s)
{
    LAPACK_larfg_body(s)
}
LAPACK_larfg(d)
{
    LAPACK_larfg_body(d)
}
LAPACK_larfg(c)
{
    LAPACK_larfg_body(c)
}
LAPACK_larfg(z)
{
    LAPACK_larfg_body(z)
}

#define LAPACK_larfgp(prefix)                                           \
  int F77_ ## prefix ## larfgp( int *n,                                 \
                                PREFIX2LAPACK_TYPEDEF(prefix)* chi,     \
                                PREFIX2LAPACK_TYPEDEF(prefix)* x2,      \
                                int* inc_x2,                            \
                                PREFIX2LAPACK_TYPEDEF(prefix)* tau )

LAPACK_larfgp(s)
{
    LAPACK_larfg_body(s)
}
LAPACK_larfgp(d)
{
    LAPACK_larfg_body(d)
}
LAPACK_larfgp(c)
{
    LAPACK_larfg_body(c)
}
LAPACK_larfgp(z)
{
    LAPACK_larfg_body(z)
}


#endif
