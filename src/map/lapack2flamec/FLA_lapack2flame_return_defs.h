/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

#ifndef FLA_LAPACK2FLAME_RETURN_DEFS_H
#define FLA_LAPACK2FLAME_RETURN_DEFS_H

// --- LAPACK return values ----

#define LAPACK_SUCCESS 512
#define LAPACK_FAILURE 312
#define LAPACK_QUICK_RETURN 212
#define LAPACK_QUERY_RETURN 112

#define LAPACK_RETURN_CHECK( r_check )                                  \
  {                                                                     \
  int r_val = r_check;                                                  \
  switch ( r_val )                                                      \
    {                                                                   \
    case LAPACK_FAILURE:      return FLA_FAILURE;                       \
    case LAPACK_QUERY_RETURN: ;                                         \
    case LAPACK_QUICK_RETURN: return 0;                                 \
    case LAPACK_SUCCESS: ;                                              \
    default: ;                                                          \
      if ( r_val > 0 ) { ; }                                            \
      else             { FLA_Check_error_code( FLA_LAPAC2FLAME_INVALID_RETURN ); } \
    }                                                                   \
  }

extern int lsame_(char *, char *);
extern int xerbla_(char *, int *);
extern int ilaenv_(int *, char *, char *, int *, int *, int *, int *);


#endif
