/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/
#include <math.h>
#include "FLAME.h"
#include "Chol_prototypes.h"

#define AA(i,j) buff_A[ (j)*ldim_A + (i) ]

FLA_Error REF_Chol( FLA_Obj A )
{
  int n, i, j, k, ldim_A;

  double *buff_A;

  n = FLA_Obj_length( A );

  ldim_A = FLA_Obj_col_stride( A );

  buff_A = (double *) FLA_Obj_buffer_at_view( A );

  for ( j=0; j<n; j++ ){
    AA( j,j ) = sqrt( AA( j,j ) );
    /* a21 = a21 / alpha11 */
    for ( i=j+1; i<n; i++ )
      AA( i,j ) = AA( i,j )  / AA( j,j );
    /* A22 = A22 - a21 * a21'  (lower triangular part only) */
    for ( k=j+1; k<n; k++ )
      for ( i=k; i<n; i++ )
        AA( i,k ) = AA( i,k ) - AA( i,j ) * AA( j,k );
  }
}

