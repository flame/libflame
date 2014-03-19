/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "FLAME.h"

FLA_Error FLA_Shift_pivots_to( FLA_Pivot_type ptype, FLA_Obj p )
{
  int  m_p, n_p;
  int* buff_p;
  int  i;

  if ( FLA_Check_error_level() >= FLA_MIN_ERROR_CHECKING )
    FLA_Shift_pivots_to_check( ptype, p );

  m_p    = FLA_Obj_length( p );
  n_p    = FLA_Obj_width( p );
  buff_p = FLA_INT_PTR( p );

  if ( m_p < 1 || n_p < 1 ) return FLA_SUCCESS;

  if ( ptype == FLA_LAPACK_PIVOTS )
  {
    // Shift FLAME pivots to LAPACK pivots.
    for ( i = 0; i < m_p; i++ )
      buff_p[ i ] += i + 1;
  }
  else
  {
    // Otherwise, shift LAPACK pivots back to FLAME.
    for ( i = 0; i < m_p; i++ )
      buff_p[ i ] -= i + 1;
  }

  return FLA_SUCCESS;
}

