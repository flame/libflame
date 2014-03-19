/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

conj1_t bl1_proj_trans1_to_conj( trans1_t trans )
{
	if ( bl1_does_conj( trans ) ) return BLIS1_CONJUGATE;
	else                          return BLIS1_NO_CONJUGATE;
}

