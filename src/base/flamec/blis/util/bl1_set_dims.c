/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

/*
*     Modifications Copyright (c) 2023 Advanced Micro Devices, Inc.  All rights reserved.
*/
#include "blis1.h"
#if FLA_ENABLE_AOCL_BLAS
#include "blis.h"
#endif

void bl1_set_dims_with_trans( trans1_t trans, integer m, integer n, integer* m_new, integer* n_new )
{
	if ( bl1_does_trans( trans ) )
	{
		*m_new = n;
		*n_new = m;
	}
	else
	{
		*m_new = m;
		*n_new = n;
	}
}

void bl1_set_dim_with_side( side1_t side, integer m, integer n, integer* dim_new )
{
	if ( bl1_is_left( side ) )
	{
		*dim_new = m;
	}
	else // if ( bl1_is_right( side ) )
	{
		*dim_new = n;
	}
}

