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

int bl1_does_trans( trans1_t trans )
{
	return ( trans == BLIS1_TRANSPOSE ||
	         trans == BLIS1_CONJ_TRANSPOSE );
}

int bl1_does_notrans( trans1_t trans )
{
	return ( trans == BLIS1_NO_TRANSPOSE ||
	         trans == BLIS1_CONJ_NO_TRANSPOSE );
}

int bl1_does_conj( trans1_t trans )
{
	return ( trans == BLIS1_CONJ_NO_TRANSPOSE ||
	         trans == BLIS1_CONJ_TRANSPOSE );
}

