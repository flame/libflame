/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#include "blis1.h"

void bl1_screate_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, float* a_save, integer a_rs_save, integer a_cs_save, float** a, integer* a_rs, integer* a_cs )
{
	integer dim_a;

	// Choose the dimension of the matrix based on the side parameter.
	if ( bl1_is_left( side ) ) dim_a = m;
	else                       dim_a = n;

	// Call the simple version with chosen dimensions.
	bl1_screate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      a,      a_rs,      a_cs );
}

void bl1_dcreate_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, double* a_save, integer a_rs_save, integer a_cs_save, double** a, integer* a_rs, integer* a_cs )
{
	integer dim_a;

	// Choose the dimension of the matrix based on the side parameter.
	if ( bl1_is_left( side ) ) dim_a = m;
	else                       dim_a = n;

	// Call the simple version with chosen dimensions.
	bl1_dcreate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      a,      a_rs,      a_cs );
}

void bl1_ccreate_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, scomplex* a_save, integer a_rs_save, integer a_cs_save, scomplex** a, integer* a_rs, integer* a_cs )
{
	integer dim_a;

	// Choose the dimension of the matrix based on the side parameter.
	if ( bl1_is_left( side ) ) dim_a = m;
	else                       dim_a = n;

	// Call the simple version with chosen dimensions.
	bl1_ccreate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      a,      a_rs,      a_cs );
}

void bl1_zcreate_contigmsr( side1_t side, uplo1_t uplo, integer m, integer n, dcomplex* a_save, integer a_rs_save, integer a_cs_save, dcomplex** a, integer* a_rs, integer* a_cs )
{
	integer dim_a;

	// Choose the dimension of the matrix based on the side parameter.
	if ( bl1_is_left( side ) ) dim_a = m;
	else                       dim_a = n;

	// Call the simple version with chosen dimensions.
	bl1_zcreate_contigmr( uplo,
	                      dim_a,
	                      dim_a,
	                      a_save, a_rs_save, a_cs_save,
	                      a,      a_rs,      a_cs );
}

