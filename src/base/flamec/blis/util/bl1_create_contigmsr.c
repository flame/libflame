
#include "blis1.h"

void bl1_screate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, float* a_save, int a_rs_save, int a_cs_save, float** a, int* a_rs, int* a_cs )
{
	int dim_a;

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

void bl1_dcreate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, double* a_save, int a_rs_save, int a_cs_save, double** a, int* a_rs, int* a_cs )
{
	int dim_a;

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

void bl1_ccreate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, scomplex* a_save, int a_rs_save, int a_cs_save, scomplex** a, int* a_rs, int* a_cs )
{
	int dim_a;

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

void bl1_zcreate_contigmsr( side1_t side, uplo1_t uplo, int m, int n, dcomplex* a_save, int a_rs_save, int a_cs_save, dcomplex** a, int* a_rs, int* a_cs )
{
	int dim_a;

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

