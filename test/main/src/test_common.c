/*
	Copyright (c) 2022 Advanced Micro Devices, Inc.  All rights reserved.
*/

#include "test_libflame.h"
#include "test_common.h"


/* create matrix of given datatype*/
void create_matrix(integer datatype, integer m, integer n, void **A)
{
	*A = NULL;

	switch(datatype)
	{
		case INT:
			*A = (integer *)malloc(m * n * sizeof(integer));
			break;

		case FLOAT:
			*A = (float *)malloc(m * n * sizeof(float));
			break;

		case DOUBLE:
			*A = (double *)malloc(m * n * sizeof(double));
			break;

		case COMPLEX:
			*A = (scomplex *)malloc(m * n * sizeof(scomplex));
			break;

		case DOUBLE_COMPLEX:
			*A = (dcomplex *)malloc(m * n * sizeof(dcomplex));
			break;
	}

	if(*A == NULL)
	{
		fprintf( stderr, "malloc() returned NULL pointer\n");
		abort();
	}

	return;
}


void create_realtype_matrix(integer datatype, integer m, integer n, void **A)
{
	*A = NULL;

	if(datatype == FLOAT || datatype == COMPLEX)
		*A = (float *)malloc(m * n * sizeof(float));
	else
		*A = (double *)malloc(m * n * sizeof(double));

	if(*A == NULL)
	{
		fprintf( stderr, "malloc() returned NULL pointer\n");
		abort();
	}

	return;
}


/* free matrix */
void free_matrix(void *A)
{

	if(!A)
		return;

	free(A);

}


/* Initialize matrix with random values */
void rand_matrix(integer datatype, integer m, integer n, void *A)
{

	switch( datatype )
	{
		case FLOAT:
		{
			float *buff_A = ( float * ) A;
			rand_matrix_s( buff_A, m, n, m );
			break;
		}

		case DOUBLE:
		{
			double *buff_A = ( double * ) A;
			rand_matrix_d( buff_A, m, n, m );
			break;
		}

		case COMPLEX:
		{
			scomplex *buff_A = ( scomplex * ) A;
			rand_matrix_c( buff_A, m, n, m );
			break;
		}

		case DOUBLE_COMPLEX:
		{
			dcomplex *buff_A = ( dcomplex * ) A;
			rand_matrix_z( buff_A, m, n, m );
			break;
		}
	}

	return;
}


void rand_matrix_s( float * A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[i * LDA + j] = SRAND();
		}
	}

	return;
}

void rand_matrix_d( double * A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[i * LDA + j] = DRAND();
		}
	}

	return;
}

void rand_matrix_c( scomplex * A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[i * LDA + j].real = SRAND();
			A[i * LDA + j].imag = SRAND();
		}
	}

	return;
}

void rand_matrix_z( dcomplex * A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[i * LDA + j].real = DRAND();
			A[i * LDA + j].imag = DRAND();
		}
	}

	return;
}


/* Initialize symmetric matrix with random values */
void rand_sym_matrix_s( float *A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < M; j++ )
		{
			A[i * LDA + j] = SRAND();
		}
	}

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < M; j++ )
		{
			A[j * LDA + i] = A[i * LDA + j];
		}
	}

	return;
}

void rand_sym_matrix_d( double * A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < M; j++ )
		{
			A[i * LDA + j] = DRAND();
		}
	}

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < M; j++ )
		{
			A[j * LDA + i] = A[i * LDA + j];
		}
	}

	return;
}

void rand_sym_matrix_c( scomplex * A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < M; j++ )
		{
			A[i * LDA + j].real = SRAND();
			A[i * LDA + j].imag = SRAND();
		}
	}

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < M; j++ )
		{
			A[j * LDA + i].real = A[i * LDA + j].real;
			A[j * LDA + i].imag = A[i * LDA + j].imag;
		}
	}

	return;
}

void rand_sym_matrix_z( dcomplex * A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < M; j++ )
		{
			A[i * LDA + j].real = DRAND();
			A[i * LDA + j].imag = DRAND();
		}
	}

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < M; j++ )
		{
			A[j * LDA + i].real = A[i * LDA + j].real;
			A[j * LDA + i].imag = A[i * LDA + j].imag;
		}
	}

	return;
}


/* Copy a matrix */
void copy_matrix(integer datatype, integer m, integer n, void *A, void *B)
{

	switch( datatype )
	{
		case INT:
		{
			integer *buff_A = (integer *)A;
			integer *buff_B = (integer *)B;

			copy_matrix_i( buff_A, buff_B, m, n, m, n);

			break;
		}

		case FLOAT:
		{
			float *buff_A = (float *)A;
			float *buff_B = (float *)B;

			copy_matrix_s( buff_A, buff_B, m, n, m, n);

			break;
		}

		case DOUBLE:
		{
			double *buff_A = (double *)A;
			double *buff_B = (double *)B;

			copy_matrix_d( buff_A, buff_B, m, n, m, n);

			break;
		}

		case COMPLEX:
		{
			scomplex *buff_A = (scomplex *)A;
			scomplex *buff_B = (scomplex *)B;

			copy_matrix_c( buff_A, buff_B, m, n, m, n);

			break;
		}

		case DOUBLE_COMPLEX:
		{

			dcomplex *buff_A = (dcomplex *)A;
			dcomplex *buff_B = (dcomplex *)B;

			copy_matrix_z( buff_A, buff_B, m, n, m, n);

			break;
		}
	}

	return;
}


void copy_realtype_matrix(integer datatype, integer m, integer n, void *A, void *B)
{

	if(datatype == FLOAT || datatype == COMPLEX)
	{
		float *buff_A = (float *)A;
		float *buff_B = (float *)B;

		copy_matrix_s( buff_A, buff_B, m, n, m, n);
	}
	else
	{
		double *buff_A = (double *)A;
		double *buff_B = (double *)B;

		copy_matrix_d( buff_A, buff_B, m, n, m, n);
	}

	return;
}


void copy_matrix_i( integer *sM, integer *dM, integer M, integer N, integer LDS, integer LDD )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			dM[ i * LDD + j ] = sM[ i * LDS + j ];
		}
	}

	return;
}

void copy_matrix_s( float *sM, float *dM, integer M, integer N, integer LDS, integer LDD )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			dM[ i * LDD + j ] = sM[ i * LDS + j ];
		}
	}

	return;
}

void copy_matrix_d( double *sM, double *dM, integer M, integer N, integer LDS, integer LDD )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			dM[ i * LDD + j ] = sM[ i * LDS + j ];
		}
	}

	return;
}

void copy_matrix_c( scomplex *sM, scomplex *dM, integer M, integer N, integer LDS, integer LDD )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			dM[ i * LDD + j ].real = sM[ i * LDS + j ].real;
			dM[ i * LDD + j ].imag = sM[ i * LDS + j ].imag;
		}
	}

	return;
}

void copy_matrix_z( dcomplex *sM, dcomplex *dM, integer M, integer N, integer LDS, integer LDD )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			dM[ i * LDD + j ].real = sM[ i * LDS + j ].real;
			dM[ i * LDD + j ].imag = sM[ i * LDS + j ].imag;
		}
	}

	return;
}


/* Pack a symmetric matrix in column first order */
void pack_matrix_lt_s( float *A, float *B, integer N, integer LDA )
{
	integer i, j;
	float *bptr = B;

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < N; j++ )
		{
			*bptr++ = A[ i * LDA + j ];
		}
	}

	return;
}

void pack_matrix_lt_d( double *A, double *B, integer N, integer LDA )
{
	integer i, j;
	double *bptr = B;

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < N; j++ )
		{
			*bptr++ = A[ i * LDA + j ];
		}
	}

	return;
}

void pack_matrix_lt_c( scomplex *A, scomplex *B, integer N, integer LDA )
{
	integer i, j;
	scomplex *bptr = B;

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < N; j++ )
		{
			bptr->real = A[ i * LDA + j ].real;
			bptr->imag = A[ i * LDA + j ].imag;
			bptr++;
		}
	}

	return;
}

void pack_matrix_lt_z( dcomplex *A, dcomplex *B, integer N, integer LDA )
{
	integer i, j;
	dcomplex *bptr = B;

	for( i = 0; i < N; i++ )
	{
		for( j = i; j < N; j++ )
		{
			bptr->real = A[ i * LDA + j ].real;
			bptr->imag = A[ i * LDA + j ].imag;
			bptr++;
		}
	}

	return;
}


/* Initialize a matrix with zeros */
void reset_matrix_s( float *A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[ i * LDA + j ] = 0.f;
		}
	}

	return;
}

void reset_matrix_d( double *A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[ i * LDA + j ] = 0.;
		}
	}

	return;
}

void reset_matrix_c( scomplex *A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[ i * LDA + j ].real = 0.f;
			A[ i * LDA + j ].imag = 0.f;
		}
	}

	return;
}

void reset_matrix_z( dcomplex *A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[ i * LDA + j ].real = 0.;
			A[ i * LDA + j ].imag = 0.;
		}
	}

	return;
}


/* Set a matrix to identity */
void set_identity_s( float *A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[ i * LDA + j ] = 0.f;
		}
		A[ i * LDA + i ] = 1.0f;
	}

	return;
}

void set_identity_d( double *A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[ i * LDA + j ] = 0.;
		}
		A[ i * LDA + i ] = 1.0;
	}

	return;
}

void set_identity_c( scomplex *A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[ i * LDA + j ].real = 0.f;
			A[ i * LDA + j ].imag = 0.f;
		}
		A[ i * LDA + i ].real = 1.0f;
		A[ i * LDA + i ].imag = 0.0f;
	}

	return;
}

void set_identity_z( dcomplex *A, integer M, integer N, integer LDA )
{
	integer i, j;

	for( i = 0; i < N; i++ )
	{
		for( j = 0; j < M; j++ )
		{
			A[ i * LDA + j ].real = 0.;
			A[ i * LDA + j ].imag = 0.;
		}
		A[ i * LDA + i ].real = 1.0;
		A[ i * LDA + i ].imag = 0.0;
	}

	return;
}

void z_div_t(dcomplex *cp, dcomplex *ap, dcomplex *bp)
{
	dcomplex a = *ap;
	dcomplex b = *bp;
	double temp;

	temp = b.real * b.real + b.imag * b.imag;
	if(!temp)
	{
		fprintf( stderr, "z_div_t : temp is zero. Abort\n");
		abort();
	}

	cp->real = ( a.real * b.real + a.imag * b.imag ) / temp;
	cp->imag = ( a.imag * b.real - a.real * b.imag ) / temp;
}


/* Division of complex types */
void c_div_t(scomplex *cp, scomplex *ap, scomplex *bp)
{
	scomplex a = *ap;
	scomplex b = *bp;
	float temp;

	temp = b.real * b.real + b.imag * b.imag;
	if(!temp)
	{
		fprintf( stderr, "z_div_t : temp is zero. Abort\n");
		abort();
	}

	cp->real = ( a.real * b.real + a.imag * b.imag ) / temp;
	cp->imag = ( a.imag * b.real - a.real * b.imag ) / temp;
}

/* work value calculation */
integer get_work_value( integer datatype, void *work )
{

	if(!work)
		return 0;

	if( datatype == FLOAT || datatype == COMPLEX )
		return (integer) (*(float*)work);
	else
		return (integer) (*(double*)work);
}


void diagmv( integer datatype, integer m, integer n, void* x, integer incx, void* a, integer a_rs, integer a_cs )
{
	integer inca, lda;
	integer n_iter;
	integer n_elem;
	integer j;

	if(m == 0 || n == 0)
		return;

	// Initialize with optimal values for column-major storage.
	inca   = a_rs;
	lda    = a_cs;
	n_iter = n;
	n_elem = m;

	switch(datatype)
	{
		case FLOAT:
		{
			float *a_begin;
			for ( j = 0; j < n_iter; j++ )
			{
				a_begin = (float *)a + j*lda;
				scalv( datatype, n_elem, x, incx, a_begin, inca );
			}
			break;
		}

		case DOUBLE:
		{
			double *a_begin;
			for ( j = 0; j < n_iter; j++ )
			{
				a_begin = (double *)a + j*lda;
				scalv( datatype, n_elem, x, incx, a_begin, inca );
			}
			break;
		}

		case COMPLEX:
		{
			scomplex *a_begin;
			for ( j = 0; j < n_iter; j++ )
			{
				a_begin = (scomplex *)a + j*lda;
				scalv( datatype, n_elem, x, incx, a_begin, inca );
			}
			break;
		}

		case DOUBLE_COMPLEX:
		{
			dcomplex *a_begin;
			for ( j = 0; j < n_iter; j++ )
			{
				a_begin = (dcomplex *)a + j*lda;
				scalv( datatype, n_elem, x, incx, a_begin, inca );
			}
			break;
		}
	}
}

void scalv( integer datatype, integer n, void* x, integer incx, void* y, integer incy )
{
	integer i;

	switch(datatype)
	{
		case FLOAT:
		{
			float *chi, *psi;
			for ( i = 0; i < n; ++i )
			{
				chi = (float *)x + i*incx;
				psi = (float *)y + i*incy;

				(*psi) = (*chi) * (*psi);
			}
			break;
		}

		case DOUBLE:
		{
			double *chi, *psi;
			for ( i = 0; i < n; ++i )
			{
				chi = (double *)x + i*incx;
				psi = (double *)y + i*incy;

				(*psi) = (*chi) * (*psi);
			}
			break;
		}

		case COMPLEX:
		{
			float *chi;
			scomplex *psi;

			for ( i = 0; i < n; ++i )
			{
				chi = (float *)x + i*incx;
				psi = (scomplex *)y + i*incy;

				psi->real = (*chi) * (psi)->real;
				psi->imag = (*chi) * (psi)->imag;
			}
			break;
		}

		case DOUBLE_COMPLEX:
		{
			double *chi;
			dcomplex *psi;
			for ( i = 0; i < n; ++i )
			{
				chi = (double *)x + i*incx;
				psi = (dcomplex *)y + i*incy;

				psi->real = (*chi) * (psi)->real;
				psi->imag = (*chi) * (psi)->imag;
			}
			break;
		}
	}
}