/*
    Copyright (c) 2020 Advanced Micro Devices, Inc.  All rights reserved.
    Sep 24, 2020
*/

#include "FLAME.h"
#include "test_libflame.h"
#include "test_common.h"

/* Initialize matrix with random values */

void rand_matrix_s( float * A, int M, int N, int LDA )
{
   int i, j;

   for( i = 0; i < N; i++ )
   {
      for( j = 0; j < M; j++ )
      {
         A[i * LDA + j] = SRAND();
      }
   }

   return;
}

void rand_matrix_d( double * A, int M, int N, int LDA )
{
   int i, j;

   for( i = 0; i < N; i++ )
   {
      for( j = 0; j < M; j++ )
      {
         A[i * LDA + j] = DRAND();
      }
   }

   return;
}

void rand_matrix_c( scomplex * A, int M, int N, int LDA )
{
   int i, j;

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

void rand_matrix_z( dcomplex * A, int M, int N, int LDA )
{
   int i, j;

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

void rand_sym_matrix_s( float * A, int M, int N, int LDA )
{
   int i, j;

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

void rand_sym_matrix_d( double * A, int M, int N, int LDA )
{
   int i, j;

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

void rand_sym_matrix_c( scomplex * A, int M, int N, int LDA )
{
   int i, j;

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

void rand_sym_matrix_z( dcomplex * A, int M, int N, int LDA )
{
   int i, j;

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

void copy_matrix_s( float *sM, float *dM, int M, int N, int LDS, int LDD )
{
   int i, j;

   for( i = 0; i < N; i++ )
   {
      for( j = 0; j < M; j++ )
      {
         dM[ i * LDD + j ] = sM[ i * LDS + j ];
      }
   }

   return;
}

void copy_matrix_d( double *sM, double *dM, int M, int N, int LDS, int LDD )
{
   int i, j;

   for( i = 0; i < N; i++ )
   {
      for( j = 0; j < M; j++ )
      {
         dM[ i * LDD + j ] = sM[ i * LDS + j ];
      }
   }

   return;
}

void copy_matrix_c( scomplex *sM, scomplex *dM, int M, int N, int LDS, int LDD )
{
   int i, j;

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

void copy_matrix_z( dcomplex *sM, dcomplex *dM, int M, int N, int LDS, int LDD )
{
   int i, j;

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

void pack_matrix_lt_s( float *A, float *B, int N, int LDA )
{
   int i, j;
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

void pack_matrix_lt_d( double *A, double *B, int N, int LDA )
{
   int i, j;
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

void pack_matrix_lt_c( scomplex *A, scomplex *B, int N, int LDA )
{
   int i, j;
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

void pack_matrix_lt_z( dcomplex *A, dcomplex *B, int N, int LDA )
{
   int i, j;
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

void reset_matrix_s( float *A, int M, int N, int LDA )
{
   int i, j;

   for( i = 0; i < N; i++ )
   {
      for( j = 0; j < M; j++ )
      {
         A[ i * LDA + j ] = 0.f;
      }
   }

   return;
}

void reset_matrix_d( double *A, int M, int N, int LDA )
{
   int i, j;

   for( i = 0; i < N; i++ )
   {
      for( j = 0; j < M; j++ )
      {
         A[ i * LDA + j ] = 0.;
      }
   }

   return;
}

void reset_matrix_c( scomplex *A, int M, int N, int LDA )
{
   int i, j;

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

void reset_matrix_z( dcomplex *A, int M, int N, int LDA )
{
   int i, j;

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

void set_identity_s( float *A, int M, int N, int LDA )
{
   int i, j;

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

void set_identity_d( double *A, int M, int N, int LDA )
{
   int i, j;

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

void set_identity_c( scomplex *A, int M, int N, int LDA )
{
   int i, j;

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

void set_identity_z( dcomplex *A, int M, int N, int LDA )
{
   int i, j;

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
   cp->real = ( a.real * b.real + a.imag * b.imag ) / temp;
   cp->imag = ( a.imag * b.real - a.real * b.imag ) / temp;
}

