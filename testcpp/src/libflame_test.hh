/******************************************************************************
* Copyright (c) 2019 - present Advanced Micro Devices, Inc. All rights reserved.
*
* Permission is hereby granted, free of charge, to any person obtaining a copy
* of this software and associated documentation files (the "Software"), to deal
* in the Software without restriction, including without limitation the rights
* to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
* copies of the Software, and to permit persons to whom the Software is
* furnished to do so, subject to the following conditions:
*
* The above copyright notice and this permission notice shall be included in
* all copies or substantial portions of the Software.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
* AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
* LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
* OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
* THE SOFTWARE.
*******************************************************************************/

/*! @file libflame_test.hh
 *  libflame_test.hh defines all the common functions which is required
 *  to validate CPP template interfaces for all libflame modules
 *  */
#ifndef LIBFLAME_TEST_HH
#define LIBFLAME_TEST_HH

#include "libflame_interface.hh"
//#include "FLAME.h"
//#define FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES 1

#define FLA_ENABLE_EXTERNAL_LAPACK_INTERFACES

using namespace std;
FLA_Error geqpf_C( FLA_Obj A, int *jpvt, FLA_Obj t, FLA_Obj rwork );

template< typename T >
T* typecast( FLA_Obj buff, int dataType)
{
  if(dataType == FLA_FLOAT)
    return  (T * ) FLA_FLOAT_PTR( buff );
  else if(dataType == FLA_DOUBLE)
    return  (T * ) FLA_DOUBLE_PTR( buff );
  if(dataType == FLA_COMPLEX)
    return  (T * ) FLA_COMPLEX_PTR( buff );
  else if(dataType == FLA_DOUBLE_COMPLEX)
    return  (T * ) FLA_DOUBLE_COMPLEX_PTR( buff );
  else if(dataType == FLA_INT)
    return  (T * ) FLA_INT_PTR( buff );
}
template< typename T >
int getDatatype()
{
  if (is_same<float, T>::value)
    return FLA_FLOAT;
  if (is_same<double, T>::value)
    return FLA_DOUBLE;
  if (is_same<lapack_complex_float, T>::value)
    return FLA_COMPLEX;
  if (is_same<lapack_complex_double, T>::value)
    return FLA_DOUBLE_COMPLEX;
  if (is_same<int, T>::value)
    return FLA_INT;
}
template< typename T >
double computeDiff(int size, T *Out, T *Out_ref)
{
  int j;
  T diff = 0;

  for ( j = 0; j < size; j ++ ) {
   diff += abs(Out_ref[j] - Out[j]) ;
  }
  return diff;
}

/*Compute difference of two buffers of complex number type */
template< typename T >
double computeErrorComplex(int size, T *Out, T *Out_ref)
{
  int j;
  double diff = 0;
  for ( j = 0; j < size; j ++ ) {
 // if((isnan(creal(Out_ref[j])) && isnan(creal(Out_ref[j])))
 //   || (isnan(cimag(Out_ref[j])) && isnan(cimag(Out_ref[j]))))
    {
      diff += abs(creal(Out_ref[j]) - creal(Out[j])) ;
      diff += abs(cimag(Out_ref[j]) - cimag(Out[j])) ;
    }
  }
   return diff;
}

FILE *fp;
   //fp = fopen("test/in","a+");
   //print(fp,m*m,aInBuff, aRefBuff);
   //fclose(fp);
   //  fp = fopen("test/out","a+");
   //print(fp,m*m,aInBuff, aRefBuff);
   //fclose(fp);
   //  fp = fopen("test/out1","a+");
   //print(fp,m*m,aInBuff, aInBuff);
   //fclose(fp);
   //  fp = fopen("test/out2","a+");
   //print(fp,m*m,aRefBuff, aRefBuff);
   //fclose(fp);

void print(FILE *fp, int size, int *Out, int *Out_ref)
{
  fprintf(fp, "\n****Starts*****\n");
  for (int j = 0; j < size; j ++ ) {
   fprintf(fp,"%d %d\t",Out[j], Out_ref[j]);
  }
  fprintf(fp,"\n****Ends*****\n");
}

void print(FILE *fp, int size, float *Out, float *Out_ref)
{
  fprintf(fp, "\n****Starts*****\n");
  for (int j = 0; j < size; j ++ ) {
   fprintf(fp,"%e %e\t",Out[j], Out_ref[j]);
  }
  fprintf(fp,"\n****Ends*****\n");
}

void print(FILE *fp, int size, double *Out, double *Out_ref)
{
  fprintf(fp, "\n****Starts*****\n");
  for (int j = 0; j < size; j ++ ) {
   fprintf(fp,"%e %e\t",Out[j], Out_ref[j]);
  }
  fprintf(fp,"\n****Ends*****\n");
}

void print(FILE *fp,int size, lapack_complex_float *Out, lapack_complex_float *Out_ref)
{
  int j;

  for ( j = 0; j < size; j ++ ) {
    fprintf(fp," %e %e %e %e \n",creal(Out_ref[j]) , creal(Out[j]),cimag(Out_ref[j]) , cimag(Out[j])) ;
  }
  fprintf(fp,"\n");
}

void print(FILE *fp, int size, lapack_complex_double *Out, lapack_complex_double *Out_ref)
{
  int j;

  for ( j = 0; j < size; j ++ ) {
    fprintf(fp," %e %e %e %e \n",creal(Out_ref[j]) , creal(Out[j]),cimag(Out_ref[j]) , cimag(Out[j])) ;
  }
  fprintf(fp,"\n");
}

// Compute difference of two buffers with Template datatype
template< typename T >
double computeError(int n, int m, T *out, T *out_ref)
{
  double retDiff = 0;
  int dataType = getDatatype<T>();
  if( dataType == FLA_FLOAT)
    retDiff = computeDiff<float>(n*m, (float*)out, (float*)out_ref);
  else if (dataType == FLA_DOUBLE)
    retDiff = computeDiff<double>(n*m, (double *)out, (double *)out_ref);
  else if (dataType == FLA_INT)
    retDiff = computeDiff<int>(n*m, (int *)out, (int *)out_ref );
  else if( dataType == FLA_COMPLEX || dataType == FLA_COMPLEX )
    retDiff = computeErrorComplex<lapack_complex_float>(n*m, (lapack_complex_float *)out, (lapack_complex_float *)out_ref);
  else if (dataType == FLA_DOUBLE_COMPLEX)
    retDiff = computeErrorComplex<lapack_complex_double>(n*m, (lapack_complex_double *)out, (lapack_complex_double *)out_ref);

  return retDiff;
}

// Allocate memory and initialise memory with random values
void allocate_init_buffer(int *&aIn, int *&aRef, int size)
{
  aIn =  new int [size];
  aRef = new int [size];
  for(int i = 0; i < size; i++)
  {
    aIn[i] = ((int) rand() / ((int) RAND_MAX / 2.0)) - 1.0;
    aRef[i] =  aIn[i] ;
  }
}

void allocate_init_buffer(float *&aIn, float *&aRef, int size)
{
  aIn =  new float [size];
  aRef = new float [size];
  for(int i = 0; i < size; i++)
  {
    aIn[i] = ((float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;
    aRef[i] =  aIn[i] ;
  }

}
void allocate_init_buffer(double *&aIn, double *&aRef, int size)
{
  aIn =  new double [size];
  aRef = new double [size];
  for(int i = 0; i < size; i++)
  {
    aIn[i] = ( (double) rand() / ((double) RAND_MAX / 2.0)) - 1.0;
    aRef[i] =  aIn[i] ;
  }
}

/*Allocate memory and initialise memory with random values for complex number*/
void allocate_init_buffer(lapack_complex_float *&Ain, lapack_complex_float *&Aref, int size)
{
  Ain =  new lapack_complex_float [size];
  Aref = (lapack_complex_float *) malloc(size * sizeof(lapack_complex_float));
  for(int i = 0; i < size; i++)
  {
    float real = ( (float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;
    float imag = ( (float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;
    Ain[i] =  lapack_make_complex_float( real, imag);
    Aref[i] = Ain[i];
  }
}

void allocate_init_buffer(lapack_complex_double *&Ain, lapack_complex_double *&Aref, int size)
{
  Ain =  new lapack_complex_double [size];
  Aref = (lapack_complex_double *) malloc(size * sizeof(lapack_complex_double));
  for(int i = 0; i < size; i++)
  {
    double real = ( (double) rand() / ((double) RAND_MAX / 2.0)) - 1.0;
    double imag = ( (double) rand() / ((double) RAND_MAX / 2.0)) - 1.0;
    Ain[i] =  lapack_make_complex_float( real, imag);
    Aref[i] = Ain[i];
  }
}

void allocate_init_buffer(float *&aIn, float *&aRef, int size, int value)
{
  aIn =  new float [size];
  aRef = new float [size];
  for(int i = 0; i < size; i++)
  {
    aIn[i] = value;
    aRef[i]  = value;
  }

}
void allocate_init_buffer(double *&aIn, double *&aRef, int size, int value)
{
  aIn =  new double [size];
  aRef = new double [size];
  for(int i = 0; i < size; i++)
  {
    aIn[i] = value;
    aRef[i]  = value;
  }
}

#endif        //  #ifndef LIBFLAME_TEST_HH