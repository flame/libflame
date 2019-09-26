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
using namespace std;

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
int getDataType()
{
  if (is_same<float, T>::value)
    return FLA_FLOAT;
  if (is_same<double, T>::value)
    return FLA_DOUBLE;
  if (is_same<scomplex, T>::value)
    return FLA_COMPLEX;
  if (is_same<dcomplex, T>::value)
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
    diff += abs(Out_ref[j].real - Out[j].real) ;
    diff += abs(Out_ref[j].imag - Out[j].imag) ;
  }
  return diff;
}

// Compute difference of two buffers with Template datatype
template< typename T >
double computeError(int n, int m, T *out, T *out_ref)
{
  double retDiff = 0;
  int dataType = getDataType<T>();
  if( dataType == FLA_FLOAT)
    retDiff = computeDiff<float>(n*m, (float*)out, (float*)out_ref);
  else if (dataType == FLA_DOUBLE)
    retDiff = computeDiff<double>(n*m, (double *)out, (double *)out_ref);
  else if (dataType == FLA_INT)
    retDiff = computeDiff<int>(n*m, (int *)out, (int *)out_ref );
  else if( dataType == FLA_COMPLEX)
    retDiff = computeErrorComplex<scomplex>(n*m, (scomplex *)out, (scomplex *)out_ref);
  else if (dataType == FLA_DOUBLE_COMPLEX)
    retDiff = computeErrorComplex<dcomplex>(n*m, (dcomplex *)out, (dcomplex *)out_ref);

  return retDiff;
}

// Allocate memory and initialise memory with random values
void allocate_init_buffer(int *&aIn, int *&aRef, int size)
{
  aIn =  new int [size];
  aRef = (int *) malloc(size * sizeof(int));
  for(int i = 0; i < size; i++)
  {
    aIn[i] = ( (int) rand() / ((int) RAND_MAX / 2.0)) - 1.0;
    aRef[i] =  aIn[i] ;
  }
}
void allocate_init_buffer(float *&aIn, float *&aRef, int size)
{
  aIn =  new float [size];
  aRef = (float *) malloc(size * sizeof(float));
  for(int i = 0; i < size; i++)
  {
    aIn[i] = ( (float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;
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

// Allocate memory and initialise memory with random values for complex number
void allocate_init_buffer(scomplex *&aIn, scomplex *&aRef, int size)
{
  aIn = new scomplex [size];
  aRef = new scomplex [size];
  for(int i = 0; i < size; i++)
  {
    aIn[i].real = ( (float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;
    aIn[i].imag = ( (float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;
    aRef[i].real = aIn[i].real;
    aRef[i].imag = aIn[i].imag;
  }
}
void allocate_init_buffer(dcomplex *&aIn, dcomplex *&aRef, int size)
{
  aIn =  new dcomplex [size];
  aRef = (dcomplex *) malloc(size * sizeof(dcomplex));
  for(int i = 0; i < size; i++)
  {
    aIn[i].real = ( (float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;
    aIn[i].imag = ( (float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;
    aRef[i].real = aIn[i].real;
    aRef[i].imag = aIn[i].imag;
  }
}

#endif        //  #ifndef LIBFLAME_TEST_HH