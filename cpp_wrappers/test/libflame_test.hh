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

/*Compute difference of two buffers with Template datatype */
template< typename T >
void computeError(int size, T *Out, T *Out_ref)
{
  int j;
  T diff = 0;
  
  for ( j = 0; j < size; j ++ ) {
    diff += abs(Out_ref[j] - Out[j]) ;
  }
  if(diff)
  {
    printf( "Failure Diff = %E\n", diff);
  }
  else{
    printf( "Success\n");
  }
  
}

/*Compute difference of two buffers of complex number type */
template< typename T >
void computeErrorComplex(int size, T *Out, T *Out_ref)
{
  int j;
  double diff = 0;
  for ( j = 0; j < size; j ++ ) {
    diff += abs(Out_ref[j].real - Out[j].real) ;
    diff += abs(Out_ref[j].imag - Out[j].imag) ;
  }
  if(diff)
  {
    printf( "Failure Diff = %E\n", diff);
  }
  else{
    printf( "Success\n");
  }
  
}

/*Allocate memory and initialise memory with random values*/
template< typename T >
void allocate_init_buffer(T *&Ain, T *&Aref, int size)
{
  Ain =  new T [size];
  Aref = (T *) malloc(size * sizeof(T));
  for(int i = 0; i < size; i++)
  {
    Ain[i] = ( (T) rand() / ((T) RAND_MAX / 2.0)) - 1.0;
    Aref[i] =  Ain[i] ; 
  }
}

/*Allocate memory and initialise memory with random values for complex number*/
template< typename T >
void allocate_init_buffer_complex(T *&Ain, T *&Aref, int size)
{
  Ain =  new T [size];
  Aref = (T *) malloc(size * sizeof(T));
  for(int i = 0; i < size; i++)
  {
    Ain[i].real = ( (float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;          
    Ain[i].imag = ( (float) rand() / ((float) RAND_MAX / 2.0)) - 1.0;          
    Aref[i].real = Ain[i].real;
    Aref[i].imag = Ain[i].imag;
  }
}

#endif        //  #ifndef LIBFLAME_TEST_HH