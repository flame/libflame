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

/*! @file liblame_test.cc
 *  libflame_test.cc Test application to validate CPP template interfaces
 *  for all libflame modules
 *  */

#include "libflame_test.hh"

template< typename T >
void getrf_test()
{
  int rval;
  int m = 1000; 
  int n = 2000;  

  int *piv_buff_ref = (int *) malloc(m * sizeof(int));
  int *piv_buff =  new int [m];
  
  srand (time(NULL));
  if(is_same<float, T>::value)
  {
	float *Ain,  *Aref;
	/*Allocate memory for in and out buffer*/
	allocate_init_buffer(Ain, Aref, m*n);
	/*Call C and CPP function*/
    sgetrf_(&m, &n, Aref, &m, piv_buff_ref, &rval);
    libflame::getrf(&m, &n, Ain, &m, piv_buff, &rval);
	/*Compute difference of CPP interface against C ref function*/
	printf("F77_sgetrf() Output comparison : ");
	computeError(m*n, Ain, Aref);
	printf("F77_sgetrf() Pivot comparison : ");
	computeError<int>(m, piv_buff, piv_buff_ref);
  }
  else if(is_same<double, T>::value)
  {
	double *Ain,  *Aref;
	allocate_init_buffer(Ain, Aref, m*n);
    dgetrf_(&m, &n, Aref, &m, piv_buff_ref, &rval);
    libflame::getrf(&m, &n, Ain, &m, piv_buff, &rval);
	printf("F77_dgetrf() Output comparison : ");
	computeError(m*n, Ain, Aref);
	printf("F77_dgetrf() Pivot comparison : ");
	computeError<int>(m, piv_buff, piv_buff_ref);
  }
  else if(is_same<std::complex<float>, T>::value)
  {
	scomplex *Ain,  *Aref;
	allocate_init_buffer_complex(Ain, Aref, m*n);
    cgetrf_(&m, &n, Aref, &m, piv_buff_ref, &rval);
    libflame::getrf(&m, &n, Ain, &m, piv_buff, &rval);
	printf("F77_cgetrf() Output comparison : ");
	computeErrorComplex(m*n, Ain, Aref);
	printf("F77_cgetrf() Pivot comparison : ");
	computeError<int>(m, piv_buff, piv_buff_ref);
  }
  else if(is_same<std::complex<double>, T>::value)
  {
	dcomplex *Ain,  *Aref;
	allocate_init_buffer_complex(Ain, Aref, m*n);
    zgetrf_(&m, &n, Aref, &m, piv_buff_ref, &rval);
    libflame::getrf(&m, &n, Ain, &m, piv_buff, &rval);
	printf("F77_zgetrf() Output comparison : ");
	computeErrorComplex(m*n, Ain, Aref);
	printf("F77_zgetrf() Pivot comparison : ");
	computeError<int>(m, piv_buff, piv_buff_ref);
  }
 }
 
int main(int argc, char *argv[])
{  
	getrf_test<float>();
	getrf_test<double>();
	getrf_test<std::complex<float>>();
	getrf_test<std::complex<double>>();
}
