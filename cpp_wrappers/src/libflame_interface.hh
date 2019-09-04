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

/*! @file libflame_interface.hh
 *  libflame_interface.hh defines all the Libflame CPP templated public interfaces
 *  */
#ifndef LIBFLAME_INTERFACE_HH
#define LIBFLAME_INTERFACE_HH

#include "libflame.hh"
namespace libflame{

/*! @brief LU factorization of a general M-by-N matrix A
 *         using partial pivoting with row interchanges. 
 * 
 *  @details
 * \b Purpose:
 * \verbatim
     LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges. 
     The factorization has the form
        A = P * L * U
     where P is a permutation matrix, L is lower triangular with unit diagonal elements (lower
     trapezoidal if M > N), and U is upper triangular (upper trapezoidal if M < N).
     
     This is the right-looking Level 3 BLAS version of the algorithm.  
\endverbatim
 * @param[in] M
M is int* \n
The number of rows of the matrix A.  M >= 0.
 * @param[in] N
N is int* \n
The number of columns of the matrix A.  N >= 0.
 * @param[in,out] A
A is float|double|complex|complex double array, dimension (LDA,N) \n
On entry, the M-by-N matrix to be factored. \n
On exit, the factors L and U from the factorization
A = P*L*U; the unit diagonal elements of L are not stored.
 * @param[in] LDA
LDA is int* \n
The leading dimension of the matrix A, LDA >= max(1,M)
 * @param[out] IPIV
IPIV is int array, dimension (min(M,N)) \n
The pivot indices; for 1 <= i <= min(M,N), row i of the
matrix was interchanged with row IPIV(i).
 * @param[out] INFO
INFO is int* \n
= 0:  successful exit \n
< 0:  if INFO = -i, the i-th argument had an illegal value \n
\> 0:  if INFO = i, U(i,i) is exactly zero. The factorization
      has been completed, but the factor U is exactly
      singular, and division by zero will occur if it is used
      to solve a system of equations.
 *  */ 
template< typename T >
void getrf(
     int* m, int* n, T* buff_A, int* ldim_A, int* buff_p, int* info )
{
  libflame::libflame_getrf(m, n, buff_A, ldim_A, buff_p, info);
}

}  // namespace libflame
#endif        //  #ifndef LIBFLAME_INTERFACE_HH


