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

/*! @file libflame.hh
 *  libflame.hh defines all the overloaded CPP functions to be invoked from 
 *  template interfaces
 *  */
#ifndef LIBFLAME_HH
#define LIBFLAME_HH

#include  <complex>
extern "C" {
#include "FLAME.h"
#include "../../src/map/lapack2flamec/FLA_lapack2flame_util_defs.h"
}

namespace libflame{
  // ==========================================================================
  // getrf()
  // --------------------------------------------------------------------------

  inline void
  libflame_getrf(
     int* m, int* n, lapack_stype * buff_A, int* ldim_A, int* buff_p, int* info )
  {
//      F77_sgetrf(m, n, buff_A, ldim_A, buff_p, info );
    sgetrf_(m, n, buff_A, ldim_A, buff_p, info);
      
  }

  inline void
  libflame_getrf(
    int* m, int* n, lapack_dtype* buff_A, int* ldim_A, int* buff_p, int* i)
  {
    dgetrf_(m, n, buff_A, ldim_A, buff_p, i);
  }

  inline void
  libflame_getrf(
     int* m, int* n, lapack_ctype* buff_A, int* ldim_A, int* buff_p, int* info )
  {
    cgetrf_(m, n, buff_A, ldim_A, buff_p, info );
  }

  inline void
  libflame_getrf(
    int* m, int* n, lapack_ztype* buff_A, int* ldim_A, int* buff_p, int* info )
  {
    zgetrf_(m, n, buff_A, ldim_A, buff_p, info );
  }

}//namespace libflame

#endif        //  #ifndef LIBFLAME_HH

