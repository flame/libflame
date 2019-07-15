/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef FLA_EXTERN_DEFS_H
#define FLA_EXTERN_DEFS_H

extern __thread FLA_Obj FLA_THREE;
extern __thread FLA_Obj FLA_TWO;
extern __thread FLA_Obj FLA_ONE;
extern __thread FLA_Obj FLA_ONE_HALF;
extern __thread FLA_Obj FLA_ZERO;
extern __thread FLA_Obj FLA_MINUS_ONE_HALF;
extern __thread FLA_Obj FLA_MINUS_ONE;
extern __thread FLA_Obj FLA_MINUS_TWO;
extern __thread FLA_Obj FLA_MINUS_THREE;

extern __thread FLA_Obj FLA_EPSILON;
extern __thread FLA_Obj FLA_SAFE_MIN;
extern __thread FLA_Obj FLA_SAFE_MIN_SQUARE;
extern __thread FLA_Obj FLA_SAFE_INV_MIN;
extern __thread FLA_Obj FLA_SAFE_INV_MIN_SQUARE;
extern __thread FLA_Obj FLA_UNDERFLOW_THRES;
extern __thread FLA_Obj FLA_OVERFLOW_THRES;
extern __thread FLA_Obj FLA_UNDERFLOW_SQUARE_THRES;
extern __thread FLA_Obj FLA_OVERFLOW_SQUARE_THRES;

extern __thread const float    fzero;
extern __thread const double   dzero;
extern __thread const scomplex czero;
extern __thread const dcomplex zzero;

#endif

