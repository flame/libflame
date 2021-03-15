/*

    Copyright (C) 2014, The University of Texas at Austin

    This file is part of libflame and is available under the 3-Clause
    BSD license, which can be found in the LICENSE file at the top-level
    directory, or at http://opensource.org/licenses/BSD-3-Clause

*/

#ifndef FLA_EXTERN_DEFS_H
#define FLA_EXTERN_DEFS_H

extern TLS_CLASS_SPEC FLA_Obj FLA_THREE;
extern TLS_CLASS_SPEC FLA_Obj FLA_TWO;
extern TLS_CLASS_SPEC FLA_Obj FLA_ONE;
extern TLS_CLASS_SPEC FLA_Obj FLA_ONE_HALF;
extern TLS_CLASS_SPEC FLA_Obj FLA_ZERO;
extern TLS_CLASS_SPEC FLA_Obj FLA_MINUS_ONE_HALF;
extern TLS_CLASS_SPEC FLA_Obj FLA_MINUS_ONE;
extern TLS_CLASS_SPEC FLA_Obj FLA_MINUS_TWO;
extern TLS_CLASS_SPEC FLA_Obj FLA_MINUS_THREE;

extern TLS_CLASS_SPEC FLA_Obj FLA_EPSILON;
extern TLS_CLASS_SPEC FLA_Obj FLA_SAFE_MIN;
extern TLS_CLASS_SPEC FLA_Obj FLA_SAFE_MIN_SQUARE;
extern TLS_CLASS_SPEC FLA_Obj FLA_SAFE_INV_MIN;
extern TLS_CLASS_SPEC FLA_Obj FLA_SAFE_INV_MIN_SQUARE;
extern TLS_CLASS_SPEC FLA_Obj FLA_UNDERFLOW_THRES;
extern TLS_CLASS_SPEC FLA_Obj FLA_OVERFLOW_THRES;
extern TLS_CLASS_SPEC FLA_Obj FLA_UNDERFLOW_SQUARE_THRES;
extern TLS_CLASS_SPEC FLA_Obj FLA_OVERFLOW_SQUARE_THRES;

extern TLS_CLASS_SPEC const float    fzero;
extern TLS_CLASS_SPEC const double   dzero;
extern TLS_CLASS_SPEC const scomplex czero;
extern TLS_CLASS_SPEC const dcomplex zzero;

#endif

