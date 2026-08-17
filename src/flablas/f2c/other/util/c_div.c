/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
#include "FLA_f2c.h"
 void c_div(scomplex *cp, scomplex *ap, scomplex *bp) {
 scomplex a = *ap;
 scomplex b = *bp;
 real temp;
 temp = b.real * b.real + b.imag * b.imag;
 cp->real = ( a.real * b.real + a.imag * b.imag ) / temp;
 cp->imag = ( a.imag * b.real - a.real * b.imag ) / temp;
 }
 
