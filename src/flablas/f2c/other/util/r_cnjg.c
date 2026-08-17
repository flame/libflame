/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
#include "FLA_f2c.h"
 void r_cnjg(scomplex *dest, scomplex *src) {
 dest->real = src->real ;
 dest->imag = -(src->imag);
 }
 
