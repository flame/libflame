/*
 * Copyright (C) 2025-2026, Advanced Micro Devices, Inc. All rights reserved.
 */
#include "FLA_f2c.h"
 void d_cnjg(dcomplex *dest, dcomplex *src) {
 dest->real = src->real ;
 dest->imag = -(src->imag);
 }
 
