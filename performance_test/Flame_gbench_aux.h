#ifndef _FLAME_GBENCH_AUX_
#define _FLAME_GBENCH_AUX_

#include <cstdio>
#include <iostream>
#include <memory>
#include <math.h>
#include <stdio.h>
#include "lapack.h"


void  Flame_gbench_init_buffer_rand( double *buf1, int size);

void  Flame_gbench_init_buffer_rand( float *buf1, int size);

void  Flame_gbench_init_buffer_rand( lapack_complex_float *buf1, int size);

void  Flame_gbench_init_buffer_rand( lapack_complex_double *buf1, int size);

void  Flame_gbench_buffer_memset( float *buf1, int size, float value);

void  Flame_gbench_buffer_memset( double *buf1, int size, double value);

void  Flame_gbench_buffer_memset( lapack_complex_double *buf1, int size, double rvalue, double ivalue);

void  Flame_gbench_buffer_memset( lapack_complex_float *buf1, int size, float rvalue, float ivalue);

#endif // _FLAME_GBENCH_AUX_