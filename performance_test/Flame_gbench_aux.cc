
#include "Flame_gbench_aux.h"

void  Flame_gbench_init_buffer_rand( double *buf1, int size)
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
    }
}

void  Flame_gbench_init_buffer_rand( float *buf1, int size)
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
    }
}

void  Flame_gbench_init_buffer_rand( lapack_complex_float *buf1, int size)
{ 
   int j;
   float stemp_r, stemp_i;
   
    for( j = 0; j < size; j++ ) {
       stemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
       stemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
       //buf1[j] = lapack_make_complex_float( stemp_r, stemp_i);
	   buf1[j] = stemp_r + stemp_i*_Complex_I;

    }
    
}

void  Flame_gbench_init_buffer_rand( lapack_complex_double *buf1, int size)
{ 
   int j;
   double dtemp_r, dtemp_i;
   
    for( j = 0; j < size; j++ ) {
       dtemp_r =(((double) rand()) / ((double) RAND_MAX) - 0.5);
       dtemp_i =(((double) rand()) / ((double) RAND_MAX) - 0.5);
       //buf1[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
       buf1[j] = dtemp_r + dtemp_i*_Complex_I;
    }
    
}

void  Flame_gbench_buffer_memset( float *buf1, int size, float value)
{ 
    for( int j = 0; j < size; j++ ) {
       buf1[j] = value;
    }
}

void  Flame_gbench_buffer_memset( double *buf1, int size, double value)
{ 
    for( int j = 0; j < size; j++ ) {
       buf1[j] = value;
    }
}

void  Flame_gbench_buffer_memset( lapack_complex_double *buf1, int size, double rvalue, double ivalue)
{ 
    for( int j = 0; j < size; j++ ) {
       buf1[j] = rvalue + ivalue * _Complex_I;
    }
}

void  Flame_gbench_buffer_memset( lapack_complex_float *buf1, int size, float rvalue, float ivalue)
{ 
    for( int j = 0; j < size; j++ ) {
       buf1[j] = rvalue + ivalue * _Complex_I;
    }
}
