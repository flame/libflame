#include <lapacke.h>
#include "lapacke_gtest_helper.h"

/* helper routine: printing a matrix */
void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );

        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", mat[i*ldm+j] );
                printf( "\n" );
        }
}


/* helper routine: printing a matrix */
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n, double* mat, lapack_int ldm ) {
        lapack_int i, j;
        printf( "\n %s\n", desc );

        for( i = 0; i < m; i++ ) {
                for( j = 0; j < n; j++ ) printf( " %6.2f", mat[i+j*ldm] );
                printf( "\n" );
        }
}

/* helper routine: printing a vector of integers */
void print_vector( char* desc, lapack_int n, lapack_int* vec ) {
        lapack_int j;
        printf( "\n %s\n", desc );
        for( j = 0; j < n; j++ ) printf( " %6i", vec[j] );
        printf( "\n" );
}

/*Compute difference of two buffers of double type */
double computeDiff_d(int size, double *Out, double *Out_ref)
{
  int j;
  double diff = 0;
  double temp;
  for ( j = 0; j < size; j ++ ) {
    //diff += fabs(Out_ref[j] - Out[j]) ;
#if LAPACKE_TEST_VERBOSE
    temp = Out_ref[j] - Out[j];
    if (temp >10.0)
        printf("\n difference exceedes 10.0 at idx:%d ref: %lf  LF:%lf",
                    j, Out_ref[j], Out[j]);
#endif
    diff += (Out_ref[j] - Out[j]) * (Out_ref[j] - Out[j]) ;
  }
  return SQRT_DOUBLE(diff);
}

/*Compute difference of two buffers of float type */
float computeDiff_s(int size, float *Out, float *Out_ref)
{
  int j;
  float diff = 0;

  for ( j = 0; j < size; j ++ ) {
    //diff += fabs(Out_ref[j] - Out[j]) ;
    diff += (Out_ref[j] - Out[j]) * (Out_ref[j] - Out[j]);
#if LAPACKE_TEST_VERBOSE
    if ((fabs(Out_ref[j] - Out[j])) > 0.00001)
       printf("\n diff error exceeded threshold: at index:%d  diff:%f ref: %f  actual: %f\n",
   j, fabs(Out_ref[j] - Out[j]), Out_ref[j], Out[j] );
#endif

  }
  //return diff;
  return SQRT_FLOAT(diff);
}

/*Compute difference of two buffers of single precision complex type */
float computeDiff_c(int size, lapack_complex_float *Out, lapack_complex_float *Out_ref)
{
  int j;
  float diff = 0, temp;

  for ( j = 0; j < size; j ++ ) {
    temp = fabs( lapack_complex_float_real(Out_ref[j]) - lapack_complex_float_real(Out[j])) ;
    diff += temp *temp;
#if LAPACKE_TEST_VERBOSE
    if (fabs( lapack_complex_float_real(Out_ref[j]) - lapack_complex_float_real(Out[j]))> FLOAT_DIFF_THRESHOLD)
       printf("\n diff error exceeded threshold for real part @computeDiff_c : at index:%d  diff:%f \n" 
        , j,fabs( lapack_complex_float_real(Out_ref[j]) - lapack_complex_float_real(Out[j])));
#endif

    temp = fabs( lapack_complex_float_imag(Out_ref[j]) - lapack_complex_float_imag(Out[j])) ;
    diff += temp*temp;
#if LAPACKE_TEST_VERBOSE
    if (fabs( lapack_complex_float_imag(Out_ref[j]) - lapack_complex_float_imag(Out[j]))> FLOAT_DIFF_THRESHOLD)
       printf("\n diff error exceeded threshold for imag part @computeDiff_c : at index:%d  diff:%f \n" 
        ,j, fabs( lapack_complex_float_imag(Out_ref[j]) - lapack_complex_float_imag(Out[j])));
#endif 

  }
  //return diff;
  return SQRT_FLOAT(diff);
}

/*Compute difference of two buffers of double precision complex type */
double  computeDiff_z(int size, lapack_complex_double *Out, lapack_complex_double *Out_ref)
{
  int j;
  double diff = 0, dtemp = 0;

  for ( j = 0; j < size; j ++ ) {
    dtemp = fabs( lapack_complex_double_real(Out_ref[j]) - lapack_complex_double_real(Out[j])) ;
    diff += dtemp*dtemp;

#if LAPACKE_TEST_VERBOSE
    if (dtemp > DOUBLE_DIFF_THRESHOLD)
       printf("\n diff error exceeded threshold for real part @computeDiff_z : at index:%d  diff:%lf \n", j, dtemp);
#endif
    
    dtemp = fabs( lapack_complex_double_imag(Out_ref[j]) - lapack_complex_double_imag(Out[j]));
    diff += dtemp*dtemp;
#if LAPACKE_TEST_VERBOSE
    if (dtemp > DOUBLE_DIFF_THRESHOLD)
       printf("\n diff error exceeded threshold for imag part @computeDiff_z : at index:%d diff:%lf \n", j, dtemp);
#endif 
  }
  //return diff;
  return SQRT_DOUBLE(diff);
}

/*Compute difference of two buffers of int type */
int computeDiff_i(int size, int *Out, int *Out_ref)
{
  int j;
  int diff = 0;

  for ( j = 0; j < size; j ++ ) {
    diff += fabs(Out_ref[j] - Out[j]) ;
#if LAPACKE_TEST_VERBOSE
    if ((fabs(Out_ref[j] - Out[j])) > 0)
       printf("\n diff error exceeded threshold: at index:%d  diff:%f \n", j, fabs(Out_ref[j] - Out[j]));
#endif   
  }
  return diff;
}

void  lapacke_gtest_init_double_buffer_pair_rand( double *buf1, double *buf2, int size) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       //buf1[j] = (double) rand();
       buf1[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       buf2[j] = buf1[j];
       //printf("%lf \n",buf1[j]);
    }
}

void  lapacke_gtest_init_double_buffer_rand( double *buf1, int size) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
    }
}

void  lapacke_gtest_init_float_buffer_rand( float *buf1, int size) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
    }
}

void  lapacke_gtest_init_float_buffer_pair_rand( float *buf1, float *buf2, int size) 
{ 
   int j;
   int temp;
    for( j = 0; j < size; j++ ) {
	   //temp = rand() %3;
       buf1[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       buf2[j] = buf1[j];
    }
}

void  lapacke_gtest_init_int_buffer_pair_rand( int *buf1, int *buf2, int size)
{
    int j;
     for( j = 0; j < size; j++ ) {
        buf1[j] = ((int) rand()) / ((int) RAND_MAX);
        buf2[j] = buf1[j];
     }
}


void  lapacke_gtest_init_scomplex_buffer_pair_rand( lapack_complex_float *buf1, lapack_complex_float *buf2, int size) 
{ 
   int j;
   float stemp_r, stemp_i;
   
    for( j = 0; j < size; j++ ) {
       stemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
       stemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
       buf1[j] = lapack_make_complex_float( stemp_r, stemp_i);
       buf2[j] = lapack_make_complex_float( stemp_r, stemp_i);
    }
    
}

void  lapacke_gtest_init_scomplex_buffer_rand( lapack_complex_float *buf1, int size) 
{ 
   int j;
   float stemp_r, stemp_i;
   
    for( j = 0; j < size; j++ ) {
       stemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
       stemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
       buf1[j] = lapack_make_complex_float( stemp_r, stemp_i);
    }
    
}

void  lapacke_gtest_init_dcomplex_buffer_pair_rand( lapack_complex_double *buf1, lapack_complex_double *buf2, int size) 
{ 
   int j;
   double dtemp_r, dtemp_i;
   
    for( j = 0; j < size; j++ ) {
       dtemp_r =(((double) rand()) / ((double) RAND_MAX) - 0.5);
       dtemp_i =(((double) rand()) / ((double) RAND_MAX) - 0.5);
       buf1[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
       buf2[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
    }
    
}

void  lapacke_gtest_init_dcomplex_buffer_rand( lapack_complex_double *buf1, int size) 
{ 
   int j;
   double dtemp_r, dtemp_i;
   
    for( j = 0; j < size; j++ ) {
       dtemp_r =(((double) rand()) / ((double) RAND_MAX) - 0.5);
       dtemp_i =(((double) rand()) / ((double) RAND_MAX) - 0.5);
       buf1[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
    }
    
}

void  lapacke_gtest_init_double_buffer_pair_with_constant( double *buf1, double *buf2, int size, double constant) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = constant;
       buf2[j] = constant;
    }
}

void  lapacke_gtest_init_double_buffer_with_constant( double *buf1, int size, double constant) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = constant;
    }
}

void  lapacke_gtest_init_float_buffer_pair_with_constant( float *buf1, float *buf2, int size, float constant) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = constant;
       buf2[j] = constant;
    }
}

void  lapacke_gtest_init_float_buffer_with_constant( float *buf1, int size, float constant) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = constant;
    }
}

void  lapacke_gtest_init_int_buffer_pair_with_constant( int *buf1, int *buf2, int size, int constant) 
{
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = constant;
       buf2[j] = constant;
    }
}

void  lapacke_gtest_init_char_buffer_pair_with_constant( char *buf1, char *buf2, int size, char constant) 
{
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = constant;
       buf2[j] = constant;
    }
}

void  lapacke_gtest_init_scomplex_buffer_pair_with_constant( lapack_complex_float *buf1, lapack_complex_float *buf2, int size, float constant) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = lapack_make_complex_float( constant, constant);
       buf2[j] = lapack_make_complex_float( constant, constant);
    }
}

void  lapacke_gtest_init_scomplex_buffer_with_constant( lapack_complex_float *buf1, int size, float constant) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = lapack_make_complex_float( constant, constant);
    }
}

void  lapacke_gtest_init_dcomplex_buffer_pair_with_constant( lapack_complex_double *buf1, lapack_complex_double *buf2, int size, double constant) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = lapack_make_complex_double( constant, constant);
       buf2[j] = lapack_make_complex_double( constant, constant);
    }
}

void  lapacke_gtest_init_dcomplex_buffer_with_constant( lapack_complex_double *buf1, int size, double constant) 
{ 
   int j;
    for( j = 0; j < size; j++ ) {
       buf1[j] = lapack_make_complex_double( constant, constant);
    }
}


void  lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( double *buf1,
                                  double *buf2, int row, int col, char type) 
{ 
   int i, j;
   double temp;
   int m = (row<col)? row:col;


   switch ( type)
   {
       // Symmetric matrix
       case 'H' :
       case 'S' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = ((double) rand()) / ((double) RAND_MAX) - 0.5;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = temp;
                    buf2[i+m*j] = temp;
                    buf1[j+m*i] = temp;
                    buf2[j+m*i] = temp;
                }
            }
            break;
       case 'L' : // Lower triangular matrix
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = ((double) rand()) / ((double) RAND_MAX) - 0.5;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = 0.0;
                    buf2[i+m*j] = 0.0;
                    buf1[j+m*i] = temp;
                    buf2[j+m*i] = temp;
                }
            }
            break;
       case 'U' :  // upper triangular matrix
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = ((double) rand()) / ((double) RAND_MAX) - 0.5;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = temp;
                    buf2[i+m*j] = temp;
                    buf1[j+m*i] = 0.0;
                    buf2[j+m*i] = 0.0;
                }
            }
            break;
       case 'I' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = 0.0;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = temp;
                    buf2[i+m*j] = temp;
                    buf1[j+m*i] = temp;
                    buf2[j+m*i] = temp;
                }
            }
            break;
       default: 
           break;
   }
}

void  lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( float *buf1,
                                  float *buf2, int row, int col, char type)
{ 
   int i, j;
   float temp;
   int m = (row<col)? row:col;


   switch ( type)
   {
       // Symmetric matrix
       case 'S' :
       case 'H' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = ((float) rand()) / ((float) RAND_MAX) - 0.5;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = temp;
                    buf2[i+m*j] = temp;
                    buf1[j+m*i] = temp;
                    buf2[j+m*i] = temp;
                }
            }
            break;
       case 'L' :  // Lower triangular matrix
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = ((float) rand()) / ((float) RAND_MAX) - 0.5;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = 0.0;
                    buf2[i+m*j] = 0.0;
                    buf1[j+m*i] = temp;
                    buf2[j+m*i] = temp;
                }
            }
            break;
       case 'U' :  // upper triangular matrix
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = ((float) rand()) / ((float) RAND_MAX) - 0.5;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = temp;
                    buf2[i+m*j] = temp;
                    buf1[j+m*i] = 0.0;
                    buf2[j+m*i] = 0.0;
                }
            }
            break;
       case 'I' :  // upper triangular matrix
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = 0.0;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = temp;
                    buf2[i+m*j] = temp;
                    buf1[j+m*i] = temp;
                    buf2[j+m*i] = temp;
                }
            }
            break;
       default: 
           break;
   }
}

void  lapacke_gtest_init_float_buffer_rand_custom_matrix( float *buf1,
                                  int row, int col, char type)
{ 
   int i, j;
   float temp;
   int m = (row<col)? row:col;


   switch ( type)
   {
       // Symmetric matrix
       case 'S' :
       case 'H' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = ((float) rand()) / ((float) RAND_MAX) - 0.5;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = temp;
                    buf1[j+m*i] = temp;
                }
            }
            break;
       case 'L' :  // Lower triangular matrix
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = ((float) rand()) / ((float) RAND_MAX) - 0.5;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = 0.0;
                    buf1[j+m*i] = temp;
                }
            }
            break;
       case 'U' :  // upper triangular matrix
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    temp = ((float) rand()) / ((float) RAND_MAX) - 0.5;
                    if (i ==j)
                        temp =1.0;
                    buf1[i+m*j] = temp;
                    buf1[j+m*i] = 0.0;
                }
            }
            break;
       default: 
           break;
   }
}


void  lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( lapack_complex_float *buf1,
                     lapack_complex_float *buf2,  int row, int col, char type) 
{ 
   float stemp_r, stemp_i;
   int i, j;
   float temp;
   int m = (row<col)? row:col;

   switch ( type)
   {
       // Symmetric matrix
       case 'S' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    stemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
                    stemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
                    buf1[j+m*i] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf2[j+m*i] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf1[i+m*j] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf2[i+m*j] = lapack_make_complex_float(  stemp_r, stemp_i);
                }
            }
            break;
       
       // Hermitian matrix
       case 'H' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    stemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
                    stemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
                    buf1[j+m*i] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf2[j+m*i] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf1[i+m*j] = lapack_make_complex_float( stemp_r, (stemp_i* -1));
                    buf2[i+m*j] = lapack_make_complex_float(  stemp_r, (stemp_i * -1));
                }
            }
            break;
        // Lower traingular matrix
        case 'L' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    stemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
                    stemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
                    buf1[j+m*i] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf2[j+m*i] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf1[i+m*j] = lapack_make_complex_float( 0.0, 0.0);
                    buf2[i+m*j] = lapack_make_complex_float(  0.0, 0.0);
                }
            }           
            break;
        // Upper triangular
       case 'U' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    stemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
                    stemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
                    buf1[i+m*j] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf2[i+m*j] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf1[j+m*i] = lapack_make_complex_float( 0.0, 0.0);
                    buf2[j+m*i] = lapack_make_complex_float(  0.0, 0.0);
                }
            }
       case 'I' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    stemp_r = 0.0;
                    stemp_i = 0.0;
					if ( i == j)
					{
						stemp_r = 1.0;
						stemp_i = 1.0;
					}
					buf1[j+m*i] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf2[j+m*i] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf1[i+m*j] = lapack_make_complex_float( stemp_r, stemp_i);
                    buf2[i+m*j] = lapack_make_complex_float(  stemp_r, stemp_i);
                }
            }
            break;
            break;
       default: 
           break;
   }
    
}

void  lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( lapack_complex_double *buf1,
               lapack_complex_double *buf2,  int row, int col, char type)
{ 
   double stemp_r, stemp_i;
   int i, j;
   double temp;
   int m = (row<col)? row:col;

   switch ( type)
   {
       // Symmetric matrix
       case 'S' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    stemp_r =(((double) rand()) / ((double) RAND_MAX) - 0.5);
                    stemp_i =(((double) rand()) / ((double) RAND_MAX) - 0.5);
                    buf1[j+m*i] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf2[j+m*i] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf1[i+m*j] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf2[i+m*j] = lapack_make_complex_double(  stemp_r, stemp_i);
                }
            }
            break;
       
       // Hermitian matrix
       case 'H' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    stemp_r =(((double) rand()) / ((double) RAND_MAX) - 0.5);
                    stemp_i =(((double) rand()) / ((double) RAND_MAX) - 0.5);
                    buf1[j+m*i] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf2[j+m*i] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf1[i+m*j] = lapack_make_complex_double( stemp_r, (stemp_i* -1));
                    buf2[i+m*j] = lapack_make_complex_double(  stemp_r, (stemp_i * -1));
                }
            }
            break;
        // Lower traingular matrix
        case 'L' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    stemp_r =(((double) rand()) / ((double) RAND_MAX) - 0.5);
                    stemp_i =(((double) rand()) / ((double) RAND_MAX) - 0.5);
                    buf1[j+m*i] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf2[j+m*i] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf1[i+m*j] = lapack_make_complex_double( 0.0, 0.0);
                    buf2[i+m*j] = lapack_make_complex_double(  0.0, 0.0);
                }
            }           
            break;
        // Upper triangular
       case 'U' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    stemp_r =(((double) rand()) / ((double) RAND_MAX) - 0.5);
                    stemp_i =(((double) rand()) / ((double) RAND_MAX) - 0.5);
                    buf1[i+m*j] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf2[i+m*j] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf1[j+m*i] = lapack_make_complex_double( 0.0, 0.0);
                    buf2[j+m*i] = lapack_make_complex_double(  0.0, 0.0);
                }
            }
            break;
       case 'I' :
            for( i = 0; i < m; i++ ) {
                for (j = i; j< m; j++) {
                    stemp_r = 0.0;
                    stemp_i = 0.0;
					if (i == j)
					{
						stemp_r = 1.0;
						stemp_i = 1.0;
					}
                    buf1[j+m*i] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf2[j+m*i] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf1[i+m*j] = lapack_make_complex_double( stemp_r, stemp_i);
                    buf2[i+m*j] = lapack_make_complex_double(  stemp_r, stemp_i);
                }
            }
            break;
       default: 
           break;
   }
    
}

/* create pivot buffer with random values within the range 0 to n-1. **/
void  lapacke_gtest_init_pivot_buffer_pair_rand( int *ipiv, int *ipivref, int size)
{
    int j;
    int n = size;

    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = j;
       ipivref[j] = ipiv[j];
    }
    for(int j=0; j < n; j++) {
       int randIndex = (int) (rand() %  n);
       int tmp = ipiv[j];
       ipiv[j] = ipiv[randIndex];
       ipiv[randIndex] = tmp;
    }
    for( j = 0; j < n; j++ ) {
       ipivref[j] = ipiv[j];
    }

}

