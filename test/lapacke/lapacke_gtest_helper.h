#ifndef _LAPACKE_EXAMPLE_HELPER_
#define _LAPACKE_EXAMPLE_HELPER_

#include <math.h>
#include <complex.h>

#define FLOAT_DIFF_THRESHOLD (10.0f)
#define DOUBLE_DIFF_THRESHOLD (10.0)
//#define LAPACKE_TEST_VERBOSE 0
#define SQRT_DOUBLE sqrt
#define SQRT_FLOAT  sqrtf

inline void  lapacke_gtest_alloc_double_buffer_pair( double **buf1, double **buf2, int siz)
{
   int i;
   *buf1 = (double *)malloc(siz*sizeof(double));
   *buf2 = (double *)malloc(siz*sizeof(double));
   for (i =0 ; i<siz; i++){
      (*buf1)[i] = 0.0;
	  (*buf2)[i] = 0.0;
   }
}

inline void  lapacke_gtest_alloc_float_buffer_pair( float **buf1, float **buf2, int siz)
{
   int i;
   *buf1 = (float *)malloc(siz*sizeof(float));
   *buf2 = (float *)malloc(siz*sizeof(float));
   for (i =0 ; i<siz; i++){
      (*buf1)[i] = 0.0;
	  (*buf2)[i] = 0.0;
   }

}

inline void  lapacke_gtest_alloc_lapack_scomplex_buffer_pair( lapack_complex_float **buf1, lapack_complex_float **buf2, int siz)
{
   int i;
   *buf1 = (lapack_complex_float *)malloc(siz*sizeof(lapack_complex_float));
   *buf2 = (lapack_complex_float *)malloc(siz*sizeof(lapack_complex_float));
   for (i =0 ; i<siz; i++){
       (*buf1)[i] = lapack_make_complex_float( 0.0, 0.0);
       (*buf2)[i] = lapack_make_complex_float( 0.0, 0.0);
   }
}

inline void  lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( lapack_complex_double **buf1, lapack_complex_double **buf2, int siz)
{
   int i;
   *buf1 = (lapack_complex_double *)malloc(siz*sizeof(lapack_complex_double));
   *buf2 = (lapack_complex_double *)malloc(siz*sizeof(lapack_complex_double));
   for (i =0 ; i<siz; i++){
       (*buf1)[i] = lapack_make_complex_double( 0.0, 0.0);
       (*buf2)[i] = lapack_make_complex_double( 0.0, 0.0);
   }
}

inline void  lapacke_gtest_alloc_int_buffer_pair( int **buf1, int **buf2, int siz)
{
   *buf1 = (int *)malloc(siz*sizeof(int));
   *buf2 = (int *)malloc(siz*sizeof(int));
}

inline void  lapacke_gtest_alloc_char_buffer_pair( char **buf1, char **buf2, int siz)
{
   *buf1 = (char *)malloc(siz*sizeof(char));
   *buf2 = (char *)malloc(siz*sizeof(char));
}


void  lapacke_gtest_init_double_buffer_pair_rand( double *buf1, double *buf2, int size);
void  lapacke_gtest_init_double_buffer_rand( double *buf1, int size);

void  lapacke_gtest_init_double_buffer_pair_with_constant( double *buf1, double *buf2, 
                                                            int size, double constant);

void  lapacke_gtest_init_double_buffer_with_constant( double *buf1, int size, double constant);


void  lapacke_gtest_init_float_buffer_pair_rand( float *buf1, float *buf2, int size);
void  lapacke_gtest_init_float_buffer_rand( float *buf1, int size);

void  lapacke_gtest_init_float_buffer_pair_with_constant( float *buf1, float *buf2, 
                                                            int size, float constant);

void  lapacke_gtest_init_float_buffer_with_constant( float *buf1, int size, float constant);

void  lapacke_gtest_init_int_buffer_pair_with_constant( int *buf1, int *buf2,
                                                            int size, int constant);

void  lapacke_gtest_init_int_buffer_pair_rand( int *buf1, int *buf2, int size);

void  lapacke_gtest_init_char_buffer_pair_with_constant( char *buf1, char *buf2, int size, char constant);

void  lapacke_gtest_init_scomplex_buffer_pair_rand( lapack_complex_float *buf1, 
													lapack_complex_float *buf2, 
													int size);	

void  lapacke_gtest_init_scomplex_buffer_rand( lapack_complex_float *buf1, 
													int size);	

void  lapacke_gtest_init_scomplex_buffer_pair_with_constant( lapack_complex_float *buf1,
															 lapack_complex_float *buf2,
															 int size,
															 float constant) ;

void  lapacke_gtest_init_scomplex_buffer_with_constant( lapack_complex_float *buf1,
													int size, float constant);

void  lapacke_gtest_init_dcomplex_buffer_pair_rand( lapack_complex_double *buf1,
													lapack_complex_double *buf2,
													int size);

void  lapacke_gtest_init_dcomplex_buffer_rand( lapack_complex_double *buf1,
													int size);

void  lapacke_gtest_init_dcomplex_buffer_pair_with_constant( lapack_complex_double *buf1,
													lapack_complex_double *buf2,
													int size,
													double constant);

void  lapacke_gtest_init_dcomplex_buffer_with_constant( lapack_complex_double *buf1, 
													int size, double constant);

// Routines for Custom matrix generation
void  lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( double *buf1,
                                  double *buf2, int row, int col, char type);

void  lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( float *buf1,
                                  float *buf2, int row, int col, char type);

void  lapacke_gtest_init_float_buffer_rand_custom_matrix( float *buf1,
                                  int row, int col, char type);

void  lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( lapack_complex_float *buf1,
                     lapack_complex_float *buf2,  int row, int col, char type) ;
                     
void  lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( lapack_complex_double *buf1,
               lapack_complex_double *buf2,  int row, int col, char type);


void print_matrix_rowmajor( char* desc, lapack_int m, lapack_int n,
                        double* mat, lapack_int ldm );
void print_matrix_colmajor( char* desc, lapack_int m, lapack_int n,
                        double* mat, lapack_int ldm );
void print_vector( char* desc, lapack_int n, lapack_int* vec );

/*Compute difference of two buffers of double type */
double computeDiff_d(int size, double *Out, double *Out_ref);

/*Compute difference of two buffers of float type */
float computeDiff_s(int size, float *Out, float *Out_ref);

/*Compute difference of two buffers of single precision complex type */
float computeDiff_c(int size, lapack_complex_float *Out, lapack_complex_float *Out_ref);

/*Compute difference of two buffers of double precision complex type */
double computeDiff_z(int size, lapack_complex_double *Out, lapack_complex_double *Out_ref);

/*Compute difference of two buffers of int type */
int computeDiff_i(int size, int *Out, int *Out_ref);

/* create pivot buffer with random values within the range 0 to n-1. **/
void  lapacke_gtest_init_pivot_buffer_pair_rand( int *ipiv, int *ipivref, int size);

#endif /* _LAPACKE_EXAMPLE_HELPER_*/
