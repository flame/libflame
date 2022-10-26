//#include "lapacke_utils.h"
#include "lapacke_gtest_lin_common.hh"

/* lin_solver_double_parameters  class constructor definition */
lin_solver_double_parameters::lin_solver_double_parameters (int matrix_layout_i,
   char uplo_i, int n_i, int nrhs_i, int lda_i, int ldb_i)
{
	int j;

	matrix_layout = matrix_layout_i;
	uplo = uplo_i;
	n = n_i;
	nrhs = nrhs_i;
	lda = lda_i;
	ldb = ldb_i;
	
/* Memory allocation of the buffers */
        a = (double *)malloc(n*n*sizeof(double));
	aref = (double *)malloc(n*n*sizeof(double));
   
        ipiv    = (int *)malloc(n*sizeof(int));
	ipivref = (int *)malloc(n*sizeof(int));
  
	b    = (double *)malloc((ldb*nrhs)*sizeof(double));
	bref = (double *)malloc((ldb*nrhs)*sizeof(double));

	
	if( (a==NULL) || (aref==NULL) || (b==NULL) || (bref==NULL)|| \
	    (ipiv==NULL) || (ipivref==NULL) ){

       printf(" lin_solver_double_parameters object: malloc error. Exiting...\n");
       free (a   );
       free (aref);
       free (b   );
       free (bref);
       free (ipiv   );
       free (ipivref);
	   exit (-1);
	}

    /* Initialization of input Buffers */
	  
    for( j = 0; j < n*n; j++ ) {
       a[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       aref[j] = a[j];
    }

    for( j = 0; j < (ldb*nrhs); j++ ) {
       b[j] = ((double) rand()) / ((double) RAND_MAX) - 0.5;
       bref[j] = b[j];
    }
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
	
	info = 0;
	inforef = 0;
	
}

/* lin_solver_double_parameters  class Destructor definition **/
lin_solver_double_parameters:: ~lin_solver_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
       printf(" lin_solver_double_parameters object: destructor invoked. \n");
#endif
       free (a   );
       free (aref);
       free (b   );
       free (bref);
       free (ipiv   );
       free (ipivref);
}

/* lin_solver_float_parameters  class constructor definition */
lin_solver_float_parameters::lin_solver_float_parameters (int matrix_layout_i,
   char uplo_i, int n_i, int nrhs_i, int lda_i, int ldb_i)
{
	int j;

	matrix_layout = matrix_layout_i;
	uplo = uplo_i;
	n = n_i;
	nrhs = nrhs_i;
	lda = lda_i;
	ldb = ldb_i;
	
/* Memory allocation of the buffers */
        a = (float *)malloc(n*n*sizeof(float));
	aref = (float *)malloc(n*n*sizeof(float));
   
        ipiv    = (int *)malloc(n*sizeof(int));
	ipivref = (int *)malloc(n*sizeof(int));
  
	b    = (float *)malloc((ldb*nrhs)*sizeof(float));
	bref = (float *)malloc((ldb*nrhs)*sizeof(float));

	
	if( (a==NULL) || (aref==NULL) || (b==NULL) || (bref==NULL)|| \
	    (ipiv==NULL) || (ipivref==NULL) ){

       printf(" lin_solver_float_parameters object: malloc error. Exiting...\n");
       free (a   );
       free (aref);
       free (b   );
       free (bref);
       free (ipiv   );
       free (ipivref);
	   exit (-1);
	}

    /* Initialization of input Buffers */
	  
    for( j = 0; j < n*n; j++ ) {
       a[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       aref[j] = a[j];
    }

    for( j = 0; j < (ldb*nrhs); j++ ) {
       b[j] = ((float) rand()) / ((float) RAND_MAX) - 0.5;
       bref[j] = b[j];
    }
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
	
	info = 0;
	inforef = 0;
	
}

/* Destructor definition **/
lin_solver_float_parameters:: ~lin_solver_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
       printf(" lin_solver_float_parameters object: destructor invoked. \n");
#endif
       free (a   );
       free (aref);
       free (b   );
       free (bref);
       free (ipiv   );
       free (ipivref);
}

/* lin_solver_scomplex_parameters  class constructor definition */
lin_solver_scomplex_parameters::lin_solver_scomplex_parameters (int matrix_layout_i,
   char uplo_i, int n_i, int nrhs_i, int lda_i, int ldb_i)
{
	int j;

	matrix_layout = matrix_layout_i;
	uplo = uplo_i;
	n = n_i;
	nrhs = nrhs_i;
	lda = lda_i;
	ldb = ldb_i;
	
/* Memory allocation of the buffers */
    a = (lapack_complex_float *)malloc(n*n*sizeof(lapack_complex_float));
	aref = (lapack_complex_float *)malloc(n*n*sizeof(lapack_complex_float));
   
    ipiv    = (int *)malloc(n*sizeof(int));
	ipivref = (int *)malloc(n*sizeof(int));
  
	b    = (lapack_complex_float *)malloc((ldb*nrhs)*sizeof(lapack_complex_float));
	bref = (lapack_complex_float *)malloc((ldb*nrhs)*sizeof(lapack_complex_float));

	
	if( (a==NULL) || (aref==NULL) || (b==NULL) || (bref==NULL)|| \
	    (ipiv==NULL) || (ipivref==NULL) ){

       printf(" lin_solver_scomplex_parameters object: malloc error. Exiting...\n");
       free (a   );
       free (aref);
       free (b   );
       free (bref);
       free (ipiv   );
       free (ipivref);
	   exit (-1);
	}

    /* Initialization of input Buffers */
	  
    for( j = 0; j < n*n; j++ ) {
	   float dtemp_r, dtemp_i;
	   dtemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   dtemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   a[j] = lapack_make_complex_float( dtemp_r, dtemp_i);
       aref[j] = lapack_make_complex_float( dtemp_r, dtemp_i);
    }

    for( j = 0; j < (ldb*nrhs); j++ ) {
	   float dtemp_r, dtemp_i;
	   dtemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   dtemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   b[j] = lapack_make_complex_float( dtemp_r, dtemp_i);
       bref[j] = lapack_make_complex_float( dtemp_r, dtemp_i);
    }
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
	
	info = 0;
	inforef = 0;
	
}

/* Destructor definition: lin_solver_scomplex_parameters  class **/
lin_solver_scomplex_parameters:: ~lin_solver_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
       printf(" lin_solver_scomplex_parameters object: destructor invoked. \n");
#endif
       free (a   );
       free (aref);
       free (b   );
       free (bref);
       free (ipiv   );
       free (ipivref);
}

/* 'lin_solver_dcomplex_parameters'  class constructor definition */
lin_solver_dcomplex_parameters::lin_solver_dcomplex_parameters (int matrix_layout_i,
   char uplo_i, int n_i, int nrhs_i, int lda_i, int ldb_i)
{
	int j;

	matrix_layout = matrix_layout_i;
	uplo = uplo_i;
	n = n_i;
	nrhs = nrhs_i;
	lda = lda_i;
	ldb = ldb_i;
	
/* Memory allocation of the buffers */
    a = (lapack_complex_double *)malloc(n*n*sizeof(lapack_complex_double));
	aref = (lapack_complex_double *)malloc(n*n*sizeof(lapack_complex_double));
   
    ipiv    = (int *)malloc(n*sizeof(int));
	ipivref = (int *)malloc(n*sizeof(int));
  
	b    = (lapack_complex_double *)malloc((ldb*nrhs)*sizeof(lapack_complex_double));
	bref = (lapack_complex_double *)malloc((ldb*nrhs)*sizeof(lapack_complex_double));

	
	if( (a==NULL) || (aref==NULL) || (b==NULL) || (bref==NULL)|| \
	    (ipiv==NULL) || (ipivref==NULL) ){

       printf(" lin_solver_dcomplex_parameters object: malloc error. Exiting...\n");
       free (a   );
       free (aref);
       free (b   );
       free (bref);
       free (ipiv   );
       free (ipivref);
	   exit (-1);
	}

    /* Initialization of input Buffers */
	  
    for( j = 0; j < n*n; j++ ) {
	   float dtemp_r, dtemp_i;
	   dtemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   dtemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   a[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
       aref[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
    }

    for( j = 0; j < (ldb*nrhs); j++ ) {
	   float dtemp_r, dtemp_i;
	   dtemp_r =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   dtemp_i =(((float) rand()) / ((float) RAND_MAX) - 0.5);
	   b[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
       bref[j] = lapack_make_complex_double( dtemp_r, dtemp_i);
    }
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
	
	info = 0;
	inforef = 0;
	
}

/* Destructor definition **/
lin_solver_dcomplex_parameters:: ~lin_solver_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
       printf(" lin_solver_dcomplex_parameters object: destructor invoked. \n");
#endif
       free (a   );
       free (aref);
       free (b   );
       free (bref);
       free (ipiv   );
       free (ipivref);
}
