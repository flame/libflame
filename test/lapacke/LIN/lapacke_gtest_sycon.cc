#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#include <math.h>
#include "lapacke.h"

#define sycon_free() \
       free (a   ); \
       free (aref); \
       free (ipiv   ); \
       free (ipivref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

class sycon_double_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;   
   char uplo;
   int 	n;
   double  *a, *aref;
   int 	lda;
   int     *ipiv, *ipivref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      sycon_double_parameters (int matrix_layout, char uplo, int n);
      ~sycon_double_parameters ();

}; /* end of sycon_double_parameters  class definition */

/* Destructor definition **/
sycon_double_parameters:: ~sycon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sycon_double_parameters object: destructor invoked. \n");
#endif
   sycon_free();
}

sycon_double_parameters::sycon_double_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &a,  &aref,  (lda*n));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       sycon_free();
       printf(" sycon_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, (lda*n));
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class dsycon_test  : public  ::testing::Test {
public:
   sycon_double_parameters  *dsycon_obj;
   void SetUp();
   void TearDown () { delete dsycon_obj; }
};

void dsycon_test::SetUp(){

   /* LAPACKE DSYCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dsycon) ( int matrix_layout, char uplo, 
                  lapack_int n, const double* a, lapack_int lda, 
			      const lapack_int* ipiv, double anorm, double* rcond );
						   
   Fptr_NL_LAPACKE_dsycon DSYCON;
   double diff;
   void *hModule, *dModule;

	
   dsycon_obj = new sycon_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   dsycon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   DSYCON = (Fptr_NL_LAPACKE_dsycon)dlsym(hModule, "LAPACKE_dsycon");
   if (NULL == DSYCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dsycon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   dsycon_obj->info =  LAPACKE_dsytrf (dsycon_obj->matrix_layout, dsycon_obj->uplo, 
                                 dsycon_obj->n,
                                 dsycon_obj->a, 
								 dsycon_obj->lda, 
					             dsycon_obj->ipiv );

   dsycon_obj->inforef =  LAPACKE_dsytrf  (dsycon_obj->matrix_layout, dsycon_obj->uplo, 
                                 dsycon_obj->n,
                                 dsycon_obj->aref, 
								 dsycon_obj->lda, 
					             dsycon_obj->ipivref );

   /* Compute libflame's Lapacke o/p  */
   dsycon_obj->info    = LAPACKE_dsycon( dsycon_obj->matrix_layout, dsycon_obj->uplo, 
                                        dsycon_obj->n, 
                                        (const double * )dsycon_obj->a, dsycon_obj->lda, 
										(const lapack_int *)dsycon_obj->ipiv, 
										dsycon_obj->anorm, &dsycon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   dsycon_obj->inforef = DSYCON( dsycon_obj->matrix_layout, dsycon_obj->uplo, 
                                dsycon_obj->n,  
								(const double * )dsycon_obj->aref, 
								dsycon_obj->lda, 
								(const lapack_int *)dsycon_obj->ipivref, 
								dsycon_obj->anorm, &dsycon_obj->rcondref);

	if( dsycon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dsycon is wrong\n", 
		         dsycon_obj->info );
	}
	if( dsycon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dsycon is wrong\n", 
		dsycon_obj->inforef );
	}

    printf("dsycon_obj->rcond: %lf   dsycon_obj->rcondref: %lf \n", dsycon_obj->rcond, dsycon_obj->rcondref);
	EXPECT_NEAR(0.0, abs(dsycon_obj->rcond - dsycon_obj->rcondref), DOUBLE_DIFF_THRESHOLD);
}

TEST(dsycon_test, dsycon_1) {}
TEST(dsycon_test, dsycon_2) {}
TEST(dsycon_test, dsycon_3) {}
TEST(dsycon_test, dsycon_4) {}

class sycon_float_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   
   char uplo;
   int 	n;
   
   float  *a, *aref;
   int 	lda;
   int     *ipiv, *ipivref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      sycon_float_parameters (int matrix_layout, char uplo, int n);
      ~sycon_float_parameters ();

}; /* end of sycon_float_parameters  class definition */

/* Destructor definition **/
sycon_float_parameters:: ~sycon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sycon_float_parameters object: destructor invoked. \n");
#endif
   sycon_free();
}

sycon_float_parameters::sycon_float_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &a,  &aref,  (lda*n));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       sycon_free();
       printf(" sycon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, (lda*n));
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class ssycon_test  : public  ::testing::Test {
public:
   sycon_float_parameters  *ssycon_obj;
   void SetUp();
   void TearDown () { delete ssycon_obj; }
};

void ssycon_test::SetUp(){

   /* LAPACKE SSYCON prototype */
   typedef int (*Fptr_NL_LAPACKE_ssycon) ( int matrix_layout, char uplo, 
                                           lapack_int n, const float* a, 
										   lapack_int lda, 
										   const lapack_int* ipiv, 
										   float anorm, float* rcond );
						   
   Fptr_NL_LAPACKE_ssycon SSYCON;
   float diff;
   void *hModule, *dModule;

	
   ssycon_obj = new sycon_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   ssycon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   SSYCON = (Fptr_NL_LAPACKE_ssycon)dlsym(hModule, "LAPACKE_ssycon");
   if (NULL == SSYCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dsycon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   ssycon_obj->info =  LAPACKE_ssytrf  (ssycon_obj->matrix_layout, ssycon_obj->uplo, 
                                 ssycon_obj->n,
                                 ssycon_obj->a, 
								 ssycon_obj->lda, 
					             ssycon_obj->ipiv );

   ssycon_obj->inforef =  LAPACKE_ssytrf  (ssycon_obj->matrix_layout, ssycon_obj->uplo, 
                                 ssycon_obj->n,
                                 ssycon_obj->aref, 
								 ssycon_obj->lda, 
					             ssycon_obj->ipivref );
   /* Compute libflame's Lapacke o/p  */
   ssycon_obj->info    = LAPACKE_ssycon( ssycon_obj->matrix_layout, ssycon_obj->uplo, 
                                        ssycon_obj->n, 
                                        (const float*)ssycon_obj->a, ssycon_obj->lda, 
										(const lapack_int *)ssycon_obj->ipiv, 
										ssycon_obj->anorm, &ssycon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   ssycon_obj->inforef = SSYCON( ssycon_obj->matrix_layout, ssycon_obj->uplo, 
                                ssycon_obj->n,  
								(const float*)ssycon_obj->aref, ssycon_obj->lda,
								(const lapack_int *)ssycon_obj->ipivref, 
								ssycon_obj->anorm, &ssycon_obj->rcondref);
	if( ssycon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dsycon is wrong\n", 
		         ssycon_obj->info );
	}
	if( ssycon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dsycon is wrong\n", 
		ssycon_obj->inforef );
	}

    printf("ssycon_obj->rcond: %lf   ssycon_obj->rcondref: %lf \n", ssycon_obj->rcond, ssycon_obj->rcondref);
	EXPECT_NEAR(0.0, abs(ssycon_obj->rcond - ssycon_obj->rcondref), DOUBLE_DIFF_THRESHOLD);
}

TEST(ssycon_test, ssycon_1) {}
TEST(ssycon_test, ssycon_2) {}
TEST(ssycon_test, ssycon_3) {}
TEST(ssycon_test, ssycon_4) {}

class sycon_scomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   lapack_complex_float  *a, *aref;
   int 	lda;
   int     *ipiv, *ipivref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      sycon_scomplex_parameters (int matrix_layout, char uplo, int n);
      ~sycon_scomplex_parameters ();

}; /* end of sycon_float_parameters  class definition */

/* Destructor definition **/
sycon_scomplex_parameters:: ~sycon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sycon_float_parameters object: destructor invoked. \n");
#endif
   sycon_free();
}

sycon_scomplex_parameters::sycon_scomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a,  &aref,  (lda*n));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       sycon_free();
       printf(" sycon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (lda*n));
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class csycon_test  : public  ::testing::Test {
public:
   sycon_scomplex_parameters  *csycon_obj;
   void SetUp();
   void TearDown () { delete csycon_obj; }
};

void csycon_test::SetUp(){

   /* LAPACKE CSYCON prototype */
   typedef int (*Fptr_NL_LAPACKE_csycon) ( int matrix_layout, char uplo, 
             lapack_int n, const lapack_complex_float* a, lapack_int lda, 
			 const lapack_int* ipiv, float anorm, float* rcond );
						   
   Fptr_NL_LAPACKE_csycon CSYCON;
   float diff;
   void *hModule, *dModule;

	
   csycon_obj = new sycon_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   csycon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   CSYCON = (Fptr_NL_LAPACKE_csycon)dlsym(hModule, "LAPACKE_csycon");
   if (NULL == CSYCON)
   {
   	  printf("Could not get the symbol -LAPACKE_csycon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   csycon_obj->info =  LAPACKE_csytrf  (csycon_obj->matrix_layout, csycon_obj->uplo, 
                                 csycon_obj->n,
                                 csycon_obj->a, 
								 csycon_obj->lda, 
					             csycon_obj->ipiv );

   csycon_obj->inforef =  LAPACKE_csytrf  (csycon_obj->matrix_layout, csycon_obj->uplo, 
                                 csycon_obj->n,
                                 csycon_obj->aref, 
								 csycon_obj->lda, 
					             csycon_obj->ipivref );

   /* Compute libflame's Lapacke o/p  */  
   csycon_obj->info    = LAPACKE_csycon( csycon_obj->matrix_layout, csycon_obj->uplo, 
                                        csycon_obj->n,
                                        (const lapack_complex_float*)csycon_obj->a, 
										csycon_obj->lda, 
					                    (const lapack_int *)csycon_obj->ipiv, csycon_obj->anorm,
			                            &csycon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   csycon_obj->inforef = CSYCON( csycon_obj->matrix_layout, csycon_obj->uplo, 
                                csycon_obj->n, 
				                (const lapack_complex_float*)csycon_obj->aref, 
								csycon_obj->lda, 
			                    (const lapack_int *)csycon_obj->ipivref, csycon_obj->anorm,
				                &csycon_obj->rcondref);

	if( csycon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_csycon is wrong\n", 
		         csycon_obj->info );
	}
	if( csycon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_csycon is wrong\n", 
		csycon_obj->inforef );
	}

    printf("csycon_obj->rcond: %lf   csycon_obj->rcondref: %lf \n", csycon_obj->rcond, csycon_obj->rcondref);
	EXPECT_NEAR(0.0, abs(csycon_obj->rcond - csycon_obj->rcondref), DOUBLE_DIFF_THRESHOLD);

}

TEST(csycon_test, csycon_1) {}
TEST(csycon_test, csycon_2) {}
TEST(csycon_test, csycon_3) {}
TEST(csycon_test, csycon_4) {}

class sycon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   
   lapack_complex_double  *a, *aref;
   int 	lda;
   int     *ipiv, *ipivref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      sycon_dcomplex_parameters (int matrix_layout, char uplo, int n);
      ~sycon_dcomplex_parameters ();

}; /* end of sycon_float_parameters  class definition */

/* Destructor definition **/
sycon_dcomplex_parameters:: ~sycon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sycon_float_parameters object: destructor invoked. \n");
#endif
   sycon_free();
}

sycon_dcomplex_parameters::sycon_dcomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a,  &aref,  (lda*n));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       sycon_free();
       printf(" sycon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, (lda*n));
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class zsycon_test  : public  ::testing::Test {
public:
   sycon_dcomplex_parameters  *zsycon_obj;
   void SetUp();
   void TearDown () { delete zsycon_obj; }
};

void zsycon_test::SetUp(){

   /* LAPACKE ZSYCON prototype */
   typedef int (*Fptr_NL_LAPACKE_zsycon) ( int matrix_layout, char uplo, 
                          lapack_int n, const lapack_complex_double* a, 
						  lapack_int lda, const lapack_int* ipiv, 
						  double anorm, double* rcond );
						   
   Fptr_NL_LAPACKE_zsycon ZSYCON;
   double diff;
   void *hModule, *dModule;

	
   zsycon_obj = new sycon_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   zsycon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   ZSYCON = (Fptr_NL_LAPACKE_zsycon)dlsym(hModule, "LAPACKE_zsycon");
   if (NULL == ZSYCON)
   {
   	  printf("Could not get the symbol -LAPACKE_zsycon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   zsycon_obj->info =  LAPACKE_zsytrf  (zsycon_obj->matrix_layout, zsycon_obj->uplo, 
                                 zsycon_obj->n,
                                 zsycon_obj->a, 
								 zsycon_obj->lda, 
					             zsycon_obj->ipiv );

   zsycon_obj->inforef =  LAPACKE_zsytrf  (zsycon_obj->matrix_layout, zsycon_obj->uplo, 
                                 zsycon_obj->n,
                                 zsycon_obj->aref, 
								 zsycon_obj->lda, 
					             zsycon_obj->ipivref );

   /* Compute libflame's Lapacke o/p  */  
   zsycon_obj->info    = LAPACKE_zsycon( zsycon_obj->matrix_layout, zsycon_obj->uplo, 
                                        zsycon_obj->n,
                                        (const lapack_complex_double*)zsycon_obj->a, 
										zsycon_obj->lda, 
					                    (const lapack_int *)zsycon_obj->ipiv, zsycon_obj->anorm,
			                            &zsycon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   zsycon_obj->inforef = ZSYCON( zsycon_obj->matrix_layout, zsycon_obj->uplo, 
                                zsycon_obj->n, 
				                (const lapack_complex_double*)zsycon_obj->aref, 
								zsycon_obj->lda, 
			                    (const lapack_int *)zsycon_obj->ipivref, zsycon_obj->anorm,
				                &zsycon_obj->rcondref);

	if( zsycon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_csycon is wrong\n", 
		         zsycon_obj->info );
	}
	if( zsycon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_csycon is wrong\n", 
		zsycon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
    printf("zsycon_obj->rcond: %lf   zsycon_obj->rcondref: %lf \n", zsycon_obj->rcond, zsycon_obj->rcondref);
	EXPECT_NEAR(0.0, abs(zsycon_obj->rcond - zsycon_obj->rcondref), DOUBLE_DIFF_THRESHOLD);

}

TEST(zsycon_test, zsycon_1) {}
TEST(zsycon_test, zsycon_2) {}
TEST(zsycon_test, zsycon_3) {}
TEST(zsycon_test, zsycon_4) {}
