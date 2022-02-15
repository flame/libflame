#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#include <math.h>
#include "lapacke.h"

#define sycon_3_free() \
       free (a   ); \
       free (aref); \
       free (e   ); \
       free (eref); \
       free (ipiv   ); \
       free (ipivref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

class sycon_3_double_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   double  *a, *aref;
   double  *e, *eref;
   int 	lda;
   int     *ipiv, *ipivref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      sycon_3_double_parameters (int matrix_layout, char uplo, int n);
      ~sycon_3_double_parameters ();

}; /* end of sycon_3_double_parameters  class definition */

/* Destructor definition **/
sycon_3_double_parameters:: ~sycon_3_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sycon_3_double_parameters object: destructor invoked. \n");
#endif
   sycon_3_free();
}

sycon_3_double_parameters::sycon_3_double_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &a,  &aref,  (lda*n));
   lapacke_gtest_alloc_double_buffer_pair( &e,  &eref,  n);
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (e==NULL) || (eref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       sycon_3_free();
       printf(" sycon_3_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, (lda*n));
    lapacke_gtest_init_double_buffer_pair_rand( e, eref, n);
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
	   e[j] = 0;
	   eref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class dsycon_3_test  : public  ::testing::Test {
public:
   sycon_3_double_parameters  *dsycon_3_obj;
   void SetUp();
   void TearDown () { delete dsycon_3_obj; }
};

void dsycon_3_test::SetUp(){

   /* LAPACKE DSYCON_3 prototype */
   typedef int (*Fptr_NL_LAPACKE_dsycon_3) ( int matrix_layout, char uplo, 
             lapack_int n, const double* a, lapack_int lda,
             const double* e,
			 const lapack_int* ipiv, double anorm, double* rcond );
						   
   Fptr_NL_LAPACKE_dsycon_3 DSYCON_3;
   double diff;
   void *hModule, *dModule;

	
   dsycon_3_obj = new sycon_3_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   dsycon_3_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   DSYCON_3 = (Fptr_NL_LAPACKE_dsycon_3)dlsym(hModule, "LAPACKE_dsycon_3");
   if (NULL == DSYCON_3)
   {
   	  printf("Could not get the symbol -LAPACKE_dsycon_3- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   dsycon_3_obj->info =  LAPACKE_dsytrf_rk  (dsycon_3_obj->matrix_layout, dsycon_3_obj->uplo, 
                                 dsycon_3_obj->n,
                                 dsycon_3_obj->a, 
								 dsycon_3_obj->lda, 
                                 dsycon_3_obj->e, 
					             dsycon_3_obj->ipiv );

   dsycon_3_obj->inforef =  LAPACKE_dsytrf_rk  (dsycon_3_obj->matrix_layout, dsycon_3_obj->uplo, 
                                 dsycon_3_obj->n,
                                 dsycon_3_obj->aref, 
								 dsycon_3_obj->lda, 
                                 dsycon_3_obj->eref, 
					             dsycon_3_obj->ipivref );

   /* Compute libflame's Lapacke o/p  */  
   dsycon_3_obj->info    = LAPACKE_dsycon_3( dsycon_3_obj->matrix_layout, dsycon_3_obj->uplo, 
                                 dsycon_3_obj->n,
                                 (const double *)dsycon_3_obj->a, 
								 dsycon_3_obj->lda, 
                                 (const double *)dsycon_3_obj->e, 
					             (const lapack_int *)dsycon_3_obj->ipiv, 
								 dsycon_3_obj->anorm, &dsycon_3_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   dsycon_3_obj->inforef = DSYCON_3( dsycon_3_obj->matrix_layout, dsycon_3_obj->uplo, 
                                dsycon_3_obj->n, 
				                (const double *)dsycon_3_obj->aref, 
								dsycon_3_obj->lda,
								(const double* )dsycon_3_obj->eref, 
			                    (const lapack_int *)dsycon_3_obj->ipivref, 
								dsycon_3_obj->anorm, &dsycon_3_obj->rcondref);

	if( dsycon_3_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dsycon_3 is wrong\n", 
		         dsycon_3_obj->info );
	}
	if( dsycon_3_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dsycon_3 is wrong\n", 
		dsycon_3_obj->inforef );
	}

    printf("dsycon_3_obj->rcond: %lf   dsycon_3_obj->rcondref: %lf \n", dsycon_3_obj->rcond, dsycon_3_obj->rcondref);
	EXPECT_NEAR(0.0, abs(dsycon_3_obj->rcond - dsycon_3_obj->rcondref), dsycon_3_obj->threshold);

}

TEST(dsycon_3_test, dsycon_3_1) {}
TEST(dsycon_3_test, dsycon_3_2) {}
TEST(dsycon_3_test, dsycon_3_3) {}
TEST(dsycon_3_test, dsycon_3_4) {}

class sycon_3_float_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   float  *a, *aref;
   float  *e, *eref;
   int 	lda;
   int     *ipiv, *ipivref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      sycon_3_float_parameters (int matrix_layout, char uplo, int n);
      ~sycon_3_float_parameters ();

}; /* end of sycon_3_float_parameters  class definition */

/* Destructor definition **/
sycon_3_float_parameters:: ~sycon_3_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sycon_3_float_parameters object: destructor invoked. \n");
#endif
   sycon_3_free();
}

sycon_3_float_parameters::sycon_3_float_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &a,  &aref,  (lda*n));
   lapacke_gtest_alloc_float_buffer_pair( &e,  &eref,  n);
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (e==NULL) || (eref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       sycon_3_free();
       printf(" sycon_3_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, (lda*n));
    lapacke_gtest_init_float_buffer_pair_rand( e, eref, n);
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
	   e[j] = 0;
	   eref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class ssycon_3_test  : public  ::testing::Test {
public:
   sycon_3_float_parameters  *ssycon_3_obj;
   void SetUp();
   void TearDown () { delete ssycon_3_obj; }
};

void ssycon_3_test::SetUp(){

   /* LAPACKE SSYCON_3 prototype */
   typedef int (*Fptr_NL_LAPACKE_ssycon_3) ( int matrix_layout, char uplo, 
             lapack_int n, const float* a, lapack_int lda,
             const float* e,
			 const lapack_int* ipiv, float anorm, float* rcond );
						   
   Fptr_NL_LAPACKE_ssycon_3 SSYCON_3;
   float diff;
   void *hModule, *dModule;

	
   ssycon_3_obj = new sycon_3_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   ssycon_3_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   SSYCON_3 = (Fptr_NL_LAPACKE_ssycon_3)dlsym(hModule, "LAPACKE_ssycon_3");
   if (NULL == SSYCON_3)
   {
   	  printf("Could not get the symbol -LAPACKE_ssycon_3- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   ssycon_3_obj->info =  LAPACKE_ssytrf_rk  (ssycon_3_obj->matrix_layout, ssycon_3_obj->uplo, 
                                 ssycon_3_obj->n,
                                 ssycon_3_obj->a, 
								 ssycon_3_obj->lda, 
                                 ssycon_3_obj->e, 
					             ssycon_3_obj->ipiv );

   ssycon_3_obj->inforef =  LAPACKE_ssytrf_rk  (ssycon_3_obj->matrix_layout, ssycon_3_obj->uplo, 
                                 ssycon_3_obj->n,
                                 ssycon_3_obj->aref, 
								 ssycon_3_obj->lda, 
                                 ssycon_3_obj->eref, 
					             ssycon_3_obj->ipivref );

   /* Compute libflame's Lapacke o/p  */  
   ssycon_3_obj->info    = LAPACKE_ssycon_3( ssycon_3_obj->matrix_layout, ssycon_3_obj->uplo, 
                                 ssycon_3_obj->n,
                                 (const float *)ssycon_3_obj->a, 
								 ssycon_3_obj->lda, 
                                 (const float *)ssycon_3_obj->e, 
					             (const lapack_int *)ssycon_3_obj->ipiv, 
								 ssycon_3_obj->anorm, &ssycon_3_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   ssycon_3_obj->inforef = SSYCON_3( ssycon_3_obj->matrix_layout, ssycon_3_obj->uplo, 
                                ssycon_3_obj->n, 
				                (const float *)ssycon_3_obj->aref, 
								ssycon_3_obj->lda,
								(const float* )ssycon_3_obj->eref, 
			                    (const lapack_int *)ssycon_3_obj->ipivref, 
								ssycon_3_obj->anorm, &ssycon_3_obj->rcondref);

	if( ssycon_3_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_ssycon_3 is wrong\n", 
		         ssycon_3_obj->info );
	}
	if( ssycon_3_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_ssycon_3 is wrong\n", 
		ssycon_3_obj->inforef );
	}

    printf("ssycon_3_obj->rcond: %lf   ssycon_3_obj->rcondref: %lf \n", ssycon_3_obj->rcond, ssycon_3_obj->rcondref);
	EXPECT_NEAR(0.0, abs(ssycon_3_obj->rcond - ssycon_3_obj->rcondref), ssycon_3_obj->threshold);

}

TEST(ssycon_3_test, ssycon_3_1) {}
TEST(ssycon_3_test, ssycon_3_2) {}
TEST(ssycon_3_test, ssycon_3_3) {}
TEST(ssycon_3_test, ssycon_3_4) {}


class sycon_3_scomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   lapack_complex_float  *a, *aref;
   lapack_complex_float  *e, *eref;
   int 	lda;
   int     *ipiv, *ipivref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      sycon_3_scomplex_parameters (int matrix_layout, char uplo, int n);
      ~sycon_3_scomplex_parameters ();

}; /* end of sycon_3_float_parameters  class definition */

/* Destructor definition **/
sycon_3_scomplex_parameters:: ~sycon_3_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sycon_3_float_parameters object: destructor invoked. \n");
#endif
   sycon_3_free();
}

sycon_3_scomplex_parameters::sycon_3_scomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a,  &aref,  (lda*n));
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e,  &eref,  n);
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (e==NULL) || (eref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       sycon_3_free();
       printf(" sycon_3_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (lda*n));
    lapacke_gtest_init_scomplex_buffer_pair_rand( e, eref, n);
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
	   e[j] = 0;
	   eref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class csycon_3_test  : public  ::testing::Test {
public:
   sycon_3_scomplex_parameters  *csycon_3_obj;
   void SetUp();
   void TearDown () { delete csycon_3_obj; }
};

void csycon_3_test::SetUp(){

   /* LAPACKE CSYCON_3 prototype */
   typedef int (*Fptr_NL_LAPACKE_csycon_3) ( int matrix_layout, char uplo, 
             lapack_int n, const lapack_complex_float* a, lapack_int lda,
             const lapack_complex_float* e,
			 const lapack_int* ipiv, float anorm, float* rcond );
						   
   Fptr_NL_LAPACKE_csycon_3 CSYCON_3;
   float diff;
   void *hModule, *dModule;

	
   csycon_3_obj = new sycon_3_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   csycon_3_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   CSYCON_3 = (Fptr_NL_LAPACKE_csycon_3)dlsym(hModule, "LAPACKE_csycon_3");
   if (NULL == CSYCON_3)
   {
   	  printf("Could not get the symbol -LAPACKE_csycon_3- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   csycon_3_obj->info =  LAPACKE_csytrf_rk  (csycon_3_obj->matrix_layout, csycon_3_obj->uplo, 
                                 csycon_3_obj->n,
                                 csycon_3_obj->a, 
								 csycon_3_obj->lda, 
                                 csycon_3_obj->e, 
					             csycon_3_obj->ipiv );

   csycon_3_obj->inforef =  LAPACKE_csytrf_rk  (csycon_3_obj->matrix_layout, csycon_3_obj->uplo, 
                                 csycon_3_obj->n,
                                 csycon_3_obj->aref, 
								 csycon_3_obj->lda, 
                                 csycon_3_obj->eref, 
					             csycon_3_obj->ipivref );

   /* Compute libflame's Lapacke o/p  */  
   csycon_3_obj->info    = LAPACKE_csycon_3( csycon_3_obj->matrix_layout, csycon_3_obj->uplo, 
                                 csycon_3_obj->n,
                                 (const lapack_complex_float *)csycon_3_obj->a, 
								 csycon_3_obj->lda, 
                                 (const lapack_complex_float *)csycon_3_obj->e, 
					             (const lapack_int *)csycon_3_obj->ipiv, 
								 csycon_3_obj->anorm, &csycon_3_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   csycon_3_obj->inforef = CSYCON_3( csycon_3_obj->matrix_layout, csycon_3_obj->uplo, 
                                csycon_3_obj->n, 
				                (const lapack_complex_float *)csycon_3_obj->aref, 
								csycon_3_obj->lda,
								(const lapack_complex_float* )csycon_3_obj->eref, 
			                    (const lapack_int *)csycon_3_obj->ipivref, 
								csycon_3_obj->anorm, &csycon_3_obj->rcondref);

	if( csycon_3_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_csycon_3 is wrong\n", 
		         csycon_3_obj->info );
	}
	if( csycon_3_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_csycon_3 is wrong\n", 
		csycon_3_obj->inforef );
	}

    printf("csycon_3_obj->rcond: %lf   csycon_3_obj->rcondref: %lf \n", csycon_3_obj->rcond, csycon_3_obj->rcondref);
	EXPECT_NEAR(0.0, abs(csycon_3_obj->rcond - csycon_3_obj->rcondref), csycon_3_obj->threshold);

}

TEST(csycon_3_test, csycon_3_1) {}
TEST(csycon_3_test, csycon_3_2) {}
TEST(csycon_3_test, csycon_3_3) {}
TEST(csycon_3_test, csycon_3_4) {}

class sycon_3_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   lapack_complex_double  *a, *aref;
   lapack_complex_double  *e, *eref;
   int 	lda;
   int     *ipiv, *ipivref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      sycon_3_dcomplex_parameters (int matrix_layout, char uplo, int n);
      ~sycon_3_dcomplex_parameters ();

}; /* end of sycon_3_double_parameters  class definition */

/* Destructor definition **/
sycon_3_dcomplex_parameters:: ~sycon_3_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sycon_3_double_parameters object: destructor invoked. \n");
#endif
   sycon_3_free();
}

sycon_3_dcomplex_parameters::sycon_3_dcomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a,  &aref,  (lda*n));
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e,  &eref,  n);
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (e==NULL) || (eref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       sycon_3_free();
       printf(" sycon_3_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, (lda*n));
    lapacke_gtest_init_dcomplex_buffer_pair_rand( e, eref, n);
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
	   e[j] = 0;
	   eref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class zsycon_3_test  : public  ::testing::Test {
public:
   sycon_3_dcomplex_parameters  *zsycon_3_obj;
   void SetUp();
   void TearDown () { delete zsycon_3_obj; }
};

void zsycon_3_test::SetUp(){

   /* LAPACKE ZSYCON_3 prototype */
   typedef int (*Fptr_NL_LAPACKE_zsycon_3) ( int matrix_layout, char uplo, 
             lapack_int n, const lapack_complex_double* a, lapack_int lda,
             const lapack_complex_double* e,
			 const lapack_int* ipiv, double anorm, double* rcond );
						   
   Fptr_NL_LAPACKE_zsycon_3 ZSYCON_3;
   double diff;
   void *hModule, *dModule;

	
   zsycon_3_obj = new sycon_3_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   zsycon_3_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   ZSYCON_3 = (Fptr_NL_LAPACKE_zsycon_3)dlsym(hModule, "LAPACKE_zsycon_3");
   if (NULL == ZSYCON_3)
   {
   	  printf("Could not get the symbol -LAPACKE_zsycon_3- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   zsycon_3_obj->info =  LAPACKE_zsytrf_rk  (zsycon_3_obj->matrix_layout, zsycon_3_obj->uplo, 
                                 zsycon_3_obj->n,
                                 zsycon_3_obj->a, 
								 zsycon_3_obj->lda, 
                                 zsycon_3_obj->e, 
					             zsycon_3_obj->ipiv );

   zsycon_3_obj->inforef =  LAPACKE_zsytrf_rk  (zsycon_3_obj->matrix_layout, zsycon_3_obj->uplo, 
                                 zsycon_3_obj->n,
                                 zsycon_3_obj->aref, 
								 zsycon_3_obj->lda, 
                                 zsycon_3_obj->eref, 
					             zsycon_3_obj->ipivref );

   /* Compute libflame's Lapacke o/p  */  
   zsycon_3_obj->info    = LAPACKE_zsycon_3( zsycon_3_obj->matrix_layout, zsycon_3_obj->uplo, 
                                 zsycon_3_obj->n,
                                 (const lapack_complex_double *)zsycon_3_obj->a, 
								 zsycon_3_obj->lda, 
                                 (const lapack_complex_double *)zsycon_3_obj->e, 
					             (const lapack_int *)zsycon_3_obj->ipiv, 
								 zsycon_3_obj->anorm, &zsycon_3_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   zsycon_3_obj->inforef = ZSYCON_3( zsycon_3_obj->matrix_layout, zsycon_3_obj->uplo, 
                                zsycon_3_obj->n, 
				                (const lapack_complex_double *)zsycon_3_obj->aref, 
								zsycon_3_obj->lda,
								(const lapack_complex_double* )zsycon_3_obj->eref, 
			                    (const lapack_int *)zsycon_3_obj->ipivref, 
								zsycon_3_obj->anorm, &zsycon_3_obj->rcondref);

	if( zsycon_3_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_zsycon_3 is wrong\n", 
		         zsycon_3_obj->info );
	}
	if( zsycon_3_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_zsycon_3 is wrong\n", 
		zsycon_3_obj->inforef );
	}

    printf("zsycon_3_obj->rcond: %lf   zsycon_3_obj->rcondref: %lf \n", zsycon_3_obj->rcond, zsycon_3_obj->rcondref);
	EXPECT_NEAR(0.0, abs(zsycon_3_obj->rcond - zsycon_3_obj->rcondref), zsycon_3_obj->threshold);

}

TEST(zsycon_3_test, zsycon_3_1) {}
TEST(zsycon_3_test, zsycon_3_2) {}
TEST(zsycon_3_test, zsycon_3_3) {}
TEST(zsycon_3_test, zsycon_3_4) {}
