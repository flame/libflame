#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#include <math.h>
#include "lapacke.h"

#define gtcon_free() \
       free (dl   ); \
       free (dlref); \
       free (d   ); \
       free (dref); \
       free (du   ); \
       free (duref); \
       free (du2   ); \
       free (du2ref); \
       free (ipiv   ); \
       free (ipivref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

class gtcon_double_parameters{
   public:
   /* input params to the API **/
   char norm;
   int 	n;
   double  *dl, *dlref;
   double  *d, *dref;
   double  *du, *duref;
   double  *du2, *du2ref;
   int     *ipiv, *ipivref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   float threshold;

   public: 
      gtcon_double_parameters (char norm, int n );
      ~gtcon_double_parameters ();

}; /* end of gtcon_double_parameters  class definition */

/* Destructor definition **/
gtcon_double_parameters:: ~gtcon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtcon_double_parameters object: destructor invoked. \n");
#endif
   gtcon_free();
}

gtcon_double_parameters::gtcon_double_parameters (char norm_i, int n_i)
{
   int j;
   
   norm = norm_i;
   n = n_i;
   anorm = (double)n/2.0;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &d,  &dref,  n);
   lapacke_gtest_alloc_double_buffer_pair( &dl,  &dlref,  n-1);
   lapacke_gtest_alloc_double_buffer_pair( &du,  &duref,  n-1);
   lapacke_gtest_alloc_double_buffer_pair( &du2,  &du2ref,  n-2);
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (d==NULL) || (dref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (du==NULL) || (duref==NULL) ||  \
        (du2==NULL) || (du2ref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       gtcon_free();
       printf(" gtcon_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_double_buffer_pair_rand( dl, dlref, n-1);
    lapacke_gtest_init_double_buffer_pair_rand( du, duref, n-1);
    lapacke_gtest_init_double_buffer_pair_rand( du2, du2ref, n-2);

    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }

	info = 0;
	inforef = 0;
}


/*  Test fixture class definition */
class dgtcon_test  : public  ::testing::Test {
public:
   gtcon_double_parameters  *dgtcon_obj;
   void SetUp();
   void TearDown () { delete dgtcon_obj; }
};

void dgtcon_test::SetUp() {

   /* LAPACKE DGTCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dgtcon) ( char norm, lapack_int n, const double *dl,
                          const double *d, const double *du, const double *du2, 
				          const lapack_int *ipiv, double anorm, double *rcond );
						   
   Fptr_NL_LAPACKE_dgtcon DGTCON;
   double diff;
   void *hModule, *dModule;

   dgtcon_obj = new gtcon_double_parameters (
					lin_solver_paramslist[idx].norm_gbcon,
                    lin_solver_paramslist[idx].n );

   dgtcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   DGTCON = (Fptr_NL_LAPACKE_dgtcon)dlsym(hModule, "LAPACKE_dgtcon");
   if (NULL == DGTCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dgtcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   dgtcon_obj->info    = LAPACKE_dgttrf( dgtcon_obj->n, 
                                         dgtcon_obj->dl,
                                         dgtcon_obj->d, 
										 dgtcon_obj->du, 
										 dgtcon_obj->du2,
										dgtcon_obj->ipiv );

   dgtcon_obj->inforef = LAPACKE_dgttrf( dgtcon_obj->n, 
                                         dgtcon_obj->dlref,
                                         dgtcon_obj->dref, 
										 dgtcon_obj->duref, 
										 dgtcon_obj->du2ref,
										dgtcon_obj->ipivref );
   
   /* Compute libflame's Lapacke o/p  */
   dgtcon_obj->info    = LAPACKE_dgtcon( dgtcon_obj->norm, dgtcon_obj->n, 
                                        (const double *)dgtcon_obj->dl,
                                        (const double *)dgtcon_obj->d, 
										(const double *)dgtcon_obj->du, 
										(const double *)dgtcon_obj->du2,
										dgtcon_obj->ipiv, dgtcon_obj->anorm,
						                &dgtcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
    dgtcon_obj->inforef = DGTCON( dgtcon_obj->norm, dgtcon_obj->n, 
	                            (const double *)dgtcon_obj->dlref, 
                                (const double *)dgtcon_obj->dref, 
								(const double *)dgtcon_obj->duref, 
								(const double *)dgtcon_obj->du2ref, 
								dgtcon_obj->ipivref, dgtcon_obj->anorm,
								&dgtcon_obj->rcondref);

	if( dgtcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument to libflame LAPACKE_dgtcon call is wrong\n", 
		         dgtcon_obj->info );
	}
	if( dgtcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument to Netlib LAPACKE_dgtcon call is wrong\n", 
		dgtcon_obj->inforef );
	}

#if LAPACKE_TEST_VERBOSE
    printf("dgtcon_obj->rcond: %lf   dgtcon_obj->rcondref: %lf \n", dgtcon_obj->rcond, dgtcon_obj->rcondref);
#endif
	EXPECT_NEAR(0.0, abs(dgtcon_obj->rcond - dgtcon_obj->rcondref), dgtcon_obj->threshold);
}

TEST(dgtcon_test, dgtcon1) {}
TEST(dgtcon_test, dgtcon2) {}
TEST(dgtcon_test, dgtcon3) {}
TEST(dgtcon_test, dgtcon4) {}

class gtcon_float_parameters{
   public:
   /* input params to the API **/
   char norm;
   int 	n;
   float  *dl, *dlref;
   float  *d, *dref;
   float  *du, *duref;
   float  *du2, *du2ref;
   int     *ipiv, *ipivref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   float threshold;

   public: 
      gtcon_float_parameters (char norm, int n );
      ~gtcon_float_parameters ();

}; /* end of gtcon_float_parameters  class definition */

/* Destructor definition **/
gtcon_float_parameters:: ~gtcon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtcon_float_parameters object: destructor invoked. \n");
#endif
   gtcon_free();
}

gtcon_float_parameters::gtcon_float_parameters (char norm_i, int n_i)
{
   int j;
   
   norm = norm_i;
   n = n_i;
   anorm = (float)n/2.0;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &d,  &dref,  n);
   lapacke_gtest_alloc_float_buffer_pair( &dl,  &dlref,  n-1);
   lapacke_gtest_alloc_float_buffer_pair( &du,  &duref,  n-1);
   lapacke_gtest_alloc_float_buffer_pair( &du2,  &du2ref,  n-2);
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (d==NULL) || (dref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (du==NULL) || (duref==NULL) ||  \
        (du2==NULL) || (du2ref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       gtcon_free();
       printf(" gtcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_float_buffer_pair_rand( dl, dlref, n-1);
    lapacke_gtest_init_float_buffer_pair_rand( du, duref, n-1);
    lapacke_gtest_init_float_buffer_pair_rand( du2, du2ref, n-2);

    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }

	info = 0;
	inforef = 0;
}


/*  Test fixture class definition */
class sgtcon_test  : public  ::testing::Test {
public:
   gtcon_float_parameters  *sgtcon_obj;
   void SetUp();
   void TearDown () { delete sgtcon_obj; }
};

void sgtcon_test::SetUp() {

   /* LAPACKE sgtCON prototype */
   typedef int (*Fptr_NL_LAPACKE_sgtcon) ( char norm, lapack_int n, const float *dl,
                          const float *d, const float *du, const float *du2, 
				          const lapack_int *ipiv, float anorm, float *rcond );
						   
   Fptr_NL_LAPACKE_sgtcon SGTCON;
   float diff;
   void *hModule, *dModule;

   sgtcon_obj = new gtcon_float_parameters (
					lin_solver_paramslist[idx].norm_gbcon,
                    lin_solver_paramslist[idx].n );

   sgtcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   SGTCON = (Fptr_NL_LAPACKE_sgtcon)dlsym(hModule, "LAPACKE_sgtcon");
   if (NULL == SGTCON)
   {
   	  printf("Could not get the symbol -LAPACKE_sgtcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   sgtcon_obj->info    = LAPACKE_sgttrf( sgtcon_obj->n, 
                                         sgtcon_obj->dl,
                                         sgtcon_obj->d, 
										 sgtcon_obj->du, 
										 sgtcon_obj->du2,
										sgtcon_obj->ipiv );

   sgtcon_obj->inforef = LAPACKE_sgttrf( sgtcon_obj->n, 
                                         sgtcon_obj->dlref,
                                         sgtcon_obj->dref, 
										 sgtcon_obj->duref, 
										 sgtcon_obj->du2ref,
										sgtcon_obj->ipivref );
   
   /* Compute libflame's Lapacke o/p  */
   sgtcon_obj->info    = LAPACKE_sgtcon( sgtcon_obj->norm, sgtcon_obj->n, 
                                        (const float *)sgtcon_obj->dl,
                                        (const float *)sgtcon_obj->d, 
										(const float *)sgtcon_obj->du, 
										(const float *)sgtcon_obj->du2,
										sgtcon_obj->ipiv, sgtcon_obj->anorm,
						                &sgtcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
    sgtcon_obj->inforef = SGTCON( sgtcon_obj->norm, sgtcon_obj->n, 
	                            (const float *)sgtcon_obj->dlref, 
                                (const float *)sgtcon_obj->dref, 
								(const float *)sgtcon_obj->duref, 
								(const float *)sgtcon_obj->du2ref, 
								sgtcon_obj->ipivref, sgtcon_obj->anorm,
								&sgtcon_obj->rcondref);

	if( sgtcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument to libflame LAPACKE_sgtcon call is wrong\n", 
		         sgtcon_obj->info );
	}
	if( sgtcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument to Netlib LAPACKE_sgtcon call is wrong\n", 
		sgtcon_obj->inforef );
	}

#if LAPACKE_TEST_VERBOSE
    printf("sgtcon_obj->rcond: %lf   sgtcon_obj->rcondref: %lf \n", sgtcon_obj->rcond, sgtcon_obj->rcondref);
#endif
	EXPECT_NEAR(0.0, abs(sgtcon_obj->rcond - sgtcon_obj->rcondref), sgtcon_obj->threshold);
}

TEST(sgtcon_test, sgtcon1) {}
TEST(sgtcon_test, sgtcon2) {}
TEST(sgtcon_test, sgtcon3) {}
TEST(sgtcon_test, sgtcon4) {}

class gtcon_scomplex_parameters{
   public:
   /* input params to the API **/
   char norm;
   int 	n;
   lapack_complex_float  *dl, *dlref;
   lapack_complex_float  *d, *dref;
   lapack_complex_float  *du, *duref;
   lapack_complex_float  *du2, *du2ref;
   int     *ipiv, *ipivref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   float threshold;

   public: 
      gtcon_scomplex_parameters (char norm, int n );
      ~gtcon_scomplex_parameters ();

}; /* end of gtcon_scomplex_parameters  class definition */

/* Destructor definition **/
gtcon_scomplex_parameters:: ~gtcon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtcon_scomplex_parameters object: destructor invoked. \n");
#endif
   gtcon_free();
}

gtcon_scomplex_parameters::gtcon_scomplex_parameters (char norm_i, int n_i)
{
   int j;
   
   norm = norm_i;
   n = n_i;
   anorm = (float)n/2.0;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &d,  &dref,  n);
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &dl,  &dlref,  n-1);
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &du,  &duref,  n-1);
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &du2,  &du2ref,  n-2);
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (d==NULL) || (dref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (du==NULL) || (duref==NULL) ||  \
        (du2==NULL) || (du2ref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       gtcon_free();
       printf(" gtcon_scomplex_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( dl, dlref, n-1);
    lapacke_gtest_init_scomplex_buffer_pair_rand( du, duref, n-1);
    lapacke_gtest_init_scomplex_buffer_pair_rand( du2, du2ref, n-2);

    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }

	info = 0;
	inforef = 0;
}


/*  Test fixture class definition */
class cgtcon_test  : public  ::testing::Test {
public:
   gtcon_scomplex_parameters  *cgtcon_obj;
   void SetUp();
   void TearDown () { delete cgtcon_obj; }
};

void cgtcon_test::SetUp() {

   /* LAPACKE sgtCON prototype */
   typedef int (*Fptr_NL_LAPACKE_cgtcon) ( char norm, lapack_int n, const lapack_complex_float *dl,
                          const lapack_complex_float *d, const lapack_complex_float *du, const lapack_complex_float *du2, 
				          const lapack_int *ipiv, float anorm, float *rcond );
						   
   Fptr_NL_LAPACKE_cgtcon CGTCON;
   lapack_complex_float diff;
   void *hModule, *dModule;

   cgtcon_obj = new gtcon_scomplex_parameters (
					lin_solver_paramslist[idx].norm_gbcon,
                    lin_solver_paramslist[idx].n );

   cgtcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   CGTCON = (Fptr_NL_LAPACKE_cgtcon)dlsym(hModule, "LAPACKE_cgtcon");
   if (NULL == CGTCON)
   {
   	  printf("Could not get the symbol -LAPACKE_cgtcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   cgtcon_obj->info    = LAPACKE_cgttrf( cgtcon_obj->n, 
                                         cgtcon_obj->dl,
                                         cgtcon_obj->d, 
										 cgtcon_obj->du, 
										 cgtcon_obj->du2,
										cgtcon_obj->ipiv );

   cgtcon_obj->inforef = LAPACKE_cgttrf( cgtcon_obj->n, 
                                         cgtcon_obj->dlref,
                                         cgtcon_obj->dref, 
										 cgtcon_obj->duref, 
										 cgtcon_obj->du2ref,
										cgtcon_obj->ipivref );
   
   /* Compute libflame's Lapacke o/p  */
   cgtcon_obj->info    = LAPACKE_cgtcon( cgtcon_obj->norm, cgtcon_obj->n, 
                                        (const lapack_complex_float *)cgtcon_obj->dl,
                                        (const lapack_complex_float *)cgtcon_obj->d, 
										(const lapack_complex_float *)cgtcon_obj->du, 
										(const lapack_complex_float *)cgtcon_obj->du2,
										cgtcon_obj->ipiv, cgtcon_obj->anorm,
						                &cgtcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
    cgtcon_obj->inforef = CGTCON( cgtcon_obj->norm, cgtcon_obj->n, 
	                            (const lapack_complex_float *)cgtcon_obj->dlref, 
                                (const lapack_complex_float *)cgtcon_obj->dref, 
								(const lapack_complex_float *)cgtcon_obj->duref, 
								(const lapack_complex_float *)cgtcon_obj->du2ref, 
								cgtcon_obj->ipivref, cgtcon_obj->anorm,
								&cgtcon_obj->rcondref);

	if( cgtcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument to libflame LAPACKE_cgtcon call is wrong\n", 
		         cgtcon_obj->info );
	}
	if( cgtcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument to Netlib LAPACKE_cgtcon call is wrong\n", 
		cgtcon_obj->inforef );
	}

#if LAPACKE_TEST_VERBOSE
    printf("cgtcon_obj->rcond: %lf   cgtcon_obj->rcondref: %lf \n", cgtcon_obj->rcond, cgtcon_obj->rcondref);
#endif
	EXPECT_NEAR(0.0, abs(cgtcon_obj->rcond - cgtcon_obj->rcondref), cgtcon_obj->threshold);
}

TEST(cgtcon_test, cgtcon1) {}
TEST(cgtcon_test, cgtcon2) {}
TEST(cgtcon_test, cgtcon3) {}
TEST(cgtcon_test, cgtcon4) {}

class gtcon_dcomplex_parameters{
   public:
   /* input params to the API **/
   char norm;
   int 	n;
   lapack_complex_double  *dl, *dlref;
   lapack_complex_double  *d, *dref;
   lapack_complex_double  *du, *duref;
   lapack_complex_double  *du2, *du2ref;
   int     *ipiv, *ipivref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   float threshold;

   public: 
      gtcon_dcomplex_parameters (char norm, int n );
      ~gtcon_dcomplex_parameters ();

}; /* end of gtcon_dcomplex_parameters  class definition */

/* Destructor definition **/
gtcon_dcomplex_parameters:: ~gtcon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gtcon_dcomplex_parameters object: destructor invoked. \n");
#endif
   gtcon_free();
}

gtcon_dcomplex_parameters::gtcon_dcomplex_parameters (char norm_i, int n_i)
{
   int j;
   
   norm = norm_i;
   n = n_i;
   anorm = (double)n/2.0;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &d,  &dref,  n);
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &dl,  &dlref,  n-1);
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &du,  &duref,  n-1);
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &du2,  &du2ref,  n-2);
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (d==NULL) || (dref==NULL) ||  \
        (dl==NULL) || (dlref==NULL) ||  \
        (du==NULL) || (duref==NULL) ||  \
        (du2==NULL) || (du2ref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       gtcon_free();
       printf(" gtcon_dcomplex_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( dl, dlref, n-1);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( du, duref, n-1);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( du2, du2ref, n-2);

    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }

	info = 0;
	inforef = 0;
}


/*  Test fixture class definition */
class zgtcon_test  : public  ::testing::Test {
public:
   gtcon_dcomplex_parameters  *zgtcon_obj;
   void SetUp();
   void TearDown () { delete zgtcon_obj; }
};

void zgtcon_test::SetUp() {

   /* LAPACKE sgtCON prototype */
   typedef int (*Fptr_NL_LAPACKE_zgtcon) ( char norm, lapack_int n, const lapack_complex_double *dl,
                          const lapack_complex_double *d, const lapack_complex_double *du, const lapack_complex_double *du2, 
				          const lapack_int *ipiv, double anorm, double *rcond );
						   
   Fptr_NL_LAPACKE_zgtcon ZGTCON;
   lapack_complex_double diff;
   void *hModule, *dModule;

   zgtcon_obj = new gtcon_dcomplex_parameters (
					lin_solver_paramslist[idx].norm_gbcon,
                    lin_solver_paramslist[idx].n );

   zgtcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   ZGTCON = (Fptr_NL_LAPACKE_zgtcon)dlsym(hModule, "LAPACKE_zgtcon");
   if (NULL == ZGTCON)
   {
   	  printf("Could not get the symbol -LAPACKE_zgtcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   zgtcon_obj->info    = LAPACKE_zgttrf( zgtcon_obj->n, 
                                         zgtcon_obj->dl,
                                         zgtcon_obj->d, 
										 zgtcon_obj->du, 
										 zgtcon_obj->du2,
										zgtcon_obj->ipiv );

   zgtcon_obj->inforef = LAPACKE_zgttrf( zgtcon_obj->n, 
                                         zgtcon_obj->dlref,
                                         zgtcon_obj->dref, 
										 zgtcon_obj->duref, 
										 zgtcon_obj->du2ref,
										zgtcon_obj->ipivref );
   
   /* Compute libflame's Lapacke o/p  */
   zgtcon_obj->info    = LAPACKE_zgtcon( zgtcon_obj->norm, zgtcon_obj->n, 
                                        (const lapack_complex_double *)zgtcon_obj->dl,
                                        (const lapack_complex_double *)zgtcon_obj->d, 
										(const lapack_complex_double *)zgtcon_obj->du, 
										(const lapack_complex_double *)zgtcon_obj->du2,
										zgtcon_obj->ipiv, zgtcon_obj->anorm,
						                &zgtcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
    zgtcon_obj->inforef = ZGTCON( zgtcon_obj->norm, zgtcon_obj->n, 
	                            (const lapack_complex_double *)zgtcon_obj->dlref, 
                                (const lapack_complex_double *)zgtcon_obj->dref, 
								(const lapack_complex_double *)zgtcon_obj->duref, 
								(const lapack_complex_double *)zgtcon_obj->du2ref, 
								zgtcon_obj->ipivref, zgtcon_obj->anorm,
								&zgtcon_obj->rcondref);

	if( zgtcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument to libflame LAPACKE_zgtcon call is wrong\n", 
		         zgtcon_obj->info );
	}
	if( zgtcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument to Netlib LAPACKE_zgtcon call is wrong\n", 
		zgtcon_obj->inforef );
	}

#if LAPACKE_TEST_VERBOSE
    printf("zgtcon_obj->rcond: %lf   zgtcon_obj->rcondref: %lf \n", zgtcon_obj->rcond, zgtcon_obj->rcondref);
#endif
	EXPECT_NEAR(0.0, abs(zgtcon_obj->rcond - zgtcon_obj->rcondref), zgtcon_obj->threshold);
}

TEST(zgtcon_test, zgtcon1) {}
TEST(zgtcon_test, zgtcon2) {}
TEST(zgtcon_test, zgtcon3) {}
TEST(zgtcon_test, zgtcon4) {}
