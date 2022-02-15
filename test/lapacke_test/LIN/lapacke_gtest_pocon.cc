#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#include <math.h>
#include "lapacke.h"

#define pocon_free() \
       free (a   ); \
       free (aref);

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

class pocon_double_parameters{
   public:
   float threshold;
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   double  *a, *aref;
   int 	lda;
   double anorm;
   double rcond, rcondref;   int info, inforef;
   public: 
      pocon_double_parameters (int matrix_layout, char uplo, int n);
      ~pocon_double_parameters ();

}; /* end of pocon_double_parameters  class definition */

/* Destructor definition **/
pocon_double_parameters:: ~pocon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pocon_double_parameters object: destructor invoked. \n");
#endif
   pocon_free();
}

pocon_double_parameters::pocon_double_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &a,  &aref,  (lda*n));
	
	if( (a==NULL) || (aref==NULL) ){
       pocon_free();
       printf(" pocon_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, (lda*n));
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class dpocon_test  : public  ::testing::Test {
public:
   pocon_double_parameters  *dpocon_obj;
   void SetUp();
   void TearDown () { delete dpocon_obj; }
};

void dpocon_test::SetUp(){

   /* LAPACKE DPOCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dpocon) ( int matrix_layout, char uplo, 
                          lapack_int n, const double* a, lapack_int lda, 
						  double anorm, double* rcond  );
						  
   Fptr_NL_LAPACKE_dpocon DPOCON;
   double diff;
   void *hModule, *dModule;

   dpocon_obj = new pocon_double_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   dpocon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   DPOCON = (Fptr_NL_LAPACKE_dpocon)dlsym(hModule, "LAPACKE_dpocon");
   if (NULL == DPOCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dpocon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   dpocon_obj->info    = LAPACKE_dpocon( dpocon_obj->matrix_layout, dpocon_obj->uplo, 
                                        dpocon_obj->n, (const double*)dpocon_obj->a, dpocon_obj->lda, 
										dpocon_obj->anorm, &dpocon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   dpocon_obj->inforef = DPOCON( dpocon_obj->matrix_layout, dpocon_obj->uplo, 
                                dpocon_obj->n, (const double*)dpocon_obj->aref, dpocon_obj->lda, 
								dpocon_obj->anorm, &dpocon_obj->rcondref);
 	if( dpocon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dpocon is wrong\n", 
		         dpocon_obj->info );
	}
	if( dpocon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dpocon is wrong\n", 
		dpocon_obj->inforef );
	}

	/* Compute Difference between Netlib ref and libflame o/ps */
	
	diff =  computeDiff_d( dpocon_obj->lda*dpocon_obj->n, dpocon_obj->a, dpocon_obj->aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
#if LAPACKE_TEST_VERBOSE
    printf("dpocon_obj->rcond: %lf   dpocon_obj->rcondref: %lf \n", dpocon_obj->rcond, dpocon_obj->rcondref);
	printf("diff: %lf  \n", diff);
#endif
	EXPECT_NEAR(0.0, abs(dpocon_obj->rcond - dpocon_obj->rcondref), dpocon_obj->threshold);
}

TEST(pocon, dpocon1) {}
TEST(pocon, dpocon2) {}
TEST(pocon, dpocon3) {}
TEST(pocon, dpocon4) {}

class pocon_float_parameters{
   public:
   float threshold;
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   float  *a, *aref;
   int 	lda;
   float anorm;
   float rcond, rcondref;   int info, inforef;
   public: 
      pocon_float_parameters (int matrix_layout, char uplo, int n);
      ~pocon_float_parameters ();

}; /* end of pocon_float_parameters  class definition */

/* Destructor definition **/
pocon_float_parameters:: ~pocon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pocon_float_parameters object: destructor invoked. \n");
#endif
   pocon_free();
}

pocon_float_parameters::pocon_float_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n_i;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &a,  &aref,  (lda*n));
	
	if( (a==NULL) || (aref==NULL) ){
       pocon_free();
       printf(" pocon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, (lda*n));
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class spocon_test  : public  ::testing::Test {
public:
   pocon_float_parameters  *spocon_obj;
   void SetUp();
   void TearDown () { delete spocon_obj; }
};

void spocon_test::SetUp(){

   /* LAPACKE SPOCON prototype */
   typedef int (*Fptr_NL_LAPACKE_spocon) ( int matrix_layout, char uplo, 
                          lapack_int n, const float* a, lapack_int lda, 
						  float anorm, float* rcond  );
						  
   Fptr_NL_LAPACKE_spocon SPOCON;
   float diff;
   void *hModule, *dModule;

   spocon_obj = new pocon_float_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   spocon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   SPOCON = (Fptr_NL_LAPACKE_spocon)dlsym(hModule, "LAPACKE_spocon");
   if (NULL == SPOCON)
   {
   	  printf("Could not get the symbol -LAPACKE_spocon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   spocon_obj->info    = LAPACKE_spocon( spocon_obj->matrix_layout, spocon_obj->uplo, 
                                        spocon_obj->n, (const float*)spocon_obj->a, spocon_obj->lda, 
										spocon_obj->anorm, &spocon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   spocon_obj->inforef = SPOCON( spocon_obj->matrix_layout, spocon_obj->uplo, 
                                spocon_obj->n, (const float*)spocon_obj->aref, spocon_obj->lda, 
								spocon_obj->anorm, &spocon_obj->rcondref);
 	if( spocon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_spocon is wrong\n", 
		         spocon_obj->info );
	}
	if( spocon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_spocon is wrong\n", 
		spocon_obj->inforef );
	}

	/* Compute Difference between Netlib ref and libflame o/ps */
	
	diff =  computeDiff_s( spocon_obj->lda*spocon_obj->n, spocon_obj->a, spocon_obj->aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
#if LAPACKE_TEST_VERBOSE
    printf("spocon_obj->rcond: %lf   spocon_obj->rcondref: %lf \n", spocon_obj->rcond, spocon_obj->rcondref);
	printf("diff: %lf  \n", diff);
#endif
	EXPECT_NEAR(0.0, abs(spocon_obj->rcond - spocon_obj->rcondref), spocon_obj->threshold);
	EXPECT_NEAR(0.0, diff, spocon_obj->threshold);
}

TEST(pocon, spocon1) {}
TEST(pocon, spocon2) {}
TEST(pocon, spocon3) {}
TEST(pocon, spocon4) {}

class pocon_scomplex_parameters{
   public:
   float threshold;
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   lapack_complex_float  *a, *aref;
   int 	lda;
   float anorm;
   float rcond, rcondref;   int info, inforef;
   public: 
      pocon_scomplex_parameters (int matrix_layout, char uplo, int n);
      ~pocon_scomplex_parameters ();

}; /* end of pocon_scomplex_parameters  class definition */

/* Destructor definition **/
pocon_scomplex_parameters:: ~pocon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pocon_scomplex_parameters object: destructor invoked. \n");
#endif
   pocon_free();
}

pocon_scomplex_parameters::pocon_scomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n_i;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a,  &aref,  (lda*n));
	
	if( (a==NULL) || (aref==NULL) ){
       pocon_free();
       printf(" pocon_scomplex_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (lda*n));
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class cpocon_test  : public  ::testing::Test {
public:
   pocon_scomplex_parameters  *cpocon_obj;
   void SetUp();
   void TearDown () { delete cpocon_obj; }
};

void cpocon_test::SetUp(){

   /* LAPACKE CPOCON prototype */
   typedef int (*Fptr_NL_LAPACKE_cpocon) ( int matrix_layout, char uplo, lapack_int n,
                      const lapack_complex_float *a, lapack_int lda,
					  float anorm, float *rcond );
						  
   Fptr_NL_LAPACKE_cpocon CPOCON;
   float diff;
   void *hModule, *dModule;

   cpocon_obj = new pocon_scomplex_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   cpocon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   CPOCON = (Fptr_NL_LAPACKE_cpocon)dlsym(hModule, "LAPACKE_cpocon");
   if (NULL == CPOCON)
   {
   	  printf("Could not get the symbol -LAPACKE_cpocon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   cpocon_obj->info    = LAPACKE_cpocon( cpocon_obj->matrix_layout, cpocon_obj->uplo, 
                                        cpocon_obj->n, (const lapack_complex_float*)cpocon_obj->a,
										cpocon_obj->lda, cpocon_obj->anorm, &cpocon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   cpocon_obj->inforef = CPOCON( cpocon_obj->matrix_layout, cpocon_obj->uplo, 
                                cpocon_obj->n, (const lapack_complex_float*)cpocon_obj->aref,
								cpocon_obj->lda, cpocon_obj->anorm, &cpocon_obj->rcondref);

 	if( cpocon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_cpocon is wrong\n", 
		         cpocon_obj->info );
	}
	if( cpocon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_cpocon is wrong\n", 
		cpocon_obj->inforef );
	}

	/* Compute Difference between Netlib ref and libflame o/ps */
	
	diff =  computeDiff_c( cpocon_obj->lda*cpocon_obj->n, cpocon_obj->a, cpocon_obj->aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
#if LAPACKE_TEST_VERBOSE
    printf("cpocon_obj->rcond: %lf   cpocon_obj->rcondref: %lf \n", cpocon_obj->rcond, cpocon_obj->rcondref);
	printf("diff: %lf  \n", diff);
#endif
	EXPECT_NEAR(0.0, abs(cpocon_obj->rcond - cpocon_obj->rcondref), cpocon_obj->threshold);
	EXPECT_NEAR(0.0, diff, cpocon_obj->threshold);
}

TEST(pocon, cpocon1) {}
TEST(pocon, cpocon2) {}
TEST(pocon, cpocon3) {}
TEST(pocon, cpocon4) {}

class pocon_dcomplex_parameters{
   public:
   double threshold;
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   lapack_complex_double  *a, *aref;
   int 	lda;
   double anorm;
   double rcond, rcondref;   int info, inforef;
   public: 
      pocon_dcomplex_parameters (int matrix_layout, char uplo, int n);
      ~pocon_dcomplex_parameters ();

}; /* end of pocon_dcomplex_parameters  class definition */

/* Destructor definition **/
pocon_dcomplex_parameters:: ~pocon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pocon_dcomplex_parameters object: destructor invoked. \n");
#endif
   pocon_free();
}

pocon_dcomplex_parameters::pocon_dcomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a,  &aref,  (lda*n));
	
	if( (a==NULL) || (aref==NULL) ){
       pocon_free();
       printf(" pocon_dcomplex_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, (lda*n));
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class zpocon_test  : public  ::testing::Test {
public:
   pocon_dcomplex_parameters  *zpocon_obj;
   void SetUp();
   void TearDown () { delete zpocon_obj; }
};

void zpocon_test::SetUp(){

   /* LAPACKE CPOCON prototype */
   typedef int (*Fptr_NL_LAPACKE_zpocon) ( int matrix_layout, char uplo, lapack_int n,
                      const lapack_complex_double *a, lapack_int lda,
					  double anorm, double *rcond );
						  
   Fptr_NL_LAPACKE_zpocon CPOCON;
   double diff;
   void *hModule, *dModule;

   zpocon_obj = new pocon_dcomplex_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   zpocon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   CPOCON = (Fptr_NL_LAPACKE_zpocon)dlsym(hModule, "LAPACKE_zpocon");
   if (NULL == CPOCON)
   {
   	  printf("Could not get the symbol -LAPACKE_zpocon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   zpocon_obj->info    = LAPACKE_zpocon( zpocon_obj->matrix_layout, zpocon_obj->uplo, 
                                        zpocon_obj->n, (const lapack_complex_double*)zpocon_obj->a,
										zpocon_obj->lda, zpocon_obj->anorm, &zpocon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   zpocon_obj->inforef = CPOCON( zpocon_obj->matrix_layout, zpocon_obj->uplo, 
                                zpocon_obj->n, (const lapack_complex_double*)zpocon_obj->aref,
								zpocon_obj->lda, zpocon_obj->anorm, &zpocon_obj->rcondref);

 	if( zpocon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_zpocon is wrong\n", 
		         zpocon_obj->info );
	}
	if( zpocon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_zpocon is wrong\n", 
		zpocon_obj->inforef );
	}

	/* Compute Difference between Netlib ref and libflame o/ps */
	
	diff =  computeDiff_z( zpocon_obj->lda*zpocon_obj->n, zpocon_obj->a, zpocon_obj->aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
#if LAPACKE_TEST_VERBOSE
    printf("zpocon_obj->rcond: %lf   zpocon_obj->rcondref: %lf \n", zpocon_obj->rcond, zpocon_obj->rcondref);
	printf("diff: %lf  \n", diff);
#endif
	EXPECT_NEAR(0.0, abs(zpocon_obj->rcond - zpocon_obj->rcondref), zpocon_obj->threshold);
	EXPECT_NEAR(0.0, diff, zpocon_obj->threshold);
}

TEST(pocon, zpocon1) {}
TEST(pocon, zpocon2) {}
TEST(pocon, zpocon3) {}
TEST(pocon, zpocon4) {}
