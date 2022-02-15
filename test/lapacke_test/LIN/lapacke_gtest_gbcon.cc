#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#include <math.h>
#include "lapacke.h"

#define gbcon_free() \
       free (ab   ); \
       free (abref); \
       free (ipiv   ); \
       free (ipivref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

class gbcon_double_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;   
   char norm;
   int 	n;
   int 	kl, ku;
   double  *ab, *abref;
   int 	ldab;
   int     *ipiv, *ipivref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      gbcon_double_parameters (int matrix_layout, char norm, int n,
                            int kl, int ku, int ldab );
      ~gbcon_double_parameters ();

}; /* end of gbcon_double_parameters  class definition */

/* Destructor definition **/
gbcon_double_parameters:: ~gbcon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbcon_double_parameters object: destructor invoked. \n");
#endif
   gbcon_free();
}

gbcon_double_parameters::gbcon_double_parameters (int matrix_layout_i, 
                             char norm_i, int n_i, int kl_i, int ku_i,
                             int ldab_i){
   int j;
   
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   n = n_i;
   kl = kl_i;
   ku = ku_i;
   ldab = ldab_i;
   anorm = (float)n/2.0;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &ab,  &abref,  (ldab*n));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (ab==NULL) || (abref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       gbcon_free();
       printf(" gbcon_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ab, abref, (ldab*n));
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class dgbcon_test  : public  ::testing::Test {
public:
   gbcon_double_parameters  *dgbcon_obj;
   void SetUp();
   void TearDown () { delete dgbcon_obj; }
};

void dgbcon_test::SetUp() {

   /* LAPACKE DGBCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dgbcon) ( int matrix_layout, char norm, lapack_int n,
                                          lapack_int kl, lapack_int ku, const double *ab, 
										  lapack_int ldab, const lapack_int *ipiv, 
										  double anorm, double *rcond );
						   
   Fptr_NL_LAPACKE_dgbcon DGBCON;
   double diff;
   void *hModule, *dModule;

   dgbcon_obj = new gbcon_double_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl_gbcon,
						   lin_solver_paramslist[idx].ku_gbcon,
						   lin_solver_paramslist[idx].ldab_gbcon );

   dgbcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
 
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   DGBCON = (Fptr_NL_LAPACKE_dgbcon)dlsym(hModule, "LAPACKE_dgbcon");
   if (NULL == DGBCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dgbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   /* Invoke the gbtrf API to get the factored contents & pivot buffer computation */
   dgbcon_obj->info    = LAPACKE_dgbtrf( dgbcon_obj->matrix_layout, dgbcon_obj->n, 
                                        dgbcon_obj->n, dgbcon_obj->kl, dgbcon_obj->ku,
                                        dgbcon_obj->ab, dgbcon_obj->ldab, 
										dgbcon_obj->ipiv);

   dgbcon_obj->inforef = LAPACKE_dgbtrf( dgbcon_obj->matrix_layout, dgbcon_obj->n, 
                                        dgbcon_obj->n, dgbcon_obj->kl, dgbcon_obj->ku,
                                        dgbcon_obj->abref, dgbcon_obj->ldab, 
										dgbcon_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */
   dgbcon_obj->info    = LAPACKE_dgbcon( dgbcon_obj->matrix_layout, dgbcon_obj->norm, 
                                        dgbcon_obj->n, dgbcon_obj->kl, dgbcon_obj->ku,
                                        dgbcon_obj->ab, dgbcon_obj->ldab, 
										dgbcon_obj->ipiv, dgbcon_obj->anorm,
						                &dgbcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   dgbcon_obj->inforef = DGBCON( dgbcon_obj->matrix_layout, dgbcon_obj->norm, 
                                dgbcon_obj->n, dgbcon_obj->kl, dgbcon_obj->ku, 
								dgbcon_obj->abref, dgbcon_obj->ldab, 
								dgbcon_obj->ipivref, dgbcon_obj->anorm,
								&dgbcon_obj->rcondref);

	if( dgbcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dgbcon is wrong\n", 
		         dgbcon_obj->info );
	}
	if( dgbcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dgbcon is wrong\n", 
		dgbcon_obj->inforef );
	}

	/* Compute Difference in libflame and the  reference o/p buffer */	
	diff =  computeDiff_d( dgbcon_obj->ldab*dgbcon_obj->n, dgbcon_obj->ab, dgbcon_obj->abref );
	if ( diff > dgbcon_obj->threshold)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
    printf("dgbcon_obj->rcond: %lf   dgbcon_obj->rcondref: %lf \n", dgbcon_obj->rcond, dgbcon_obj->rcondref);
	EXPECT_NEAR(0.0, abs(dgbcon_obj->rcond - dgbcon_obj->rcondref), dgbcon_obj->threshold);
}

TEST(dgbcon_test, dgbcon1) {}
TEST(dgbcon_test, dgbcon2) {}
TEST(dgbcon_test, dgbcon3) {}
TEST(dgbcon_test, dgbcon4) {}

class gbcon_float_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;   
   char norm;
   int 	n;
   int 	kl, ku;
   float  *ab, *abref;
   int 	ldab;
   int     *ipiv, *ipivref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      gbcon_float_parameters (int matrix_layout, char norm, int n,
                            int kl, int ku, int ldab );
      ~gbcon_float_parameters ();

}; /* end of gbcon_float_parameters  class definition */

/* Destructor definition **/
gbcon_float_parameters:: ~gbcon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbcon_float_parameters object: destructor invoked. \n");
#endif
   gbcon_free();
}

gbcon_float_parameters::gbcon_float_parameters (int matrix_layout_i, 
                             char norm_i, int n_i, int kl_i, int ku_i,
                             int ldab_i){
   int j;
   
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   n = n_i;
   kl = kl_i;
   ku = ku_i;
   ldab = ldab_i;
   anorm = (float)n/2.0;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &ab,  &abref,  (ldab*n));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (ab==NULL) || (abref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       gbcon_free();
       printf(" gbcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ab, abref, (ldab*n));
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class sgbcon_test  : public  ::testing::Test {
public:
   gbcon_float_parameters  *sgbcon_obj;
   void SetUp();
   void TearDown () { delete sgbcon_obj; }
};

void sgbcon_test::SetUp() {

   /* LAPACKE sgbcon prototype */
   typedef int (*Fptr_NL_LAPACKE_sgbcon) ( int matrix_layout, char norm, lapack_int n,
                                          lapack_int kl, lapack_int ku, const float *ab, 
										  lapack_int ldab, const lapack_int *ipiv, 
										  float anorm, float *rcond );
						   
   Fptr_NL_LAPACKE_sgbcon sgbcon;
   float diff;
   void *hModule, *dModule;

   sgbcon_obj = new gbcon_float_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl_gbcon,
						   lin_solver_paramslist[idx].ku_gbcon,
						   lin_solver_paramslist[idx].ldab_gbcon );

   sgbcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
 
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   sgbcon = (Fptr_NL_LAPACKE_sgbcon)dlsym(hModule, "LAPACKE_sgbcon");
   if (NULL == sgbcon)
   {
   	  printf("Could not get the symbol -LAPACKE_sgbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   /* Invoke the gbtrf API to get the factored contents & pivot buffer computation */
   sgbcon_obj->info    = LAPACKE_sgbtrf( sgbcon_obj->matrix_layout, sgbcon_obj->n, 
                                        sgbcon_obj->n, sgbcon_obj->kl, sgbcon_obj->ku,
                                        sgbcon_obj->ab, sgbcon_obj->ldab, 
										sgbcon_obj->ipiv);

   sgbcon_obj->inforef = LAPACKE_sgbtrf( sgbcon_obj->matrix_layout, sgbcon_obj->n, 
                                        sgbcon_obj->n, sgbcon_obj->kl, sgbcon_obj->ku,
                                        sgbcon_obj->abref, sgbcon_obj->ldab, 
										sgbcon_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */
   sgbcon_obj->info    = LAPACKE_sgbcon( sgbcon_obj->matrix_layout, sgbcon_obj->norm, 
                                        sgbcon_obj->n, sgbcon_obj->kl, sgbcon_obj->ku,
                                        sgbcon_obj->ab, sgbcon_obj->ldab, 
										sgbcon_obj->ipiv, sgbcon_obj->anorm,
						                &sgbcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   sgbcon_obj->inforef = sgbcon( sgbcon_obj->matrix_layout, sgbcon_obj->norm, 
                                sgbcon_obj->n, sgbcon_obj->kl, sgbcon_obj->ku, 
								sgbcon_obj->abref, sgbcon_obj->ldab, 
								sgbcon_obj->ipivref, sgbcon_obj->anorm,
								&sgbcon_obj->rcondref);

	if( sgbcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_sgbcon is wrong\n", 
		         sgbcon_obj->info );
	}
	if( sgbcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_sgbcon is wrong\n", 
		sgbcon_obj->inforef );
	}

	/* Compute Difference in libflame and the  reference o/p buffer */	
	diff =  computeDiff_s( sgbcon_obj->ldab*sgbcon_obj->n, sgbcon_obj->ab, sgbcon_obj->abref );
	if ( diff > sgbcon_obj->threshold)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
    printf("sgbcon_obj->rcond: %lf   sgbcon_obj->rcondref: %lf \n", sgbcon_obj->rcond, sgbcon_obj->rcondref);
	EXPECT_NEAR(0.0, abs(sgbcon_obj->rcond - sgbcon_obj->rcondref), sgbcon_obj->threshold);
}

TEST(sgbcon_test, sgbcon1) {}
TEST(sgbcon_test, sgbcon2) {}
TEST(sgbcon_test, sgbcon3) {}
TEST(sgbcon_test, sgbcon4) {}


class gbcon_scomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;   
   char norm;
   int 	n;
   int 	kl, ku;
   lapack_complex_float  *ab, *abref;
   int 	ldab;
   int     *ipiv, *ipivref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      gbcon_scomplex_parameters (int matrix_layout, char norm, int n,
                            int kl, int ku, int ldab );
      ~gbcon_scomplex_parameters ();

}; /* end of gbcon_scomplex_parameters  class definition */

/* Destructor definition **/
gbcon_scomplex_parameters:: ~gbcon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbcon_scomplex_parameters object: destructor invoked. \n");
#endif
   gbcon_free();
}

gbcon_scomplex_parameters::gbcon_scomplex_parameters (int matrix_layout_i, 
                             char norm_i, int n_i, int kl_i, int ku_i,
                             int ldab_i){
   int j;
   
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   n = n_i;
   kl = kl_i;
   ku = ku_i;
   ldab = ldab_i;
   anorm = (float)n/2.0;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ab,  &abref,  (ldab*n));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (ab==NULL) || (abref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       gbcon_free();
       printf(" gbcon_scomplex_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ab, abref, (ldab*n));
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class cgbcon_test  : public  ::testing::Test {
public:
   gbcon_scomplex_parameters  *cgbcon_obj;
   void SetUp();
   void TearDown () { delete cgbcon_obj; }
};

void cgbcon_test::SetUp() {

   /* LAPACKE cgbcon prototype */
   typedef int (*Fptr_NL_LAPACKE_cgbcon) ( int matrix_layout, char norm, lapack_int n,
                                          lapack_int kl, lapack_int ku, const lapack_complex_float *ab, 
										  lapack_int ldab, const lapack_int *ipiv, 
										  float anorm, float *rcond );
						   
   Fptr_NL_LAPACKE_cgbcon cgbcon;
   float diff;
   void *hModule, *dModule;

   cgbcon_obj = new gbcon_scomplex_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl_gbcon,
						   lin_solver_paramslist[idx].ku_gbcon,
						   lin_solver_paramslist[idx].ldab_gbcon );

   cgbcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
 
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   cgbcon = (Fptr_NL_LAPACKE_cgbcon)dlsym(hModule, "LAPACKE_cgbcon");
   if (NULL == cgbcon)
   {
   	  printf("Could not get the symbol -LAPACKE_cgbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   /* Invoke the gbtrf API to get the factored contents & pivot buffer computation */
   cgbcon_obj->info    = LAPACKE_cgbtrf( cgbcon_obj->matrix_layout, cgbcon_obj->n, 
                                        cgbcon_obj->n, cgbcon_obj->kl, cgbcon_obj->ku,
                                        cgbcon_obj->ab, cgbcon_obj->ldab, 
										cgbcon_obj->ipiv);

   cgbcon_obj->inforef = LAPACKE_cgbtrf( cgbcon_obj->matrix_layout, cgbcon_obj->n, 
                                        cgbcon_obj->n, cgbcon_obj->kl, cgbcon_obj->ku,
                                        cgbcon_obj->abref, cgbcon_obj->ldab, 
										cgbcon_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */
   cgbcon_obj->info    = LAPACKE_cgbcon( cgbcon_obj->matrix_layout, cgbcon_obj->norm, 
                                        cgbcon_obj->n, cgbcon_obj->kl, cgbcon_obj->ku,
                                        cgbcon_obj->ab, cgbcon_obj->ldab, 
										cgbcon_obj->ipiv, cgbcon_obj->anorm,
						                &cgbcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   cgbcon_obj->inforef = cgbcon( cgbcon_obj->matrix_layout, cgbcon_obj->norm, 
                                cgbcon_obj->n, cgbcon_obj->kl, cgbcon_obj->ku, 
								cgbcon_obj->abref, cgbcon_obj->ldab, 
								cgbcon_obj->ipivref, cgbcon_obj->anorm,
								&cgbcon_obj->rcondref);

	if( cgbcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_cgbcon is wrong\n", 
		         cgbcon_obj->info );
	}
	if( cgbcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_cgbcon is wrong\n", 
		cgbcon_obj->inforef );
	}

	/* Compute Difference in libflame and the  reference o/p buffer */	
	diff =  computeDiff_c( cgbcon_obj->ldab*cgbcon_obj->n, cgbcon_obj->ab, cgbcon_obj->abref );
	if ( diff > cgbcon_obj->threshold)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
    printf("cgbcon_obj->rcond: %lf   cgbcon_obj->rcondref: %lf \n", cgbcon_obj->rcond, cgbcon_obj->rcondref);
	EXPECT_NEAR(0.0, abs(cgbcon_obj->rcond - cgbcon_obj->rcondref), cgbcon_obj->threshold);
	EXPECT_NEAR(0.0, diff, cgbcon_obj->threshold);
}

TEST(cgbcon_test, cgbcon1) {}
TEST(cgbcon_test, cgbcon2) {}
TEST(cgbcon_test, cgbcon3) {}
TEST(cgbcon_test, cgbcon4) {}

class gbcon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;   
   char norm;
   int 	n;
   int 	kl, ku;
   lapack_complex_double  *ab, *abref;
   int 	ldab;
   int     *ipiv, *ipivref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      gbcon_dcomplex_parameters (int matrix_layout, char norm, int n,
                            int kl, int ku, int ldab );
      ~gbcon_dcomplex_parameters ();

}; /* end of gbcon_dcomplex_parameters  class definition */

/* Destructor definition **/
gbcon_dcomplex_parameters:: ~gbcon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gbcon_dcomplex_parameters object: destructor invoked. \n");
#endif
   gbcon_free();
}

gbcon_dcomplex_parameters::gbcon_dcomplex_parameters (int matrix_layout_i, 
                             char norm_i, int n_i, int kl_i, int ku_i,
                             int ldab_i){
   int j;
   
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   n = n_i;
   kl = kl_i;
   ku = ku_i;
   ldab = ldab_i;
   anorm = (double)n/2.0;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ab,  &abref,  (ldab*n));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (ab==NULL) || (abref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       gbcon_free();
       printf(" gbcon_dcomplex_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ab, abref, (ldab*n));
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class zgbcon_test  : public  ::testing::Test {
public:
   gbcon_dcomplex_parameters  *zgbcon_obj;
   void SetUp();
   void TearDown () { delete zgbcon_obj; }
};

void zgbcon_test::SetUp() {

   /* LAPACKE zgbcon prototype */
   typedef int (*Fptr_NL_LAPACKE_zgbcon) ( int matrix_layout, char norm, lapack_int n,
                                          lapack_int kl, lapack_int ku, const lapack_complex_double *ab, 
										  lapack_int ldab, const lapack_int *ipiv, 
										  double anorm, double *rcond );
						   
   Fptr_NL_LAPACKE_zgbcon zgbcon;
   double diff;
   void *hModule, *dModule;

   zgbcon_obj = new gbcon_dcomplex_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl_gbcon,
						   lin_solver_paramslist[idx].ku_gbcon,
						   lin_solver_paramslist[idx].ldab_gbcon );

   zgbcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
 
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   zgbcon = (Fptr_NL_LAPACKE_zgbcon)dlsym(hModule, "LAPACKE_zgbcon");
   if (NULL == zgbcon)
   {
   	  printf("Could not get the symbol -LAPACKE_zgbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   /* Invoke the gbtrf API to get the factored contents & pivot buffer computation */
   zgbcon_obj->info    = LAPACKE_zgbtrf( zgbcon_obj->matrix_layout, zgbcon_obj->n, 
                                        zgbcon_obj->n, zgbcon_obj->kl, zgbcon_obj->ku,
                                        zgbcon_obj->ab, zgbcon_obj->ldab, 
										zgbcon_obj->ipiv);

   zgbcon_obj->inforef = LAPACKE_zgbtrf( zgbcon_obj->matrix_layout, zgbcon_obj->n, 
                                        zgbcon_obj->n, zgbcon_obj->kl, zgbcon_obj->ku,
                                        zgbcon_obj->abref, zgbcon_obj->ldab, 
										zgbcon_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */
   zgbcon_obj->info    = LAPACKE_zgbcon( zgbcon_obj->matrix_layout, zgbcon_obj->norm, 
                                        zgbcon_obj->n, zgbcon_obj->kl, zgbcon_obj->ku,
                                        zgbcon_obj->ab, zgbcon_obj->ldab, 
										zgbcon_obj->ipiv, zgbcon_obj->anorm,
						                &zgbcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   zgbcon_obj->inforef = zgbcon( zgbcon_obj->matrix_layout, zgbcon_obj->norm, 
                                zgbcon_obj->n, zgbcon_obj->kl, zgbcon_obj->ku, 
								zgbcon_obj->abref, zgbcon_obj->ldab, 
								zgbcon_obj->ipivref, zgbcon_obj->anorm,
								&zgbcon_obj->rcondref);

	if( zgbcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_zgbcon is wrong\n", 
		         zgbcon_obj->info );
	}
	if( zgbcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_zgbcon is wrong\n", 
		zgbcon_obj->inforef );
	}

	/* Compute Difference in libflame and the  reference o/p buffer */	
	diff =  computeDiff_z( zgbcon_obj->ldab*zgbcon_obj->n, zgbcon_obj->ab, zgbcon_obj->abref );
	if ( diff > zgbcon_obj->threshold)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
    printf("zgbcon_obj->rcond: %lf   zgbcon_obj->rcondref: %lf \n", zgbcon_obj->rcond, zgbcon_obj->rcondref);
	EXPECT_NEAR(0.0, abs(zgbcon_obj->rcond - zgbcon_obj->rcondref), zgbcon_obj->threshold);
	EXPECT_NEAR(0.0, diff, zgbcon_obj->threshold);
}

TEST(zgbcon_test, zgbcon1) {}
TEST(zgbcon_test, zgbcon2) {}
TEST(zgbcon_test, zgbcon3) {}
TEST(zgbcon_test, zgbcon4) {}
