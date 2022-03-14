#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#include <math.h>
#include "lapacke.h"

#define hecon_3_free() \
       free (a   ); \
       free (aref); \
       free (e   ); \
       free (eref); \
       free (ipiv   ); \
       free (ipivref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin hecon_3_scomplex_parameters  class definition */
class hecon_3_scomplex_parameters{
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
      hecon_3_scomplex_parameters (int matrix_layout, char uplo, int n);
      ~hecon_3_scomplex_parameters ();

}; /* end of hecon_3_float_parameters  class definition */

/* Destructor definition **/
hecon_3_scomplex_parameters:: ~hecon_3_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hecon_3_float_parameters object: destructor invoked. \n");
#endif
   hecon_3_free();
}

hecon_3_scomplex_parameters::hecon_3_scomplex_parameters (int matrix_layout_i, 
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
       hecon_3_free();
       printf(" hecon_3_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    //lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (lda*n));
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(a, aref,
                     n, lda, uplo);
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
class checon_3_test  : public  ::testing::Test {
public:
   hecon_3_scomplex_parameters  *checon_3_obj;
   void SetUp();
   void TearDown () { delete checon_3_obj; }
};

void checon_3_test::SetUp(){

   /* LAPACKE CHECON_3 prototype */
   typedef int (*Fptr_NL_LAPACKE_checon_3) ( int matrix_layout, char uplo, 
             lapack_int n, const lapack_complex_float* a, lapack_int lda,
             const lapack_complex_float* e,
			 const lapack_int* ipiv, float anorm, float* rcond );
						   
   Fptr_NL_LAPACKE_checon_3 CHECON_3;
   float diff;
   void *hModule, *dModule;

	
   checon_3_obj = new hecon_3_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   checon_3_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
 
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   CHECON_3 = (Fptr_NL_LAPACKE_checon_3)dlsym(hModule, "LAPACKE_checon_3");
   if (NULL == CHECON_3)
   {
   	  printf("Could not get the symbol -LAPACKE_checon_3- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   checon_3_obj->info =  LAPACKE_chetrf_rk (checon_3_obj->matrix_layout, checon_3_obj->uplo, 
                                 checon_3_obj->n,
                                 checon_3_obj->a, 
								 checon_3_obj->lda, 
                                 checon_3_obj->e, 
					             checon_3_obj->ipiv );

   checon_3_obj->inforef =  LAPACKE_chetrf_rk (checon_3_obj->matrix_layout, checon_3_obj->uplo, 
                                 checon_3_obj->n,
                                 checon_3_obj->aref, 
								 checon_3_obj->lda, 
                                 checon_3_obj->eref, 
					             checon_3_obj->ipivref );

   /* Compute libflame's Lapacke o/p  */  
   checon_3_obj->info    = LAPACKE_checon_3( checon_3_obj->matrix_layout, checon_3_obj->uplo, 
                                 checon_3_obj->n,
                                 (const lapack_complex_float *)checon_3_obj->a, 
								 checon_3_obj->lda, 
                                 (const lapack_complex_float *)checon_3_obj->e, 
					             (const lapack_int *)checon_3_obj->ipiv, 
								 checon_3_obj->anorm, &checon_3_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   checon_3_obj->inforef = CHECON_3( checon_3_obj->matrix_layout, checon_3_obj->uplo, 
                                checon_3_obj->n, 
				                (const lapack_complex_float *)checon_3_obj->aref, 
								checon_3_obj->lda,
								(const lapack_complex_float* )checon_3_obj->eref, 
			                    (const lapack_int *)checon_3_obj->ipivref, 
								checon_3_obj->anorm, &checon_3_obj->rcondref);

	if( checon_3_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_checon_3 is wrong\n", 
		         checon_3_obj->info );
	}
	if( checon_3_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_checon_3 is wrong\n", 
		checon_3_obj->inforef );
	}

    printf("checon_3_obj->rcond: %lf   checon_3_obj->rcondref: %lf \n", checon_3_obj->rcond, checon_3_obj->rcondref);
	EXPECT_NEAR(0.0, abs(checon_3_obj->rcond - checon_3_obj->rcondref), checon_3_obj->threshold);

}

TEST(checon_3_test, checon_3_1) {}
TEST(checon_3_test, checon_3_2) {}
TEST(checon_3_test, checon_3_3) {}
TEST(checon_3_test, checon_3_4) {}


class hecon_3_dcomplex_parameters{
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
      hecon_3_dcomplex_parameters (int matrix_layout, char uplo, int n);
      ~hecon_3_dcomplex_parameters ();

}; /* end of hecon_3_float_parameters  class definition */

/* Destructor definition **/
hecon_3_dcomplex_parameters:: ~hecon_3_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hecon_3_float_parameters object: destructor invoked. \n");
#endif
   hecon_3_free();
}

hecon_3_dcomplex_parameters::hecon_3_dcomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
    
    
   lda = n;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a,  &aref,  (lda*n));
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e,  &eref, n);
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (e==NULL) || (eref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       hecon_3_free();
       printf(" hecon_3_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    //lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (lda*n));
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(a, aref,
                     n, lda, uplo);
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
class zhecon_3_test  : public  ::testing::Test {
public:
   hecon_3_dcomplex_parameters  *zhecon_3_obj;
   void SetUp();
   void TearDown () { delete zhecon_3_obj; }
};

void zhecon_3_test::SetUp(){

   /* LAPACKE ZHECON_3 prototype */
   typedef int (*Fptr_NL_LAPACKE_zhecon_3) ( int matrix_layout, char uplo, 
                          lapack_int n, const lapack_complex_double* a, 
						  lapack_int lda, const lapack_complex_double* e,
						  const lapack_int* ipiv, 
						  double anorm, double* rcond );
						   
   Fptr_NL_LAPACKE_zhecon_3 ZHECON_3;
   double diff;
   void *hModule, *dModule;

	
   zhecon_3_obj = new hecon_3_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   zhecon_3_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
 
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   ZHECON_3 = (Fptr_NL_LAPACKE_zhecon_3)dlsym(hModule, "LAPACKE_zhecon_3");
   if (NULL == ZHECON_3)
   {
   	  printf("Could not get the symbol -LAPACKE_zhecon_3- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   zhecon_3_obj->info =  LAPACKE_zhetrf_rk (zhecon_3_obj->matrix_layout, zhecon_3_obj->uplo, 
                                 zhecon_3_obj->n,
                                 zhecon_3_obj->a, 
								 zhecon_3_obj->lda, 
                                 zhecon_3_obj->e, 
					             zhecon_3_obj->ipiv );

   zhecon_3_obj->inforef =  LAPACKE_zhetrf_rk (zhecon_3_obj->matrix_layout, zhecon_3_obj->uplo, 
                                 zhecon_3_obj->n,
                                 zhecon_3_obj->aref, 
								 zhecon_3_obj->lda, 
                                 zhecon_3_obj->eref, 
					             zhecon_3_obj->ipivref );

   /* Compute libflame's Lapacke o/p  */  
   zhecon_3_obj->info    = LAPACKE_zhecon_3( zhecon_3_obj->matrix_layout, zhecon_3_obj->uplo, 
                                        zhecon_3_obj->n,
                                        (const lapack_complex_double*)zhecon_3_obj->a, 
										zhecon_3_obj->lda,
                                        (const lapack_complex_double*)zhecon_3_obj->e,										
					                    (const lapack_int *)zhecon_3_obj->ipiv, 
										zhecon_3_obj->anorm, &zhecon_3_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   zhecon_3_obj->inforef = ZHECON_3( zhecon_3_obj->matrix_layout, zhecon_3_obj->uplo, 
                                zhecon_3_obj->n, 
				                (const lapack_complex_double*)zhecon_3_obj->aref, 
								zhecon_3_obj->lda, 
								(const lapack_complex_double*)zhecon_3_obj->eref,
			                    (const lapack_int *)zhecon_3_obj->ipivref, 
								zhecon_3_obj->anorm, &zhecon_3_obj->rcondref);

	if( zhecon_3_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_zhecon_3 is wrong\n", 
		         zhecon_3_obj->info );
	}
	if( zhecon_3_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_zhecon_3 is wrong\n", 
		zhecon_3_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
    printf("zhecon_3_obj->rcond: %lf   zhecon_3_obj->rcondref: %lf \n", zhecon_3_obj->rcond, zhecon_3_obj->rcondref);
	EXPECT_NEAR(0.0, abs(zhecon_3_obj->rcond - zhecon_3_obj->rcondref), zhecon_3_obj->threshold);

}

TEST(zhecon_3_test, zhecon_3_1) {}
TEST(zhecon_3_test, zhecon_3_2) {}
TEST(zhecon_3_test, zhecon_3_3) {}
TEST(zhecon_3_test, zhecon_3_4) {}

