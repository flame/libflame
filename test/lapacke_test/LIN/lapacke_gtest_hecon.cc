#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#include <math.h>
#include "lapacke.h"
#include <stdio.h>
#define hecon_free() \
       free (a   ); \
       free (aref); \
       free (ipiv   ); \
       free (ipivref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin hecon_scomplex_parameters  class definition */
class hecon_scomplex_parameters{
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
      hecon_scomplex_parameters (int matrix_layout, char uplo, int n);
      ~hecon_scomplex_parameters ();

}; /* end of hecon_float_parameters  class definition */

/* Destructor definition **/
hecon_scomplex_parameters:: ~hecon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hecon_float_parameters object: destructor invoked. \n");
#endif
   hecon_free();
}

hecon_scomplex_parameters::hecon_scomplex_parameters (int matrix_layout_i, 
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
       hecon_free();
       printf(" hecon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    //lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (lda*n));
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(a, aref,
                     n, lda, uplo);
					 
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }

    /* create pivot buffer with random values within the range 0 to n-1. 
    for(int j=0; j < n; j++) {
       int randIndex = (int) (rand() %  n);
       int tmp = ipiv[j];
       ipiv[j] = ipiv[randIndex];
       ipiv[randIndex] = tmp;

    }
    for( j = 0; j < n; j++ ) {
       ipivref[j] = ipiv[j];
    }**/
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class checon_test  : public  ::testing::Test {
public:
   hecon_scomplex_parameters  *checon_obj;
   void SetUp();
   void TearDown () { delete checon_obj; }
};

void checon_test::SetUp(){

   /* LAPACKE CHECON prototype */
   typedef int (*Fptr_NL_LAPACKE_checon) ( int matrix_layout, char uplo, 
             lapack_int n, const lapack_complex_float* a, lapack_int lda, 
			 const lapack_int* ipiv, float anorm, float* rcond );
						   
   Fptr_NL_LAPACKE_checon CHECON;
   float diff;
   void *hModule, *dModule;
	
   checon_obj = new hecon_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   checon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
 
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   CHECON = (Fptr_NL_LAPACKE_checon)dlsym(hModule, "LAPACKE_checon");
   if (NULL == CHECON)
   {
   	  printf("Could not get the symbol -LAPACKE_checon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   LAPACKE_chetrf ( checon_obj->matrix_layout, checon_obj->uplo, 
                                        checon_obj->n,
                                        checon_obj->a, 
										checon_obj->lda, 
					                    checon_obj->ipiv);

   LAPACKE_chetrf ( checon_obj->matrix_layout, checon_obj->uplo, 
                                        checon_obj->n,
                                        checon_obj->aref, 
										checon_obj->lda, 
					                    checon_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */  
   checon_obj->info    = LAPACKE_checon( checon_obj->matrix_layout, checon_obj->uplo, 
                                        checon_obj->n,
                                        (const lapack_complex_float*)checon_obj->a, 
										checon_obj->lda, 
					                    (const lapack_int *)checon_obj->ipiv, checon_obj->anorm,
			                            &checon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   checon_obj->inforef = CHECON( checon_obj->matrix_layout, checon_obj->uplo, 
                                checon_obj->n, 
				                (const lapack_complex_float*)checon_obj->aref, 
								checon_obj->lda, 
			                    (const lapack_int *)checon_obj->ipivref, checon_obj->anorm,
				                &checon_obj->rcondref);

	if( checon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_checon is wrong\n", 
		         checon_obj->info );
	}
	if( checon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_checon is wrong\n", 
		checon_obj->inforef );
	}

    printf("checon_obj->rcond: %lf   checon_obj->rcondref: %lf \n", checon_obj->rcond, checon_obj->rcondref);
	EXPECT_NEAR(0.0, abs(checon_obj->rcond - checon_obj->rcondref), checon_obj->threshold);
}

TEST(checon_test, checon1) {}
TEST(checon_test, checon2) {}
TEST(checon_test, checon3) {}
TEST(checon_test, checon4) {}

class hecon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   
   lapack_complex_double  *a, *aref;
   int 	lda;
   int  *ipiv, *ipivref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      hecon_dcomplex_parameters (int matrix_layout, char uplo, int n);
      ~hecon_dcomplex_parameters ();

}; /* end of hecon_float_parameters  class definition */

/* Destructor definition **/
hecon_dcomplex_parameters:: ~hecon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hecon_float_parameters object: destructor invoked. \n");
#endif
   hecon_free();
}

hecon_dcomplex_parameters::hecon_dcomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a,  &aref,  (lda*n));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       hecon_free();
       printf(" hecon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(a, aref,
                     n, lda, uplo);
					 
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class zhecon_test  : public  ::testing::Test {
public:
   hecon_dcomplex_parameters  *zhecon_obj;
   void SetUp();
   void TearDown () { delete zhecon_obj; }
};

void zhecon_test::SetUp(){

   /* LAPACKE ZHECON prototype */
   typedef int (*Fptr_NL_LAPACKE_zhecon) ( int matrix_layout, char uplo, 
                          lapack_int n, const lapack_complex_double* a, 
						  lapack_int lda, const lapack_int* ipiv, 
						  double anorm, double* rcond );
						   
   Fptr_NL_LAPACKE_zhecon ZHECON;
   float diff;
   void *hModule, *dModule;

	
   zhecon_obj = new hecon_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   zhecon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   ZHECON = (Fptr_NL_LAPACKE_zhecon)dlsym(hModule, "LAPACKE_zhecon");
   if (NULL == ZHECON)
   {
   	  printf("Could not get the symbol -LAPACKE_zhecon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   LAPACKE_zhetrf ( zhecon_obj->matrix_layout, zhecon_obj->uplo, 
                                        zhecon_obj->n,
                                        zhecon_obj->a, 
										zhecon_obj->lda, 
					                    zhecon_obj->ipiv);

   LAPACKE_zhetrf ( zhecon_obj->matrix_layout, zhecon_obj->uplo, 
                                        zhecon_obj->n,
                                        zhecon_obj->aref, 
										zhecon_obj->lda, 
					                    zhecon_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */  
   zhecon_obj->info    = LAPACKE_zhecon( zhecon_obj->matrix_layout, zhecon_obj->uplo, 
                                        zhecon_obj->n,
                                        (const lapack_complex_double*)zhecon_obj->a, 
										zhecon_obj->lda, 
					                    (const lapack_int *)zhecon_obj->ipiv, zhecon_obj->anorm,
			                            &zhecon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   zhecon_obj->inforef = ZHECON( zhecon_obj->matrix_layout, zhecon_obj->uplo, 
                                zhecon_obj->n, 
				                (const lapack_complex_double*)zhecon_obj->aref, 
								zhecon_obj->lda, 
			                    (const lapack_int *)zhecon_obj->ipivref, zhecon_obj->anorm,
				                &zhecon_obj->rcondref);

	if( zhecon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_zhecon is wrong\n", 
		         zhecon_obj->info );
	}
	if( zhecon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_zhecon is wrong\n", 
		zhecon_obj->inforef );
	}

	/* Compute Difference in the o/ps */
    printf("zhecon_obj->rcond: %lf   zhecon_obj->rcondref: %lf \n", zhecon_obj->rcond, zhecon_obj->rcondref);

    /* validation of o/p */    
	EXPECT_NEAR(0.0, abs(zhecon_obj->rcond - zhecon_obj->rcondref), zhecon_obj->threshold);

}

TEST(zhecon_test, zhecon1) {}
TEST(zhecon_test, zhecon2) {}
TEST(zhecon_test, zhecon3) {}
TEST(zhecon_test, zhecon4) {}
