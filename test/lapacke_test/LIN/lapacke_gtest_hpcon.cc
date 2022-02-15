#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#include <math.h>
#include "lapacke.h"

#define hpcon_free() \
       free (a   ); \
       free (aref); \
       free (ipiv   ); \
       free (ipivref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin hpcon_scomplex_parameters  class definition */
class hpcon_scomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   lapack_complex_float  *a,*aref;
   int 	lda;
   int     *ipiv,*ipivref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      hpcon_scomplex_parameters (int matrix_layout, char uplo, int n);
      ~hpcon_scomplex_parameters ();

}; /* end of hpcon_float_parameters  class definition */

/* Destructor definition **/
hpcon_scomplex_parameters:: ~hpcon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hpcon_float_parameters object: destructor invoked. \n");
#endif
   hpcon_free();
}

hpcon_scomplex_parameters::hpcon_scomplex_parameters (int matrix_layout_i,
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*(n+1)/2));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       hpcon_free();
       printf(" hpcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,(n*(n+1)/2));
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class chpcon_test  : public  ::testing::Test {
public:
   hpcon_scomplex_parameters  *chpcon_obj;
   void SetUp();
   void TearDown () { delete chpcon_obj; }
};

void chpcon_test::SetUp(){

   /* LAPACKE CHPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_chpcon) ( int matrix_layout,char uplo,
             lapack_int n,const lapack_complex_float* a,
			 const lapack_int* ipiv,float anorm,float* rcond );
						   Fptr_NL_LAPACKE_chpcon CHPCON;
   float diff;
   void *hModule,*dModule;
	
   chpcon_obj = new hpcon_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   chpcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   CHPCON = (Fptr_NL_LAPACKE_chpcon)dlsym(hModule,"LAPACKE_chpcon");
   if (NULL == CHPCON)
   {
   	  printf("Could not get the symbol -LAPACKE_chpcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }

   LAPACKE_chptrf ( chpcon_obj->matrix_layout, chpcon_obj->uplo, 
                                        chpcon_obj->n,
                                        chpcon_obj->a, 
					                    chpcon_obj->ipiv);

   LAPACKE_chptrf ( chpcon_obj->matrix_layout, chpcon_obj->uplo, 
                                        chpcon_obj->n,
                                        chpcon_obj->aref, 
					                    chpcon_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */  
   chpcon_obj->info    = LAPACKE_chpcon( chpcon_obj->matrix_layout,chpcon_obj->uplo,
                                        chpcon_obj->n,
                                        (const lapack_complex_float*)chpcon_obj->a,
					                    (const lapack_int *)chpcon_obj->ipiv,chpcon_obj->anorm,
			                            &chpcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   chpcon_obj->inforef = CHPCON( chpcon_obj->matrix_layout,chpcon_obj->uplo,
                                chpcon_obj->n,
				                (const lapack_complex_float*)chpcon_obj->aref,
			                    (const lapack_int *)chpcon_obj->ipivref,chpcon_obj->anorm,
				                &chpcon_obj->rcondref);

	if( chpcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_chpcon is wrong\n",
		         chpcon_obj->info );
	}
	if( chpcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_chpcon is wrong\n",
		chpcon_obj->inforef );
	}

#if LAPACKE_TEST_VERBOSE
    printf("chpcon_obj->rcond: %lf   chpcon_obj->rcondref: %lf \n",chpcon_obj->rcond,chpcon_obj->rcondref);
#endif
	EXPECT_NEAR(0.0,abs(chpcon_obj->rcond - chpcon_obj->rcondref),DOUBLE_DIFF_THRESHOLD);

}

TEST(chpcon_test, chpcon1) {}
TEST(chpcon_test, chpcon2) {}
TEST(chpcon_test, chpcon3) {}
TEST(chpcon_test, chpcon4) {}

class hpcon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char uplo;
   int 	n;
   lapack_complex_double  *a,*aref;
   int 	lda;
   int     *ipiv,*ipivref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   float threshold;
   
   public: 
      hpcon_dcomplex_parameters (int matrix_layout, char uplo, int n);
      ~hpcon_dcomplex_parameters ();

}; /* end of hpcon_float_parameters  class definition */

/* Destructor definition **/
hpcon_dcomplex_parameters:: ~hpcon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hpcon_float_parameters object: destructor invoked. \n");
#endif
   hpcon_free();
}

hpcon_dcomplex_parameters::hpcon_dcomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   lda = n;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*(n+1)/2));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv,&ipivref,n);
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       hpcon_free();
       printf(" hpcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,(n*(n+1)/2));
    /* create pivot buffer with random values within the range 0 to n-1. **/
    for( j = 0; j < n; j++ ) {
       ipiv[j] = 0;
       ipivref[j] = 0;
    }
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class zhpcon_test  : public  ::testing::Test {
public:
   hpcon_dcomplex_parameters  *zhpcon_obj;
   void SetUp();
   void TearDown () { delete zhpcon_obj; }
};

void zhpcon_test::SetUp(){

   /* LAPACKE ZHPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_zhpcon) ( int matrix_layout,char uplo,
                          lapack_int n,const lapack_complex_double* a,
						  const lapack_int* ipiv,
						  double anorm,double* rcond );
						   Fptr_NL_LAPACKE_zhpcon ZHPCON;
   double diff;
   void *hModule,*dModule;

	
   zhpcon_obj = new hpcon_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   zhpcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);
   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   ZHPCON = (Fptr_NL_LAPACKE_zhpcon)dlsym(hModule,"LAPACKE_zhpcon");
   if (NULL == ZHPCON)
   {
   	  printf("Could not get the symbol -LAPACKE_zhpcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   LAPACKE_zhptrf ( zhpcon_obj->matrix_layout, zhpcon_obj->uplo, 
                                        zhpcon_obj->n,
                                        zhpcon_obj->a, 
					                    zhpcon_obj->ipiv);

   LAPACKE_zhptrf ( zhpcon_obj->matrix_layout, zhpcon_obj->uplo, 
                                        zhpcon_obj->n,
                                        zhpcon_obj->aref, 
					                    zhpcon_obj->ipivref);

   /* Compute libflame's Lapacke o/p  */  
   zhpcon_obj->info    = LAPACKE_zhpcon( zhpcon_obj->matrix_layout,zhpcon_obj->uplo,
                                        zhpcon_obj->n,
                                        (const lapack_complex_double*)zhpcon_obj->a,
					                    (const lapack_int *)zhpcon_obj->ipiv,zhpcon_obj->anorm,
			                            &zhpcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   zhpcon_obj->inforef = ZHPCON( zhpcon_obj->matrix_layout,zhpcon_obj->uplo,
                                zhpcon_obj->n,
				                (const lapack_complex_double*)zhpcon_obj->aref,
			                    (const lapack_int *)zhpcon_obj->ipivref,zhpcon_obj->anorm,
				                &zhpcon_obj->rcondref);

	if( zhpcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_zhpcon is wrong\n",
		         zhpcon_obj->info );
	}
	if( zhpcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_zhpcon is wrong\n",
		zhpcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("zhpcon_obj->rcond: %lf   zhpcon_obj->rcondref: %lf \n",zhpcon_obj->rcond,zhpcon_obj->rcondref);
#endif
	EXPECT_NEAR(0.0,abs(zhpcon_obj->rcond - zhpcon_obj->rcondref),DOUBLE_DIFF_THRESHOLD);

}

TEST(zhpcon_test, zhpcon1) {}
TEST(zhpcon_test, zhpcon2) {}
TEST(zhpcon_test, zhpcon3) {}
TEST(zhpcon_test, zhpcon4) {}
