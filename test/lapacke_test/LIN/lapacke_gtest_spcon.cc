#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#include <math.h>
#include "lapacke.h"

#define spcon_free() \
       free (a   ); \
       free (aref); \
       free (ipiv   ); \
       free (ipivref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

class spcon_double_parameters{
   public:
   float threshold;
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
   
   
   public: 
      spcon_double_parameters (int matrix_layout, char uplo, int n );
      ~spcon_double_parameters ();

}; /* end of spcon_double_parameters  class definition */

/* Destructor definition **/
spcon_double_parameters:: ~spcon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spcon_double_parameters object: destructor invoked. \n");
#endif
   spcon_free();
}

spcon_double_parameters::spcon_double_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;    
    
   lda = n_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &a,  &aref,  (n*(n+1)/2));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       spcon_free();
       printf(" spcon_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, (n*(n+1)/2));
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class dspcon_test  : public  ::testing::Test {
public:
   spcon_double_parameters  *dspcon_obj;
   void SetUp();
   void TearDown () { delete dspcon_obj; }
};

void dspcon_test::SetUp(){


   /* LAPACKE DSPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dspcon) ( int matrix_layout, char uplo, 
                  lapack_int n, const double* a, 
			      const lapack_int* ipiv, double anorm, double* rcond );
						   
   Fptr_NL_LAPACKE_dspcon DSPCON;
   double diff;
   void *hModule, *dModule;

   dspcon_obj = new spcon_double_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   dspcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   DSPCON = (Fptr_NL_LAPACKE_dspcon)dlsym(hModule, "LAPACKE_dspcon");
   if (NULL == DSPCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dspcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   /* call SPTRF before calling SPCON  */
    dspcon_obj->info = LAPACKE_dsptrf( dspcon_obj->matrix_layout,dspcon_obj->uplo,
                                                     dspcon_obj->n,dspcon_obj->a,
                                                               dspcon_obj->ipiv);

    dspcon_obj->info = LAPACKE_dsptrf( dspcon_obj->matrix_layout,dspcon_obj->uplo,
                                                     dspcon_obj->n,dspcon_obj->aref,
                                                               dspcon_obj->ipivref);
   
   /* Compute libflame's Lapacke o/p  */
   dspcon_obj->info    = LAPACKE_dspcon( dspcon_obj->matrix_layout, dspcon_obj->uplo, 
                                        dspcon_obj->n, 
                                        (const double * )dspcon_obj->a,
										(const lapack_int *)dspcon_obj->ipiv, 
										dspcon_obj->anorm, &dspcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   dspcon_obj->inforef = DSPCON( dspcon_obj->matrix_layout, dspcon_obj->uplo, 
                                dspcon_obj->n,  
								(const double * )dspcon_obj->aref, 
								(const lapack_int *)dspcon_obj->ipivref, 
								dspcon_obj->anorm, &dspcon_obj->rcondref);

	if( dspcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dspcon is wrong\n", 
		         dspcon_obj->info );
	}
	if( dspcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dspcon is wrong\n", 
		dspcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("dspcon_obj->rcond: %lf   dspcon_obj->rcondref: %lf \n", dspcon_obj->rcond, dspcon_obj->rcondref);
#endif
	EXPECT_NEAR(0.0, abs(dspcon_obj->rcond - dspcon_obj->rcondref), DOUBLE_DIFF_THRESHOLD);
}

TEST(spcon, dspcon1) {}
TEST(spcon, dspcon2) {}
TEST(spcon, dspcon3) {}
TEST(spcon, dspcon4) {}

class spcon_float_parameters{
   public:
   float threshold;
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
   
   
   public: 
      spcon_float_parameters (int matrix_layout, char uplo, int n );
      ~spcon_float_parameters ();

}; /* end of spcon_float_parameters  class definition */

/* Destructor definition **/
spcon_float_parameters:: ~spcon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spcon_float_parameters object: destructor invoked. \n");
#endif
   spcon_free();
}

spcon_float_parameters::spcon_float_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;    
    
   lda = n_i;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &a,  &aref,  (n*(n+1)/2));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       spcon_free();
       printf(" spcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, (n*(n+1)/2));
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class sspcon_test  : public  ::testing::Test {
public:
   spcon_float_parameters  *sspcon_obj;
   void SetUp();
   void TearDown () { delete sspcon_obj; }
};

void sspcon_test::SetUp(){


   /* LAPACKE DSPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_sspcon) ( int matrix_layout, char uplo, 
                  lapack_int n, const float* a, 
			      const lapack_int* ipiv, float anorm, float* rcond );
						   
   Fptr_NL_LAPACKE_sspcon DSPCON;
   float diff;
   void *hModule, *dModule;

   sspcon_obj = new spcon_float_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   sspcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   DSPCON = (Fptr_NL_LAPACKE_sspcon)dlsym(hModule, "LAPACKE_sspcon");
   if (NULL == DSPCON)
   {
   	  printf("Could not get the symbol -LAPACKE_sspcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   /* call SPTRF before calling SPCON  */
    sspcon_obj->info = LAPACKE_ssptrf( sspcon_obj->matrix_layout,sspcon_obj->uplo,
                                                     sspcon_obj->n,sspcon_obj->a,
                                                               sspcon_obj->ipiv);

    sspcon_obj->info = LAPACKE_ssptrf( sspcon_obj->matrix_layout,sspcon_obj->uplo,
                                                     sspcon_obj->n,sspcon_obj->aref,
                                                               sspcon_obj->ipivref);
   
   /* Compute libflame's Lapacke o/p  */
   sspcon_obj->info    = LAPACKE_sspcon( sspcon_obj->matrix_layout, sspcon_obj->uplo, 
                                        sspcon_obj->n, 
                                        (const float * )sspcon_obj->a,
										(const lapack_int *)sspcon_obj->ipiv, 
										sspcon_obj->anorm, &sspcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   sspcon_obj->inforef = DSPCON( sspcon_obj->matrix_layout, sspcon_obj->uplo, 
                                sspcon_obj->n,  
								(const float * )sspcon_obj->aref, 
								(const lapack_int *)sspcon_obj->ipivref, 
								sspcon_obj->anorm, &sspcon_obj->rcondref);

	if( sspcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_sspcon is wrong\n", 
		         sspcon_obj->info );
	}
	if( sspcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_sspcon is wrong\n", 
		sspcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("sspcon_obj->rcond: %lf   sspcon_obj->rcondref: %lf \n", sspcon_obj->rcond, sspcon_obj->rcondref);
#endif
	EXPECT_NEAR(0.0, abs(sspcon_obj->rcond - sspcon_obj->rcondref), DOUBLE_DIFF_THRESHOLD);
}

TEST(spcon, sspcon1) {}
TEST(spcon, sspcon2) {}
TEST(spcon, sspcon3) {}
TEST(spcon, sspcon4) {}

class cpcon_scomplex_parameters{
   public:
   float threshold;
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
   
   
   public: 
      cpcon_scomplex_parameters (int matrix_layout, char uplo, int n );
      ~cpcon_scomplex_parameters ();

}; /* end of cpcon_scomplex_parameters  class definition */

/* Destructor definition **/
cpcon_scomplex_parameters:: ~cpcon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" cpcon_scomplex_parameters object: destructor invoked. \n");
#endif
   spcon_free();
}

cpcon_scomplex_parameters::cpcon_scomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;    
    
   lda = n_i;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a,  &aref,  (n*(n+1)/2));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       spcon_free();
       printf(" cpcon_scomplex_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (n*(n+1)/2));
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class cspcon_test  : public  ::testing::Test {
public:
   cpcon_scomplex_parameters  *cspcon_obj;
   void SetUp();
   void TearDown () { delete cspcon_obj; }
};

void cspcon_test::SetUp(){


   /* LAPACKE CSPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_cspcon) ( int matrix_layout, char uplo, 
                  lapack_int n, const lapack_complex_float* a, 
			      const lapack_int* ipiv, float anorm, float* rcond );
						   
   Fptr_NL_LAPACKE_cspcon CSPCON;
   lapack_complex_float diff;
   void *hModule, *dModule;

   cspcon_obj = new cpcon_scomplex_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   cspcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   CSPCON = (Fptr_NL_LAPACKE_cspcon)dlsym(hModule, "LAPACKE_cspcon");
   if (NULL == CSPCON)
   {
   	  printf("Could not get the symbol -LAPACKE_cspcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   /* call SPTRF before calling SPCON  */
    cspcon_obj->info = LAPACKE_csptrf( cspcon_obj->matrix_layout,cspcon_obj->uplo,
                                                     cspcon_obj->n,cspcon_obj->a,
                                                               cspcon_obj->ipiv);

    cspcon_obj->info = LAPACKE_csptrf( cspcon_obj->matrix_layout,cspcon_obj->uplo,
                                                     cspcon_obj->n,cspcon_obj->aref,
                                                               cspcon_obj->ipivref);
   
   /* Compute libflame's Lapacke o/p  */
   cspcon_obj->info    = LAPACKE_cspcon( cspcon_obj->matrix_layout, cspcon_obj->uplo, 
                                        cspcon_obj->n, 
                                        (const lapack_complex_float * )cspcon_obj->a,
										(const lapack_int *)cspcon_obj->ipiv, 
										cspcon_obj->anorm, &cspcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   cspcon_obj->inforef = CSPCON( cspcon_obj->matrix_layout, cspcon_obj->uplo, 
                                cspcon_obj->n,  
								(const lapack_complex_float * )cspcon_obj->aref, 
								(const lapack_int *)cspcon_obj->ipivref, 
								cspcon_obj->anorm, &cspcon_obj->rcondref);

	if( cspcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_cspcon is wrong\n", 
		         cspcon_obj->info );
	}
	if( cspcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_cspcon is wrong\n", 
		cspcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("cspcon_obj->rcond: %lf   cspcon_obj->rcondref: %lf \n", cspcon_obj->rcond, cspcon_obj->rcondref);
#endif
	EXPECT_NEAR(0.0, abs(cspcon_obj->rcond - cspcon_obj->rcondref), DOUBLE_DIFF_THRESHOLD);
}

TEST(cpcon, cspcon1) {}
TEST(cpcon, cspcon2) {}
TEST(cpcon, cspcon3) {}
TEST(cpcon, cspcon4) {}

class spcon_dcomplex_parameters{
   public:
   float threshold;
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
   
   
   public: 
      spcon_dcomplex_parameters (int matrix_layout, char uplo, int n );
      ~spcon_dcomplex_parameters ();

}; /* end of spcon_dcomplex_parameters  class definition */

/* Destructor definition **/
spcon_dcomplex_parameters:: ~spcon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spcon_dcomplex_parameters object: destructor invoked. \n");
#endif
   spcon_free();
}

spcon_dcomplex_parameters::spcon_dcomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i ){
   
   int j;
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;    
    
   lda = n_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a,  &aref,  (n*(n+1)/2));
   lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
	
	if( (a==NULL) || (aref==NULL) ||  \
	    (ipiv==NULL) || (ipivref==NULL) ){
       spcon_free();
       printf(" spcon_dcomplex_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, (n*(n+1)/2));
	
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class zspcon_test  : public  ::testing::Test {
public:
   spcon_dcomplex_parameters  *zspcon_obj;
   void SetUp();
   void TearDown () { delete zspcon_obj; }
};

void zspcon_test::SetUp(){


   /* LAPACKE ZSPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_zspcon) ( int matrix_layout, char uplo, 
                  lapack_int n, const lapack_complex_double* a, 
			      const lapack_int* ipiv, double anorm, double* rcond );
						   
   Fptr_NL_LAPACKE_zspcon ZSPCON;
   lapack_complex_double diff;
   void *hModule, *dModule;

   zspcon_obj = new spcon_dcomplex_parameters (lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n );

   zspcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   ZSPCON = (Fptr_NL_LAPACKE_zspcon)dlsym(hModule, "LAPACKE_zspcon");
   if (NULL == ZSPCON)
   {
   	  printf("Could not get the symbol -LAPACKE_zspcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   
   /* call SPTRF before calling SPCON  */
    zspcon_obj->info = LAPACKE_zsptrf( zspcon_obj->matrix_layout,zspcon_obj->uplo,
                                                     zspcon_obj->n,zspcon_obj->a,
                                                               zspcon_obj->ipiv);

    zspcon_obj->info = LAPACKE_zsptrf( zspcon_obj->matrix_layout,zspcon_obj->uplo,
                                                     zspcon_obj->n,zspcon_obj->aref,
                                                               zspcon_obj->ipivref);
   
   /* Compute libflame's Lapacke o/p  */
   zspcon_obj->info    = LAPACKE_zspcon( zspcon_obj->matrix_layout, zspcon_obj->uplo, 
                                        zspcon_obj->n, 
                                        (const lapack_complex_double * )zspcon_obj->a,
										(const lapack_int *)zspcon_obj->ipiv, 
										zspcon_obj->anorm, &zspcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   zspcon_obj->inforef = ZSPCON( zspcon_obj->matrix_layout, zspcon_obj->uplo, 
                                zspcon_obj->n,  
								(const lapack_complex_double * )zspcon_obj->aref, 
								(const lapack_int *)zspcon_obj->ipivref, 
								zspcon_obj->anorm, &zspcon_obj->rcondref);

	if( zspcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_zspcon is wrong\n", 
		         zspcon_obj->info );
	}
	if( zspcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_zspcon is wrong\n", 
		zspcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("zspcon_obj->rcond: %lf   zspcon_obj->rcondref: %lf \n", zspcon_obj->rcond, zspcon_obj->rcondref);
#endif
	EXPECT_NEAR(0.0, abs(zspcon_obj->rcond - zspcon_obj->rcondref), DOUBLE_DIFF_THRESHOLD);
}

TEST(spcon, zspcon1) {}
TEST(spcon, zspcon2) {}
TEST(spcon, zspcon3) {}
TEST(spcon, zspcon4) {}
