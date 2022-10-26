#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#include <math.h>
#include "lapacke.h"

#define ptcon_free() \
       free (d   ); \
       free (dref); \
       free (e   ); \
       free (eref);

class ptcon_double_parameters{
   public:
   /* input params to the API **/
   int 	n;
   double  *d, *dref;
   double  *e, *eref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   public: 
      ptcon_double_parameters (int n, double anorm );
      ~ptcon_double_parameters ();

}; /* end of ptcon_double_parameters  class definition */

/* Destructor definition **/
ptcon_double_parameters:: ~ptcon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptcon_double_parameters object: destructor invoked. \n");
#endif
   ptcon_free();
}

ptcon_double_parameters::ptcon_double_parameters ( int n_i, double anormi ){
   n = n_i;
   anorm = anormi;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &d,  &dref,  n);
   lapacke_gtest_alloc_double_buffer_pair( &e,  &eref,  n-1);
	if( (d==NULL) || (dref==NULL) ||  \
	    (e==NULL) || (eref==NULL) ){
	   ptcon_free();
       printf(" ptcon_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_double_buffer_pair_rand( e, eref, n-1 );
	
	info = 0;
	inforef = 0;
}

TEST(ptcon, dptcon1) {

   /* LAPACKE DPTCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dptcon) ( lapack_int n, const double* d,
                          const double* e, double anorm, double* rcond  );
						  
   Fptr_NL_LAPACKE_dptcon DPTCON;
   double diff;
   void *hModule, *dModule;

	
   ptcon_double_parameters dptcon_obj( 512, 121.0 );
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   DPTCON = (Fptr_NL_LAPACKE_dptcon)dlsym(hModule, "LAPACKE_dptcon");
   if (NULL == DPTCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dptcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   dptcon_obj.info    = LAPACKE_dptcon( dptcon_obj.n, (const double*)dptcon_obj.d, 
                                        (const double*)dptcon_obj.e, 
										dptcon_obj.anorm, &dptcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   dptcon_obj.inforef = DPTCON(  dptcon_obj.n, (const double*)dptcon_obj.dref, 
                                        (const double*)dptcon_obj.eref, 
										dptcon_obj.anorm, &dptcon_obj.rcondref);
 	if( dptcon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dptcon is wrong\n", 
		         dptcon_obj.info );
	}
	if( dptcon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dptcon is wrong\n", 
		dptcon_obj.inforef );
	}

#if LAPACKE_TEST_VERBOSE
    printf("dptcon_obj.rcond: %lf   dptcon_obj.rcondref: %lf \n", dptcon_obj.rcond, dptcon_obj.rcondref);
#endif
	EXPECT_NEAR(0.0, abs(dptcon_obj.rcond - dptcon_obj.rcondref), DOUBLE_DIFF_THRESHOLD);
}



class ptcon_float_parameters{
   public:
   /* input params to the API **/
   int 	n;
   float  *d, *dref;
   float  *e, *eref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   public: 
      ptcon_float_parameters (int n, float anorm  );
      ~ptcon_float_parameters ();

}; /* end of ptcon_float_parameters  class definition */

/* Destructor definition **/
ptcon_float_parameters:: ~ptcon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptcon_float_parameters object: destructor invoked. \n");
#endif
   ptcon_free();
}

ptcon_float_parameters::ptcon_float_parameters ( int n_i, float anormi ){
   n = n_i;
   anorm = anormi;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &d,  &dref,  n);
   lapacke_gtest_alloc_float_buffer_pair( &e,  &eref,  n-1);
	if( (d==NULL) || (dref==NULL) ||  \
	    (e==NULL) || (eref==NULL) ){
	   ptcon_free();
       printf(" ptcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_float_buffer_pair_rand( e, eref, n-1 );
	
	info = 0;
	inforef = 0;
}

TEST(ptcon, sptcon1) {

   /* LAPACKE SPTCON prototype */
   typedef int (*Fptr_NL_LAPACKE_sptcon) ( lapack_int n, const float* d,
                          const float* e, float anorm, float* rcond  );
						  
   Fptr_NL_LAPACKE_sptcon SPTCON;
   float diff;
   void *hModule, *dModule;

	
   ptcon_float_parameters sptcon_obj( 512, 121.0 );
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   SPTCON = (Fptr_NL_LAPACKE_sptcon)dlsym(hModule, "LAPACKE_sptcon");
   if (NULL == SPTCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dptcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   sptcon_obj.info    = LAPACKE_sptcon( sptcon_obj.n, (const float*)sptcon_obj.d, 
                                        (const float*)sptcon_obj.e, 
										sptcon_obj.anorm, &sptcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   sptcon_obj.inforef = SPTCON(  sptcon_obj.n, (const float*)sptcon_obj.dref, 
                                        (const float*)sptcon_obj.eref, 
										sptcon_obj.anorm, &sptcon_obj.rcondref);
 	if( sptcon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dptcon is wrong\n", 
		         sptcon_obj.info );
	}
	if( sptcon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dptcon is wrong\n", 
		sptcon_obj.inforef );
	}

#if LAPACKE_TEST_VERBOSE
    printf("sptcon_obj.rcond: %lf   sptcon_obj.rcondref: %lf \n", sptcon_obj.rcond, sptcon_obj.rcondref);
#endif
	EXPECT_NEAR(0.0, fabs(sptcon_obj.rcond - sptcon_obj.rcondref), FLOAT_DIFF_THRESHOLD);
}

class ptcon_scomplex_parameters{
   public:
   /* input params to the API **/
   int 	n;
   float  *d, *dref;
   lapack_complex_float  *e, *eref;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   public: 
      ptcon_scomplex_parameters (int n, float anorm  );
      ~ptcon_scomplex_parameters ();

}; /* end of ptcon_scomplex_parameters  class definition */

/* Destructor definition **/
ptcon_scomplex_parameters:: ~ptcon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptcon_scomplex_parameters object: destructor invoked. \n");
#endif
   ptcon_free();
}

ptcon_scomplex_parameters::ptcon_scomplex_parameters ( int n_i, float anormi ){
   n = n_i;
   anorm = anormi;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &d,  &dref,  n);
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e,  &eref,  n-1);
	if( (d==NULL) || (dref==NULL) ||  \
	    (e==NULL) || (eref==NULL) ){
	   ptcon_free();
       printf(" ptcon_scomplex_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( e, eref, n-1 );
	
	info = 0;
	inforef = 0;
}

TEST(ptcon, cptcon1) {

   /* LAPACKE CPTCON prototype */
   typedef int (*Fptr_NL_LAPACKE_cptcon) ( lapack_int n, const float* d, 
            const lapack_complex_float* e, float anorm, float* rcond  );
						  
   Fptr_NL_LAPACKE_cptcon CPTCON;
   float diff;
   void *hModule, *dModule;

	
   ptcon_scomplex_parameters cptcon_obj( 256, 12.87 );
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   CPTCON = (Fptr_NL_LAPACKE_cptcon)dlsym(hModule, "LAPACKE_cptcon");
   if (NULL == CPTCON)
   {
   	  printf("Could not get the symbol -LAPACKE_cptcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   cptcon_obj.info    = LAPACKE_cptcon( cptcon_obj.n, (const float*)cptcon_obj.d, 
                                        (lapack_complex_float*)cptcon_obj.e, 
										cptcon_obj.anorm, &cptcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   cptcon_obj.inforef = CPTCON(  cptcon_obj.n, (const float*)cptcon_obj.dref, 
                                        (const lapack_complex_float*)cptcon_obj.eref, 
										cptcon_obj.anorm, &cptcon_obj.rcondref);
 	if( cptcon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dptcon is wrong\n", 
		         cptcon_obj.info );
	}
	if( cptcon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dptcon is wrong\n", 
		cptcon_obj.inforef );
	}

#if LAPACKE_TEST_VERBOSE
    printf("cptcon_obj.rcond: %lf   cptcon_obj.rcondref: %lf \n", cptcon_obj.rcond, cptcon_obj.rcondref);
#endif
	EXPECT_NEAR(0.0, fabs(cptcon_obj.rcond - cptcon_obj.rcondref), FLOAT_DIFF_THRESHOLD);
}


class ptcon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	n;
   double  *d, *dref;
   lapack_complex_double  *e, *eref;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   public: 
      ptcon_dcomplex_parameters (int n, double anorm  );
      ~ptcon_dcomplex_parameters ();

}; /* end of ptcon_dcomplex_parameters  class definition */

/* Destructor definition **/
ptcon_dcomplex_parameters:: ~ptcon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptcon_dcomplex_parameters object: destructor invoked. \n");
#endif
   ptcon_free();
}

ptcon_dcomplex_parameters::ptcon_dcomplex_parameters ( int n_i, double anormi ){
   n = n_i;
   anorm = anormi;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &d,  &dref,  n);
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e,  &eref,  n-1);
	if( (d==NULL) || (dref==NULL) ||  \
	    (e==NULL) || (eref==NULL) ){
	   ptcon_free();
       printf(" ptcon_dcomplex_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( e, eref, n-1 );
	
	info = 0;
	inforef = 0;
}

TEST(ptcon, zptcon1) {

   /* LAPACKE ZPTCON prototype */
   typedef int (*Fptr_NL_LAPACKE_zptcon) ( lapack_int n, const double* d, 
            const lapack_complex_double* e, double anorm, double* rcond  );
						  
   Fptr_NL_LAPACKE_zptcon ZPTCON;
   double diff;
   void *hModule, *dModule;

	
   ptcon_dcomplex_parameters zptcon_obj( 256, 12.87 );
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   ZPTCON = (Fptr_NL_LAPACKE_zptcon)dlsym(hModule, "LAPACKE_zptcon");
   if (NULL == ZPTCON)
   {
   	  printf("Could not get the symbol -LAPACKE_zptcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   zptcon_obj.info    = LAPACKE_zptcon( zptcon_obj.n, (const double*)zptcon_obj.d, 
                                        (lapack_complex_double*)zptcon_obj.e, 
										zptcon_obj.anorm, &zptcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   zptcon_obj.inforef = ZPTCON(  zptcon_obj.n, (const double*)zptcon_obj.dref, 
                                        (lapack_complex_double*)zptcon_obj.eref, 
										zptcon_obj.anorm, &zptcon_obj.rcondref);
 	if( zptcon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dptcon is wrong\n", 
		         zptcon_obj.info );
	}
	if( zptcon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dptcon is wrong\n", 
		zptcon_obj.inforef );
	}

#if LAPACKE_TEST_VERBOSE
    printf("zptcon_obj.rcond: %lf   zptcon_obj.rcondref: %lf \n", zptcon_obj.rcond, zptcon_obj.rcondref);
#endif
	EXPECT_NEAR(0.0, fabs(zptcon_obj.rcond - zptcon_obj.rcondref), FLOAT_DIFF_THRESHOLD);
}
