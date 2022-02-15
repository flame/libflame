#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#include <math.h>
#include "lapacke.h"

#define gecon_free() \
       free (a   ); \
       free (aref);

class gecon_double_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   int 	n;
   double  *a, *aref;
   int 	lda;
   double anorm;
   double rcond, rcondref;   int info, inforef;
   public: 
      gecon_double_parameters (int matrix_layout, char norm, int n,
                               int lda, double anorm );
      ~gecon_double_parameters ();

}; /* end of gecon_double_parameters  class definition */

/* Destructor definition **/
gecon_double_parameters:: ~gecon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gecon_double_parameters object: destructor invoked. \n");
#endif
   gecon_free();
}

gecon_double_parameters::gecon_double_parameters (int matrix_layout_i, 
                             char norm_i, int n_i, int lda_i, 
							 double anorm_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   n = n_i;
   lda = lda_i;
   anorm = anorm_i;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &a,  &aref,  (lda*n));
	
	if( (a==NULL) || (aref==NULL) ){
       gecon_free();
       printf(" gecon_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, (lda*n));
	
	info = 0;
	inforef = 0;
}

TEST(gecon, dgecon1) {

   /* LAPACKE DGECON prototype */
   typedef int (*Fptr_NL_LAPACKE_dgecon) ( int matrix_layout, char norm, 
                          lapack_int n, const double* a, lapack_int lda, 
						  double anorm, double* rcond  );
						  
   Fptr_NL_LAPACKE_dgecon DGECON;
   double diff;
   void *hModule, *dModule;

	
   gecon_double_parameters dgecon_obj(LAPACK_COL_MAJOR, 
                             'I', 512, 512, 121.0 );
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   DGECON = (Fptr_NL_LAPACKE_dgecon)dlsym(hModule, "LAPACKE_dgecon");
   if (NULL == DGECON)
   {
   	  printf("Could not get the symbol -LAPACKE_dgecon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   dgecon_obj.info    = LAPACKE_dgecon( dgecon_obj.matrix_layout, dgecon_obj.norm, 
                                        dgecon_obj.n, (const double*)dgecon_obj.a, dgecon_obj.lda, 
										dgecon_obj.anorm, &dgecon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   dgecon_obj.inforef = DGECON( dgecon_obj.matrix_layout, dgecon_obj.norm, 
                                dgecon_obj.n, (const double*)dgecon_obj.aref, dgecon_obj.lda, 
								dgecon_obj.anorm, &dgecon_obj.rcondref);
 	if( dgecon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dgecon is wrong\n", 
		         dgecon_obj.info );
	}
	if( dgecon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dgecon is wrong\n", 
		dgecon_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
	
	diff =  computeDiff_d( dgecon_obj.lda*dgecon_obj.n, dgecon_obj.a, dgecon_obj.aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
	/*
	FLA_Nrm2_external( dgecon_obj.A, dgecon_obj.norm );
    FLA_Nrm2_external( dgecon_obj.Aref, dgecon_obj. normref );
	diff = FLA_abs( dgecon_obj. norm,  dgecon_obj. normref);
	*/
#if LAPACKE_TEST_VERBOSE
    printf("dgecon_obj.rcond: %lf   dgecon_obj.rcondref: %lf \n", dgecon_obj.rcond, dgecon_obj.rcondref);
	printf("diff: %lf  \n", diff);
#endif
	EXPECT_NEAR(0.0, abs(dgecon_obj.rcond - dgecon_obj.rcondref), DOUBLE_DIFF_THRESHOLD);
}



class gecon_float_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   int 	n;
   int 	kl, ku;
   float  *a, *aref;
   int 	lda;
   float anorm;
   float rcond, rcondref;   int info, inforef;
   public: 
      gecon_float_parameters (int matrix_layout, char norm, int n,
                            int kl, int ku, int lda, float anorm );
      ~gecon_float_parameters ();

}; /* end of gecon_float_parameters  class definition */

/* Destructor definition **/
gecon_float_parameters:: ~gecon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gecon_float_parameters object: destructor invoked. \n");
#endif
   gecon_free();
}

gecon_float_parameters::gecon_float_parameters (int matrix_layout_i, 
                             char norm_i, int n_i, int kl_i, int ku_i,
                             int lda_i, float anorm_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   n = n_i;
   kl = kl_i;
   ku = ku_i;
   lda = lda_i;
   anorm = anorm_i;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &a,  &aref,  (lda*n));
	
	if( (a==NULL) || (aref==NULL)){
       gecon_free();
       printf(" gecon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, (lda*n));
	info = 0;
	inforef = 0;
}

TEST(gecon, sgecon1) {

   /* LAPACKE SGECON prototype */
   typedef int (*Fptr_NL_LAPACKE_sgecon) ( int matrix_layout, char norm, lapack_int n,
                                          const float *a, lapack_int lda, 
										  float anorm, float *rcond );
						   Fptr_NL_LAPACKE_sgecon SGECON;
   float diff;
   void *hModule, *dModule;

	
   gecon_float_parameters sgecon_obj(LAPACK_COL_MAJOR, 
                             'O', 1024, 500, 400,
                             1450, 0.86 );
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   SGECON = (Fptr_NL_LAPACKE_sgecon)dlsym(hModule, "LAPACKE_dgecon");
   if (NULL == SGECON)
   {
   	  printf("Could not get the symbol -LAPACKE_dgecon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   sgecon_obj.info    = LAPACKE_sgecon( sgecon_obj.matrix_layout, sgecon_obj.norm, 
                                        sgecon_obj.n, 
                                        sgecon_obj.a, sgecon_obj.lda, 
										sgecon_obj.anorm,
						                &sgecon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   sgecon_obj.inforef = SGECON( sgecon_obj.matrix_layout, sgecon_obj.norm, 
                                sgecon_obj.n, 
								sgecon_obj.aref, sgecon_obj.lda, 
								sgecon_obj.anorm,
								&sgecon_obj.rcondref);

	if( sgecon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with libflame LAPACKE_sgecon is wrong\n", 
		         sgecon_obj.info );
	}
	if( sgecon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_sgecon is wrong\n", 
		sgecon_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
	
	diff =  computeDiff_s( sgecon_obj.lda*sgecon_obj.n, sgecon_obj.a, sgecon_obj.aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
	/*
	FLA_Nrm2_external( sgecon_obj.A, sgecon_obj.norm );
    FLA_Nrm2_external( sgecon_obj.Aref, sgecon_obj. normref );
	diff = FLA_abs( sgecon_obj. norm,  sgecon_obj. normref);
	*/
#if LAPACKE_TEST_VERBOSE
    printf("sgecon_obj.rcond: %lf   sgecon_obj.rcondref: %lf \n", sgecon_obj.rcond, sgecon_obj.rcondref);
	printf("diff: %f  \n", diff);
#endif
   EXPECT_NEAR(0.0, abs(sgecon_obj.rcond - sgecon_obj.rcondref), DOUBLE_DIFF_THRESHOLD);
}


class gecon_scomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   int 	n;
   int 	kl, ku;
   lapack_complex_float  *a, *aref;
   int 	lda;
   float anorm;
   float rcond, rcondref;   int info, inforef;
   public: 
      gecon_scomplex_parameters (int matrix_layout, char norm, int n,
                            int kl, int ku, int lda, float anorm );
      ~gecon_scomplex_parameters ();

}; /* end of gecon_float_parameters  class definition */

/* Destructor definition **/
gecon_scomplex_parameters:: ~gecon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gecon_float_parameters object: destructor invoked. \n");
#endif
   gecon_free();
}

gecon_scomplex_parameters::gecon_scomplex_parameters (int matrix_layout_i, 
                             char norm_i, int n_i, int kl_i, int ku_i,
                             int lda_i, float anorm_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   n = n_i;
   kl = kl_i;
   ku = ku_i;
   lda = lda_i;
   anorm = anorm_i;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a,  &aref,  (lda*n));
	
	if( (a==NULL) || (aref==NULL)){
       gecon_free();
       printf(" gecon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (lda*n));
	
	info = 0;
	inforef = 0;
}

TEST(gecon, cgecon1) {

   /* LAPACKE CGECON prototype */
   typedef int (*Fptr_NL_LAPACKE_cgecon) ( int matrix_layout, char norm, lapack_int n,
                      const lapack_complex_float *a, 
					  lapack_int lda,
					  float anorm, float *rcond );
						   Fptr_NL_LAPACKE_cgecon CGECON;
   float diff;
   void *hModule, *dModule;

	
   gecon_scomplex_parameters cgecon_obj(LAPACK_COL_MAJOR, 
                             '1', 256, 125, 400,
                             1450, 11.0 );
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   CGECON = (Fptr_NL_LAPACKE_cgecon)dlsym(hModule, "LAPACKE_cgecon");
   if (NULL == CGECON)
   {
   	  printf("Could not get the symbol -LAPACKE_cgecon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */  
   cgecon_obj.info    = LAPACKE_cgecon( cgecon_obj.matrix_layout, cgecon_obj.norm, 
                                        cgecon_obj.n,
                                        cgecon_obj.a, cgecon_obj.lda, 
					                    cgecon_obj.anorm,
			                            &cgecon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   cgecon_obj.inforef = CGECON( cgecon_obj.matrix_layout, cgecon_obj.norm, 
                                cgecon_obj.n,
				                cgecon_obj.aref, cgecon_obj.lda, 
			                    cgecon_obj.anorm,
				                &cgecon_obj.rcondref);

	if( cgecon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_cgecon is wrong\n", 
		         cgecon_obj.info );
	}
	if( cgecon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_cgecon is wrong\n", 
		cgecon_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
	
	diff =  computeDiff_c( cgecon_obj.lda*cgecon_obj.n, cgecon_obj.a, cgecon_obj.aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
	/*
	FLA_Nrm2_external( cgecon_obj.A, cgecon_obj.norm );
    FLA_Nrm2_external( cgecon_obj.Aref, cgecon_obj. normref );
	diff = FLA_abs( cgecon_obj. norm,  cgecon_obj. normref);
	*/
#if LAPACKE_TEST_VERBOSE
    printf("cgecon_obj.rcond: %lf   cgecon_obj.rcondref: %lf \n", cgecon_obj.rcond, cgecon_obj.rcondref);
	printf("diff: %f  \n", diff);
#endif
	EXPECT_NEAR(0.0, abs(cgecon_obj.rcond - cgecon_obj.rcondref), DOUBLE_DIFF_THRESHOLD);

}

class gecon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   int 	n;
   int 	kl, ku;
   lapack_complex_double  *a, *aref;
   int 	lda;
   double anorm;
   double rcond, rcondref;   int info, inforef;
   public: 
      gecon_dcomplex_parameters (int matrix_layout, char norm, int n,
                            int kl, int ku, int lda, double anorm );
      ~gecon_dcomplex_parameters ();

}; /* end of gecon_float_parameters  class definition */

/* Destructor definition **/
gecon_dcomplex_parameters:: ~gecon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" gecon_float_parameters object: destructor invoked. \n");
#endif
   gecon_free();
}

gecon_dcomplex_parameters::gecon_dcomplex_parameters (int matrix_layout_i, 
                             char norm_i, int n_i, int kl_i, int ku_i,
                             int lda_i, double anorm_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   n = n_i;
   kl = kl_i;
   ku = ku_i;
   lda = lda_i;
   anorm = anorm_i;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a,  &aref,  (lda*n));
	
	if( (a==NULL) || (aref==NULL) ){
       gecon_free();
       printf(" gecon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, (lda*n));
	
	info = 0;
	inforef = 0;
}

TEST(gecon, zgecon1) {

   /* LAPACKE ZGECON prototype */
   typedef int (*Fptr_NL_LAPACKE_zgecon) ( int matrix_layout, char norm, lapack_int n,
                      const lapack_complex_double *a, 
					  lapack_int lda,
					  double anorm, double *rcond );
						   Fptr_NL_LAPACKE_zgecon ZGECON;
   double diff;
   void *hModule, *dModule;

	
   gecon_dcomplex_parameters zgecon_obj(LAPACK_COL_MAJOR, 
                             'O', 512, 500, 400,
                             1450, 167.0 );
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   ZGECON = (Fptr_NL_LAPACKE_zgecon)dlsym(hModule, "LAPACKE_zgecon");
   if (NULL == ZGECON)
   {
   	  printf("Could not get the symbol -LAPACKE_zgecon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */  
   zgecon_obj.info    = LAPACKE_zgecon( zgecon_obj.matrix_layout, zgecon_obj.norm, 
                                        zgecon_obj.n,
                                        zgecon_obj.a, zgecon_obj.lda, 
					                    zgecon_obj.anorm,
			                            &zgecon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   zgecon_obj.inforef = ZGECON( zgecon_obj.matrix_layout, zgecon_obj.norm, 
                                zgecon_obj.n,
				                zgecon_obj.aref, zgecon_obj.lda, 
			                    zgecon_obj.anorm,
				                &zgecon_obj.rcondref);

	if( zgecon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_cgecon is wrong\n", 
		         zgecon_obj.info );
	}
	if( zgecon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_cgecon is wrong\n", 
		zgecon_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
	
	diff =  computeDiff_z( zgecon_obj.lda*zgecon_obj.n, zgecon_obj.a, zgecon_obj.aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	
	/*
	FLA_Nrm2_external( zgecon_obj.A, zgecon_obj.norm );
    FLA_Nrm2_external( zgecon_obj.Aref, zgecon_obj. normref );
	diff = FLA_abs( zgecon_obj. norm,  zgecon_obj. normref);
	*/
#if LAPACKE_TEST_VERBOSE
    printf("zgecon_obj.rcond: %lf   zgecon_obj.rcondref: %lf \n", zgecon_obj.rcond, zgecon_obj.rcondref);
	printf("diff: %lf  \n", diff);
#endif
	EXPECT_NEAR(0.0, abs(zgecon_obj.rcond - zgecon_obj.rcondref), DOUBLE_DIFF_THRESHOLD);

}


