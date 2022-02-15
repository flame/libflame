#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#include <math.h>
#include "lapacke.h"

#define pbcon_free() \
       free (ab   ); \
       free (abref)

class pbcon_double_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   
   char uplo;
   int 	n;
   int 	kd, ku;
   double  *ab, *abref;
   int 	ldab;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   
   
   public: 
      pbcon_double_parameters (int matrix_layout, char uplo, int n,
                            int kd, int ku, int ldab, double anorm );
      ~pbcon_double_parameters ();

}; /* end of pbcon_double_parameters  class definition */

/* Destructor definition **/
pbcon_double_parameters:: ~pbcon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbcon_double_parameters object: destructor invoked. \n");
#endif
   pbcon_free();
}

pbcon_double_parameters::pbcon_double_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i, int kd_i, int ku_i,
                             int ldab_i, double anorm_i ){
   int j;
   
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   kd = kd_i;
   ku = ku_i;
   ldab = ldab_i;
   anorm = anorm_i;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &ab,  &abref,  (ldab*n));
	
	if( (ab==NULL) || (abref==NULL)){
       pbcon_free();
       printf(" pbcon_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ab, abref, (ldab*n));
	
	info = 0;
	inforef = 0;
}

TEST(pbcon, dpbcon1) {

   /* LAPACKE DPBCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dpbcon) ( int matrix_layout, char uplo, int n,
                                           int kd, const double* ab, int ldab, 
										   double anorm, double* rcond );
						   
   Fptr_NL_LAPACKE_dpbcon DPBCON;
   double diff;
   void *hModule, *dModule;

	
   pbcon_double_parameters dpbcon_obj(LAPACK_COL_MAJOR, 
                             'U', 512, 500, 400,
                             1450, 1.0 );
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   DPBCON = (Fptr_NL_LAPACKE_dpbcon)dlsym(hModule, "LAPACKE_dpbcon");
   if (NULL == DPBCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dpbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   dpbcon_obj.info    = LAPACKE_dpbcon( dpbcon_obj.matrix_layout, dpbcon_obj.uplo,
                                        dpbcon_obj.n, dpbcon_obj.kd,
                                        dpbcon_obj.ab, dpbcon_obj.ldab, 
										dpbcon_obj.anorm, &dpbcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   dpbcon_obj.inforef = DPBCON( dpbcon_obj.matrix_layout, dpbcon_obj.uplo, 
                                dpbcon_obj.n, dpbcon_obj.kd, 
								dpbcon_obj.abref, dpbcon_obj.ldab, 
								dpbcon_obj.anorm, &dpbcon_obj.rcondref);

	if( dpbcon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dpbcon is wrong\n", 
		         dpbcon_obj.info );
	}
	if( dpbcon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dpbcon is wrong\n", 
		dpbcon_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("dpbcon_obj.rcond: %lf   dpbcon_obj.rcondref: %lf \n", dpbcon_obj.rcond, dpbcon_obj.rcondref);
#endif
	EXPECT_NEAR(0.0, abs(dpbcon_obj.rcond - dpbcon_obj.rcondref), DOUBLE_DIFF_THRESHOLD);
}



class pbcon_float_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   
   char uplo;
   int 	n;
   int 	kd, ku;
   float  *ab, *abref;
   int 	ldab;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   
   
   public: 
      pbcon_float_parameters (int matrix_layout, char uplo, int n,
                            int kd, int ku, int ldab, float anorm );
      ~pbcon_float_parameters ();

}; /* end of pbcon_float_parameters  class definition */

/* Destructor definition **/
pbcon_float_parameters:: ~pbcon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbcon_float_parameters object: destructor invoked. \n");
#endif
   pbcon_free();
}

pbcon_float_parameters::pbcon_float_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i, int kd_i, int ku_i,
                             int ldab_i, float anorm_i ){
   int j;
   
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   kd = kd_i;
   ku = ku_i;
   ldab = ldab_i;
   anorm = anorm_i;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &ab,  &abref,  (ldab*n));
	
	if( (ab==NULL) || (abref==NULL) ){
       pbcon_free();
       printf(" pbcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ab, abref, (ldab*n));
	
	info = 0;
	inforef = 0;
}

TEST(pbcon, spbcon1) {

   /* LAPACKE SPBCON prototype */
   typedef int (*Fptr_NL_LAPACKE_spbcon) ( int matrix_layout, char uplo, 
                                          lapack_int n, lapack_int kd, 
										  const float* ab, lapack_int ldab, 
										  float anorm, float* rcond);
						   
   Fptr_NL_LAPACKE_spbcon SPBCON;
   float diff;
   void *hModule, *dModule;

	
   pbcon_float_parameters spbcon_obj(LAPACK_ROW_MAJOR, 
                             'L', 512, 500, 400,
                             1450, 1.0 );
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   SPBCON = (Fptr_NL_LAPACKE_spbcon)dlsym(hModule, "LAPACKE_spbcon");
   if (NULL == SPBCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dpbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   spbcon_obj.info    = LAPACKE_spbcon( spbcon_obj.matrix_layout, spbcon_obj.uplo, 
                                        spbcon_obj.n, spbcon_obj.kd,
                                        spbcon_obj.ab, spbcon_obj.ldab, 
										spbcon_obj.anorm, &spbcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   spbcon_obj.inforef = SPBCON( spbcon_obj.matrix_layout, spbcon_obj.uplo, 
                                spbcon_obj.n, spbcon_obj.kd,
								spbcon_obj.abref, spbcon_obj.ldab, 
								spbcon_obj.anorm, &spbcon_obj.rcondref);
	if( spbcon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dpbcon is wrong\n", 
		         spbcon_obj.info );
	}
	if( spbcon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dpbcon is wrong\n", 
		spbcon_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("spbcon_obj.rcond: %lf   spbcon_obj.rcondref: %lf \n", spbcon_obj.rcond, spbcon_obj.rcondref);
#endif
	EXPECT_NEAR(0.0, abs(spbcon_obj.rcond - spbcon_obj.rcondref), DOUBLE_DIFF_THRESHOLD);
}


class pbcon_scomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   
   char uplo;
   int 	n;
   int 	kd, ku;
   lapack_complex_float  *ab, *abref;
   int 	ldab;
   float anorm;
   float rcond, rcondref;   
   int info, inforef;
   
   
   public: 
      pbcon_scomplex_parameters (int matrix_layout, char uplo, int n,
                            int kd, int ku, int ldab, float anorm );
      ~pbcon_scomplex_parameters ();

}; /* end of pbcon_float_parameters  class definition */

/* Destructor definition **/
pbcon_scomplex_parameters:: ~pbcon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbcon_float_parameters object: destructor invoked. \n");
#endif
   pbcon_free();
}

pbcon_scomplex_parameters::pbcon_scomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i, int kd_i, int ku_i,
                             int ldab_i, float anorm_i ){
   int j;
   
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   kd = kd_i;
   ku = ku_i;
   ldab = ldab_i;
   anorm = anorm_i;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ab,  &abref,  (ldab*n));
	
	if( (ab==NULL) || (abref==NULL) ){
       pbcon_free();
       printf(" pbcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ab, abref, (ldab*n));

	info = 0;
	inforef = 0;
}

TEST(pbcon, cpbcon1) {

   /* LAPACKE CPBCON prototype */
   typedef int (*Fptr_NL_LAPACKE_cpbcon) ( int matrix_layout, char uplo, 
             lapack_int n, lapack_int kd, const lapack_complex_float* ab, 
			 lapack_int ldab, float anorm, float* rcond );
						   
   Fptr_NL_LAPACKE_cpbcon CPBCON;
   float diff;
   void *hModule, *dModule;

	
   pbcon_scomplex_parameters cpbcon_obj(LAPACK_ROW_MAJOR, 
                             'U', 512, 500, 400,
                             1450, 1.0 );
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   CPBCON = (Fptr_NL_LAPACKE_cpbcon)dlsym(hModule, "LAPACKE_cpbcon");
   if (NULL == CPBCON)
   {
   	  printf("Could not get the symbol -LAPACKE_cpbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */  
   cpbcon_obj.info    = LAPACKE_cpbcon( cpbcon_obj.matrix_layout, cpbcon_obj.uplo, 
                                        cpbcon_obj.n, cpbcon_obj.kd,
                                        cpbcon_obj.ab, cpbcon_obj.ldab, 
					                    cpbcon_obj.anorm, &cpbcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   cpbcon_obj.inforef = CPBCON( cpbcon_obj.matrix_layout, cpbcon_obj.uplo, 
                                cpbcon_obj.n, cpbcon_obj.kd,
				                cpbcon_obj.abref, cpbcon_obj.ldab, 
			                    cpbcon_obj.anorm, &cpbcon_obj.rcondref);

	if( cpbcon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_cpbcon is wrong\n", 
		         cpbcon_obj.info );
	}
	if( cpbcon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_cpbcon is wrong\n", 
		cpbcon_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("cpbcon_obj.rcond: %lf   cpbcon_obj.rcondref: %lf \n", cpbcon_obj.rcond, cpbcon_obj.rcondref);
#endif
	EXPECT_NEAR(0.0, abs(cpbcon_obj.rcond - cpbcon_obj.rcondref), DOUBLE_DIFF_THRESHOLD);

}

class pbcon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   
   char uplo;
   int 	n;
   int 	kd, ku;
   lapack_complex_double  *ab, *abref;
   int 	ldab;
   double anorm;
   double rcond, rcondref;   
   int info, inforef;
   
   
   public: 
      pbcon_dcomplex_parameters (int matrix_layout, char uplo, int n,
                            int kd, int ku, int ldab, double anorm );
      ~pbcon_dcomplex_parameters ();

}; /* end of pbcon_float_parameters  class definition */

/* Destructor definition **/
pbcon_dcomplex_parameters:: ~pbcon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbcon_float_parameters object: destructor invoked. \n");
#endif
   pbcon_free();
}

pbcon_dcomplex_parameters::pbcon_dcomplex_parameters (int matrix_layout_i, 
                             char uplo_i, int n_i, int kd_i, int ku_i,
                             int ldab_i, double anorm_i ){
   int j;
   
   matrix_layout = matrix_layout_i;
   uplo = uplo_i;
   n = n_i;
   kd = kd_i;
   ku = ku_i;
   ldab = ldab_i;
   anorm = anorm_i;

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ab,  &abref,  (ldab*n));
	
	if( (ab==NULL) || (abref==NULL) ){
       pbcon_free();
       printf(" pbcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ab, abref, (ldab*n));
	
	info = 0;
	inforef = 0;
}

TEST(pbcon, zpbcon1) {

   /* LAPACKE ZPBCON prototype */
   typedef int (*Fptr_NL_LAPACKE_zpbcon) ( int matrix_layout, char uplo, 
           lapack_int n, lapack_int kd, const lapack_complex_double* ab, 
           lapack_int ldab, double anorm, double* rcond );
						   
   Fptr_NL_LAPACKE_zpbcon ZPBCON;
   double diff;
   void *hModule, *dModule;

	
   pbcon_dcomplex_parameters zpbcon_obj(LAPACK_COL_MAJOR, 
                             'L', 512, 500, 400,
                             1450, 1.0 );
   
   dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
   
   hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
   
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   
   ZPBCON = (Fptr_NL_LAPACKE_zpbcon)dlsym(hModule, "LAPACKE_zpbcon");
   if (NULL == ZPBCON)
   {
   	  printf("Could not get the symbol -LAPACKE_zpbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */  
   zpbcon_obj.info    = LAPACKE_zpbcon( zpbcon_obj.matrix_layout, zpbcon_obj.uplo, 
                                        zpbcon_obj.n, zpbcon_obj.kd,
                                        zpbcon_obj.ab, zpbcon_obj.ldab, 
					                    zpbcon_obj.anorm, &zpbcon_obj.rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */	
   zpbcon_obj.inforef = ZPBCON( zpbcon_obj.matrix_layout, zpbcon_obj.uplo, 
                                zpbcon_obj.n, zpbcon_obj.kd,
				                zpbcon_obj.abref, zpbcon_obj.ldab, 
			                    zpbcon_obj.anorm, &zpbcon_obj.rcondref);

	if( zpbcon_obj.info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_cpbcon is wrong\n", 
		         zpbcon_obj.info );
	}
	if( zpbcon_obj.inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_cpbcon is wrong\n", 
		zpbcon_obj.inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("zpbcon_obj.rcond: %lf   zpbcon_obj.rcondref: %lf \n", zpbcon_obj.rcond, zpbcon_obj.rcondref);
#endif
	EXPECT_NEAR(0.0, abs(zpbcon_obj.rcond - zpbcon_obj.rcondref), DOUBLE_DIFF_THRESHOLD);

}


