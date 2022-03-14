#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"
#include <math.h>
#include "lapacke.h"

#define tbcon_free() \
       free (ab   ); \
       free (abref);

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

class tbcon_double_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   char uplo,diag;
   int 	n;
   int kd;
   double  *ab,*abref;
   int 	ldab;
   double anorm;
   double rcond,rcondref;   
   int info,inforef;
   float threshold;
   public: 
      tbcon_double_parameters (int matrix_layout,char norm,char uplo,char diag,int n,
                               int kd,int ldab);
      ~tbcon_double_parameters ();

}; /* end of tbcon_double_parameters  class definition */

/* Destructor definition **/
tbcon_double_parameters:: ~tbcon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbcon_double_parameters object: destructor invoked. \n");
#endif
   tbcon_free();
}

tbcon_double_parameters::tbcon_double_parameters (int matrix_layout_i,
                             char norm_i,char uplo_i,char diag_i,
							 int n_i,int kd_i,int ldab_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   kd = kd_i;
   ldab = ldab_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &ab, &abref, (ldab*n));
	if( (ab==NULL) || (abref==NULL) ){
       tbcon_free();
       printf(" tbcon_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( ab,abref,(ldab*n));
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class dtbcon_test  : public  ::testing::Test {
public:
   tbcon_double_parameters  *dtbcon_obj;
   void SetUp();
   void TearDown () { delete dtbcon_obj; }
};

void dtbcon_test::SetUp(){

   /* LAPACKE DTBCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dtbcon) ( int matrix_layout,char norm,
                      char uplo,char diag,lapack_int n,lapack_int kd,
					  const double* ab,lapack_int ldab,double* rcond );
   Fptr_NL_LAPACKE_dtbcon DTBCON;
   double diff;
   void *hModule,*dModule;

   dtbcon_obj = new tbcon_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].kd,
						   lin_solver_paramslist[idx].ldab_gbcon );

   dtbcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   DTBCON = (Fptr_NL_LAPACKE_dtbcon)dlsym(hModule,"LAPACKE_dtbcon");
   if (NULL == DTBCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dtbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   dtbcon_obj->info    = LAPACKE_dtbcon( dtbcon_obj->matrix_layout,dtbcon_obj->norm,
                             dtbcon_obj->uplo,dtbcon_obj->diag,dtbcon_obj->n,
							 dtbcon_obj->kd,(const double*)dtbcon_obj->ab,
							 dtbcon_obj->ldab,&dtbcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   dtbcon_obj->inforef = DTBCON( dtbcon_obj->matrix_layout,dtbcon_obj->norm,
                                dtbcon_obj->uplo,dtbcon_obj->diag,
                                dtbcon_obj->n,dtbcon_obj->kd,
								(const double*)dtbcon_obj->abref,
								dtbcon_obj->ldab, &dtbcon_obj->rcondref);

 	if( dtbcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dtbcon is wrong\n",
		         dtbcon_obj->info );
	}
	if( dtbcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dtbcon is wrong\n",
		dtbcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("dtbcon_obj->rcond: %lf   dtbcon_obj->rcondref: %lf \n",
	                       dtbcon_obj->rcond,dtbcon_obj->rcondref);
	printf("diff: %lf  \n",diff);
#endif
	EXPECT_NEAR(0.0,abs(dtbcon_obj->rcond - dtbcon_obj->rcondref),DOUBLE_DIFF_THRESHOLD);
}

TEST(dtbcon_test, dtbcon_1) {}
TEST(dtbcon_test, dtbcon_2) {}
TEST(dtbcon_test, dtbcon_3) {}
TEST(dtbcon_test, dtbcon_4) {}
class tbcon_float_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   char uplo,diag;
   int 	n;
   int 	kd;
   float  *ab,*abref;
   int 	ldab;
   float anorm;
   float rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      tbcon_float_parameters (int matrix_layout,char norm,char uplo,
	                          char diag,int n,int kd,
							  int ldab );
      ~tbcon_float_parameters ();

}; /* end of tbcon_float_parameters  class definition */

/* Destructor definition **/
tbcon_float_parameters:: ~tbcon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbcon_float_parameters object: destructor invoked. \n");
#endif
   tbcon_free();
}

tbcon_float_parameters::tbcon_float_parameters (int matrix_layout_i,
                             char norm_i, char uplo_i,char diag_i,
							 int n_i,int kd_i,
                             int ldab_i){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   kd = kd_i;
   ldab = ldab_i;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &ab, &abref, (ldab*n));
	if( (ab==NULL) || (abref==NULL)){
       tbcon_free();
       printf(" tbcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( ab,abref,(ldab*n));
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class stbcon_test  : public  ::testing::Test {
public:
   tbcon_float_parameters  *stbcon_obj;
   void SetUp();
   void TearDown () { delete stbcon_obj; }
};

void stbcon_test::SetUp(){

   /* LAPACKE STBCON prototype */
   typedef int (*Fptr_NL_LAPACKE_stbcon) ( int matrix_layout,char norm,
                      char uplo,char diag,lapack_int n,lapack_int kd,
					  const float* ab,lapack_int ldab,float* rcond);

   Fptr_NL_LAPACKE_stbcon STBCON;
   float diff;
   void *hModule,*dModule;

   stbcon_obj = new tbcon_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].kd,
						   lin_solver_paramslist[idx].ldab_gbcon );

   stbcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   STBCON = (Fptr_NL_LAPACKE_stbcon)dlsym(hModule,"LAPACKE_dtbcon");
   if (NULL == STBCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dtbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   stbcon_obj->info    = LAPACKE_stbcon( stbcon_obj->matrix_layout,stbcon_obj->norm,
                                        stbcon_obj->uplo,stbcon_obj->diag,
                                        stbcon_obj->n,stbcon_obj->kd,
                                        stbcon_obj->ab,stbcon_obj->ldab,
						                &stbcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   stbcon_obj->inforef = STBCON( stbcon_obj->matrix_layout,stbcon_obj->norm,
                                stbcon_obj->uplo,stbcon_obj->diag,
								stbcon_obj->n,stbcon_obj->kd,
								stbcon_obj->abref,stbcon_obj->ldab,
								&stbcon_obj->rcondref);

	if( stbcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with libflame LAPACKE_stbcon is wrong\n",
		         stbcon_obj->info );
	}
	if( stbcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_stbcon is wrong\n",
		stbcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("stbcon_obj->rcond: %lf   stbcon_obj->rcondref: %lf \n",
	                       stbcon_obj->rcond,stbcon_obj->rcondref);
	printf("diff: %f  \n",diff);
#endif
   EXPECT_NEAR(0.0,abs(stbcon_obj->rcond - stbcon_obj->rcondref),DOUBLE_DIFF_THRESHOLD);
}

TEST(stbcon_test, stbcon_1) {}
TEST(stbcon_test, stbcon_2) {}
TEST(stbcon_test, stbcon_3) {}
TEST(stbcon_test, stbcon_4) {}

class tbcon_scomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   char uplo,diag;
   int 	n;
   int 	kd;
   lapack_complex_float  *ab,*abref;
   int 	ldab;
   float anorm;
   float rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      tbcon_scomplex_parameters (int matrix_layout,char norm,char uplo,
	                             char diag,int n,int kd,
								 int ldab);
      ~tbcon_scomplex_parameters ();

}; /* end of tbcon_float_parameters  class definition */

/* Destructor definition **/
tbcon_scomplex_parameters:: ~tbcon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbcon_float_parameters object: destructor invoked. \n");
#endif
   tbcon_free();
}

tbcon_scomplex_parameters::tbcon_scomplex_parameters (int matrix_layout_i,
                             char norm_i, char uplo_i,char diag_i,
							 int n_i,int kd_i,
                             int ldab_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   kd = kd_i;
   ldab = ldab_i;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &ab, &abref, (ldab*n));
	if( (ab==NULL) || (abref==NULL)){
       tbcon_free();
       printf(" tbcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( ab,abref,(ldab*n));
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class ctbcon_test  : public  ::testing::Test {
public:
   tbcon_scomplex_parameters  *ctbcon_obj;
   void SetUp();
   void TearDown () { delete ctbcon_obj; }
};

void ctbcon_test::SetUp(){

   /* LAPACKE CTBCON prototype */
   typedef int (*Fptr_NL_LAPACKE_ctbcon) ( int matrix_layout,char norm,
                      char uplo,char diag,lapack_int n,lapack_int kd,
					  const lapack_complex_float* ab,lapack_int ldab,
					  float* rcond );

   Fptr_NL_LAPACKE_ctbcon CTBCON;
   float diff;
   void *hModule,*dModule;

   ctbcon_obj = new tbcon_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].kd,
						   lin_solver_paramslist[idx].ldab_gbcon );

   ctbcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   CTBCON = (Fptr_NL_LAPACKE_ctbcon)dlsym(hModule,"LAPACKE_ctbcon");
   if (NULL == CTBCON)
   {
   	  printf("Could not get the symbol -LAPACKE_ctbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */  
   ctbcon_obj->info    = LAPACKE_ctbcon( ctbcon_obj->matrix_layout,ctbcon_obj->norm,
                                        ctbcon_obj->uplo,ctbcon_obj->diag,
                                        ctbcon_obj->n,ctbcon_obj->kd,
                                        ctbcon_obj->ab,ctbcon_obj->ldab,
			                            &ctbcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   ctbcon_obj->inforef = CTBCON( ctbcon_obj->matrix_layout,ctbcon_obj->norm,
                                ctbcon_obj->uplo,ctbcon_obj->diag,
                                ctbcon_obj->n,ctbcon_obj->kd,
				                ctbcon_obj->abref,ctbcon_obj->ldab,
				                &ctbcon_obj->rcondref);

	if( ctbcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_ctbcon is wrong\n",
		         ctbcon_obj->info );
	}
	if( ctbcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_ctbcon is wrong\n",
		ctbcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("ctbcon_obj->rcond: %lf   ctbcon_obj->rcondref: %lf \n",
	                       ctbcon_obj->rcond,ctbcon_obj->rcondref);
	printf("diff: %f  \n",diff);
#endif
	EXPECT_NEAR(0.0,abs(ctbcon_obj->rcond - ctbcon_obj->rcondref),
	                                      DOUBLE_DIFF_THRESHOLD);

}

TEST(ctbcon_test, ctbcon_1) {}
TEST(ctbcon_test, ctbcon_2) {}
TEST(ctbcon_test, ctbcon_3) {}
TEST(ctbcon_test, ctbcon_4) {}

class tbcon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   char uplo,diag;
   int 	n;
   int 	kd;
   lapack_complex_double  *ab,*abref;
   int 	ldab;
   double anorm;
   double rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      tbcon_dcomplex_parameters (int matrix_layout,char norm,char uplo,
	                             char diag,int n,int kd,
								 int ldab );
      ~tbcon_dcomplex_parameters ();

}; /* end of tbcon_float_parameters  class definition */

/* Destructor definition **/
tbcon_dcomplex_parameters:: ~tbcon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbcon_float_parameters object: destructor invoked. \n");
#endif
   tbcon_free();
}

tbcon_dcomplex_parameters::tbcon_dcomplex_parameters (int matrix_layout_i,
                             char norm_i,char uplo_i,char diag_i,
							 int n_i,int kd_i,
                             int ldab_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   kd = kd_i;
   ldab = ldab_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &ab, &abref, (ldab*n));
	if( (ab==NULL) || (abref==NULL) ){
       tbcon_free();
       printf(" tbcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( ab,abref,(ldab*n));
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class ztbcon_test  : public  ::testing::Test {
public:
   tbcon_dcomplex_parameters  *ztbcon_obj;
   void SetUp();
   void TearDown () { delete ztbcon_obj; }
};

void ztbcon_test::SetUp(){

   /* LAPACKE ZTBCON prototype */
   typedef int (*Fptr_NL_LAPACKE_ztbcon) ( int matrix_layout,char norm,
                      char uplo,char diag,lapack_int n,lapack_int kd,
					  const lapack_complex_double* ab,lapack_int ldab,
					  double* rcond );

   Fptr_NL_LAPACKE_ztbcon ZTBCON;
   double diff;
   void *hModule,*dModule;

   ztbcon_obj = new tbcon_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].kd,
						   lin_solver_paramslist[idx].ldab_gbcon );

   ztbcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   ZTBCON = (Fptr_NL_LAPACKE_ztbcon)dlsym(hModule,"LAPACKE_ztbcon");
   if (NULL == ZTBCON)
   {
   	  printf("Could not get the symbol -LAPACKE_ztbcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */  
   ztbcon_obj->info    = LAPACKE_ztbcon( ztbcon_obj->matrix_layout,ztbcon_obj->norm,
                                        ztbcon_obj->uplo,ztbcon_obj->diag,
                                        ztbcon_obj->n,ztbcon_obj->kd,
                                        ztbcon_obj->ab,ztbcon_obj->ldab,
			                            &ztbcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   ztbcon_obj->inforef = ZTBCON( ztbcon_obj->matrix_layout,ztbcon_obj->norm,
                                ztbcon_obj->uplo,ztbcon_obj->diag,
                                ztbcon_obj->n,ztbcon_obj->kd,
				                ztbcon_obj->abref,ztbcon_obj->ldab,
				                &ztbcon_obj->rcondref);

	if( ztbcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_ctbcon is wrong\n",
		         ztbcon_obj->info );
	}
	if( ztbcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_ctbcon is wrong\n",
		ztbcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("ztbcon_obj->rcond: %lf   ztbcon_obj->rcondref: %lf \n",ztbcon_obj->rcond,ztbcon_obj->rcondref);
	printf("diff: %lf  \n",diff);
#endif
	EXPECT_NEAR(0.0,abs(ztbcon_obj->rcond - ztbcon_obj->rcondref),DOUBLE_DIFF_THRESHOLD);

}

TEST(ztbcon_test, ztbcon_1) {}
TEST(ztbcon_test, ztbcon_2) {}
TEST(ztbcon_test, ztbcon_3) {}
TEST(ztbcon_test, ztbcon_4) {}
