#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#include <math.h>
#include "lapacke.h"

#define trcon_free() \
       free (a   ); \
       free (aref);

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

class trcon_double_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   char uplo,diag;
   int 	n;
   double  *a,*aref;
   int 	lda;
   double anorm;
   double rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      trcon_double_parameters (int matrix_layout,char norm,char uplo,char diag,int n);
      ~trcon_double_parameters ();

}; /* end of trcon_double_parameters  class definition */

/* Destructor definition **/
trcon_double_parameters:: ~trcon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trcon_double_parameters object: destructor invoked. \n");
#endif
   trcon_free();
}

trcon_double_parameters::trcon_double_parameters (int matrix_layout_i,
                             char norm_i,char uplo_i,char diag_i,
                             int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   lda = n_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (lda*n));
	if( (a==NULL) || (aref==NULL) ){
       trcon_free();
       printf(" trcon_double_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,(lda*n));
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class dtrcon_test  : public  ::testing::Test {
public:
   trcon_double_parameters  *dtrcon_obj;
   void SetUp();
   void TearDown () { delete dtrcon_obj; }
};

void dtrcon_test::SetUp(){

   /* LAPACKE DTRCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dtrcon) ( int matrix_layout,char norm,
                    char uplo,char diag,lapack_int n,const double* a,
					lapack_int lda,double* rcond  );
						  
   Fptr_NL_LAPACKE_dtrcon DTRCON;
   double diff;
   void *hModule,*dModule;

   dtrcon_obj = new trcon_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n );

   dtrcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);

   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   DTRCON = (Fptr_NL_LAPACKE_dtrcon)dlsym(hModule,"LAPACKE_dtrcon");
   if (NULL == DTRCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dtrcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   dtrcon_obj->info    = LAPACKE_dtrcon( dtrcon_obj->matrix_layout,dtrcon_obj->norm,
                             dtrcon_obj->uplo,dtrcon_obj->diag,dtrcon_obj->n,
							 (const double*)dtrcon_obj->a,dtrcon_obj->lda,
							  &dtrcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   dtrcon_obj->inforef = DTRCON( dtrcon_obj->matrix_layout,dtrcon_obj->norm,
                                dtrcon_obj->uplo,dtrcon_obj->diag,
                                dtrcon_obj->n,(const double*)dtrcon_obj->aref,
								dtrcon_obj->lda, &dtrcon_obj->rcondref);

 	if( dtrcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_dtrcon is wrong\n",
		         dtrcon_obj->info );
	}
	if( dtrcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_dtrcon is wrong\n",
		dtrcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_d( dtrcon_obj->lda*dtrcon_obj->n,dtrcon_obj->a,dtrcon_obj->aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	/*
	FLA_Nrm2_external( dtrcon_obj->A,dtrcon_obj->norm );
    FLA_Nrm2_external( dtrcon_obj->Aref,dtrcon_obj-> normref );
	diff = FLA_abs( dtrcon_obj-> norm, dtrcon_obj-> normref);
	*/
#if LAPACKE_TEST_VERBOSE
    printf("dtrcon_obj->rcond: %lf   dtrcon_obj->rcondref: %lf \n",
	                       dtrcon_obj->rcond,dtrcon_obj->rcondref);
	printf("diff: %lf  \n",diff);
#endif
	EXPECT_NEAR(0.0,abs(dtrcon_obj->rcond - dtrcon_obj->rcondref),dtrcon_obj->threshold);
}

TEST(dtrcon_test, dtrcon_1) {}
TEST(dtrcon_test, dtrcon_2) {}
TEST(dtrcon_test, dtrcon_3) {}
TEST(dtrcon_test, dtrcon_4) {}
class trcon_float_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   char uplo,diag;
   int 	n;
   float  *a,*aref;
   int 	lda;
   float anorm;
   float rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      trcon_float_parameters (int matrix_layout,char norm,char uplo,char diag,int n);
      ~trcon_float_parameters ();

}; /* end of trcon_float_parameters  class definition */

/* Destructor definition **/
trcon_float_parameters:: ~trcon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trcon_float_parameters object: destructor invoked. \n");
#endif
   trcon_free();
}

trcon_float_parameters::trcon_float_parameters (int matrix_layout_i,
                             char norm_i, char uplo_i,char diag_i,
                             int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   lda = n_i;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (lda*n));
	if( (a==NULL) || (aref==NULL)){
       trcon_free();
       printf(" trcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,(lda*n));
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class strcon_test  : public  ::testing::Test {
public:
   trcon_float_parameters  *strcon_obj;
   void SetUp();
   void TearDown () { delete strcon_obj; }
};

void strcon_test::SetUp(){

   /* LAPACKE STRCON prototype */
   typedef int (*Fptr_NL_LAPACKE_strcon) ( int matrix_layout,char norm,
                     char uplo,char diag,lapack_int n,const float* a,
					 lapack_int lda,float* rcond );

   Fptr_NL_LAPACKE_strcon STRCON;
   float diff;
   void *hModule,*dModule;

   strcon_obj = new trcon_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n );

   strcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   STRCON = (Fptr_NL_LAPACKE_strcon)dlsym(hModule,"LAPACKE_strcon");
   if (NULL == STRCON)
   {
   	  printf("Could not get the symbol -LAPACKE_dtrcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute the reference o/p by invoking Netlib-Lapack's API */
   strcon_obj->inforef = STRCON( strcon_obj->matrix_layout,strcon_obj->norm,
                                strcon_obj->uplo,strcon_obj->diag,                               strcon_obj->n,
				strcon_obj->aref,strcon_obj->lda,
				&strcon_obj->rcondref);

   /* Compute libflame's Lapacke o/p  */
   strcon_obj->info    = LAPACKE_strcon( strcon_obj->matrix_layout,strcon_obj->norm,
                                        strcon_obj->uplo,strcon_obj->diag,
                                        strcon_obj->n,
                                        strcon_obj->a,strcon_obj->lda,
				        &strcon_obj->rcond);

	if( strcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with libflame LAPACKE_strcon is wrong\n",
		         strcon_obj->info );
	}
	if( strcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_strcon is wrong\n",
		strcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_s( strcon_obj->lda*strcon_obj->n,strcon_obj->a,strcon_obj->aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	/*
	FLA_Nrm2_external( strcon_obj->A,strcon_obj->norm );
    FLA_Nrm2_external( strcon_obj->Aref,strcon_obj-> normref );
	diff = FLA_abs( strcon_obj-> norm, strcon_obj-> normref);
	*/
#if LAPACKE_TEST_VERBOSE
    printf("strcon_obj->rcond: %lf   strcon_obj->rcondref: %lf \n",
	                       strcon_obj->rcond,strcon_obj->rcondref);
	printf("diff: %f  \n",diff);
#endif
   EXPECT_NEAR(0.0,abs(strcon_obj->rcond - strcon_obj->rcondref),DOUBLE_DIFF_THRESHOLD);
}

TEST(strcon_test, strcon_1) {}
TEST(strcon_test, strcon_2) {}
TEST(strcon_test, strcon_3) {}
TEST(strcon_test, strcon_4) {}

class trcon_scomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   char uplo,diag;
   int 	n;
   lapack_complex_float  *a,*aref;
   int 	lda;
   float anorm;
   float rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      trcon_scomplex_parameters (int matrix_layout,char norm,char uplo,char diag,int n);
      ~trcon_scomplex_parameters ();

}; /* end of trcon_float_parameters  class definition */

/* Destructor definition **/
trcon_scomplex_parameters:: ~trcon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trcon_float_parameters object: destructor invoked. \n");
#endif
   trcon_free();
}

trcon_scomplex_parameters::trcon_scomplex_parameters (int matrix_layout_i,
                             char norm_i, char uplo_i,char diag_i,
                             int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   lda = n_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation
/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (lda*n));
	if( (a==NULL) || (aref==NULL)){
       trcon_free();
       printf(" trcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,(lda*n));
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class ctrcon_test  : public  ::testing::Test {
public:
   trcon_scomplex_parameters  *ctrcon_obj;
   void SetUp();
   void TearDown () { delete ctrcon_obj; }
};

void ctrcon_test::SetUp(){

   /* LAPACKE CTRCON prototype */
   typedef int (*Fptr_NL_LAPACKE_ctrcon) ( int matrix_layout,char norm,
        char uplo,char diag,lapack_int n,const lapack_complex_float* a,
		lapack_int lda,float* rcond );

   Fptr_NL_LAPACKE_ctrcon CTRCON;
   float diff;
   void *hModule,*dModule;

   ctrcon_obj = new trcon_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n );

   ctrcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   CTRCON = (Fptr_NL_LAPACKE_ctrcon)dlsym(hModule,"LAPACKE_ctrcon");
   if (NULL == CTRCON)
   {
   	  printf("Could not get the symbol -LAPACKE_ctrcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */  
   ctrcon_obj->info    = LAPACKE_ctrcon( ctrcon_obj->matrix_layout,ctrcon_obj->norm,
                                        ctrcon_obj->uplo,ctrcon_obj->diag,
                                        ctrcon_obj->n,
                                        ctrcon_obj->a,ctrcon_obj->lda,
			                            &ctrcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   ctrcon_obj->inforef = CTRCON( ctrcon_obj->matrix_layout,ctrcon_obj->norm,
                                ctrcon_obj->uplo,ctrcon_obj->diag,
                                ctrcon_obj->n,
				                ctrcon_obj->aref,ctrcon_obj->lda,
				                &ctrcon_obj->rcondref);

	if( ctrcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_ctrcon is wrong\n",
		         ctrcon_obj->info );
	}
	if( ctrcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_ctrcon is wrong\n",
		ctrcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_c( ctrcon_obj->lda*ctrcon_obj->n,ctrcon_obj->a,ctrcon_obj->aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	/*
	FLA_Nrm2_external( ctrcon_obj->A,ctrcon_obj->norm );
    FLA_Nrm2_external( ctrcon_obj->Aref,ctrcon_obj-> normref );
	diff = FLA_abs( ctrcon_obj-> norm, ctrcon_obj-> normref);
	*/
#if LAPACKE_TEST_VERBOSE
    printf("ctrcon_obj->rcond: %lf   ctrcon_obj->rcondref: %lf \n",
	                       ctrcon_obj->rcond,ctrcon_obj->rcondref);
	printf("diff: %f  \n",diff);
#endif
	EXPECT_NEAR(0.0,abs(ctrcon_obj->rcond - ctrcon_obj->rcondref),
	                                      DOUBLE_DIFF_THRESHOLD);

}

TEST(ctrcon_test, ctrcon_1) {}
TEST(ctrcon_test, ctrcon_2) {}
TEST(ctrcon_test, ctrcon_3) {}
TEST(ctrcon_test, ctrcon_4) {}

class trcon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int 	matrix_layout;
   char norm;
   char uplo,diag;
   int 	n;
   lapack_complex_double  *a,*aref;
   int 	lda;
   double anorm;
   double rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      trcon_dcomplex_parameters (int matrix_layout,char norm,char uplo,char diag,int n);
      ~trcon_dcomplex_parameters ();

}; /* end of trcon_float_parameters  class definition */

/* Destructor definition **/
trcon_dcomplex_parameters:: ~trcon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trcon_float_parameters object: destructor invoked. \n");
#endif
   trcon_free();
}

trcon_dcomplex_parameters::trcon_dcomplex_parameters (int matrix_layout_i,
                             char norm_i,char uplo_i,char diag_i,
                             int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   lda = n_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation
/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (lda*n));
	if( (a==NULL) || (aref==NULL) ){
       trcon_free();
       printf(" trcon_float_parameters object: malloc error. Exiting...\n");
	   exit(-1);
	}

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,(lda*n));
	info = 0;
	inforef = 0;
}

/*  Test fixture class definition */
class ztrcon_test  : public  ::testing::Test {
public:
   trcon_dcomplex_parameters  *ztrcon_obj;
   void SetUp();
   void TearDown () { delete ztrcon_obj; }
};

void ztrcon_test::SetUp(){

   /* LAPACKE ZTRCON prototype */
   typedef int (*Fptr_NL_LAPACKE_ztrcon) ( int matrix_layout,char norm,
      char uplo,char diag,lapack_int n,const lapack_complex_double* a,
	  lapack_int lda,double* rcond );

   Fptr_NL_LAPACKE_ztrcon ZTRCON;
   double diff;
   void *hModule,*dModule;

   ztrcon_obj = new trcon_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n );

   ztrcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
   	   printf("Load Library failed. Exiting ....\n");
   	   exit (0);
   }
   ZTRCON = (Fptr_NL_LAPACKE_ztrcon)dlsym(hModule,"LAPACKE_ztrcon");
   if (NULL == ZTRCON)
   {
   	  printf("Could not get the symbol -LAPACKE_ztrcon- . Exiting...\n");
   	  dlclose(hModule);
   	  dlclose(dModule);
   	  exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */  
   ztrcon_obj->info    = LAPACKE_ztrcon( ztrcon_obj->matrix_layout,ztrcon_obj->norm,
                                        ztrcon_obj->uplo,ztrcon_obj->diag,
                                        ztrcon_obj->n,
                                        ztrcon_obj->a,ztrcon_obj->lda,
			                            &ztrcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   ztrcon_obj->inforef = ZTRCON( ztrcon_obj->matrix_layout,ztrcon_obj->norm,
                                ztrcon_obj->uplo,ztrcon_obj->diag,
                                ztrcon_obj->n,
				                ztrcon_obj->aref,ztrcon_obj->lda,
				                &ztrcon_obj->rcondref);

	if( ztrcon_obj->info < 0 ) {
		printf( "\n warning: The i:%d th argument with LAPACKE_ztrcon is wrong\n",
		         ztrcon_obj->info );
	}
	if( ztrcon_obj->inforef < 0 ) {
		printf( "The i:%d th argument with Netlib LAPACKE_ztrcon is wrong\n",
		ztrcon_obj->inforef );
	}

	/* Compute Difference in C and CPP buffer */
	diff =  computeDiff_z( ztrcon_obj->lda*ztrcon_obj->n,ztrcon_obj->a,ztrcon_obj->aref );
	if ( diff > DOUBLE_DIFF_THRESHOLD)
	{
	   printf("\n Warning: factored band matrix AB differences exist.");
	}
	/*
	FLA_Nrm2_external( ztrcon_obj->A,ztrcon_obj->norm );
    FLA_Nrm2_external( ztrcon_obj->Aref,ztrcon_obj-> normref );
	diff = FLA_abs( ztrcon_obj-> norm, ztrcon_obj-> normref);
	*/
#if LAPACKE_TEST_VERBOSE
    printf("ztrcon_obj->rcond: %lf   ztrcon_obj->rcondref: %lf \n",ztrcon_obj->rcond,ztrcon_obj->rcondref);
	printf("diff: %lf  \n",diff);
#endif
	EXPECT_NEAR(0.0,abs(ztrcon_obj->rcond - ztrcon_obj->rcondref),DOUBLE_DIFF_THRESHOLD);

}

TEST(ztrcon_test, ztrcon_1) {}
TEST(ztrcon_test, ztrcon_2) {}
TEST(ztrcon_test, ztrcon_3) {}
TEST(ztrcon_test, ztrcon_4) {}
