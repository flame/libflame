#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"
#include <math.h>
#include "lapacke.h"

#define tpcon_free() \
       free (a   ); \
       free (aref);

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

class tpcon_double_parameters{
   public:
   /* input params to the API **/
   int  matrix_layout;
   char norm;
   char uplo,diag;
   int  n;
   double  *a,*aref;
   double anorm;
   double rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      tpcon_double_parameters (int matrix_layout,char norm,char uplo,char diag,int n);
      ~tpcon_double_parameters ();

}; /* end of tpcon_double_parameters  class definition */

/* Destructor definition **/
tpcon_double_parameters:: ~tpcon_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tpcon_double_parameters object: destructor invoked. \n");
#endif
   tpcon_free();
}

tpcon_double_parameters::tpcon_double_parameters (int matrix_layout_i,
                             char norm_i,char uplo_i,char diag_i,
                             int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*(n+1)/2));
    if( (a==NULL) || (aref==NULL) ){
       tpcon_free();
       printf(" tpcon_double_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a,aref,(n*(n+1)/2));
    info = 0;
    inforef = 0;
}

/*  Test fixture class definition */
class dtpcon_test  : public  ::testing::Test {
public:
   tpcon_double_parameters  *dtpcon_obj;
   void SetUp();
   void TearDown () { delete dtpcon_obj; }
};

void dtpcon_test::SetUp(){

   /* LAPACKE DTPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_dtpcon) ( int matrix_layout,char norm,
                    char uplo,char diag,lapack_int n,const double* a,
                     double* rcond  );
                          
   Fptr_NL_LAPACKE_dtpcon DTPCON;
   double diff;
   void *hModule,*dModule;

   dtpcon_obj = new tpcon_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n );

   dtpcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
   DTPCON = (Fptr_NL_LAPACKE_dtpcon)dlsym(hModule,"LAPACKE_dtpcon");
   if (NULL == DTPCON)
   {
      printf("Could not get the symbol -LAPACKE_dtpcon- . Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   dtpcon_obj->info    = LAPACKE_dtpcon( dtpcon_obj->matrix_layout,dtpcon_obj->norm,
                             dtpcon_obj->uplo,dtpcon_obj->diag,dtpcon_obj->n,
                             (const double*)dtpcon_obj->a,
                              &dtpcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   dtpcon_obj->inforef = DTPCON( dtpcon_obj->matrix_layout,dtpcon_obj->norm,
                                dtpcon_obj->uplo,dtpcon_obj->diag,
                                dtpcon_obj->n,(const double*)dtpcon_obj->aref,
                                &dtpcon_obj->rcondref);

    if( dtpcon_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with LAPACKE_dtpcon is wrong\n",
                 dtpcon_obj->info );
    }
    if( dtpcon_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtpcon is wrong\n",
        dtpcon_obj->inforef );
    }

    /* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("dtpcon_obj->rcond: %lf   dtpcon_obj->rcondref: %lf \n",
                           dtpcon_obj->rcond,dtpcon_obj->rcondref);
    printf("diff: %lf  \n",diff);
#endif
    EXPECT_NEAR(0.0,abs(dtpcon_obj->rcond - dtpcon_obj->rcondref),dtpcon_obj->threshold);
}

TEST(dtpcon_test, dtpcon_1) {}
TEST(dtpcon_test, dtpcon_2) {}
TEST(dtpcon_test, dtpcon_3) {}
TEST(dtpcon_test, dtpcon_4) {}
class tpcon_float_parameters{
   public:
   /* input params to the API **/
   int  matrix_layout;
   char norm;
   char uplo,diag;
   int  n;
   float  *a,*aref;
   float anorm;
   float rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      tpcon_float_parameters (int matrix_layout,char norm,char uplo,char diag,int n);
      ~tpcon_float_parameters ();

}; /* end of tpcon_float_parameters  class definition */

/* Destructor definition **/
tpcon_float_parameters:: ~tpcon_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tpcon_float_parameters object: destructor invoked. \n");
#endif
   tpcon_free();
}

tpcon_float_parameters::tpcon_float_parameters (int matrix_layout_i,
                             char norm_i,char uplo_i,char diag_i,
                             int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*(n+1)/2));
    if( (a==NULL) || (aref==NULL) ){
       tpcon_free();
       printf(" tpcon_float_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a,aref,(n*(n+1)/2));
    info = 0;
    inforef = 0;
}

/*  Test fixture class definition */
class stpcon_test  : public  ::testing::Test {
public:
   tpcon_float_parameters  *stpcon_obj;
   void SetUp();
   void TearDown () { delete stpcon_obj; }
};

void stpcon_test::SetUp(){

   /* LAPACKE STPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_stpcon) ( int matrix_layout,char norm,
                    char uplo,char diag,lapack_int n,const float* a,
                     float* rcond  );
                          
   Fptr_NL_LAPACKE_stpcon STPCON;
   float diff;
   void *hModule,*dModule;

   stpcon_obj = new tpcon_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n );

   stpcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
   STPCON = (Fptr_NL_LAPACKE_stpcon)dlsym(hModule,"LAPACKE_stpcon");
   if (NULL == STPCON)
   {
      printf("Could not get the symbol -LAPACKE_stpcon- . Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   stpcon_obj->info    = LAPACKE_stpcon( stpcon_obj->matrix_layout,stpcon_obj->norm,
                             stpcon_obj->uplo,stpcon_obj->diag,stpcon_obj->n,
                             (const float*)stpcon_obj->a,
                              &stpcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   stpcon_obj->inforef = STPCON( stpcon_obj->matrix_layout,stpcon_obj->norm,
                                stpcon_obj->uplo,stpcon_obj->diag,
                                stpcon_obj->n,(const float*)stpcon_obj->aref,
                                &stpcon_obj->rcondref);

    if( stpcon_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with LAPACKE_stpcon is wrong\n",
                 stpcon_obj->info );
    }
    if( stpcon_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stpcon is wrong\n",
        stpcon_obj->inforef );
    }

    /* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("stpcon_obj->rcond: %lf   stpcon_obj->rcondref: %lf \n",
                           stpcon_obj->rcond,stpcon_obj->rcondref);
    printf("diff: %lf  \n",diff);
#endif
    EXPECT_NEAR(0.0,abs(stpcon_obj->rcond - stpcon_obj->rcondref),stpcon_obj->threshold);
}

TEST(stpcon_test, stpcon_1) {}
TEST(stpcon_test, stpcon_2) {}
TEST(stpcon_test, stpcon_3) {}
TEST(stpcon_test, stpcon_4) {}

class tpcon_scomplex_parameters{
   public:
   /* input params to the API **/
   int  matrix_layout;
   char norm;
   char uplo,diag;
   int  n;
   lapack_complex_float  *a,*aref;
   float anorm;
   float rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      tpcon_scomplex_parameters (int matrix_layout,char norm,char uplo,char diag,int n);
      ~tpcon_scomplex_parameters ();

}; /* end of tpcon_scomplex_parameters  class definition */

/* Destructor definition **/
tpcon_scomplex_parameters:: ~tpcon_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tpcon_scomplex_parameters object: destructor invoked. \n");
#endif
   tpcon_free();
}

tpcon_scomplex_parameters::tpcon_scomplex_parameters (int matrix_layout_i,
                             char norm_i,char uplo_i,char diag_i,
                             int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   anorm = (float)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*(n+1)/2));
    if( (a==NULL) || (aref==NULL) ){
       tpcon_free();
       printf(" tpcon_scomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a,aref,(n*(n+1)/2));
    info = 0;
    inforef = 0;
}

/*  Test fixture class definition */
class ctpcon_test  : public  ::testing::Test {
public:
   tpcon_scomplex_parameters  *ctpcon_obj;
   void SetUp();
   void TearDown () { delete ctpcon_obj; }
};

void ctpcon_test::SetUp(){

   /* LAPACKE CTPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_ctpcon) ( int matrix_layout,char norm,
                    char uplo,char diag,lapack_int n, const lapack_complex_float* a,
                     float* rcond  );
                          
   Fptr_NL_LAPACKE_ctpcon CTPCON;
   float diff;
   void *hModule,*dModule;

   ctpcon_obj = new tpcon_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n );

   ctpcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
   CTPCON = (Fptr_NL_LAPACKE_ctpcon)dlsym(hModule,"LAPACKE_ctpcon");
   if (NULL == CTPCON)
   {
      printf("Could not get the symbol -LAPACKE_ctpcon- . Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   ctpcon_obj->info    = LAPACKE_ctpcon( ctpcon_obj->matrix_layout,ctpcon_obj->norm,
                             ctpcon_obj->uplo,ctpcon_obj->diag,ctpcon_obj->n,
                             (const lapack_complex_float*)ctpcon_obj->a,
                              &ctpcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   ctpcon_obj->inforef = CTPCON( ctpcon_obj->matrix_layout,ctpcon_obj->norm,
                                ctpcon_obj->uplo,ctpcon_obj->diag,
                                ctpcon_obj->n,(const lapack_complex_float*)ctpcon_obj->aref,
                                &ctpcon_obj->rcondref);

    if( ctpcon_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with LAPACKE_ctpcon is wrong\n",
                 ctpcon_obj->info );
    }
    if( ctpcon_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctpcon is wrong\n",
        ctpcon_obj->inforef );
    }

    /* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("ctpcon_obj->rcond: %lf   ctpcon_obj->rcondref: %lf \n",
                           ctpcon_obj->rcond,ctpcon_obj->rcondref);
    printf("diff: %lf  \n",diff);
#endif
    EXPECT_NEAR(0.0,abs(ctpcon_obj->rcond - ctpcon_obj->rcondref),ctpcon_obj->threshold);
}

TEST(ctpcon_test, ctpcon_1) {}
TEST(ctpcon_test, ctpcon_2) {}
TEST(ctpcon_test, ctpcon_3) {}
TEST(ctpcon_test, ctpcon_4) {}

class tpcon_dcomplex_parameters{
   public:
   /* input params to the API **/
   int  matrix_layout;
   char norm;
   char uplo,diag;
   int  n;
   lapack_complex_double  *a,*aref;
   double anorm;
   double rcond,rcondref;   int info,inforef;
   float threshold;
   public: 
      tpcon_dcomplex_parameters (int matrix_layout,char norm,char uplo,char diag,int n);
      ~tpcon_dcomplex_parameters ();

}; /* end of tpcon_dcomplex_parameters  class definition */

/* Destructor definition **/
tpcon_dcomplex_parameters:: ~tpcon_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tpcon_dcomplex_parameters object: destructor invoked. \n");
#endif
   tpcon_free();
}

tpcon_dcomplex_parameters::tpcon_dcomplex_parameters (int matrix_layout_i,
                             char norm_i,char uplo_i,char diag_i,
                             int n_i ){
   int j;
   matrix_layout = matrix_layout_i;
   norm = norm_i;
   uplo = uplo_i;
   diag = diag_i;
   n = n_i;
   anorm = (double)n/2.0; //TODO: replace with L1 norm calculation

/* Memory allocation of the buffers */
   lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*(n+1)/2));
    if( (a==NULL) || (aref==NULL) ){
       tpcon_free();
       printf(" tpcon_dcomplex_parameters object: malloc error. Exiting...\n");
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a,aref,(n*(n+1)/2));
    info = 0;
    inforef = 0;
}

/*  Test fixture class definition */
class ztpcon_test  : public  ::testing::Test {
public:
   tpcon_dcomplex_parameters  *ztpcon_obj;
   void SetUp();
   void TearDown () { delete ztpcon_obj; }
};

void ztpcon_test::SetUp(){

   /* LAPACKE ZTPCON prototype */
   typedef int (*Fptr_NL_LAPACKE_ztpcon) ( int matrix_layout,char norm,
                    char uplo,char diag,lapack_int n, const lapack_complex_double* a,
                     double* rcond  );
                          
   Fptr_NL_LAPACKE_ztpcon ZTPCON;
   double diff;
   void *hModule,*dModule;

   ztpcon_obj = new tpcon_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].norm_gbcon,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
						   lin_solver_paramslist[idx].n );

   ztpcon_obj->threshold = lin_solver_paramslist[idx].solver_threhold;
   idx = Circular_Increment_Index(idx);

   dModule = dlopen(NETLIB_BLAS_LIB,RTLD_NOW | RTLD_GLOBAL);
   hModule = dlopen(NETLIB_LAPACKE_LIB,RTLD_NOW);
   if ((NULL == hModule) || (NULL == dModule) )
   {
       printf("Load Library failed. Exiting ....\n");
       exit (0);
   }
   ZTPCON = (Fptr_NL_LAPACKE_ztpcon)dlsym(hModule,"LAPACKE_ztpcon");
   if (NULL == ZTPCON)
   {
      printf("Could not get the symbol -LAPACKE_ztpcon- . Exiting...\n");
      dlclose(hModule);
      dlclose(dModule);
      exit (-1);
   }
   /* Compute libflame's Lapacke o/p  */
   ztpcon_obj->info    = LAPACKE_ztpcon( ztpcon_obj->matrix_layout,ztpcon_obj->norm,
                             ztpcon_obj->uplo,ztpcon_obj->diag,ztpcon_obj->n,
                             (const lapack_complex_double*)ztpcon_obj->a,
                              &ztpcon_obj->rcond);

   /* Compute the reference o/p by invoking Netlib-Lapack's API */   
   ztpcon_obj->inforef = ZTPCON( ztpcon_obj->matrix_layout,ztpcon_obj->norm,
                                ztpcon_obj->uplo,ztpcon_obj->diag,
                                ztpcon_obj->n,(const lapack_complex_double*)ztpcon_obj->aref,
                                &ztpcon_obj->rcondref);

    if( ztpcon_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with LAPACKE_ztpcon is wrong\n",
                 ztpcon_obj->info );
    }
    if( ztpcon_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztpcon is wrong\n",
        ztpcon_obj->inforef );
    }

    /* Compute Difference in C and CPP buffer */
#if LAPACKE_TEST_VERBOSE
    printf("ztpcon_obj->rcond: %lf   ztpcon_obj->rcondref: %lf \n",
                           ztpcon_obj->rcond,ztpcon_obj->rcondref);
    printf("diff: %lf  \n",diff);
#endif
    EXPECT_NEAR(0.0,abs(ztpcon_obj->rcond - ztpcon_obj->rcondref),ztpcon_obj->threshold);
}

TEST(ztpcon_test, ztpcon_1) {}
TEST(ztpcon_test, ztpcon_2) {}
TEST(ztpcon_test, ztpcon_3) {}
TEST(ztpcon_test, ztpcon_4) {}
