#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define sytri2_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (ipiv != NULL) free (ipiv); \
  if (ipivref != NULL)free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin sytri2_float_parameters  class definition */
class sytri2_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Input/ Output parameters */
      float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri2_float_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri2_float_parameters (); 
};  /* end of sytri2_float_parameters  class definition */


/* Constructor sytri2_float_parameters definition */
sytri2_float_parameters:: sytri2_float_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytri2 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       sytri2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri2_float_parameters:: ~sytri2_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri2_float_parameters object: destructor invoked. \n");
#endif
   sytri2_free();
}

//  Test fixture class definition
class ssytri2_test  : public  ::testing::Test {
public:
   sytri2_float_parameters  *ssytri2_obj;
   void SetUp();  
   void TearDown () { delete ssytri2_obj; }
};


void ssytri2_test::SetUp(){

    /* LAPACKE SSYTRI2 prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytri2) ( int matrix_layout, char uplo,
                                lapack_int n,  float * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_ssytri2 SSYTRI2;

    typedef int (*Fptr_NL_LAPACKE_ssytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, float *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_ssytrf SSYTRF;

    ssytri2_obj = new sytri2_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    ssytri2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytri2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytri2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytri2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SSYTRI2 = (Fptr_NL_LAPACKE_ssytri2)dlsym(ssytri2_obj->hModule, "LAPACKE_ssytri2");
    ASSERT_TRUE(SSYTRI2 != NULL) << "failed to get the Netlib LAPACKE_ssytri2 symbol";
    
    SSYTRF = (Fptr_NL_LAPACKE_ssytrf)dlsym(ssytri2_obj->hModule,"LAPACKE_ssytrf");
    ASSERT_TRUE(SSYTRF != NULL) << "failed to get the Netlib LAPACKE_ssytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri2 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    ssytri2_obj->inforef = SSYTRF(   ssytri2_obj->matrix_layout,
                                    ssytri2_obj->uplo,
                                    ssytri2_obj->n,
                                    ssytri2_obj->aref,
                                    ssytri2_obj->lda,
                                    ssytri2_obj->ipivref);

    ssytri2_obj->inforef = SSYTRI2(   ssytri2_obj->matrix_layout,
                                    ssytri2_obj->uplo, 
                                    ssytri2_obj->n,
                                    ssytri2_obj->aref,
                                    ssytri2_obj->lda,
                                    ssytri2_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    ssytri2_obj->info  = LAPACKE_ssytrf( ssytri2_obj->matrix_layout,
                                        ssytri2_obj->uplo,
                                        ssytri2_obj->n,
                                        ssytri2_obj->a,
                                        ssytri2_obj->lda,
                                        ssytri2_obj->ipiv);

    ssytri2_obj->info = LAPACKE_ssytri2(  ssytri2_obj->matrix_layout, 
                                        ssytri2_obj->uplo,
                                        ssytri2_obj->n, 
                                        ssytri2_obj->a, 
                                        ssytri2_obj->lda,
                                        ssytri2_obj->ipiv);

    if( ssytri2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssytri2 is wrong\n", ssytri2_obj->info );
    }
    if( ssytri2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytri2 is wrong\n", 
        ssytri2_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ssytri2_obj->diff =  computeDiff_s( (ssytri2_obj->n)*(ssytri2_obj->lda),
                           ssytri2_obj->a, ssytri2_obj->aref );
}

TEST_F(ssytri2_test, ssytri21) {
    EXPECT_NEAR(0.0, ssytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri2_test, ssytri22) {
    EXPECT_NEAR(0.0, ssytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri2_test, ssytri23) {
    EXPECT_NEAR(0.0, ssytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri2_test, ssytri24) {
    EXPECT_NEAR(0.0, ssytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytri2_double_parameters  class definition */
class sytri2_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Input/ Output parameters */
      double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri2_double_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri2_double_parameters (); 
};  /* end of sytri2_double_parameters  class definition */


/* Constructor sytri2_double_parameters definition */
sytri2_double_parameters:: sytri2_double_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytri2 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       sytri2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri2_double_parameters:: ~sytri2_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri2_double_parameters object: destructor invoked. \n");
#endif
   sytri2_free();
}

//  Test fixture class definition
class dsytri2_test  : public  ::testing::Test {
public:
   sytri2_double_parameters  *dsytri2_obj;
   void SetUp();  
   void TearDown () { delete dsytri2_obj; }
};


void dsytri2_test::SetUp(){

    /* LAPACKE DSYTRI2 prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytri2) ( int matrix_layout, char uplo,
                                lapack_int n,  double * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_dsytri2 DSYTRI2;

    typedef int (*Fptr_NL_LAPACKE_dsytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, double *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_dsytrf DSYTRF;

    dsytri2_obj = new sytri2_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    dsytri2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytri2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytri2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytri2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DSYTRI2 = (Fptr_NL_LAPACKE_dsytri2)dlsym(dsytri2_obj->hModule, "LAPACKE_dsytri2");
    ASSERT_TRUE(DSYTRI2 != NULL) << "failed to get the Netlib LAPACKE_dsytri2 symbol";
    
    DSYTRF = (Fptr_NL_LAPACKE_dsytrf)dlsym(dsytri2_obj->hModule,"LAPACKE_dsytrf");
    ASSERT_TRUE(DSYTRF != NULL) << "failed to get the Netlib LAPACKE_dsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri2 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dsytri2_obj->inforef = DSYTRF(   dsytri2_obj->matrix_layout,
                                    dsytri2_obj->uplo,
                                    dsytri2_obj->n,
                                    dsytri2_obj->aref,
                                    dsytri2_obj->lda,
                                    dsytri2_obj->ipivref);

    dsytri2_obj->inforef = DSYTRI2(   dsytri2_obj->matrix_layout,
                                    dsytri2_obj->uplo, 
                                    dsytri2_obj->n,
                                    dsytri2_obj->aref,
                                    dsytri2_obj->lda,
                                    dsytri2_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    dsytri2_obj->info  = LAPACKE_dsytrf( dsytri2_obj->matrix_layout,
                                        dsytri2_obj->uplo,
                                        dsytri2_obj->n,
                                        dsytri2_obj->a,
                                        dsytri2_obj->lda,
                                        dsytri2_obj->ipiv);

    dsytri2_obj->info = LAPACKE_dsytri2(  dsytri2_obj->matrix_layout, 
                                        dsytri2_obj->uplo,
                                        dsytri2_obj->n, 
                                        dsytri2_obj->a, 
                                        dsytri2_obj->lda,
                                        dsytri2_obj->ipiv);

    if( dsytri2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsytri2 is wrong\n", dsytri2_obj->info );
    }
    if( dsytri2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytri2 is wrong\n", 
        dsytri2_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dsytri2_obj->diff =  computeDiff_d( (dsytri2_obj->n)*(dsytri2_obj->lda),
                           dsytri2_obj->a, dsytri2_obj->aref );
}

TEST_F(dsytri2_test, dsytri21) {
    EXPECT_NEAR(0.0, dsytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri2_test, dsytri22) {
    EXPECT_NEAR(0.0, dsytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri2_test, dsytri23) {
    EXPECT_NEAR(0.0, dsytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri2_test, dsytri24) {
    EXPECT_NEAR(0.0, dsytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytri2_scomplex_parameters  class definition */
class sytri2_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri2_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri2_scomplex_parameters (); 
};  /* end of sytri2_scomplex_parameters  class definition */


/* Constructor sytri2_scomplex_parameters definition */
sytri2_scomplex_parameters:: sytri2_scomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytri2 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       sytri2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri2_scomplex_parameters:: ~sytri2_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri2_scomplex_parameters object: destructor invoked. \n");
#endif
   sytri2_free();
}

//  Test fixture class definition
class csytri2_test  : public  ::testing::Test {
public:
   sytri2_scomplex_parameters  *csytri2_obj;
   void SetUp();  
   void TearDown () { delete csytri2_obj; }
};


void csytri2_test::SetUp(){

    /* LAPACKE CSYTRI2 prototype */
    typedef int (*Fptr_NL_LAPACKE_csytri2) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_float * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_csytri2 CSYTRI2;

    typedef int (*Fptr_NL_LAPACKE_csytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_csytrf CSYTRF;

    csytri2_obj = new sytri2_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    csytri2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytri2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytri2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytri2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CSYTRI2 = (Fptr_NL_LAPACKE_csytri2)dlsym(csytri2_obj->hModule, "LAPACKE_csytri2");
    ASSERT_TRUE(CSYTRI2 != NULL) << "failed to get the Netlib LAPACKE_csytri2 symbol";
    
    CSYTRF = (Fptr_NL_LAPACKE_csytrf)dlsym(csytri2_obj->hModule,"LAPACKE_csytrf");
    ASSERT_TRUE(CSYTRF != NULL) << "failed to get the Netlib LAPACKE_csytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri2 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    csytri2_obj->inforef = CSYTRF(   csytri2_obj->matrix_layout,
                                    csytri2_obj->uplo,
                                    csytri2_obj->n,
                                    csytri2_obj->aref,
                                    csytri2_obj->lda,
                                    csytri2_obj->ipivref);

    csytri2_obj->inforef = CSYTRI2(   csytri2_obj->matrix_layout,
                                    csytri2_obj->uplo, 
                                    csytri2_obj->n,
                                    csytri2_obj->aref,
                                    csytri2_obj->lda,
                                    csytri2_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    csytri2_obj->info  = LAPACKE_csytrf( csytri2_obj->matrix_layout,
                                        csytri2_obj->uplo,
                                        csytri2_obj->n,
                                        csytri2_obj->a,
                                        csytri2_obj->lda,
                                        csytri2_obj->ipiv);

    csytri2_obj->info = LAPACKE_csytri2(  csytri2_obj->matrix_layout, 
                                        csytri2_obj->uplo,
                                        csytri2_obj->n, 
                                        csytri2_obj->a, 
                                        csytri2_obj->lda,
                                        csytri2_obj->ipiv);

    if( csytri2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csytri2 is wrong\n", csytri2_obj->info );
    }
    if( csytri2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytri2 is wrong\n", 
        csytri2_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    csytri2_obj->diff =  computeDiff_c( (csytri2_obj->n)*(csytri2_obj->lda),
                           csytri2_obj->a, csytri2_obj->aref );
}

TEST_F(csytri2_test, csytri21) {
    EXPECT_NEAR(0.0, csytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri2_test, csytri22) {
    EXPECT_NEAR(0.0, csytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri2_test, csytri23) {
    EXPECT_NEAR(0.0, csytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri2_test, csytri24) {
    EXPECT_NEAR(0.0, csytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytri2_dcomplex_parameters  class definition */
class sytri2_dcomplex_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Input/ Output parameters */
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri2_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri2_dcomplex_parameters (); 
};  /* end of sytri2_dcomplex_parameters  class definition */


/* Constructor sytri2_dcomplex_parameters definition */
sytri2_dcomplex_parameters:: sytri2_dcomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytri2 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       sytri2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri2_dcomplex_parameters:: ~sytri2_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri2_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytri2_free();
}

//  Test fixture class definition
class zsytri2_test  : public  ::testing::Test {
public:
   sytri2_dcomplex_parameters  *zsytri2_obj;
   void SetUp();  
   void TearDown () { delete zsytri2_obj; }
};


void zsytri2_test::SetUp(){

    /* LAPACKE ZSYTRI2 prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytri2) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_double * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_zsytri2 ZSYTRI2;

    typedef int (*Fptr_NL_LAPACKE_zsytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_zsytrf ZSYTRF;

    zsytri2_obj = new sytri2_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    zsytri2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytri2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytri2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytri2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZSYTRI2 = (Fptr_NL_LAPACKE_zsytri2)dlsym(zsytri2_obj->hModule, "LAPACKE_zsytri2");
    ASSERT_TRUE(ZSYTRI2 != NULL) << "failed to get the Netlib LAPACKE_zsytri2 symbol";
    
    ZSYTRF = (Fptr_NL_LAPACKE_zsytrf)dlsym(zsytri2_obj->hModule,"LAPACKE_zsytrf");
    ASSERT_TRUE(ZSYTRF != NULL) << "failed to get the Netlib LAPACKE_zsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri2 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zsytri2_obj->inforef = ZSYTRF(   zsytri2_obj->matrix_layout,
                                    zsytri2_obj->uplo,
                                    zsytri2_obj->n,
                                    zsytri2_obj->aref,
                                    zsytri2_obj->lda,
                                    zsytri2_obj->ipivref);

    zsytri2_obj->inforef = ZSYTRI2(   zsytri2_obj->matrix_layout,
                                    zsytri2_obj->uplo, 
                                    zsytri2_obj->n,
                                    zsytri2_obj->aref,
                                    zsytri2_obj->lda,
                                    zsytri2_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    zsytri2_obj->info  = LAPACKE_zsytrf( zsytri2_obj->matrix_layout,
                                        zsytri2_obj->uplo,
                                        zsytri2_obj->n,
                                        zsytri2_obj->a,
                                        zsytri2_obj->lda,
                                        zsytri2_obj->ipiv);

    zsytri2_obj->info = LAPACKE_zsytri2(  zsytri2_obj->matrix_layout, 
                                        zsytri2_obj->uplo,
                                        zsytri2_obj->n, 
                                        zsytri2_obj->a, 
                                        zsytri2_obj->lda,
                                        zsytri2_obj->ipiv);

    if( zsytri2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsytri2 is wrong\n", zsytri2_obj->info );
    }
    if( zsytri2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytri2 is wrong\n", 
        zsytri2_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zsytri2_obj->diff =  computeDiff_z( (zsytri2_obj->n)*(zsytri2_obj->lda),
                           zsytri2_obj->a, zsytri2_obj->aref );
}

TEST_F(zsytri2_test, zsytri21) {
    EXPECT_NEAR(0.0, zsytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri2_test, zsytri22) {
    EXPECT_NEAR(0.0, zsytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri2_test, zsytri23) {
    EXPECT_NEAR(0.0, zsytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri2_test, zsytri24) {
    EXPECT_NEAR(0.0, zsytri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

