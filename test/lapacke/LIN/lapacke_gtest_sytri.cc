#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define sytri_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (ipiv != NULL) free (ipiv); \
  if (ipivref != NULL)free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin sytri_float_parameters  class definition */
class sytri_float_parameters{
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
      sytri_float_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri_float_parameters (); 
};  /* end of sytri_float_parameters  class definition */


/* Constructor sytri_float_parameters definition */
sytri_float_parameters:: sytri_float_parameters ( int matrix_layout_i, 
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
   printf(" \n sytri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       sytri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri_float_parameters:: ~sytri_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri_float_parameters object: destructor invoked. \n");
#endif
   sytri_free();
}

//  Test fixture class definition
class ssytri_test  : public  ::testing::Test {
public:
   sytri_float_parameters  *ssytri_obj;
   void SetUp();  
   void TearDown () { delete ssytri_obj; }
};


void ssytri_test::SetUp(){

    /* LAPACKE SSYTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytri) ( int matrix_layout, char uplo,
                                lapack_int n,  float * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_ssytri SSYTRI;

    typedef int (*Fptr_NL_LAPACKE_ssytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, float *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_ssytrf SSYTRF;

    ssytri_obj = new sytri_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    ssytri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SSYTRI = (Fptr_NL_LAPACKE_ssytri)dlsym(ssytri_obj->hModule, "LAPACKE_ssytri");
    ASSERT_TRUE(SSYTRI != NULL) << "failed to get the Netlib LAPACKE_ssytri symbol";
    
    SSYTRF = (Fptr_NL_LAPACKE_ssytrf)dlsym(ssytri_obj->hModule,"LAPACKE_ssytrf");
    ASSERT_TRUE(SSYTRF != NULL) << "failed to get the Netlib LAPACKE_ssytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    ssytri_obj->inforef = SSYTRF(   ssytri_obj->matrix_layout,
                                    ssytri_obj->uplo,
                                    ssytri_obj->n,
                                    ssytri_obj->aref,
                                    ssytri_obj->lda,
                                    ssytri_obj->ipivref);

    ssytri_obj->inforef = SSYTRI(   ssytri_obj->matrix_layout,
                                    ssytri_obj->uplo, 
                                    ssytri_obj->n,
                                    ssytri_obj->aref,
                                    ssytri_obj->lda,
                                    ssytri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    ssytri_obj->info  = LAPACKE_ssytrf( ssytri_obj->matrix_layout,
                                        ssytri_obj->uplo,
                                        ssytri_obj->n,
                                        ssytri_obj->a,
                                        ssytri_obj->lda,
                                        ssytri_obj->ipiv);

    ssytri_obj->info = LAPACKE_ssytri(  ssytri_obj->matrix_layout, 
                                        ssytri_obj->uplo,
                                        ssytri_obj->n, 
                                        ssytri_obj->a, 
                                        ssytri_obj->lda,
                                        ssytri_obj->ipiv);

    if( ssytri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssytri is wrong\n", ssytri_obj->info );
    }
    if( ssytri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytri is wrong\n", 
        ssytri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ssytri_obj->diff =  computeDiff_s( (ssytri_obj->n)*(ssytri_obj->lda),
                           ssytri_obj->a, ssytri_obj->aref );
}

TEST_F(ssytri_test, ssytri1) {
    EXPECT_NEAR(0.0, ssytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri_test, ssytri2) {
    EXPECT_NEAR(0.0, ssytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri_test, ssytri3) {
    EXPECT_NEAR(0.0, ssytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri_test, ssytri4) {
    EXPECT_NEAR(0.0, ssytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytri_double_parameters  class definition */
class sytri_double_parameters{
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
      sytri_double_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri_double_parameters (); 
};  /* end of sytri_double_parameters  class definition */


/* Constructor sytri_double_parameters definition */
sytri_double_parameters:: sytri_double_parameters ( int matrix_layout_i, 
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
   printf(" \n sytri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       sytri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri_double_parameters:: ~sytri_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri_double_parameters object: destructor invoked. \n");
#endif
   sytri_free();
}

//  Test fixture class definition
class dsytri_test  : public  ::testing::Test {
public:
   sytri_double_parameters  *dsytri_obj;
   void SetUp();  
   void TearDown () { delete dsytri_obj; }
};


void dsytri_test::SetUp(){

    /* LAPACKE DSYTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytri) ( int matrix_layout, char uplo,
                                lapack_int n,  double * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_dsytri DSYTRI;

    typedef int (*Fptr_NL_LAPACKE_dsytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, double *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_dsytrf DSYTRF;

    dsytri_obj = new sytri_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    dsytri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DSYTRI = (Fptr_NL_LAPACKE_dsytri)dlsym(dsytri_obj->hModule, "LAPACKE_dsytri");
    ASSERT_TRUE(DSYTRI != NULL) << "failed to get the Netlib LAPACKE_dsytri symbol";
    
    DSYTRF = (Fptr_NL_LAPACKE_dsytrf)dlsym(dsytri_obj->hModule,"LAPACKE_dsytrf");
    ASSERT_TRUE(DSYTRF != NULL) << "failed to get the Netlib LAPACKE_dsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dsytri_obj->inforef = DSYTRF(   dsytri_obj->matrix_layout,
                                    dsytri_obj->uplo,
                                    dsytri_obj->n,
                                    dsytri_obj->aref,
                                    dsytri_obj->lda,
                                    dsytri_obj->ipivref);

    dsytri_obj->inforef = DSYTRI(   dsytri_obj->matrix_layout,
                                    dsytri_obj->uplo, 
                                    dsytri_obj->n,
                                    dsytri_obj->aref,
                                    dsytri_obj->lda,
                                    dsytri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    dsytri_obj->info  = LAPACKE_dsytrf( dsytri_obj->matrix_layout,
                                        dsytri_obj->uplo,
                                        dsytri_obj->n,
                                        dsytri_obj->a,
                                        dsytri_obj->lda,
                                        dsytri_obj->ipiv);

    dsytri_obj->info = LAPACKE_dsytri(  dsytri_obj->matrix_layout, 
                                        dsytri_obj->uplo,
                                        dsytri_obj->n, 
                                        dsytri_obj->a, 
                                        dsytri_obj->lda,
                                        dsytri_obj->ipiv);

    if( dsytri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsytri is wrong\n", dsytri_obj->info );
    }
    if( dsytri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytri is wrong\n", 
        dsytri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dsytri_obj->diff =  computeDiff_d( (dsytri_obj->n)*(dsytri_obj->lda),
                           dsytri_obj->a, dsytri_obj->aref );
}

TEST_F(dsytri_test, dsytri1) {
    EXPECT_NEAR(0.0, dsytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri_test, dsytri2) {
    EXPECT_NEAR(0.0, dsytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri_test, dsytri3) {
    EXPECT_NEAR(0.0, dsytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri_test, dsytri4) {
    EXPECT_NEAR(0.0, dsytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytri_scomplex_parameters  class definition */
class sytri_scomplex_parameters{
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
      sytri_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri_scomplex_parameters (); 
};  /* end of sytri_scomplex_parameters  class definition */


/* Constructor sytri_scomplex_parameters definition */
sytri_scomplex_parameters:: sytri_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n sytri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       sytri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri_scomplex_parameters:: ~sytri_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri_scomplex_parameters object: destructor invoked. \n");
#endif
   sytri_free();
}

//  Test fixture class definition
class csytri_test  : public  ::testing::Test {
public:
   sytri_scomplex_parameters  *csytri_obj;
   void SetUp();  
   void TearDown () { delete csytri_obj; }
};


void csytri_test::SetUp(){

    /* LAPACKE CSYTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_csytri) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_float * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_csytri CSYTRI;

    typedef int (*Fptr_NL_LAPACKE_csytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_csytrf CSYTRF;

    csytri_obj = new sytri_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    csytri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CSYTRI = (Fptr_NL_LAPACKE_csytri)dlsym(csytri_obj->hModule, "LAPACKE_csytri");
    ASSERT_TRUE(CSYTRI != NULL) << "failed to get the Netlib LAPACKE_csytri symbol";
    
    CSYTRF = (Fptr_NL_LAPACKE_csytrf)dlsym(csytri_obj->hModule,"LAPACKE_csytrf");
    ASSERT_TRUE(CSYTRF != NULL) << "failed to get the Netlib LAPACKE_csytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    csytri_obj->inforef = CSYTRF(   csytri_obj->matrix_layout,
                                    csytri_obj->uplo,
                                    csytri_obj->n,
                                    csytri_obj->aref,
                                    csytri_obj->lda,
                                    csytri_obj->ipivref);

    csytri_obj->inforef = CSYTRI(   csytri_obj->matrix_layout,
                                    csytri_obj->uplo, 
                                    csytri_obj->n,
                                    csytri_obj->aref,
                                    csytri_obj->lda,
                                    csytri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    csytri_obj->info  = LAPACKE_csytrf( csytri_obj->matrix_layout,
                                        csytri_obj->uplo,
                                        csytri_obj->n,
                                        csytri_obj->a,
                                        csytri_obj->lda,
                                        csytri_obj->ipiv);

    csytri_obj->info = LAPACKE_csytri(  csytri_obj->matrix_layout, 
                                        csytri_obj->uplo,
                                        csytri_obj->n, 
                                        csytri_obj->a, 
                                        csytri_obj->lda,
                                        csytri_obj->ipiv);

    if( csytri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csytri is wrong\n", csytri_obj->info );
    }
    if( csytri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytri is wrong\n", 
        csytri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    csytri_obj->diff =  computeDiff_c( (csytri_obj->n)*(csytri_obj->lda),
                           csytri_obj->a, csytri_obj->aref );
}

TEST_F(csytri_test, csytri1) {
    EXPECT_NEAR(0.0, csytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri_test, csytri2) {
    EXPECT_NEAR(0.0, csytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri_test, csytri3) {
    EXPECT_NEAR(0.0, csytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri_test, csytri4) {
    EXPECT_NEAR(0.0, csytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytri_dcomplex_parameters  class definition */
class sytri_dcomplex_parameters{
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
      sytri_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri_dcomplex_parameters (); 
};  /* end of sytri_dcomplex_parameters  class definition */


/* Constructor sytri_dcomplex_parameters definition */
sytri_dcomplex_parameters:: sytri_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n sytri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       sytri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri_dcomplex_parameters:: ~sytri_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytri_free();
}

//  Test fixture class definition
class zsytri_test  : public  ::testing::Test {
public:
   sytri_dcomplex_parameters  *zsytri_obj;
   void SetUp();  
   void TearDown () { delete zsytri_obj; }
};


void zsytri_test::SetUp(){

    /* LAPACKE ZSYTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytri) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_double * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_zsytri ZSYTRI;

    typedef int (*Fptr_NL_LAPACKE_zsytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_zsytrf ZSYTRF;

    zsytri_obj = new sytri_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    zsytri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZSYTRI = (Fptr_NL_LAPACKE_zsytri)dlsym(zsytri_obj->hModule, "LAPACKE_zsytri");
    ASSERT_TRUE(ZSYTRI != NULL) << "failed to get the Netlib LAPACKE_zsytri symbol";
    
    ZSYTRF = (Fptr_NL_LAPACKE_zsytrf)dlsym(zsytri_obj->hModule,"LAPACKE_zsytrf");
    ASSERT_TRUE(ZSYTRF != NULL) << "failed to get the Netlib LAPACKE_zsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zsytri_obj->inforef = ZSYTRF(   zsytri_obj->matrix_layout,
                                    zsytri_obj->uplo,
                                    zsytri_obj->n,
                                    zsytri_obj->aref,
                                    zsytri_obj->lda,
                                    zsytri_obj->ipivref);

    zsytri_obj->inforef = ZSYTRI(   zsytri_obj->matrix_layout,
                                    zsytri_obj->uplo, 
                                    zsytri_obj->n,
                                    zsytri_obj->aref,
                                    zsytri_obj->lda,
                                    zsytri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    zsytri_obj->info  = LAPACKE_zsytrf( zsytri_obj->matrix_layout,
                                        zsytri_obj->uplo,
                                        zsytri_obj->n,
                                        zsytri_obj->a,
                                        zsytri_obj->lda,
                                        zsytri_obj->ipiv);

    zsytri_obj->info = LAPACKE_zsytri(  zsytri_obj->matrix_layout, 
                                        zsytri_obj->uplo,
                                        zsytri_obj->n, 
                                        zsytri_obj->a, 
                                        zsytri_obj->lda,
                                        zsytri_obj->ipiv);

    if( zsytri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsytri is wrong\n", zsytri_obj->info );
    }
    if( zsytri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytri is wrong\n", 
        zsytri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zsytri_obj->diff =  computeDiff_z( (zsytri_obj->n)*(zsytri_obj->lda),
                           zsytri_obj->a, zsytri_obj->aref );
}

TEST_F(zsytri_test, zsytri1) {
    EXPECT_NEAR(0.0, zsytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri_test, zsytri2) {
    EXPECT_NEAR(0.0, zsytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri_test, zsytri3) {
    EXPECT_NEAR(0.0, zsytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri_test, zsytri4) {
    EXPECT_NEAR(0.0, zsytri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

