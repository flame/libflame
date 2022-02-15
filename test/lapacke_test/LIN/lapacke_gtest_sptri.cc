#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define sptri_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (ipiv != NULL) free (ipiv); \
  if (ipivref != NULL)free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin sptri_float_parameters  class definition */
class sptri_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Input/ Output parameters */
      float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sptri_float_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i);
              
      ~sptri_float_parameters (); 
};  /* end of sptri_float_parameters  class definition */


/* Constructor sptri_float_parameters definition */
sptri_float_parameters:: sptri_float_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sptri Double:  n: %d, Uplo: %c \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       sptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sptri_float_parameters:: ~sptri_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptri_float_parameters object: destructor invoked. \n");
#endif
   sptri_free();
}

//  Test fixture class definition
class ssptri_test  : public  ::testing::Test {
public:
   sptri_float_parameters  *ssptri_obj;
   void SetUp();  
   void TearDown () { delete ssptri_obj; }
};


void ssptri_test::SetUp(){

    /* LAPACKE SSPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_ssptri) ( int matrix_layout, char uplo,
                                lapack_int n,  float * a,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_ssptri SSPTRI;

    typedef int (*Fptr_NL_LAPACKE_ssptrf) ( int matrix_layout ,char uplo,
                               lapack_int n, float *a,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_ssptrf SSPTRF;

    ssptri_obj = new sptri_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);

    ssptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SSPTRI = (Fptr_NL_LAPACKE_ssptri)dlsym(ssptri_obj->hModule, "LAPACKE_ssptri");
    ASSERT_TRUE(SSPTRI != NULL) << "failed to get the Netlib LAPACKE_ssptri symbol";
    
    SSPTRF = (Fptr_NL_LAPACKE_ssptrf)dlsym(ssptri_obj->hModule,"LAPACKE_ssptrf");
    ASSERT_TRUE(SSPTRF != NULL) << "failed to get the Netlib LAPACKE_ssptrf symbol";

    /* Pre condition: need to call sptrf - before calling sptri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    ssptri_obj->inforef = SSPTRF(   ssptri_obj->matrix_layout,
                                    ssptri_obj->uplo,
                                    ssptri_obj->n,
                                    ssptri_obj->aref,
                                    ssptri_obj->ipivref);

    ssptri_obj->inforef = SSPTRI(   ssptri_obj->matrix_layout,
                                    ssptri_obj->uplo, 
                                    ssptri_obj->n,
                                    ssptri_obj->aref,
                                    ssptri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    ssptri_obj->info  = LAPACKE_ssptrf( ssptri_obj->matrix_layout,
                                        ssptri_obj->uplo,
                                        ssptri_obj->n,
                                        ssptri_obj->a,
                                        ssptri_obj->ipiv);

    ssptri_obj->info = LAPACKE_ssptri(  ssptri_obj->matrix_layout, 
                                        ssptri_obj->uplo,
                                        ssptri_obj->n, 
                                        ssptri_obj->a, 
                                        ssptri_obj->ipiv);

    if( ssptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssptri is wrong\n", ssptri_obj->info );
    }
    if( ssptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssptri is wrong\n", 
        ssptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ssptri_obj->diff =  computeDiff_s( (ssptri_obj->n)*(ssptri_obj->n),
                           ssptri_obj->a, ssptri_obj->aref );
}

TEST_F(ssptri_test, ssptri1) {
    EXPECT_NEAR(0.0, ssptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssptri_test, ssptri2) {
    EXPECT_NEAR(0.0, ssptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssptri_test, ssptri3) {
    EXPECT_NEAR(0.0, ssptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssptri_test, ssptri4) {
    EXPECT_NEAR(0.0, ssptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sptri_double_parameters  class definition */
class sptri_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Input/ Output parameters */
      double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sptri_double_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i);
              
      ~sptri_double_parameters (); 
};  /* end of sptri_double_parameters  class definition */


/* Constructor sptri_double_parameters definition */
sptri_double_parameters:: sptri_double_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sptri Double:  n: %d, Uplo: %c \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       sptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sptri_double_parameters:: ~sptri_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptri_double_parameters object: destructor invoked. \n");
#endif
   sptri_free();
}

//  Test fixture class definition
class dsptri_test  : public  ::testing::Test {
public:
   sptri_double_parameters  *dsptri_obj;
   void SetUp();  
   void TearDown () { delete dsptri_obj; }
};


void dsptri_test::SetUp(){

    /* LAPACKE DSPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_dsptri) ( int matrix_layout, char uplo,
                                lapack_int n,  double * a,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_dsptri DSPTRI;

    typedef int (*Fptr_NL_LAPACKE_dsptrf) ( int matrix_layout ,char uplo,
                               lapack_int n, double *a,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_dsptrf DSPTRF;

    dsptri_obj = new sptri_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);

    dsptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DSPTRI = (Fptr_NL_LAPACKE_dsptri)dlsym(dsptri_obj->hModule, "LAPACKE_dsptri");
    ASSERT_TRUE(DSPTRI != NULL) << "failed to get the Netlib LAPACKE_dsptri symbol";
    
    DSPTRF = (Fptr_NL_LAPACKE_dsptrf)dlsym(dsptri_obj->hModule,"LAPACKE_dsptrf");
    ASSERT_TRUE(DSPTRF != NULL) << "failed to get the Netlib LAPACKE_dsptrf symbol";

    /* Pre condition: need to call sptrf - before calling sptri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dsptri_obj->inforef = DSPTRF(   dsptri_obj->matrix_layout,
                                    dsptri_obj->uplo,
                                    dsptri_obj->n,
                                    dsptri_obj->aref,
                                    dsptri_obj->ipivref);

    dsptri_obj->inforef = DSPTRI(   dsptri_obj->matrix_layout,
                                    dsptri_obj->uplo, 
                                    dsptri_obj->n,
                                    dsptri_obj->aref,
                                    dsptri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    dsptri_obj->info  = LAPACKE_dsptrf( dsptri_obj->matrix_layout,
                                        dsptri_obj->uplo,
                                        dsptri_obj->n,
                                        dsptri_obj->a,
                                        dsptri_obj->ipiv);

    dsptri_obj->info = LAPACKE_dsptri(  dsptri_obj->matrix_layout, 
                                        dsptri_obj->uplo,
                                        dsptri_obj->n, 
                                        dsptri_obj->a, 
                                        dsptri_obj->ipiv);

    if( dsptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsptri is wrong\n", dsptri_obj->info );
    }
    if( dsptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsptri is wrong\n", 
        dsptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dsptri_obj->diff =  computeDiff_d( (dsptri_obj->n)*(dsptri_obj->n),
                           dsptri_obj->a, dsptri_obj->aref );
}

TEST_F(dsptri_test, dsptri1) {
    EXPECT_NEAR(0.0, dsptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsptri_test, dsptri2) {
    EXPECT_NEAR(0.0, dsptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsptri_test, dsptri3) {
    EXPECT_NEAR(0.0, dsptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsptri_test, dsptri4) {
    EXPECT_NEAR(0.0, dsptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sptri_scomplex_parameters  class definition */
class sptri_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sptri_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i);
              
      ~sptri_scomplex_parameters (); 
};  /* end of sptri_scomplex_parameters  class definition */


/* Constructor sptri_scomplex_parameters definition */
sptri_scomplex_parameters:: sptri_scomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sptri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       sptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sptri_scomplex_parameters:: ~sptri_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptri_scomplex_parameters object: destructor invoked. \n");
#endif
   sptri_free();
}

//  Test fixture class definition
class csptri_test  : public  ::testing::Test {
public:
   sptri_scomplex_parameters  *csptri_obj;
   void SetUp();  
   void TearDown () { delete csptri_obj; }
};


void csptri_test::SetUp(){

    /* LAPACKE CSPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_csptri) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_float * a, 
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_csptri CSPTRI;

    typedef int (*Fptr_NL_LAPACKE_csptrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, 
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_csptrf CSPTRF;

    csptri_obj = new sptri_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);

    csptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CSPTRI = (Fptr_NL_LAPACKE_csptri)dlsym(csptri_obj->hModule, "LAPACKE_csptri");
    ASSERT_TRUE(CSPTRI != NULL) << "failed to get the Netlib LAPACKE_csptri symbol";
    
    CSPTRF = (Fptr_NL_LAPACKE_csptrf)dlsym(csptri_obj->hModule,"LAPACKE_csptrf");
    ASSERT_TRUE(CSPTRF != NULL) << "failed to get the Netlib LAPACKE_csptrf symbol";

    /* Pre condition: need to call sptrf - before calling sptri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    csptri_obj->inforef = CSPTRF(   csptri_obj->matrix_layout,
                                    csptri_obj->uplo,
                                    csptri_obj->n,
                                    csptri_obj->aref,
                                    csptri_obj->ipivref);

    csptri_obj->inforef = CSPTRI(   csptri_obj->matrix_layout,
                                    csptri_obj->uplo, 
                                    csptri_obj->n,
                                    csptri_obj->aref,
                                    csptri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    csptri_obj->info  = LAPACKE_csptrf( csptri_obj->matrix_layout,
                                        csptri_obj->uplo,
                                        csptri_obj->n,
                                        csptri_obj->a,
                                        csptri_obj->ipiv);

    csptri_obj->info = LAPACKE_csptri(  csptri_obj->matrix_layout, 
                                        csptri_obj->uplo,
                                        csptri_obj->n, 
                                        csptri_obj->a, 
                                        csptri_obj->ipiv);

    if( csptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csptri is wrong\n", csptri_obj->info );
    }
    if( csptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csptri is wrong\n", 
        csptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    csptri_obj->diff =  computeDiff_c( (csptri_obj->n)*(csptri_obj->n),
                           csptri_obj->a, csptri_obj->aref );
}

TEST_F(csptri_test, csptri1) {
    EXPECT_NEAR(0.0, csptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csptri_test, csptri2) {
    EXPECT_NEAR(0.0, csptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csptri_test, csptri3) {
    EXPECT_NEAR(0.0, csptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csptri_test, csptri4) {
    EXPECT_NEAR(0.0, csptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sptri_dcomplex_parameters  class definition */
class sptri_dcomplex_parameters{
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
      sptri_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i);
              
      ~sptri_dcomplex_parameters (); 
};  /* end of sptri_dcomplex_parameters  class definition */


/* Constructor sptri_dcomplex_parameters definition */
sptri_dcomplex_parameters:: sptri_dcomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sptri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       sptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sptri_dcomplex_parameters:: ~sptri_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptri_dcomplex_parameters object: destructor invoked. \n");
#endif
   sptri_free();
}

//  Test fixture class definition
class zsptri_test  : public  ::testing::Test {
public:
   sptri_dcomplex_parameters  *zsptri_obj;
   void SetUp();  
   void TearDown () { delete zsptri_obj; }
};


void zsptri_test::SetUp(){

    /* LAPACKE ZSPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_zsptri) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_double * a, 
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_zsptri ZSPTRI;

    typedef int (*Fptr_NL_LAPACKE_zsptrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, 
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_zsptrf ZSPTRF;

    zsptri_obj = new sptri_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);

    zsptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZSPTRI = (Fptr_NL_LAPACKE_zsptri)dlsym(zsptri_obj->hModule, "LAPACKE_zsptri");
    ASSERT_TRUE(ZSPTRI != NULL) << "failed to get the Netlib LAPACKE_zsptri symbol";
    
    ZSPTRF = (Fptr_NL_LAPACKE_zsptrf)dlsym(zsptri_obj->hModule,"LAPACKE_zsptrf");
    ASSERT_TRUE(ZSPTRF != NULL) << "failed to get the Netlib LAPACKE_zsptrf symbol";

    /* Pre condition: need to call sptrf - before calling sptri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zsptri_obj->inforef = ZSPTRF(   zsptri_obj->matrix_layout,
                                    zsptri_obj->uplo,
                                    zsptri_obj->n,
                                    zsptri_obj->aref,
                                    zsptri_obj->ipivref);

    zsptri_obj->inforef = ZSPTRI(   zsptri_obj->matrix_layout,
                                    zsptri_obj->uplo, 
                                    zsptri_obj->n,
                                    zsptri_obj->aref,
                                    zsptri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    zsptri_obj->info  = LAPACKE_zsptrf( zsptri_obj->matrix_layout,
                                        zsptri_obj->uplo,
                                        zsptri_obj->n,
                                        zsptri_obj->a,
                                        zsptri_obj->ipiv);

    zsptri_obj->info = LAPACKE_zsptri(  zsptri_obj->matrix_layout, 
                                        zsptri_obj->uplo,
                                        zsptri_obj->n, 
                                        zsptri_obj->a, 
                                        zsptri_obj->ipiv);

    if( zsptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsptri is wrong\n", zsptri_obj->info );
    }
    if( zsptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsptri is wrong\n", 
        zsptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zsptri_obj->diff =  computeDiff_z( (zsptri_obj->n)*(zsptri_obj->n),
                           zsptri_obj->a, zsptri_obj->aref );
}

TEST_F(zsptri_test, zsptri1) {
    EXPECT_NEAR(0.0, zsptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsptri_test, zsptri2) {
    EXPECT_NEAR(0.0, zsptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsptri_test, zsptri3) {
    EXPECT_NEAR(0.0, zsptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsptri_test, zsptri4) {
    EXPECT_NEAR(0.0, zsptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

