#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define sytri2x_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (ipiv != NULL) free (ipiv); \
  if (ipivref != NULL)free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin sytri2x_float_parameters  class definition */
class sytri2x_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int nb; // block size
      /* Input/ Output parameters */
      float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri2x_float_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i,
                                lapack_int nb_i);
              
      ~sytri2x_float_parameters (); 
};  /* end of sytri2x_float_parameters  class definition */


/* Constructor sytri2x_float_parameters definition */
sytri2x_float_parameters:: sytri2x_float_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i,
                       lapack_int nb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    nb = nb_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytri2x Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       sytri2x_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri2x_float_parameters:: ~sytri2x_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri2x_float_parameters object: destructor invoked. \n");
#endif
   sytri2x_free();
}

//  Test fixture class definition
class ssytri2x_test  : public  ::testing::Test {
public:
   sytri2x_float_parameters  *ssytri2x_obj;
   void SetUp();  
   void TearDown () { delete ssytri2x_obj; }
};


void ssytri2x_test::SetUp(){

    /* LAPACKE SSYTRI2X prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytri2x) ( int matrix_layout, char uplo,
                                lapack_int n,  float * a, lapack_int lda,
                                const lapack_int *ipiv, lapack_int nb  );

    Fptr_NL_LAPACKE_ssytri2x SSYTRI2X;

    typedef int (*Fptr_NL_LAPACKE_ssytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, float *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_ssytrf SSYTRF;

    ssytri2x_obj = new sytri2x_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].ku);

    idx = Circular_Increment_Index(idx);

    ssytri2x_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytri2x_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytri2x_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytri2x_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SSYTRI2X = (Fptr_NL_LAPACKE_ssytri2x)dlsym(ssytri2x_obj->hModule, "LAPACKE_ssytri2x");
    ASSERT_TRUE(SSYTRI2X != NULL) << "failed to get the Netlib LAPACKE_ssytri2x symbol";
    
    SSYTRF = (Fptr_NL_LAPACKE_ssytrf)dlsym(ssytri2x_obj->hModule,"LAPACKE_ssytrf");
    ASSERT_TRUE(SSYTRF != NULL) << "failed to get the Netlib LAPACKE_ssytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri2x function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    ssytri2x_obj->inforef = SSYTRF(   ssytri2x_obj->matrix_layout,
                                    ssytri2x_obj->uplo,
                                    ssytri2x_obj->n,
                                    ssytri2x_obj->aref,
                                    ssytri2x_obj->lda,
                                    ssytri2x_obj->ipivref);

    ssytri2x_obj->inforef = SSYTRI2X(   ssytri2x_obj->matrix_layout,
                                    ssytri2x_obj->uplo, 
                                    ssytri2x_obj->n,
                                    ssytri2x_obj->aref,
                                    ssytri2x_obj->lda,
                                    ssytri2x_obj->ipivref,
                                    ssytri2x_obj->nb);
    /* Compute libflame's Lapacke o/p  */
    ssytri2x_obj->info  = LAPACKE_ssytrf( ssytri2x_obj->matrix_layout,
                                        ssytri2x_obj->uplo,
                                        ssytri2x_obj->n,
                                        ssytri2x_obj->a,
                                        ssytri2x_obj->lda,
                                        ssytri2x_obj->ipiv);

    ssytri2x_obj->info = LAPACKE_ssytri2x(  ssytri2x_obj->matrix_layout, 
                                        ssytri2x_obj->uplo,
                                        ssytri2x_obj->n, 
                                        ssytri2x_obj->a, 
                                        ssytri2x_obj->lda,
                                        ssytri2x_obj->ipiv,
                                        ssytri2x_obj->nb );

    if( ssytri2x_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssytri2x is wrong\n", ssytri2x_obj->info );
    }
    if( ssytri2x_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytri2x is wrong\n", 
        ssytri2x_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ssytri2x_obj->diff =  computeDiff_s( (ssytri2x_obj->n)*(ssytri2x_obj->lda),
                           ssytri2x_obj->a, ssytri2x_obj->aref );
}

TEST_F(ssytri2x_test, ssytri2x1) {
    EXPECT_NEAR(0.0, ssytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri2x_test, ssytri2x2) {
    EXPECT_NEAR(0.0, ssytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri2x_test, ssytri2x3) {
    EXPECT_NEAR(0.0, ssytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri2x_test, ssytri2x4) {
    EXPECT_NEAR(0.0, ssytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytri2x_double_parameters  class definition */
class sytri2x_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int nb; // block size

      /* Input/ Output parameters */
      double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri2x_double_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i,
                                lapack_int nb_i);
      ~sytri2x_double_parameters (); 
};  /* end of sytri2x_double_parameters  class definition */


/* Constructor sytri2x_double_parameters definition */
sytri2x_double_parameters:: sytri2x_double_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nb_i){
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    nb = nb_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytri2x Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       sytri2x_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri2x_double_parameters:: ~sytri2x_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri2x_double_parameters object: destructor invoked. \n");
#endif
   sytri2x_free();
}

//  Test fixture class definition
class dsytri2x_test  : public  ::testing::Test {
public:
   sytri2x_double_parameters  *dsytri2x_obj;
   void SetUp();  
   void TearDown () { delete dsytri2x_obj; }
};


void dsytri2x_test::SetUp(){

    /* LAPACKE DSYTRI2X prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytri2x) ( int matrix_layout, char uplo,
                                lapack_int n,  double * a, lapack_int lda,
                                const lapack_int *ipiv, lapack_int nb  );

    Fptr_NL_LAPACKE_dsytri2x DSYTRI2X;

    typedef int (*Fptr_NL_LAPACKE_dsytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, double *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_dsytrf DSYTRF;

    dsytri2x_obj = new sytri2x_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].ku);
    idx = Circular_Increment_Index(idx);

    dsytri2x_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytri2x_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytri2x_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytri2x_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DSYTRI2X = (Fptr_NL_LAPACKE_dsytri2x)dlsym(dsytri2x_obj->hModule, "LAPACKE_dsytri2x");
    ASSERT_TRUE(DSYTRI2X != NULL) << "failed to get the Netlib LAPACKE_dsytri2x symbol";
    
    DSYTRF = (Fptr_NL_LAPACKE_dsytrf)dlsym(dsytri2x_obj->hModule,"LAPACKE_dsytrf");
    ASSERT_TRUE(DSYTRF != NULL) << "failed to get the Netlib LAPACKE_dsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri2x function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dsytri2x_obj->inforef = DSYTRF(   dsytri2x_obj->matrix_layout,
                                    dsytri2x_obj->uplo,
                                    dsytri2x_obj->n,
                                    dsytri2x_obj->aref,
                                    dsytri2x_obj->lda,
                                    dsytri2x_obj->ipivref);

    dsytri2x_obj->inforef = DSYTRI2X(   dsytri2x_obj->matrix_layout,
                                    dsytri2x_obj->uplo, 
                                    dsytri2x_obj->n,
                                    dsytri2x_obj->aref,
                                    dsytri2x_obj->lda,
                                    dsytri2x_obj->ipivref,
                                    dsytri2x_obj->nb );
                        

    /* Compute libflame's Lapacke o/p  */
    dsytri2x_obj->info  = LAPACKE_dsytrf( dsytri2x_obj->matrix_layout,
                                        dsytri2x_obj->uplo,
                                        dsytri2x_obj->n,
                                        dsytri2x_obj->a,
                                        dsytri2x_obj->lda,
                                        dsytri2x_obj->ipiv);

    dsytri2x_obj->info = LAPACKE_dsytri2x(  dsytri2x_obj->matrix_layout, 
                                        dsytri2x_obj->uplo,
                                        dsytri2x_obj->n, 
                                        dsytri2x_obj->a, 
                                        dsytri2x_obj->lda,
                                        dsytri2x_obj->ipiv,
                                        dsytri2x_obj->nb );

    if( dsytri2x_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsytri2x is wrong\n", dsytri2x_obj->info );
    }
    if( dsytri2x_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytri2x is wrong\n", 
        dsytri2x_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dsytri2x_obj->diff =  computeDiff_d( (dsytri2x_obj->n)*(dsytri2x_obj->lda),
                           dsytri2x_obj->a, dsytri2x_obj->aref );
}

TEST_F(dsytri2x_test, dsytri2x1) {
    EXPECT_NEAR(0.0, dsytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri2x_test, dsytri2x2) {
    EXPECT_NEAR(0.0, dsytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri2x_test, dsytri2x3) {
    EXPECT_NEAR(0.0, dsytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri2x_test, dsytri2x4) {
    EXPECT_NEAR(0.0, dsytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytri2x_scomplex_parameters  class definition */
class sytri2x_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int nb; // block size

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri2x_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i,
                                lapack_int nb_i);
      ~sytri2x_scomplex_parameters (); 
};  /* end of sytri2x_scomplex_parameters  class definition */


/* Constructor sytri2x_scomplex_parameters definition */
sytri2x_scomplex_parameters:: sytri2x_scomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nb_i){
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    nb = nb_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytri2x Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       sytri2x_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri2x_scomplex_parameters:: ~sytri2x_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri2x_scomplex_parameters object: destructor invoked. \n");
#endif
   sytri2x_free();
}

//  Test fixture class definition
class csytri2x_test  : public  ::testing::Test {
public:
   sytri2x_scomplex_parameters  *csytri2x_obj;
   void SetUp();  
   void TearDown () { delete csytri2x_obj; }
};


void csytri2x_test::SetUp(){

    /* LAPACKE CSYTRI2X prototype */
    typedef int (*Fptr_NL_LAPACKE_csytri2x) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_float * a, lapack_int lda,
                                const lapack_int *ipiv, lapack_int nb  );

    Fptr_NL_LAPACKE_csytri2x CSYTRI2X;

    typedef int (*Fptr_NL_LAPACKE_csytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_csytrf CSYTRF;

    csytri2x_obj = new sytri2x_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].ku);
    idx = Circular_Increment_Index(idx);

    csytri2x_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytri2x_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytri2x_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytri2x_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CSYTRI2X = (Fptr_NL_LAPACKE_csytri2x)dlsym(csytri2x_obj->hModule, "LAPACKE_csytri2x");
    ASSERT_TRUE(CSYTRI2X != NULL) << "failed to get the Netlib LAPACKE_csytri2x symbol";
    
    CSYTRF = (Fptr_NL_LAPACKE_csytrf)dlsym(csytri2x_obj->hModule,"LAPACKE_csytrf");
    ASSERT_TRUE(CSYTRF != NULL) << "failed to get the Netlib LAPACKE_csytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri2x function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    csytri2x_obj->inforef = CSYTRF(   csytri2x_obj->matrix_layout,
                                    csytri2x_obj->uplo,
                                    csytri2x_obj->n,
                                    csytri2x_obj->aref,
                                    csytri2x_obj->lda,
                                    csytri2x_obj->ipivref);

    csytri2x_obj->inforef = CSYTRI2X(   csytri2x_obj->matrix_layout,
                                    csytri2x_obj->uplo, 
                                    csytri2x_obj->n,
                                    csytri2x_obj->aref,
                                    csytri2x_obj->lda,
                                    csytri2x_obj->ipivref,
                                    csytri2x_obj->nb );
                        

    /* Compute libflame's Lapacke o/p  */
    csytri2x_obj->info  = LAPACKE_csytrf( csytri2x_obj->matrix_layout,
                                        csytri2x_obj->uplo,
                                        csytri2x_obj->n,
                                        csytri2x_obj->a,
                                        csytri2x_obj->lda,
                                        csytri2x_obj->ipiv);

    csytri2x_obj->info = LAPACKE_csytri2x(  csytri2x_obj->matrix_layout, 
                                        csytri2x_obj->uplo,
                                        csytri2x_obj->n, 
                                        csytri2x_obj->a, 
                                        csytri2x_obj->lda,
                                        csytri2x_obj->ipiv,
                                        csytri2x_obj->nb );

    if( csytri2x_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csytri2x is wrong\n", csytri2x_obj->info );
    }
    if( csytri2x_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytri2x is wrong\n", 
        csytri2x_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    csytri2x_obj->diff =  computeDiff_c( (csytri2x_obj->n)*(csytri2x_obj->lda),
                           csytri2x_obj->a, csytri2x_obj->aref );
}

TEST_F(csytri2x_test, csytri2x1) {
    EXPECT_NEAR(0.0, csytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri2x_test, csytri2x2) {
    EXPECT_NEAR(0.0, csytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri2x_test, csytri2x3) {
    EXPECT_NEAR(0.0, csytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri2x_test, csytri2x4) {
    EXPECT_NEAR(0.0, csytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytri2x_dcomplex_parameters  class definition */
class sytri2x_dcomplex_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int nb; // block size

      /* Input/ Output parameters */
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri2x_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i,
                                lapack_int nb_i);
      ~sytri2x_dcomplex_parameters (); 
};  /* end of sytri2x_dcomplex_parameters  class definition */


/* Constructor sytri2x_dcomplex_parameters definition */
sytri2x_dcomplex_parameters:: sytri2x_dcomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nb_i){
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    nb = nb_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytri2x Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       sytri2x_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

sytri2x_dcomplex_parameters:: ~sytri2x_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri2x_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytri2x_free();
}

//  Test fixture class definition
class zsytri2x_test  : public  ::testing::Test {
public:
   sytri2x_dcomplex_parameters  *zsytri2x_obj;
   void SetUp();  
   void TearDown () { delete zsytri2x_obj; }
};


void zsytri2x_test::SetUp(){

    /* LAPACKE ZSYTRI2X prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytri2x) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_double * a, lapack_int lda,
                                const lapack_int *ipiv, lapack_int nb  );

    Fptr_NL_LAPACKE_zsytri2x ZSYTRI2X;

    typedef int (*Fptr_NL_LAPACKE_zsytrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_zsytrf ZSYTRF;

    zsytri2x_obj = new sytri2x_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].ku);
    idx = Circular_Increment_Index(idx);

    zsytri2x_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytri2x_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytri2x_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytri2x_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZSYTRI2X = (Fptr_NL_LAPACKE_zsytri2x)dlsym(zsytri2x_obj->hModule, "LAPACKE_zsytri2x");
    ASSERT_TRUE(ZSYTRI2X != NULL) << "failed to get the Netlib LAPACKE_zsytri2x symbol";
    
    ZSYTRF = (Fptr_NL_LAPACKE_zsytrf)dlsym(zsytri2x_obj->hModule,"LAPACKE_zsytrf");
    ASSERT_TRUE(ZSYTRF != NULL) << "failed to get the Netlib LAPACKE_zsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytri2x function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zsytri2x_obj->inforef = ZSYTRF(   zsytri2x_obj->matrix_layout,
                                    zsytri2x_obj->uplo,
                                    zsytri2x_obj->n,
                                    zsytri2x_obj->aref,
                                    zsytri2x_obj->lda,
                                    zsytri2x_obj->ipivref);

    zsytri2x_obj->inforef = ZSYTRI2X(   zsytri2x_obj->matrix_layout,
                                    zsytri2x_obj->uplo, 
                                    zsytri2x_obj->n,
                                    zsytri2x_obj->aref,
                                    zsytri2x_obj->lda,
                                    zsytri2x_obj->ipivref,
                                    zsytri2x_obj->nb );
                        

    /* Compute libflame's Lapacke o/p  */
    zsytri2x_obj->info  = LAPACKE_zsytrf( zsytri2x_obj->matrix_layout,
                                        zsytri2x_obj->uplo,
                                        zsytri2x_obj->n,
                                        zsytri2x_obj->a,
                                        zsytri2x_obj->lda,
                                        zsytri2x_obj->ipiv);

    zsytri2x_obj->info = LAPACKE_zsytri2x(  zsytri2x_obj->matrix_layout, 
                                        zsytri2x_obj->uplo,
                                        zsytri2x_obj->n, 
                                        zsytri2x_obj->a, 
                                        zsytri2x_obj->lda,
                                        zsytri2x_obj->ipiv,
                                        zsytri2x_obj->nb );

    if( zsytri2x_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsytri2x is wrong\n", zsytri2x_obj->info );
    }
    if( zsytri2x_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytri2x is wrong\n", 
        zsytri2x_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zsytri2x_obj->diff =  computeDiff_z( (zsytri2x_obj->n)*(zsytri2x_obj->lda),
                           zsytri2x_obj->a, zsytri2x_obj->aref );
}

TEST_F(zsytri2x_test, zsytri2x1) {
    EXPECT_NEAR(0.0, zsytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri2x_test, zsytri2x2) {
    EXPECT_NEAR(0.0, zsytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri2x_test, zsytri2x3) {
    EXPECT_NEAR(0.0, zsytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri2x_test, zsytri2x4) {
    EXPECT_NEAR(0.0, zsytri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

