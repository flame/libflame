#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define sytri_3_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (e != NULL)    free (e   ); \
  if (eref != NULL) free (eref); \
  if (ipiv != NULL) free (ipiv); \
  if (ipivref != NULL)free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin sytri_3_float_parameters  class definition */
class sytri_3_float_parameters{
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
      float *e, *eref; // superdiagonal (or subdiagonal) elements of
      //  the symmetric block diagonal matrix D
      
      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri_3_float_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri_3_float_parameters (); 
};  /* end of sytri_3_float_parameters  class definition */


/* Constructor sytri_3_float_parameters definition */
sytri_3_float_parameters:: sytri_3_float_parameters ( int matrix_layout_i, 
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
   printf(" \n sytri_3 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_float_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       sytri_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    lapacke_gtest_init_float_buffer_pair_rand( e, eref, n);    
   } /* end of Constructor  */

sytri_3_float_parameters:: ~sytri_3_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri_3_float_parameters object: destructor invoked. \n");
#endif
   sytri_3_free();
}

//  Test fixture class definition
class ssytri_3_test  : public  ::testing::Test {
public:
   sytri_3_float_parameters  *ssytri_3_obj;
   void SetUp();  
   void TearDown () { delete ssytri_3_obj; }
};


void ssytri_3_test::SetUp(){

    /* LAPACKE SSYTRI_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytri_3) ( int matrix_layout, char uplo,
                                lapack_int n,  float *a, lapack_int lda,
                                const float *e, const lapack_int *ipiv  );
    Fptr_NL_LAPACKE_ssytri_3 SSYTRI_3;

    typedef int (*Fptr_NL_LAPACKE_ssytrf_rk) ( int matrix_layout ,char uplo,
                               lapack_int n, float *a, lapack_int lda,
                                float *e, lapack_int *ipiv );
    Fptr_NL_LAPACKE_ssytrf_rk SSYTRF_RK;

    ssytri_3_obj = new sytri_3_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    ssytri_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytri_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytri_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytri_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SSYTRI_3 = (Fptr_NL_LAPACKE_ssytri_3)dlsym(ssytri_3_obj->hModule, "LAPACKE_ssytri_3");
    ASSERT_TRUE(SSYTRI_3 != NULL) << "failed to get the Netlib LAPACKE_ssytri_3 symbol";
    
    SSYTRF_RK = (Fptr_NL_LAPACKE_ssytrf_rk)dlsym(ssytri_3_obj->hModule,"LAPACKE_ssytrf_rk");
    ASSERT_TRUE(SSYTRF_RK != NULL) << "failed to get the Netlib LAPACKE_ssytrf_rk symbol";

    /* Pre condition: need to call sytrf - before calling sytri_3 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    ssytri_3_obj->inforef = SSYTRF_RK(   ssytri_3_obj->matrix_layout,
                                    ssytri_3_obj->uplo,
                                    ssytri_3_obj->n,
                                    ssytri_3_obj->aref,
                                    ssytri_3_obj->lda,
                                    ssytri_3_obj->eref,
                                    ssytri_3_obj->ipivref);

    ssytri_3_obj->inforef = SSYTRI_3(   ssytri_3_obj->matrix_layout,
                                    ssytri_3_obj->uplo, 
                                    ssytri_3_obj->n,
                                    ssytri_3_obj->aref,
                                    ssytri_3_obj->lda,
                                    ssytri_3_obj->eref,
                                    ssytri_3_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    ssytri_3_obj->info  = LAPACKE_ssytrf_rk( ssytri_3_obj->matrix_layout,
                                        ssytri_3_obj->uplo,
                                        ssytri_3_obj->n,
                                        ssytri_3_obj->a,
                                        ssytri_3_obj->lda,
                                        ssytri_3_obj->e,
                                        ssytri_3_obj->ipiv);

    ssytri_3_obj->info = LAPACKE_ssytri_3(  ssytri_3_obj->matrix_layout, 
                                        ssytri_3_obj->uplo,
                                        ssytri_3_obj->n, 
                                        ssytri_3_obj->a, 
                                        ssytri_3_obj->lda,
                                        ssytri_3_obj->e,
                                        ssytri_3_obj->ipiv);

    if( ssytri_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ssytri_3 is wrong\n", ssytri_3_obj->info );
    }
    if( ssytri_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytri_3 is wrong\n", 
        ssytri_3_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ssytri_3_obj->diff =  computeDiff_s( (ssytri_3_obj->n)*(ssytri_3_obj->lda),
                           ssytri_3_obj->a, ssytri_3_obj->aref );
}

TEST_F(ssytri_3_test, ssytri_31) {
    EXPECT_NEAR(0.0, ssytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri_3_test, ssytri_32) {
    EXPECT_NEAR(0.0, ssytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri_3_test, ssytri_33) {
    EXPECT_NEAR(0.0, ssytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytri_3_test, ssytri_34) {
    EXPECT_NEAR(0.0, ssytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
/* Begin sytri_3_double_parameters  class definition */
class sytri_3_double_parameters{
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
      double *e, *eref; // superdiagonal (or subdiagonal) elements of
      //  the symmetric block diagonal matrix D
      
      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri_3_double_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri_3_double_parameters (); 
};  /* end of sytri_3_double_parameters  class definition */


/* Constructor sytri_3_double_parameters definition */
sytri_3_double_parameters:: sytri_3_double_parameters ( int matrix_layout_i, 
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
   printf(" \n sytri_3 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_double_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       sytri_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    lapacke_gtest_init_double_buffer_pair_rand( e, eref, n);    
   } /* end of Constructor  */

sytri_3_double_parameters:: ~sytri_3_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri_3_double_parameters object: destructor invoked. \n");
#endif
   sytri_3_free();
}

//  Test fixture class definition
class dsytri_3_test  : public  ::testing::Test {
public:
   sytri_3_double_parameters  *dsytri_3_obj;
   void SetUp();  
   void TearDown () { delete dsytri_3_obj; }
};


void dsytri_3_test::SetUp(){

    /* LAPACKE DSYTRI_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytri_3) ( int matrix_layout, char uplo,
                                lapack_int n,  double *a, lapack_int lda,
                                const double *e, const lapack_int *ipiv  );
    Fptr_NL_LAPACKE_dsytri_3 DSYTRI_3;

    typedef int (*Fptr_NL_LAPACKE_dsytrf_rk) ( int matrix_layout ,char uplo,
                               lapack_int n, double *a, lapack_int lda,
                                double *e, lapack_int *ipiv );
    Fptr_NL_LAPACKE_dsytrf_rk DSYTRF_RK;

    dsytri_3_obj = new sytri_3_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    dsytri_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytri_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytri_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytri_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DSYTRI_3 = (Fptr_NL_LAPACKE_dsytri_3)dlsym(dsytri_3_obj->hModule, "LAPACKE_dsytri_3");
    ASSERT_TRUE(DSYTRI_3 != NULL) << "failed to get the Netlib LAPACKE_dsytri_3 symbol";
    
    DSYTRF_RK = (Fptr_NL_LAPACKE_dsytrf_rk)dlsym(dsytri_3_obj->hModule,"LAPACKE_dsytrf_rk");
    ASSERT_TRUE(DSYTRF_RK != NULL) << "failed to get the Netlib LAPACKE_dsytrf_rk symbol";

    /* Pre condition: need to call sytrf - before calling sytri_3 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dsytri_3_obj->inforef = DSYTRF_RK(   dsytri_3_obj->matrix_layout,
                                    dsytri_3_obj->uplo,
                                    dsytri_3_obj->n,
                                    dsytri_3_obj->aref,
                                    dsytri_3_obj->lda,
                                    dsytri_3_obj->eref,
                                    dsytri_3_obj->ipivref);

    dsytri_3_obj->inforef = DSYTRI_3(   dsytri_3_obj->matrix_layout,
                                    dsytri_3_obj->uplo, 
                                    dsytri_3_obj->n,
                                    dsytri_3_obj->aref,
                                    dsytri_3_obj->lda,
                                    dsytri_3_obj->eref,
                                    dsytri_3_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    dsytri_3_obj->info  = LAPACKE_dsytrf_rk( dsytri_3_obj->matrix_layout,
                                        dsytri_3_obj->uplo,
                                        dsytri_3_obj->n,
                                        dsytri_3_obj->a,
                                        dsytri_3_obj->lda,
                                        dsytri_3_obj->e,
                                        dsytri_3_obj->ipiv);

    dsytri_3_obj->info = LAPACKE_dsytri_3(  dsytri_3_obj->matrix_layout, 
                                        dsytri_3_obj->uplo,
                                        dsytri_3_obj->n, 
                                        dsytri_3_obj->a, 
                                        dsytri_3_obj->lda,
                                        dsytri_3_obj->e,
                                        dsytri_3_obj->ipiv);

    if( dsytri_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dsytri_3 is wrong\n", dsytri_3_obj->info );
    }
    if( dsytri_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytri_3 is wrong\n", 
        dsytri_3_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dsytri_3_obj->diff =  computeDiff_d( (dsytri_3_obj->n)*(dsytri_3_obj->lda),
                           dsytri_3_obj->a, dsytri_3_obj->aref );
}

TEST_F(dsytri_3_test, dsytri_31) {
    EXPECT_NEAR(0.0, dsytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri_3_test, dsytri_32) {
    EXPECT_NEAR(0.0, dsytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri_3_test, dsytri_33) {
    EXPECT_NEAR(0.0, dsytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytri_3_test, dsytri_34) {
    EXPECT_NEAR(0.0, dsytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytri_3_scomplex_parameters  class definition */
class sytri_3_scomplex_parameters{
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
      lapack_complex_float *a, *aref; // The array ab contains the matrix A
      lapack_complex_float *e, *eref; // superdiagonal (or subdiagonal) elements of
      //  the symmetric block diagonal matrix D
      
      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri_3_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri_3_scomplex_parameters (); 
};  /* end of sytri_3_scomplex_parameters  class definition */


/* Constructor sytri_3_scomplex_parameters definition */
sytri_3_scomplex_parameters:: sytri_3_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n sytri_3 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       sytri_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( e, eref, n);    
   } /* end of Constructor  */

sytri_3_scomplex_parameters:: ~sytri_3_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri_3_scomplex_parameters object: destructor invoked. \n");
#endif
   sytri_3_free();
}

//  Test fixture class definition
class csytri_3_test  : public  ::testing::Test {
public:
   sytri_3_scomplex_parameters  *csytri_3_obj;
   void SetUp();  
   void TearDown () { delete csytri_3_obj; }
};


void csytri_3_test::SetUp(){

    /* LAPACKE CSYTRI_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_csytri_3) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_float *a, lapack_int lda,
                                const lapack_complex_float *e, const lapack_int *ipiv  );
    Fptr_NL_LAPACKE_csytri_3 CSYTRI_3;

    typedef int (*Fptr_NL_LAPACKE_csytrf_rk) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda,
                                lapack_complex_float *e, lapack_int *ipiv );
    Fptr_NL_LAPACKE_csytrf_rk CSYTRF_RK;

    csytri_3_obj = new sytri_3_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    csytri_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytri_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytri_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytri_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CSYTRI_3 = (Fptr_NL_LAPACKE_csytri_3)dlsym(csytri_3_obj->hModule, "LAPACKE_csytri_3");
    ASSERT_TRUE(CSYTRI_3 != NULL) << "failed to get the Netlib LAPACKE_csytri_3 symbol";
    
    CSYTRF_RK = (Fptr_NL_LAPACKE_csytrf_rk)dlsym(csytri_3_obj->hModule,"LAPACKE_csytrf_rk");
    ASSERT_TRUE(CSYTRF_RK != NULL) << "failed to get the Netlib LAPACKE_csytrf_rk symbol";

    /* Pre condition: need to call sytrf - before calling sytri_3 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    csytri_3_obj->inforef = CSYTRF_RK(   csytri_3_obj->matrix_layout,
                                    csytri_3_obj->uplo,
                                    csytri_3_obj->n,
                                    csytri_3_obj->aref,
                                    csytri_3_obj->lda,
                                    csytri_3_obj->eref,
                                    csytri_3_obj->ipivref);

    csytri_3_obj->inforef = CSYTRI_3(   csytri_3_obj->matrix_layout,
                                    csytri_3_obj->uplo, 
                                    csytri_3_obj->n,
                                    csytri_3_obj->aref,
                                    csytri_3_obj->lda,
                                    csytri_3_obj->eref,
                                    csytri_3_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    csytri_3_obj->info  = LAPACKE_csytrf_rk( csytri_3_obj->matrix_layout,
                                        csytri_3_obj->uplo,
                                        csytri_3_obj->n,
                                        csytri_3_obj->a,
                                        csytri_3_obj->lda,
                                        csytri_3_obj->e,
                                        csytri_3_obj->ipiv);

    csytri_3_obj->info = LAPACKE_csytri_3(  csytri_3_obj->matrix_layout, 
                                        csytri_3_obj->uplo,
                                        csytri_3_obj->n, 
                                        csytri_3_obj->a, 
                                        csytri_3_obj->lda,
                                        csytri_3_obj->e,
                                        csytri_3_obj->ipiv);

    if( csytri_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csytri_3 is wrong\n", csytri_3_obj->info );
    }
    if( csytri_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytri_3 is wrong\n", 
        csytri_3_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    csytri_3_obj->diff =  computeDiff_c( (csytri_3_obj->n)*(csytri_3_obj->lda),
                           csytri_3_obj->a, csytri_3_obj->aref );
}

TEST_F(csytri_3_test, csytri_31) {
    EXPECT_NEAR(0.0, csytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri_3_test, csytri_32) {
    EXPECT_NEAR(0.0, csytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri_3_test, csytri_33) {
    EXPECT_NEAR(0.0, csytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytri_3_test, csytri_34) {
    EXPECT_NEAR(0.0, csytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sytri_3_dcomplex_parameters  class definition */
class sytri_3_dcomplex_parameters{
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
      lapack_complex_double *a, *aref; // The array ab contains the matrix A
      lapack_complex_double *e, *eref; // superdiagonal (or subdiagonal) elements of
      //  the symmetric block diagonal matrix D
      
      /* Return Values */
      lapack_int info, inforef;

   public: 
      sytri_3_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~sytri_3_dcomplex_parameters (); 
};  /* end of sytri_3_dcomplex_parameters  class definition */


/* Constructor sytri_3_dcomplex_parameters definition */
sytri_3_dcomplex_parameters:: sytri_3_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n sytri_3 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       sytri_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( e, eref, n);    
   } /* end of Constructor  */

sytri_3_dcomplex_parameters:: ~sytri_3_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytri_3_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytri_3_free();
}

//  Test fixture class definition
class zsytri_3_test  : public  ::testing::Test {
public:
   sytri_3_dcomplex_parameters  *zsytri_3_obj;
   void SetUp();  
   void TearDown () { delete zsytri_3_obj; }
};


void zsytri_3_test::SetUp(){

    /* LAPACKE ZSYTRI_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytri_3) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_double *a, lapack_int lda,
                                const lapack_complex_double *e, const lapack_int *ipiv  );
    Fptr_NL_LAPACKE_zsytri_3 ZSYTRI_3;

    typedef int (*Fptr_NL_LAPACKE_zsytrf_rk) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda,
                                lapack_complex_double *e, lapack_int *ipiv );
    Fptr_NL_LAPACKE_zsytrf_rk ZSYTRF_RK;

    zsytri_3_obj = new sytri_3_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    zsytri_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytri_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytri_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytri_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZSYTRI_3 = (Fptr_NL_LAPACKE_zsytri_3)dlsym(zsytri_3_obj->hModule, "LAPACKE_zsytri_3");
    ASSERT_TRUE(ZSYTRI_3 != NULL) << "failed to get the Netlib LAPACKE_zsytri_3 symbol";
    
    ZSYTRF_RK = (Fptr_NL_LAPACKE_zsytrf_rk)dlsym(zsytri_3_obj->hModule,"LAPACKE_zsytrf_rk");
    ASSERT_TRUE(ZSYTRF_RK != NULL) << "failed to get the Netlib LAPACKE_zsytrf_rk symbol";

    /* Pre condition: need to call sytrf - before calling sytri_3 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zsytri_3_obj->inforef = ZSYTRF_RK(   zsytri_3_obj->matrix_layout,
                                    zsytri_3_obj->uplo,
                                    zsytri_3_obj->n,
                                    zsytri_3_obj->aref,
                                    zsytri_3_obj->lda,
                                    zsytri_3_obj->eref,
                                    zsytri_3_obj->ipivref);

    zsytri_3_obj->inforef = ZSYTRI_3(   zsytri_3_obj->matrix_layout,
                                    zsytri_3_obj->uplo, 
                                    zsytri_3_obj->n,
                                    zsytri_3_obj->aref,
                                    zsytri_3_obj->lda,
                                    zsytri_3_obj->eref,
                                    zsytri_3_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    zsytri_3_obj->info  = LAPACKE_zsytrf_rk( zsytri_3_obj->matrix_layout,
                                        zsytri_3_obj->uplo,
                                        zsytri_3_obj->n,
                                        zsytri_3_obj->a,
                                        zsytri_3_obj->lda,
                                        zsytri_3_obj->e,
                                        zsytri_3_obj->ipiv);

    zsytri_3_obj->info = LAPACKE_zsytri_3(  zsytri_3_obj->matrix_layout, 
                                        zsytri_3_obj->uplo,
                                        zsytri_3_obj->n, 
                                        zsytri_3_obj->a, 
                                        zsytri_3_obj->lda,
                                        zsytri_3_obj->e,
                                        zsytri_3_obj->ipiv);

    if( zsytri_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsytri_3 is wrong\n", zsytri_3_obj->info );
    }
    if( zsytri_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytri_3 is wrong\n", 
        zsytri_3_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zsytri_3_obj->diff =  computeDiff_z( (zsytri_3_obj->n)*(zsytri_3_obj->lda),
                           zsytri_3_obj->a, zsytri_3_obj->aref );
}

TEST_F(zsytri_3_test, zsytri_31) {
    EXPECT_NEAR(0.0, zsytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri_3_test, zsytri_32) {
    EXPECT_NEAR(0.0, zsytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri_3_test, zsytri_33) {
    EXPECT_NEAR(0.0, zsytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytri_3_test, zsytri_34) {
    EXPECT_NEAR(0.0, zsytri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


