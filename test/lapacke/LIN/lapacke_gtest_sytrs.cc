#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define LAPACKE_TEST_VERBOSE (1)
#define sytrs_free() \
       if (ipiv != NULL) free (ipiv); \
       if (bref != NULL) free (bref); \
       if (b != NULL)    free (b   ); \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (ipivref != NULL)free (ipivref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin sytrs_double_parameters  class definition */
class sytrs_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sytrs_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_double_parameters ();
};  /* end of sytrs_double_parameters  class definition */


/* Constructor sytrs_double_parameters definition */
sytrs_double_parameters:: sytrs_double_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytrs Double:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif


    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_double_parameters object: malloc error.";
       sytrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs_double_parameters:: ~sytrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_double_parameters object: destructor invoked. \n");
#endif
   sytrs_free();
}


//  Test fixture class definition
class dsytrs_test  : public  ::testing::Test {
public:
   sytrs_double_parameters  *dsytrs_obj;
   void SetUp();  
   void TearDown () { delete dsytrs_obj; }
};


void dsytrs_test::SetUp(){

    /* LAPACKE DSYTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs, const double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                                            double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dsytrs DSYTRS;

     /* LAPACKE DSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf) ( int matrix_layout, char uplo,
             lapack_int n, double* a, lapack_int lda, lapack_int* ipiv );

    Fptr_NL_LAPACKE_dsytrf DSYTRF;

    dsytrs_obj = new  sytrs_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dsytrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYTRS = (Fptr_NL_LAPACKE_dsytrs)dlsym(dsytrs_obj->hModule, "LAPACKE_dsytrs");
    ASSERT_TRUE(DSYTRS != NULL) << "failed to get the Netlib LAPACKE_dsytrs symbol";

    DSYTRF = (Fptr_NL_LAPACKE_dsytrf)dlsym(dsytrs_obj->hModule,"LAPACKE_dsytrf");
    ASSERT_TRUE(DSYTRF != NULL) << "failed to get the Netlib LAPACKE_dsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dsytrs_obj->inforef = DSYTRF( dsytrs_obj->matrix_layout,
                            dsytrs_obj->uplo, dsytrs_obj->n,
                                     dsytrs_obj->aref,
                      dsytrs_obj->lda, dsytrs_obj->ipivref);

    dsytrs_obj->inforef = DSYTRS( dsytrs_obj->matrix_layout,
                                  dsytrs_obj->uplo, dsytrs_obj->n,
                                  dsytrs_obj->nrhs,
                                  (const double *)dsytrs_obj->aref,
                              dsytrs_obj->lda, dsytrs_obj->ipivref,
                                dsytrs_obj->bref, dsytrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dsytrs_obj->info = LAPACKE_dsytrf( dsytrs_obj->matrix_layout,
                                 dsytrs_obj->uplo, dsytrs_obj->n,
                                     dsytrs_obj->a,
                               dsytrs_obj->lda, dsytrs_obj->ipiv);

    dsytrs_obj->info = LAPACKE_dsytrs( dsytrs_obj->matrix_layout,
                dsytrs_obj->uplo, dsytrs_obj->n, dsytrs_obj->nrhs,
                                  (const double *)dsytrs_obj->a,
                               dsytrs_obj->lda, dsytrs_obj->ipiv,
                                 dsytrs_obj->b, dsytrs_obj->ldb );


    if( dsytrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrs is wrong\n",
                    dsytrs_obj->info );
    }
    if( dsytrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrs is wrong\n",
        dsytrs_obj->inforef );
    }
}

TEST_F(dsytrs_test, dsytrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_obj->b_bufsize,
                           dsytrs_obj->b, dsytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsytrs_test, dsytrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_obj->b_bufsize,
                           dsytrs_obj->b, dsytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsytrs_test, dsytrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_obj->b_bufsize,
                           dsytrs_obj->b, dsytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsytrs_test, dsytrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_obj->b_bufsize,
                           dsytrs_obj->b, dsytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrs_float_parameters  class definition */
class sytrs_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sytrs_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~sytrs_float_parameters ();
};  /* end of sytrs_float_parameters  class definition */


/* Constructor sytrs_float_parameters definition */
sytrs_float_parameters:: sytrs_float_parameters ( int matrix_layout_i,
                       char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytrs float:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (lda*n));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_double_parameters object: malloc error.";
       sytrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

sytrs_float_parameters:: ~sytrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_float_parameters object: destructor invoked. \n");
#endif
   sytrs_free();
}


//  Test fixture class definition
class ssytrs_test  : public  ::testing::Test {
public:
   sytrs_float_parameters  *ssytrs_obj;
   void SetUp();  
   void TearDown () { delete ssytrs_obj; }
};


void ssytrs_test::SetUp(){

    /* LAPACKE SSYTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs, const float * a,
                                  lapack_int lda, const lapack_int * ipiv,
                                            float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_ssytrs SSYTRS;

     /* LAPACKE SSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf) ( int matrix_layout, char uplo,
                lapack_int n, float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_ssytrf SSYTRF;

    ssytrs_obj = new  sytrs_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    ssytrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYTRS = (Fptr_NL_LAPACKE_ssytrs)dlsym(ssytrs_obj->hModule, "LAPACKE_ssytrs");
    ASSERT_TRUE(SSYTRS != NULL) << "failed to get the Netlib LAPACKE_ssytrs symbol";

    SSYTRF = (Fptr_NL_LAPACKE_ssytrf)dlsym(ssytrs_obj->hModule,"LAPACKE_ssytrf");
    ASSERT_TRUE(SSYTRF != NULL) << "failed to get the Netlib LAPACKE_ssytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    ssytrs_obj->inforef = SSYTRF( ssytrs_obj->matrix_layout,
                                    ssytrs_obj->uplo, ssytrs_obj->n,
                                     ssytrs_obj->aref,
                               ssytrs_obj->lda, ssytrs_obj->ipivref);

    ssytrs_obj->inforef = SSYTRS( ssytrs_obj->matrix_layout,
                                  ssytrs_obj->uplo, ssytrs_obj->n,
                                  ssytrs_obj->nrhs,
                                  (const float *)ssytrs_obj->aref,
                                  ssytrs_obj->lda, ssytrs_obj->ipivref,
                                  ssytrs_obj->bref, ssytrs_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssytrs_obj->info = LAPACKE_ssytrf( ssytrs_obj->matrix_layout,
                                 ssytrs_obj->uplo, ssytrs_obj->n,
                                     ssytrs_obj->a,
                               ssytrs_obj->lda, ssytrs_obj->ipiv);

    ssytrs_obj->info = LAPACKE_ssytrs( ssytrs_obj->matrix_layout,
                ssytrs_obj->uplo, ssytrs_obj->n, ssytrs_obj->nrhs,
                                  (const float *)ssytrs_obj->a,
                               ssytrs_obj->lda, ssytrs_obj->ipiv,
                                 ssytrs_obj->b, ssytrs_obj->ldb );
    if( ssytrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrs is wrong\n",
                    ssytrs_obj->info );
    }
    if( ssytrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrs is wrong\n",
        ssytrs_obj->inforef );
    }
}

TEST_F(ssytrs_test, ssytrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs_obj->b_bufsize,
                           ssytrs_obj->b, ssytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytrs_test, ssytrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs_obj->b_bufsize,
                           ssytrs_obj->b, ssytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytrs_test, ssytrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs_obj->b_bufsize,
                           ssytrs_obj->b, ssytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytrs_test, ssytrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs_obj->b_bufsize,
                           ssytrs_obj->b, ssytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sytrs_scomplex_parameters  class definition */
class sytrs_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sytrs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_scomplex_parameters ();
};  /* end of sytrs_scomplex_parameters  class definition */


/* Constructor sytrs_scomplex_parameters definition */
sytrs_scomplex_parameters:: sytrs_scomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytrs scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_scomplex_parameters object: malloc error.";
       sytrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs_scomplex_parameters:: ~sytrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrs_free();
}


//  Test fixture class definition
class csytrs_test  : public  ::testing::Test {
public:
   sytrs_scomplex_parameters  *csytrs_obj;
   void SetUp();  
   void TearDown () { delete csytrs_obj; }
};


void csytrs_test::SetUp(){

    /* LAPACKE CSYTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_float * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_csytrs CSYTRS;

     /* LAPACKE CSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf) ( int matrix_layout, char uplo,
    lapack_int n, lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_csytrf CSYTRF;


    csytrs_obj = new  sytrs_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    csytrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYTRS = (Fptr_NL_LAPACKE_csytrs)dlsym(csytrs_obj->hModule, "LAPACKE_csytrs");
    ASSERT_TRUE(CSYTRS != NULL) << "failed to get the Netlib LAPACKE_csytrs symbol";

    CSYTRF = (Fptr_NL_LAPACKE_csytrf)dlsym(csytrs_obj->hModule,"LAPACKE_csytrf");
    ASSERT_TRUE(CSYTRF != NULL) << "failed to get the Netlib LAPACKE_csytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    csytrs_obj->inforef = CSYTRF( csytrs_obj->matrix_layout,
                            csytrs_obj->uplo, csytrs_obj->n,
                                     csytrs_obj->aref,
                      csytrs_obj->lda, csytrs_obj->ipivref);

    csytrs_obj->inforef = CSYTRS( csytrs_obj->matrix_layout,
                                  csytrs_obj->uplo, csytrs_obj->n,
                                  csytrs_obj->nrhs,
                                  (const lapack_complex_float *)csytrs_obj->aref,
                                  csytrs_obj->lda, csytrs_obj->ipivref,
                                  csytrs_obj->bref, csytrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    csytrs_obj->info = LAPACKE_csytrf( csytrs_obj->matrix_layout,
                                 csytrs_obj->uplo, csytrs_obj->n,
                                     csytrs_obj->a,
                               csytrs_obj->lda, csytrs_obj->ipiv);

    csytrs_obj->info = LAPACKE_csytrs( csytrs_obj->matrix_layout,
                csytrs_obj->uplo, csytrs_obj->n, csytrs_obj->nrhs,
                                  (const lapack_complex_float *)csytrs_obj->a,
                               csytrs_obj->lda, csytrs_obj->ipiv,
                                 csytrs_obj->b, csytrs_obj->ldb );


    if( csytrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrs is wrong\n",
                    csytrs_obj->info );
    }
    if( csytrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrs is wrong\n",
        csytrs_obj->inforef );
    }
}

TEST_F(csytrs_test, csytrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_obj->b_bufsize,
                           csytrs_obj->b, csytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs_test, csytrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_obj->b_bufsize,
                           csytrs_obj->b, csytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs_test, csytrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_obj->b_bufsize,
                           csytrs_obj->b, csytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs_test, csytrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_obj->b_bufsize,
                           csytrs_obj->b, csytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sytrs_dcomplex_parameters  class definition */
class sytrs_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sytrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_dcomplex_parameters ();
};  /* end of sytrs_dcomplex_parameters  class definition */


/* Constructor sytrs_dcomplex_parameters definition */
sytrs_dcomplex_parameters:: sytrs_dcomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, lapack_int lda_i,
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sytrs DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, lda, ldb, nrhs);
#endif

    if(matrix_layout==LAPACK_COL_MAJOR){
        b_bufsize = ldb*nrhs;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        b_bufsize = ldb*n;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_dcomplex_parameters object: malloc error.";
       sytrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs_dcomplex_parameters:: ~sytrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrs_free();
}


//  Test fixture class definition
class zsytrs_test  : public  ::testing::Test {
public:
   sytrs_dcomplex_parameters  *zsytrs_obj;
   void SetUp();  
   void TearDown () { delete zsytrs_obj; }
};


void zsytrs_test::SetUp(){

    /* LAPACKE ZSYTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zsytrs ZSYTRS;

     /* LAPACKE ZSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zsytrf ZSYTRF;


    zsytrs_obj = new  sytrs_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zsytrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYTRS = (Fptr_NL_LAPACKE_zsytrs)dlsym(zsytrs_obj->hModule, "LAPACKE_zsytrs");
    ASSERT_TRUE(ZSYTRS != NULL) << "failed to get the Netlib LAPACKE_zsytrs symbol";

    ZSYTRF = (Fptr_NL_LAPACKE_zsytrf)dlsym(zsytrs_obj->hModule,"LAPACKE_zsytrf");
    ASSERT_TRUE(ZSYTRF != NULL) << "failed to get the Netlib LAPACKE_zsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zsytrs_obj->inforef = ZSYTRF( zsytrs_obj->matrix_layout,
                                    zsytrs_obj->uplo, zsytrs_obj->n,
                                     zsytrs_obj->aref,
                               zsytrs_obj->lda, zsytrs_obj->ipivref);

    zsytrs_obj->inforef = ZSYTRS( zsytrs_obj->matrix_layout,
                                  zsytrs_obj->uplo, zsytrs_obj->n,
                                  zsytrs_obj->nrhs,
                                  (const lapack_complex_double *)zsytrs_obj->aref,
                                  zsytrs_obj->lda, zsytrs_obj->ipivref,
                                  zsytrs_obj->bref, zsytrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zsytrs_obj->info = LAPACKE_zsytrf( zsytrs_obj->matrix_layout,
                                 zsytrs_obj->uplo, zsytrs_obj->n,
                                     zsytrs_obj->a,
                               zsytrs_obj->lda, zsytrs_obj->ipiv);

    zsytrs_obj->info = LAPACKE_zsytrs( zsytrs_obj->matrix_layout,
                zsytrs_obj->uplo, zsytrs_obj->n, zsytrs_obj->nrhs,
                                  (const lapack_complex_double *)zsytrs_obj->a,
                               zsytrs_obj->lda, zsytrs_obj->ipiv,
                                 zsytrs_obj->b, zsytrs_obj->ldb );


    if( zsytrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrs is wrong\n",
                    zsytrs_obj->info );
    }
    if( zsytrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrs is wrong\n",
        zsytrs_obj->inforef );
    }
}

TEST_F(zsytrs_test, zsytrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_obj->b_bufsize,
                           zsytrs_obj->b, zsytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsytrs_test, zsytrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_obj->b_bufsize,
                           zsytrs_obj->b, zsytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsytrs_test, zsytrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_obj->b_bufsize,
                           zsytrs_obj->b, zsytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsytrs_test, zsytrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_obj->b_bufsize,
                           zsytrs_obj->b, zsytrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
