#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define sytrs_rook_free() \
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


/* Begin sytrs_rook_double_parameters  class definition */
class sytrs_rook_double_parameters{
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
      sytrs_rook_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_rook_double_parameters ();
};  /* end of sytrs_rook_double_parameters  class definition */


/* Constructor sytrs_rook_double_parameters definition */
sytrs_rook_double_parameters:: sytrs_rook_double_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_rook Double:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sytrs_rook_double_parameters object: malloc error.";
       sytrs_rook_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs_rook_double_parameters:: ~sytrs_rook_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_rook_double_parameters object: destructor invoked. \n");
#endif
   sytrs_rook_free();
}


//  Test fixture class definition
class dsytrs_rook_test  : public  ::testing::Test {
public:
   sytrs_rook_double_parameters  *dsytrs_rook_obj;
   void SetUp();  
   void TearDown () { delete dsytrs_rook_obj; }
};


void dsytrs_rook_test::SetUp(){

    /* LAPACKE DSYTRS_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrs_rook) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs, const double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                                            double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dsytrs_rook DSYTRS_ROOK;

     /* LAPACKE DSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf) ( int matrix_layout, char uplo,
             lapack_int n, double* a, lapack_int lda, lapack_int* ipiv );

    Fptr_NL_LAPACKE_dsytrf DSYTRF;

    dsytrs_rook_obj = new  sytrs_rook_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dsytrs_rook_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytrs_rook_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytrs_rook_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytrs_rook_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYTRS_ROOK = (Fptr_NL_LAPACKE_dsytrs_rook)dlsym(dsytrs_rook_obj->hModule, "LAPACKE_dsytrs_rook");
    ASSERT_TRUE(DSYTRS_ROOK != NULL) << "failed to get the Netlib LAPACKE_dsytrs_rook symbol";

    DSYTRF = (Fptr_NL_LAPACKE_dsytrf)dlsym(dsytrs_rook_obj->hModule,"LAPACKE_dsytrf");
    ASSERT_TRUE(DSYTRF != NULL) << "failed to get the Netlib LAPACKE_dsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs_rook function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dsytrs_rook_obj->inforef = DSYTRF( dsytrs_rook_obj->matrix_layout,
                            dsytrs_rook_obj->uplo, dsytrs_rook_obj->n,
                                     dsytrs_rook_obj->aref,
                      dsytrs_rook_obj->lda, dsytrs_rook_obj->ipivref);

    dsytrs_rook_obj->inforef = DSYTRS_ROOK( dsytrs_rook_obj->matrix_layout,
                                  dsytrs_rook_obj->uplo, dsytrs_rook_obj->n,
                                  dsytrs_rook_obj->nrhs,
                                  (const double *)dsytrs_rook_obj->aref,
                              dsytrs_rook_obj->lda, dsytrs_rook_obj->ipivref,
                                dsytrs_rook_obj->bref, dsytrs_rook_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dsytrs_rook_obj->info = LAPACKE_dsytrf( dsytrs_rook_obj->matrix_layout,
                                 dsytrs_rook_obj->uplo, dsytrs_rook_obj->n,
                                     dsytrs_rook_obj->a,
                               dsytrs_rook_obj->lda, dsytrs_rook_obj->ipiv);

    dsytrs_rook_obj->info = LAPACKE_dsytrs_rook( dsytrs_rook_obj->matrix_layout,
                dsytrs_rook_obj->uplo, dsytrs_rook_obj->n, dsytrs_rook_obj->nrhs,
                                  (const double *)dsytrs_rook_obj->a,
                               dsytrs_rook_obj->lda, dsytrs_rook_obj->ipiv,
                                 dsytrs_rook_obj->b, dsytrs_rook_obj->ldb );


    if( dsytrs_rook_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrs_rook is wrong\n",
                    dsytrs_rook_obj->info );
    }
    if( dsytrs_rook_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrs_rook is wrong\n",
        dsytrs_rook_obj->inforef );
    }
}

TEST_F(dsytrs_rook_test, dsytrs_rook1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_rook_obj->b_bufsize,
                           dsytrs_rook_obj->b, dsytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsytrs_rook_test, dsytrs_rook2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_rook_obj->b_bufsize,
                           dsytrs_rook_obj->b, dsytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsytrs_rook_test, dsytrs_rook3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_rook_obj->b_bufsize,
                           dsytrs_rook_obj->b, dsytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsytrs_rook_test, dsytrs_rook4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_rook_obj->b_bufsize,
                           dsytrs_rook_obj->b, dsytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrs_rook_float_parameters  class definition */
class sytrs_rook_float_parameters{

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
      sytrs_rook_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~sytrs_rook_float_parameters ();
};  /* end of sytrs_rook_float_parameters  class definition */


/* Constructor sytrs_rook_float_parameters definition */
sytrs_rook_float_parameters:: sytrs_rook_float_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_rook float:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sytrs_rook_double_parameters object: malloc error.";
       sytrs_rook_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

sytrs_rook_float_parameters:: ~sytrs_rook_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_rook_float_parameters object: destructor invoked. \n");
#endif
   sytrs_rook_free();
}


//  Test fixture class definition
class ssytrs_rook_test  : public  ::testing::Test {
public:
   sytrs_rook_float_parameters  *ssytrs_rook_obj;
   void SetUp();  
   void TearDown () { delete ssytrs_rook_obj; }
};


void ssytrs_rook_test::SetUp(){

    /* LAPACKE SSYTRS_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrs_rook) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs, const float * a,
                                  lapack_int lda, const lapack_int * ipiv,
                                            float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_ssytrs_rook SSYTRS_ROOK;

     /* LAPACKE SSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf) ( int matrix_layout, char uplo,
                lapack_int n, float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_ssytrf SSYTRF;

    ssytrs_rook_obj = new  sytrs_rook_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    ssytrs_rook_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytrs_rook_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytrs_rook_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytrs_rook_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYTRS_ROOK = (Fptr_NL_LAPACKE_ssytrs_rook)dlsym(ssytrs_rook_obj->hModule, "LAPACKE_ssytrs_rook");
    ASSERT_TRUE(SSYTRS_ROOK != NULL) << "failed to get the Netlib LAPACKE_ssytrs_rook symbol";

    SSYTRF = (Fptr_NL_LAPACKE_ssytrf)dlsym(ssytrs_rook_obj->hModule,"LAPACKE_ssytrf");
    ASSERT_TRUE(SSYTRF != NULL) << "failed to get the Netlib LAPACKE_ssytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs_rook function */

    /* Compute the Netlib-Lapacke's reference o/p */
    ssytrs_rook_obj->inforef = SSYTRF( ssytrs_rook_obj->matrix_layout,
                                    ssytrs_rook_obj->uplo, ssytrs_rook_obj->n,
                                     ssytrs_rook_obj->aref,
                               ssytrs_rook_obj->lda, ssytrs_rook_obj->ipivref);

    ssytrs_rook_obj->inforef = SSYTRS_ROOK( ssytrs_rook_obj->matrix_layout,
                                  ssytrs_rook_obj->uplo, ssytrs_rook_obj->n,
                                  ssytrs_rook_obj->nrhs,
                                  (const float *)ssytrs_rook_obj->aref,
                                  ssytrs_rook_obj->lda, ssytrs_rook_obj->ipivref,
                                  ssytrs_rook_obj->bref, ssytrs_rook_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssytrs_rook_obj->info = LAPACKE_ssytrf( ssytrs_rook_obj->matrix_layout,
                                 ssytrs_rook_obj->uplo, ssytrs_rook_obj->n,
                                     ssytrs_rook_obj->a,
                               ssytrs_rook_obj->lda, ssytrs_rook_obj->ipiv);

    ssytrs_rook_obj->info = LAPACKE_ssytrs_rook( ssytrs_rook_obj->matrix_layout,
                ssytrs_rook_obj->uplo, ssytrs_rook_obj->n, ssytrs_rook_obj->nrhs,
                                  (const float *)ssytrs_rook_obj->a,
                               ssytrs_rook_obj->lda, ssytrs_rook_obj->ipiv,
                                 ssytrs_rook_obj->b, ssytrs_rook_obj->ldb );
    if( ssytrs_rook_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrs_rook is wrong\n",
                    ssytrs_rook_obj->info );
    }
    if( ssytrs_rook_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrs_rook is wrong\n",
        ssytrs_rook_obj->inforef );
    }
}

TEST_F(ssytrs_rook_test, ssytrs_rook1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs_rook_obj->b_bufsize,
                           ssytrs_rook_obj->b, ssytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sytrs_rook_scomplex_parameters  class definition */
class sytrs_rook_scomplex_parameters{
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
      sytrs_rook_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_rook_scomplex_parameters ();
};  /* end of sytrs_rook_scomplex_parameters  class definition */


/* Constructor sytrs_rook_scomplex_parameters definition */
sytrs_rook_scomplex_parameters:: sytrs_rook_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_rook scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sytrs_rook_scomplex_parameters object: malloc error.";
       sytrs_rook_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs_rook_scomplex_parameters:: ~sytrs_rook_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_rook_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrs_rook_free();
}


//  Test fixture class definition
class csytrs_rook_test  : public  ::testing::Test {
public:
   sytrs_rook_scomplex_parameters  *csytrs_rook_obj;
   void SetUp();  
   void TearDown () { delete csytrs_rook_obj; }
};


void csytrs_rook_test::SetUp(){

    /* LAPACKE CSYTRS_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrs_rook) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_float * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_csytrs_rook CSYTRS_ROOK;

     /* LAPACKE CSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf) ( int matrix_layout, char uplo,
    lapack_int n, lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_csytrf CSYTRF;


    csytrs_rook_obj = new  sytrs_rook_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    csytrs_rook_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytrs_rook_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytrs_rook_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytrs_rook_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYTRS_ROOK = (Fptr_NL_LAPACKE_csytrs_rook)dlsym(csytrs_rook_obj->hModule, "LAPACKE_csytrs_rook");
    ASSERT_TRUE(CSYTRS_ROOK != NULL) << "failed to get the Netlib LAPACKE_csytrs_rook symbol";

    CSYTRF = (Fptr_NL_LAPACKE_csytrf)dlsym(csytrs_rook_obj->hModule,"LAPACKE_csytrf");
    ASSERT_TRUE(CSYTRF != NULL) << "failed to get the Netlib LAPACKE_csytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs_rook function */

    /* Compute the Netlib-Lapacke's reference o/p */
    csytrs_rook_obj->inforef = CSYTRF( csytrs_rook_obj->matrix_layout,
                            csytrs_rook_obj->uplo, csytrs_rook_obj->n,
                                     csytrs_rook_obj->aref,
                      csytrs_rook_obj->lda, csytrs_rook_obj->ipivref);

    csytrs_rook_obj->inforef = CSYTRS_ROOK( csytrs_rook_obj->matrix_layout,
                                  csytrs_rook_obj->uplo, csytrs_rook_obj->n,
                                  csytrs_rook_obj->nrhs,
                                  (const lapack_complex_float *)csytrs_rook_obj->aref,
                                  csytrs_rook_obj->lda, csytrs_rook_obj->ipivref,
                                  csytrs_rook_obj->bref, csytrs_rook_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    csytrs_rook_obj->info = LAPACKE_csytrf( csytrs_rook_obj->matrix_layout,
                                 csytrs_rook_obj->uplo, csytrs_rook_obj->n,
                                     csytrs_rook_obj->a,
                               csytrs_rook_obj->lda, csytrs_rook_obj->ipiv);

    csytrs_rook_obj->info = LAPACKE_csytrs_rook( csytrs_rook_obj->matrix_layout,
                csytrs_rook_obj->uplo, csytrs_rook_obj->n, csytrs_rook_obj->nrhs,
                                  (const lapack_complex_float *)csytrs_rook_obj->a,
                               csytrs_rook_obj->lda, csytrs_rook_obj->ipiv,
                                 csytrs_rook_obj->b, csytrs_rook_obj->ldb );


    if( csytrs_rook_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrs_rook is wrong\n",
                    csytrs_rook_obj->info );
    }
    if( csytrs_rook_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrs_rook is wrong\n",
        csytrs_rook_obj->inforef );
    }
}

TEST_F(csytrs_rook_test, csytrs_rook1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_rook_obj->b_bufsize,
                           csytrs_rook_obj->b, csytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs_rook_test, csytrs_rook2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_rook_obj->b_bufsize,
                           csytrs_rook_obj->b, csytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs_rook_test, csytrs_rook3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_rook_obj->b_bufsize,
                           csytrs_rook_obj->b, csytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs_rook_test, csytrs_rook4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_rook_obj->b_bufsize,
                           csytrs_rook_obj->b, csytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sytrs_rook_dcomplex_parameters  class definition */
class sytrs_rook_dcomplex_parameters{
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
      sytrs_rook_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_rook_dcomplex_parameters ();
};  /* end of sytrs_rook_dcomplex_parameters  class definition */


/* Constructor sytrs_rook_dcomplex_parameters definition */
sytrs_rook_dcomplex_parameters:: sytrs_rook_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_rook DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sytrs_rook_dcomplex_parameters object: malloc error.";
       sytrs_rook_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs_rook_dcomplex_parameters:: ~sytrs_rook_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_rook_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrs_rook_free();
}


//  Test fixture class definition
class zsytrs_rook_test  : public  ::testing::Test {
public:
   sytrs_rook_dcomplex_parameters  *zsytrs_rook_obj;
   void SetUp();  
   void TearDown () { delete zsytrs_rook_obj; }
};


void zsytrs_rook_test::SetUp(){

    /* LAPACKE ZSYTRS_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrs_rook) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zsytrs_rook ZSYTRS_ROOK;

     /* LAPACKE ZSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zsytrf ZSYTRF;


    zsytrs_rook_obj = new  sytrs_rook_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zsytrs_rook_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytrs_rook_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytrs_rook_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytrs_rook_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYTRS_ROOK = (Fptr_NL_LAPACKE_zsytrs_rook)dlsym(zsytrs_rook_obj->hModule, "LAPACKE_zsytrs_rook");
    ASSERT_TRUE(ZSYTRS_ROOK != NULL) << "failed to get the Netlib LAPACKE_zsytrs_rook symbol";

    ZSYTRF = (Fptr_NL_LAPACKE_zsytrf)dlsym(zsytrs_rook_obj->hModule,"LAPACKE_zsytrf");
    ASSERT_TRUE(ZSYTRF != NULL) << "failed to get the Netlib LAPACKE_zsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs_rook function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zsytrs_rook_obj->inforef = ZSYTRF( zsytrs_rook_obj->matrix_layout,
                                    zsytrs_rook_obj->uplo, zsytrs_rook_obj->n,
                                     zsytrs_rook_obj->aref,
                               zsytrs_rook_obj->lda, zsytrs_rook_obj->ipivref);

    zsytrs_rook_obj->inforef = ZSYTRS_ROOK( zsytrs_rook_obj->matrix_layout,
                                  zsytrs_rook_obj->uplo, zsytrs_rook_obj->n,
                                  zsytrs_rook_obj->nrhs,
                                  (const lapack_complex_double *)zsytrs_rook_obj->aref,
                                  zsytrs_rook_obj->lda, zsytrs_rook_obj->ipivref,
                                  zsytrs_rook_obj->bref, zsytrs_rook_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zsytrs_rook_obj->info = LAPACKE_zsytrf( zsytrs_rook_obj->matrix_layout,
                                 zsytrs_rook_obj->uplo, zsytrs_rook_obj->n,
                                     zsytrs_rook_obj->a,
                               zsytrs_rook_obj->lda, zsytrs_rook_obj->ipiv);

    zsytrs_rook_obj->info = LAPACKE_zsytrs_rook( zsytrs_rook_obj->matrix_layout,
                zsytrs_rook_obj->uplo, zsytrs_rook_obj->n, zsytrs_rook_obj->nrhs,
                                  (const lapack_complex_double *)zsytrs_rook_obj->a,
                               zsytrs_rook_obj->lda, zsytrs_rook_obj->ipiv,
                                 zsytrs_rook_obj->b, zsytrs_rook_obj->ldb );


    if( zsytrs_rook_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrs_rook is wrong\n",
                    zsytrs_rook_obj->info );
    }
    if( zsytrs_rook_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrs_rook is wrong\n",
        zsytrs_rook_obj->inforef );
    }
}

TEST_F(zsytrs_rook_test, zsytrs_rook1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_rook_obj->b_bufsize,
                           zsytrs_rook_obj->b, zsytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsytrs_rook_test, zsytrs_rook2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_rook_obj->b_bufsize,
                           zsytrs_rook_obj->b, zsytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsytrs_rook_test, zsytrs_rook3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_rook_obj->b_bufsize,
                           zsytrs_rook_obj->b, zsytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsytrs_rook_test, zsytrs_rook4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_rook_obj->b_bufsize,
                           zsytrs_rook_obj->b, zsytrs_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
