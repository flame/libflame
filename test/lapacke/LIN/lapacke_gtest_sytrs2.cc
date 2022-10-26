#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define sytrs2_free() \
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


/* Begin sytrs2_double_parameters  class definition */
class sytrs2_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      sytrs2_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs2_double_parameters ();
};  /* end of sytrs2_double_parameters  class definition */


/* Constructor sytrs2_double_parameters definition */
sytrs2_double_parameters:: sytrs2_double_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs2 Double:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sytrs2_double_parameters object: malloc error.";
       sytrs2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs2_double_parameters:: ~sytrs2_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs2_double_parameters object: destructor invoked. \n");
#endif
   sytrs2_free();
}


//  Test fixture class definition
class dsytrs2_test  : public  ::testing::Test {
public:
   sytrs2_double_parameters  *dsytrs2_obj;
   void SetUp();  
   void TearDown () { delete dsytrs2_obj; }
};


void dsytrs2_test::SetUp(){

    /* LAPACKE DSYTRS2 prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrs2) ( int matrix_layout, char uplo,
                                lapack_int n, lapack_int nrhs, double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                                            double * b, lapack_int ldb );

    Fptr_NL_LAPACKE_dsytrs2 DSYTRS2;

     /* LAPACKE DSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf) ( int matrix_layout, char uplo,
             lapack_int n, double* a, lapack_int lda, lapack_int* ipiv );

    Fptr_NL_LAPACKE_dsytrf DSYTRF;

    dsytrs2_obj      = new  sytrs2_double_parameters(
                                lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dsytrs2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytrs2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytrs2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytrs2_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYTRS2 = (Fptr_NL_LAPACKE_dsytrs2)dlsym(dsytrs2_obj->hModule, "LAPACKE_dsytrs2");
    ASSERT_TRUE(DSYTRS2 != NULL) << "failed to get the Netlib LAPACKE_dsytrs2 symbol";

    DSYTRF = (Fptr_NL_LAPACKE_dsytrf)dlsym(dsytrs2_obj->hModule,"LAPACKE_dsytrf");
    ASSERT_TRUE(DSYTRF != NULL) << "failed to get the Netlib LAPACKE_dsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs2 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dsytrs2_obj->inforef = DSYTRF( dsytrs2_obj->matrix_layout,
                            dsytrs2_obj->uplo, dsytrs2_obj->n,
                                     dsytrs2_obj->aref,
                      dsytrs2_obj->lda, dsytrs2_obj->ipivref);

    dsytrs2_obj->inforef = DSYTRS2( dsytrs2_obj->matrix_layout,
                                  dsytrs2_obj->uplo, dsytrs2_obj->n,
                                  dsytrs2_obj->nrhs, dsytrs2_obj->aref,
                              dsytrs2_obj->lda, dsytrs2_obj->ipivref,
                                dsytrs2_obj->bref, dsytrs2_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dsytrs2_obj->info = LAPACKE_dsytrf( dsytrs2_obj->matrix_layout,
                 dsytrs2_obj->uplo, dsytrs2_obj->n, dsytrs2_obj->a,
                              dsytrs2_obj->lda, dsytrs2_obj->ipiv);

    dsytrs2_obj->info = LAPACKE_dsytrs2( dsytrs2_obj->matrix_layout,
               dsytrs2_obj->uplo, dsytrs2_obj->n, dsytrs2_obj->nrhs,
                dsytrs2_obj->a, dsytrs2_obj->lda, dsytrs2_obj->ipiv,
                                 dsytrs2_obj->b, dsytrs2_obj->ldb );


    if( dsytrs2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrs2 is wrong\n",
                    dsytrs2_obj->info );
    }
    if( dsytrs2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrs2 is wrong\n",
        dsytrs2_obj->inforef );
    }
}

TEST_F(dsytrs2_test, dsytrs21) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff = computeDiff_d( dsytrs2_obj->b_bufsize, dsytrs2_obj->b, dsytrs2_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytrs2_test, dsytrs22) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs2_obj->b_bufsize, dsytrs2_obj->b, dsytrs2_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytrs2_test, dsytrs23) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff = computeDiff_d( dsytrs2_obj->b_bufsize, dsytrs2_obj->b, dsytrs2_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytrs2_test, dsytrs24) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs2_obj->b_bufsize, dsytrs2_obj->b, dsytrs2_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sytrs2_float_parameters  class definition */
class sytrs2_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      sytrs2_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~sytrs2_float_parameters ();
};  /* end of sytrs2_float_parameters  class definition */


/* Constructor sytrs2_float_parameters definition */
sytrs2_float_parameters:: sytrs2_float_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs2 float:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sytrs2_double_parameters object: malloc error.";
       sytrs2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

sytrs2_float_parameters:: ~sytrs2_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs2_float_parameters object: destructor invoked. \n");
#endif
   sytrs2_free();
}


//  Test fixture class definition
class ssytrs2_test  : public  ::testing::Test {
public:
   sytrs2_float_parameters  *ssytrs2_obj;
   void SetUp();  
   void TearDown () { delete ssytrs2_obj; }
};


void ssytrs2_test::SetUp(){

    /* LAPACKE SSYTRS2 prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrs2) ( int matrix_layout, char uplo,
                                lapack_int n, lapack_int nrhs,  float * a,
                                  lapack_int lda, const lapack_int * ipiv,
                                            float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_ssytrs2 SSYTRS2;

     /* LAPACKE SSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf) ( int matrix_layout, char uplo,
                lapack_int n, float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_ssytrf SSYTRF;

    ssytrs2_obj = new  sytrs2_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    ssytrs2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytrs2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytrs2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytrs2_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYTRS2 = (Fptr_NL_LAPACKE_ssytrs2)dlsym(ssytrs2_obj->hModule, "LAPACKE_ssytrs2");
    ASSERT_TRUE(SSYTRS2 != NULL) << "failed to get the Netlib LAPACKE_ssytrs2 symbol";

    SSYTRF = (Fptr_NL_LAPACKE_ssytrf)dlsym(ssytrs2_obj->hModule,"LAPACKE_ssytrf");
    ASSERT_TRUE(SSYTRF != NULL) << "failed to get the Netlib LAPACKE_ssytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs2 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    ssytrs2_obj->inforef = SSYTRF( ssytrs2_obj->matrix_layout,
                                    ssytrs2_obj->uplo, 
                                    ssytrs2_obj->n,
                                    ssytrs2_obj->aref,
                                    ssytrs2_obj->lda,
                                    ssytrs2_obj->ipivref);

    ssytrs2_obj->inforef = SSYTRS2( ssytrs2_obj->matrix_layout,
                                  ssytrs2_obj->uplo, ssytrs2_obj->n,
                               ssytrs2_obj->nrhs, ssytrs2_obj->aref,
                             ssytrs2_obj->lda, ssytrs2_obj->ipivref,
                               ssytrs2_obj->bref, ssytrs2_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssytrs2_obj->info = LAPACKE_ssytrf( ssytrs2_obj->matrix_layout,
                                 ssytrs2_obj->uplo, ssytrs2_obj->n,
                  ssytrs2_obj->a, ssytrs2_obj->lda, ssytrs2_obj->ipiv);

    ssytrs2_obj->info = LAPACKE_ssytrs2( ssytrs2_obj->matrix_layout,
                ssytrs2_obj->uplo, ssytrs2_obj->n, ssytrs2_obj->nrhs,
                ssytrs2_obj->a, ssytrs2_obj->lda, ssytrs2_obj->ipiv,
                                 ssytrs2_obj->b, ssytrs2_obj->ldb );

    if( ssytrs2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrs2 is wrong\n",
                    ssytrs2_obj->info );
    }
    if( ssytrs2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrs2 is wrong\n",
        ssytrs2_obj->inforef );
    }
}

TEST_F(ssytrs2_test, ssytrs21) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs2_obj->b_bufsize, ssytrs2_obj->b, ssytrs2_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytrs2_test, ssytrs22) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs2_obj->b_bufsize, ssytrs2_obj->b, ssytrs2_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytrs2_test, ssytrs23) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs2_obj->b_bufsize, ssytrs2_obj->b, ssytrs2_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytrs2_test, ssytrs24) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs2_obj->b_bufsize, ssytrs2_obj->b, ssytrs2_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sytrs2_scomplex_parameters  class definition */
class sytrs2_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      sytrs2_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs2_scomplex_parameters ();
};  /* end of sytrs2_scomplex_parameters  class definition */


/* Constructor sytrs2_scomplex_parameters definition */
sytrs2_scomplex_parameters:: sytrs2_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs2 scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sytrs2_scomplex_parameters object: malloc error.";
       sytrs2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs2_scomplex_parameters:: ~sytrs2_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs2_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrs2_free();
}


//  Test fixture class definition
class csytrs2_test  : public  ::testing::Test {
public:
   sytrs2_scomplex_parameters  *csytrs2_obj;
   void SetUp();  
   void TearDown () { delete csytrs2_obj; }
};


void csytrs2_test::SetUp(){

    /* LAPACKE CSYTRS2 prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrs2)( int matrix_layout, char uplo,
                                          lapack_int n, lapack_int nrhs,
                                               lapack_complex_float * a,
                                lapack_int lda, const lapack_int * ipiv,
                             lapack_complex_float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_csytrs2 CSYTRS2;

     /* LAPACKE CSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf) ( int matrix_layout, char uplo,
    lapack_int n, lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_csytrf CSYTRF;


    csytrs2_obj = new  sytrs2_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    csytrs2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytrs2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytrs2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytrs2_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYTRS2 = (Fptr_NL_LAPACKE_csytrs2)dlsym(csytrs2_obj->hModule, "LAPACKE_csytrs2");
    ASSERT_TRUE(CSYTRS2 != NULL) << "failed to get the Netlib LAPACKE_csytrs2 symbol";

    CSYTRF = (Fptr_NL_LAPACKE_csytrf)dlsym(csytrs2_obj->hModule,"LAPACKE_csytrf");
    ASSERT_TRUE(CSYTRF != NULL) << "failed to get the Netlib LAPACKE_csytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs2 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    csytrs2_obj->inforef = CSYTRF( csytrs2_obj->matrix_layout,
                            csytrs2_obj->uplo, csytrs2_obj->n,
                                     csytrs2_obj->aref,
                      csytrs2_obj->lda, csytrs2_obj->ipivref);

    csytrs2_obj->inforef = CSYTRS2( csytrs2_obj->matrix_layout,
                                  csytrs2_obj->uplo, csytrs2_obj->n,
                                  csytrs2_obj->nrhs,
                                  csytrs2_obj->aref,
                                  csytrs2_obj->lda, csytrs2_obj->ipivref,
                                  csytrs2_obj->bref, csytrs2_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    csytrs2_obj->info = LAPACKE_csytrf( csytrs2_obj->matrix_layout,
                                 csytrs2_obj->uplo, csytrs2_obj->n,
                                     csytrs2_obj->a,
                               csytrs2_obj->lda, csytrs2_obj->ipiv);

    csytrs2_obj->info = LAPACKE_csytrs2( csytrs2_obj->matrix_layout,
                csytrs2_obj->uplo, csytrs2_obj->n, csytrs2_obj->nrhs,
                                  csytrs2_obj->a,
                               csytrs2_obj->lda, csytrs2_obj->ipiv,
                                 csytrs2_obj->b, csytrs2_obj->ldb );


    if( csytrs2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrs2 is wrong\n",
                    csytrs2_obj->info );
    }
    if( csytrs2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrs2 is wrong\n",
        csytrs2_obj->inforef );
    }
}

TEST_F(csytrs2_test, csytrs21) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs2_obj->b_bufsize,
                           csytrs2_obj->b, csytrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs2_test, csytrs22) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs2_obj->b_bufsize,
                           csytrs2_obj->b, csytrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs2_test, csytrs23) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs2_obj->b_bufsize,
                           csytrs2_obj->b, csytrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs2_test, csytrs24) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs2_obj->b_bufsize,
                           csytrs2_obj->b, csytrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sytrs2_dcomplex_parameters  class definition */
class sytrs2_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
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
      sytrs2_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs2_dcomplex_parameters ();
};  /* end of sytrs2_dcomplex_parameters  class definition */


/* Constructor sytrs2_dcomplex_parameters definition */
sytrs2_dcomplex_parameters:: sytrs2_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs2 DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sytrs2_dcomplex_parameters object: malloc error.";
       sytrs2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs2_dcomplex_parameters:: ~sytrs2_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs2_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrs2_free();
}


//  Test fixture class definition
class zsytrs2_test  : public  ::testing::Test {
public:
   sytrs2_dcomplex_parameters  *zsytrs2_obj;
   void SetUp();  
   void TearDown () { delete zsytrs2_obj; }
};


void zsytrs2_test::SetUp(){

    /* LAPACKE ZSYTRS2 prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrs2) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          lapack_complex_double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zsytrs2 ZSYTRS2;

     /* LAPACKE ZSYTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zsytrf ZSYTRF;


    zsytrs2_obj = new  sytrs2_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zsytrs2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytrs2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytrs2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytrs2_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYTRS2 = (Fptr_NL_LAPACKE_zsytrs2)dlsym(zsytrs2_obj->hModule, "LAPACKE_zsytrs2");
    ASSERT_TRUE(ZSYTRS2 != NULL) << "failed to get the Netlib LAPACKE_zsytrs2 symbol";

    ZSYTRF = (Fptr_NL_LAPACKE_zsytrf)dlsym(zsytrs2_obj->hModule,"LAPACKE_zsytrf");
    ASSERT_TRUE(ZSYTRF != NULL) << "failed to get the Netlib LAPACKE_zsytrf symbol";

    /* Pre condition: need to call sytrf - before calling sytrs2 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zsytrs2_obj->inforef = ZSYTRF( zsytrs2_obj->matrix_layout,
                                    zsytrs2_obj->uplo, zsytrs2_obj->n,
                                     zsytrs2_obj->aref,
                               zsytrs2_obj->lda, zsytrs2_obj->ipivref);

    zsytrs2_obj->inforef = ZSYTRS2( zsytrs2_obj->matrix_layout,
                                  zsytrs2_obj->uplo, zsytrs2_obj->n,
                                  zsytrs2_obj->nrhs,
                                  zsytrs2_obj->aref,
                                  zsytrs2_obj->lda, zsytrs2_obj->ipivref,
                                  zsytrs2_obj->bref, zsytrs2_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zsytrs2_obj->info = LAPACKE_zsytrf( zsytrs2_obj->matrix_layout,
                                 zsytrs2_obj->uplo, zsytrs2_obj->n,
                                     zsytrs2_obj->a,
                               zsytrs2_obj->lda, zsytrs2_obj->ipiv);

    zsytrs2_obj->info = LAPACKE_zsytrs2( zsytrs2_obj->matrix_layout,
                zsytrs2_obj->uplo, zsytrs2_obj->n, zsytrs2_obj->nrhs,
                                  zsytrs2_obj->a,
                               zsytrs2_obj->lda, zsytrs2_obj->ipiv,
                                 zsytrs2_obj->b, zsytrs2_obj->ldb );


    if( zsytrs2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrs2 is wrong\n",
                    zsytrs2_obj->info );
    }
    if( zsytrs2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrs2 is wrong\n",
        zsytrs2_obj->inforef );
    }
}

TEST_F(zsytrs2_test, zsytrs21) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs2_obj->b_bufsize,
                           zsytrs2_obj->b, zsytrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsytrs2_test, zsytrs22) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs2_obj->b_bufsize,
                           zsytrs2_obj->b, zsytrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsytrs2_test, zsytrs23) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs2_obj->b_bufsize,
                           zsytrs2_obj->b, zsytrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsytrs2_test, zsytrs24) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs2_obj->b_bufsize,
                           zsytrs2_obj->b, zsytrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
