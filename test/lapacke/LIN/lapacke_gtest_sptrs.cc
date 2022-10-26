#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define LAPACKE_TEST_VERBOSE (1)
#define sptrs_free() \
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


/* Begin sptrs_double_parameters  class definition */
class sptrs_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
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
      sptrs_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sptrs_double_parameters ();
};  /* end of sptrs_double_parameters  class definition */


/* Constructor sptrs_double_parameters definition */
sptrs_double_parameters:: sptrs_double_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sptrs Double:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
             n, uplo, ldb, nrhs);
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
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*(n+1)/2));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sptrs_double_parameters object: malloc error.";
       sptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sptrs_double_parameters:: ~sptrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptrs_double_parameters object: destructor invoked. \n");
#endif
   sptrs_free();
}


//  Test fixture class definition
class dsptrs_test  : public  ::testing::Test {
public:
   sptrs_double_parameters  *dsptrs_obj;
   void SetUp();  
   void TearDown () { delete dsptrs_obj; }
};


void dsptrs_test::SetUp(){

    /* LAPACKE DSPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dsptrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs, const double *a,
                                   const lapack_int *ipiv,
                                            double *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dsptrs DSPTRS;

     /* LAPACKE DSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dsptrf) ( int matrix_layout, char uplo,
             lapack_int n, double* a,  lapack_int* ipiv );

    Fptr_NL_LAPACKE_dsptrf DSPTRF;

    dsptrs_obj = new  sptrs_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dsptrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsptrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsptrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsptrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSPTRS = (Fptr_NL_LAPACKE_dsptrs)dlsym(dsptrs_obj->hModule, "LAPACKE_dsptrs");
    ASSERT_TRUE(DSPTRS != NULL) << "failed to get the Netlib LAPACKE_dsptrs symbol";

    DSPTRF = (Fptr_NL_LAPACKE_dsptrf)dlsym(dsptrs_obj->hModule,"LAPACKE_dsptrf");
    ASSERT_TRUE(DSPTRF != NULL) << "failed to get the Netlib LAPACKE_dsptrf symbol";

    /* Pre condition: need to call sptrf - before calling sptrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dsptrs_obj->inforef = DSPTRF( dsptrs_obj->matrix_layout,
                            dsptrs_obj->uplo, dsptrs_obj->n,
                                     dsptrs_obj->aref,
                      dsptrs_obj->ipivref);

    dsptrs_obj->inforef = DSPTRS( dsptrs_obj->matrix_layout,
                                  dsptrs_obj->uplo, dsptrs_obj->n,
                                  dsptrs_obj->nrhs,
                                  (const double *)dsptrs_obj->aref,
                              dsptrs_obj->ipivref,
                                dsptrs_obj->bref, dsptrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dsptrs_obj->info = LAPACKE_dsptrf( dsptrs_obj->matrix_layout,
                                 dsptrs_obj->uplo, dsptrs_obj->n,
                                     dsptrs_obj->a,
                               dsptrs_obj->ipiv);

    dsptrs_obj->info = LAPACKE_dsptrs( dsptrs_obj->matrix_layout,
                dsptrs_obj->uplo, dsptrs_obj->n, dsptrs_obj->nrhs,
                                  (const double *)dsptrs_obj->a,
                               dsptrs_obj->ipiv,
                                 dsptrs_obj->b, dsptrs_obj->ldb );


    if( dsptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsptrs is wrong\n",
                    dsptrs_obj->info );
    }
    if( dsptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsptrs is wrong\n",
        dsptrs_obj->inforef );
    }
}

TEST_F(dsptrs_test, dsptrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsptrs_obj->b_bufsize,
                           dsptrs_obj->b, dsptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsptrs_test, dsptrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsptrs_obj->b_bufsize,
                           dsptrs_obj->b, dsptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsptrs_test, dsptrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsptrs_obj->b_bufsize,
                           dsptrs_obj->b, dsptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsptrs_test, dsptrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsptrs_obj->b_bufsize,
                           dsptrs_obj->b, dsptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sptrs_float_parameters  class definition */
class sptrs_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
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
      sptrs_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~sptrs_float_parameters ();
};  /* end of sptrs_float_parameters  class definition */


/* Constructor sptrs_float_parameters definition */
sptrs_float_parameters:: sptrs_float_parameters ( int matrix_layout_i,
                       char uplo_i, lapack_int n_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sptrs float:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
             n, uplo, ldb, nrhs);
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
    lapacke_gtest_alloc_float_buffer_pair(  &a, &aref, (n*(n+1)/2));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sptrs_double_parameters object: malloc error.";
       sptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

sptrs_float_parameters:: ~sptrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptrs_float_parameters object: destructor invoked. \n");
#endif
   sptrs_free();
}


//  Test fixture class definition
class ssptrs_test  : public  ::testing::Test {
public:
   sptrs_float_parameters  *ssptrs_obj;
   void SetUp();  
   void TearDown () { delete ssptrs_obj; }
};


void ssptrs_test::SetUp(){

    /* LAPACKE SSPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ssptrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs, const float *a,
                                   const lapack_int *ipiv,
                                            float *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_ssptrs SSPTRS;

     /* LAPACKE SSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_ssptrf) ( int matrix_layout, char uplo,
                lapack_int n, float* a,lapack_int* ipiv );

    Fptr_NL_LAPACKE_ssptrf SSPTRF;

    ssptrs_obj = new  sptrs_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    ssptrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssptrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssptrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssptrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSPTRS = (Fptr_NL_LAPACKE_ssptrs)dlsym(ssptrs_obj->hModule, "LAPACKE_ssptrs");
    ASSERT_TRUE(SSPTRS != NULL) << "failed to get the Netlib LAPACKE_ssptrs symbol";

    SSPTRF = (Fptr_NL_LAPACKE_ssptrf)dlsym(ssptrs_obj->hModule,"LAPACKE_ssptrf");
    ASSERT_TRUE(SSPTRF != NULL) << "failed to get the Netlib LAPACKE_ssptrf symbol";

    /* Pre condition: need to call sptrf - before calling sptrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    ssptrs_obj->inforef = SSPTRF( ssptrs_obj->matrix_layout,
                                    ssptrs_obj->uplo, ssptrs_obj->n,
                                     ssptrs_obj->aref,
                                ssptrs_obj->ipivref);

    ssptrs_obj->inforef = SSPTRS( ssptrs_obj->matrix_layout,
                                  ssptrs_obj->uplo, ssptrs_obj->n,
                                  ssptrs_obj->nrhs,
                                  (const float *)ssptrs_obj->aref,
                                   ssptrs_obj->ipivref,
                                  ssptrs_obj->bref, ssptrs_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssptrs_obj->info = LAPACKE_ssptrf( ssptrs_obj->matrix_layout,
                                 ssptrs_obj->uplo, ssptrs_obj->n,
                                     ssptrs_obj->a,
                                ssptrs_obj->ipiv);

    ssptrs_obj->info = LAPACKE_ssptrs( ssptrs_obj->matrix_layout,
                ssptrs_obj->uplo, ssptrs_obj->n, ssptrs_obj->nrhs,
                                  (const float *)ssptrs_obj->a,
                                ssptrs_obj->ipiv,
                                 ssptrs_obj->b, ssptrs_obj->ldb );
    if( ssptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssptrs is wrong\n",
                    ssptrs_obj->info );
    }
    if( ssptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssptrs is wrong\n",
        ssptrs_obj->inforef );
    }
}

TEST_F(ssptrs_test, ssptrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssptrs_obj->b_bufsize,
                           ssptrs_obj->b,
						   ssptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssptrs_test, ssptrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssptrs_obj->b_bufsize,
                           ssptrs_obj->b,
						   ssptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssptrs_test, ssptrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssptrs_obj->b_bufsize,
                           ssptrs_obj->b,
						   ssptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssptrs_test, ssptrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssptrs_obj->b_bufsize,
                           ssptrs_obj->b,
						   ssptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sptrs_scomplex_parameters  class definition */
class sptrs_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
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
      sptrs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sptrs_scomplex_parameters ();
};  /* end of sptrs_scomplex_parameters  class definition */


/* Constructor sptrs_scomplex_parameters definition */
sptrs_scomplex_parameters:: sptrs_scomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sptrs scomplex:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
             n, uplo, ldb, nrhs);
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
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*(n+1)/2));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sptrs_scomplex_parameters object: malloc error.";
       sptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sptrs_scomplex_parameters:: ~sptrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptrs_scomplex_parameters object: destructor invoked. \n");
#endif
   sptrs_free();
}


//  Test fixture class definition
class csptrs_test  : public  ::testing::Test {
public:
   sptrs_scomplex_parameters  *csptrs_obj;
   void SetUp();  
   void TearDown () { delete csptrs_obj; }
};


void csptrs_test::SetUp(){

    /* LAPACKE CSPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_csptrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_float *a,
                                   const lapack_int *ipiv,
                              lapack_complex_float *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_csptrs CSPTRS;

     /* LAPACKE CSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_csptrf) ( int matrix_layout, char uplo,
    lapack_int n, lapack_complex_float* a,lapack_int* ipiv );

    Fptr_NL_LAPACKE_csptrf CSPTRF;


    csptrs_obj = new  sptrs_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    csptrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csptrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csptrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csptrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSPTRS = (Fptr_NL_LAPACKE_csptrs)dlsym(csptrs_obj->hModule, "LAPACKE_csptrs");
    ASSERT_TRUE(CSPTRS != NULL) << "failed to get the Netlib LAPACKE_csptrs symbol";

    CSPTRF = (Fptr_NL_LAPACKE_csptrf)dlsym(csptrs_obj->hModule,"LAPACKE_csptrf");
    ASSERT_TRUE(CSPTRF != NULL) << "failed to get the Netlib LAPACKE_csptrf symbol";

    /* Pre condition: need to call sptrf - before calling sptrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    csptrs_obj->inforef = CSPTRF( csptrs_obj->matrix_layout,
                            csptrs_obj->uplo, csptrs_obj->n,
                                     csptrs_obj->aref,
                       csptrs_obj->ipivref);

    csptrs_obj->inforef = CSPTRS( csptrs_obj->matrix_layout,
                                  csptrs_obj->uplo, csptrs_obj->n,
                                  csptrs_obj->nrhs,
                      (const lapack_complex_float *)csptrs_obj->aref,
                                  csptrs_obj->ipivref,
                                  csptrs_obj->bref, csptrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    csptrs_obj->info = LAPACKE_csptrf( csptrs_obj->matrix_layout,
                                 csptrs_obj->uplo, csptrs_obj->n,
                                 csptrs_obj->a, csptrs_obj->ipiv);

    csptrs_obj->info = LAPACKE_csptrs( csptrs_obj->matrix_layout,
                csptrs_obj->uplo, csptrs_obj->n, csptrs_obj->nrhs,
                      (const lapack_complex_float *)csptrs_obj->a,
                csptrs_obj->ipiv, csptrs_obj->b, csptrs_obj->ldb );


    if( csptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csptrs is wrong\n",
                    csptrs_obj->info );
    }
    if( csptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csptrs is wrong\n",
        csptrs_obj->inforef );
    }
}

TEST_F(csptrs_test, csptrs1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csptrs_obj->b_bufsize,
                           csptrs_obj->b, 
						   csptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csptrs_test, csptrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csptrs_obj->b_bufsize,
                           csptrs_obj->b, 
						   csptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csptrs_test, csptrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csptrs_obj->b_bufsize,
                           csptrs_obj->b, csptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csptrs_test, csptrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csptrs_obj->b_bufsize,
                           csptrs_obj->b, csptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sptrs_dcomplex_parameters  class definition */
class sptrs_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
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
      sptrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sptrs_dcomplex_parameters ();
};  /* end of sptrs_dcomplex_parameters  class definition */


/* Constructor sptrs_dcomplex_parameters definition */
sptrs_dcomplex_parameters:: sptrs_dcomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n sptrs DComplex:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
             n, uplo, ldb, nrhs);
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
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*(n+1)/2));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sptrs_dcomplex_parameters object: malloc error.";
       sptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sptrs_dcomplex_parameters:: ~sptrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sptrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   sptrs_free();
}


//  Test fixture class definition
class zsptrs_test  : public  ::testing::Test {
public:
   sptrs_dcomplex_parameters  *zsptrs_obj;
   void SetUp();  
   void TearDown () { delete zsptrs_obj; }
};


void zsptrs_test::SetUp(){

    /* LAPACKE ZSPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zsptrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_double *a,
                                   const lapack_int *ipiv,
                              lapack_complex_double *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zsptrs ZSPTRS;

     /* LAPACKE ZSPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zsptrf) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_double *a,lapack_int *ipiv );

    Fptr_NL_LAPACKE_zsptrf ZSPTRF;


    zsptrs_obj = new  sptrs_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zsptrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsptrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsptrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsptrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSPTRS = (Fptr_NL_LAPACKE_zsptrs)dlsym(zsptrs_obj->hModule, "LAPACKE_zsptrs");
    ASSERT_TRUE(ZSPTRS != NULL) << "failed to get the Netlib LAPACKE_zsptrs symbol";

    ZSPTRF = (Fptr_NL_LAPACKE_zsptrf)dlsym(zsptrs_obj->hModule,"LAPACKE_zsptrf");
    ASSERT_TRUE(ZSPTRF != NULL) << "failed to get the Netlib LAPACKE_zsptrf symbol";

    /* Pre condition: need to call sptrf - before calling sptrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zsptrs_obj->inforef = ZSPTRF( zsptrs_obj->matrix_layout,
                                    zsptrs_obj->uplo, zsptrs_obj->n,
                                     zsptrs_obj->aref,
                               zsptrs_obj->ipivref);

    zsptrs_obj->inforef = ZSPTRS( zsptrs_obj->matrix_layout,
                                  zsptrs_obj->uplo, zsptrs_obj->n,
                                  zsptrs_obj->nrhs,
                                  (const lapack_complex_double *)zsptrs_obj->aref,
                                  zsptrs_obj->ipivref,
                                  zsptrs_obj->bref, zsptrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zsptrs_obj->info = LAPACKE_zsptrf( zsptrs_obj->matrix_layout,
                                 zsptrs_obj->uplo, zsptrs_obj->n,
                                     zsptrs_obj->a,
                               zsptrs_obj->ipiv);

    zsptrs_obj->info = LAPACKE_zsptrs( zsptrs_obj->matrix_layout,
                zsptrs_obj->uplo, zsptrs_obj->n, zsptrs_obj->nrhs,
                                  (const lapack_complex_double *)zsptrs_obj->a,
                               zsptrs_obj->ipiv,
                                 zsptrs_obj->b, zsptrs_obj->ldb );


    if( zsptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsptrs is wrong\n",
                    zsptrs_obj->info );
    }
    if( zsptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsptrs is wrong\n",
        zsptrs_obj->inforef );
    }
}

TEST_F(zsptrs_test, zsptrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsptrs_obj->b_bufsize,
                           zsptrs_obj->b, zsptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsptrs_test, zsptrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsptrs_obj->b_bufsize,
                           zsptrs_obj->b, zsptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsptrs_test, zsptrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsptrs_obj->b_bufsize,
                           zsptrs_obj->b, zsptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsptrs_test, zsptrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsptrs_obj->b_bufsize,
                           zsptrs_obj->b, zsptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
