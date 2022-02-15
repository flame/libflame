#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define sysv_rook_free() \
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


/* Begin sysv_rook_double_parameters  class definition */
class sysv_rook_double_parameters{
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
      sysv_rook_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_rook_double_parameters ();
};  /* end of sysv_rook_double_parameters  class definition */


/* Constructor sysv_rook_double_parameters definition */
sysv_rook_double_parameters:: sysv_rook_double_parameters ( int matrix_layout_i,
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
   printf(" \n sysv_rook Double:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sysv_rook_double_parameters object: malloc error.";
       sysv_rook_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_rook_double_parameters:: ~sysv_rook_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_rook_double_parameters object: destructor invoked. \n");
#endif
   sysv_rook_free();
}


//  Test fixture class definition
class dsysv_rook_test  : public  ::testing::Test {
public:
   sysv_rook_double_parameters  *dsysv_rook_obj;
   void SetUp();  
   void TearDown () { delete dsysv_rook_obj; }
};


void dsysv_rook_test::SetUp(){

    /* LAPACKE DSYSV_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_dsysv_rook) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,  double * a,
                                  lapack_int lda,  lapack_int * ipiv,
                                            double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dsysv_rook DSYSV_ROOK;

    dsysv_rook_obj = new  sysv_rook_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dsysv_rook_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsysv_rook_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsysv_rook_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsysv_rook_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYSV_ROOK = (Fptr_NL_LAPACKE_dsysv_rook)dlsym(dsysv_rook_obj->hModule, "LAPACKE_dsysv_rook");
    ASSERT_TRUE(DSYSV_ROOK != NULL) << "failed to get the Netlib LAPACKE_dsysv_rook symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    dsysv_rook_obj->inforef = DSYSV_ROOK( dsysv_rook_obj->matrix_layout,
                                  dsysv_rook_obj->uplo, dsysv_rook_obj->n,
                                  dsysv_rook_obj->nrhs,
                                  ( double *)dsysv_rook_obj->aref,
                              dsysv_rook_obj->lda, dsysv_rook_obj->ipivref,
                                dsysv_rook_obj->bref, dsysv_rook_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dsysv_rook_obj->info = LAPACKE_dsysv_rook( dsysv_rook_obj->matrix_layout,
                dsysv_rook_obj->uplo, dsysv_rook_obj->n, dsysv_rook_obj->nrhs,
                                  ( double *)dsysv_rook_obj->a,
                               dsysv_rook_obj->lda, dsysv_rook_obj->ipiv,
                                 dsysv_rook_obj->b, dsysv_rook_obj->ldb );


    if( dsysv_rook_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsysv_rook is wrong\n",
                    dsysv_rook_obj->info );
    }
    if( dsysv_rook_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsysv_rook is wrong\n",
        dsysv_rook_obj->inforef );
    }
}

TEST_F(dsysv_rook_test, dsysv_rook1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_rook_obj->b_bufsize,
                           dsysv_rook_obj->b, dsysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsysv_rook_test, dsysv_rook2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_rook_obj->b_bufsize,
                           dsysv_rook_obj->b, dsysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsysv_rook_test, dsysv_rook3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_rook_obj->b_bufsize,
                           dsysv_rook_obj->b, dsysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsysv_rook_test, dsysv_rook4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_rook_obj->b_bufsize,
                           dsysv_rook_obj->b, dsysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sysv_rook_float_parameters  class definition */
class sysv_rook_float_parameters{

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
      sysv_rook_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~sysv_rook_float_parameters ();
};  /* end of sysv_rook_float_parameters  class definition */


/* Constructor sysv_rook_float_parameters definition */
sysv_rook_float_parameters:: sysv_rook_float_parameters ( int matrix_layout_i,
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
   printf(" \n sysv_rook float:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sysv_rook_double_parameters object: malloc error.";
       sysv_rook_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

sysv_rook_float_parameters:: ~sysv_rook_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_rook_float_parameters object: destructor invoked. \n");
#endif
   sysv_rook_free();
}


//  Test fixture class definition
class ssysv_rook_test  : public  ::testing::Test {
public:
   sysv_rook_float_parameters  *ssysv_rook_obj;
   void SetUp();  
   void TearDown () { delete ssysv_rook_obj; }
};


void ssysv_rook_test::SetUp(){

    /* LAPACKE SSYSV_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_ssysv_rook) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,  float * a,
                                  lapack_int lda,  lapack_int * ipiv,
                                            float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_ssysv_rook SSYSV_ROOK;

    ssysv_rook_obj = new  sysv_rook_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    ssysv_rook_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssysv_rook_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssysv_rook_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssysv_rook_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYSV_ROOK = (Fptr_NL_LAPACKE_ssysv_rook)dlsym(ssysv_rook_obj->hModule, "LAPACKE_ssysv_rook");
    ASSERT_TRUE(SSYSV_ROOK != NULL) << "failed to get the Netlib LAPACKE_ssysv_rook symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    ssysv_rook_obj->inforef = SSYSV_ROOK( ssysv_rook_obj->matrix_layout,
                                  ssysv_rook_obj->uplo, ssysv_rook_obj->n,
                                  ssysv_rook_obj->nrhs,
                                  ( float *)ssysv_rook_obj->aref,
                                  ssysv_rook_obj->lda, ssysv_rook_obj->ipivref,
                                  ssysv_rook_obj->bref, ssysv_rook_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssysv_rook_obj->info = LAPACKE_ssysv_rook( ssysv_rook_obj->matrix_layout,
                ssysv_rook_obj->uplo, ssysv_rook_obj->n, ssysv_rook_obj->nrhs,
                                  ( float *)ssysv_rook_obj->a,
                               ssysv_rook_obj->lda, ssysv_rook_obj->ipiv,
                                 ssysv_rook_obj->b, ssysv_rook_obj->ldb );
    if( ssysv_rook_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssysv_rook is wrong\n",
                    ssysv_rook_obj->info );
    }
    if( ssysv_rook_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssysv_rook is wrong\n",
        ssysv_rook_obj->inforef );
    }
}

TEST_F(ssysv_rook_test, ssysv_rook1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_rook_obj->b_bufsize,
                           ssysv_rook_obj->b, ssysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_rook_test, ssysv_rook2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_rook_obj->b_bufsize,
                           ssysv_rook_obj->b, ssysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_rook_test, ssysv_rook3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_rook_obj->b_bufsize,
                           ssysv_rook_obj->b, ssysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_rook_test, ssysv_rook4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_rook_obj->b_bufsize,
                           ssysv_rook_obj->b, ssysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sysv_rook_scomplex_parameters  class definition */
class sysv_rook_scomplex_parameters{
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
      sysv_rook_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_rook_scomplex_parameters ();
};  /* end of sysv_rook_scomplex_parameters  class definition */


/* Constructor sysv_rook_scomplex_parameters definition */
sysv_rook_scomplex_parameters:: sysv_rook_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sysv_rook scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sysv_rook_scomplex_parameters object: malloc error.";
       sysv_rook_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_rook_scomplex_parameters:: ~sysv_rook_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_rook_scomplex_parameters object: destructor invoked. \n");
#endif
   sysv_rook_free();
}


//  Test fixture class definition
class csysv_rook_test  : public  ::testing::Test {
public:
   sysv_rook_scomplex_parameters  *csysv_rook_obj;
   void SetUp();  
   void TearDown () { delete csysv_rook_obj; }
};


void csysv_rook_test::SetUp(){

    /* LAPACKE CSYSV_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_csysv_rook) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_float * a,
                                  lapack_int lda,  lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_csysv_rook CSYSV_ROOK;

    csysv_rook_obj = new  sysv_rook_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    csysv_rook_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csysv_rook_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csysv_rook_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csysv_rook_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYSV_ROOK = (Fptr_NL_LAPACKE_csysv_rook)dlsym(csysv_rook_obj->hModule, "LAPACKE_csysv_rook");
    ASSERT_TRUE(CSYSV_ROOK != NULL) << "failed to get the Netlib LAPACKE_csysv_rook symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    csysv_rook_obj->inforef = CSYSV_ROOK( csysv_rook_obj->matrix_layout,
                                  csysv_rook_obj->uplo, csysv_rook_obj->n,
                                  csysv_rook_obj->nrhs,
                                  ( lapack_complex_float *)csysv_rook_obj->aref,
                                  csysv_rook_obj->lda, csysv_rook_obj->ipivref,
                                  csysv_rook_obj->bref, csysv_rook_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    csysv_rook_obj->info = LAPACKE_csysv_rook( csysv_rook_obj->matrix_layout,
                csysv_rook_obj->uplo, csysv_rook_obj->n, csysv_rook_obj->nrhs,
                                  ( lapack_complex_float *)csysv_rook_obj->a,
                               csysv_rook_obj->lda, csysv_rook_obj->ipiv,
                                 csysv_rook_obj->b, csysv_rook_obj->ldb );


    if( csysv_rook_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csysv_rook is wrong\n",
                    csysv_rook_obj->info );
    }
    if( csysv_rook_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csysv_rook is wrong\n",
        csysv_rook_obj->inforef );
    }
}

TEST_F(csysv_rook_test, csysv_rook1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_rook_obj->b_bufsize,
                           csysv_rook_obj->b, csysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csysv_rook_test, csysv_rook2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_rook_obj->b_bufsize,
                           csysv_rook_obj->b, csysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csysv_rook_test, csysv_rook3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_rook_obj->b_bufsize,
                           csysv_rook_obj->b, csysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csysv_rook_test, csysv_rook4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_rook_obj->b_bufsize,
                           csysv_rook_obj->b, csysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sysv_rook_dcomplex_parameters  class definition */
class sysv_rook_dcomplex_parameters{
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
      sysv_rook_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_rook_dcomplex_parameters ();
};  /* end of sysv_rook_dcomplex_parameters  class definition */


/* Constructor sysv_rook_dcomplex_parameters definition */
sysv_rook_dcomplex_parameters:: sysv_rook_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sysv_rook DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sysv_rook_dcomplex_parameters object: malloc error.";
       sysv_rook_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_rook_dcomplex_parameters:: ~sysv_rook_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_rook_dcomplex_parameters object: destructor invoked. \n");
#endif
   sysv_rook_free();
}


//  Test fixture class definition
class zsysv_rook_test  : public  ::testing::Test {
public:
   sysv_rook_dcomplex_parameters  *zsysv_rook_obj;
   void SetUp();  
   void TearDown () { delete zsysv_rook_obj; }
};


void zsysv_rook_test::SetUp(){

    /* LAPACKE ZSYSV_ROOK prototype */
    typedef int (*Fptr_NL_LAPACKE_zsysv_rook) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_double * a,
                                  lapack_int lda,  lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zsysv_rook ZSYSV_ROOK;

    zsysv_rook_obj = new  sysv_rook_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zsysv_rook_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsysv_rook_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsysv_rook_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsysv_rook_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYSV_ROOK = (Fptr_NL_LAPACKE_zsysv_rook)dlsym(zsysv_rook_obj->hModule, "LAPACKE_zsysv_rook");
    ASSERT_TRUE(ZSYSV_ROOK != NULL) << "failed to get the Netlib LAPACKE_zsysv_rook symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zsysv_rook_obj->inforef = ZSYSV_ROOK( zsysv_rook_obj->matrix_layout,
                                  zsysv_rook_obj->uplo, zsysv_rook_obj->n,
                                  zsysv_rook_obj->nrhs,
                                  ( lapack_complex_double *)zsysv_rook_obj->aref,
                                  zsysv_rook_obj->lda, zsysv_rook_obj->ipivref,
                                  zsysv_rook_obj->bref, zsysv_rook_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zsysv_rook_obj->info = LAPACKE_zsysv_rook( zsysv_rook_obj->matrix_layout,
                zsysv_rook_obj->uplo, zsysv_rook_obj->n, zsysv_rook_obj->nrhs,
                                  ( lapack_complex_double *)zsysv_rook_obj->a,
                               zsysv_rook_obj->lda, zsysv_rook_obj->ipiv,
                                 zsysv_rook_obj->b, zsysv_rook_obj->ldb );


    if( zsysv_rook_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsysv_rook is wrong\n",
                    zsysv_rook_obj->info );
    }
    if( zsysv_rook_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsysv_rook is wrong\n",
        zsysv_rook_obj->inforef );
    }
}

TEST_F(zsysv_rook_test, zsysv_rook1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_rook_obj->b_bufsize,
                           zsysv_rook_obj->b, zsysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsysv_rook_test, zsysv_rook2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_rook_obj->b_bufsize,
                           zsysv_rook_obj->b, zsysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsysv_rook_test, zsysv_rook3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_rook_obj->b_bufsize,
                           zsysv_rook_obj->b, zsysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsysv_rook_test, zsysv_rook4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_rook_obj->b_bufsize,
                           zsysv_rook_obj->b, zsysv_rook_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
