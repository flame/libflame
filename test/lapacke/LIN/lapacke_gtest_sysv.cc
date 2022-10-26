#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define sysv_free() \
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


/* Begin sysv_double_parameters  class definition */
class sysv_double_parameters{
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
      sysv_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_double_parameters ();
};  /* end of sysv_double_parameters  class definition */


/* Constructor sysv_double_parameters definition */
sysv_double_parameters:: sysv_double_parameters ( int matrix_layout_i,
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
   printf(" \n sysv Double:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sysv_double_parameters object: malloc error.";
       sysv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_double_parameters:: ~sysv_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_double_parameters object: destructor invoked. \n");
#endif
   sysv_free();
}


//  Test fixture class definition
class dsysv_test  : public  ::testing::Test {
public:
   sysv_double_parameters  *dsysv_obj;
   void SetUp();  
   void TearDown () { delete dsysv_obj; }
};


void dsysv_test::SetUp(){

    /* LAPACKE DSYSV prototype */
    typedef int (*Fptr_NL_LAPACKE_dsysv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,  double * a,
                                  lapack_int lda,  lapack_int * ipiv,
                                            double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dsysv DSYSV;

    dsysv_obj = new  sysv_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dsysv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsysv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsysv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsysv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYSV = (Fptr_NL_LAPACKE_dsysv)dlsym(dsysv_obj->hModule, "LAPACKE_dsysv");
    ASSERT_TRUE(DSYSV != NULL) << "failed to get the Netlib LAPACKE_dsysv symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    dsysv_obj->inforef = DSYSV( dsysv_obj->matrix_layout,
                                  dsysv_obj->uplo, dsysv_obj->n,
                                  dsysv_obj->nrhs,
                                  ( double *)dsysv_obj->aref,
                              dsysv_obj->lda, dsysv_obj->ipivref,
                                dsysv_obj->bref, dsysv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dsysv_obj->info = LAPACKE_dsysv( dsysv_obj->matrix_layout,
                dsysv_obj->uplo, dsysv_obj->n, dsysv_obj->nrhs,
                                  ( double *)dsysv_obj->a,
                               dsysv_obj->lda, dsysv_obj->ipiv,
                                 dsysv_obj->b, dsysv_obj->ldb );


    if( dsysv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsysv is wrong\n",
                    dsysv_obj->info );
    }
    if( dsysv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsysv is wrong\n",
        dsysv_obj->inforef );
    }
}

TEST_F(dsysv_test, dsysv1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_obj->b_bufsize,
                           dsysv_obj->b, dsysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsysv_test, dsysv2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_obj->b_bufsize,
                           dsysv_obj->b, dsysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsysv_test, dsysv3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_obj->b_bufsize,
                           dsysv_obj->b, dsysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dsysv_test, dsysv4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsysv_obj->b_bufsize,
                           dsysv_obj->b, dsysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sysv_float_parameters  class definition */
class sysv_float_parameters{

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
      sysv_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~sysv_float_parameters ();
};  /* end of sysv_float_parameters  class definition */


/* Constructor sysv_float_parameters definition */
sysv_float_parameters:: sysv_float_parameters ( int matrix_layout_i,
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
   printf(" \n sysv float:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sysv_double_parameters object: malloc error.";
       sysv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

sysv_float_parameters:: ~sysv_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_float_parameters object: destructor invoked. \n");
#endif
   sysv_free();
}


//  Test fixture class definition
class ssysv_test  : public  ::testing::Test {
public:
   sysv_float_parameters  *ssysv_obj;
   void SetUp();  
   void TearDown () { delete ssysv_obj; }
};


void ssysv_test::SetUp(){

    /* LAPACKE SSYSV prototype */
    typedef int (*Fptr_NL_LAPACKE_ssysv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,  float * a,
                                  lapack_int lda,  lapack_int * ipiv,
                                            float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_ssysv SSYSV;

    ssysv_obj = new  sysv_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    ssysv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssysv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssysv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssysv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYSV = (Fptr_NL_LAPACKE_ssysv)dlsym(ssysv_obj->hModule, "LAPACKE_ssysv");
    ASSERT_TRUE(SSYSV != NULL) << "failed to get the Netlib LAPACKE_ssysv symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    ssysv_obj->inforef = SSYSV( ssysv_obj->matrix_layout,
                                  ssysv_obj->uplo, ssysv_obj->n,
                                  ssysv_obj->nrhs,
                                  ( float *)ssysv_obj->aref,
                                  ssysv_obj->lda, ssysv_obj->ipivref,
                                  ssysv_obj->bref, ssysv_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssysv_obj->info = LAPACKE_ssysv( ssysv_obj->matrix_layout,
                ssysv_obj->uplo, ssysv_obj->n, ssysv_obj->nrhs,
                                  ( float *)ssysv_obj->a,
                               ssysv_obj->lda, ssysv_obj->ipiv,
                                 ssysv_obj->b, ssysv_obj->ldb );
    if( ssysv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssysv is wrong\n",
                    ssysv_obj->info );
    }
    if( ssysv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssysv is wrong\n",
        ssysv_obj->inforef );
    }
}

TEST_F(ssysv_test, ssysv1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_obj->b_bufsize,
                           ssysv_obj->b, ssysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_test, ssysv2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_obj->b_bufsize,
                           ssysv_obj->b, ssysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_test, ssysv3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_obj->b_bufsize,
                           ssysv_obj->b, ssysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssysv_test, ssysv4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssysv_obj->b_bufsize,
                           ssysv_obj->b, ssysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sysv_scomplex_parameters  class definition */
class sysv_scomplex_parameters{
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
      sysv_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_scomplex_parameters ();
};  /* end of sysv_scomplex_parameters  class definition */


/* Constructor sysv_scomplex_parameters definition */
sysv_scomplex_parameters:: sysv_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sysv scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sysv_scomplex_parameters object: malloc error.";
       sysv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_scomplex_parameters:: ~sysv_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_scomplex_parameters object: destructor invoked. \n");
#endif
   sysv_free();
}


//  Test fixture class definition
class csysv_test  : public  ::testing::Test {
public:
   sysv_scomplex_parameters  *csysv_obj;
   void SetUp();  
   void TearDown () { delete csysv_obj; }
};


void csysv_test::SetUp(){

    /* LAPACKE CSYSV prototype */
    typedef int (*Fptr_NL_LAPACKE_csysv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_float * a,
                                  lapack_int lda,  lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_csysv CSYSV;

    csysv_obj = new  sysv_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    csysv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csysv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csysv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csysv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYSV = (Fptr_NL_LAPACKE_csysv)dlsym(csysv_obj->hModule, "LAPACKE_csysv");
    ASSERT_TRUE(CSYSV != NULL) << "failed to get the Netlib LAPACKE_csysv symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    csysv_obj->inforef = CSYSV( csysv_obj->matrix_layout,
                                  csysv_obj->uplo, csysv_obj->n,
                                  csysv_obj->nrhs,
                                  ( lapack_complex_float *)csysv_obj->aref,
                                  csysv_obj->lda, csysv_obj->ipivref,
                                  csysv_obj->bref, csysv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    csysv_obj->info = LAPACKE_csysv( csysv_obj->matrix_layout,
                csysv_obj->uplo, csysv_obj->n, csysv_obj->nrhs,
                                  ( lapack_complex_float *)csysv_obj->a,
                               csysv_obj->lda, csysv_obj->ipiv,
                                 csysv_obj->b, csysv_obj->ldb );


    if( csysv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csysv is wrong\n",
                    csysv_obj->info );
    }
    if( csysv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csysv is wrong\n",
        csysv_obj->inforef );
    }
}

TEST_F(csysv_test, csysv1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_obj->b_bufsize,
                           csysv_obj->b, csysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csysv_test, csysv2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_obj->b_bufsize,
                           csysv_obj->b, csysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csysv_test, csysv3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_obj->b_bufsize,
                           csysv_obj->b, csysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csysv_test, csysv4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csysv_obj->b_bufsize,
                           csysv_obj->b, csysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sysv_dcomplex_parameters  class definition */
class sysv_dcomplex_parameters{
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
      sysv_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sysv_dcomplex_parameters ();
};  /* end of sysv_dcomplex_parameters  class definition */


/* Constructor sysv_dcomplex_parameters definition */
sysv_dcomplex_parameters:: sysv_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sysv DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "sysv_dcomplex_parameters object: malloc error.";
       sysv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sysv_dcomplex_parameters:: ~sysv_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sysv_dcomplex_parameters object: destructor invoked. \n");
#endif
   sysv_free();
}


//  Test fixture class definition
class zsysv_test  : public  ::testing::Test {
public:
   sysv_dcomplex_parameters  *zsysv_obj;
   void SetUp();  
   void TearDown () { delete zsysv_obj; }
};


void zsysv_test::SetUp(){

    /* LAPACKE ZSYSV prototype */
    typedef int (*Fptr_NL_LAPACKE_zsysv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_double * a,
                                  lapack_int lda,  lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zsysv ZSYSV;

    zsysv_obj = new  sysv_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zsysv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsysv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsysv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsysv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYSV = (Fptr_NL_LAPACKE_zsysv)dlsym(zsysv_obj->hModule, "LAPACKE_zsysv");
    ASSERT_TRUE(ZSYSV != NULL) << "failed to get the Netlib LAPACKE_zsysv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zsysv_obj->inforef = ZSYSV( zsysv_obj->matrix_layout,
                                  zsysv_obj->uplo, zsysv_obj->n,
                                  zsysv_obj->nrhs,
                                  ( lapack_complex_double *)zsysv_obj->aref,
                                  zsysv_obj->lda, zsysv_obj->ipivref,
                                  zsysv_obj->bref, zsysv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zsysv_obj->info = LAPACKE_zsysv( zsysv_obj->matrix_layout,
                zsysv_obj->uplo, zsysv_obj->n, zsysv_obj->nrhs,
                                  ( lapack_complex_double *)zsysv_obj->a,
                               zsysv_obj->lda, zsysv_obj->ipiv,
                                 zsysv_obj->b, zsysv_obj->ldb );


    if( zsysv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsysv is wrong\n",
                    zsysv_obj->info );
    }
    if( zsysv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsysv is wrong\n",
        zsysv_obj->inforef );
    }
}

TEST_F(zsysv_test, zsysv1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_obj->b_bufsize,
                           zsysv_obj->b, zsysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsysv_test, zsysv2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_obj->b_bufsize,
                           zsysv_obj->b, zsysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsysv_test, zsysv3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_obj->b_bufsize,
                           zsysv_obj->b, zsysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zsysv_test, zsysv4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsysv_obj->b_bufsize,
                           zsysv_obj->b, zsysv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
