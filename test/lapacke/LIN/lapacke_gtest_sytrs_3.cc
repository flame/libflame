#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define sytrs_3_free() \
       if (ipiv != NULL) free (ipiv); \
       if (bref != NULL) free (bref); \
       if (b != NULL)    free (b   ); \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (e != NULL)    free (e   ); \
       if (eref != NULL) free (eref); \
       if (ipivref != NULL)free (ipivref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin sytrs_3_double_parameters  class definition */
class sytrs_3_double_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array 'a' contains the matrix A
      double *e, *eref; // superdiagonal (or subdiagonal) elements of
      //  the symmetric block diagonal matrix D
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sytrs_3_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_3_double_parameters ();
};  /* end of sytrs_3_double_parameters  class definition */


/* Constructor sytrs_3_double_parameters definition */
sytrs_3_double_parameters:: sytrs_3_double_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_3 Double:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
    lapacke_gtest_alloc_double_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_3_double_parameters object: malloc error.";
       sytrs_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( e, eref, n);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs_3_double_parameters:: ~sytrs_3_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_3_double_parameters object: destructor invoked. \n");
#endif
   sytrs_3_free();
}


//  Test fixture class definition
class dsytrs_3_test  : public  ::testing::Test {
public:
   sytrs_3_double_parameters  *dsytrs_3_obj;
   void SetUp();  
   void TearDown () { delete dsytrs_3_obj; }
};


void dsytrs_3_test::SetUp(){

    /* LAPACKE DSYTRS_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrs_3) ( int matrix_layout, char uplo,
                            lapack_int n, lapack_int nrhs, const double *a,
                   lapack_int lda, const double *e, const lapack_int *ipiv,
                                              double * b, lapack_int ldb );

    Fptr_NL_LAPACKE_dsytrs_3 DSYTRS_3;

     /* LAPACKE DSYTRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf_rk) ( int matrix_layout, char uplo,
     lapack_int n, double* a, lapack_int lda, double *e, lapack_int* ipiv );

    Fptr_NL_LAPACKE_dsytrf_rk DSYTRF_RK;

    dsytrs_3_obj = new  sytrs_3_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dsytrs_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytrs_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytrs_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytrs_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYTRS_3 = (Fptr_NL_LAPACKE_dsytrs_3)dlsym(dsytrs_3_obj->hModule, "LAPACKE_dsytrs_3");
    ASSERT_TRUE(DSYTRS_3 != NULL) << "failed to get the Netlib LAPACKE_dsytrs_3 symbol";

    DSYTRF_RK = (Fptr_NL_LAPACKE_dsytrf_rk)dlsym(dsytrs_3_obj->hModule,"LAPACKE_dsytrf_rk");
    ASSERT_TRUE(DSYTRF_RK != NULL) << "failed to get the Netlib LAPACKE_dsytrf_rk symbol";

    /* Pre condition: need to call sytrf_rk - before calling sytrs_3 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dsytrs_3_obj->inforef = DSYTRF_RK( dsytrs_3_obj->matrix_layout,
                                       dsytrs_3_obj->uplo,
                                       dsytrs_3_obj->n,
                                       dsytrs_3_obj->aref,
                                       dsytrs_3_obj->lda,
                                       dsytrs_3_obj->eref,
                                       dsytrs_3_obj->ipivref);

    dsytrs_3_obj->inforef = DSYTRS_3( dsytrs_3_obj->matrix_layout,
                                      dsytrs_3_obj->uplo,
                                      dsytrs_3_obj->n,
                                      dsytrs_3_obj->nrhs,
                                  (const double *)dsytrs_3_obj->aref,
                                      dsytrs_3_obj->lda,
                                  (const double *)dsytrs_3_obj->eref,
                                      dsytrs_3_obj->ipivref,
                                      dsytrs_3_obj->bref,
                                      dsytrs_3_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dsytrs_3_obj->info = LAPACKE_dsytrf_rk( dsytrs_3_obj->matrix_layout,
                                            dsytrs_3_obj->uplo,
                                            dsytrs_3_obj->n,
                                            dsytrs_3_obj->a,
                                            dsytrs_3_obj->lda,
                                            dsytrs_3_obj->e,
                                            dsytrs_3_obj->ipiv);

    dsytrs_3_obj->info = LAPACKE_dsytrs_3( dsytrs_3_obj->matrix_layout,
                                           dsytrs_3_obj->uplo,
                                           dsytrs_3_obj->n,
                                           dsytrs_3_obj->nrhs,
                                       (const double *)dsytrs_3_obj->a,
                                           dsytrs_3_obj->lda,
                                       (const double *)dsytrs_3_obj->e,
                                           dsytrs_3_obj->ipiv,
                                           dsytrs_3_obj->b,
                                           dsytrs_3_obj->ldb );

    if( dsytrs_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrs_3 is wrong\n",
                    dsytrs_3_obj->info );
    }

    if( dsytrs_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrs_3 is wrong\n",
        dsytrs_3_obj->inforef );
    }
}

TEST_F(dsytrs_3_test, dsytrs_31) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_3_obj->b_bufsize,
                           dsytrs_3_obj->b,
                           dsytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytrs_3_test, dsytrs_32) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_3_obj->b_bufsize,
                           dsytrs_3_obj->b,
                           dsytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytrs_3_test, dsytrs_33) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_3_obj->b_bufsize,
                           dsytrs_3_obj->b,
                           dsytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsytrs_3_test, dsytrs_34) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dsytrs_3_obj->b_bufsize,
                           dsytrs_3_obj->b,
                           dsytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin sytrs_3_float_parameters  class definition */
class sytrs_3_float_parameters{

   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array 'a' contains the matrix A
      float *e, *eref; // superdiagonal (or subdiagonal) elements of
      //  the symmetric block diagonal matrix D
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sytrs_3_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~sytrs_3_float_parameters ();
};  /* end of sytrs_3_float_parameters  class definition */


/* Constructor sytrs_3_float_parameters definition */
sytrs_3_float_parameters:: sytrs_3_float_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_3 float:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
    lapacke_gtest_alloc_float_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_3_double_parameters object: malloc error.";
       sytrs_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, lda*n);
    lapacke_gtest_init_float_buffer_pair_rand( e, eref, n);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

sytrs_3_float_parameters:: ~sytrs_3_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_3_float_parameters object: destructor invoked. \n");
#endif
   sytrs_3_free();
}


//  Test fixture class definition
class ssytrs_3_test  : public  ::testing::Test {
public:
   sytrs_3_float_parameters  *ssytrs_3_obj;
   void SetUp();  
   void TearDown () { delete ssytrs_3_obj; }
};


void ssytrs_3_test::SetUp(){

    /* LAPACKE SSYTRS_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrs_3) ( int matrix_layout,
                       char uplo, lapack_int n, lapack_int nrhs, 
                                const float * a, lapack_int lda, 
                        const float *e, const lapack_int * ipiv,
                                    float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_ssytrs_3 SSYTRS_3;

     /* LAPACKE SSYTRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf_rk) ( int matrix_layout, char uplo,
                                     lapack_int n, float* a, lapack_int lda, 
                                               float* e, lapack_int* ipiv );

    Fptr_NL_LAPACKE_ssytrf_rk SSYTRF_RK;

    ssytrs_3_obj = new  sytrs_3_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    ssytrs_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytrs_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytrs_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytrs_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYTRS_3 = (Fptr_NL_LAPACKE_ssytrs_3)dlsym(ssytrs_3_obj->hModule, "LAPACKE_ssytrs_3");
    ASSERT_TRUE(SSYTRS_3 != NULL) << "failed to get the Netlib LAPACKE_ssytrs_3 symbol";

    SSYTRF_RK = (Fptr_NL_LAPACKE_ssytrf_rk)dlsym(ssytrs_3_obj->hModule,"LAPACKE_ssytrf_rk");
    ASSERT_TRUE(SSYTRF_RK != NULL) << "failed to get the Netlib LAPACKE_ssytrf_rk symbol";

    /* Pre condition: need to call sytrf_rk - before calling sytrs_3 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    ssytrs_3_obj->inforef = SSYTRF_RK( ssytrs_3_obj->matrix_layout,
                                       ssytrs_3_obj->uplo,
                                       ssytrs_3_obj->n,
                                       ssytrs_3_obj->aref,
                                       ssytrs_3_obj->lda,
                                       ssytrs_3_obj->eref,
                                       ssytrs_3_obj->ipivref);

    ssytrs_3_obj->inforef = SSYTRS_3( ssytrs_3_obj->matrix_layout,
                                      ssytrs_3_obj->uplo,
                                      ssytrs_3_obj->n,
                                      ssytrs_3_obj->nrhs,
                                  (const float *)ssytrs_3_obj->aref,
                                      ssytrs_3_obj->lda,
                                  (const float *)ssytrs_3_obj->eref,
                                      ssytrs_3_obj->ipivref,
                                      ssytrs_3_obj->bref,
                                      ssytrs_3_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssytrs_3_obj->info = LAPACKE_ssytrf_rk( ssytrs_3_obj->matrix_layout,
                                            ssytrs_3_obj->uplo,
                                            ssytrs_3_obj->n,
                                            ssytrs_3_obj->a,
                                            ssytrs_3_obj->lda,
                                            ssytrs_3_obj->e,
                                            ssytrs_3_obj->ipiv);

    ssytrs_3_obj->info = LAPACKE_ssytrs_3(  ssytrs_3_obj->matrix_layout,
                                            ssytrs_3_obj->uplo,
                                            ssytrs_3_obj->n,
                                            ssytrs_3_obj->nrhs,
                                  (const float *)ssytrs_3_obj->a,
                                            ssytrs_3_obj->lda,
                                  (const float *)ssytrs_3_obj->e,
                                            ssytrs_3_obj->ipiv,
                                            ssytrs_3_obj->b,
                                            ssytrs_3_obj->ldb );
    if( ssytrs_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrs_3 is wrong\n",
                    ssytrs_3_obj->info );
    }
    if( ssytrs_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrs_3 is wrong\n",
        ssytrs_3_obj->inforef );
    }
}

TEST_F(ssytrs_3_test, ssytrs_31) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs_3_obj->b_bufsize,
                           ssytrs_3_obj->b, 
                           ssytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytrs_3_test, ssytrs_32) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs_3_obj->b_bufsize,
                           ssytrs_3_obj->b, 
                           ssytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytrs_3_test, ssytrs_33) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs_3_obj->b_bufsize,
                           ssytrs_3_obj->b, 
                           ssytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssytrs_3_test, ssytrs_34) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( ssytrs_3_obj->b_bufsize,
                           ssytrs_3_obj->b, 
                           ssytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sytrs_3_scomplex_parameters  class definition */
class sytrs_3_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A
      lapack_complex_float *e, *eref; // superdiagonal (or subdiagonal) elements of
      //  the symmetric block diagonal matrix D
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sytrs_3_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_3_scomplex_parameters ();
};  /* end of sytrs_3_scomplex_parameters  class definition */


/* Constructor sytrs_3_scomplex_parameters definition */
sytrs_3_scomplex_parameters:: sytrs_3_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_3 scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_3_scomplex_parameters object: malloc error.";
       sytrs_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( e, eref, n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs_3_scomplex_parameters:: ~sytrs_3_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_3_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrs_3_free();
}


//  Test fixture class definition
class csytrs_3_test  : public  ::testing::Test {
public:
   sytrs_3_scomplex_parameters  *csytrs_3_obj;
   void SetUp();  
   void TearDown () { delete csytrs_3_obj; }
};


void csytrs_3_test::SetUp(){

    /* LAPACKE CSYTRS_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrs_3) ( int matrix_layout, char uplo,
                                            lapack_int n, lapack_int nrhs,
                                            const lapack_complex_float *a,
                                            lapack_int lda, 
                                            const lapack_complex_float *e,
                                            const lapack_int *ipiv,
                                            lapack_complex_float *b,
                                            lapack_int ldb  );

    Fptr_NL_LAPACKE_csytrs_3 CSYTRS_3;

     /* LAPACKE CSYTRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf_rk) ( int matrix_layout, char uplo,
                      lapack_int n, lapack_complex_float *a, lapack_int lda, 
                                lapack_complex_float *e, lapack_int* ipiv );

    Fptr_NL_LAPACKE_csytrf_rk CSYTRF_RK;


    csytrs_3_obj = new  sytrs_3_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    csytrs_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytrs_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytrs_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytrs_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYTRS_3 = (Fptr_NL_LAPACKE_csytrs_3)dlsym(csytrs_3_obj->hModule, "LAPACKE_csytrs_3");
    ASSERT_TRUE(CSYTRS_3 != NULL) << "failed to get the Netlib LAPACKE_csytrs_3 symbol";

    CSYTRF_RK = (Fptr_NL_LAPACKE_csytrf_rk)dlsym(csytrs_3_obj->hModule,"LAPACKE_csytrf_rk");
    ASSERT_TRUE(CSYTRF_RK != NULL) << "failed to get the Netlib LAPACKE_csytrf_rk symbol";

    /* Pre condition: need to call sytrf_rk - before calling sytrs_3 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    csytrs_3_obj->inforef = CSYTRF_RK( csytrs_3_obj->matrix_layout,
                                     csytrs_3_obj->uplo,
                                     csytrs_3_obj->n,
                                     csytrs_3_obj->aref,
                                     csytrs_3_obj->lda,
                                     csytrs_3_obj->eref,
                                     csytrs_3_obj->ipivref );

    csytrs_3_obj->inforef = CSYTRS_3( csytrs_3_obj->matrix_layout,
                                    csytrs_3_obj->uplo, csytrs_3_obj->n,
                                    csytrs_3_obj->nrhs,
                            (const lapack_complex_float *)csytrs_3_obj->aref,
                                    csytrs_3_obj->lda,
                            (const lapack_complex_float *)csytrs_3_obj->eref,
                                    csytrs_3_obj->ipivref,
                                    csytrs_3_obj->bref, csytrs_3_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    csytrs_3_obj->info = LAPACKE_csytrf_rk( csytrs_3_obj->matrix_layout,
                                            csytrs_3_obj->uplo, csytrs_3_obj->n,
                                            csytrs_3_obj->a,
                                            csytrs_3_obj->lda,
                                            csytrs_3_obj->e,
                                            csytrs_3_obj->ipiv);

    csytrs_3_obj->info = LAPACKE_csytrs_3( csytrs_3_obj->matrix_layout,
                                           csytrs_3_obj->uplo,
                                           csytrs_3_obj->n,
                                           csytrs_3_obj->nrhs,
                                  (const lapack_complex_float *)csytrs_3_obj->a,
                                           csytrs_3_obj->lda,
                                  (const lapack_complex_float *)csytrs_3_obj->e,
                                           csytrs_3_obj->ipiv,
                                           csytrs_3_obj->b,
                                           csytrs_3_obj->ldb );


    if( csytrs_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrs_3 is wrong\n",
                    csytrs_3_obj->info );
    }
    if( csytrs_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrs_3 is wrong\n",
        csytrs_3_obj->inforef );
    }
}

TEST_F(csytrs_3_test, csytrs_31) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_3_obj->b_bufsize,
                           csytrs_3_obj->b, csytrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs_3_test, csytrs_32) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_3_obj->b_bufsize,
                           csytrs_3_obj->b, csytrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(csytrs_3_test, csytrs_33) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_3_obj->b_bufsize,
                           csytrs_3_obj->b,
                           csytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csytrs_3_test, csytrs_34) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( csytrs_3_obj->b_bufsize,
                           csytrs_3_obj->b,
                           csytrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin sytrs_3_dcomplex_parameters  class definition */
class sytrs_3_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A
      lapack_complex_double *e, *eref; // superdiagonal (or subdiagonal) elements of
      //  the symmetric block diagonal matrix D
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      sytrs_3_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~sytrs_3_dcomplex_parameters ();
};  /* end of sytrs_3_dcomplex_parameters  class definition */


/* Constructor sytrs_3_dcomplex_parameters definition */
sytrs_3_dcomplex_parameters:: sytrs_3_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n sytrs_3 DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "sytrs_3_dcomplex_parameters object: malloc error.";
       sytrs_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( e, eref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

sytrs_3_dcomplex_parameters:: ~sytrs_3_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" sytrs_3_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrs_3_free();
}


//  Test fixture class definition
class zsytrs_3_test  : public  ::testing::Test {
public:
   sytrs_3_dcomplex_parameters  *zsytrs_3_obj;
   void SetUp();  
   void TearDown () { delete zsytrs_3_obj; }
};


void zsytrs_3_test::SetUp(){

    /* LAPACKE ZSYTRS_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrs_3)( int matrix_layout, char uplo,
                                            lapack_int n, lapack_int nrhs,
                                            const lapack_complex_double * a,
                                            lapack_int lda, 
                                            const lapack_complex_double * e,
                                            const lapack_int * ipiv,
                                            lapack_complex_double *b,
                                            lapack_int ldb  );

    Fptr_NL_LAPACKE_zsytrs_3 ZSYTRS_3;

     /* LAPACKE ZSYTRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf_rk) ( int matrix_layout,char uplo,
                                               lapack_int n,
                                               lapack_complex_double* a,
                                               lapack_int lda,
                                               lapack_complex_double* e,
                                               lapack_int* ipiv );

    Fptr_NL_LAPACKE_zsytrf_rk ZSYTRF_RK;


    zsytrs_3_obj = new  sytrs_3_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zsytrs_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytrs_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytrs_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytrs_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYTRS_3 = (Fptr_NL_LAPACKE_zsytrs_3)dlsym(zsytrs_3_obj->hModule, "LAPACKE_zsytrs_3");
    ASSERT_TRUE(ZSYTRS_3 != NULL) << "failed to get the Netlib LAPACKE_zsytrs_3 symbol";

    ZSYTRF_RK = (Fptr_NL_LAPACKE_zsytrf_rk)dlsym(zsytrs_3_obj->hModule,"LAPACKE_zsytrf_rk");
    ASSERT_TRUE(ZSYTRF_RK != NULL) << "failed to get the Netlib LAPACKE_zsytrf_rk symbol";

    /* Pre condition: need to call sytrf_rk - before calling sytrs_3 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zsytrs_3_obj->inforef = ZSYTRF_RK( zsytrs_3_obj->matrix_layout,
                                       zsytrs_3_obj->uplo,
                                       zsytrs_3_obj->n,
                                       zsytrs_3_obj->aref,
                                       zsytrs_3_obj->lda,
                                       zsytrs_3_obj->eref,
                                       zsytrs_3_obj->ipivref );
                                     
    zsytrs_3_obj->inforef = ZSYTRS_3(zsytrs_3_obj->matrix_layout,
                                     zsytrs_3_obj->uplo, zsytrs_3_obj->n,
                                     zsytrs_3_obj->nrhs,
                            (const lapack_complex_double *)zsytrs_3_obj->aref,
                                     zsytrs_3_obj->lda,
                            (const lapack_complex_double *)zsytrs_3_obj->eref,
                                     zsytrs_3_obj->ipivref,
                                     zsytrs_3_obj->bref, zsytrs_3_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zsytrs_3_obj->info = LAPACKE_zsytrf_rk( zsytrs_3_obj->matrix_layout,
                                            zsytrs_3_obj->uplo, zsytrs_3_obj->n,
                                            zsytrs_3_obj->a,
                                            zsytrs_3_obj->lda,
                                            zsytrs_3_obj->e,
                                            zsytrs_3_obj->ipiv);

    zsytrs_3_obj->info = LAPACKE_zsytrs_3( zsytrs_3_obj->matrix_layout,
                                           zsytrs_3_obj->uplo,
                                           zsytrs_3_obj->n,
                                           zsytrs_3_obj->nrhs,
                                  (const lapack_complex_double *)zsytrs_3_obj->a,
                                           zsytrs_3_obj->lda,
                                  (const lapack_complex_double *)zsytrs_3_obj->e,
                                           zsytrs_3_obj->ipiv,
                                           zsytrs_3_obj->b,
                                           zsytrs_3_obj->ldb );

    if( zsytrs_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrs_3 is wrong\n",
                    zsytrs_3_obj->info );
    }
    if( zsytrs_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrs_3 is wrong\n",
        zsytrs_3_obj->inforef );
    }
}

TEST_F(zsytrs_3_test, zsytrs_31) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_3_obj->b_bufsize,
                           zsytrs_3_obj->b, 
                           zsytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytrs_3_test, zsytrs_32) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_3_obj->b_bufsize,
                           zsytrs_3_obj->b, 
                           zsytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytrs_3_test, zsytrs_33) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_3_obj->b_bufsize,
                           zsytrs_3_obj->b, 
                           zsytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsytrs_3_test, zsytrs_34) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zsytrs_3_obj->b_bufsize,
                           zsytrs_3_obj->b, 
                           zsytrs_3_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}