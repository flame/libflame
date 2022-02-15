#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define spsv_free() \
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

/* Begin spsv_double_parameters  class definition */
class spsv_double_parameters{
   public:
      int b_bufsize;
      double diff; //to capture reference and libflame o/ps' difference
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      spsv_double_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~spsv_double_parameters ();
};  /* end of spsv_double_parameters  class definition */


/* Constructor spsv_double_parameters definition */
spsv_double_parameters:: spsv_double_parameters ( int matrix_layout_i,
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
   printf(" \n spsv Double:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "spsv_double_parameters object: malloc error.";
       spsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

spsv_double_parameters:: ~spsv_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spsv_double_parameters object: destructor invoked. \n");
#endif
   spsv_free();
}


//  Test fixture class definition
class dspsv_test  : public  ::testing::Test {
public:
   spsv_double_parameters  *dspsv_obj;
   void SetUp();  
   void TearDown () { delete dspsv_obj; }
};


void dspsv_test::SetUp(){

    /* LAPACKE DSPSV prototype */
    typedef int (*Fptr_NL_LAPACKE_dspsv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,  double *a,
                                    lapack_int *ipiv,
                                            double *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dspsv DSPSV;

    dspsv_obj = new  spsv_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dspsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dspsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dspsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dspsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSPSV = (Fptr_NL_LAPACKE_dspsv)dlsym(dspsv_obj->hModule, "LAPACKE_dspsv");
    ASSERT_TRUE(DSPSV != NULL) << "failed to get the Netlib LAPACKE_dspsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    dspsv_obj->inforef = DSPSV( dspsv_obj->matrix_layout,
                                  dspsv_obj->uplo, dspsv_obj->n,
                                  dspsv_obj->nrhs,
                                  ( double *)dspsv_obj->aref,
                              dspsv_obj->ipivref,
                                dspsv_obj->bref, dspsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dspsv_obj->info = LAPACKE_dspsv( dspsv_obj->matrix_layout,
                dspsv_obj->uplo, dspsv_obj->n, dspsv_obj->nrhs,
                                  ( double *)dspsv_obj->a,
                               dspsv_obj->ipiv,
                                 dspsv_obj->b, dspsv_obj->ldb );


    if( dspsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dspsv is wrong\n",
                    dspsv_obj->info );
    }
    if( dspsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dspsv is wrong\n",
        dspsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dspsv_obj->diff =  computeDiff_d( dspsv_obj->b_bufsize,
                           dspsv_obj->b, dspsv_obj->bref );
}

TEST_F(dspsv_test, dspsv1) {
    EXPECT_NEAR(0.0, dspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dspsv_test, dspsv2) {
    EXPECT_NEAR(0.0, dspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dspsv_test, dspsv3) {
    EXPECT_NEAR(0.0, dspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dspsv_test, dspsv4) {
    EXPECT_NEAR(0.0, dspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin spsv_float_parameters  class definition */
class spsv_float_parameters{

   public:
      float diff; //to capture reference and libflame o/ps' difference
      int b_bufsize;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      spsv_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
      ~spsv_float_parameters ();
};  /* end of spsv_float_parameters  class definition */


/* Constructor spsv_float_parameters definition */
spsv_float_parameters:: spsv_float_parameters ( int matrix_layout_i,
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
   printf(" \n spsv float:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "spsv_double_parameters object: malloc error.";
       spsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

spsv_float_parameters:: ~spsv_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spsv_float_parameters object: destructor invoked. \n");
#endif
   spsv_free();
}


//  Test fixture class definition
class sspsv_test  : public  ::testing::Test {
public:
   spsv_float_parameters  *sspsv_obj;
   void SetUp();  
   void TearDown () { delete sspsv_obj; }
};


void sspsv_test::SetUp(){

    /* LAPACKE SSPSV prototype */
    typedef int (*Fptr_NL_LAPACKE_sspsv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,  float *a,
                                    lapack_int *ipiv,
                                            float *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_sspsv SSPSV;
    sspsv_obj = new  spsv_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    sspsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sspsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sspsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sspsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSPSV = (Fptr_NL_LAPACKE_sspsv)dlsym(sspsv_obj->hModule, "LAPACKE_sspsv");
    ASSERT_TRUE(SSPSV != NULL) << "failed to get the Netlib LAPACKE_sspsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    sspsv_obj->inforef = SSPSV( sspsv_obj->matrix_layout,
                                  sspsv_obj->uplo, sspsv_obj->n,
                                  sspsv_obj->nrhs,
                                  ( float *)sspsv_obj->aref,
                                   sspsv_obj->ipivref,
                                  sspsv_obj->bref, sspsv_obj->ldb);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    sspsv_obj->info = LAPACKE_sspsv( sspsv_obj->matrix_layout,
                sspsv_obj->uplo, sspsv_obj->n, sspsv_obj->nrhs,
                                  ( float *)sspsv_obj->a,
                                sspsv_obj->ipiv,
                                 sspsv_obj->b, sspsv_obj->ldb );
    if( sspsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sspsv is wrong\n",
                    sspsv_obj->info );
    }
    if( sspsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sspsv is wrong\n",
        sspsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    sspsv_obj->diff =  computeDiff_s( sspsv_obj->b_bufsize,
                           sspsv_obj->b,
                           sspsv_obj->bref );
}

TEST_F(sspsv_test, sspsv2) {
    EXPECT_NEAR(0.0, sspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sspsv_test, sspsv3) {
    EXPECT_NEAR(0.0, sspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sspsv_test, sspsv4) {
    EXPECT_NEAR(0.0, sspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin spsv_scomplex_parameters  class definition */
class spsv_scomplex_parameters{
   public:
      int b_bufsize;
      float diff; //to capture reference and libflame o/ps' difference
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      spsv_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~spsv_scomplex_parameters ();
};  /* end of spsv_scomplex_parameters  class definition */


/* Constructor spsv_scomplex_parameters definition */
spsv_scomplex_parameters:: spsv_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n spsv scomplex:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "spsv_scomplex_parameters object: malloc error.";
       spsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

spsv_scomplex_parameters:: ~spsv_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spsv_scomplex_parameters object: destructor invoked. \n");
#endif
   spsv_free();
}


//  Test fixture class definition
class cspsv_test  : public  ::testing::Test {
public:
   spsv_scomplex_parameters  *cspsv_obj;
   void SetUp();  
   void TearDown () { delete cspsv_obj; }
};


void cspsv_test::SetUp(){

    /* LAPACKE CSPSV prototype */
    typedef int (*Fptr_NL_LAPACKE_cspsv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_float *a,
                                    lapack_int *ipiv,
                              lapack_complex_float *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_cspsv CSPSV;

    cspsv_obj = new  spsv_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    cspsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cspsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cspsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cspsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSPSV = (Fptr_NL_LAPACKE_cspsv)dlsym(cspsv_obj->hModule, "LAPACKE_cspsv");
    ASSERT_TRUE(CSPSV != NULL) << "failed to get the Netlib LAPACKE_cspsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    cspsv_obj->inforef = CSPSV( cspsv_obj->matrix_layout,
                                  cspsv_obj->uplo, cspsv_obj->n,
                                  cspsv_obj->nrhs,
                      ( lapack_complex_float *)cspsv_obj->aref,
                                  cspsv_obj->ipivref,
                                  cspsv_obj->bref, cspsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    cspsv_obj->info = LAPACKE_cspsv( cspsv_obj->matrix_layout,
                cspsv_obj->uplo, cspsv_obj->n, cspsv_obj->nrhs,
                      ( lapack_complex_float *)cspsv_obj->a,
                cspsv_obj->ipiv, cspsv_obj->b, cspsv_obj->ldb );


    if( cspsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cspsv is wrong\n",
                    cspsv_obj->info );
    }
    if( cspsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cspsv is wrong\n",
        cspsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cspsv_obj->diff =  computeDiff_c( cspsv_obj->b_bufsize,
                           cspsv_obj->b, 
                           cspsv_obj->bref );
}

TEST_F(cspsv_test, cspsv1) {
    EXPECT_NEAR(0.0, cspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cspsv_test, cspsv2) {
    EXPECT_NEAR(0.0, cspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cspsv_test, cspsv3) {
    EXPECT_NEAR(0.0, cspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cspsv_test, cspsv4) {
    EXPECT_NEAR(0.0, cspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}



/* Begin spsv_dcomplex_parameters  class definition */
class spsv_dcomplex_parameters{
   public:
      int b_bufsize;
      double diff; //to capture reference and libflame o/ps' difference
     /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      spsv_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~spsv_dcomplex_parameters ();
};  /* end of spsv_dcomplex_parameters  class definition */


/* Constructor spsv_dcomplex_parameters definition */
spsv_dcomplex_parameters:: spsv_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n spsv DComplex:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "spsv_dcomplex_parameters object: malloc error.";
       spsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

spsv_dcomplex_parameters:: ~spsv_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" spsv_dcomplex_parameters object: destructor invoked. \n");
#endif
   spsv_free();
}


//  Test fixture class definition
class zspsv_test  : public  ::testing::Test {
public:
   spsv_dcomplex_parameters  *zspsv_obj;
   void SetUp();  
   void TearDown () { delete zspsv_obj; }
};


void zspsv_test::SetUp(){

    /* LAPACKE ZSPSV prototype */
    typedef int (*Fptr_NL_LAPACKE_zspsv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_double *a,
                                    lapack_int *ipiv,
                              lapack_complex_double *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zspsv ZSPSV;
    zspsv_obj = new  spsv_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zspsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zspsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zspsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zspsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSPSV = (Fptr_NL_LAPACKE_zspsv)dlsym(zspsv_obj->hModule, "LAPACKE_zspsv");
    ASSERT_TRUE(ZSPSV != NULL) << "failed to get the Netlib LAPACKE_zspsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zspsv_obj->inforef = ZSPSV( zspsv_obj->matrix_layout,
                                  zspsv_obj->uplo, zspsv_obj->n,
                                  zspsv_obj->nrhs,
                                  ( lapack_complex_double *)zspsv_obj->aref,
                                  zspsv_obj->ipivref,
                                  zspsv_obj->bref, zspsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zspsv_obj->info = LAPACKE_zspsv( zspsv_obj->matrix_layout,
                zspsv_obj->uplo, zspsv_obj->n, zspsv_obj->nrhs,
                                  ( lapack_complex_double *)zspsv_obj->a,
                               zspsv_obj->ipiv,
                                 zspsv_obj->b, zspsv_obj->ldb );


    if( zspsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zspsv is wrong\n",
                    zspsv_obj->info );
    }
    if( zspsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zspsv is wrong\n",
        zspsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zspsv_obj->diff =  computeDiff_z( zspsv_obj->b_bufsize,
                           zspsv_obj->b, zspsv_obj->bref );
}

TEST_F(zspsv_test, zspsv1) {
    EXPECT_NEAR(0.0, zspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zspsv_test, zspsv2) {
    EXPECT_NEAR(0.0, zspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zspsv_test, zspsv3) {
    EXPECT_NEAR(0.0, zspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zspsv_test, zspsv4) {
    EXPECT_NEAR(0.0, zspsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
