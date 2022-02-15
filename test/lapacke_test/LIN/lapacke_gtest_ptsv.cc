#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define ptsv_free() \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (d != NULL)    free (d   ); \
       if (dref != NULL) free (dref); \
       if (e != NULL)    free (e   ); \
       if (eref != NULL) free (eref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin ptsv_double_parameters  class definition */
class ptsv_double_parameters{
   public:
   
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      double * e, *eref; // (n - 1) multipliers that define the matrix L
      double * d, *dref; //   n diagonal elements of the upper triangular matrix U
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      ptsv_double_parameters ( int matrix_layout_i, lapack_int n_i,
                               lapack_int nrhs_i, lapack_int ldb_i);
             
      ~ptsv_double_parameters ();
};  /* end of ptsv_double_parameters  class definition */


/* Constructor ptsv_double_parameters definition */
ptsv_double_parameters:: ptsv_double_parameters ( int matrix_layout_i, 
                 lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n ptsv Double:  n: %d, trans: %c ldb: %d nrhs: %d \n",
             n, trans, ldb, nrhs);
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
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_double_buffer_pair( &e, &eref, (n-1));

    if( (b==NULL) || (bref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (d==NULL) || (dref==NULL)   ){
       ptsv_free();
       EXPECT_FALSE( true) << "ptsv_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_double_buffer_pair_rand( e, eref, (n-1));

   } /* end of Constructor  */

ptsv_double_parameters:: ~ptsv_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptsv_double_parameters object: destructor invoked. \n");
#endif
   ptsv_free();
}


//  Test fixture class definition
class dptsv_test  : public  ::testing::Test {
public:
   ptsv_double_parameters  *dptsv_obj;
   void SetUp();  
   void TearDown () { delete dptsv_obj; }
};


void dptsv_test::SetUp(){

    /* LAPACKE DPTSV prototype */
    typedef int (*Fptr_NL_LAPACKE_dptsv) ( int matrix_layout, lapack_int n,
                          lapack_int nrhs,  double* d,  double* e,
                              double* b, lapack_int ldb);

    Fptr_NL_LAPACKE_dptsv DPTSV;

    dptsv_obj = new  ptsv_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dptsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dptsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dptsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dptsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DPTSV = (Fptr_NL_LAPACKE_dptsv)dlsym(dptsv_obj->hModule, "LAPACKE_dptsv");
    ASSERT_TRUE(DPTSV != NULL) << "failed to get the Netlib LAPACKE_dptsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    dptsv_obj->inforef = DPTSV( dptsv_obj->matrix_layout,
                                  dptsv_obj->n,
                                  dptsv_obj->nrhs,
                                  ( double *)dptsv_obj->dref,
                                  ( double *)dptsv_obj->eref,
                                  dptsv_obj->bref,
                                  dptsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dptsv_obj->info = LAPACKE_dptsv( dptsv_obj->matrix_layout,
                                 dptsv_obj->n, dptsv_obj->nrhs,
                                  ( double *)dptsv_obj->d,
                                  ( double *)dptsv_obj->e,
                                  dptsv_obj->b, dptsv_obj->ldb );

    if( dptsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
                   LAPACKE_dptsv is wrong\n", dptsv_obj->info );
    }
    if( dptsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dptsv \
                            is wrong\n",  dptsv_obj->inforef );
    }
}

TEST_F(dptsv_test, dptsv1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dptsv_obj->b_bufsize,
                           dptsv_obj->b, dptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dptsv_test, dptsv2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dptsv_obj->b_bufsize,
                           dptsv_obj->b, dptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dptsv_test, dptsv3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dptsv_obj->b_bufsize,
                           dptsv_obj->b, dptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dptsv_test, dptsv4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dptsv_obj->b_bufsize,
                           dptsv_obj->b, dptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin ptsv_float_parameters  class definition */
class ptsv_float_parameters{
   public:
   
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      float * e, *eref; // (n - 1) multipliers that define the matrix L
      float * d, *dref; //   n diagonal elements of the upper triangular matrix U
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      ptsv_float_parameters ( int matrix_layout_i, lapack_int n_i,
                               lapack_int nrhs_i, lapack_int ldb_i);
             
      ~ptsv_float_parameters ();
};  /* end of ptsv_float_parameters  class definition */


/* Constructor ptsv_float_parameters definition */
ptsv_float_parameters:: ptsv_float_parameters ( int matrix_layout_i, 
               lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n ptsv float:  n: %d, trans: %c ldb: %d nrhs: %d \n",
             n, trans, ldb, nrhs);
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
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_float_buffer_pair( &e, &eref, (n-1));

    if( (b==NULL) || (bref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (d==NULL) || (dref==NULL)   ){
       ptsv_free();
       EXPECT_FALSE( true) << "ptsv_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_float_buffer_pair_rand( e, eref, (n-1));

   } /* end of Constructor  */

ptsv_float_parameters:: ~ptsv_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptsv_float_parameters object: destructor invoked. \n");
#endif
   ptsv_free();
}


//  Test fixture class definition
class sptsv_test  : public  ::testing::Test {
public:
   ptsv_float_parameters  *sptsv_obj;
   void SetUp();  
   void TearDown () { delete sptsv_obj; }
};


void sptsv_test::SetUp(){

    /* LAPACKE SPTSV prototype */
    typedef int (*Fptr_NL_LAPACKE_sptsv) ( int matrix_layout, lapack_int n,
                          lapack_int nrhs,  float* d,  float* e,
                              float* b, lapack_int ldb);

    Fptr_NL_LAPACKE_sptsv SPTSV;

    sptsv_obj = new  ptsv_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    sptsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sptsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sptsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sptsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SPTSV = (Fptr_NL_LAPACKE_sptsv)dlsym(sptsv_obj->hModule, "LAPACKE_sptsv");
    ASSERT_TRUE(SPTSV != NULL) << "failed to get the Netlib LAPACKE_sptsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    sptsv_obj->inforef = SPTSV( sptsv_obj->matrix_layout,
                                  sptsv_obj->n,
                                  sptsv_obj->nrhs,
                                  ( float *)sptsv_obj->dref,
                                  ( float *)sptsv_obj->eref,
                                  sptsv_obj->bref,
                                  sptsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    sptsv_obj->info = LAPACKE_sptsv( sptsv_obj->matrix_layout,
                                 sptsv_obj->n, sptsv_obj->nrhs,
                                  ( float *)sptsv_obj->d,
                                  ( float *)sptsv_obj->e,
                                  sptsv_obj->b, sptsv_obj->ldb );

    if( sptsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
                   LAPACKE_sptsv is wrong\n", sptsv_obj->info );
    }
    if( sptsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sptsv \
                            is wrong\n",  sptsv_obj->inforef );
    }
}

TEST_F(sptsv_test, sptsv1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sptsv_obj->b_bufsize,
                           sptsv_obj->b, sptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sptsv_test, sptsv2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sptsv_obj->b_bufsize,
                           sptsv_obj->b, sptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sptsv_test, sptsv3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sptsv_obj->b_bufsize,
                           sptsv_obj->b, sptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sptsv_test, sptsv4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sptsv_obj->b_bufsize,
                           sptsv_obj->b, sptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin ptsv_scomplex_parameters  class definition */
class ptsv_scomplex_parameters{
   public:
   
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float * e, *eref; // (n - 1) multipliers that define the matrix L
      float * d, *dref; //   n diagonal elements of the upper triangular matrix U
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      ptsv_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
             lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~ptsv_scomplex_parameters ();
};  /* end of ptsv_scomplex_parameters  class definition */


/* Constructor ptsv_scomplex_parameters definition */
ptsv_scomplex_parameters:: ptsv_scomplex_parameters ( int matrix_layout_i, 
         char uplo_i, lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    n = n_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n ptsv lapack_complex_float:  n: %d, trans: %c ldb: %d nrhs: %d \n",
             n, trans, ldb, nrhs);
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
    lapacke_gtest_alloc_float_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e, &eref, (n-1));

    if( (b==NULL) || (bref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (d==NULL) || (dref==NULL)   ){
       ptsv_free();
       EXPECT_FALSE( true) << "ptsv_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_scomplex_buffer_pair_rand( e, eref, (n-1));

   } /* end of Constructor  */

ptsv_scomplex_parameters:: ~ptsv_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptsv_scomplex_parameters object: destructor invoked. \n");
#endif
   ptsv_free();
}


//  Test fixture class definition
class cptsv_test  : public  ::testing::Test {
public:
   ptsv_scomplex_parameters  *cptsv_obj;
   void SetUp();  
   void TearDown () { delete cptsv_obj; }
};


void cptsv_test::SetUp(){

    /* LAPACKE CPTSV prototype */
    typedef int (*Fptr_NL_LAPACKE_cptsv) ( int matrix_layout, 
                           lapack_int n, lapack_int nrhs,  float* d, 
                                            lapack_complex_float* e,
                                lapack_complex_float* b, lapack_int ldb);

    Fptr_NL_LAPACKE_cptsv CPTSV;
    cptsv_obj = new  ptsv_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    cptsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cptsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cptsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cptsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CPTSV = (Fptr_NL_LAPACKE_cptsv)dlsym(cptsv_obj->hModule, "LAPACKE_cptsv");
    ASSERT_TRUE(CPTSV != NULL) << "failed to get the Netlib LAPACKE_cptsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    cptsv_obj->inforef = CPTSV( cptsv_obj->matrix_layout,
                                  
                                  cptsv_obj->n,
                                  cptsv_obj->nrhs,
                                  cptsv_obj->dref,
                                  ( lapack_complex_float *)cptsv_obj->eref,
                                  cptsv_obj->bref,
                                  cptsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    cptsv_obj->info = LAPACKE_cptsv( cptsv_obj->matrix_layout,
                cptsv_obj->n, cptsv_obj->nrhs,
                                    cptsv_obj->d,
                     ( lapack_complex_float *)cptsv_obj->e,
                                  cptsv_obj->b, cptsv_obj->ldb );

    if( cptsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
                   LAPACKE_cptsv is wrong\n", cptsv_obj->info );
    }
    if( cptsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cptsv \
                            is wrong\n",  cptsv_obj->inforef );
    }
}

TEST_F(cptsv_test, cptsv1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cptsv_obj->b_bufsize,
                           cptsv_obj->b, cptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cptsv_test, cptsv2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cptsv_obj->b_bufsize,
                           cptsv_obj->b, cptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cptsv_test, cptsv3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cptsv_obj->b_bufsize,
                           cptsv_obj->b, cptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cptsv_test, cptsv4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cptsv_obj->b_bufsize,
                           cptsv_obj->b, cptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin ptsv_dcomplex_parameters  class definition */
class ptsv_dcomplex_parameters{
   public:
   
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      char uplo; //  Must be 'U' or 'L'.
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double * e, *eref; // (n - 1) multipliers that define the matrix L
      double * d, *dref; //   n diagonal elements of the upper triangular matrix U
      void *hModule, *dModule;
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      ptsv_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
             lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~ptsv_dcomplex_parameters ();
};  /* end of ptsv_dcomplex_parameters  class definition */


/* Constructor ptsv_dcomplex_parameters definition */
ptsv_dcomplex_parameters:: ptsv_dcomplex_parameters ( int matrix_layout_i, 
         char uplo_i, lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n ptsv lapack_complex_double:  n: %d, trans: %c ldb: %d nrhs: %d \n",
             n, trans, ldb, nrhs);
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
    lapacke_gtest_alloc_double_buffer_pair( &d, &dref, n);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e, &eref, (n-1));

    if( (b==NULL) || (bref==NULL) ||  \
        (e==NULL) || (eref==NULL) ||  \
        (d==NULL) || (dref==NULL)   ){
       ptsv_free();
       EXPECT_FALSE( true) << "ptsv_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( e, eref, (n-1));

   } /* end of Constructor  */

ptsv_dcomplex_parameters:: ~ptsv_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ptsv_dcomplex_parameters object: destructor invoked. \n");
#endif
   ptsv_free();
}


//  Test fixture class definition
class zptsv_test  : public  ::testing::Test {
public:
   ptsv_dcomplex_parameters  *zptsv_obj;
   void SetUp();  
   void TearDown () { delete zptsv_obj; }
};


void zptsv_test::SetUp(){

    /* LAPACKE ZPTSV prototype */
    typedef int (*Fptr_NL_LAPACKE_zptsv) ( int matrix_layout, 
                          lapack_int n, lapack_int nrhs,  double* d, 
                                            lapack_complex_double* e,
                                lapack_complex_double* b, lapack_int ldb);

    Fptr_NL_LAPACKE_zptsv ZPTSV;

    zptsv_obj = new  ptsv_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    zptsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zptsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zptsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zptsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZPTSV = (Fptr_NL_LAPACKE_zptsv)dlsym(zptsv_obj->hModule, "LAPACKE_zptsv");
    ASSERT_TRUE(ZPTSV != NULL) << "failed to get the Netlib LAPACKE_zptsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zptsv_obj->inforef = ZPTSV( zptsv_obj->matrix_layout,
                                  zptsv_obj->n,
                                  zptsv_obj->nrhs,
                                  zptsv_obj->dref,
                                  ( lapack_complex_double *)zptsv_obj->eref,
                                  zptsv_obj->bref,
                                  zptsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zptsv_obj->info = LAPACKE_zptsv( zptsv_obj->matrix_layout,
                                 zptsv_obj->n, zptsv_obj->nrhs,
                                  zptsv_obj->d,
                    ( lapack_complex_double *)zptsv_obj->e,
                                  zptsv_obj->b, zptsv_obj->ldb );

    if( zptsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
                   LAPACKE_zptsv is wrong\n", zptsv_obj->info );
    }
    if( zptsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zptsv \
                            is wrong\n",  zptsv_obj->inforef );
    }
}

TEST_F(zptsv_test, zptsv1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zptsv_obj->b_bufsize,
                           zptsv_obj->b, zptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zptsv_test, zptsv2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zptsv_obj->b_bufsize,
                           zptsv_obj->b, zptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zptsv_test, zptsv3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zptsv_obj->b_bufsize,
                           zptsv_obj->b, zptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zptsv_test, zptsv4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zptsv_obj->b_bufsize,
                           zptsv_obj->b, zptsv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}