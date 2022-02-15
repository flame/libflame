#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define pttrs_free() \
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


/* Begin pttrs_double_parameters  class definition */
class pttrs_double_parameters{
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
      pttrs_double_parameters ( int matrix_layout_i, lapack_int n_i,
                               lapack_int nrhs_i, lapack_int ldb_i);
             
      ~pttrs_double_parameters ();
};  /* end of pttrs_double_parameters  class definition */


/* Constructor pttrs_double_parameters definition */
pttrs_double_parameters:: pttrs_double_parameters ( int matrix_layout_i, 
                 lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pttrs Double:  n: %d, trans: %c ldb: %d nrhs: %d \n",
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
       pttrs_free();
       EXPECT_FALSE( true) << "pttrs_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_double_buffer_pair_rand( e, eref, (n-1));

   } /* end of Constructor  */

pttrs_double_parameters:: ~pttrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pttrs_double_parameters object: destructor invoked. \n");
#endif
   pttrs_free();
}


//  Test fixture class definition
class dpttrs_test  : public  ::testing::Test {
public:
   pttrs_double_parameters  *dpttrs_obj;
   void SetUp();  
   void TearDown () { delete dpttrs_obj; }
};


void dpttrs_test::SetUp(){

    /* LAPACKE DPTTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dpttrs) ( int matrix_layout, lapack_int n,
                          lapack_int nrhs, const double* d, const double* e,
                              double* b, lapack_int ldb);

    Fptr_NL_LAPACKE_dpttrs DPTTRS;

     /* LAPACKE DPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_dpttrf) (lapack_int n, double* d, double *e);

    Fptr_NL_LAPACKE_dpttrf DPTTRF;

    dpttrs_obj = new  pttrs_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    dpttrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpttrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpttrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpttrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DPTTRS = (Fptr_NL_LAPACKE_dpttrs)dlsym(dpttrs_obj->hModule, "LAPACKE_dpttrs");
    ASSERT_TRUE(DPTTRS != NULL) << "failed to get the Netlib LAPACKE_dpttrs symbol";

    DPTTRF = (Fptr_NL_LAPACKE_dpttrf)dlsym(dpttrs_obj->hModule,"LAPACKE_dpttrf");
    ASSERT_TRUE(DPTTRF != NULL) << "failed to get the Netlib LAPACKE_dpttrf symbol";

    /* Pre condition: need to call pttrf - before calling pttrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    dpttrs_obj->inforef = DPTTRF( dpttrs_obj->n, dpttrs_obj->dref, dpttrs_obj->eref);

    dpttrs_obj->inforef = DPTTRS( dpttrs_obj->matrix_layout,
                                  dpttrs_obj->n,
                                  dpttrs_obj->nrhs,
                                  (const double *)dpttrs_obj->dref,
                                  (const double *)dpttrs_obj->eref,
                                  dpttrs_obj->bref,
                                  dpttrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    dpttrs_obj->info = LAPACKE_dpttrf( dpttrs_obj->n, dpttrs_obj->d,
                                     dpttrs_obj->e);

    dpttrs_obj->info = LAPACKE_dpttrs( dpttrs_obj->matrix_layout,
                                 dpttrs_obj->n, dpttrs_obj->nrhs,
                                  (const double *)dpttrs_obj->d,
                                  (const double *)dpttrs_obj->e,
                                  dpttrs_obj->b, dpttrs_obj->ldb );

    if( dpttrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
                   LAPACKE_dpttrs is wrong\n", dpttrs_obj->info );
    }
    if( dpttrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpttrs \
                            is wrong\n",  dpttrs_obj->inforef );
    }
}

TEST_F(dpttrs_test, dpttrs1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpttrs_obj->b_bufsize,
                           dpttrs_obj->b, dpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpttrs_test, dpttrs2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpttrs_obj->b_bufsize,
                           dpttrs_obj->b, dpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpttrs_test, dpttrs3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpttrs_obj->b_bufsize,
                           dpttrs_obj->b, dpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpttrs_test, dpttrs4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpttrs_obj->b_bufsize,
                           dpttrs_obj->b, dpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pttrs_float_parameters  class definition */
class pttrs_float_parameters{
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
      pttrs_float_parameters ( int matrix_layout_i, lapack_int n_i,
                               lapack_int nrhs_i, lapack_int ldb_i);
             
      ~pttrs_float_parameters ();
};  /* end of pttrs_float_parameters  class definition */


/* Constructor pttrs_float_parameters definition */
pttrs_float_parameters:: pttrs_float_parameters ( int matrix_layout_i, 
               lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pttrs float:  n: %d, trans: %c ldb: %d nrhs: %d \n",
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
       pttrs_free();
       EXPECT_FALSE( true) << "pttrs_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_float_buffer_pair_rand( e, eref, (n-1));

   } /* end of Constructor  */

pttrs_float_parameters:: ~pttrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pttrs_float_parameters object: destructor invoked. \n");
#endif
   pttrs_free();
}


//  Test fixture class definition
class spttrs_test  : public  ::testing::Test {
public:
   pttrs_float_parameters  *spttrs_obj;
   void SetUp();  
   void TearDown () { delete spttrs_obj; }
};


void spttrs_test::SetUp(){

    /* LAPACKE SPTTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_spttrs) ( int matrix_layout, lapack_int n,
                          lapack_int nrhs, const float* d, const float* e,
                              float* b, lapack_int ldb);

    Fptr_NL_LAPACKE_spttrs SPTTRS;

     /* LAPACKE SPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_spttrf) (lapack_int n, float* d, float *e);

    Fptr_NL_LAPACKE_spttrf SPTTRF;

    spttrs_obj = new  pttrs_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    spttrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spttrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spttrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spttrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SPTTRS = (Fptr_NL_LAPACKE_spttrs)dlsym(spttrs_obj->hModule, "LAPACKE_spttrs");
    ASSERT_TRUE(SPTTRS != NULL) << "failed to get the Netlib LAPACKE_spttrs symbol";

    SPTTRF = (Fptr_NL_LAPACKE_spttrf)dlsym(spttrs_obj->hModule,"LAPACKE_spttrf");
    ASSERT_TRUE(SPTTRF != NULL) << "failed to get the Netlib LAPACKE_spttrf symbol";

    /* Pre condition: need to call pttrf - before calling pttrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    spttrs_obj->inforef = SPTTRF( spttrs_obj->n, spttrs_obj->dref, spttrs_obj->eref);

    spttrs_obj->inforef = SPTTRS( spttrs_obj->matrix_layout,
                                  spttrs_obj->n,
                                  spttrs_obj->nrhs,
                                  (const float *)spttrs_obj->dref,
                                  (const float *)spttrs_obj->eref,
                                  spttrs_obj->bref,
                                  spttrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    spttrs_obj->info = LAPACKE_spttrf( spttrs_obj->n, spttrs_obj->d,
                                     spttrs_obj->e);

    spttrs_obj->info = LAPACKE_spttrs( spttrs_obj->matrix_layout,
                                 spttrs_obj->n, spttrs_obj->nrhs,
                                  (const float *)spttrs_obj->d,
                                  (const float *)spttrs_obj->e,
                                  spttrs_obj->b, spttrs_obj->ldb );

    if( spttrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
                   LAPACKE_spttrs is wrong\n", spttrs_obj->info );
    }
    if( spttrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spttrs \
                            is wrong\n",  spttrs_obj->inforef );
    }
}

TEST_F(spttrs_test, spttrs1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spttrs_obj->b_bufsize,
                           spttrs_obj->b, spttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spttrs_test, spttrs2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spttrs_obj->b_bufsize,
                           spttrs_obj->b, spttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spttrs_test, spttrs3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spttrs_obj->b_bufsize,
                           spttrs_obj->b, spttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spttrs_test, spttrs4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spttrs_obj->b_bufsize,
                           spttrs_obj->b, spttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pttrs_scomplex_parameters  class definition */
class pttrs_scomplex_parameters{
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
      pttrs_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
             lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~pttrs_scomplex_parameters ();
};  /* end of pttrs_scomplex_parameters  class definition */


/* Constructor pttrs_scomplex_parameters definition */
pttrs_scomplex_parameters:: pttrs_scomplex_parameters ( int matrix_layout_i, 
         char uplo_i, lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    uplo = uplo_i;
    n = n_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pttrs lapack_complex_float:  n: %d, trans: %c ldb: %d nrhs: %d \n",
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
       pttrs_free();
       EXPECT_FALSE( true) << "pttrs_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_scomplex_buffer_pair_rand( e, eref, (n-1));

   } /* end of Constructor  */

pttrs_scomplex_parameters:: ~pttrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pttrs_scomplex_parameters object: destructor invoked. \n");
#endif
   pttrs_free();
}


//  Test fixture class definition
class cpttrs_test  : public  ::testing::Test {
public:
   pttrs_scomplex_parameters  *cpttrs_obj;
   void SetUp();  
   void TearDown () { delete cpttrs_obj; }
};


void cpttrs_test::SetUp(){

    /* LAPACKE CPTTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_cpttrs) ( int matrix_layout, char uplo,
                           lapack_int n, lapack_int nrhs, const float* d, 
                                           const lapack_complex_float* e,
                                lapack_complex_float* b, lapack_int ldb);

    Fptr_NL_LAPACKE_cpttrs CPTTRS;

     /* LAPACKE CPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_cpttrf) (lapack_int n, 
                     float* d, lapack_complex_float *e);

    Fptr_NL_LAPACKE_cpttrf CPTTRF;

    cpttrs_obj = new  pttrs_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    cpttrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpttrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpttrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpttrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CPTTRS = (Fptr_NL_LAPACKE_cpttrs)dlsym(cpttrs_obj->hModule, "LAPACKE_cpttrs");
    ASSERT_TRUE(CPTTRS != NULL) << "failed to get the Netlib LAPACKE_cpttrs symbol";

    CPTTRF = (Fptr_NL_LAPACKE_cpttrf)dlsym(cpttrs_obj->hModule,"LAPACKE_cpttrf");
    ASSERT_TRUE(CPTTRF != NULL) << "failed to get the Netlib LAPACKE_cpttrf symbol";

    /* Pre condition: need to call pttrf - before calling pttrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    cpttrs_obj->inforef = CPTTRF( cpttrs_obj->n, cpttrs_obj->dref, cpttrs_obj->eref);

    cpttrs_obj->inforef = CPTTRS( cpttrs_obj->matrix_layout,
                                  cpttrs_obj->uplo,
                                  cpttrs_obj->n,
                                  cpttrs_obj->nrhs,
                                  (const float *)cpttrs_obj->dref,
                                  (const lapack_complex_float *)cpttrs_obj->eref,
                                  cpttrs_obj->bref,
                                  cpttrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    cpttrs_obj->info = LAPACKE_cpttrf( cpttrs_obj->n, cpttrs_obj->d,
                                     cpttrs_obj->e);

    cpttrs_obj->info = LAPACKE_cpttrs( cpttrs_obj->matrix_layout,
               cpttrs_obj->uplo, cpttrs_obj->n, cpttrs_obj->nrhs,
                                    (const float *)cpttrs_obj->d,
                     (const lapack_complex_float *)cpttrs_obj->e,
                                  cpttrs_obj->b, cpttrs_obj->ldb );

    if( cpttrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
                   LAPACKE_cpttrs is wrong\n", cpttrs_obj->info );
    }
    if( cpttrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpttrs \
                            is wrong\n",  cpttrs_obj->inforef );
    }
}

TEST_F(cpttrs_test, cpttrs1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpttrs_obj->b_bufsize,
                           cpttrs_obj->b, cpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpttrs_test, cpttrs2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpttrs_obj->b_bufsize,
                           cpttrs_obj->b, cpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpttrs_test, cpttrs3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpttrs_obj->b_bufsize,
                           cpttrs_obj->b, cpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpttrs_test, cpttrs4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpttrs_obj->b_bufsize,
                           cpttrs_obj->b, cpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pttrs_dcomplex_parameters  class definition */
class pttrs_dcomplex_parameters{
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
      pttrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
             lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~pttrs_dcomplex_parameters ();
};  /* end of pttrs_dcomplex_parameters  class definition */


/* Constructor pttrs_dcomplex_parameters definition */
pttrs_dcomplex_parameters:: pttrs_dcomplex_parameters ( int matrix_layout_i, 
         char uplo_i, lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pttrs lapack_complex_double:  n: %d, trans: %c ldb: %d nrhs: %d \n",
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
       pttrs_free();
       EXPECT_FALSE( true) << "pttrs_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( e, eref, (n-1));

   } /* end of Constructor  */

pttrs_dcomplex_parameters:: ~pttrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pttrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   pttrs_free();
}


//  Test fixture class definition
class zpttrs_test  : public  ::testing::Test {
public:
   pttrs_dcomplex_parameters  *zpttrs_obj;
   void SetUp();  
   void TearDown () { delete zpttrs_obj; }
};


void zpttrs_test::SetUp(){

    /* LAPACKE ZPTTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zpttrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs, const double* d, 
                                           const lapack_complex_double* e,
                                lapack_complex_double* b, lapack_int ldb);

    Fptr_NL_LAPACKE_zpttrs ZPTTRS;

     /* LAPACKE ZPTTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zpttrf) (lapack_int n, 
                     double* d, lapack_complex_double *e);

    Fptr_NL_LAPACKE_zpttrf ZPTTRF;

    zpttrs_obj = new  pttrs_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);

    zpttrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpttrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpttrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpttrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZPTTRS = (Fptr_NL_LAPACKE_zpttrs)dlsym(zpttrs_obj->hModule, "LAPACKE_zpttrs");
    ASSERT_TRUE(ZPTTRS != NULL) << "failed to get the Netlib LAPACKE_zpttrs symbol";

    ZPTTRF = (Fptr_NL_LAPACKE_zpttrf)dlsym(zpttrs_obj->hModule,"LAPACKE_zpttrf");
    ASSERT_TRUE(ZPTTRF != NULL) << "failed to get the Netlib LAPACKE_zpttrf symbol";

    /* Pre condition: need to call pttrf - before calling pttrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zpttrs_obj->inforef = ZPTTRF( zpttrs_obj->n, zpttrs_obj->dref, zpttrs_obj->eref);

    zpttrs_obj->inforef = ZPTTRS( zpttrs_obj->matrix_layout,
                                  zpttrs_obj->uplo,
                                  zpttrs_obj->n,
                                  zpttrs_obj->nrhs,
                                  (const double *)zpttrs_obj->dref,
                                  (const lapack_complex_double *)zpttrs_obj->eref,
                                  zpttrs_obj->bref,
                                  zpttrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zpttrs_obj->info = LAPACKE_zpttrf( zpttrs_obj->n, zpttrs_obj->d,
                                     zpttrs_obj->e);

    zpttrs_obj->info = LAPACKE_zpttrs( zpttrs_obj->matrix_layout,
                                               zpttrs_obj->uplo,
                                 zpttrs_obj->n, zpttrs_obj->nrhs,
                                  (const double *)zpttrs_obj->d,
                    (const lapack_complex_double *)zpttrs_obj->e,
                                  zpttrs_obj->b, zpttrs_obj->ldb );

    if( zpttrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
                   LAPACKE_zpttrs is wrong\n", zpttrs_obj->info );
    }
    if( zpttrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpttrs \
                            is wrong\n",  zpttrs_obj->inforef );
    }
}

TEST_F(zpttrs_test, zpttrs1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpttrs_obj->b_bufsize,
                           zpttrs_obj->b, zpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpttrs_test, zpttrs2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpttrs_obj->b_bufsize,
                           zpttrs_obj->b, zpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpttrs_test, zpttrs3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpttrs_obj->b_bufsize,
                           zpttrs_obj->b, zpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpttrs_test, zpttrs4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpttrs_obj->b_bufsize,
                           zpttrs_obj->b, zpttrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}