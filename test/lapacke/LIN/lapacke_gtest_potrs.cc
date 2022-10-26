#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define potrs_free() \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin potrs_double_parameters  class definition */
class potrs_double_parameters{
   public:
      int b_bufsize;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      potrs_double_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~potrs_double_parameters (); 
};  /* end of potrs_double_parameters  class definition */


/* Constructor potrs_double_parameters definition */
potrs_double_parameters:: potrs_double_parameters ( int matrix_layout_i, 
                         char uplo_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n potrs Double:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, ldab, ldb, nrhs);
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

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       potrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

potrs_double_parameters:: ~potrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrs_double_parameters object: destructor invoked. \n");
#endif
   potrs_free();
}

//  Test fixture class definition
class dpotrs_test  : public  ::testing::Test {
public:
   potrs_double_parameters  *dpotrs_obj;
   void SetUp();  
   void TearDown () { delete dpotrs_obj; }
};


void dpotrs_test::SetUp(){

    /* LAPACKE DPOTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dpotrs) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs, const double * a,
                          lapack_int lda, double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dpotrs DPOTRS;

    typedef int (*Fptr_NL_LAPACKE_dpotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, double *a, lapack_int lda );
    Fptr_NL_LAPACKE_dpotrf DPOTRF;

    void *hModule, *dModule;

    dpotrs_obj = new potrs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";


    DPOTRS = (Fptr_NL_LAPACKE_dpotrs)dlsym(hModule, "LAPACKE_dpotrs");
    ASSERT_TRUE(DPOTRS != NULL) << "failed to get the Netlib LAPACKE_dpotrs symbol";
    
    DPOTRF = (Fptr_NL_LAPACKE_dpotrf)dlsym(hModule,"LAPACKE_dpotrf");
    ASSERT_TRUE(DPOTRF != NULL) << "failed to get the Netlib LAPACKE_dpotrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dpotrs_obj->inforef = DPOTRF( dpotrs_obj->matrix_layout,dpotrs_obj->uplo,
                            dpotrs_obj->n,dpotrs_obj->aref,dpotrs_obj->lda);

    dpotrs_obj->inforef = DPOTRS( dpotrs_obj->matrix_layout, dpotrs_obj->uplo, 
                             dpotrs_obj->n, dpotrs_obj->nrhs,
                             (const double *)dpotrs_obj->aref, 
                          dpotrs_obj->lda, dpotrs_obj->bref, dpotrs_obj->ldb);
                          

    /* Compute libflame's Lapacke o/p  */
    dpotrs_obj->info  = LAPACKE_dpotrf( dpotrs_obj->matrix_layout,
                                 dpotrs_obj->uplo,dpotrs_obj->n,
                                 dpotrs_obj->a,dpotrs_obj->lda);

    dpotrs_obj->info = LAPACKE_dpotrs( dpotrs_obj->matrix_layout, 
               dpotrs_obj->uplo, dpotrs_obj->n, dpotrs_obj->nrhs, 
                                   (const double *)dpotrs_obj->a, 
               dpotrs_obj->lda, dpotrs_obj->b, dpotrs_obj->ldb );

    if( dpotrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dpotrs is wrong\n", dpotrs_obj->info );
    }
    if( dpotrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpotrs is wrong\n", 
        dpotrs_obj->inforef );
    }
       if( hModule != NULL) dlclose(hModule);
       if(dModule != NULL) dlclose(dModule);
}

TEST_F(dpotrs_test, dpotrs1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpotrs_obj->b_bufsize, dpotrs_obj->b,
                                            dpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpotrs_test, dpotrs2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpotrs_obj->b_bufsize, dpotrs_obj->b,
                                            dpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpotrs_test, dpotrs3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpotrs_obj->b_bufsize, dpotrs_obj->b,
                                            dpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpotrs_test, dpotrs4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dpotrs_obj->b_bufsize, dpotrs_obj->b,
                                            dpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin potrs_float_parameters  class definition */
class potrs_float_parameters{
   public:
      int b_bufsize;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      potrs_float_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~potrs_float_parameters (); 
};  /* end of potrs_float_parameters  class definition */


/* Constructor potrs_float_parameters definition */
potrs_float_parameters:: potrs_float_parameters ( int matrix_layout_i, 
                         char uplo_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n potrs float:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, ldab, ldb, nrhs);
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
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*lda));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       potrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

potrs_float_parameters:: ~potrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrs_float_parameters object: destructor invoked. \n");
#endif
   potrs_free();
}

//  Test fixture class definition
class spotrs_test  : public  ::testing::Test {
public:
   potrs_float_parameters  *spotrs_obj;
   void SetUp();  
   void TearDown () { delete spotrs_obj; }
};


void spotrs_test::SetUp(){

    /* LAPACKE SPOTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_spotrs) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs, const float * a,
                          lapack_int lda, float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_spotrs SPOTRS;

    typedef int (*Fptr_NL_LAPACKE_spotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, float *a, lapack_int lda );
    Fptr_NL_LAPACKE_spotrf SPOTRF;

    void *hModule, *dModule;

    spotrs_obj = new potrs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";


    SPOTRS = (Fptr_NL_LAPACKE_spotrs)dlsym(hModule, "LAPACKE_spotrs");
    ASSERT_TRUE(SPOTRS != NULL) << "failed to get the Netlib LAPACKE_spotrs symbol";
    
    SPOTRF = (Fptr_NL_LAPACKE_spotrf)dlsym(hModule,"LAPACKE_spotrf");
    ASSERT_TRUE(SPOTRF != NULL) << "failed to get the Netlib LAPACKE_spotrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    spotrs_obj->inforef = SPOTRF( spotrs_obj->matrix_layout,spotrs_obj->uplo,
                            spotrs_obj->n,spotrs_obj->aref,spotrs_obj->lda);

    spotrs_obj->inforef = SPOTRS( spotrs_obj->matrix_layout, spotrs_obj->uplo, 
                             spotrs_obj->n, spotrs_obj->nrhs,
                             (const float *)spotrs_obj->aref, 
                          spotrs_obj->lda, spotrs_obj->bref, spotrs_obj->ldb);
                          

    /* Compute libflame's Lapacke o/p  */
    spotrs_obj->info  = LAPACKE_spotrf( spotrs_obj->matrix_layout,
                                 spotrs_obj->uplo,spotrs_obj->n,
                                 spotrs_obj->a,spotrs_obj->lda);

    spotrs_obj->info = LAPACKE_spotrs( spotrs_obj->matrix_layout, 
               spotrs_obj->uplo, spotrs_obj->n, spotrs_obj->nrhs, 
                                   (const float *)spotrs_obj->a, 
               spotrs_obj->lda, spotrs_obj->b, spotrs_obj->ldb );

    if( spotrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_spotrs is wrong\n", spotrs_obj->info );
    }
    if( spotrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spotrs is wrong\n", 
        spotrs_obj->inforef );
    }
       if( hModule != NULL) dlclose(hModule);
       if(dModule != NULL) dlclose(dModule);
}

TEST_F(spotrs_test, spotrs1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spotrs_obj->b_bufsize, spotrs_obj->b,
                                            spotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spotrs_test, spotrs2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spotrs_obj->b_bufsize, spotrs_obj->b,
                                            spotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spotrs_test, spotrs3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spotrs_obj->b_bufsize, spotrs_obj->b,
                                            spotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spotrs_test, spotrs4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( spotrs_obj->b_bufsize, spotrs_obj->b,
                                            spotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin potrs_scomplex_parameters  class definition */
class potrs_scomplex_parameters{
   public:
      int b_bufsize;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      potrs_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~potrs_scomplex_parameters (); 
};  /* end of potrs_scomplex_parameters  class definition */


/* Constructor potrs_scomplex_parameters definition */
potrs_scomplex_parameters:: potrs_scomplex_parameters ( int matrix_layout_i, 
                         char uplo_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n potrs lapack_complex_float:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, ldab, ldb, nrhs);
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

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       potrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

potrs_scomplex_parameters:: ~potrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrs_scomplex_parameters object: destructor invoked. \n");
#endif
   potrs_free();
}

//  Test fixture class definition
class cpotrs_test  : public  ::testing::Test {
public:
   potrs_scomplex_parameters  *cpotrs_obj;
   void SetUp();  
   void TearDown () { delete cpotrs_obj; }
};


void cpotrs_test::SetUp(){

    /* LAPACKE CPOTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_cpotrs) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs, const lapack_complex_float * a,
                          lapack_int lda, lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_cpotrs CPOTRS;

    typedef int (*Fptr_NL_LAPACKE_cpotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda );
    Fptr_NL_LAPACKE_cpotrf CPOTRF;

    void *hModule, *dModule;

    cpotrs_obj = new potrs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";


    CPOTRS = (Fptr_NL_LAPACKE_cpotrs)dlsym(hModule, "LAPACKE_cpotrs");
    ASSERT_TRUE(CPOTRS != NULL) << "failed to get the Netlib LAPACKE_cpotrs symbol";
    
    CPOTRF = (Fptr_NL_LAPACKE_cpotrf)dlsym(hModule,"LAPACKE_cpotrf");
    ASSERT_TRUE(CPOTRF != NULL) << "failed to get the Netlib LAPACKE_cpotrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    cpotrs_obj->inforef = CPOTRF( cpotrs_obj->matrix_layout,cpotrs_obj->uplo,
                            cpotrs_obj->n,cpotrs_obj->aref,cpotrs_obj->lda);

    cpotrs_obj->inforef = CPOTRS( cpotrs_obj->matrix_layout, cpotrs_obj->uplo, 
                             cpotrs_obj->n, cpotrs_obj->nrhs,
                             (const lapack_complex_float *)cpotrs_obj->aref, 
                          cpotrs_obj->lda, cpotrs_obj->bref, cpotrs_obj->ldb);
                          

    /* Compute libflame's Lapacke o/p  */
    cpotrs_obj->info  = LAPACKE_cpotrf( cpotrs_obj->matrix_layout,
                                 cpotrs_obj->uplo,cpotrs_obj->n,
                                 cpotrs_obj->a,cpotrs_obj->lda);

    cpotrs_obj->info = LAPACKE_cpotrs( cpotrs_obj->matrix_layout, 
               cpotrs_obj->uplo, cpotrs_obj->n, cpotrs_obj->nrhs, 
                                   (const lapack_complex_float *)cpotrs_obj->a, 
               cpotrs_obj->lda, cpotrs_obj->b, cpotrs_obj->ldb );

    if( cpotrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpotrs is wrong\n", cpotrs_obj->info );
    }
    if( cpotrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpotrs is wrong\n", 
        cpotrs_obj->inforef );
    }
       if( hModule != NULL) dlclose(hModule);
       if(dModule != NULL) dlclose(dModule);
}

TEST_F(cpotrs_test, cpotrs1) {
    float diff;    
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpotrs_obj->b_bufsize, cpotrs_obj->b,
                                            cpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpotrs_test, cpotrs2) {
    float diff;    
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpotrs_obj->b_bufsize, cpotrs_obj->b,
                                            cpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpotrs_test, cpotrs3) {
    float diff;    
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpotrs_obj->b_bufsize, cpotrs_obj->b,
                                            cpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpotrs_test, cpotrs4) {
    float diff;    
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cpotrs_obj->b_bufsize, cpotrs_obj->b,
                                            cpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin potrs_dcomplex_parameters  class definition */
class potrs_dcomplex_parameters{
   public:
      int b_bufsize;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      potrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~potrs_dcomplex_parameters (); 
};  /* end of potrs_dcomplex_parameters  class definition */


/* Constructor potrs_dcomplex_parameters definition */
potrs_dcomplex_parameters:: potrs_dcomplex_parameters ( int matrix_layout_i, 
                         char uplo_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n potrs lapack_complex_double:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
             n, uplo, ldab, ldb, nrhs);
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

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       potrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

potrs_dcomplex_parameters:: ~potrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   potrs_free();
}

//  Test fixture class definition
class zpotrs_test  : public  ::testing::Test {
public:
   potrs_dcomplex_parameters  *zpotrs_obj;
   void SetUp();  
   void TearDown () { delete zpotrs_obj; }
};


void zpotrs_test::SetUp(){

    /* LAPACKE ZPOTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zpotrs) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs, const lapack_complex_double * a,
                          lapack_int lda, lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zpotrs ZPOTRS;

    typedef int (*Fptr_NL_LAPACKE_zpotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda );
    Fptr_NL_LAPACKE_zpotrf ZPOTRF;

    void *hModule, *dModule;

    zpotrs_obj = new potrs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(hModule != NULL) << "Netlib lapacke handle NULL";


    ZPOTRS = (Fptr_NL_LAPACKE_zpotrs)dlsym(hModule, "LAPACKE_zpotrs");
    ASSERT_TRUE(ZPOTRS != NULL) << "failed to get the Netlib LAPACKE_zpotrs symbol";
    
    ZPOTRF = (Fptr_NL_LAPACKE_zpotrf)dlsym(hModule,"LAPACKE_zpotrf");
    ASSERT_TRUE(ZPOTRF != NULL) << "failed to get the Netlib LAPACKE_zpotrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zpotrs_obj->inforef = ZPOTRF( zpotrs_obj->matrix_layout,zpotrs_obj->uplo,
                            zpotrs_obj->n,zpotrs_obj->aref,zpotrs_obj->lda);

    zpotrs_obj->inforef = ZPOTRS( zpotrs_obj->matrix_layout, zpotrs_obj->uplo, 
                             zpotrs_obj->n, zpotrs_obj->nrhs,
                             (const lapack_complex_double *)zpotrs_obj->aref, 
                          zpotrs_obj->lda, zpotrs_obj->bref, zpotrs_obj->ldb);
                          

    /* Compute libflame's Lapacke o/p  */
    zpotrs_obj->info  = LAPACKE_zpotrf( zpotrs_obj->matrix_layout,
                                 zpotrs_obj->uplo,zpotrs_obj->n,
                                 zpotrs_obj->a,zpotrs_obj->lda);

    zpotrs_obj->info = LAPACKE_zpotrs( zpotrs_obj->matrix_layout, 
               zpotrs_obj->uplo, zpotrs_obj->n, zpotrs_obj->nrhs, 
                    (const lapack_complex_double *)zpotrs_obj->a, 
               zpotrs_obj->lda, zpotrs_obj->b, zpotrs_obj->ldb );

    if( zpotrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zpotrs is wrong\n", zpotrs_obj->info );
    }
    if( zpotrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpotrs is wrong\n", 
        zpotrs_obj->inforef );
    }
       if( hModule != NULL) dlclose(hModule);
       if(dModule != NULL) dlclose(dModule);
}

TEST_F(zpotrs_test, zpotrs1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpotrs_obj->b_bufsize, zpotrs_obj->b,
                                            zpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpotrs_test, zpotrs2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpotrs_obj->b_bufsize, zpotrs_obj->b,
                                            zpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpotrs_test, zpotrs3) {
    double diff;        
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpotrs_obj->b_bufsize, zpotrs_obj->b,
                                            zpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpotrs_test, zpotrs4) {
    double diff;        
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zpotrs_obj->b_bufsize, zpotrs_obj->b,
                                            zpotrs_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
