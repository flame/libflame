#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define posv_free() \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin posv_double_parameters  class definition */
class posv_double_parameters{
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
      posv_double_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~posv_double_parameters (); 
};  /* end of posv_double_parameters  class definition */


/* Constructor posv_double_parameters definition */
posv_double_parameters:: posv_double_parameters ( int matrix_layout_i, 
                         char uplo_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n posv Double:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       posv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

posv_double_parameters:: ~posv_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" posv_double_parameters object: destructor invoked. \n");
#endif
   posv_free();
}

//  Test fixture class definition
class dposv_test  : public  ::testing::Test {
public:
   posv_double_parameters  *dposv_obj;
   void SetUp();  
   void TearDown () { delete dposv_obj; }
};


void dposv_test::SetUp(){

    /* LAPACKE DPOSV prototype */
    typedef int (*Fptr_NL_LAPACKE_dposv) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs,  double * a,
                          lapack_int lda, double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_dposv DPOSV;

    typedef int (*Fptr_NL_LAPACKE_dpotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, double *a, lapack_int lda );
    Fptr_NL_LAPACKE_dpotrf DPOTRF;

    void *hModule, *dModule;

    dposv_obj = new posv_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
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


    DPOSV = (Fptr_NL_LAPACKE_dposv)dlsym(hModule, "LAPACKE_dposv");
    ASSERT_TRUE(DPOSV != NULL) << "failed to get the Netlib LAPACKE_dposv symbol";
    
    DPOTRF = (Fptr_NL_LAPACKE_dpotrf)dlsym(hModule,"LAPACKE_dpotrf");
    ASSERT_TRUE(DPOTRF != NULL) << "failed to get the Netlib LAPACKE_dpotrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the reference o/p by invoking Netlib-LapackE's API */
#if 0
    dposv_obj->inforef = DPOTRF( dposv_obj->matrix_layout,dposv_obj->uplo,
                            dposv_obj->n,dposv_obj->aref,dposv_obj->lda);
#endif
    dposv_obj->inforef = DPOSV( dposv_obj->matrix_layout, dposv_obj->uplo, 
                             dposv_obj->n, dposv_obj->nrhs,
                             ( double *)dposv_obj->aref, 
                          dposv_obj->lda, dposv_obj->bref, dposv_obj->ldb);
                          

    /* Compute libflame's Lapacke o/p  */
#if 0
    dposv_obj->info  = LAPACKE_dpotrf( dposv_obj->matrix_layout,
                                 dposv_obj->uplo,dposv_obj->n,
                                 dposv_obj->a,dposv_obj->lda);
#endif
    dposv_obj->info = LAPACKE_dposv( dposv_obj->matrix_layout, 
               dposv_obj->uplo, dposv_obj->n, dposv_obj->nrhs, 
                                   ( double *)dposv_obj->a, 
               dposv_obj->lda, dposv_obj->b, dposv_obj->ldb );

    if( dposv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dposv is wrong\n", dposv_obj->info );
    }
    if( dposv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dposv is wrong\n", 
        dposv_obj->inforef );
    }
       if( hModule != NULL) dlclose(hModule);
       if(dModule != NULL) dlclose(dModule);
}

TEST_F(dposv_test, dposv1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dposv_obj->b_bufsize, dposv_obj->b,
                                            dposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dposv_test, dposv2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dposv_obj->b_bufsize, dposv_obj->b,
                                            dposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dposv_test, dposv3) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dposv_obj->b_bufsize, dposv_obj->b,
                                            dposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dposv_test, dposv4) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_d( dposv_obj->b_bufsize, dposv_obj->b,
                                            dposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin posv_float_parameters  class definition */
class posv_float_parameters{
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
      posv_float_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~posv_float_parameters (); 
};  /* end of posv_float_parameters  class definition */


/* Constructor posv_float_parameters definition */
posv_float_parameters:: posv_float_parameters ( int matrix_layout_i, 
                         char uplo_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n posv float:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       posv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

posv_float_parameters:: ~posv_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" posv_float_parameters object: destructor invoked. \n");
#endif
   posv_free();
}

//  Test fixture class definition
class sposv_test  : public  ::testing::Test {
public:
   posv_float_parameters  *sposv_obj;
   void SetUp();  
   void TearDown () { delete sposv_obj; }
};


void sposv_test::SetUp(){

    /* LAPACKE SPOSV prototype */
    typedef int (*Fptr_NL_LAPACKE_sposv) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs,  float * a,
                          lapack_int lda, float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_sposv SPOSV;

    typedef int (*Fptr_NL_LAPACKE_spotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, float *a, lapack_int lda );
    Fptr_NL_LAPACKE_spotrf SPOTRF;

    void *hModule, *dModule;

    sposv_obj = new posv_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
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


    SPOSV = (Fptr_NL_LAPACKE_sposv)dlsym(hModule, "LAPACKE_sposv");
    ASSERT_TRUE(SPOSV != NULL) << "failed to get the Netlib LAPACKE_sposv symbol";
    
    SPOTRF = (Fptr_NL_LAPACKE_spotrf)dlsym(hModule,"LAPACKE_spotrf");
    ASSERT_TRUE(SPOTRF != NULL) << "failed to get the Netlib LAPACKE_spotrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
#if 0
    sposv_obj->inforef = SPOTRF( sposv_obj->matrix_layout,sposv_obj->uplo,
                            sposv_obj->n,sposv_obj->aref,sposv_obj->lda);
#endif
    sposv_obj->inforef = SPOSV( sposv_obj->matrix_layout, sposv_obj->uplo, 
                             sposv_obj->n, sposv_obj->nrhs,
                             ( float *)sposv_obj->aref, 
                          sposv_obj->lda, sposv_obj->bref, sposv_obj->ldb);
                          

    /* Compute libflame's Lapacke o/p  */
#if 0
    sposv_obj->info  = LAPACKE_spotrf( sposv_obj->matrix_layout,
                                 sposv_obj->uplo,sposv_obj->n,
                                 sposv_obj->a,sposv_obj->lda);
#endif

    sposv_obj->info = LAPACKE_sposv( sposv_obj->matrix_layout, 
               sposv_obj->uplo, sposv_obj->n, sposv_obj->nrhs, 
                                   ( float *)sposv_obj->a, 
               sposv_obj->lda, sposv_obj->b, sposv_obj->ldb );

    if( sposv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sposv is wrong\n", sposv_obj->info );
    }
    if( sposv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sposv is wrong\n", 
        sposv_obj->inforef );
    }
       if( hModule != NULL) dlclose(hModule);
       if(dModule != NULL) dlclose(dModule);
}

TEST_F(sposv_test, sposv1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sposv_obj->b_bufsize, sposv_obj->b,
                                            sposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sposv_test, sposv2) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sposv_obj->b_bufsize, sposv_obj->b,
                                            sposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sposv_test, sposv3) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sposv_obj->b_bufsize, sposv_obj->b,
                                            sposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sposv_test, sposv4) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_s( sposv_obj->b_bufsize, sposv_obj->b,
                                            sposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin posv_scomplex_parameters  class definition */
class posv_scomplex_parameters{
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
      posv_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~posv_scomplex_parameters (); 
};  /* end of posv_scomplex_parameters  class definition */


/* Constructor posv_scomplex_parameters definition */
posv_scomplex_parameters:: posv_scomplex_parameters ( int matrix_layout_i, 
                         char uplo_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n posv lapack_complex_float:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       posv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

posv_scomplex_parameters:: ~posv_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" posv_scomplex_parameters object: destructor invoked. \n");
#endif
   posv_free();
}

//  Test fixture class definition
class cposv_test  : public  ::testing::Test {
public:
   posv_scomplex_parameters  *cposv_obj;
   void SetUp();  
   void TearDown () { delete cposv_obj; }
};


void cposv_test::SetUp(){

    /* LAPACKE CPOSV prototype */
    typedef int (*Fptr_NL_LAPACKE_cposv) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs,  lapack_complex_float * a,
                          lapack_int lda, lapack_complex_float * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_cposv CPOSV;

    typedef int (*Fptr_NL_LAPACKE_cpotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda );
    Fptr_NL_LAPACKE_cpotrf CPOTRF;

    void *hModule, *dModule;

    cposv_obj = new posv_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
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


    CPOSV = (Fptr_NL_LAPACKE_cposv)dlsym(hModule, "LAPACKE_cposv");
    ASSERT_TRUE(CPOSV != NULL) << "failed to get the Netlib LAPACKE_cposv symbol";
    
    CPOTRF = (Fptr_NL_LAPACKE_cpotrf)dlsym(hModule,"LAPACKE_cpotrf");
    ASSERT_TRUE(CPOTRF != NULL) << "failed to get the Netlib LAPACKE_cpotrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
#if 0
    cposv_obj->inforef = CPOTRF( cposv_obj->matrix_layout,cposv_obj->uplo,
                            cposv_obj->n,cposv_obj->aref,cposv_obj->lda);
#endif
    cposv_obj->inforef = CPOSV( cposv_obj->matrix_layout, cposv_obj->uplo, 
                             cposv_obj->n, cposv_obj->nrhs,
                             ( lapack_complex_float *)cposv_obj->aref, 
                          cposv_obj->lda, cposv_obj->bref, cposv_obj->ldb);
                          

    /* Compute libflame's Lapacke o/p  */
#if 0
    cposv_obj->info  = LAPACKE_cpotrf( cposv_obj->matrix_layout,
                                 cposv_obj->uplo,cposv_obj->n,
                                 cposv_obj->a,cposv_obj->lda);
#endif
    cposv_obj->info = LAPACKE_cposv( cposv_obj->matrix_layout, 
               cposv_obj->uplo, cposv_obj->n, cposv_obj->nrhs, 
                                   ( lapack_complex_float *)cposv_obj->a, 
               cposv_obj->lda, cposv_obj->b, cposv_obj->ldb );

    if( cposv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cposv is wrong\n", cposv_obj->info );
    }
    if( cposv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cposv is wrong\n", 
        cposv_obj->inforef );
    }
       if( hModule != NULL) dlclose(hModule);
       if(dModule != NULL) dlclose(dModule);
}

TEST_F(cposv_test, cposv1) {
    float diff;    
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cposv_obj->b_bufsize, cposv_obj->b,
                                            cposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cposv_test, cposv2) {
    float diff;    
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cposv_obj->b_bufsize, cposv_obj->b,
                                            cposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cposv_test, cposv3) {
    float diff;    
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cposv_obj->b_bufsize, cposv_obj->b,
                                            cposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cposv_test, cposv4) {
    float diff;    
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( cposv_obj->b_bufsize, cposv_obj->b,
                                            cposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin posv_dcomplex_parameters  class definition */
class posv_dcomplex_parameters{
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
      posv_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
              
      ~posv_dcomplex_parameters (); 
};  /* end of posv_dcomplex_parameters  class definition */


/* Constructor posv_dcomplex_parameters definition */
posv_dcomplex_parameters:: posv_dcomplex_parameters ( int matrix_layout_i, 
                         char uplo_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n posv lapack_complex_double:  n: %d, Uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       posv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

posv_dcomplex_parameters:: ~posv_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" posv_dcomplex_parameters object: destructor invoked. \n");
#endif
   posv_free();
}

//  Test fixture class definition
class zposv_test  : public  ::testing::Test {
public:
   posv_dcomplex_parameters  *zposv_obj;
   void SetUp();  
   void TearDown () { delete zposv_obj; }
};


void zposv_test::SetUp(){

    /* LAPACKE ZPOSV prototype */
    typedef int (*Fptr_NL_LAPACKE_zposv) (int matrix_layout, char uplo,
                        lapack_int n, lapack_int nrhs,  lapack_complex_double * a,
                          lapack_int lda, lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zposv ZPOSV;

    typedef int (*Fptr_NL_LAPACKE_zpotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda );
    Fptr_NL_LAPACKE_zpotrf ZPOTRF;

    void *hModule, *dModule;

    zposv_obj = new posv_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
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


    ZPOSV = (Fptr_NL_LAPACKE_zposv)dlsym(hModule, "LAPACKE_zposv");
    ASSERT_TRUE(ZPOSV != NULL) << "failed to get the Netlib LAPACKE_zposv symbol";
    
    ZPOTRF = (Fptr_NL_LAPACKE_zpotrf)dlsym(hModule,"LAPACKE_zpotrf");
    ASSERT_TRUE(ZPOTRF != NULL) << "failed to get the Netlib LAPACKE_zpotrf symbol";

    /* Pre condition: need to call pbtrf - before calling pbtrs function */

    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
#if 0
    zposv_obj->inforef = ZPOTRF( zposv_obj->matrix_layout,zposv_obj->uplo,
                            zposv_obj->n,zposv_obj->aref,zposv_obj->lda);
#endif
    zposv_obj->inforef = ZPOSV( zposv_obj->matrix_layout, zposv_obj->uplo, 
                             zposv_obj->n, zposv_obj->nrhs,
                             ( lapack_complex_double *)zposv_obj->aref, 
                          zposv_obj->lda, zposv_obj->bref, zposv_obj->ldb);
                          

    /* Compute libflame's Lapacke o/p  */
#if 0
    zposv_obj->info  = LAPACKE_zpotrf( zposv_obj->matrix_layout,
                                 zposv_obj->uplo,zposv_obj->n,
                                 zposv_obj->a,zposv_obj->lda);
#endif
    zposv_obj->info = LAPACKE_zposv( zposv_obj->matrix_layout, 
               zposv_obj->uplo, zposv_obj->n, zposv_obj->nrhs, 
                    ( lapack_complex_double *)zposv_obj->a, 
               zposv_obj->lda, zposv_obj->b, zposv_obj->ldb );

    if( zposv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zposv is wrong\n", zposv_obj->info );
    }
    if( zposv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zposv is wrong\n", 
        zposv_obj->inforef );
    }
       if( hModule != NULL) dlclose(hModule);
       if(dModule != NULL) dlclose(dModule);
}

TEST_F(zposv_test, zposv1) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zposv_obj->b_bufsize, zposv_obj->b,
                                            zposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zposv_test, zposv2) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zposv_obj->b_bufsize, zposv_obj->b,
                                            zposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zposv_test, zposv3) {
    double diff;        
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zposv_obj->b_bufsize, zposv_obj->b,
                                            zposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zposv_test, zposv4) {
    double diff;        
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zposv_obj->b_bufsize, zposv_obj->b,
                                            zposv_obj->bref );
    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
