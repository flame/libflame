#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define trtrs_free() \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin trtrs_double_parameters  class definition */
class trtrs_double_parameters{
   public:
      int b_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
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
      trtrs_double_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i);
              
      ~trtrs_double_parameters (); 
};  /* end of trtrs_double_parameters  class definition */


/* Constructor trtrs_double_parameters definition */
trtrs_double_parameters:: trtrs_double_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i, lapack_int lda_i, 
                                 lapack_int nrhs_i, lapack_int ldb_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n trtrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
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
       trtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
//    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

trtrs_double_parameters:: ~trtrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trtrs_double_parameters object: destructor invoked. \n");
#endif
   trtrs_free();
}

//  Test fixture class definition
class dtrtrs_test  : public  ::testing::Test {
public:
   trtrs_double_parameters  *dtrtrs_obj;
   void SetUp();  
   void TearDown () { delete dtrtrs_obj; }
};


void dtrtrs_test::SetUp(){

    /* LAPACKE DTRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dtrtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const double *a, 
                                    lapack_int lda, double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_dtrtrs DTRTRS;

    dtrtrs_obj = new trtrs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    dtrtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtrtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtrtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtrtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DTRTRS = (Fptr_NL_LAPACKE_dtrtrs)dlsym(dtrtrs_obj->hModule, "LAPACKE_dtrtrs");
    ASSERT_TRUE(DTRTRS != NULL) << "failed to get the Netlib LAPACKE_dtrtrs symbol";
    

    dtrtrs_obj->inforef = DTRTRS( dtrtrs_obj->matrix_layout, dtrtrs_obj->uplo,
                                  dtrtrs_obj->trans, dtrtrs_obj->diag, 
                                  dtrtrs_obj->n, dtrtrs_obj->nrhs,
                                  (const double *)dtrtrs_obj->aref, 
                                  dtrtrs_obj->lda, dtrtrs_obj->bref,
                                  dtrtrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    dtrtrs_obj->info = LAPACKE_dtrtrs( dtrtrs_obj->matrix_layout, dtrtrs_obj->uplo,
                                  dtrtrs_obj->trans, dtrtrs_obj->diag, 
                                  dtrtrs_obj->n, dtrtrs_obj->nrhs,
                                  (const double *)dtrtrs_obj->a, 
                                  dtrtrs_obj->lda, dtrtrs_obj->b,
                                  dtrtrs_obj->ldb );

    if( dtrtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtrtrs is wrong\n", dtrtrs_obj->info );
    }
    if( dtrtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtrtrs is wrong\n", 
        dtrtrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtrtrs_obj->diff =  computeDiff_d( dtrtrs_obj->b_bufsize, 
                dtrtrs_obj->b, dtrtrs_obj->bref );

}

TEST_F(dtrtrs_test, dtrtrs1) {
    EXPECT_NEAR(0.0, dtrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrtrs_test, dtrtrs2) {
    EXPECT_NEAR(0.0, dtrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrtrs_test, dtrtrs3) {
    EXPECT_NEAR(0.0, dtrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrtrs_test, dtrtrs4) {
    EXPECT_NEAR(0.0, dtrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}



/* Begin trtrs_float_parameters  class definition */
class trtrs_float_parameters{
   public:
      int b_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
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
      trtrs_float_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i);
              
      ~trtrs_float_parameters (); 
};  /* end of trtrs_float_parameters  class definition */


/* Constructor trtrs_float_parameters definition */
trtrs_float_parameters:: trtrs_float_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i, lapack_int lda_i, 
                                 lapack_int nrhs_i, lapack_int ldb_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n trtrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
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
       trtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

trtrs_float_parameters:: ~trtrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trtrs_float_parameters object: destructor invoked. \n");
#endif
   trtrs_free();
}

//  Test fixture class definition
class strtrs_test  : public  ::testing::Test {
public:
   trtrs_float_parameters  *strtrs_obj;
   void SetUp();  
   void TearDown () { delete strtrs_obj; }
};


void strtrs_test::SetUp(){

    /* LAPACKE STRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_strtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const float *a, 
                                    lapack_int lda, float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_strtrs STRTRS;

    strtrs_obj = new trtrs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    strtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    strtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(strtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(strtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    STRTRS = (Fptr_NL_LAPACKE_strtrs)dlsym(strtrs_obj->hModule, "LAPACKE_strtrs");
    ASSERT_TRUE(STRTRS != NULL) << "failed to get the Netlib LAPACKE_strtrs symbol";
    

    strtrs_obj->inforef = STRTRS( strtrs_obj->matrix_layout, strtrs_obj->uplo,
                                  strtrs_obj->trans, strtrs_obj->diag, 
                                  strtrs_obj->n, strtrs_obj->nrhs,
                                  (const float *)strtrs_obj->aref, 
                                  strtrs_obj->lda, strtrs_obj->bref,
                                  strtrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    strtrs_obj->info = LAPACKE_strtrs( strtrs_obj->matrix_layout, strtrs_obj->uplo,
                                  strtrs_obj->trans, strtrs_obj->diag, 
                                  strtrs_obj->n, strtrs_obj->nrhs,
                                  (const float *)strtrs_obj->a, 
                                  strtrs_obj->lda, strtrs_obj->b,
                                  strtrs_obj->ldb );

    if( strtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_strtrs is wrong\n", strtrs_obj->info );
    }
    if( strtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_strtrs is wrong\n", 
        strtrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    strtrs_obj->diff =  computeDiff_s( strtrs_obj->b_bufsize, 
                strtrs_obj->b, strtrs_obj->bref );

}

TEST_F(strtrs_test, strtrs1) {
    EXPECT_NEAR(0.0, strtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strtrs_test, strtrs2) {
    EXPECT_NEAR(0.0, strtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strtrs_test, strtrs3) {
    EXPECT_NEAR(0.0, strtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strtrs_test, strtrs4) {
    EXPECT_NEAR(0.0, strtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin trtrs_scomplex_parameters  class definition */
class trtrs_scomplex_parameters{
   public:
      int b_bufsize;
      float diff; // capture difference between ref o/p & libflame lapacke o/p.
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
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
      trtrs_scomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i);
              
      ~trtrs_scomplex_parameters (); 
};  /* end of trtrs_scomplex_parameters  class definition */


/* Constructor trtrs_scomplex_parameters definition */
trtrs_scomplex_parameters:: trtrs_scomplex_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i, lapack_int lda_i, 
                                 lapack_int nrhs_i, lapack_int ldb_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n trtrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
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
       trtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

trtrs_scomplex_parameters:: ~trtrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trtrs_scomplex_parameters object: destructor invoked. \n");
#endif
   trtrs_free();
}

//  Test fixture class definition
class ctrtrs_test  : public  ::testing::Test {
public:
   trtrs_scomplex_parameters  *ctrtrs_obj;
   void SetUp();  
   void TearDown () { delete ctrtrs_obj; }
};


void ctrtrs_test::SetUp(){

    /* LAPACKE CTRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ctrtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const lapack_complex_float *a, 
                                    lapack_int lda, lapack_complex_float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ctrtrs CTRTRS;

    ctrtrs_obj = new trtrs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    ctrtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctrtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctrtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctrtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CTRTRS = (Fptr_NL_LAPACKE_ctrtrs)dlsym(ctrtrs_obj->hModule, "LAPACKE_ctrtrs");
    ASSERT_TRUE(CTRTRS != NULL) << "failed to get the Netlib LAPACKE_ctrtrs symbol";
    

    ctrtrs_obj->inforef = CTRTRS( ctrtrs_obj->matrix_layout, ctrtrs_obj->uplo,
                                  ctrtrs_obj->trans, ctrtrs_obj->diag, 
                                  ctrtrs_obj->n, ctrtrs_obj->nrhs,
                                  (const lapack_complex_float *)ctrtrs_obj->aref, 
                                  ctrtrs_obj->lda, ctrtrs_obj->bref,
                                  ctrtrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    ctrtrs_obj->info = LAPACKE_ctrtrs( ctrtrs_obj->matrix_layout, ctrtrs_obj->uplo,
                                  ctrtrs_obj->trans, ctrtrs_obj->diag, 
                                  ctrtrs_obj->n, ctrtrs_obj->nrhs,
                                  (const lapack_complex_float *)ctrtrs_obj->a, 
                                  ctrtrs_obj->lda, ctrtrs_obj->b,
                                  ctrtrs_obj->ldb );

    if( ctrtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctrtrs is wrong\n", ctrtrs_obj->info );
    }
    if( ctrtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctrtrs is wrong\n", 
        ctrtrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctrtrs_obj->diff =  computeDiff_c( ctrtrs_obj->b_bufsize, 
                ctrtrs_obj->b, ctrtrs_obj->bref );

}

TEST_F(ctrtrs_test, ctrtrs1) {
    EXPECT_NEAR(0.0, ctrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrtrs_test, ctrtrs2) {
    EXPECT_NEAR(0.0, ctrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrtrs_test, ctrtrs3) {
    EXPECT_NEAR(0.0, ctrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrtrs_test, ctrtrs4) {
    EXPECT_NEAR(0.0, ctrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin trtrs_dcomplex_parameters  class definition */
class trtrs_dcomplex_parameters{
   public:
      int b_bufsize;
      double diff; // capture difference between ref o/p & libflame lapacke o/p.
      void *hModule, *dModule;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char trans; //  Must be 'N' , 'T' or 'C'.
      char diag; //  Must be 'U' or 'N'.
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
      trtrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i);
              
      ~trtrs_dcomplex_parameters (); 
};  /* end of trtrs_dcomplex_parameters  class definition */


/* Constructor trtrs_dcomplex_parameters definition */
trtrs_dcomplex_parameters:: trtrs_dcomplex_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i, lapack_int lda_i, 
                                 lapack_int nrhs_i, lapack_int ldb_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n trtrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
    ldb: %d nrhs: %d \n",  n, uplo, trans, diag, lda, ldb, nrhs);
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
       trtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

trtrs_dcomplex_parameters:: ~trtrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trtrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   trtrs_free();
}

//  Test fixture class definition
class ztrtrs_test  : public  ::testing::Test {
public:
   trtrs_dcomplex_parameters  *ztrtrs_obj;
   void SetUp();  
   void TearDown () { delete ztrtrs_obj; }
};


void ztrtrs_test::SetUp(){

    /* LAPACKE ZTRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ztrtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const lapack_complex_double *a, 
                                    lapack_int lda, lapack_complex_double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ztrtrs ZTRTRS;

    ztrtrs_obj = new trtrs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    ztrtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztrtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztrtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztrtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZTRTRS = (Fptr_NL_LAPACKE_ztrtrs)dlsym(ztrtrs_obj->hModule, "LAPACKE_ztrtrs");
    ASSERT_TRUE(ZTRTRS != NULL) << "failed to get the Netlib LAPACKE_ztrtrs symbol";
    

    ztrtrs_obj->inforef = ZTRTRS( ztrtrs_obj->matrix_layout, ztrtrs_obj->uplo,
                                  ztrtrs_obj->trans, ztrtrs_obj->diag, 
                                  ztrtrs_obj->n, ztrtrs_obj->nrhs,
                                  (const lapack_complex_double *)ztrtrs_obj->aref, 
                                  ztrtrs_obj->lda, ztrtrs_obj->bref,
                                  ztrtrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    ztrtrs_obj->info = LAPACKE_ztrtrs( ztrtrs_obj->matrix_layout, ztrtrs_obj->uplo,
                                  ztrtrs_obj->trans, ztrtrs_obj->diag, 
                                  ztrtrs_obj->n, ztrtrs_obj->nrhs,
                                  (const lapack_complex_double *)ztrtrs_obj->a, 
                                  ztrtrs_obj->lda, ztrtrs_obj->b,
                                  ztrtrs_obj->ldb );

    if( ztrtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztrtrs is wrong\n", ztrtrs_obj->info );
    }
    if( ztrtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztrtrs is wrong\n", 
        ztrtrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztrtrs_obj->diff =  computeDiff_z( ztrtrs_obj->b_bufsize, 
                ztrtrs_obj->b, ztrtrs_obj->bref );

}

TEST_F(ztrtrs_test, ztrtrs1) {
    EXPECT_NEAR(0.0, ztrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrtrs_test, ztrtrs2) {
    EXPECT_NEAR(0.0, ztrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrtrs_test, ztrtrs3) {
    EXPECT_NEAR(0.0, ztrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrtrs_test, ztrtrs4) {
    EXPECT_NEAR(0.0, ztrtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

