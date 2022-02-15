#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define tptrs_free() \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin tptrs_double_parameters  class definition */
class tptrs_double_parameters{
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
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tptrs_double_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i,  
                                lapack_int nrhs_i, lapack_int ldb_i);
              
      ~tptrs_double_parameters (); 
};  /* end of tptrs_double_parameters  class definition */


/* Constructor tptrs_double_parameters definition */
tptrs_double_parameters:: tptrs_double_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i, 
                                 lapack_int nrhs_i, lapack_int ldb_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tptrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
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
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*(n+1)/2));
    lapacke_gtest_alloc_double_buffer_pair( &b, &bref, b_bufsize);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       tptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, n, (n+1)/2, uplo);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

tptrs_double_parameters:: ~tptrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tptrs_double_parameters object: destructor invoked. \n");
#endif
   tptrs_free();
}

//  Test fixture class definition
class dtptrs_test  : public  ::testing::Test {
public:
   tptrs_double_parameters  *dtptrs_obj;
   void SetUp();  
   void TearDown () { delete dtptrs_obj; }
};


void dtptrs_test::SetUp(){

    /* LAPACKE DTPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dtptrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const double *a, 
                                    double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_dtptrs DTPTRS;

    dtptrs_obj = new tptrs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    dtptrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtptrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtptrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtptrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DTPTRS = (Fptr_NL_LAPACKE_dtptrs)dlsym(dtptrs_obj->hModule, "LAPACKE_dtptrs");
    ASSERT_TRUE(DTPTRS != NULL) << "failed to get the Netlib LAPACKE_dtptrs symbol";
    

    dtptrs_obj->inforef = DTPTRS( dtptrs_obj->matrix_layout, dtptrs_obj->uplo,
                                  dtptrs_obj->trans, dtptrs_obj->diag, 
                                  dtptrs_obj->n, dtptrs_obj->nrhs,
                                  (const double *)dtptrs_obj->aref, 
                                  dtptrs_obj->bref,
                                  dtptrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    dtptrs_obj->info = LAPACKE_dtptrs( dtptrs_obj->matrix_layout, dtptrs_obj->uplo,
                                  dtptrs_obj->trans, dtptrs_obj->diag, 
                                  dtptrs_obj->n, dtptrs_obj->nrhs,
                                  (const double *)dtptrs_obj->a, 
                                  dtptrs_obj->b,
                                  dtptrs_obj->ldb );

    if( dtptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtptrs is wrong\n", dtptrs_obj->info );
    }
    if( dtptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtptrs is wrong\n", 
        dtptrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtptrs_obj->diff =  computeDiff_d( dtptrs_obj->b_bufsize, 
                dtptrs_obj->b, dtptrs_obj->bref );

}

TEST_F(dtptrs_test, dtptrs1) {
    EXPECT_NEAR(0.0, dtptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtptrs_test, dtptrs2) {
    EXPECT_NEAR(0.0, dtptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtptrs_test, dtptrs3) {
    EXPECT_NEAR(0.0, dtptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtptrs_test, dtptrs4) {
    EXPECT_NEAR(0.0, dtptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin tptrs_float_parameters  class definition */
class tptrs_float_parameters{
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
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tptrs_float_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i,  
                                lapack_int nrhs_i, lapack_int ldb_i);
              
      ~tptrs_float_parameters (); 
};  /* end of tptrs_float_parameters  class definition */


/* Constructor tptrs_float_parameters definition */
tptrs_float_parameters:: tptrs_float_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i,  
                                 lapack_int nrhs_i, lapack_int ldb_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tptrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
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
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*(n+1)/2));
    lapacke_gtest_alloc_float_buffer_pair( &b, &bref, b_bufsize);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       tptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, n, (n+1)/2, uplo);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

tptrs_float_parameters:: ~tptrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tptrs_float_parameters object: destructor invoked. \n");
#endif
   tptrs_free();
}

//  Test fixture class definition
class stptrs_test  : public  ::testing::Test {
public:
   tptrs_float_parameters  *stptrs_obj;
   void SetUp();  
   void TearDown () { delete stptrs_obj; }
};


void stptrs_test::SetUp(){

    /* LAPACKE STPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_stptrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const float *a, 
                                    float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_stptrs STPTRS;

    stptrs_obj = new tptrs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    stptrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stptrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stptrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stptrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    STPTRS = (Fptr_NL_LAPACKE_stptrs)dlsym(stptrs_obj->hModule, "LAPACKE_stptrs");
    ASSERT_TRUE(STPTRS != NULL) << "failed to get the Netlib LAPACKE_stptrs symbol";
    

    stptrs_obj->inforef = STPTRS( stptrs_obj->matrix_layout, stptrs_obj->uplo,
                                  stptrs_obj->trans, stptrs_obj->diag, 
                                  stptrs_obj->n, stptrs_obj->nrhs,
                                  (const float *)stptrs_obj->aref, 
                                  stptrs_obj->bref,
                                  stptrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    stptrs_obj->info = LAPACKE_stptrs( stptrs_obj->matrix_layout, stptrs_obj->uplo,
                                  stptrs_obj->trans, stptrs_obj->diag, 
                                  stptrs_obj->n, stptrs_obj->nrhs,
                                  (const float *)stptrs_obj->a, 
                                  stptrs_obj->b,
                                  stptrs_obj->ldb );

    if( stptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stptrs is wrong\n", stptrs_obj->info );
    }
    if( stptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stptrs is wrong\n", 
        stptrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stptrs_obj->diff =  computeDiff_s( stptrs_obj->b_bufsize, 
                stptrs_obj->b, stptrs_obj->bref );

}

TEST_F(stptrs_test, stptrs1) {
    EXPECT_NEAR(0.0, stptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stptrs_test, stptrs2) {
    EXPECT_NEAR(0.0, stptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stptrs_test, stptrs3) {
    EXPECT_NEAR(0.0, stptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stptrs_test, stptrs4) {
    EXPECT_NEAR(0.0, stptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin tptrs_scomplex_parameters  class definition */
class tptrs_scomplex_parameters{
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
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tptrs_scomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i,  
                                lapack_int nrhs_i, lapack_int ldb_i);
              
      ~tptrs_scomplex_parameters (); 
};  /* end of tptrs_scomplex_parameters  class definition */


/* Constructor tptrs_scomplex_parameters definition */
tptrs_scomplex_parameters:: tptrs_scomplex_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i,  
                                 lapack_int nrhs_i, lapack_int ldb_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tptrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
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
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*(n+1)/2));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &b, &bref, b_bufsize);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       tptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, n, (n+1)/2, uplo);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

tptrs_scomplex_parameters:: ~tptrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tptrs_scomplex_parameters object: destructor invoked. \n");
#endif
   tptrs_free();
}

//  Test fixture class definition
class ctptrs_test  : public  ::testing::Test {
public:
   tptrs_scomplex_parameters  *ctptrs_obj;
   void SetUp();  
   void TearDown () { delete ctptrs_obj; }
};


void ctptrs_test::SetUp(){

    /* LAPACKE CTPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ctptrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const lapack_complex_float *a, 
                                    lapack_complex_float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ctptrs CTPTRS;

    ctptrs_obj = new tptrs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    ctptrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctptrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctptrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctptrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CTPTRS = (Fptr_NL_LAPACKE_ctptrs)dlsym(ctptrs_obj->hModule, "LAPACKE_ctptrs");
    ASSERT_TRUE(CTPTRS != NULL) << "failed to get the Netlib LAPACKE_ctptrs symbol";
    

    ctptrs_obj->inforef = CTPTRS( ctptrs_obj->matrix_layout, ctptrs_obj->uplo,
                                  ctptrs_obj->trans, ctptrs_obj->diag, 
                                  ctptrs_obj->n, ctptrs_obj->nrhs,
                                  (const lapack_complex_float *)ctptrs_obj->aref, 
                                  ctptrs_obj->bref,
                                  ctptrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    ctptrs_obj->info = LAPACKE_ctptrs( ctptrs_obj->matrix_layout, ctptrs_obj->uplo,
                                  ctptrs_obj->trans, ctptrs_obj->diag, 
                                  ctptrs_obj->n, ctptrs_obj->nrhs,
                                  (const lapack_complex_float *)ctptrs_obj->a, 
                                  ctptrs_obj->b,
                                  ctptrs_obj->ldb );

    if( ctptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctptrs is wrong\n", ctptrs_obj->info );
    }
    if( ctptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctptrs is wrong\n", 
        ctptrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctptrs_obj->diff =  computeDiff_c( ctptrs_obj->b_bufsize, 
                ctptrs_obj->b, ctptrs_obj->bref );

}

TEST_F(ctptrs_test, ctptrs1) {
    EXPECT_NEAR(0.0, ctptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctptrs_test, ctptrs2) {
    EXPECT_NEAR(0.0, ctptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctptrs_test, ctptrs3) {
    EXPECT_NEAR(0.0, ctptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctptrs_test, ctptrs4) {
    EXPECT_NEAR(0.0, ctptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin tptrs_dcomplex_parameters  class definition */
class tptrs_dcomplex_parameters{
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
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tptrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i,  
                                lapack_int nrhs_i, lapack_int ldb_i);
              
      ~tptrs_dcomplex_parameters (); 
};  /* end of tptrs_dcomplex_parameters  class definition */


/* Constructor tptrs_dcomplex_parameters definition */
tptrs_dcomplex_parameters:: tptrs_dcomplex_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i,  
                                 lapack_int nrhs_i, lapack_int ldb_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tptrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
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
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*(n+1)/2));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL)   ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       tptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, n, (n+1)/2, uplo);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

tptrs_dcomplex_parameters:: ~tptrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tptrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   tptrs_free();
}

//  Test fixture class definition
class ztptrs_test  : public  ::testing::Test {
public:
   tptrs_dcomplex_parameters  *ztptrs_obj;
   void SetUp();  
   void TearDown () { delete ztptrs_obj; }
};


void ztptrs_test::SetUp(){

    /* LAPACKE ZTPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ztptrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int nrhs, const lapack_complex_double *a, 
                                    lapack_complex_double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ztptrs ZTPTRS;

    ztptrs_obj = new tptrs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    ztptrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztptrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztptrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztptrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZTPTRS = (Fptr_NL_LAPACKE_ztptrs)dlsym(ztptrs_obj->hModule, "LAPACKE_ztptrs");
    ASSERT_TRUE(ZTPTRS != NULL) << "failed to get the Netlib LAPACKE_ztptrs symbol";
    

    ztptrs_obj->inforef = ZTPTRS( ztptrs_obj->matrix_layout, ztptrs_obj->uplo,
                                  ztptrs_obj->trans, ztptrs_obj->diag, 
                                  ztptrs_obj->n, ztptrs_obj->nrhs,
                                  (const lapack_complex_double *)ztptrs_obj->aref, 
                                  ztptrs_obj->bref,
                                  ztptrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    ztptrs_obj->info = LAPACKE_ztptrs( ztptrs_obj->matrix_layout, ztptrs_obj->uplo,
                                  ztptrs_obj->trans, ztptrs_obj->diag, 
                                  ztptrs_obj->n, ztptrs_obj->nrhs,
                                  (const lapack_complex_double *)ztptrs_obj->a, 
                                  ztptrs_obj->b,
                                  ztptrs_obj->ldb );

    if( ztptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztptrs is wrong\n", ztptrs_obj->info );
    }
    if( ztptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztptrs is wrong\n", 
        ztptrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztptrs_obj->diff =  computeDiff_z( ztptrs_obj->b_bufsize, 
                ztptrs_obj->b, ztptrs_obj->bref );

}

TEST_F(ztptrs_test, ztptrs1) {
    EXPECT_NEAR(0.0, ztptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztptrs_test, ztptrs2) {
    EXPECT_NEAR(0.0, ztptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztptrs_test, ztptrs3) {
    EXPECT_NEAR(0.0, ztptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztptrs_test, ztptrs4) {
    EXPECT_NEAR(0.0, ztptrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
