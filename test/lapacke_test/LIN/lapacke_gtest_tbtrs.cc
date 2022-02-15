#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define tbtrs_free() \
  if (b != NULL)    free (b   ); \
  if (bref != NULL) free (bref); \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if( hModule != NULL) dlclose(hModule); \
  if(dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin tbtrs_double_parameters  class definition */
class tbtrs_double_parameters{
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
      lapack_int kd; // The number of superdiagonals or subdiagonals in A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tbtrs_double_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int lda_i, 
                            lapack_int nrhs_i, lapack_int ldb_i, lapack_int kd_i);
              
      ~tbtrs_double_parameters (); 
};  /* end of tbtrs_double_parameters  class definition */


/* Constructor tbtrs_double_parameters definition */
tbtrs_double_parameters:: tbtrs_double_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i, lapack_int lda_i, 
                        lapack_int nrhs_i, lapack_int ldb_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
	kd = kd_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tbtrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
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
       tbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
//    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_double_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

tbtrs_double_parameters:: ~tbtrs_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbtrs_double_parameters object: destructor invoked. \n");
#endif
   tbtrs_free();
}

//  Test fixture class definition
class dtbtrs_test  : public  ::testing::Test {
public:
   tbtrs_double_parameters  *dtbtrs_obj;
   void SetUp();  
   void TearDown () { delete dtbtrs_obj; }
};


void dtbtrs_test::SetUp(){

    /* LAPACKE DTRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_dtbtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                               lapack_int kd, lapack_int nrhs, const double *a, 
                                    lapack_int lda, double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_dtbtrs DTRTRS;

    dtbtrs_obj = new tbtrs_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb,
                           lin_solver_paramslist[idx].kd );

    idx = Circular_Increment_Index(idx);

    dtbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DTRTRS = (Fptr_NL_LAPACKE_dtbtrs)dlsym(dtbtrs_obj->hModule, "LAPACKE_dtbtrs");
    ASSERT_TRUE(DTRTRS != NULL) << "failed to get the Netlib LAPACKE_dtbtrs symbol";
    

    dtbtrs_obj->inforef = DTRTRS( dtbtrs_obj->matrix_layout, dtbtrs_obj->uplo,
                                  dtbtrs_obj->trans, dtbtrs_obj->diag, 
                                dtbtrs_obj->n, dtbtrs_obj->kd, dtbtrs_obj->nrhs,
                                  (const double *)dtbtrs_obj->aref, 
                                  dtbtrs_obj->lda, dtbtrs_obj->bref,
                                  dtbtrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    dtbtrs_obj->info = LAPACKE_dtbtrs( dtbtrs_obj->matrix_layout, dtbtrs_obj->uplo,
                                  dtbtrs_obj->trans, dtbtrs_obj->diag, 
                                  dtbtrs_obj->n, dtbtrs_obj->kd, dtbtrs_obj->nrhs,
                                  (const double *)dtbtrs_obj->a, 
                                  dtbtrs_obj->lda, dtbtrs_obj->b,
                                  dtbtrs_obj->ldb );

    if( dtbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtbtrs is wrong\n", dtbtrs_obj->info );
    }
    if( dtbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtbtrs is wrong\n", 
        dtbtrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtbtrs_obj->diff =  computeDiff_d( dtbtrs_obj->b_bufsize, 
                dtbtrs_obj->b, dtbtrs_obj->bref );

}

TEST_F(dtbtrs_test, dtbtrs1) {
    EXPECT_NEAR(0.0, dtbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtbtrs_test, dtbtrs2) {
    EXPECT_NEAR(0.0, dtbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtbtrs_test, dtbtrs3) {
    EXPECT_NEAR(0.0, dtbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtbtrs_test, dtbtrs4) {
    EXPECT_NEAR(0.0, dtbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}



/* Begin tbtrs_float_parameters  class definition */
class tbtrs_float_parameters{
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
      lapack_int kd; // The number of superdiagonals or subdiagonals in A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tbtrs_float_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i, lapack_int kd_i);
              
      ~tbtrs_float_parameters (); 
};  /* end of tbtrs_float_parameters  class definition */


/* Constructor tbtrs_float_parameters definition */
tbtrs_float_parameters:: tbtrs_float_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i, lapack_int lda_i, 
                                 lapack_int nrhs_i, lapack_int ldb_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
	kd = kd_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tbtrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
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
       tbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_float_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

tbtrs_float_parameters:: ~tbtrs_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbtrs_float_parameters object: destructor invoked. \n");
#endif
   tbtrs_free();
}

//  Test fixture class definition
class stbtrs_test  : public  ::testing::Test {
public:
   tbtrs_float_parameters  *stbtrs_obj;
   void SetUp();  
   void TearDown () { delete stbtrs_obj; }
};


void stbtrs_test::SetUp(){

    /* LAPACKE STRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_stbtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int kd, lapack_int nrhs, const float *a, 
                                    lapack_int lda, float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_stbtrs STRTRS;

    stbtrs_obj = new tbtrs_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb,
                           lin_solver_paramslist[idx].kd );

    idx = Circular_Increment_Index(idx);

    stbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    STRTRS = (Fptr_NL_LAPACKE_stbtrs)dlsym(stbtrs_obj->hModule, "LAPACKE_stbtrs");
    ASSERT_TRUE(STRTRS != NULL) << "failed to get the Netlib LAPACKE_stbtrs symbol";
    

    stbtrs_obj->inforef = STRTRS( stbtrs_obj->matrix_layout, stbtrs_obj->uplo,
                                  stbtrs_obj->trans, stbtrs_obj->diag, 
                                  stbtrs_obj->n, stbtrs_obj->kd, stbtrs_obj->nrhs,
                                  (const float *)stbtrs_obj->aref, 
                                  stbtrs_obj->lda, stbtrs_obj->bref,
                                  stbtrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    stbtrs_obj->info = LAPACKE_stbtrs( stbtrs_obj->matrix_layout, stbtrs_obj->uplo,
                                  stbtrs_obj->trans, stbtrs_obj->diag, 
                                  stbtrs_obj->n, stbtrs_obj->kd, stbtrs_obj->nrhs,
                                  (const float *)stbtrs_obj->a, 
                                  stbtrs_obj->lda, stbtrs_obj->b,
                                  stbtrs_obj->ldb );

    if( stbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stbtrs is wrong\n", stbtrs_obj->info );
    }
    if( stbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stbtrs is wrong\n", 
        stbtrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stbtrs_obj->diff =  computeDiff_s( stbtrs_obj->b_bufsize, 
                stbtrs_obj->b, stbtrs_obj->bref );

}

TEST_F(stbtrs_test, stbtrs1) {
    EXPECT_NEAR(0.0, stbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stbtrs_test, stbtrs2) {
    EXPECT_NEAR(0.0, stbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stbtrs_test, stbtrs3) {
    EXPECT_NEAR(0.0, stbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stbtrs_test, stbtrs4) {
    EXPECT_NEAR(0.0, stbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin tbtrs_scomplex_parameters  class definition */
class tbtrs_scomplex_parameters{
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
      lapack_int kd; // The number of superdiagonals or subdiagonals in A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tbtrs_scomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i, lapack_int kd_i);
              
      ~tbtrs_scomplex_parameters (); 
};  /* end of tbtrs_scomplex_parameters  class definition */


/* Constructor tbtrs_scomplex_parameters definition */
tbtrs_scomplex_parameters:: tbtrs_scomplex_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i, lapack_int lda_i, 
                                 lapack_int nrhs_i, lapack_int ldb_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
	kd = kd_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tbtrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
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
       tbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

tbtrs_scomplex_parameters:: ~tbtrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbtrs_scomplex_parameters object: destructor invoked. \n");
#endif
   tbtrs_free();
}

//  Test fixture class definition
class ctbtrs_test  : public  ::testing::Test {
public:
   tbtrs_scomplex_parameters  *ctbtrs_obj;
   void SetUp();  
   void TearDown () { delete ctbtrs_obj; }
};


void ctbtrs_test::SetUp(){

    /* LAPACKE CTRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ctbtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int kd, lapack_int nrhs,
									const lapack_complex_float *a, 
                                    lapack_int lda, lapack_complex_float *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ctbtrs CTRTRS;

    ctbtrs_obj = new tbtrs_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb,
                           lin_solver_paramslist[idx].kd );

    idx = Circular_Increment_Index(idx);

    ctbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CTRTRS = (Fptr_NL_LAPACKE_ctbtrs)dlsym(ctbtrs_obj->hModule, "LAPACKE_ctbtrs");
    ASSERT_TRUE(CTRTRS != NULL) << "failed to get the Netlib LAPACKE_ctbtrs symbol";
    

    ctbtrs_obj->inforef = CTRTRS( ctbtrs_obj->matrix_layout, ctbtrs_obj->uplo,
                                  ctbtrs_obj->trans, ctbtrs_obj->diag, 
                                  ctbtrs_obj->n, ctbtrs_obj->kd, ctbtrs_obj->nrhs,
                                  (const lapack_complex_float *)ctbtrs_obj->aref, 
                                  ctbtrs_obj->lda, ctbtrs_obj->bref,
                                  ctbtrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    ctbtrs_obj->info = LAPACKE_ctbtrs( ctbtrs_obj->matrix_layout, ctbtrs_obj->uplo,
                                  ctbtrs_obj->trans, ctbtrs_obj->diag, 
                                  ctbtrs_obj->n, ctbtrs_obj->kd, ctbtrs_obj->nrhs,
                                  (const lapack_complex_float *)ctbtrs_obj->a, 
                                  ctbtrs_obj->lda, ctbtrs_obj->b,
                                  ctbtrs_obj->ldb );

    if( ctbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctbtrs is wrong\n", ctbtrs_obj->info );
    }
    if( ctbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctbtrs is wrong\n", 
        ctbtrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctbtrs_obj->diff =  computeDiff_c( ctbtrs_obj->b_bufsize, 
                ctbtrs_obj->b, ctbtrs_obj->bref );

}

TEST_F(ctbtrs_test, ctbtrs1) {
    EXPECT_NEAR(0.0, ctbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctbtrs_test, ctbtrs2) {
    EXPECT_NEAR(0.0, ctbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctbtrs_test, ctbtrs3) {
    EXPECT_NEAR(0.0, ctbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctbtrs_test, ctbtrs4) {
    EXPECT_NEAR(0.0, ctbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin tbtrs_dcomplex_parameters  class definition */
class tbtrs_dcomplex_parameters{
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
      lapack_int kd; // The number of superdiagonals or subdiagonals in A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tbtrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i, char trans_i,
                                char diag_i, lapack_int n_i, lapack_int lda_i, 
                                lapack_int nrhs_i, lapack_int ldb_i, lapack_int kd_i);
              
      ~tbtrs_dcomplex_parameters (); 
};  /* end of tbtrs_dcomplex_parameters  class definition */


/* Constructor tbtrs_dcomplex_parameters definition */
tbtrs_dcomplex_parameters:: tbtrs_dcomplex_parameters ( int matrix_layout_i, 
                                 char uplo_i, char trans_i, char diag_i, 
                                       lapack_int n_i, lapack_int lda_i, 
                                 lapack_int nrhs_i, lapack_int ldb_i, lapack_int kd_i) {
                                    
    matrix_layout = matrix_layout_i;
    n = n_i;
	kd = kd_i;
    uplo = uplo_i;
    trans = trans_i;
    diag = diag_i;
    lda = lda_i;
    ldb = ldb_i;
    nrhs = nrhs_i;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tbtrs Double:  n: %d, Uplo: %c trans: %c diag: %c lda: %d  \
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
       tbtrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
    
   } /* end of Constructor  */

tbtrs_dcomplex_parameters:: ~tbtrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tbtrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   tbtrs_free();
}

//  Test fixture class definition
class ztbtrs_test  : public  ::testing::Test {
public:
   tbtrs_dcomplex_parameters  *ztbtrs_obj;
   void SetUp();  
   void TearDown () { delete ztbtrs_obj; }
};


void ztbtrs_test::SetUp(){

    /* LAPACKE ZTRTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_ztbtrs) (int matrix_layout, char uplo,
                                    char trans, char diag, lapack_int n,
                                    lapack_int kd, lapack_int nrhs,
									const lapack_complex_double *a, 
                                    lapack_int lda, lapack_complex_double *b, lapack_int ldb);

    Fptr_NL_LAPACKE_ztbtrs ZTRTRS;

    ztbtrs_obj = new tbtrs_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
                           lin_solver_paramslist[idx].nrhs,
                           lin_solver_paramslist[idx].ldb,
                           lin_solver_paramslist[idx].kd );

    idx = Circular_Increment_Index(idx);

    ztbtrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztbtrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztbtrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztbtrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZTRTRS = (Fptr_NL_LAPACKE_ztbtrs)dlsym(ztbtrs_obj->hModule, "LAPACKE_ztbtrs");
    ASSERT_TRUE(ZTRTRS != NULL) << "failed to get the Netlib LAPACKE_ztbtrs symbol";
    

    ztbtrs_obj->inforef = ZTRTRS( ztbtrs_obj->matrix_layout, ztbtrs_obj->uplo,
                                  ztbtrs_obj->trans, ztbtrs_obj->diag, 
                                  ztbtrs_obj->n, ztbtrs_obj->kd, ztbtrs_obj->nrhs,
                                  (const lapack_complex_double *)ztbtrs_obj->aref, 
                                  ztbtrs_obj->lda, ztbtrs_obj->bref,
                                  ztbtrs_obj->ldb );

    /* Compute libflame's Lapacke o/p  */
    ztbtrs_obj->info = LAPACKE_ztbtrs( ztbtrs_obj->matrix_layout, ztbtrs_obj->uplo,
                                  ztbtrs_obj->trans, ztbtrs_obj->diag, 
                                  ztbtrs_obj->n, ztbtrs_obj->kd, ztbtrs_obj->nrhs,
                                  (const lapack_complex_double *)ztbtrs_obj->a, 
                                  ztbtrs_obj->lda, ztbtrs_obj->b,
                                  ztbtrs_obj->ldb );

    if( ztbtrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztbtrs is wrong\n", ztbtrs_obj->info );
    }
    if( ztbtrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztbtrs is wrong\n", 
        ztbtrs_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztbtrs_obj->diff =  computeDiff_z( ztbtrs_obj->b_bufsize, 
                ztbtrs_obj->b, ztbtrs_obj->bref );

}

TEST_F(ztbtrs_test, ztbtrs1) {
    EXPECT_NEAR(0.0, ztbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztbtrs_test, ztbtrs2) {
    EXPECT_NEAR(0.0, ztbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztbtrs_test, ztbtrs3) {
    EXPECT_NEAR(0.0, ztbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztbtrs_test, ztbtrs4) {
    EXPECT_NEAR(0.0, ztbtrs_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

