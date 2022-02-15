#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define hesv_aa_free() \
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

/* Begin hesv_aa_scomplex_parameters  class definition */
class hesv_aa_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      int b_bufsize;
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

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      hesv_aa_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hesv_aa_scomplex_parameters ();
};  /* end of hesv_aa_scomplex_parameters  class definition */


/* Constructor hesv_aa_scomplex_parameters definition */
hesv_aa_scomplex_parameters:: hesv_aa_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hesv_aa scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hesv_aa_scomplex_parameters object: malloc error.";
       hesv_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hesv_aa_scomplex_parameters:: ~hesv_aa_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hesv_aa_scomplex_parameters object: destructor invoked. \n");
#endif
   hesv_aa_free();
}


//  Test fixture class definition
class chesv_aa_test  : public  ::testing::Test {
public:
   hesv_aa_scomplex_parameters  *chesv_aa_obj;
   void SetUp();  
   void TearDown () { delete chesv_aa_obj; }
};


void chesv_aa_test::SetUp(){

    /* LAPACKE CHESV_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_chesv_aa) ( int matrix_layout, char uplo,
                                           lapack_int n, lapack_int nrhs,
                                           lapack_complex_float * a,
                                 lapack_int lda,  lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_chesv_aa CHESV_AA;

    chesv_aa_obj = new  hesv_aa_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    chesv_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chesv_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chesv_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chesv_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHESV_AA = (Fptr_NL_LAPACKE_chesv_aa)dlsym(chesv_aa_obj->hModule, "LAPACKE_chesv_aa");
    ASSERT_TRUE(CHESV_AA != NULL) << "failed to get the Netlib LAPACKE_chesv_aa symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    chesv_aa_obj->inforef = CHESV_AA( chesv_aa_obj->matrix_layout,
                                  chesv_aa_obj->uplo, chesv_aa_obj->n,
                                  chesv_aa_obj->nrhs,
                                  ( lapack_complex_float *)chesv_aa_obj->aref,
                                  chesv_aa_obj->lda, chesv_aa_obj->ipivref,
                                  chesv_aa_obj->bref, chesv_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    chesv_aa_obj->info = LAPACKE_chesv_aa( chesv_aa_obj->matrix_layout,
                chesv_aa_obj->uplo, chesv_aa_obj->n, chesv_aa_obj->nrhs,
                                  ( lapack_complex_float *)chesv_aa_obj->a,
                               chesv_aa_obj->lda, chesv_aa_obj->ipiv,
                                 chesv_aa_obj->b, chesv_aa_obj->ldb );


    if( chesv_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chesv_aa is wrong\n",
                    chesv_aa_obj->info );
    }
    if( chesv_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chesv_aa is wrong\n",
        chesv_aa_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    chesv_aa_obj->diff =  computeDiff_c( chesv_aa_obj->b_bufsize,
                           chesv_aa_obj->b, chesv_aa_obj->bref );
}

TEST_F(chesv_aa_test, chesv_aa1) {
    EXPECT_NEAR(0.0, chesv_aa_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chesv_aa_test, chesv_aa2) {
    EXPECT_NEAR(0.0, chesv_aa_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chesv_aa_test, chesv_aa3) {
    EXPECT_NEAR(0.0, chesv_aa_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chesv_aa_test, chesv_aa4) {
    EXPECT_NEAR(0.0, chesv_aa_obj->diff, LAPACKE_GTEST_THRESHOLD);
}



/* Begin hesv_aa_dcomplex_parameters  class definition */
class hesv_aa_dcomplex_parameters{
   public:
      int b_bufsize;
      double diff; //to capture reference and libflame o/ps' difference
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

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      hesv_aa_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hesv_aa_dcomplex_parameters ();
};  /* end of hesv_aa_dcomplex_parameters  class definition */


/* Constructor hesv_aa_dcomplex_parameters definition */
hesv_aa_dcomplex_parameters:: hesv_aa_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hesv_aa DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hesv_aa_dcomplex_parameters object: malloc error.";
       hesv_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hesv_aa_dcomplex_parameters:: ~hesv_aa_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hesv_aa_dcomplex_parameters object: destructor invoked. \n");
#endif
   hesv_aa_free();
}


//  Test fixture class definition
class zhesv_aa_test  : public  ::testing::Test {
public:
   hesv_aa_dcomplex_parameters  *zhesv_aa_obj;
   void SetUp();  
   void TearDown () { delete zhesv_aa_obj; }
};


void zhesv_aa_test::SetUp(){

    /* LAPACKE ZHESV_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_zhesv_aa) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_double * a,
                                  lapack_int lda,  lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zhesv_aa ZHESV_AA;

    zhesv_aa_obj = new  hesv_aa_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zhesv_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhesv_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhesv_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhesv_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHESV_AA = (Fptr_NL_LAPACKE_zhesv_aa)dlsym(zhesv_aa_obj->hModule, "LAPACKE_zhesv_aa");
    ASSERT_TRUE(ZHESV_AA != NULL) << "failed to get the Netlib LAPACKE_zhesv_aa symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zhesv_aa_obj->inforef = ZHESV_AA( zhesv_aa_obj->matrix_layout,
                                  zhesv_aa_obj->uplo, zhesv_aa_obj->n,
                                  zhesv_aa_obj->nrhs,
                                  ( lapack_complex_double *)zhesv_aa_obj->aref,
                                  zhesv_aa_obj->lda, zhesv_aa_obj->ipivref,
                                  zhesv_aa_obj->bref, zhesv_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zhesv_aa_obj->info = LAPACKE_zhesv_aa( zhesv_aa_obj->matrix_layout,
                zhesv_aa_obj->uplo, zhesv_aa_obj->n, zhesv_aa_obj->nrhs,
                                  ( lapack_complex_double *)zhesv_aa_obj->a,
                               zhesv_aa_obj->lda, zhesv_aa_obj->ipiv,
                                 zhesv_aa_obj->b, zhesv_aa_obj->ldb );

    if( zhesv_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhesv_aa is wrong\n",
                    zhesv_aa_obj->info );
    }
    if( zhesv_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhesv_aa is wrong\n",
        zhesv_aa_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zhesv_aa_obj->diff =  computeDiff_z( zhesv_aa_obj->b_bufsize,
                           zhesv_aa_obj->b, zhesv_aa_obj->bref );
}

TEST_F(zhesv_aa_test, zhesv_aa1) {
    EXPECT_NEAR(0.0, zhesv_aa_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhesv_aa_test, zhesv_aa2) {
    EXPECT_NEAR(0.0, zhesv_aa_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhesv_aa_test, zhesv_aa3) {
    EXPECT_NEAR(0.0, zhesv_aa_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhesv_aa_test, zhesv_aa4) {
    EXPECT_NEAR(0.0, zhesv_aa_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
