#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define hesv_free() \
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

/* Begin hesv_scomplex_parameters  class definition */
class hesv_scomplex_parameters{
   public:
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
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      hesv_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hesv_scomplex_parameters ();
};  /* end of hesv_scomplex_parameters  class definition */


/* Constructor hesv_scomplex_parameters definition */
hesv_scomplex_parameters:: hesv_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hesv scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hesv_scomplex_parameters object: malloc error.";
       hesv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hesv_scomplex_parameters:: ~hesv_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hesv_scomplex_parameters object: destructor invoked. \n");
#endif
   hesv_free();
}


//  Test fixture class definition
class chesv_test  : public  ::testing::Test {
public:
   hesv_scomplex_parameters  *chesv_obj;
   void SetUp();  
   void TearDown () { delete chesv_obj; }
};


void chesv_test::SetUp(){

    /* LAPACKE CHESV prototype */
    typedef int (*Fptr_NL_LAPACKE_chesv) ( int matrix_layout, char uplo,
                                           lapack_int n, lapack_int nrhs,
                                           lapack_complex_float * a,
                                 lapack_int lda,  lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_chesv CHESV;

     /* LAPACKE CHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf) ( int matrix_layout, char uplo,
    lapack_int n, lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_chetrf CHETRF;


    chesv_obj = new  hesv_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    chesv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chesv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chesv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chesv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHESV = (Fptr_NL_LAPACKE_chesv)dlsym(chesv_obj->hModule, "LAPACKE_chesv");
    ASSERT_TRUE(CHESV != NULL) << "failed to get the Netlib LAPACKE_chesv symbol";

    CHETRF = (Fptr_NL_LAPACKE_chetrf)dlsym(chesv_obj->hModule,"LAPACKE_chetrf");
    ASSERT_TRUE(CHETRF != NULL) << "failed to get the Netlib LAPACKE_chetrf symbol";

    /* Pre condition: need to call hetrf - before calling hesv function */

    /* Compute the Netlib-Lapacke's reference o/p */
    chesv_obj->inforef = CHETRF( chesv_obj->matrix_layout,
                            chesv_obj->uplo, chesv_obj->n,
                                           chesv_obj->aref,
                      chesv_obj->lda, chesv_obj->ipivref);

    chesv_obj->inforef = CHESV( chesv_obj->matrix_layout,
                                  chesv_obj->uplo, chesv_obj->n,
                                  chesv_obj->nrhs,
                                  ( lapack_complex_float *)chesv_obj->aref,
                                  chesv_obj->lda, chesv_obj->ipivref,
                                  chesv_obj->bref, chesv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    chesv_obj->info = LAPACKE_chetrf( chesv_obj->matrix_layout,
                                 chesv_obj->uplo, chesv_obj->n,
                                     chesv_obj->a,
                               chesv_obj->lda, chesv_obj->ipiv);

    chesv_obj->info = LAPACKE_chesv( chesv_obj->matrix_layout,
                chesv_obj->uplo, chesv_obj->n, chesv_obj->nrhs,
                                  ( lapack_complex_float *)chesv_obj->a,
                               chesv_obj->lda, chesv_obj->ipiv,
                                 chesv_obj->b, chesv_obj->ldb );


    if( chesv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chesv is wrong\n",
                    chesv_obj->info );
    }
    if( chesv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chesv is wrong\n",
        chesv_obj->inforef );
    }
}

TEST_F(chesv_test, chesv1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chesv_obj->b_bufsize,
                           chesv_obj->b, chesv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chesv_test, chesv2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chesv_obj->b_bufsize,
                           chesv_obj->b, chesv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chesv_test, chesv3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chesv_obj->b_bufsize,
                           chesv_obj->b, chesv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chesv_test, chesv4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chesv_obj->b_bufsize,
                           chesv_obj->b, chesv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin hesv_dcomplex_parameters  class definition */
class hesv_dcomplex_parameters{
   public:
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
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      hesv_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hesv_dcomplex_parameters ();
};  /* end of hesv_dcomplex_parameters  class definition */


/* Constructor hesv_dcomplex_parameters definition */
hesv_dcomplex_parameters:: hesv_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hesv DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hesv_dcomplex_parameters object: malloc error.";
       hesv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hesv_dcomplex_parameters:: ~hesv_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hesv_dcomplex_parameters object: destructor invoked. \n");
#endif
   hesv_free();
}


//  Test fixture class definition
class zhesv_test  : public  ::testing::Test {
public:
   hesv_dcomplex_parameters  *zhesv_obj;
   void SetUp();  
   void TearDown () { delete zhesv_obj; }
};


void zhesv_test::SetUp(){

    /* LAPACKE ZHESV prototype */
    typedef int (*Fptr_NL_LAPACKE_zhesv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_double * a,
                                  lapack_int lda,  lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zhesv ZHESV;

     /* LAPACKE ZHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zhetrf ZHETRF;


    zhesv_obj = new  hesv_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zhesv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhesv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhesv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhesv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHESV = (Fptr_NL_LAPACKE_zhesv)dlsym(zhesv_obj->hModule, "LAPACKE_zhesv");
    ASSERT_TRUE(ZHESV != NULL) << "failed to get the Netlib LAPACKE_zhesv symbol";

    ZHETRF = (Fptr_NL_LAPACKE_zhetrf)dlsym(zhesv_obj->hModule,"LAPACKE_zhetrf");
    ASSERT_TRUE(ZHETRF != NULL) << "failed to get the Netlib LAPACKE_zhetrf symbol";

    /* Pre condition: need to call hetrf - before calling hesv function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zhesv_obj->inforef = ZHETRF( zhesv_obj->matrix_layout,
                                    zhesv_obj->uplo, zhesv_obj->n,
                                     zhesv_obj->aref,
                               zhesv_obj->lda, zhesv_obj->ipivref);

    zhesv_obj->inforef = ZHESV( zhesv_obj->matrix_layout,
                                  zhesv_obj->uplo, zhesv_obj->n,
                                  zhesv_obj->nrhs,
                                  ( lapack_complex_double *)zhesv_obj->aref,
                                  zhesv_obj->lda, zhesv_obj->ipivref,
                                  zhesv_obj->bref, zhesv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zhesv_obj->info = LAPACKE_zhetrf( zhesv_obj->matrix_layout,
                                 zhesv_obj->uplo, zhesv_obj->n,
                                     zhesv_obj->a,
                               zhesv_obj->lda, zhesv_obj->ipiv);

    zhesv_obj->info = LAPACKE_zhesv( zhesv_obj->matrix_layout,
                zhesv_obj->uplo, zhesv_obj->n, zhesv_obj->nrhs,
                                  ( lapack_complex_double *)zhesv_obj->a,
                               zhesv_obj->lda, zhesv_obj->ipiv,
                                 zhesv_obj->b, zhesv_obj->ldb );


    if( zhesv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhesv is wrong\n",
                    zhesv_obj->info );
    }
    if( zhesv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhesv is wrong\n",
        zhesv_obj->inforef );
    }
}

TEST_F(zhesv_test, zhesv1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhesv_obj->b_bufsize,
                           zhesv_obj->b, zhesv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhesv_test, zhesv2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhesv_obj->b_bufsize,
                           zhesv_obj->b, zhesv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhesv_test, zhesv3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhesv_obj->b_bufsize,
                           zhesv_obj->b, zhesv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhesv_test, zhesv4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhesv_obj->b_bufsize,
                           zhesv_obj->b, zhesv_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
