#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define hesv_rk_free() \
       if (ipiv != NULL) free (ipiv); \
       if (bref != NULL) free (bref); \
       if (b != NULL)    free (b   ); \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (e != NULL)    free (e   ); \
       if (eref != NULL) free (eref); \
       if (ipivref != NULL)free (ipivref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin hesv_rk_scomplex_parameters  class definition */
class hesv_rk_scomplex_parameters{
   public:
      int b_bufsize;
      float diff; //to capture reference and libflame o/ps' difference
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A
      lapack_complex_float *e, *eref; //the superdiagonal (or subdiagonal) elements
      // of the Hermitian block diagonal matrix D with 1-by-1 or 2-by-2 diagonal blocks.
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      hesv_rk_scomplex_parameters( int matrix_layout_i, char uplo_i,
                                   lapack_int n_i, lapack_int lda_i,
                               lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hesv_rk_scomplex_parameters ();
};  /* end of hesv_rk_scomplex_parameters  class definition */


/* Constructor hesv_rk_scomplex_parameters definition */
hesv_rk_scomplex_parameters:: hesv_rk_scomplex_parameters( int matrix_layout_i,
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
   printf(" \n hesv_rk scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "hesv_rk_scomplex_parameters object: malloc error.";
       hesv_rk_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_scomplex_buffer_pair_rand( e, eref, n );

   } /* end of Constructor  */

hesv_rk_scomplex_parameters:: ~hesv_rk_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hesv_rk_scomplex_parameters object: destructor invoked. \n");
#endif
   hesv_rk_free();
}


//  Test fixture class definition
class chesv_rk_test  : public  ::testing::Test {
public:
   hesv_rk_scomplex_parameters  *chesv_rk_obj;
   void SetUp();  
   void TearDown () { delete chesv_rk_obj; }
};


void chesv_rk_test::SetUp(){

    /* LAPACKE CHESV_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_chesv_rk) ( int matrix_layout, char uplo,
                                           lapack_int n, lapack_int nrhs,
                                            lapack_complex_float* a,
                           lapack_int lda,  lapack_complex_float *e,
                                                  lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_chesv_rk CHESV_RK;
    chesv_rk_obj = new  hesv_rk_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    chesv_rk_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chesv_rk_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chesv_rk_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chesv_rk_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHESV_RK = (Fptr_NL_LAPACKE_chesv_rk)dlsym(chesv_rk_obj->hModule, "LAPACKE_chesv_rk");
    ASSERT_TRUE(CHESV_RK != NULL) << "failed to get the Netlib LAPACKE_chesv_rk symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    chesv_rk_obj->inforef = CHESV_RK( chesv_rk_obj->matrix_layout,
                                      chesv_rk_obj->uplo,
                                      chesv_rk_obj->n,
                                      chesv_rk_obj->nrhs,
                        ( lapack_complex_float *)chesv_rk_obj->aref,
                                      chesv_rk_obj->lda, 
                        ( lapack_complex_float *)chesv_rk_obj->eref,
                                      chesv_rk_obj->ipivref,
                                      chesv_rk_obj->bref,
                                      chesv_rk_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    chesv_rk_obj->info = LAPACKE_chesv_rk( chesv_rk_obj->matrix_layout,
                                           chesv_rk_obj->uplo, 
                                           chesv_rk_obj->n, 
                                           chesv_rk_obj->nrhs,
                              ( lapack_complex_float *)chesv_rk_obj->a,
                                           chesv_rk_obj->lda,
                              ( lapack_complex_float *)chesv_rk_obj->e,
                                           chesv_rk_obj->ipiv, 
                                           chesv_rk_obj->b,
                                           chesv_rk_obj->ldb );


    if( chesv_rk_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chesv_rk is wrong\n",
                    chesv_rk_obj->info );
    }
    if( chesv_rk_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chesv_rk is wrong\n",
        chesv_rk_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    chesv_rk_obj->diff =  computeDiff_c( chesv_rk_obj->b_bufsize,
                           chesv_rk_obj->b, chesv_rk_obj->bref );
}

TEST_F(chesv_rk_test, chesv_rk1) {
    EXPECT_NEAR(0.0, chesv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chesv_rk_test, chesv_rk2) {
    EXPECT_NEAR(0.0, chesv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chesv_rk_test, chesv_rk3) {
    EXPECT_NEAR(0.0, chesv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chesv_rk_test, chesv_rk4) {
    EXPECT_NEAR(0.0, chesv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin hesv_rk_dcomplex_parameters  class definition */
class hesv_rk_dcomplex_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      int b_bufsize;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A
      lapack_complex_double *e, *eref; //the superdiagonal (or subdiagonal) elements
      // of the Hermitian block diagonal matrix D with 1-by-1 or 2-by-2 diagonal blocks.
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      hesv_rk_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                    lapack_int n_i, lapack_int lda_i,
                                    lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hesv_rk_dcomplex_parameters ();
};  /* end of hesv_rk_dcomplex_parameters  class definition */


/* Constructor hesv_rk_dcomplex_parameters definition */
hesv_rk_dcomplex_parameters:: hesv_rk_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hesv_rk DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &b, &bref, b_bufsize);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "hesv_rk_dcomplex_parameters object: malloc error.";
       hesv_rk_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( e, eref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hesv_rk_dcomplex_parameters:: ~hesv_rk_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hesv_rk_dcomplex_parameters object: destructor invoked. \n");
#endif
   hesv_rk_free();
}


//  Test fixture class definition
class zhesv_rk_test  : public  ::testing::Test {
public:
   hesv_rk_dcomplex_parameters  *zhesv_rk_obj;
   void SetUp();  
   void TearDown () { delete zhesv_rk_obj; }
};


void zhesv_rk_test::SetUp(){

    /* LAPACKE ZHESV_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_zhesv_rk)(int matrix_layout, char uplo,
                                            lapack_int n, lapack_int nrhs,
                                             lapack_complex_double * a,
                                            lapack_int lda, 
                                             lapack_complex_double * e,
                                             lapack_int * ipiv,
                                            lapack_complex_double * b,
                                            lapack_int ldb );

    Fptr_NL_LAPACKE_zhesv_rk ZHESV_RK;

     /* LAPACKE ZHETRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf_rk) ( int matrix_layout,char uplo ,lapack_int n,
   lapack_complex_double* a,lapack_int lda, lapack_complex_double* e, lapack_int* ipiv );

    Fptr_NL_LAPACKE_zhetrf_rk ZHETRF_RK;


    zhesv_rk_obj = new  hesv_rk_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    zhesv_rk_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhesv_rk_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhesv_rk_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhesv_rk_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHESV_RK = (Fptr_NL_LAPACKE_zhesv_rk)dlsym(zhesv_rk_obj->hModule, "LAPACKE_zhesv_rk");
    ASSERT_TRUE(ZHESV_RK != NULL) << "failed to get the Netlib LAPACKE_zhesv_rk symbol";
    /* Compute the Netlib-Lapacke's reference o/p */
    zhesv_rk_obj->inforef = ZHESV_RK( zhesv_rk_obj->matrix_layout,
                                      zhesv_rk_obj->uplo, 
                                      zhesv_rk_obj->n,
                                      zhesv_rk_obj->nrhs,
                          ( lapack_complex_double *)zhesv_rk_obj->aref,
                                      zhesv_rk_obj->lda,
                          ( lapack_complex_double *)zhesv_rk_obj->eref,
                                      zhesv_rk_obj->ipivref,
                                      zhesv_rk_obj->bref,
                                      zhesv_rk_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zhesv_rk_obj->info = LAPACKE_zhesv_rk(zhesv_rk_obj->matrix_layout,
                                          zhesv_rk_obj->uplo,
                                          zhesv_rk_obj->n,
                                          zhesv_rk_obj->nrhs,
                                  ( lapack_complex_double *)zhesv_rk_obj->a,
                                          zhesv_rk_obj->lda,
                                  ( lapack_complex_double *)zhesv_rk_obj->e,
                                          zhesv_rk_obj->ipiv,
                                          zhesv_rk_obj->b,
                                          zhesv_rk_obj->ldb );


    if( zhesv_rk_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhesv_rk is wrong\n",
                    zhesv_rk_obj->info );
    }
    if( zhesv_rk_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhesv_rk is wrong\n",
        zhesv_rk_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zhesv_rk_obj->diff =  computeDiff_z( zhesv_rk_obj->b_bufsize,
                           zhesv_rk_obj->b,
                           zhesv_rk_obj->bref );
}

TEST_F(zhesv_rk_test, zhesv_rk1) {
    EXPECT_NEAR(0.0, zhesv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhesv_rk_test, zhesv_rk2) {
    EXPECT_NEAR(0.0, zhesv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhesv_rk_test, zhesv_rk3) {
    EXPECT_NEAR(0.0, zhesv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhesv_rk_test, zhesv_rk4) {
    EXPECT_NEAR(0.0, zhesv_rk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

