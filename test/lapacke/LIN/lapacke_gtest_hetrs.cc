#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define hetrs_free() \
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


/* Begin hetrs_scomplex_parameters  class definition */
class hetrs_scomplex_parameters{
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
      hetrs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hetrs_scomplex_parameters ();
};  /* end of hetrs_scomplex_parameters  class definition */


/* Constructor hetrs_scomplex_parameters definition */
hetrs_scomplex_parameters:: hetrs_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hetrs scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hetrs_scomplex_parameters object: malloc error.";
       hetrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hetrs_scomplex_parameters:: ~hetrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrs_scomplex_parameters object: destructor invoked. \n");
#endif
   hetrs_free();
}


//  Test fixture class definition
class chetrs_test  : public  ::testing::Test {
public:
   hetrs_scomplex_parameters  *chetrs_obj;
   void SetUp();  
   void TearDown () { delete chetrs_obj; }
};


void chetrs_test::SetUp(){

    /* LAPACKE CHETRS prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrs) ( int matrix_layout, char uplo,
                                           lapack_int n, lapack_int nrhs,
                                          const lapack_complex_float * a,
                                 lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_chetrs CHETRS;

     /* LAPACKE CHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf) ( int matrix_layout, char uplo,
    lapack_int n, lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_chetrf CHETRF;


    chetrs_obj = new  hetrs_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    chetrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHETRS = (Fptr_NL_LAPACKE_chetrs)dlsym(chetrs_obj->hModule, "LAPACKE_chetrs");
    ASSERT_TRUE(CHETRS != NULL) << "failed to get the Netlib LAPACKE_chetrs symbol";

    CHETRF = (Fptr_NL_LAPACKE_chetrf)dlsym(chetrs_obj->hModule,"LAPACKE_chetrf");
    ASSERT_TRUE(CHETRF != NULL) << "failed to get the Netlib LAPACKE_chetrf symbol";

    /* Pre condition: need to call hetrf - before calling hetrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    chetrs_obj->inforef = CHETRF( chetrs_obj->matrix_layout,
                            chetrs_obj->uplo, chetrs_obj->n,
                                           chetrs_obj->aref,
                      chetrs_obj->lda, chetrs_obj->ipivref);

    chetrs_obj->inforef = CHETRS( chetrs_obj->matrix_layout,
                                  chetrs_obj->uplo, chetrs_obj->n,
                                  chetrs_obj->nrhs,
                                  (const lapack_complex_float *)chetrs_obj->aref,
                                  chetrs_obj->lda, chetrs_obj->ipivref,
                                  chetrs_obj->bref, chetrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    chetrs_obj->info = LAPACKE_chetrf( chetrs_obj->matrix_layout,
                                 chetrs_obj->uplo, chetrs_obj->n,
                                     chetrs_obj->a,
                               chetrs_obj->lda, chetrs_obj->ipiv);

    chetrs_obj->info = LAPACKE_chetrs( chetrs_obj->matrix_layout,
                chetrs_obj->uplo, chetrs_obj->n, chetrs_obj->nrhs,
                                  (const lapack_complex_float *)chetrs_obj->a,
                               chetrs_obj->lda, chetrs_obj->ipiv,
                                 chetrs_obj->b, chetrs_obj->ldb );


    if( chetrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrs is wrong\n",
                    chetrs_obj->info );
    }
    if( chetrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrs is wrong\n",
        chetrs_obj->inforef );
    }
}

TEST_F(chetrs_test, chetrs1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_obj->b_bufsize,
                           chetrs_obj->b, chetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs_test, chetrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_obj->b_bufsize,
                           chetrs_obj->b, chetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs_test, chetrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_obj->b_bufsize,
                           chetrs_obj->b, chetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs_test, chetrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_obj->b_bufsize,
                           chetrs_obj->b, chetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin hetrs_dcomplex_parameters  class definition */
class hetrs_dcomplex_parameters{
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
      hetrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hetrs_dcomplex_parameters ();
};  /* end of hetrs_dcomplex_parameters  class definition */


/* Constructor hetrs_dcomplex_parameters definition */
hetrs_dcomplex_parameters:: hetrs_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hetrs DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hetrs_dcomplex_parameters object: malloc error.";
       hetrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hetrs_dcomplex_parameters:: ~hetrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetrs_free();
}


//  Test fixture class definition
class zhetrs_test  : public  ::testing::Test {
public:
   hetrs_dcomplex_parameters  *zhetrs_obj;
   void SetUp();  
   void TearDown () { delete zhetrs_obj; }
};


void zhetrs_test::SetUp(){

    /* LAPACKE ZHETRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zhetrs ZHETRS;

     /* LAPACKE ZHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zhetrf ZHETRF;


    zhetrs_obj = new  hetrs_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zhetrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHETRS = (Fptr_NL_LAPACKE_zhetrs)dlsym(zhetrs_obj->hModule, "LAPACKE_zhetrs");
    ASSERT_TRUE(ZHETRS != NULL) << "failed to get the Netlib LAPACKE_zhetrs symbol";

    ZHETRF = (Fptr_NL_LAPACKE_zhetrf)dlsym(zhetrs_obj->hModule,"LAPACKE_zhetrf");
    ASSERT_TRUE(ZHETRF != NULL) << "failed to get the Netlib LAPACKE_zhetrf symbol";

    /* Pre condition: need to call hetrf - before calling hetrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zhetrs_obj->inforef = ZHETRF( zhetrs_obj->matrix_layout,
                                    zhetrs_obj->uplo, zhetrs_obj->n,
                                     zhetrs_obj->aref,
                               zhetrs_obj->lda, zhetrs_obj->ipivref);

    zhetrs_obj->inforef = ZHETRS( zhetrs_obj->matrix_layout,
                                  zhetrs_obj->uplo, zhetrs_obj->n,
                                  zhetrs_obj->nrhs,
                                  (const lapack_complex_double *)zhetrs_obj->aref,
                                  zhetrs_obj->lda, zhetrs_obj->ipivref,
                                  zhetrs_obj->bref, zhetrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zhetrs_obj->info = LAPACKE_zhetrf( zhetrs_obj->matrix_layout,
                                 zhetrs_obj->uplo, zhetrs_obj->n,
                                     zhetrs_obj->a,
                               zhetrs_obj->lda, zhetrs_obj->ipiv);

    zhetrs_obj->info = LAPACKE_zhetrs( zhetrs_obj->matrix_layout,
                zhetrs_obj->uplo, zhetrs_obj->n, zhetrs_obj->nrhs,
                                  (const lapack_complex_double *)zhetrs_obj->a,
                               zhetrs_obj->lda, zhetrs_obj->ipiv,
                                 zhetrs_obj->b, zhetrs_obj->ldb );


    if( zhetrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrs is wrong\n",
                    zhetrs_obj->info );
    }
    if( zhetrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrs is wrong\n",
        zhetrs_obj->inforef );
    }
}

TEST_F(zhetrs_test, zhetrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_obj->b_bufsize,
                           zhetrs_obj->b, zhetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhetrs_test, zhetrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_obj->b_bufsize,
                           zhetrs_obj->b, zhetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhetrs_test, zhetrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_obj->b_bufsize,
                           zhetrs_obj->b, zhetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhetrs_test, zhetrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_obj->b_bufsize,
                           zhetrs_obj->b, zhetrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
