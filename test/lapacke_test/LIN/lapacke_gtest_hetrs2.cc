#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define hetrs2_free() \
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


/* Begin hetrs2_scomplex_parameters  class definition */
class hetrs2_scomplex_parameters{
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
      hetrs2_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hetrs2_scomplex_parameters ();
};  /* end of hetrs2_scomplex_parameters  class definition */


/* Constructor hetrs2_scomplex_parameters definition */
hetrs2_scomplex_parameters:: hetrs2_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hetrs2 scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hetrs2_scomplex_parameters object: malloc error.";
       hetrs2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hetrs2_scomplex_parameters:: ~hetrs2_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrs2_scomplex_parameters object: destructor invoked. \n");
#endif
   hetrs2_free();
}


//  Test fixture class definition
class chetrs2_test  : public  ::testing::Test {
public:
   hetrs2_scomplex_parameters  *chetrs2_obj;
   void SetUp();  
   void TearDown () { delete chetrs2_obj; }
};


void chetrs2_test::SetUp(){

    /* LAPACKE CHETRS2 prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrs2) ( int matrix_layout, char uplo,
                                           lapack_int n, lapack_int nrhs,
                                          const lapack_complex_float * a,
                                 lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_chetrs2 CHETRS2;

     /* LAPACKE CHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf) ( int matrix_layout, char uplo,
    lapack_int n, lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_chetrf CHETRF;


    chetrs2_obj = new  hetrs2_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    chetrs2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetrs2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetrs2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetrs2_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHETRS2 = (Fptr_NL_LAPACKE_chetrs2)dlsym(chetrs2_obj->hModule, "LAPACKE_chetrs2");
    ASSERT_TRUE(CHETRS2 != NULL) << "failed to get the Netlib LAPACKE_chetrs2 symbol";

    CHETRF = (Fptr_NL_LAPACKE_chetrf)dlsym(chetrs2_obj->hModule,"LAPACKE_chetrf");
    ASSERT_TRUE(CHETRF != NULL) << "failed to get the Netlib LAPACKE_chetrf symbol";

    /* Pre condition: need to call hetrf - before calling hetrs2 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    chetrs2_obj->inforef = CHETRF( chetrs2_obj->matrix_layout,
                            chetrs2_obj->uplo, chetrs2_obj->n,
                                           chetrs2_obj->aref,
                      chetrs2_obj->lda, chetrs2_obj->ipivref);

    chetrs2_obj->inforef = CHETRS2( chetrs2_obj->matrix_layout,
                                  chetrs2_obj->uplo, chetrs2_obj->n,
                                  chetrs2_obj->nrhs,
                                  (const lapack_complex_float *)chetrs2_obj->aref,
                                  chetrs2_obj->lda, chetrs2_obj->ipivref,
                                  chetrs2_obj->bref, chetrs2_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    chetrs2_obj->info = LAPACKE_chetrf( chetrs2_obj->matrix_layout,
                                 chetrs2_obj->uplo, chetrs2_obj->n,
                                     chetrs2_obj->a,
                               chetrs2_obj->lda, chetrs2_obj->ipiv);

    chetrs2_obj->info = LAPACKE_chetrs2( chetrs2_obj->matrix_layout,
                chetrs2_obj->uplo, chetrs2_obj->n, chetrs2_obj->nrhs,
                                  (const lapack_complex_float *)chetrs2_obj->a,
                               chetrs2_obj->lda, chetrs2_obj->ipiv,
                                 chetrs2_obj->b, chetrs2_obj->ldb );


    if( chetrs2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrs2 is wrong\n",
                    chetrs2_obj->info );
    }
    if( chetrs2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrs2 is wrong\n",
        chetrs2_obj->inforef );
    }
}

TEST_F(chetrs2_test, chetrs21) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs2_obj->b_bufsize,
                           chetrs2_obj->b, chetrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs2_test, chetrs22) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs2_obj->b_bufsize,
                           chetrs2_obj->b, chetrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs2_test, chetrs23) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs2_obj->b_bufsize,
                           chetrs2_obj->b, chetrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs2_test, chetrs24) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs2_obj->b_bufsize,
                           chetrs2_obj->b, chetrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin hetrs2_dcomplex_parameters  class definition */
class hetrs2_dcomplex_parameters{
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
      hetrs2_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hetrs2_dcomplex_parameters ();
};  /* end of hetrs2_dcomplex_parameters  class definition */


/* Constructor hetrs2_dcomplex_parameters definition */
hetrs2_dcomplex_parameters:: hetrs2_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hetrs2 DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hetrs2_dcomplex_parameters object: malloc error.";
       hetrs2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hetrs2_dcomplex_parameters:: ~hetrs2_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrs2_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetrs2_free();
}


//  Test fixture class definition
class zhetrs2_test  : public  ::testing::Test {
public:
   hetrs2_dcomplex_parameters  *zhetrs2_obj;
   void SetUp();  
   void TearDown () { delete zhetrs2_obj; }
};


void zhetrs2_test::SetUp(){

    /* LAPACKE ZHETRS2 prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrs2) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zhetrs2 ZHETRS2;

     /* LAPACKE ZHETRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zhetrf ZHETRF;


    zhetrs2_obj = new  hetrs2_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zhetrs2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetrs2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetrs2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetrs2_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHETRS2 = (Fptr_NL_LAPACKE_zhetrs2)dlsym(zhetrs2_obj->hModule, "LAPACKE_zhetrs2");
    ASSERT_TRUE(ZHETRS2 != NULL) << "failed to get the Netlib LAPACKE_zhetrs2 symbol";

    ZHETRF = (Fptr_NL_LAPACKE_zhetrf)dlsym(zhetrs2_obj->hModule,"LAPACKE_zhetrf");
    ASSERT_TRUE(ZHETRF != NULL) << "failed to get the Netlib LAPACKE_zhetrf symbol";

    /* Pre condition: need to call hetrf - before calling hetrs2 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zhetrs2_obj->inforef = ZHETRF( zhetrs2_obj->matrix_layout,
                                    zhetrs2_obj->uplo, zhetrs2_obj->n,
                                     zhetrs2_obj->aref,
                               zhetrs2_obj->lda, zhetrs2_obj->ipivref);

    zhetrs2_obj->inforef = ZHETRS2( zhetrs2_obj->matrix_layout,
                                  zhetrs2_obj->uplo, zhetrs2_obj->n,
                                  zhetrs2_obj->nrhs,
                                  (const lapack_complex_double *)zhetrs2_obj->aref,
                                  zhetrs2_obj->lda, zhetrs2_obj->ipivref,
                                  zhetrs2_obj->bref, zhetrs2_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zhetrs2_obj->info = LAPACKE_zhetrf( zhetrs2_obj->matrix_layout,
                                 zhetrs2_obj->uplo, zhetrs2_obj->n,
                                     zhetrs2_obj->a,
                               zhetrs2_obj->lda, zhetrs2_obj->ipiv);

    zhetrs2_obj->info = LAPACKE_zhetrs2( zhetrs2_obj->matrix_layout,
                zhetrs2_obj->uplo, zhetrs2_obj->n, zhetrs2_obj->nrhs,
                                  (const lapack_complex_double *)zhetrs2_obj->a,
                               zhetrs2_obj->lda, zhetrs2_obj->ipiv,
                                 zhetrs2_obj->b, zhetrs2_obj->ldb );


    if( zhetrs2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrs2 is wrong\n",
                    zhetrs2_obj->info );
    }
    if( zhetrs2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrs2 is wrong\n",
        zhetrs2_obj->inforef );
    }
}

TEST_F(zhetrs2_test, zhetrs21) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs2_obj->b_bufsize,
                           zhetrs2_obj->b, zhetrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhetrs2_test, zhetrs22) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs2_obj->b_bufsize,
                           zhetrs2_obj->b, zhetrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhetrs2_test, zhetrs23) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs2_obj->b_bufsize,
                           zhetrs2_obj->b, zhetrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhetrs2_test, zhetrs24) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs2_obj->b_bufsize,
                           zhetrs2_obj->b, zhetrs2_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
