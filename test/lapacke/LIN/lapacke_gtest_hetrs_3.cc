#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define hetrs_3_free() \
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


/* Begin hetrs_3_scomplex_parameters  class definition */
class hetrs_3_scomplex_parameters{
   public:
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
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      hetrs_3_scomplex_parameters( int matrix_layout_i, char uplo_i,
                                   lapack_int n_i, lapack_int lda_i,
                               lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hetrs_3_scomplex_parameters ();
};  /* end of hetrs_3_scomplex_parameters  class definition */


/* Constructor hetrs_3_scomplex_parameters definition */
hetrs_3_scomplex_parameters:: hetrs_3_scomplex_parameters( int matrix_layout_i,
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
   printf(" \n hetrs_3 scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hetrs_3_scomplex_parameters object: malloc error.";
       hetrs_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
    lapacke_gtest_init_scomplex_buffer_pair_rand( e, eref, n );

   } /* end of Constructor  */

hetrs_3_scomplex_parameters:: ~hetrs_3_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrs_3_scomplex_parameters object: destructor invoked. \n");
#endif
   hetrs_3_free();
}


//  Test fixture class definition
class chetrs_3_test  : public  ::testing::Test {
public:
   hetrs_3_scomplex_parameters  *chetrs_3_obj;
   void SetUp();  
   void TearDown () { delete chetrs_3_obj; }
};


void chetrs_3_test::SetUp(){

    /* LAPACKE CHETRS_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrs_3) ( int matrix_layout, char uplo,
                                           lapack_int n, lapack_int nrhs,
                                           const lapack_complex_float* a,
                           lapack_int lda, const lapack_complex_float *e,
                                                 const lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_chetrs_3 CHETRS_3;

     /* LAPACKE CHETRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf_rk) ( int matrix_layout, char uplo,
                    lapack_int n, lapack_complex_float* a,lapack_int lda, 
                             lapack_complex_float* e, lapack_int* ipiv );

    Fptr_NL_LAPACKE_chetrf_rk CHETRF_RK;


    chetrs_3_obj = new  hetrs_3_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    chetrs_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetrs_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetrs_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetrs_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHETRS_3 = (Fptr_NL_LAPACKE_chetrs_3)dlsym(chetrs_3_obj->hModule, "LAPACKE_chetrs_3");
    ASSERT_TRUE(CHETRS_3 != NULL) << "failed to get the Netlib LAPACKE_chetrs_3 symbol";

    CHETRF_RK = (Fptr_NL_LAPACKE_chetrf_rk)dlsym(chetrs_3_obj->hModule,"LAPACKE_chetrf_rk");
    ASSERT_TRUE(CHETRF_RK != NULL) << "failed to get the Netlib LAPACKE_chetrf_rk symbol";

    /* Pre condition: need to call hetrf_rk - before calling hetrs_3 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    chetrs_3_obj->inforef = CHETRF_RK( chetrs_3_obj->matrix_layout,
                                       chetrs_3_obj->uplo,
                                       chetrs_3_obj->n,
                                       chetrs_3_obj->aref,
                                       chetrs_3_obj->lda,
                                       chetrs_3_obj->eref, 
                                       chetrs_3_obj->ipivref);

    chetrs_3_obj->inforef = CHETRS_3( chetrs_3_obj->matrix_layout,
                                      chetrs_3_obj->uplo,
                                      chetrs_3_obj->n,
                                      chetrs_3_obj->nrhs,
                        (const lapack_complex_float *)chetrs_3_obj->aref,
                                      chetrs_3_obj->lda, 
                        (const lapack_complex_float *)chetrs_3_obj->eref,
                                      chetrs_3_obj->ipivref,
                                      chetrs_3_obj->bref,
                                      chetrs_3_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    chetrs_3_obj->info = LAPACKE_chetrf_rk( chetrs_3_obj->matrix_layout,
                                            chetrs_3_obj->uplo,
                                            chetrs_3_obj->n,
                                            chetrs_3_obj->a,
                                            chetrs_3_obj->lda,
                                            chetrs_3_obj->e,
                                            chetrs_3_obj->ipiv);

    chetrs_3_obj->info = LAPACKE_chetrs_3( chetrs_3_obj->matrix_layout,
                                           chetrs_3_obj->uplo, 
                                           chetrs_3_obj->n, 
                                           chetrs_3_obj->nrhs,
                              (const lapack_complex_float *)chetrs_3_obj->a,
                                           chetrs_3_obj->lda,
                              (const lapack_complex_float *)chetrs_3_obj->e,
                                           chetrs_3_obj->ipiv, 
                                           chetrs_3_obj->b,
                                           chetrs_3_obj->ldb );


    if( chetrs_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrs_3 is wrong\n",
                    chetrs_3_obj->info );
    }
    if( chetrs_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrs_3 is wrong\n",
        chetrs_3_obj->inforef );
    }
}

TEST_F(chetrs_3_test, chetrs_31) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_3_obj->b_bufsize,
                           chetrs_3_obj->b, chetrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs_3_test, chetrs_32) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_3_obj->b_bufsize,
                           chetrs_3_obj->b, chetrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs_3_test, chetrs_33) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_3_obj->b_bufsize,
                           chetrs_3_obj->b, chetrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs_3_test, chetrs_34) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_3_obj->b_bufsize,
                           chetrs_3_obj->b, chetrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin hetrs_3_dcomplex_parameters  class definition */
class hetrs_3_dcomplex_parameters{
   public:
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
      int b_bufsize;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      hetrs_3_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                    lapack_int n_i, lapack_int lda_i,
                                    lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hetrs_3_dcomplex_parameters ();
};  /* end of hetrs_3_dcomplex_parameters  class definition */


/* Constructor hetrs_3_dcomplex_parameters definition */
hetrs_3_dcomplex_parameters:: hetrs_3_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hetrs_3 DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hetrs_3_dcomplex_parameters object: malloc error.";
       hetrs_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( e, eref, n);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hetrs_3_dcomplex_parameters:: ~hetrs_3_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrs_3_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetrs_3_free();
}


//  Test fixture class definition
class zhetrs_3_test  : public  ::testing::Test {
public:
   hetrs_3_dcomplex_parameters  *zhetrs_3_obj;
   void SetUp();  
   void TearDown () { delete zhetrs_3_obj; }
};


void zhetrs_3_test::SetUp(){

    /* LAPACKE ZHETRS_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrs_3)(int matrix_layout, char uplo,
                                            lapack_int n, lapack_int nrhs,
                                            const lapack_complex_double * a,
                                            lapack_int lda, 
                                            const lapack_complex_double * e,
                                            const lapack_int * ipiv,
                                            lapack_complex_double * b,
                                            lapack_int ldb );

    Fptr_NL_LAPACKE_zhetrs_3 ZHETRS_3;

     /* LAPACKE ZHETRF_RK prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf_rk) ( int matrix_layout,char uplo ,lapack_int n,
   lapack_complex_double* a,lapack_int lda, lapack_complex_double* e, lapack_int* ipiv );

    Fptr_NL_LAPACKE_zhetrf_rk ZHETRF_RK;


    zhetrs_3_obj = new  hetrs_3_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );

    idx = Circular_Increment_Index(idx);

    zhetrs_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetrs_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetrs_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetrs_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHETRS_3 = (Fptr_NL_LAPACKE_zhetrs_3)dlsym(zhetrs_3_obj->hModule, "LAPACKE_zhetrs_3");
    ASSERT_TRUE(ZHETRS_3 != NULL) << "failed to get the Netlib LAPACKE_zhetrs_3 symbol";

    ZHETRF_RK = (Fptr_NL_LAPACKE_zhetrf_rk)dlsym(zhetrs_3_obj->hModule,"LAPACKE_zhetrf_rk");
    ASSERT_TRUE(ZHETRF_RK != NULL) << "failed to get the Netlib LAPACKE_zhetrf_rk symbol";

    /* Pre condition: need to call hetrf_rk - before calling hetrs_3 function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zhetrs_3_obj->inforef = ZHETRF_RK( zhetrs_3_obj->matrix_layout,
                                       zhetrs_3_obj->uplo, 
                                       zhetrs_3_obj->n,
                                       zhetrs_3_obj->aref,
                                       zhetrs_3_obj->lda, 
                                       zhetrs_3_obj->eref,
                                       zhetrs_3_obj->ipivref );

    zhetrs_3_obj->inforef = ZHETRS_3( zhetrs_3_obj->matrix_layout,
                                      zhetrs_3_obj->uplo, 
                                      zhetrs_3_obj->n,
                                      zhetrs_3_obj->nrhs,
                          (const lapack_complex_double *)zhetrs_3_obj->aref,
                                      zhetrs_3_obj->lda,
                          (const lapack_complex_double *)zhetrs_3_obj->eref,
                                      zhetrs_3_obj->ipivref,
                                      zhetrs_3_obj->bref,
                                      zhetrs_3_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zhetrs_3_obj->info = LAPACKE_zhetrf_rk(zhetrs_3_obj->matrix_layout,
                                           zhetrs_3_obj->uplo,
                                           zhetrs_3_obj->n,
                                           zhetrs_3_obj->a,
                                           zhetrs_3_obj->lda,
                                           zhetrs_3_obj->e,
                                           zhetrs_3_obj->ipiv);

    zhetrs_3_obj->info = LAPACKE_zhetrs_3(zhetrs_3_obj->matrix_layout,
                                          zhetrs_3_obj->uplo,
                                          zhetrs_3_obj->n,
                                          zhetrs_3_obj->nrhs,
                                  (const lapack_complex_double *)zhetrs_3_obj->a,
                                          zhetrs_3_obj->lda,
                                  (const lapack_complex_double *)zhetrs_3_obj->e,
                                          zhetrs_3_obj->ipiv,
                                          zhetrs_3_obj->b,
                                          zhetrs_3_obj->ldb );


    if( zhetrs_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrs_3 is wrong\n",
                    zhetrs_3_obj->info );
    }
    if( zhetrs_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrs_3 is wrong\n",
        zhetrs_3_obj->inforef );
    }
}

TEST_F(zhetrs_3_test, zhetrs_31) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_3_obj->b_bufsize,
                           zhetrs_3_obj->b,
                           zhetrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetrs_3_test, zhetrs_32) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_3_obj->b_bufsize,
                           zhetrs_3_obj->b,
                           zhetrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetrs_3_test, zhetrs_33) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_3_obj->b_bufsize,
                           zhetrs_3_obj->b,
                           zhetrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetrs_3_test, zhetrs_34) {
    double diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_3_obj->b_bufsize,
                           zhetrs_3_obj->b,
                           zhetrs_3_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
