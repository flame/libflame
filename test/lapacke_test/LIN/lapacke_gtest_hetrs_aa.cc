#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define hetrs_aa_free() \
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


/* Begin hetrs_aa_scomplex_parameters  class definition */
class hetrs_aa_scomplex_parameters{
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
      hetrs_aa_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hetrs_aa_scomplex_parameters ();
};  /* end of hetrs_aa_scomplex_parameters  class definition */


/* Constructor hetrs_aa_scomplex_parameters definition */
hetrs_aa_scomplex_parameters:: hetrs_aa_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hetrs_aa scomplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hetrs_aa_scomplex_parameters object: malloc error.";
       hetrs_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hetrs_aa_scomplex_parameters:: ~hetrs_aa_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrs_aa_scomplex_parameters object: destructor invoked. \n");
#endif
   hetrs_aa_free();
}


//  Test fixture class definition
class chetrs_aa_test  : public  ::testing::Test {
public:
   hetrs_aa_scomplex_parameters  *chetrs_aa_obj;
   void SetUp();  
   void TearDown () { delete chetrs_aa_obj; }
};


void chetrs_aa_test::SetUp(){

    /* LAPACKE CHETRS_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrs_aa) ( int matrix_layout, char uplo,
                                           lapack_int n, lapack_int nrhs,
                                          const lapack_complex_float * a,
                                 lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_float * b, lapack_int ldb );

    Fptr_NL_LAPACKE_chetrs_aa CHETRS_AA;

     /* LAPACKE CHETRF_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf_aa) ( int matrix_layout, char uplo,
    lapack_int n, lapack_complex_float* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_chetrf_aa CHETRF_AA;


    chetrs_aa_obj = new  hetrs_aa_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    chetrs_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetrs_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetrs_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetrs_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHETRS_AA = (Fptr_NL_LAPACKE_chetrs_aa)dlsym(chetrs_aa_obj->hModule, "LAPACKE_chetrs_aa");
    ASSERT_TRUE(CHETRS_AA != NULL) << "failed to get the Netlib LAPACKE_chetrs_aa symbol";

    CHETRF_AA = (Fptr_NL_LAPACKE_chetrf_aa)dlsym(chetrs_aa_obj->hModule,"LAPACKE_chetrf_aa");
    ASSERT_TRUE(CHETRF_AA != NULL) << "failed to get the Netlib LAPACKE_chetrf_aa symbol";

    /* Pre condition: need to call hetrf_aa - before calling hetrs_aa function */

    /* Compute the Netlib-Lapacke's reference o/p */
    chetrs_aa_obj->inforef = CHETRF_AA( chetrs_aa_obj->matrix_layout,
                            chetrs_aa_obj->uplo, chetrs_aa_obj->n,
                                           chetrs_aa_obj->aref,
                      chetrs_aa_obj->lda, chetrs_aa_obj->ipivref);

    chetrs_aa_obj->inforef = CHETRS_AA( chetrs_aa_obj->matrix_layout,
                                  chetrs_aa_obj->uplo, chetrs_aa_obj->n,
                                  chetrs_aa_obj->nrhs,
                                  (const lapack_complex_float *)chetrs_aa_obj->aref,
                                  chetrs_aa_obj->lda, chetrs_aa_obj->ipivref,
                                  chetrs_aa_obj->bref, chetrs_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    chetrs_aa_obj->info = LAPACKE_chetrf_aa( chetrs_aa_obj->matrix_layout,
                                 chetrs_aa_obj->uplo, chetrs_aa_obj->n,
                                     chetrs_aa_obj->a,
                               chetrs_aa_obj->lda, chetrs_aa_obj->ipiv);

    chetrs_aa_obj->info = LAPACKE_chetrs_aa( chetrs_aa_obj->matrix_layout,
                chetrs_aa_obj->uplo, chetrs_aa_obj->n, chetrs_aa_obj->nrhs,
                                  (const lapack_complex_float *)chetrs_aa_obj->a,
                               chetrs_aa_obj->lda, chetrs_aa_obj->ipiv,
                                 chetrs_aa_obj->b, chetrs_aa_obj->ldb );


    if( chetrs_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrs_aa is wrong\n",
                    chetrs_aa_obj->info );
    }
    if( chetrs_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrs_aa is wrong\n",
        chetrs_aa_obj->inforef );
    }
}

TEST_F(chetrs_aa_test, chetrs_aa1) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_aa_obj->b_bufsize,
                           chetrs_aa_obj->b, chetrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs_aa_test, chetrs_aa2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_aa_obj->b_bufsize,
                           chetrs_aa_obj->b, chetrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs_aa_test, chetrs_aa3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_aa_obj->b_bufsize,
                           chetrs_aa_obj->b, chetrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chetrs_aa_test, chetrs_aa4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chetrs_aa_obj->b_bufsize,
                           chetrs_aa_obj->b, chetrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin hetrs_aa_dcomplex_parameters  class definition */
class hetrs_aa_dcomplex_parameters{
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
      hetrs_aa_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, lapack_int lda_i,
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hetrs_aa_dcomplex_parameters ();
};  /* end of hetrs_aa_dcomplex_parameters  class definition */


/* Constructor hetrs_aa_dcomplex_parameters definition */
hetrs_aa_dcomplex_parameters:: hetrs_aa_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hetrs_aa DComplex:  n: %d, uplo: %c lda: %d ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hetrs_aa_dcomplex_parameters object: malloc error.";
       hetrs_aa_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*lda);
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hetrs_aa_dcomplex_parameters:: ~hetrs_aa_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetrs_aa_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetrs_aa_free();
}


//  Test fixture class definition
class zhetrs_aa_test  : public  ::testing::Test {
public:
   hetrs_aa_dcomplex_parameters  *zhetrs_aa_obj;
   void SetUp();  
   void TearDown () { delete zhetrs_aa_obj; }
};


void zhetrs_aa_test::SetUp(){

    /* LAPACKE ZHETRS_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrs_aa) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_double * a,
                                  lapack_int lda, const lapack_int * ipiv,
                              lapack_complex_double * b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zhetrs_aa ZHETRS_AA;

     /* LAPACKE ZHETRF_AA prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf_aa) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_double* a,lapack_int lda,lapack_int* ipiv );

    Fptr_NL_LAPACKE_zhetrf_aa ZHETRF_AA;


    zhetrs_aa_obj = new  hetrs_aa_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].lda,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zhetrs_aa_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetrs_aa_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetrs_aa_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetrs_aa_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHETRS_AA = (Fptr_NL_LAPACKE_zhetrs_aa)dlsym(zhetrs_aa_obj->hModule, "LAPACKE_zhetrs_aa");
    ASSERT_TRUE(ZHETRS_AA != NULL) << "failed to get the Netlib LAPACKE_zhetrs_aa symbol";

    ZHETRF_AA = (Fptr_NL_LAPACKE_zhetrf_aa)dlsym(zhetrs_aa_obj->hModule,"LAPACKE_zhetrf_aa");
    ASSERT_TRUE(ZHETRF_AA != NULL) << "failed to get the Netlib LAPACKE_zhetrf_aa symbol";

    /* Pre condition: need to call hetrf_aa - before calling hetrs_aa function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zhetrs_aa_obj->inforef = ZHETRF_AA( zhetrs_aa_obj->matrix_layout,
                                    zhetrs_aa_obj->uplo, zhetrs_aa_obj->n,
                                    zhetrs_aa_obj->aref,
                               zhetrs_aa_obj->lda, zhetrs_aa_obj->ipivref);

    zhetrs_aa_obj->inforef = ZHETRS_AA( zhetrs_aa_obj->matrix_layout,
                                  zhetrs_aa_obj->uplo, zhetrs_aa_obj->n,
                                  zhetrs_aa_obj->nrhs,
                          (const lapack_complex_double *)zhetrs_aa_obj->aref,
                                  zhetrs_aa_obj->lda, zhetrs_aa_obj->ipivref,
                                  zhetrs_aa_obj->bref, zhetrs_aa_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zhetrs_aa_obj->info = LAPACKE_zhetrf_aa( zhetrs_aa_obj->matrix_layout,
                                    zhetrs_aa_obj->uplo, zhetrs_aa_obj->n,
                                    zhetrs_aa_obj->a,
                                 zhetrs_aa_obj->lda, zhetrs_aa_obj->ipiv);

    zhetrs_aa_obj->info = LAPACKE_zhetrs_aa( zhetrs_aa_obj->matrix_layout,
                zhetrs_aa_obj->uplo, zhetrs_aa_obj->n, zhetrs_aa_obj->nrhs,
                           (const lapack_complex_double *)zhetrs_aa_obj->a,
                               zhetrs_aa_obj->lda, zhetrs_aa_obj->ipiv,
                                 zhetrs_aa_obj->b, zhetrs_aa_obj->ldb );


    if( zhetrs_aa_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrs_aa is wrong\n",
                    zhetrs_aa_obj->info );
    }
    if( zhetrs_aa_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrs_aa is wrong\n",
        zhetrs_aa_obj->inforef );
    }
}

TEST_F(zhetrs_aa_test, zhetrs_aa1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_aa_obj->b_bufsize,
                           zhetrs_aa_obj->b, zhetrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhetrs_aa_test, zhetrs_aa2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_aa_obj->b_bufsize,
                           zhetrs_aa_obj->b, zhetrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhetrs_aa_test, zhetrs_aa3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_aa_obj->b_bufsize,
                           zhetrs_aa_obj->b, zhetrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhetrs_aa_test, zhetrs_aa4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhetrs_aa_obj->b_bufsize,
                           zhetrs_aa_obj->b, zhetrs_aa_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
