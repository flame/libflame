#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define LAPACKE_TEST_VERBOSE (1)
#define hptrs_free() \
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

/* Begin hptrs_scomplex_parameters  class definition */
class hptrs_scomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
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
      hptrs_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hptrs_scomplex_parameters ();
};  /* end of hptrs_scomplex_parameters  class definition */


/* Constructor hptrs_scomplex_parameters definition */
hptrs_scomplex_parameters:: hptrs_scomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n hptrs scomplex:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
             n, uplo, ldb, nrhs);
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
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "hptrs_scomplex_parameters object: malloc error.";
       hptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hptrs_scomplex_parameters:: ~hptrs_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hptrs_scomplex_parameters object: destructor invoked. \n");
#endif
   hptrs_free();
}


//  Test fixture class definition
class chptrs_test  : public  ::testing::Test {
public:
   hptrs_scomplex_parameters  *chptrs_obj;
   void SetUp();  
   void TearDown () { delete chptrs_obj; }
};


void chptrs_test::SetUp(){

    /* LAPACKE CHPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_chptrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_float *a,
                                   const lapack_int *ipiv,
                              lapack_complex_float *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_chptrs CHPTRS;

     /* LAPACKE CHPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_chptrf) ( int matrix_layout, char uplo,
    lapack_int n, lapack_complex_float* a,lapack_int* ipiv );

    Fptr_NL_LAPACKE_chptrf CHPTRF;


    chptrs_obj = new  hptrs_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    chptrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chptrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chptrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chptrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHPTRS = (Fptr_NL_LAPACKE_chptrs)dlsym(chptrs_obj->hModule, "LAPACKE_chptrs");
    ASSERT_TRUE(CHPTRS != NULL) << "failed to get the Netlib LAPACKE_chptrs symbol";

    CHPTRF = (Fptr_NL_LAPACKE_chptrf)dlsym(chptrs_obj->hModule,"LAPACKE_chptrf");
    ASSERT_TRUE(CHPTRF != NULL) << "failed to get the Netlib LAPACKE_chptrf symbol";

    /* Pre condition: need to call hptrf - before calling hptrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    chptrs_obj->inforef = CHPTRF( chptrs_obj->matrix_layout,
                            chptrs_obj->uplo, chptrs_obj->n,
                                     chptrs_obj->aref,
                       chptrs_obj->ipivref);

    chptrs_obj->inforef = CHPTRS( chptrs_obj->matrix_layout,
                                  chptrs_obj->uplo, chptrs_obj->n,
                                  chptrs_obj->nrhs,
                      (const lapack_complex_float *)chptrs_obj->aref,
                                  chptrs_obj->ipivref,
                                  chptrs_obj->bref, chptrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    chptrs_obj->info = LAPACKE_chptrf( chptrs_obj->matrix_layout,
                                 chptrs_obj->uplo, chptrs_obj->n,
                                 chptrs_obj->a, chptrs_obj->ipiv);

    chptrs_obj->info = LAPACKE_chptrs( chptrs_obj->matrix_layout,
                chptrs_obj->uplo, chptrs_obj->n, chptrs_obj->nrhs,
                      (const lapack_complex_float *)chptrs_obj->a,
                chptrs_obj->ipiv, chptrs_obj->b, chptrs_obj->ldb );


    if( chptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chptrs is wrong\n",
                    chptrs_obj->info );
    }
    if( chptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chptrs is wrong\n",
        chptrs_obj->inforef );
    }
}

TEST_F(chptrs_test, chptrs1) {
    float diff;
    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chptrs_obj->b_bufsize,
                           chptrs_obj->b, 
						   chptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chptrs_test, chptrs2) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chptrs_obj->b_bufsize,
                           chptrs_obj->b, 
						   chptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chptrs_test, chptrs3) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chptrs_obj->b_bufsize,
                           chptrs_obj->b, chptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(chptrs_test, chptrs4) {
    float diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_c( chptrs_obj->b_bufsize,
                           chptrs_obj->b, chptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin hptrs_dcomplex_parameters  class definition */
class hptrs_dcomplex_parameters{
   public:
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
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
      hptrs_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hptrs_dcomplex_parameters ();
};  /* end of hptrs_dcomplex_parameters  class definition */


/* Constructor hptrs_dcomplex_parameters definition */
hptrs_dcomplex_parameters:: hptrs_dcomplex_parameters ( int matrix_layout_i,
                         char uplo_i, lapack_int n_i, 
                                lapack_int nrhs_i, lapack_int ldb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    ldb = ldb_i;
    nrhs = nrhs_i;

    hModule = NULL;
    dModule = NULL;
#if LAPACKE_TEST_VERBOSE
   printf(" \n hptrs DComplex:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
             n, uplo, ldb, nrhs);
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
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (b==NULL) || (bref==NULL) ||  \
        (ipiv==NULL) || (ipivref==NULL) ){
       EXPECT_FALSE( true) << "hptrs_dcomplex_parameters object: malloc error.";
       hptrs_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hptrs_dcomplex_parameters:: ~hptrs_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hptrs_dcomplex_parameters object: destructor invoked. \n");
#endif
   hptrs_free();
}


//  Test fixture class definition
class zhptrs_test  : public  ::testing::Test {
public:
   hptrs_dcomplex_parameters  *zhptrs_obj;
   void SetUp();  
   void TearDown () { delete zhptrs_obj; }
};


void zhptrs_test::SetUp(){

    /* LAPACKE ZHPTRS prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptrs) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                          const lapack_complex_double *a,
                                   const lapack_int *ipiv,
                              lapack_complex_double *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zhptrs ZHPTRS;

     /* LAPACKE ZHPTRF prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptrf) ( int matrix_layout,char uplo ,lapack_int n,
                                    lapack_complex_double *a,lapack_int *ipiv );

    Fptr_NL_LAPACKE_zhptrf ZHPTRF;


    zhptrs_obj = new  hptrs_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zhptrs_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhptrs_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhptrs_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhptrs_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHPTRS = (Fptr_NL_LAPACKE_zhptrs)dlsym(zhptrs_obj->hModule, "LAPACKE_zhptrs");
    ASSERT_TRUE(ZHPTRS != NULL) << "failed to get the Netlib LAPACKE_zhptrs symbol";

    ZHPTRF = (Fptr_NL_LAPACKE_zhptrf)dlsym(zhptrs_obj->hModule,"LAPACKE_zhptrf");
    ASSERT_TRUE(ZHPTRF != NULL) << "failed to get the Netlib LAPACKE_zhptrf symbol";

    /* Pre condition: need to call hptrf - before calling hptrs function */

    /* Compute the Netlib-Lapacke's reference o/p */
    zhptrs_obj->inforef = ZHPTRF( zhptrs_obj->matrix_layout,
                                    zhptrs_obj->uplo, zhptrs_obj->n,
                                     zhptrs_obj->aref,
                               zhptrs_obj->ipivref);

    zhptrs_obj->inforef = ZHPTRS( zhptrs_obj->matrix_layout,
                                  zhptrs_obj->uplo, zhptrs_obj->n,
                                  zhptrs_obj->nrhs,
                                  (const lapack_complex_double *)zhptrs_obj->aref,
                                  zhptrs_obj->ipivref,
                                  zhptrs_obj->bref, zhptrs_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zhptrs_obj->info = LAPACKE_zhptrf( zhptrs_obj->matrix_layout,
                                 zhptrs_obj->uplo, zhptrs_obj->n,
                                     zhptrs_obj->a,
                               zhptrs_obj->ipiv);

    zhptrs_obj->info = LAPACKE_zhptrs( zhptrs_obj->matrix_layout,
                zhptrs_obj->uplo, zhptrs_obj->n, zhptrs_obj->nrhs,
                                  (const lapack_complex_double *)zhptrs_obj->a,
                               zhptrs_obj->ipiv,
                                 zhptrs_obj->b, zhptrs_obj->ldb );


    if( zhptrs_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhptrs is wrong\n",
                    zhptrs_obj->info );
    }
    if( zhptrs_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhptrs is wrong\n",
        zhptrs_obj->inforef );
    }
}

TEST_F(zhptrs_test, zhptrs1) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhptrs_obj->b_bufsize,
                           zhptrs_obj->b, zhptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhptrs_test, zhptrs2) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhptrs_obj->b_bufsize,
                           zhptrs_obj->b, zhptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhptrs_test, zhptrs3) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhptrs_obj->b_bufsize,
                           zhptrs_obj->b, zhptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(zhptrs_test, zhptrs4) {
    double diff;

    /* Compute Difference between libflame and Netlib o/ps  */
    diff =  computeDiff_z( zhptrs_obj->b_bufsize,
                           zhptrs_obj->b, zhptrs_obj->bref );

    EXPECT_NEAR(0.0, diff, LAPACKE_GTEST_THRESHOLD);
}
