#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define LAPACKE_TEST_VERBOSE (1)
#define hpsv_free() \
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

/* Begin hpsv_scomplex_parameters  class definition */
class hpsv_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      int b_bufsize;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      lapack_complex_float *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      hpsv_scomplex_parameters ( int matrix_layout_i, char uplo_i,
                       lapack_int n_i, lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hpsv_scomplex_parameters ();
};  /* end of hpsv_scomplex_parameters  class definition */


/* Constructor hpsv_scomplex_parameters definition */
hpsv_scomplex_parameters:: hpsv_scomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hpsv scomplex:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hpsv_scomplex_parameters object: malloc error.";
       hpsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, b_bufsize);
   

   } /* end of Constructor  */

hpsv_scomplex_parameters:: ~hpsv_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hpsv_scomplex_parameters object: destructor invoked. \n");
#endif
   hpsv_free();
}


//  Test fixture class definition
class chpsv_test  : public  ::testing::Test {
public:
   hpsv_scomplex_parameters  *chpsv_obj;
   void SetUp();  
   void TearDown () { delete chpsv_obj; }
};


void chpsv_test::SetUp(){

    /* LAPACKE CHPSV prototype */
    typedef int (*Fptr_NL_LAPACKE_chpsv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_float *a,
                                    lapack_int *ipiv,
                              lapack_complex_float *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_chpsv CHPSV;

    chpsv_obj = new  hpsv_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    chpsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chpsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chpsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chpsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHPSV = (Fptr_NL_LAPACKE_chpsv)dlsym(chpsv_obj->hModule, "LAPACKE_chpsv");
    ASSERT_TRUE(CHPSV != NULL) << "failed to get the Netlib LAPACKE_chpsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    chpsv_obj->inforef = CHPSV( chpsv_obj->matrix_layout,
                                  chpsv_obj->uplo, chpsv_obj->n,
                                  chpsv_obj->nrhs,
                      ( lapack_complex_float *)chpsv_obj->aref,
                                  chpsv_obj->ipivref,
                                  chpsv_obj->bref, chpsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    chpsv_obj->info = LAPACKE_chpsv( chpsv_obj->matrix_layout,
                chpsv_obj->uplo, chpsv_obj->n, chpsv_obj->nrhs,
                      ( lapack_complex_float *)chpsv_obj->a,
                chpsv_obj->ipiv, chpsv_obj->b, chpsv_obj->ldb );


    if( chpsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chpsv is wrong\n",
                    chpsv_obj->info );
    }
    if( chpsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chpsv is wrong\n",
        chpsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    chpsv_obj->diff =  computeDiff_c( chpsv_obj->b_bufsize,
                           chpsv_obj->b, 
						   chpsv_obj->bref );
}

TEST_F(chpsv_test, chpsv1) {
    EXPECT_NEAR(0.0, chpsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpsv_test, chpsv2) {
    EXPECT_NEAR(0.0, chpsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpsv_test, chpsv3) {
    EXPECT_NEAR(0.0, chpsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chpsv_test, chpsv4) {
    EXPECT_NEAR(0.0, chpsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin hpsv_dcomplex_parameters  class definition */
class hpsv_dcomplex_parameters{
   public:
      int b_bufsize;
      double diff; //to capture reference and libflame o/ps' difference
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A
      lapack_int *ipiv, *ipivref; // The pivot indices
      void *hModule, *dModule;

      /* Input/ Output parameters */
      lapack_complex_double *b, *bref; //right-hand sides for the systems of equations.

      /* Return Values */
      lapack_int info, inforef;

   public:
      hpsv_dcomplex_parameters ( int matrix_layout_i, char uplo_i,
                                lapack_int n_i, 
                            lapack_int nrhs_i, lapack_int ldb_i);
             
      ~hpsv_dcomplex_parameters ();
};  /* end of hpsv_dcomplex_parameters  class definition */


/* Constructor hpsv_dcomplex_parameters definition */
hpsv_dcomplex_parameters:: hpsv_dcomplex_parameters ( int matrix_layout_i,
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
   printf(" \n hpsv DComplex:  n: %d, uplo: %c  ldb: %d nrhs: %d \n",
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
       EXPECT_FALSE( true) << "hpsv_dcomplex_parameters object: malloc error.";
       hpsv_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, (n*(n+1)/2));
    lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, b_bufsize);

   } /* end of Constructor  */

hpsv_dcomplex_parameters:: ~hpsv_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hpsv_dcomplex_parameters object: destructor invoked. \n");
#endif
   hpsv_free();
}


//  Test fixture class definition
class zhpsv_test  : public  ::testing::Test {
public:
   hpsv_dcomplex_parameters  *zhpsv_obj;
   void SetUp();  
   void TearDown () { delete zhpsv_obj; }
};


void zhpsv_test::SetUp(){

    /* LAPACKE ZHPSV prototype */
    typedef int (*Fptr_NL_LAPACKE_zhpsv) ( int matrix_layout, char uplo,
                          lapack_int n, lapack_int nrhs,
                           lapack_complex_double *a,
                                    lapack_int *ipiv,
                              lapack_complex_double *b, lapack_int ldb  );

    Fptr_NL_LAPACKE_zhpsv ZHPSV;
    zhpsv_obj = new  hpsv_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
                                         lin_solver_paramslist[idx].nrhs,
                                         lin_solver_paramslist[idx].ldb );
    idx = Circular_Increment_Index(idx);


    zhpsv_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhpsv_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhpsv_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhpsv_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHPSV = (Fptr_NL_LAPACKE_zhpsv)dlsym(zhpsv_obj->hModule, "LAPACKE_zhpsv");
    ASSERT_TRUE(ZHPSV != NULL) << "failed to get the Netlib LAPACKE_zhpsv symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zhpsv_obj->inforef = ZHPSV( zhpsv_obj->matrix_layout,
                                  zhpsv_obj->uplo, zhpsv_obj->n,
                                  zhpsv_obj->nrhs,
                                  zhpsv_obj->aref,
                                  zhpsv_obj->ipivref,
                                  zhpsv_obj->bref, zhpsv_obj->ldb);

    /* Compute the Libflme lapacke o/p  */
    zhpsv_obj->info = LAPACKE_zhpsv( zhpsv_obj->matrix_layout,
									 zhpsv_obj->uplo, zhpsv_obj->n,
									 zhpsv_obj->nrhs, zhpsv_obj->a,
									 zhpsv_obj->ipiv,
                                     zhpsv_obj->b, zhpsv_obj->ldb );

    if( zhpsv_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhpsv is wrong\n",
                    zhpsv_obj->info );
    }
    if( zhpsv_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhpsv is wrong\n",
        zhpsv_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zhpsv_obj->diff =  computeDiff_z( zhpsv_obj->b_bufsize,
                           zhpsv_obj->b, zhpsv_obj->bref );
}

TEST_F(zhpsv_test, zhpsv1) {
    EXPECT_NEAR(0.0, zhpsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpsv_test, zhpsv2) {
    EXPECT_NEAR(0.0, zhpsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpsv_test, zhpsv3) {
    EXPECT_NEAR(0.0, zhpsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhpsv_test, zhpsv4) {
    EXPECT_NEAR(0.0, zhpsv_obj->diff, LAPACKE_GTEST_THRESHOLD);
}