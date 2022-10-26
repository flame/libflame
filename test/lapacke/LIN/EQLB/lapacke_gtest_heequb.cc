#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define heequb_free() \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (s != NULL)    free (s   ); \
       if (sref != NULL) free (sref); \
       if( hModule != NULL) dlclose(hModule); \
       if( dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin heequb_scomplex_parameters  class definition */
class heequb_scomplex_parameters{
   public:
      int a_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // order of matrix A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_complex_float *a, *aref; // the n-by-n symmetric / Hermitian positive definite matrix
      char uplo; //  Must be 'U' or 'L'.

      /* Output parameters */
      float * s, *sref; // colum scale factors, , array if size 'n'
      float scond, scondref; // ratio of the smallest c[i] to the largest c[i]
      float amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      heequb_scomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i, char uplo_i);
             
      ~heequb_scomplex_parameters ();
};  /* end of heequb_scomplex_parameters  class definition */


/* Constructor heequb_scomplex_parameters definition */
heequb_scomplex_parameters:: heequb_scomplex_parameters ( int matrix_layout_i, 
     lapack_int n_i, char uplo_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n heequb scomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = n*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       heequb_free();
       EXPECT_FALSE( true) << "heequb_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref,
                                                             n, n, 'S');
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

heequb_scomplex_parameters:: ~heequb_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" heequb_scomplex_parameters object: destructor invoked. \n");
#endif
   heequb_free();
}


//  Test fixture class definition
class cheequb_test  : public  ::testing::Test {
public:
   heequb_scomplex_parameters  *cheequb_obj;
   void SetUp();  
   void TearDown () { delete cheequb_obj; }
};


void cheequb_test::SetUp(){

    /* LAPACKE CHEEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_cheequb) ( int matrix_layout, char uplo, 
                            lapack_int n, const lapack_complex_float *a, 
                 lapack_int lda, float *s, float *scond, float *amax  );

    Fptr_NL_LAPACKE_cheequb CHEEQUB;

    cheequb_obj = new  heequb_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n, 
                                         lin_solver_paramslist[idx].Uplo);
    idx = Circular_Increment_Index(idx);

    cheequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cheequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cheequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cheequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHEEQUB = (Fptr_NL_LAPACKE_cheequb)dlsym(cheequb_obj->hModule, "LAPACKE_cheequb");
    ASSERT_TRUE(CHEEQUB != NULL) << "failed to get the Netlib LAPACKE_cheequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    cheequb_obj->inforef = CHEEQUB( cheequb_obj->matrix_layout,
                                  cheequb_obj->uplo,
                                  cheequb_obj->n,
                                  cheequb_obj->aref,
                                  cheequb_obj->lda,
                                  cheequb_obj->sref,
                                  &cheequb_obj->scondref,
                                  &cheequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    cheequb_obj->info = LAPACKE_cheequb( cheequb_obj->matrix_layout,
                                  cheequb_obj->uplo,
                                  cheequb_obj->n,
                                  cheequb_obj->a,
                                  cheequb_obj->lda,
                                  cheequb_obj->s,
                                  &cheequb_obj->scond,
                                  &cheequb_obj->amax);

    if( cheequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cheequb \
        is wrong\n", cheequb_obj->info );
    }
    if( cheequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cheequb is wrong\n",
        cheequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cheequb_obj->diff =  computeDiff_s(  cheequb_obj->n,
                                        cheequb_obj->s,
                                        cheequb_obj->sref );
}

TEST_F(cheequb_test, cheequb1) {
    EXPECT_NEAR(0.0, cheequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cheequb_obj->scond, cheequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cheequb_obj->amax, cheequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheequb_test, cheequb2) {
    EXPECT_NEAR(0.0, cheequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cheequb_obj->scond, cheequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cheequb_obj->amax, cheequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheequb_test, cheequb3) {
    EXPECT_NEAR(0.0, cheequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cheequb_obj->scond, cheequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cheequb_obj->amax, cheequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cheequb_test, cheequb4) {
    EXPECT_NEAR(0.0, cheequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cheequb_obj->scond, cheequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cheequb_obj->amax, cheequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}


/* Begin heequb_dcomplex_parameters  class definition */
class heequb_dcomplex_parameters{
   public:
      int a_bufsize;
      double diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // order of matrix A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_complex_double *a, *aref; // the n-by-n symmetric / Hermitian positive definite matrix
      char uplo; //  Must be 'U' or 'L'.

      /* Output parameters */
      double * s, *sref; // colum scale factors, , array if size 'n'
      double scond, scondref; // ratio of the smallest c[i] to the largest c[i]
      double amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      heequb_dcomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i, char uplo_i);
             
      ~heequb_dcomplex_parameters ();
};  /* end of heequb_dcomplex_parameters  class definition */


/* Constructor heequb_dcomplex_parameters definition */
heequb_dcomplex_parameters:: heequb_dcomplex_parameters ( int matrix_layout_i, 
     lapack_int n_i, char uplo_i) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n heequb dcomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = n*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       heequb_free();
       EXPECT_FALSE( true) << "heequb_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref,
                                                             n, n, 'S');
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

heequb_dcomplex_parameters:: ~heequb_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" heequb_dcomplex_parameters object: destructor invoked. \n");
#endif
   heequb_free();
}


//  Test fixture class definition
class zheequb_test  : public  ::testing::Test {
public:
   heequb_dcomplex_parameters  *zheequb_obj;
   void SetUp();  
   void TearDown () { delete zheequb_obj; }
};


void zheequb_test::SetUp(){

    /* LAPACKE ZHEEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_zheequb) ( int matrix_layout, char uplo, 
                            lapack_int n, const lapack_complex_double *a, 
            lapack_int lda, double *s, double *scond, double *amax );

    Fptr_NL_LAPACKE_zheequb ZHEEQUB;

    zheequb_obj = new  heequb_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,  
                                         lin_solver_paramslist[idx].Uplo);
    idx = Circular_Increment_Index(idx);

    zheequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zheequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zheequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zheequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHEEQUB = (Fptr_NL_LAPACKE_zheequb)dlsym(zheequb_obj->hModule, "LAPACKE_zheequb");
    ASSERT_TRUE(ZHEEQUB != NULL) << "failed to get the Netlib LAPACKE_zheequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zheequb_obj->inforef = ZHEEQUB( zheequb_obj->matrix_layout,
                                  zheequb_obj->uplo,
                                  zheequb_obj->n,
                                  zheequb_obj->aref,
                                  zheequb_obj->lda,
                                  zheequb_obj->sref,
                                  &zheequb_obj->scondref,
                                  &zheequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    zheequb_obj->info = LAPACKE_zheequb( zheequb_obj->matrix_layout,
                                  zheequb_obj->uplo,
                                  zheequb_obj->n,
                                  zheequb_obj->a,
                                  zheequb_obj->lda,
                                  zheequb_obj->s,
                                  &zheequb_obj->scond,
                                  &zheequb_obj->amax);

    if( zheequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zheequb \
        is wrong\n", zheequb_obj->info );
    }
    if( zheequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zheequb is wrong\n",
        zheequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zheequb_obj->diff +=  computeDiff_d(  zheequb_obj->n,
                                        zheequb_obj->s,
                                        zheequb_obj->sref );
}

TEST_F(zheequb_test, zheequb1) {
    EXPECT_NEAR(0.0, zheequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zheequb_obj->scond, zheequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zheequb_obj->amax, zheequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheequb_test, zheequb2) {
    EXPECT_NEAR(0.0, zheequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zheequb_obj->scond, zheequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zheequb_obj->amax, zheequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheequb_test, zheequb3) {
    EXPECT_NEAR(0.0, zheequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zheequb_obj->scond, zheequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zheequb_obj->amax, zheequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zheequb_test, zheequb4) {
    EXPECT_NEAR(0.0, zheequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zheequb_obj->scond, zheequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zheequb_obj->amax, zheequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}
