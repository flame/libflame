#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define syequb_free() \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (s != NULL)    free (s   ); \
       if (sref != NULL) free (sref); \
       if( hModule != NULL) dlclose(hModule); \
       if( dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin syequb_float_parameters  class definition */
class syequb_float_parameters{
   public:
      int a_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // order of matrix A
      float *a, *aref; // the n-by-n symmetric / Hermitian positive definite matrix
      lapack_int lda;  //  leading dimension of 'a'
      char uplo; //  Must be 'U' or 'L'.


      /* Output parameters */
      float *s, *sref; // scale factors, of size 'n'
      float scond, scondref; // ratio of the smallest s[i] to the largest s[i]
      float amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      syequb_float_parameters ( int matrix_layout_i, 
              lapack_int n_i, char uplo_i);
             
      ~syequb_float_parameters ();
};  /* end of syequb_float_parameters  class definition */


/* Constructor syequb_float_parameters definition */
syequb_float_parameters:: syequb_float_parameters ( int matrix_layout_i, 
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
   printf(" \n syequb float:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = n*n;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       syequb_free();
       EXPECT_FALSE( true) << "syequb_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref,
                                                             n, n, 'S');

    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

syequb_float_parameters:: ~syequb_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" syequb_float_parameters object: destructor invoked. \n");
#endif
   syequb_free();
}


//  Test fixture class definition
class ssyequb_test  : public  ::testing::Test {
public:
   syequb_float_parameters  *ssyequb_obj;
   void SetUp();  
   void TearDown () { delete ssyequb_obj; }
};


void ssyequb_test::SetUp(){

    /* LAPACKE SSYEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_ssyequb) (   int matrix_layout, char uplo,
                    lapack_int n, const float *a, lapack_int lda, float* s,
                                    float* scond, float* amax  );

    Fptr_NL_LAPACKE_ssyequb SSYEQUB;

    ssyequb_obj = new  syequb_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n, 
                                         lin_solver_paramslist[idx].Uplo);
    idx = Circular_Increment_Index(idx);

    ssyequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssyequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssyequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssyequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYEQUB = (Fptr_NL_LAPACKE_ssyequb)dlsym(ssyequb_obj->hModule, "LAPACKE_ssyequb");
    ASSERT_TRUE(SSYEQUB != NULL) << "failed to get the Netlib LAPACKE_ssyequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    ssyequb_obj->inforef = SSYEQUB( ssyequb_obj->matrix_layout,
                                  ssyequb_obj->uplo,
                                  ssyequb_obj->n,
                                  ssyequb_obj->aref,
                                  ssyequb_obj->lda,
                                  ssyequb_obj->sref,
                                  &ssyequb_obj->scondref,
                                  &ssyequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    ssyequb_obj->info = LAPACKE_ssyequb( ssyequb_obj->matrix_layout,
                                  ssyequb_obj->uplo,
                                  ssyequb_obj->n,
                                  ssyequb_obj->a,
                                  ssyequb_obj->lda,
                                  ssyequb_obj->s,
                                  &ssyequb_obj->scond,
                                  &ssyequb_obj->amax);

    if( ssyequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssyequb \
        is wrong\n", ssyequb_obj->info );
    }
    if( ssyequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssyequb is wrong\n",
        ssyequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ssyequb_obj->diff =  computeDiff_s(  ssyequb_obj->n,
                                        ssyequb_obj->s,
                                        ssyequb_obj->sref );
}

TEST_F(ssyequb_test, ssyequb1) {
    EXPECT_NEAR(0.0, ssyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(ssyequb_obj->scond, ssyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(ssyequb_obj->amax, ssyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssyequb_test, ssyequb2) {
    EXPECT_NEAR(0.0, ssyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(ssyequb_obj->scond, ssyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(ssyequb_obj->amax, ssyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssyequb_test, ssyequb3) {
    EXPECT_NEAR(0.0, ssyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(ssyequb_obj->scond, ssyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(ssyequb_obj->amax, ssyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ssyequb_test, ssyequb4) {
    EXPECT_NEAR(0.0, ssyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(ssyequb_obj->scond, ssyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(ssyequb_obj->amax, ssyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin syequb_double_parameters  class definition */
class syequb_double_parameters{
   public:
      int a_bufsize;
      double diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // order of matrix A
      lapack_int lda;  //  leading dimension of 'a'
      double *a, *aref; // the n-by-n symmetric / Hermitian positive definite matrix
      char uplo; //  Must be 'U' or 'L'.

      /* Output parameters */
      double * s, *sref; // colum scale factors, , array if size 'n'
      double scond, scondref; // ratio of the smallest c[i] to the largest c[i]
      double amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      syequb_double_parameters ( int matrix_layout_i, 
              lapack_int n_i, char uplo_i);
             
      ~syequb_double_parameters ();
};  /* end of syequb_double_parameters  class definition */


/* Constructor syequb_double_parameters definition */
syequb_double_parameters:: syequb_double_parameters ( int matrix_layout_i, 
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
   printf(" \n syequb double:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize =n*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       syequb_free();
       EXPECT_FALSE( true) << "syequb_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref,
                                                             n, n, 'S');
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

syequb_double_parameters:: ~syequb_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" syequb_double_parameters object: destructor invoked. \n");
#endif
   syequb_free();
}


//  Test fixture class definition
class dsyequb_test  : public  ::testing::Test {
public:
   syequb_double_parameters  *dsyequb_obj;
   void SetUp();  
   void TearDown () { delete dsyequb_obj; }
};


void dsyequb_test::SetUp(){

    /* LAPACKE DSYEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_dsyequb) ( int matrix_layout,char uplo, 
                 lapack_int n, const double *a,  lapack_int lda,
                 double *s, double *scond, double *amax  );

    Fptr_NL_LAPACKE_dsyequb DSYEQUB;

    dsyequb_obj = new  syequb_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n, 
                                         lin_solver_paramslist[idx].Uplo);
    idx = Circular_Increment_Index(idx);

    dsyequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsyequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsyequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsyequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYEQUB = (Fptr_NL_LAPACKE_dsyequb)dlsym(dsyequb_obj->hModule, "LAPACKE_dsyequb");
    ASSERT_TRUE(DSYEQUB != NULL) << "failed to get the Netlib LAPACKE_dsyequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    dsyequb_obj->inforef = DSYEQUB( dsyequb_obj->matrix_layout,
                                  dsyequb_obj->uplo,
                                  dsyequb_obj->n,
                                  dsyequb_obj->aref,
                                  dsyequb_obj->lda,
                                  dsyequb_obj->sref,
                                  &dsyequb_obj->scondref,
                                  &dsyequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    dsyequb_obj->info = LAPACKE_dsyequb( dsyequb_obj->matrix_layout,
                                  dsyequb_obj->uplo,
                                  dsyequb_obj->n,
                                  dsyequb_obj->a,
                                  dsyequb_obj->lda,
                                  dsyequb_obj->s,
                                  &dsyequb_obj->scond,
                                  &dsyequb_obj->amax);

    if( dsyequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsyequb \
        is wrong\n", dsyequb_obj->info );
    }
    if( dsyequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsyequb is wrong\n",
        dsyequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dsyequb_obj->diff =  computeDiff_d(  dsyequb_obj->n,
                                        dsyequb_obj->s,
                                        dsyequb_obj->sref );
}

TEST_F(dsyequb_test, dsyequb1) {
    EXPECT_NEAR(0.0, dsyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dsyequb_obj->scond, dsyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dsyequb_obj->amax, dsyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsyequb_test, dsyequb2) {
    EXPECT_NEAR(0.0, dsyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dsyequb_obj->scond, dsyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dsyequb_obj->amax, dsyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsyequb_test, dsyequb3) {
    EXPECT_NEAR(0.0, dsyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dsyequb_obj->scond, dsyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dsyequb_obj->amax, dsyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dsyequb_test, dsyequb4) {
    EXPECT_NEAR(0.0, dsyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dsyequb_obj->scond, dsyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dsyequb_obj->amax, dsyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin syequb_scomplex_parameters  class definition */
class syequb_scomplex_parameters{
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
      syequb_scomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i, char uplo_i);
             
      ~syequb_scomplex_parameters ();
};  /* end of syequb_scomplex_parameters  class definition */


/* Constructor syequb_scomplex_parameters definition */
syequb_scomplex_parameters:: syequb_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n syequb scomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = n*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       syequb_free();
       EXPECT_FALSE( true) << "syequb_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref,
                                                             n, n, 'S');
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

syequb_scomplex_parameters:: ~syequb_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" syequb_scomplex_parameters object: destructor invoked. \n");
#endif
   syequb_free();
}


//  Test fixture class definition
class csyequb_test  : public  ::testing::Test {
public:
   syequb_scomplex_parameters  *csyequb_obj;
   void SetUp();  
   void TearDown () { delete csyequb_obj; }
};


void csyequb_test::SetUp(){

    /* LAPACKE CSYEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_csyequb) ( int matrix_layout, char uplo, 
                            lapack_int n, const lapack_complex_float *a, 
                 lapack_int lda, float *s, float *scond, float *amax  );

    Fptr_NL_LAPACKE_csyequb CSYEQUB;

    csyequb_obj = new  syequb_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n, 
                                         lin_solver_paramslist[idx].Uplo);
    idx = Circular_Increment_Index(idx);

    csyequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csyequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csyequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csyequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYEQUB = (Fptr_NL_LAPACKE_csyequb)dlsym(csyequb_obj->hModule, "LAPACKE_csyequb");
    ASSERT_TRUE(CSYEQUB != NULL) << "failed to get the Netlib LAPACKE_csyequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    csyequb_obj->inforef = CSYEQUB( csyequb_obj->matrix_layout,
                                  csyequb_obj->uplo,
                                  csyequb_obj->n,
                                  csyequb_obj->aref,
                                  csyequb_obj->lda,
                                  csyequb_obj->sref,
                                  &csyequb_obj->scondref,
                                  &csyequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    csyequb_obj->info = LAPACKE_csyequb( csyequb_obj->matrix_layout,
                                  csyequb_obj->uplo,
                                  csyequb_obj->n,
                                  csyequb_obj->a,
                                  csyequb_obj->lda,
                                  csyequb_obj->s,
                                  &csyequb_obj->scond,
                                  &csyequb_obj->amax);

    if( csyequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csyequb \
        is wrong\n", csyequb_obj->info );
    }
    if( csyequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csyequb is wrong\n",
        csyequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    csyequb_obj->diff =  computeDiff_s(  csyequb_obj->n,
                                        csyequb_obj->s,
                                        csyequb_obj->sref );
}

TEST_F(csyequb_test, csyequb1) {
    EXPECT_NEAR(0.0, csyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(csyequb_obj->scond, csyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(csyequb_obj->amax, csyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyequb_test, csyequb2) {
    EXPECT_NEAR(0.0, csyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(csyequb_obj->scond, csyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(csyequb_obj->amax, csyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyequb_test, csyequb3) {
    EXPECT_NEAR(0.0, csyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(csyequb_obj->scond, csyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(csyequb_obj->amax, csyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyequb_test, csyequb4) {
    EXPECT_NEAR(0.0, csyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(csyequb_obj->scond, csyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(csyequb_obj->amax, csyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}


/* Begin syequb_dcomplex_parameters  class definition */
class syequb_dcomplex_parameters{
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
      syequb_dcomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i, char uplo_i);
             
      ~syequb_dcomplex_parameters ();
};  /* end of syequb_dcomplex_parameters  class definition */


/* Constructor syequb_dcomplex_parameters definition */
syequb_dcomplex_parameters:: syequb_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n syequb dcomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = n*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       syequb_free();
       EXPECT_FALSE( true) << "syequb_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref,
                                                             n, n, 'S');
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

syequb_dcomplex_parameters:: ~syequb_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" syequb_dcomplex_parameters object: destructor invoked. \n");
#endif
   syequb_free();
}


//  Test fixture class definition
class zsyequb_test  : public  ::testing::Test {
public:
   syequb_dcomplex_parameters  *zsyequb_obj;
   void SetUp();  
   void TearDown () { delete zsyequb_obj; }
};


void zsyequb_test::SetUp(){

    /* LAPACKE ZSYEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_zsyequb) ( int matrix_layout, char uplo, 
                            lapack_int n, const lapack_complex_double *a, 
            lapack_int lda, double *s, double *scond, double *amax );

    Fptr_NL_LAPACKE_zsyequb ZSYEQUB;

    zsyequb_obj = new  syequb_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,  
                                         lin_solver_paramslist[idx].Uplo);
    idx = Circular_Increment_Index(idx);

    zsyequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsyequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsyequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsyequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYEQUB = (Fptr_NL_LAPACKE_zsyequb)dlsym(zsyequb_obj->hModule, "LAPACKE_zsyequb");
    ASSERT_TRUE(ZSYEQUB != NULL) << "failed to get the Netlib LAPACKE_zsyequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zsyequb_obj->inforef = ZSYEQUB( zsyequb_obj->matrix_layout,
                                  zsyequb_obj->uplo,
                                  zsyequb_obj->n,
                                  zsyequb_obj->aref,
                                  zsyequb_obj->lda,
                                  zsyequb_obj->sref,
                                  &zsyequb_obj->scondref,
                                  &zsyequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    zsyequb_obj->info = LAPACKE_zsyequb( zsyequb_obj->matrix_layout,
                                  zsyequb_obj->uplo,
                                  zsyequb_obj->n,
                                  zsyequb_obj->a,
                                  zsyequb_obj->lda,
                                  zsyequb_obj->s,
                                  &zsyequb_obj->scond,
                                  &zsyequb_obj->amax);

    if( zsyequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsyequb \
        is wrong\n", zsyequb_obj->info );
    }
    if( zsyequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsyequb is wrong\n",
        zsyequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zsyequb_obj->diff +=  computeDiff_d(  zsyequb_obj->n,
                                        zsyequb_obj->s,
                                        zsyequb_obj->sref );
}

TEST_F(zsyequb_test, zsyequb1) {
    EXPECT_NEAR(0.0, zsyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zsyequb_obj->scond, zsyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zsyequb_obj->amax, zsyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyequb_test, zsyequb2) {
    EXPECT_NEAR(0.0, zsyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zsyequb_obj->scond, zsyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zsyequb_obj->amax, zsyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyequb_test, zsyequb3) {
    EXPECT_NEAR(0.0, zsyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zsyequb_obj->scond, zsyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zsyequb_obj->amax, zsyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyequb_test, zsyequb4) {
    EXPECT_NEAR(0.0, zsyequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zsyequb_obj->scond, zsyequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zsyequb_obj->amax, zsyequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}
