#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define ppequ_free() \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (s != NULL)    free (s   ); \
       if (sref != NULL) free (sref); \
       if( hModule != NULL) dlclose(hModule); \
       if( dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ppequ_float_parameters  class definition */
class ppequ_float_parameters{
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
      ppequ_float_parameters ( int matrix_layout_i, 
              lapack_int n_i, char uplo_i);
             
      ~ppequ_float_parameters ();
};  /* end of ppequ_float_parameters  class definition */


/* Constructor ppequ_float_parameters definition */
ppequ_float_parameters:: ppequ_float_parameters ( int matrix_layout_i, 
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
   printf(" \n ppequ float:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = ((n+1)*n)/2;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       ppequ_free();
       EXPECT_FALSE( true) << "ppequ_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, a_bufsize);

    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

ppequ_float_parameters:: ~ppequ_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppequ_float_parameters object: destructor invoked. \n");
#endif
   ppequ_free();
}


//  Test fixture class definition
class sppequ_test  : public  ::testing::Test {
public:
   ppequ_float_parameters  *sppequ_obj;
   void SetUp();  
   void TearDown () { delete sppequ_obj; }
};


void sppequ_test::SetUp(){

    /* LAPACKE SPPEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_sppequ) (   int matrix_layout, char uplo,
                    lapack_int n, const float *a, float* s,
                                    float* scond, float* amax  );

    Fptr_NL_LAPACKE_sppequ SPPEQU;

    sppequ_obj = new  ppequ_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n, 
                                         lin_solver_paramslist[idx].Uplo);
    idx = Circular_Increment_Index(idx);

    sppequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sppequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sppequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sppequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SPPEQU = (Fptr_NL_LAPACKE_sppequ)dlsym(sppequ_obj->hModule, "LAPACKE_sppequ");
    ASSERT_TRUE(SPPEQU != NULL) << "failed to get the Netlib LAPACKE_sppequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    sppequ_obj->inforef = SPPEQU( sppequ_obj->matrix_layout,
                                  sppequ_obj->uplo,
                                  sppequ_obj->n,
                                  sppequ_obj->aref,
                                  sppequ_obj->sref,
                                  &sppequ_obj->scondref,
                                  &sppequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    sppequ_obj->info = LAPACKE_sppequ( sppequ_obj->matrix_layout,
                                  sppequ_obj->uplo,
                                  sppequ_obj->n,
                                  sppequ_obj->a,
                                  sppequ_obj->s,
                                  &sppequ_obj->scond,
                                  &sppequ_obj->amax);

    if( sppequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_sppequ \
        is wrong\n", sppequ_obj->info );
    }
    if( sppequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sppequ is wrong\n",
        sppequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    sppequ_obj->diff =  computeDiff_s(  sppequ_obj->n,
                                        sppequ_obj->s,
                                        sppequ_obj->sref );
}

TEST_F(sppequ_test, sppequ1) {
    EXPECT_NEAR(0.0, sppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sppequ_obj->scond, sppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sppequ_obj->amax, sppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sppequ_test, sppequ2) {
    EXPECT_NEAR(0.0, sppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sppequ_obj->scond, sppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sppequ_obj->amax, sppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sppequ_test, sppequ3) {
    EXPECT_NEAR(0.0, sppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sppequ_obj->scond, sppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sppequ_obj->amax, sppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sppequ_test, sppequ4) {
    EXPECT_NEAR(0.0, sppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sppequ_obj->scond, sppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(sppequ_obj->amax, sppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin ppequ_double_parameters  class definition */
class ppequ_double_parameters{
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
      ppequ_double_parameters ( int matrix_layout_i, 
              lapack_int n_i, char uplo_i);
             
      ~ppequ_double_parameters ();
};  /* end of ppequ_double_parameters  class definition */


/* Constructor ppequ_double_parameters definition */
ppequ_double_parameters:: ppequ_double_parameters ( int matrix_layout_i, 
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
   printf(" \n ppequ double:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = ((n+1)*n)/2;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       ppequ_free();
       EXPECT_FALSE( true) << "ppequ_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

ppequ_double_parameters:: ~ppequ_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppequ_double_parameters object: destructor invoked. \n");
#endif
   ppequ_free();
}


//  Test fixture class definition
class dppequ_test  : public  ::testing::Test {
public:
   ppequ_double_parameters  *dppequ_obj;
   void SetUp();  
   void TearDown () { delete dppequ_obj; }
};


void dppequ_test::SetUp(){

    /* LAPACKE DPPEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_dppequ) ( int matrix_layout,char uplo, 
                 lapack_int n, const double *a,
                 double *s, double *scond, double *amax  );

    Fptr_NL_LAPACKE_dppequ DPPEQU;

    dppequ_obj = new  ppequ_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n, 
                                         lin_solver_paramslist[idx].Uplo);
    idx = Circular_Increment_Index(idx);

    dppequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dppequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dppequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dppequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DPPEQU = (Fptr_NL_LAPACKE_dppequ)dlsym(dppequ_obj->hModule, "LAPACKE_dppequ");
    ASSERT_TRUE(DPPEQU != NULL) << "failed to get the Netlib LAPACKE_dppequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    dppequ_obj->inforef = DPPEQU( dppequ_obj->matrix_layout,
                                  dppequ_obj->uplo,
                                  dppequ_obj->n,
                                  dppequ_obj->aref,
                                  dppequ_obj->sref,
                                  &dppequ_obj->scondref,
                                  &dppequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    dppequ_obj->info = LAPACKE_dppequ( dppequ_obj->matrix_layout,
                                  dppequ_obj->uplo,
                                  dppequ_obj->n,
                                  dppequ_obj->a,
                                  dppequ_obj->s,
                                  &dppequ_obj->scond,
                                  &dppequ_obj->amax);

    if( dppequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dppequ \
        is wrong\n", dppequ_obj->info );
    }
    if( dppequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dppequ is wrong\n",
        dppequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dppequ_obj->diff =  computeDiff_d(  dppequ_obj->n,
                                        dppequ_obj->s,
                                        dppequ_obj->sref );
}

TEST_F(dppequ_test, dppequ1) {
    EXPECT_NEAR(0.0, dppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dppequ_obj->scond, dppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dppequ_obj->amax, dppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dppequ_test, dppequ2) {
    EXPECT_NEAR(0.0, dppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dppequ_obj->scond, dppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dppequ_obj->amax, dppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dppequ_test, dppequ3) {
    EXPECT_NEAR(0.0, dppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dppequ_obj->scond, dppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dppequ_obj->amax, dppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dppequ_test, dppequ4) {
    EXPECT_NEAR(0.0, dppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dppequ_obj->scond, dppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dppequ_obj->amax, dppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin ppequ_scomplex_parameters  class definition */
class ppequ_scomplex_parameters{
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
      ppequ_scomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i, char uplo_i);
             
      ~ppequ_scomplex_parameters ();
};  /* end of ppequ_scomplex_parameters  class definition */


/* Constructor ppequ_scomplex_parameters definition */
ppequ_scomplex_parameters:: ppequ_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n ppequ scomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = (n*(n+1))/2;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       ppequ_free();
       EXPECT_FALSE( true) << "ppequ_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

ppequ_scomplex_parameters:: ~ppequ_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppequ_scomplex_parameters object: destructor invoked. \n");
#endif
   ppequ_free();
}


//  Test fixture class definition
class cppequ_test  : public  ::testing::Test {
public:
   ppequ_scomplex_parameters  *cppequ_obj;
   void SetUp();  
   void TearDown () { delete cppequ_obj; }
};


void cppequ_test::SetUp(){

    /* LAPACKE CPPEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_cppequ) ( int matrix_layout, char uplo, 
                            lapack_int n, const lapack_complex_float *a, 
                               float *s, float *scond, float *amax  );

    Fptr_NL_LAPACKE_cppequ CPPEQU;

    cppequ_obj = new  ppequ_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n, 
                                         lin_solver_paramslist[idx].Uplo);
    idx = Circular_Increment_Index(idx);

    cppequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cppequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cppequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cppequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CPPEQU = (Fptr_NL_LAPACKE_cppequ)dlsym(cppequ_obj->hModule, "LAPACKE_cppequ");
    ASSERT_TRUE(CPPEQU != NULL) << "failed to get the Netlib LAPACKE_cppequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    cppequ_obj->inforef = CPPEQU( cppequ_obj->matrix_layout,
                                  cppequ_obj->uplo,
                                  cppequ_obj->n,
                                  cppequ_obj->aref,
                                  cppequ_obj->sref,
                                  &cppequ_obj->scondref,
                                  &cppequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    cppequ_obj->info = LAPACKE_cppequ( cppequ_obj->matrix_layout,
                                  cppequ_obj->uplo,
                                  cppequ_obj->n,
                                  cppequ_obj->a,
                                  cppequ_obj->s,
                                  &cppequ_obj->scond,
                                  &cppequ_obj->amax);

    if( cppequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cppequ \
        is wrong\n", cppequ_obj->info );
    }
    if( cppequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cppequ is wrong\n",
        cppequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cppequ_obj->diff =  computeDiff_s(  cppequ_obj->n,
                                        cppequ_obj->s,
                                        cppequ_obj->sref );
}

TEST_F(cppequ_test, cppequ1) {
    EXPECT_NEAR(0.0, cppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cppequ_obj->scond, cppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cppequ_obj->amax, cppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cppequ_test, cppequ2) {
    EXPECT_NEAR(0.0, cppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cppequ_obj->scond, cppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cppequ_obj->amax, cppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cppequ_test, cppequ3) {
    EXPECT_NEAR(0.0, cppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cppequ_obj->scond, cppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cppequ_obj->amax, cppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cppequ_test, cppequ4) {
    EXPECT_NEAR(0.0, cppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cppequ_obj->scond, cppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cppequ_obj->amax, cppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}


/* Begin ppequ_dcomplex_parameters  class definition */
class ppequ_dcomplex_parameters{
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
      ppequ_dcomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i, char uplo_i);
             
      ~ppequ_dcomplex_parameters ();
};  /* end of ppequ_dcomplex_parameters  class definition */


/* Constructor ppequ_dcomplex_parameters definition */
ppequ_dcomplex_parameters:: ppequ_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n ppequ dcomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = (n*(n+1))/2;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       ppequ_free();
       EXPECT_FALSE( true) << "ppequ_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, a_bufsize);
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

ppequ_dcomplex_parameters:: ~ppequ_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ppequ_dcomplex_parameters object: destructor invoked. \n");
#endif
   ppequ_free();
}


//  Test fixture class definition
class zppequ_test  : public  ::testing::Test {
public:
   ppequ_dcomplex_parameters  *zppequ_obj;
   void SetUp();  
   void TearDown () { delete zppequ_obj; }
};


void zppequ_test::SetUp(){

    /* LAPACKE ZPPEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_zppequ) ( int matrix_layout, char uplo, 
                            lapack_int n, const lapack_complex_double *a, 
                                double *s, double *scond, double *amax );

    Fptr_NL_LAPACKE_zppequ ZPPEQU;

    zppequ_obj = new  ppequ_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,  
                                         lin_solver_paramslist[idx].Uplo);
    idx = Circular_Increment_Index(idx);

    zppequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zppequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zppequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zppequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZPPEQU = (Fptr_NL_LAPACKE_zppequ)dlsym(zppequ_obj->hModule, "LAPACKE_zppequ");
    ASSERT_TRUE(ZPPEQU != NULL) << "failed to get the Netlib LAPACKE_zppequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zppequ_obj->inforef = ZPPEQU( zppequ_obj->matrix_layout,
                                  zppequ_obj->uplo,
                                  zppequ_obj->n,
                                  zppequ_obj->aref,
                                  zppequ_obj->sref,
                                  &zppequ_obj->scondref,
                                  &zppequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    zppequ_obj->info = LAPACKE_zppequ( zppequ_obj->matrix_layout,
                                  zppequ_obj->uplo,
                                  zppequ_obj->n,
                                  zppequ_obj->a,
                                  zppequ_obj->s,
                                  &zppequ_obj->scond,
                                  &zppequ_obj->amax);

    if( zppequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zppequ \
        is wrong\n", zppequ_obj->info );
    }
    if( zppequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zppequ is wrong\n",
        zppequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zppequ_obj->diff +=  computeDiff_d(  zppequ_obj->n,
                                        zppequ_obj->s,
                                        zppequ_obj->sref );
}

TEST_F(zppequ_test, zppequ1) {
    EXPECT_NEAR(0.0, zppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zppequ_obj->scond, zppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zppequ_obj->amax, zppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zppequ_test, zppequ2) {
    EXPECT_NEAR(0.0, zppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zppequ_obj->scond, zppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zppequ_obj->amax, zppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zppequ_test, zppequ3) {
    EXPECT_NEAR(0.0, zppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zppequ_obj->scond, zppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zppequ_obj->amax, zppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zppequ_test, zppequ4) {
    EXPECT_NEAR(0.0, zppequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zppequ_obj->scond, zppequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zppequ_obj->amax, zppequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}
