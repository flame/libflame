#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define poequ_free() \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (s != NULL)    free (s   ); \
       if (sref != NULL) free (sref); \
       if( hModule != NULL) dlclose(hModule); \
       if( dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin poequ_float_parameters  class definition */
class poequ_float_parameters{
   public:
      int a_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // order of matrix A
      float *a, *aref; // the n-by-n symmetric / Hermitian positive definite matrix
      lapack_int lda;  //  leading dimension of 'a'

      /* Output parameters */
      float *s, *sref; // scale factors, of size 'n'
      float scond, scondref; // ratio of the smallest s[i] to the largest s[i]
      float amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      poequ_float_parameters ( int matrix_layout_i, 
              lapack_int n_i);
             
      ~poequ_float_parameters ();
};  /* end of poequ_float_parameters  class definition */


/* Constructor poequ_float_parameters definition */
poequ_float_parameters:: poequ_float_parameters ( int matrix_layout_i, 
      lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n poequ float:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = lda*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       poequ_free();
       EXPECT_FALSE( true) << "poequ_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, n, n, 'S');
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

poequ_float_parameters:: ~poequ_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" poequ_float_parameters object: destructor invoked. \n");
#endif
   poequ_free();
}


//  Test fixture class definition
class spoequ_test  : public  ::testing::Test {
public:
   poequ_float_parameters  *spoequ_obj;
   void SetUp();  
   void TearDown () { delete spoequ_obj; }
};


void spoequ_test::SetUp(){

    /* LAPACKE SPOEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_spoequ) ( int matrix_layout, lapack_int n,
                                    const float* a, lapack_int lda, float* s,
                                    float* scond, float* amax  );

    Fptr_NL_LAPACKE_spoequ SPOEQU;

    spoequ_obj = new  poequ_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    spoequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spoequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spoequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spoequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SPOEQU = (Fptr_NL_LAPACKE_spoequ)dlsym(spoequ_obj->hModule, "LAPACKE_spoequ");
    ASSERT_TRUE(SPOEQU != NULL) << "failed to get the Netlib LAPACKE_spoequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    spoequ_obj->inforef = SPOEQU( spoequ_obj->matrix_layout,
                                  spoequ_obj->n,
                                  spoequ_obj->aref,
                                  spoequ_obj->lda,
                                  spoequ_obj->sref,
                                  &spoequ_obj->scondref,
                                  &spoequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    spoequ_obj->info = LAPACKE_spoequ( spoequ_obj->matrix_layout,
                                  spoequ_obj->n,
                                  spoequ_obj->a,
                                  spoequ_obj->lda,
                                  spoequ_obj->s,
                                  &spoequ_obj->scond,
                                  &spoequ_obj->amax);

    if( spoequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spoequ \
        is wrong\n", spoequ_obj->info );
    }
    if( spoequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spoequ is wrong\n",
        spoequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    spoequ_obj->diff =  computeDiff_s(  spoequ_obj->n,
                                        spoequ_obj->s,
                                        spoequ_obj->sref );
}

TEST_F(spoequ_test, spoequ1) {
    EXPECT_NEAR(0.0, spoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequ_obj->scond, spoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequ_obj->amax, spoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spoequ_test, spoequ2) {
    EXPECT_NEAR(0.0, spoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequ_obj->scond, spoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequ_obj->amax, spoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spoequ_test, spoequ3) {
    EXPECT_NEAR(0.0, spoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequ_obj->scond, spoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequ_obj->amax, spoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spoequ_test, spoequ4) {
    EXPECT_NEAR(0.0, spoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequ_obj->scond, spoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequ_obj->amax, spoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin poequ_double_parameters  class definition */
class poequ_double_parameters{
   public:
      int a_bufsize;
      double diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // order of matrix A
      lapack_int lda;  //  leading dimension of 'a'
      double *a, *aref; // the n-by-n symmetric / Hermitian positive definite matrix

      /* Output parameters */
      double * s, *sref; // colum scale factors, , array if size 'n'
      double scond, scondref; // ratio of the smallest c[i] to the largest c[i]
      double amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      poequ_double_parameters ( int matrix_layout_i, 
              lapack_int n_i);
             
      ~poequ_double_parameters ();
};  /* end of poequ_double_parameters  class definition */


/* Constructor poequ_double_parameters definition */
poequ_double_parameters:: poequ_double_parameters ( int matrix_layout_i, 
      lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n poequ double:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = lda*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       poequ_free();
       EXPECT_FALSE( true) << "poequ_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, n, n, 'S');
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

poequ_double_parameters:: ~poequ_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" poequ_double_parameters object: destructor invoked. \n");
#endif
   poequ_free();
}


//  Test fixture class definition
class dpoequ_test  : public  ::testing::Test {
public:
   poequ_double_parameters  *dpoequ_obj;
   void SetUp();  
   void TearDown () { delete dpoequ_obj; }
};


void dpoequ_test::SetUp(){

    /* LAPACKE DPOEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_dpoequ) ( int matrix_layout, 
                 lapack_int n, const double *a, lapack_int lda,
                 double *s, double *scond, double *amax  );

    Fptr_NL_LAPACKE_dpoequ DPOEQU;

    dpoequ_obj = new  poequ_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    dpoequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpoequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpoequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpoequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DPOEQU = (Fptr_NL_LAPACKE_dpoequ)dlsym(dpoequ_obj->hModule, "LAPACKE_dpoequ");
    ASSERT_TRUE(DPOEQU != NULL) << "failed to get the Netlib LAPACKE_dpoequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    dpoequ_obj->inforef = DPOEQU( dpoequ_obj->matrix_layout,
                                  dpoequ_obj->n,
                                  dpoequ_obj->aref,
                                  dpoequ_obj->lda,
                                  dpoequ_obj->sref,
                                  &dpoequ_obj->scondref,
                                  &dpoequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    dpoequ_obj->info = LAPACKE_dpoequ( dpoequ_obj->matrix_layout,
                                  dpoequ_obj->n,
                                  dpoequ_obj->a,
                                  dpoequ_obj->lda,
                                  dpoequ_obj->s,
                                  &dpoequ_obj->scond,
                                  &dpoequ_obj->amax);

    if( dpoequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpoequ \
        is wrong\n", dpoequ_obj->info );
    }
    if( dpoequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpoequ is wrong\n",
        dpoequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dpoequ_obj->diff =  computeDiff_d(  dpoequ_obj->n,
                                        dpoequ_obj->s,
                                        dpoequ_obj->sref );
}

TEST_F(dpoequ_test, dpoequ1) {
    EXPECT_NEAR(0.0, dpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequ_obj->scond, dpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequ_obj->amax, dpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpoequ_test, dpoequ2) {
    EXPECT_NEAR(0.0, dpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequ_obj->scond, dpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequ_obj->amax, dpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpoequ_test, dpoequ3) {
    EXPECT_NEAR(0.0, dpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequ_obj->scond, dpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequ_obj->amax, dpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpoequ_test, dpoequ4) {
    EXPECT_NEAR(0.0, dpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequ_obj->scond, dpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequ_obj->amax, dpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin poequ_scomplex_parameters  class definition */
class poequ_scomplex_parameters{
   public:
      int a_bufsize;
      float diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // order of matrix A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_complex_float *a, *aref; // the n-by-n symmetric / Hermitian positive definite matrix

      /* Output parameters */
      float * s, *sref; // colum scale factors, , array if size 'n'
      float scond, scondref; // ratio of the smallest c[i] to the largest c[i]
      float amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      poequ_scomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i);
             
      ~poequ_scomplex_parameters ();
};  /* end of poequ_scomplex_parameters  class definition */


/* Constructor poequ_scomplex_parameters definition */
poequ_scomplex_parameters:: poequ_scomplex_parameters ( int matrix_layout_i, 
      lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n poequ scomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = lda*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       poequ_free();
       EXPECT_FALSE( true) << "poequ_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, n, n,
                                            lin_solver_paramslist[idx].symm);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

poequ_scomplex_parameters:: ~poequ_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" poequ_scomplex_parameters object: destructor invoked. \n");
#endif
   poequ_free();
}


//  Test fixture class definition
class cpoequ_test  : public  ::testing::Test {
public:
   poequ_scomplex_parameters  *cpoequ_obj;
   void SetUp();  
   void TearDown () { delete cpoequ_obj; }
};


void cpoequ_test::SetUp(){

    /* LAPACKE CPOEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_cpoequ) ( int matrix_layout, lapack_int n,
                               const lapack_complex_float *a, 
                               lapack_int lda, float *s, 
                               float *scond, float *amax  );

    Fptr_NL_LAPACKE_cpoequ CPOEQU;

    cpoequ_obj = new  poequ_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    cpoequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpoequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpoequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpoequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CPOEQU = (Fptr_NL_LAPACKE_cpoequ)dlsym(cpoequ_obj->hModule, "LAPACKE_cpoequ");
    ASSERT_TRUE(CPOEQU != NULL) << "failed to get the Netlib LAPACKE_cpoequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    cpoequ_obj->inforef = CPOEQU( cpoequ_obj->matrix_layout,
                                  cpoequ_obj->n,
                                  cpoequ_obj->aref,
                                  cpoequ_obj->lda,
                                  cpoequ_obj->sref,
                                  &cpoequ_obj->scondref,
                                  &cpoequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    cpoequ_obj->info = LAPACKE_cpoequ( cpoequ_obj->matrix_layout,
                                  cpoequ_obj->n,
                                  cpoequ_obj->a,
                                  cpoequ_obj->lda,
                                  cpoequ_obj->s,
                                  &cpoequ_obj->scond,
                                  &cpoequ_obj->amax);

    if( cpoequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpoequ \
        is wrong\n", cpoequ_obj->info );
    }
    if( cpoequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpoequ is wrong\n",
        cpoequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cpoequ_obj->diff =  computeDiff_s(  cpoequ_obj->n,
                                        cpoequ_obj->s,
                                        cpoequ_obj->sref );
}

TEST_F(cpoequ_test, cpoequ1) {
    EXPECT_NEAR(0.0, cpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequ_obj->scond, cpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequ_obj->amax, cpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpoequ_test, cpoequ2) {
    EXPECT_NEAR(0.0, cpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequ_obj->scond, cpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequ_obj->amax, cpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpoequ_test, cpoequ3) {
    EXPECT_NEAR(0.0, cpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequ_obj->scond, cpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequ_obj->amax, cpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpoequ_test, cpoequ4) {
    EXPECT_NEAR(0.0, cpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequ_obj->scond, cpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequ_obj->amax, cpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}


/* Begin poequ_dcomplex_parameters  class definition */
class poequ_dcomplex_parameters{
   public:
      int a_bufsize;
      double diff; // to capture the netlib, libflame o/ps
      void *hModule, *dModule;
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // order of matrix A
      lapack_int lda;  //  leading dimension of 'a'
      lapack_complex_double *a, *aref; // the n-by-n symmetric / Hermitian positive definite matrix

      /* Output parameters */
      double * s, *sref; // colum scale factors, , array if size 'n'
      double scond, scondref; // ratio of the smallest c[i] to the largest c[i]
      double amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      poequ_dcomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i);
             
      ~poequ_dcomplex_parameters ();
};  /* end of poequ_dcomplex_parameters  class definition */


/* Constructor poequ_dcomplex_parameters definition */
poequ_dcomplex_parameters:: poequ_dcomplex_parameters ( int matrix_layout_i, 
      lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n poequ dcomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = lda*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       poequ_free();
       EXPECT_FALSE( true) << "poequ_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, n, n,
                                            lin_solver_paramslist[idx].symm);
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

poequ_dcomplex_parameters:: ~poequ_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" poequ_dcomplex_parameters object: destructor invoked. \n");
#endif
   poequ_free();
}


//  Test fixture class definition
class zpoequ_test  : public  ::testing::Test {
public:
   poequ_dcomplex_parameters  *zpoequ_obj;
   void SetUp();  
   void TearDown () { delete zpoequ_obj; }
};


void zpoequ_test::SetUp(){

    /* LAPACKE ZPOEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_zpoequ) ( int matrix_layout, lapack_int n,
                               const lapack_complex_double *a, 
                               lapack_int lda, double *s, 
                               double *scond, double *amax  );

    Fptr_NL_LAPACKE_zpoequ ZPOEQU;

    zpoequ_obj = new  poequ_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    zpoequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpoequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpoequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpoequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZPOEQU = (Fptr_NL_LAPACKE_zpoequ)dlsym(zpoequ_obj->hModule, "LAPACKE_zpoequ");
    ASSERT_TRUE(ZPOEQU != NULL) << "failed to get the Netlib LAPACKE_zpoequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zpoequ_obj->inforef = ZPOEQU( zpoequ_obj->matrix_layout,
                                  zpoequ_obj->n,
                                  zpoequ_obj->aref,
                                  zpoequ_obj->lda,
                                  zpoequ_obj->sref,
                                  &zpoequ_obj->scondref,
                                  &zpoequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    zpoequ_obj->info = LAPACKE_zpoequ( zpoequ_obj->matrix_layout,
                                  zpoequ_obj->n,
                                  zpoequ_obj->a,
                                  zpoequ_obj->lda,
                                  zpoequ_obj->s,
                                  &zpoequ_obj->scond,
                                  &zpoequ_obj->amax);

    if( zpoequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpoequ \
        is wrong\n", zpoequ_obj->info );
    }
    if( zpoequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpoequ is wrong\n",
        zpoequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zpoequ_obj->diff +=  computeDiff_d(  zpoequ_obj->n,
                                        zpoequ_obj->s,
                                        zpoequ_obj->sref );
}

TEST_F(zpoequ_test, zpoequ1) {
    EXPECT_NEAR(0.0, zpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequ_obj->scond, zpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequ_obj->amax, zpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpoequ_test, zpoequ2) {
    EXPECT_NEAR(0.0, zpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequ_obj->scond, zpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequ_obj->amax, zpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpoequ_test, zpoequ3) {
    EXPECT_NEAR(0.0, zpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequ_obj->scond, zpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequ_obj->amax, zpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpoequ_test, zpoequ4) {
    EXPECT_NEAR(0.0, zpoequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequ_obj->scond, zpoequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequ_obj->amax, zpoequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}
