#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define poequb_free() \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (s != NULL)    free (s   ); \
       if (sref != NULL) free (sref); \
       if( hModule != NULL) dlclose(hModule); \
       if( dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin poequb_float_parameters  class definition */
class poequb_float_parameters{
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
      poequb_float_parameters ( int matrix_layout_i, 
              lapack_int n_i);
             
      ~poequb_float_parameters ();
};  /* end of poequb_float_parameters  class definition */


/* Constructor poequb_float_parameters definition */
poequb_float_parameters:: poequb_float_parameters ( int matrix_layout_i, 
      lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n poequb float:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = lda*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       poequb_free();
       EXPECT_FALSE( true) << "poequb_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, n, n, 'S');
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

poequb_float_parameters:: ~poequb_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" poequb_float_parameters object: destructor invoked. \n");
#endif
   poequb_free();
}


//  Test fixture class definition
class spoequb_test  : public  ::testing::Test {
public:
   poequb_float_parameters  *spoequb_obj;
   void SetUp();  
   void TearDown () { delete spoequb_obj; }
};


void spoequb_test::SetUp(){

    /* LAPACKE SPOEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_spoequb) ( int matrix_layout, lapack_int n,
                                    const float* a, lapack_int lda, float* s,
                                    float* scond, float* amax  );

    Fptr_NL_LAPACKE_spoequb SPOEQUB;

    spoequb_obj = new  poequb_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    spoequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spoequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spoequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spoequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SPOEQUB = (Fptr_NL_LAPACKE_spoequb)dlsym(spoequb_obj->hModule, "LAPACKE_spoequb");
    ASSERT_TRUE(SPOEQUB != NULL) << "failed to get the Netlib LAPACKE_spoequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    spoequb_obj->inforef = SPOEQUB( spoequb_obj->matrix_layout,
                                  spoequb_obj->n,
                                  spoequb_obj->aref,
                                  spoequb_obj->lda,
                                  spoequb_obj->sref,
                                  &spoequb_obj->scondref,
                                  &spoequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    spoequb_obj->info = LAPACKE_spoequb( spoequb_obj->matrix_layout,
                                  spoequb_obj->n,
                                  spoequb_obj->a,
                                  spoequb_obj->lda,
                                  spoequb_obj->s,
                                  &spoequb_obj->scond,
                                  &spoequb_obj->amax);

    if( spoequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spoequb \
        is wrong\n", spoequb_obj->info );
    }
    if( spoequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spoequb is wrong\n",
        spoequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    spoequb_obj->diff =  computeDiff_s(  spoequb_obj->n,
                                        spoequb_obj->s,
                                        spoequb_obj->sref );
}

TEST_F(spoequb_test, spoequb1) {
    EXPECT_NEAR(0.0, spoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequb_obj->scond, spoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequb_obj->amax, spoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spoequb_test, spoequb2) {
    EXPECT_NEAR(0.0, spoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequb_obj->scond, spoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequb_obj->amax, spoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spoequb_test, spoequb3) {
    EXPECT_NEAR(0.0, spoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequb_obj->scond, spoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequb_obj->amax, spoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spoequb_test, spoequb4) {
    EXPECT_NEAR(0.0, spoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequb_obj->scond, spoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spoequb_obj->amax, spoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin poequb_double_parameters  class definition */
class poequb_double_parameters{
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
      poequb_double_parameters ( int matrix_layout_i, 
              lapack_int n_i);
             
      ~poequb_double_parameters ();
};  /* end of poequb_double_parameters  class definition */


/* Constructor poequb_double_parameters definition */
poequb_double_parameters:: poequb_double_parameters ( int matrix_layout_i, 
      lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n poequb double:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = lda*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       poequb_free();
       EXPECT_FALSE( true) << "poequb_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, n, n, 'S');
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

poequb_double_parameters:: ~poequb_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" poequb_double_parameters object: destructor invoked. \n");
#endif
   poequb_free();
}


//  Test fixture class definition
class dpoequb_test  : public  ::testing::Test {
public:
   poequb_double_parameters  *dpoequb_obj;
   void SetUp();  
   void TearDown () { delete dpoequb_obj; }
};


void dpoequb_test::SetUp(){

    /* LAPACKE DPOEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_dpoequb) ( int matrix_layout, 
                 lapack_int n, const double *a, lapack_int lda,
                 double *s, double *scond, double *amax  );

    Fptr_NL_LAPACKE_dpoequb DPOEQUB;

    dpoequb_obj = new  poequb_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    dpoequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpoequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpoequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpoequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DPOEQUB = (Fptr_NL_LAPACKE_dpoequb)dlsym(dpoequb_obj->hModule, "LAPACKE_dpoequb");
    ASSERT_TRUE(DPOEQUB != NULL) << "failed to get the Netlib LAPACKE_dpoequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    dpoequb_obj->inforef = DPOEQUB( dpoequb_obj->matrix_layout,
                                  dpoequb_obj->n,
                                  dpoequb_obj->aref,
                                  dpoequb_obj->lda,
                                  dpoequb_obj->sref,
                                  &dpoequb_obj->scondref,
                                  &dpoequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    dpoequb_obj->info = LAPACKE_dpoequb( dpoequb_obj->matrix_layout,
                                  dpoequb_obj->n,
                                  dpoequb_obj->a,
                                  dpoequb_obj->lda,
                                  dpoequb_obj->s,
                                  &dpoequb_obj->scond,
                                  &dpoequb_obj->amax);

    if( dpoequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpoequb \
        is wrong\n", dpoequb_obj->info );
    }
    if( dpoequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpoequb is wrong\n",
        dpoequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dpoequb_obj->diff =  computeDiff_d(  dpoequb_obj->n,
                                        dpoequb_obj->s,
                                        dpoequb_obj->sref );
}

TEST_F(dpoequb_test, dpoequb1) {
    EXPECT_NEAR(0.0, dpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequb_obj->scond, dpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequb_obj->amax, dpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpoequb_test, dpoequb2) {
    EXPECT_NEAR(0.0, dpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequb_obj->scond, dpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequb_obj->amax, dpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpoequb_test, dpoequb3) {
    EXPECT_NEAR(0.0, dpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequb_obj->scond, dpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequb_obj->amax, dpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpoequb_test, dpoequb4) {
    EXPECT_NEAR(0.0, dpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequb_obj->scond, dpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpoequb_obj->amax, dpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin poequb_scomplex_parameters  class definition */
class poequb_scomplex_parameters{
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
      poequb_scomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i);
             
      ~poequb_scomplex_parameters ();
};  /* end of poequb_scomplex_parameters  class definition */


/* Constructor poequb_scomplex_parameters definition */
poequb_scomplex_parameters:: poequb_scomplex_parameters ( int matrix_layout_i, 
      lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n poequb scomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = lda*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       poequb_free();
       EXPECT_FALSE( true) << "poequb_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, n, n,
                                            lin_solver_paramslist[idx].symm);
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

poequb_scomplex_parameters:: ~poequb_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" poequb_scomplex_parameters object: destructor invoked. \n");
#endif
   poequb_free();
}


//  Test fixture class definition
class cpoequb_test  : public  ::testing::Test {
public:
   poequb_scomplex_parameters  *cpoequb_obj;
   void SetUp();  
   void TearDown () { delete cpoequb_obj; }
};


void cpoequb_test::SetUp(){

    /* LAPACKE CPOEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_cpoequb) ( int matrix_layout, lapack_int n,
                               const lapack_complex_float *a, 
                               lapack_int lda, float *s, 
                               float *scond, float *amax  );

    Fptr_NL_LAPACKE_cpoequb CPOEQUB;

    cpoequb_obj = new  poequb_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    cpoequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpoequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpoequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpoequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CPOEQUB = (Fptr_NL_LAPACKE_cpoequb)dlsym(cpoequb_obj->hModule, "LAPACKE_cpoequb");
    ASSERT_TRUE(CPOEQUB != NULL) << "failed to get the Netlib LAPACKE_cpoequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    cpoequb_obj->inforef = CPOEQUB( cpoequb_obj->matrix_layout,
                                  cpoequb_obj->n,
                                  cpoequb_obj->aref,
                                  cpoequb_obj->lda,
                                  cpoequb_obj->sref,
                                  &cpoequb_obj->scondref,
                                  &cpoequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    cpoequb_obj->info = LAPACKE_cpoequb( cpoequb_obj->matrix_layout,
                                  cpoequb_obj->n,
                                  cpoequb_obj->a,
                                  cpoequb_obj->lda,
                                  cpoequb_obj->s,
                                  &cpoequb_obj->scond,
                                  &cpoequb_obj->amax);

    if( cpoequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpoequb \
        is wrong\n", cpoequb_obj->info );
    }
    if( cpoequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpoequb is wrong\n",
        cpoequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cpoequb_obj->diff =  computeDiff_s(  cpoequb_obj->n,
                                        cpoequb_obj->s,
                                        cpoequb_obj->sref );
}

TEST_F(cpoequb_test, cpoequb1) {
    EXPECT_NEAR(0.0, cpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequb_obj->scond, cpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequb_obj->amax, cpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpoequb_test, cpoequb2) {
    EXPECT_NEAR(0.0, cpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequb_obj->scond, cpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequb_obj->amax, cpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpoequb_test, cpoequb3) {
    EXPECT_NEAR(0.0, cpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequb_obj->scond, cpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequb_obj->amax, cpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpoequb_test, cpoequb4) {
    EXPECT_NEAR(0.0, cpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequb_obj->scond, cpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpoequb_obj->amax, cpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}


/* Begin poequb_dcomplex_parameters  class definition */
class poequb_dcomplex_parameters{
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
      poequb_dcomplex_parameters ( int matrix_layout_i, 
              lapack_int n_i);
             
      ~poequb_dcomplex_parameters ();
};  /* end of poequb_dcomplex_parameters  class definition */


/* Constructor poequb_dcomplex_parameters definition */
poequb_dcomplex_parameters:: poequb_dcomplex_parameters ( int matrix_layout_i, 
      lapack_int n_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n poequb dcomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = lda*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       poequb_free();
       EXPECT_FALSE( true) << "poequb_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, n, n,
                                            lin_solver_paramslist[idx].symm);
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

poequb_dcomplex_parameters:: ~poequb_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" poequb_dcomplex_parameters object: destructor invoked. \n");
#endif
   poequb_free();
}


//  Test fixture class definition
class zpoequb_test  : public  ::testing::Test {
public:
   poequb_dcomplex_parameters  *zpoequb_obj;
   void SetUp();  
   void TearDown () { delete zpoequb_obj; }
};


void zpoequb_test::SetUp(){

    /* LAPACKE ZPOEQUB prototype */
    typedef int (*Fptr_NL_LAPACKE_zpoequb) ( int matrix_layout, lapack_int n,
                               const lapack_complex_double *a, 
                               lapack_int lda, double *s, 
                               double *scond, double *amax  );

    Fptr_NL_LAPACKE_zpoequb ZPOEQUB;

    zpoequb_obj = new  poequb_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n );
    idx = Circular_Increment_Index(idx);

    zpoequb_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpoequb_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpoequb_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpoequb_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZPOEQUB = (Fptr_NL_LAPACKE_zpoequb)dlsym(zpoequb_obj->hModule, "LAPACKE_zpoequb");
    ASSERT_TRUE(ZPOEQUB != NULL) << "failed to get the Netlib LAPACKE_zpoequb symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zpoequb_obj->inforef = ZPOEQUB( zpoequb_obj->matrix_layout,
                                  zpoequb_obj->n,
                                  zpoequb_obj->aref,
                                  zpoequb_obj->lda,
                                  zpoequb_obj->sref,
                                  &zpoequb_obj->scondref,
                                  &zpoequb_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    zpoequb_obj->info = LAPACKE_zpoequb( zpoequb_obj->matrix_layout,
                                  zpoequb_obj->n,
                                  zpoequb_obj->a,
                                  zpoequb_obj->lda,
                                  zpoequb_obj->s,
                                  &zpoequb_obj->scond,
                                  &zpoequb_obj->amax);

    if( zpoequb_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpoequb \
        is wrong\n", zpoequb_obj->info );
    }
    if( zpoequb_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpoequb is wrong\n",
        zpoequb_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zpoequb_obj->diff +=  computeDiff_d(  zpoequb_obj->n,
                                        zpoequb_obj->s,
                                        zpoequb_obj->sref );
}

TEST_F(zpoequb_test, zpoequb1) {
    EXPECT_NEAR(0.0, zpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequb_obj->scond, zpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequb_obj->amax, zpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpoequb_test, zpoequb2) {
    EXPECT_NEAR(0.0, zpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequb_obj->scond, zpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequb_obj->amax, zpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpoequb_test, zpoequb3) {
    EXPECT_NEAR(0.0, zpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequb_obj->scond, zpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequb_obj->amax, zpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpoequb_test, zpoequb4) {
    EXPECT_NEAR(0.0, zpoequb_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequb_obj->scond, zpoequb_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpoequb_obj->amax, zpoequb_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}
