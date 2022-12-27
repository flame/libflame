#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"

#define pbequ_free() \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (s != NULL)    free (s   ); \
       if (sref != NULL) free (sref); \
       if( hModule != NULL) dlclose(hModule); \
       if( dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin pbequ_float_parameters  class definition */
class pbequ_float_parameters{
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
	  lapack_int kd; //  number of superdiagonals or subdiagonals
      /* Output parameters */
      float *s, *sref; // scale factors, of size 'n'
      float scond, scondref; // ratio of the smallest s[i] to the largest s[i]
      float amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      pbequ_float_parameters ( int matrix_layout_i, 
               lapack_int n_i, char uplo_i, lapack_int kd_i );
             
      ~pbequ_float_parameters ();
};  /* end of pbequ_float_parameters  class definition */


/* Constructor pbequ_float_parameters definition */
pbequ_float_parameters:: pbequ_float_parameters ( int matrix_layout_i, 
      lapack_int n_i, char uplo_i, lapack_int kd_i ) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	kd = kd_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ fla_max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n pbequ float:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = n*n;

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       pbequ_free();
       EXPECT_FALSE( true) << "pbequ_float_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref,
                                                             n, n, 'S');

    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

pbequ_float_parameters:: ~pbequ_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbequ_float_parameters object: destructor invoked. \n");
#endif
   pbequ_free();
}


//  Test fixture class definition
class spbequ_test  : public  ::testing::Test {
public:
   pbequ_float_parameters  *spbequ_obj;
   void SetUp();  
   void TearDown () { delete spbequ_obj; }
};


void spbequ_test::SetUp(){

    /* LAPACKE SPBEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_spbequ) (   int matrix_layout, char uplo,
                    lapack_int n, lapack_int kd, const float *a, lapack_int lda, float* s,
                                    float* scond, float* amax  );

    Fptr_NL_LAPACKE_spbequ SPBEQU;

    spbequ_obj = new  pbequ_float_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n, 
                                         lin_solver_paramslist[idx].Uplo,
										 lin_solver_paramslist[idx].kd);
    idx = Circular_Increment_Index(idx);

    spbequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spbequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spbequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spbequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SPBEQU = (Fptr_NL_LAPACKE_spbequ)dlsym(spbequ_obj->hModule, "LAPACKE_spbequ");
    ASSERT_TRUE(SPBEQU != NULL) << "failed to get the Netlib LAPACKE_spbequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    spbequ_obj->inforef = SPBEQU( spbequ_obj->matrix_layout,
                                  spbequ_obj->uplo,
                                  spbequ_obj->n,
                                  spbequ_obj->kd,
                                  spbequ_obj->aref,
                                  spbequ_obj->lda,
                                  spbequ_obj->sref,
                                  &spbequ_obj->scondref,
                                  &spbequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    spbequ_obj->info = LAPACKE_spbequ( spbequ_obj->matrix_layout,
                                  spbequ_obj->uplo,
                                  spbequ_obj->n,
                                  spbequ_obj->kd,
                                  spbequ_obj->a,
                                  spbequ_obj->lda,
                                  spbequ_obj->s,
                                  &spbequ_obj->scond,
                                  &spbequ_obj->amax);

    if( spbequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_spbequ \
        is wrong\n", spbequ_obj->info );
    }
    if( spbequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spbequ is wrong\n",
        spbequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    spbequ_obj->diff =  computeDiff_s(  spbequ_obj->n,
                                        spbequ_obj->s,
                                        spbequ_obj->sref );
}

TEST_F(spbequ_test, spbequ1) {
    EXPECT_NEAR(0.0, spbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spbequ_obj->scond, spbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spbequ_obj->amax, spbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbequ_test, spbequ2) {
    EXPECT_NEAR(0.0, spbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spbequ_obj->scond, spbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spbequ_obj->amax, spbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbequ_test, spbequ3) {
    EXPECT_NEAR(0.0, spbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spbequ_obj->scond, spbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spbequ_obj->amax, spbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spbequ_test, spbequ4) {
    EXPECT_NEAR(0.0, spbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spbequ_obj->scond, spbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(spbequ_obj->amax, spbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pbequ_double_parameters  class definition */
class pbequ_double_parameters{
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
	  lapack_int kd; //  number of superdiagonals or subdiagonals
      /* Output parameters */
      double * s, *sref; // colum scale factors, , array if size 'n'
      double scond, scondref; // ratio of the smallest c[i] to the largest c[i]
      double amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      pbequ_double_parameters ( int matrix_layout_i, 
               lapack_int n_i, char uplo_i, lapack_int kd_i );
             
      ~pbequ_double_parameters ();
};  /* end of pbequ_double_parameters  class definition */


/* Constructor pbequ_double_parameters definition */
pbequ_double_parameters:: pbequ_double_parameters ( int matrix_layout_i, 
      lapack_int n_i, char uplo_i, lapack_int kd_i ) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	kd = kd_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ fla_max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n pbequ double:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize =n*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       pbequ_free();
       EXPECT_FALSE( true) << "pbequ_double_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref,
                                                             n, n, 'S');
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

pbequ_double_parameters:: ~pbequ_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbequ_double_parameters object: destructor invoked. \n");
#endif
   pbequ_free();
}


//  Test fixture class definition
class dpbequ_test  : public  ::testing::Test {
public:
   pbequ_double_parameters  *dpbequ_obj;
   void SetUp();  
   void TearDown () { delete dpbequ_obj; }
};


void dpbequ_test::SetUp(){

    /* LAPACKE DPBEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_dpbequ) ( int matrix_layout,char uplo, 
                 lapack_int n, lapack_int kd, const double *a,  lapack_int lda,
                 double *s, double *scond, double *amax  );

    Fptr_NL_LAPACKE_dpbequ DPBEQU;

    dpbequ_obj = new  pbequ_double_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n, 
                                         lin_solver_paramslist[idx].Uplo,
										 lin_solver_paramslist[idx].kd);
    idx = Circular_Increment_Index(idx);

    dpbequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpbequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpbequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpbequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DPBEQU = (Fptr_NL_LAPACKE_dpbequ)dlsym(dpbequ_obj->hModule, "LAPACKE_dpbequ");
    ASSERT_TRUE(DPBEQU != NULL) << "failed to get the Netlib LAPACKE_dpbequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    dpbequ_obj->inforef = DPBEQU( dpbequ_obj->matrix_layout,
                                  dpbequ_obj->uplo,
                                  dpbequ_obj->n,
                                  dpbequ_obj->kd,
                                  dpbequ_obj->aref,
                                  dpbequ_obj->lda,
                                  dpbequ_obj->sref,
                                  &dpbequ_obj->scondref,
                                  &dpbequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    dpbequ_obj->info = LAPACKE_dpbequ( dpbequ_obj->matrix_layout,
                                  dpbequ_obj->uplo,
                                  dpbequ_obj->n,
                                  dpbequ_obj->kd,
                                  dpbequ_obj->a,
                                  dpbequ_obj->lda,
                                  dpbequ_obj->s,
                                  &dpbequ_obj->scond,
                                  &dpbequ_obj->amax);

    if( dpbequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dpbequ \
        is wrong\n", dpbequ_obj->info );
    }
    if( dpbequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpbequ is wrong\n",
        dpbequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dpbequ_obj->diff =  computeDiff_d(  dpbequ_obj->n,
                                        dpbequ_obj->s,
                                        dpbequ_obj->sref );
}

TEST_F(dpbequ_test, dpbequ1) {
    EXPECT_NEAR(0.0, dpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpbequ_obj->scond, dpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpbequ_obj->amax, dpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbequ_test, dpbequ2) {
    EXPECT_NEAR(0.0, dpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpbequ_obj->scond, dpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpbequ_obj->amax, dpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbequ_test, dpbequ3) {
    EXPECT_NEAR(0.0, dpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpbequ_obj->scond, dpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpbequ_obj->amax, dpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpbequ_test, dpbequ4) {
    EXPECT_NEAR(0.0, dpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpbequ_obj->scond, dpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(dpbequ_obj->amax, dpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pbequ_scomplex_parameters  class definition */
class pbequ_scomplex_parameters{
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
	  lapack_int kd; //  number of superdiagonals or subdiagonals
      /* Output parameters */
      float * s, *sref; // colum scale factors, , array if size 'n'
      float scond, scondref; // ratio of the smallest c[i] to the largest c[i]
      float amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      pbequ_scomplex_parameters ( int matrix_layout_i, 
               lapack_int n_i, char uplo_i, lapack_int kd_i );
             
      ~pbequ_scomplex_parameters ();
};  /* end of pbequ_scomplex_parameters  class definition */


/* Constructor pbequ_scomplex_parameters definition */
pbequ_scomplex_parameters:: pbequ_scomplex_parameters ( int matrix_layout_i, 
      lapack_int n_i, char uplo_i, lapack_int kd_i ) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	kd = kd_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ fla_max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n pbequ scomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = n*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_float_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       pbequ_free();
       EXPECT_FALSE( true) << "pbequ_scomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref,
                                                             n, n, 'S');
    lapacke_gtest_init_float_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

pbequ_scomplex_parameters:: ~pbequ_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbequ_scomplex_parameters object: destructor invoked. \n");
#endif
   pbequ_free();
}


//  Test fixture class definition
class cpbequ_test  : public  ::testing::Test {
public:
   pbequ_scomplex_parameters  *cpbequ_obj;
   void SetUp();  
   void TearDown () { delete cpbequ_obj; }
};


void cpbequ_test::SetUp(){

    /* LAPACKE CPBEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_cpbequ) ( int matrix_layout, char uplo, 
                            lapack_int n, lapack_int kd, const lapack_complex_float *a, 
                 lapack_int lda, float *s, float *scond, float *amax  );

    Fptr_NL_LAPACKE_cpbequ CPBEQU;

    cpbequ_obj = new  pbequ_scomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n, 
                                         lin_solver_paramslist[idx].Uplo,
										 lin_solver_paramslist[idx].kd);
    idx = Circular_Increment_Index(idx);

    cpbequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpbequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpbequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpbequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CPBEQU = (Fptr_NL_LAPACKE_cpbequ)dlsym(cpbequ_obj->hModule, "LAPACKE_cpbequ");
    ASSERT_TRUE(CPBEQU != NULL) << "failed to get the Netlib LAPACKE_cpbequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    cpbequ_obj->inforef = CPBEQU( cpbequ_obj->matrix_layout,
                                  cpbequ_obj->uplo,
                                  cpbequ_obj->n,
                                  cpbequ_obj->kd,
                                  cpbequ_obj->aref,
                                  cpbequ_obj->lda,
                                  cpbequ_obj->sref,
                                  &cpbequ_obj->scondref,
                                  &cpbequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    cpbequ_obj->info = LAPACKE_cpbequ( cpbequ_obj->matrix_layout,
                                  cpbequ_obj->uplo,
                                  cpbequ_obj->n,
                                  cpbequ_obj->kd,
                                  cpbequ_obj->a,
                                  cpbequ_obj->lda,
                                  cpbequ_obj->s,
                                  &cpbequ_obj->scond,
                                  &cpbequ_obj->amax);

    if( cpbequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_cpbequ \
        is wrong\n", cpbequ_obj->info );
    }
    if( cpbequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpbequ is wrong\n",
        cpbequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cpbequ_obj->diff =  computeDiff_s(  cpbequ_obj->n,
                                        cpbequ_obj->s,
                                        cpbequ_obj->sref );
}

TEST_F(cpbequ_test, cpbequ1) {
    EXPECT_NEAR(0.0, cpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpbequ_obj->scond, cpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpbequ_obj->amax, cpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbequ_test, cpbequ2) {
    EXPECT_NEAR(0.0, cpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpbequ_obj->scond, cpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpbequ_obj->amax, cpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbequ_test, cpbequ3) {
    EXPECT_NEAR(0.0, cpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpbequ_obj->scond, cpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpbequ_obj->amax, cpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpbequ_test, cpbequ4) {
    EXPECT_NEAR(0.0, cpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpbequ_obj->scond, cpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(cpbequ_obj->amax, cpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}


/* Begin pbequ_dcomplex_parameters  class definition */
class pbequ_dcomplex_parameters{
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
	  lapack_int kd; //  number of superdiagonals or subdiagonals
      /* Output parameters */
      double * s, *sref; // colum scale factors, , array if size 'n'
      double scond, scondref; // ratio of the smallest c[i] to the largest c[i]
      double amax, amaxref; // Absolute value of the largest element of the matrix A.

      /* Return Values */
      lapack_int info, inforef;

   public:
      pbequ_dcomplex_parameters ( int matrix_layout_i, 
               lapack_int n_i, char uplo_i, lapack_int kd_i );
             
      ~pbequ_dcomplex_parameters ();
};  /* end of pbequ_dcomplex_parameters  class definition */


/* Constructor pbequ_dcomplex_parameters definition */
pbequ_dcomplex_parameters:: pbequ_dcomplex_parameters ( int matrix_layout_i, 
      lapack_int n_i, char uplo_i, lapack_int kd_i ) {
    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	kd = kd_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;

    lda = n; // as per API spec, lda≥ fla_max(1, n).

#if LAPACKE_TEST_VERBOSE
   printf(" \n pbequ dcomplex:  n: %d lda: %d \n", n, lda);
#endif
    a_bufsize = n*n;
    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, a_bufsize);
    lapacke_gtest_alloc_double_buffer_pair( &s, &sref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (s==NULL) || (sref==NULL) ){
       pbequ_free();
       EXPECT_FALSE( true) << "pbequ_dcomplex_parameters object: malloc error.";
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref,
                                                             n, n, 'S');
    lapacke_gtest_init_double_buffer_pair_with_constant(s, sref, n, 0.0);

   } /* end of Constructor  */

pbequ_dcomplex_parameters:: ~pbequ_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pbequ_dcomplex_parameters object: destructor invoked. \n");
#endif
   pbequ_free();
}


//  Test fixture class definition
class zpbequ_test  : public  ::testing::Test {
public:
   pbequ_dcomplex_parameters  *zpbequ_obj;
   void SetUp();  
   void TearDown () { delete zpbequ_obj; }
};


void zpbequ_test::SetUp(){

    /* LAPACKE ZPBEQU prototype */
    typedef int (*Fptr_NL_LAPACKE_zpbequ) ( int matrix_layout, char uplo, 
                            lapack_int n, lapack_int kd, const lapack_complex_double *a, 
            lapack_int lda, double *s, double *scond, double *amax );

    Fptr_NL_LAPACKE_zpbequ ZPBEQU;

    zpbequ_obj = new  pbequ_dcomplex_parameters(lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].n,  
                                         lin_solver_paramslist[idx].Uplo,
										 lin_solver_paramslist[idx].kd);
    idx = Circular_Increment_Index(idx);

    zpbequ_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpbequ_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpbequ_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpbequ_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZPBEQU = (Fptr_NL_LAPACKE_zpbequ)dlsym(zpbequ_obj->hModule, "LAPACKE_zpbequ");
    ASSERT_TRUE(ZPBEQU != NULL) << "failed to get the Netlib LAPACKE_zpbequ symbol";

    /* Compute the Netlib-Lapacke's reference o/p  */
    zpbequ_obj->inforef = ZPBEQU( zpbequ_obj->matrix_layout,
                                  zpbequ_obj->uplo,
                                  zpbequ_obj->n,
                                  zpbequ_obj->kd,
                                  zpbequ_obj->aref,
                                  zpbequ_obj->lda,
                                  zpbequ_obj->sref,
                                  &zpbequ_obj->scondref,
                                  &zpbequ_obj->amaxref);

    /* Compute the Libflme lapacke o/p  */
    zpbequ_obj->info = LAPACKE_zpbequ( zpbequ_obj->matrix_layout,
                                  zpbequ_obj->uplo,
                                  zpbequ_obj->n,
                                  zpbequ_obj->kd,
                                  zpbequ_obj->a,
                                  zpbequ_obj->lda,
                                  zpbequ_obj->s,
                                  &zpbequ_obj->scond,
                                  &zpbequ_obj->amax);

    if( zpbequ_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zpbequ \
        is wrong\n", zpbequ_obj->info );
    }
    if( zpbequ_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpbequ is wrong\n",
        zpbequ_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zpbequ_obj->diff +=  computeDiff_d(  zpbequ_obj->n,
                                        zpbequ_obj->s,
                                        zpbequ_obj->sref );
}

TEST_F(zpbequ_test, zpbequ1) {
    EXPECT_NEAR(0.0, zpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpbequ_obj->scond, zpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpbequ_obj->amax, zpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbequ_test, zpbequ2) {
    EXPECT_NEAR(0.0, zpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpbequ_obj->scond, zpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpbequ_obj->amax, zpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbequ_test, zpbequ3) {
    EXPECT_NEAR(0.0, zpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpbequ_obj->scond, zpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpbequ_obj->amax, zpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpbequ_test, zpbequ4) {
    EXPECT_NEAR(0.0, zpbequ_obj->diff, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpbequ_obj->scond, zpbequ_obj->scondref, LAPACKE_GTEST_THRESHOLD);
    EXPECT_NEAR(zpbequ_obj->amax, zpbequ_obj->amaxref, LAPACKE_GTEST_THRESHOLD);
}
