#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define potri_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin potri_float_parameters  class definition */
class potri_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      potri_float_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~potri_float_parameters (); 
};  /* end of potri_float_parameters  class definition */


/* Constructor potri_float_parameters definition */
potri_float_parameters:: potri_float_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n potri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       potri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

potri_float_parameters:: ~potri_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potri_float_parameters object: destructor invoked. \n");
#endif
   potri_free();
}

//  Test fixture class definition
class spotri_test  : public  ::testing::Test {
public:
   potri_float_parameters  *spotri_obj;
   void SetUp();  
   void TearDown () { delete spotri_obj; }
};


void spotri_test::SetUp(){

    /* LAPACKE SPOTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_spotri) (int matrix_layout, char uplo,
                        lapack_int n,  float * a, lapack_int lda  );

    Fptr_NL_LAPACKE_spotri SPOTRI;

    typedef int (*Fptr_NL_LAPACKE_spotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, float *a, lapack_int lda );
    Fptr_NL_LAPACKE_spotrf SPOTRF;

    spotri_obj = new potri_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    spotri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spotri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spotri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spotri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SPOTRI = (Fptr_NL_LAPACKE_spotri)dlsym(spotri_obj->hModule, "LAPACKE_spotri");
    ASSERT_TRUE(SPOTRI != NULL) << "failed to get the Netlib LAPACKE_spotri symbol";
    
    SPOTRF = (Fptr_NL_LAPACKE_spotrf)dlsym(spotri_obj->hModule,"LAPACKE_spotrf");
    ASSERT_TRUE(SPOTRF != NULL) << "failed to get the Netlib LAPACKE_spotrf symbol";

    /* Pre condition: need to call potrf - before calling potri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    spotri_obj->inforef = SPOTRF( 	spotri_obj->matrix_layout,
									spotri_obj->uplo,
									spotri_obj->n,
									spotri_obj->aref,
									spotri_obj->lda);

    spotri_obj->inforef = SPOTRI( 	spotri_obj->matrix_layout,
									spotri_obj->uplo, 
									spotri_obj->n,
									spotri_obj->aref,
									spotri_obj->lda);
                          

    /* Compute libflame's Lapacke o/p  */
    spotri_obj->info  = LAPACKE_spotrf( spotri_obj->matrix_layout,
										spotri_obj->uplo,
										spotri_obj->n,
										spotri_obj->a,
										spotri_obj->lda);

    spotri_obj->info = LAPACKE_spotri( 	spotri_obj->matrix_layout, 
										spotri_obj->uplo,
										spotri_obj->n, 
										spotri_obj->a, 
										spotri_obj->lda);

    if( spotri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_spotri is wrong\n", spotri_obj->info );
    }
    if( spotri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spotri is wrong\n", 
        spotri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    spotri_obj->diff =  computeDiff_s( (spotri_obj->n)*(spotri_obj->lda),
						   spotri_obj->a, spotri_obj->aref );
}

TEST_F(spotri_test, spotri1) {
    EXPECT_NEAR(0.0, spotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spotri_test, spotri2) {
    EXPECT_NEAR(0.0, spotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spotri_test, spotri3) {
    EXPECT_NEAR(0.0, spotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spotri_test, spotri4) {
    EXPECT_NEAR(0.0, spotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin potri_double_parameters  class definition */
class potri_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      potri_double_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~potri_double_parameters (); 
};  /* end of potri_double_parameters  class definition */


/* Constructor potri_double_parameters definition */
potri_double_parameters:: potri_double_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n potri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       potri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

potri_double_parameters:: ~potri_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potri_double_parameters object: destructor invoked. \n");
#endif
   potri_free();
}

//  Test fixture class definition
class dpotri_test  : public  ::testing::Test {
public:
   potri_double_parameters  *dpotri_obj;
   void SetUp();  
   void TearDown () { delete dpotri_obj; }
};


void dpotri_test::SetUp(){

    /* LAPACKE DPOTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_dpotri) (int matrix_layout, char uplo,
                        lapack_int n,  double * a, lapack_int lda  );

    Fptr_NL_LAPACKE_dpotri DPOTRI;

    typedef int (*Fptr_NL_LAPACKE_dpotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, double *a, lapack_int lda );
    Fptr_NL_LAPACKE_dpotrf DPOTRF;

    dpotri_obj = new potri_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    dpotri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpotri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpotri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpotri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DPOTRI = (Fptr_NL_LAPACKE_dpotri)dlsym(dpotri_obj->hModule, "LAPACKE_dpotri");
    ASSERT_TRUE(DPOTRI != NULL) << "failed to get the Netlib LAPACKE_dpotri symbol";
    
    DPOTRF = (Fptr_NL_LAPACKE_dpotrf)dlsym(dpotri_obj->hModule,"LAPACKE_dpotrf");
    ASSERT_TRUE(DPOTRF != NULL) << "failed to get the Netlib LAPACKE_dpotrf symbol";

    /* Pre condition: need to call potrf - before calling potri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dpotri_obj->inforef = DPOTRF( 	dpotri_obj->matrix_layout,
									dpotri_obj->uplo,
									dpotri_obj->n,
									dpotri_obj->aref,
									dpotri_obj->lda);

    dpotri_obj->inforef = DPOTRI( 	dpotri_obj->matrix_layout,
									dpotri_obj->uplo, 
									dpotri_obj->n,
									dpotri_obj->aref,
									dpotri_obj->lda);
                          

    /* Compute libflame's Lapacke o/p  */
    dpotri_obj->info  = LAPACKE_dpotrf( dpotri_obj->matrix_layout,
										dpotri_obj->uplo,
										dpotri_obj->n,
										dpotri_obj->a,
										dpotri_obj->lda);

    dpotri_obj->info = LAPACKE_dpotri( 	dpotri_obj->matrix_layout, 
										dpotri_obj->uplo,
										dpotri_obj->n, 
										dpotri_obj->a, 
										dpotri_obj->lda);

    if( dpotri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dpotri is wrong\n", dpotri_obj->info );
    }
    if( dpotri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpotri is wrong\n", 
        dpotri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dpotri_obj->diff =  computeDiff_d( (dpotri_obj->n)*(dpotri_obj->lda),
						   dpotri_obj->a, dpotri_obj->aref );
}

TEST_F(dpotri_test, dpotri1) {
    EXPECT_NEAR(0.0, dpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpotri_test, dpotri2) {
    EXPECT_NEAR(0.0, dpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpotri_test, dpotri3) {
    EXPECT_NEAR(0.0, dpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpotri_test, dpotri4) {
    EXPECT_NEAR(0.0, dpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin potri_lapack_complex_float_parameters  class definition */
class potri_lapack_complex_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      potri_lapack_complex_float_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~potri_lapack_complex_float_parameters (); 
};  /* end of potri_lapack_complex_float_parameters  class definition */


/* Constructor potri_lapack_complex_float_parameters definition */
potri_lapack_complex_float_parameters:: potri_lapack_complex_float_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n potri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_lapack_complex_float_parameters object: malloc error.";
       potri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

potri_lapack_complex_float_parameters:: ~potri_lapack_complex_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potri_lapack_complex_float_parameters object: destructor invoked. \n");
#endif
   potri_free();
}

//  Test fixture class definition
class cpotri_test  : public  ::testing::Test {
public:
   potri_lapack_complex_float_parameters  *cpotri_obj;
   void SetUp();  
   void TearDown () { delete cpotri_obj; }
};


void cpotri_test::SetUp(){

    /* LAPACKE CPOTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_cpotri) (int matrix_layout, char uplo,
                        lapack_int n,  lapack_complex_float * a, lapack_int lda  );

    Fptr_NL_LAPACKE_cpotri CPOTRI;

    typedef int (*Fptr_NL_LAPACKE_cpotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda );
    Fptr_NL_LAPACKE_cpotrf CPOTRF;

    cpotri_obj = new potri_lapack_complex_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    cpotri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpotri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpotri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpotri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CPOTRI = (Fptr_NL_LAPACKE_cpotri)dlsym(cpotri_obj->hModule, "LAPACKE_cpotri");
    ASSERT_TRUE(CPOTRI != NULL) << "failed to get the Netlib LAPACKE_cpotri symbol";
    
    CPOTRF = (Fptr_NL_LAPACKE_cpotrf)dlsym(cpotri_obj->hModule,"LAPACKE_cpotrf");
    ASSERT_TRUE(CPOTRF != NULL) << "failed to get the Netlib LAPACKE_cpotrf symbol";

    /* Pre condition: need to call potrf - before calling potri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    cpotri_obj->inforef = CPOTRF( 	cpotri_obj->matrix_layout,
									cpotri_obj->uplo,
									cpotri_obj->n,
									cpotri_obj->aref,
									cpotri_obj->lda);

    cpotri_obj->inforef = CPOTRI( 	cpotri_obj->matrix_layout,
									cpotri_obj->uplo, 
									cpotri_obj->n,
									cpotri_obj->aref,
									cpotri_obj->lda);
                          

    /* Compute libflame's Lapacke o/p  */
    cpotri_obj->info  = LAPACKE_cpotrf( cpotri_obj->matrix_layout,
										cpotri_obj->uplo,
										cpotri_obj->n,
										cpotri_obj->a,
										cpotri_obj->lda);

    cpotri_obj->info = LAPACKE_cpotri( 	cpotri_obj->matrix_layout, 
										cpotri_obj->uplo,
										cpotri_obj->n, 
										cpotri_obj->a, 
										cpotri_obj->lda);

    if( cpotri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpotri is wrong\n", cpotri_obj->info );
    }
    if( cpotri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpotri is wrong\n", 
        cpotri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cpotri_obj->diff =  computeDiff_c( (cpotri_obj->n)*(cpotri_obj->lda),
						   cpotri_obj->a, cpotri_obj->aref );
}

TEST_F(cpotri_test, cpotri1) {
    EXPECT_NEAR(0.0, cpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpotri_test, cpotri2) {
    EXPECT_NEAR(0.0, cpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpotri_test, cpotri3) {
    EXPECT_NEAR(0.0, cpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpotri_test, cpotri4) {
    EXPECT_NEAR(0.0, cpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin potri_lapack_complex_double_parameters  class definition */
class potri_lapack_complex_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      potri_lapack_complex_double_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~potri_lapack_complex_double_parameters (); 
};  /* end of potri_lapack_complex_double_parameters  class definition */


/* Constructor potri_lapack_complex_double_parameters definition */
potri_lapack_complex_double_parameters:: potri_lapack_complex_double_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n potri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_lapack_complex_double_parameters object: malloc error.";
       potri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

potri_lapack_complex_double_parameters:: ~potri_lapack_complex_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" potri_lapack_complex_double_parameters object: destructor invoked. \n");
#endif
   potri_free();
}

//  Test fixture class definition
class zpotri_test  : public  ::testing::Test {
public:
   potri_lapack_complex_double_parameters  *zpotri_obj;
   void SetUp();  
   void TearDown () { delete zpotri_obj; }
};


void zpotri_test::SetUp(){

    /* LAPACKE ZPOTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_zpotri) (int matrix_layout, char uplo,
                        lapack_int n,  lapack_complex_double * a, lapack_int lda  );

    Fptr_NL_LAPACKE_zpotri ZPOTRI;

    typedef int (*Fptr_NL_LAPACKE_zpotrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda );
    Fptr_NL_LAPACKE_zpotrf ZPOTRF;

    zpotri_obj = new potri_lapack_complex_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    zpotri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpotri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpotri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpotri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZPOTRI = (Fptr_NL_LAPACKE_zpotri)dlsym(zpotri_obj->hModule, "LAPACKE_zpotri");
    ASSERT_TRUE(ZPOTRI != NULL) << "failed to get the Netlib LAPACKE_zpotri symbol";
    
    ZPOTRF = (Fptr_NL_LAPACKE_zpotrf)dlsym(zpotri_obj->hModule,"LAPACKE_zpotrf");
    ASSERT_TRUE(ZPOTRF != NULL) << "failed to get the Netlib LAPACKE_zpotrf symbol";

    /* Pre condition: need to call potrf - before calling potri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zpotri_obj->inforef = ZPOTRF( 	zpotri_obj->matrix_layout,
									zpotri_obj->uplo,
									zpotri_obj->n,
									zpotri_obj->aref,
									zpotri_obj->lda);

    zpotri_obj->inforef = ZPOTRI( 	zpotri_obj->matrix_layout,
									zpotri_obj->uplo, 
									zpotri_obj->n,
									zpotri_obj->aref,
									zpotri_obj->lda);
                          

    /* Compute libflame's Lapacke o/p  */
    zpotri_obj->info  = LAPACKE_zpotrf( zpotri_obj->matrix_layout,
										zpotri_obj->uplo,
										zpotri_obj->n,
										zpotri_obj->a,
										zpotri_obj->lda);

    zpotri_obj->info = LAPACKE_zpotri( 	zpotri_obj->matrix_layout, 
										zpotri_obj->uplo,
										zpotri_obj->n, 
										zpotri_obj->a, 
										zpotri_obj->lda);

    if( zpotri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zpotri is wrong\n", zpotri_obj->info );
    }
    if( zpotri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpotri is wrong\n", 
        zpotri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zpotri_obj->diff =  computeDiff_z( (zpotri_obj->n)*(zpotri_obj->lda),
						   zpotri_obj->a, zpotri_obj->aref );
}

TEST_F(zpotri_test, zpotri1) {
    EXPECT_NEAR(0.0, zpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpotri_test, zpotri2) {
    EXPECT_NEAR(0.0, zpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpotri_test, zpotri3) {
    EXPECT_NEAR(0.0, zpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpotri_test, zpotri4) {
    EXPECT_NEAR(0.0, zpotri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
