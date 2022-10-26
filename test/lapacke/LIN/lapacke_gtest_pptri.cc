#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define pptri_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin pptri_float_parameters  class definition */
class pptri_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
      
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      char uplo; //  Must be 'U' or 'L'.

      /* Input/ Output parameters */
      float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pptri_float_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i);
              
      ~pptri_float_parameters (); 
};  /* end of pptri_float_parameters  class definition */


/* Constructor pptri_float_parameters definition */
pptri_float_parameters:: pptri_float_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pptri float:  n: %d, Uplo: %c  \n", n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       pptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

pptri_float_parameters:: ~pptri_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptri_float_parameters object: destructor invoked. \n");
#endif
   pptri_free();
}

//  Test fixture class definition
class spptri_test  : public  ::testing::Test {
public:
   pptri_float_parameters  *spptri_obj;
   void SetUp();  
   void TearDown () { delete spptri_obj; }
};


void spptri_test::SetUp(){

    /* LAPACKE SPPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_spptri) (int matrix_layout, char uplo,
                        lapack_int n,  float * ap );

    Fptr_NL_LAPACKE_spptri SPPTRI;

    typedef int (*Fptr_NL_LAPACKE_spptrf) ( int matrix_layout ,char uplo,
                               lapack_int n, float *ap );
    Fptr_NL_LAPACKE_spptrf SPPTRF;

    spptri_obj = new pptri_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);

    spptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    spptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(spptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(spptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    SPPTRI = (Fptr_NL_LAPACKE_spptri)dlsym(spptri_obj->hModule, "LAPACKE_spptri");
    ASSERT_TRUE(SPPTRI != NULL) << "failed to get the Netlib LAPACKE_spptri symbol";
    
    SPPTRF = (Fptr_NL_LAPACKE_spptrf)dlsym(spptri_obj->hModule,"LAPACKE_spptrf");
    ASSERT_TRUE(SPPTRF != NULL) << "failed to get the Netlib LAPACKE_spptrf symbol";

    /* Pre condition: need to call pptrf - before calling pptri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    spptri_obj->inforef = SPPTRF(   spptri_obj->matrix_layout,
                                    spptri_obj->uplo,
                                    spptri_obj->n,
                                    spptri_obj->aref);

    spptri_obj->inforef = SPPTRI(   spptri_obj->matrix_layout,
                                    spptri_obj->uplo, 
                                    spptri_obj->n,
                                    spptri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    spptri_obj->info  = LAPACKE_spptrf( spptri_obj->matrix_layout,
                                        spptri_obj->uplo,
                                        spptri_obj->n,
                                        spptri_obj->a);

    spptri_obj->info = LAPACKE_spptri(  spptri_obj->matrix_layout, 
                                        spptri_obj->uplo,
                                        spptri_obj->n, 
                                        spptri_obj->a);

    if( spptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_spptri is wrong\n", spptri_obj->info );
    }
    if( spptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_spptri is wrong\n", 
        spptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    spptri_obj->diff =  computeDiff_s( (spptri_obj->n)*(spptri_obj->n+1)/2,
                           spptri_obj->a, spptri_obj->aref );
}

TEST_F(spptri_test, spptri1) {
    EXPECT_NEAR(0.0, spptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spptri_test, spptri2) {
    EXPECT_NEAR(0.0, spptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spptri_test, spptri3) {
    EXPECT_NEAR(0.0, spptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(spptri_test, spptri4) {
    EXPECT_NEAR(0.0, spptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin pptri_double_parameters  class definition */
class pptri_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
      
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      char uplo; //  Must be 'U' or 'L'.

      /* Input/ Output parameters */
      double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pptri_double_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i);
              
      ~pptri_double_parameters (); 
};  /* end of pptri_double_parameters  class definition */


/* Constructor pptri_double_parameters definition */
pptri_double_parameters:: pptri_double_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pptri double:  n: %d, Uplo: %c  \n", n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       pptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

pptri_double_parameters:: ~pptri_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptri_double_parameters object: destructor invoked. \n");
#endif
   pptri_free();
}

//  Test fixture class definition
class dpptri_test  : public  ::testing::Test {
public:
   pptri_double_parameters  *dpptri_obj;
   void SetUp();  
   void TearDown () { delete dpptri_obj; }
};


void dpptri_test::SetUp(){

    /* LAPACKE DPPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_dpptri) (int matrix_layout, char uplo,
                        lapack_int n,  double * ap );

    Fptr_NL_LAPACKE_dpptri DPPTRI;

    typedef int (*Fptr_NL_LAPACKE_dpptrf) ( int matrix_layout ,char uplo,
                               lapack_int n, double *ap );
    Fptr_NL_LAPACKE_dpptrf DPPTRF;

    dpptri_obj = new pptri_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);

    dpptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dpptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dpptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dpptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DPPTRI = (Fptr_NL_LAPACKE_dpptri)dlsym(dpptri_obj->hModule, "LAPACKE_dpptri");
    ASSERT_TRUE(DPPTRI != NULL) << "failed to get the Netlib LAPACKE_dpptri symbol";
    
    DPPTRF = (Fptr_NL_LAPACKE_dpptrf)dlsym(dpptri_obj->hModule,"LAPACKE_dpptrf");
    ASSERT_TRUE(DPPTRF != NULL) << "failed to get the Netlib LAPACKE_dpptrf symbol";

    /* Pre condition: need to call pptrf - before calling pptri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dpptri_obj->inforef = DPPTRF(   dpptri_obj->matrix_layout,
                                    dpptri_obj->uplo,
                                    dpptri_obj->n,
                                    dpptri_obj->aref);

    dpptri_obj->inforef = DPPTRI(   dpptri_obj->matrix_layout,
                                    dpptri_obj->uplo, 
                                    dpptri_obj->n,
                                    dpptri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    dpptri_obj->info  = LAPACKE_dpptrf( dpptri_obj->matrix_layout,
                                        dpptri_obj->uplo,
                                        dpptri_obj->n,
                                        dpptri_obj->a);

    dpptri_obj->info = LAPACKE_dpptri(  dpptri_obj->matrix_layout, 
                                        dpptri_obj->uplo,
                                        dpptri_obj->n, 
                                        dpptri_obj->a);

    if( dpptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dpptri is wrong\n", dpptri_obj->info );
    }
    if( dpptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dpptri is wrong\n", 
        dpptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dpptri_obj->diff =  computeDiff_d( (dpptri_obj->n)*(dpptri_obj->n+1)/2,
                           dpptri_obj->a, dpptri_obj->aref );
}

TEST_F(dpptri_test, dpptri1) {
    EXPECT_NEAR(0.0, dpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpptri_test, dpptri2) {
    EXPECT_NEAR(0.0, dpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpptri_test, dpptri3) {
    EXPECT_NEAR(0.0, dpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dpptri_test, dpptri4) {
    EXPECT_NEAR(0.0, dpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pptri_scomplex_parameters  class definition */
class pptri_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
      
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      char uplo; //  Must be 'U' or 'L'.

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pptri_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i);
              
      ~pptri_scomplex_parameters (); 
};  /* end of pptri_scomplex_parameters  class definition */


/* Constructor pptri_scomplex_parameters definition */
pptri_scomplex_parameters:: pptri_scomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pptri lapack_complex_float:  n: %d, Uplo: %c  \n", n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       pptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

pptri_scomplex_parameters:: ~pptri_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptri_scomplex_parameters object: destructor invoked. \n");
#endif
   pptri_free();
}

//  Test fixture class definition
class cpptri_test  : public  ::testing::Test {
public:
   pptri_scomplex_parameters  *cpptri_obj;
   void SetUp();  
   void TearDown () { delete cpptri_obj; }
};


void cpptri_test::SetUp(){

    /* LAPACKE CPPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_cpptri) (int matrix_layout, char uplo,
                        lapack_int n,  lapack_complex_float * ap );

    Fptr_NL_LAPACKE_cpptri CPPTRI;

    typedef int (*Fptr_NL_LAPACKE_cpptrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *ap );
    Fptr_NL_LAPACKE_cpptrf CPPTRF;

    cpptri_obj = new pptri_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);

    cpptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cpptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cpptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cpptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CPPTRI = (Fptr_NL_LAPACKE_cpptri)dlsym(cpptri_obj->hModule, "LAPACKE_cpptri");
    ASSERT_TRUE(CPPTRI != NULL) << "failed to get the Netlib LAPACKE_cpptri symbol";
    
    CPPTRF = (Fptr_NL_LAPACKE_cpptrf)dlsym(cpptri_obj->hModule,"LAPACKE_cpptrf");
    ASSERT_TRUE(CPPTRF != NULL) << "failed to get the Netlib LAPACKE_cpptrf symbol";

    /* Pre condition: need to call pptrf - before calling pptri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    cpptri_obj->inforef = CPPTRF(   cpptri_obj->matrix_layout,
                                    cpptri_obj->uplo,
                                    cpptri_obj->n,
                                    cpptri_obj->aref);

    cpptri_obj->inforef = CPPTRI(   cpptri_obj->matrix_layout,
                                    cpptri_obj->uplo, 
                                    cpptri_obj->n,
                                    cpptri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    cpptri_obj->info  = LAPACKE_cpptrf( cpptri_obj->matrix_layout,
                                        cpptri_obj->uplo,
                                        cpptri_obj->n,
                                        cpptri_obj->a);

    cpptri_obj->info = LAPACKE_cpptri(  cpptri_obj->matrix_layout, 
                                        cpptri_obj->uplo,
                                        cpptri_obj->n, 
                                        cpptri_obj->a);

    if( cpptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cpptri is wrong\n", cpptri_obj->info );
    }
    if( cpptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cpptri is wrong\n", 
        cpptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    cpptri_obj->diff =  computeDiff_c( (cpptri_obj->n)*(cpptri_obj->n+1)/2,
                           cpptri_obj->a, cpptri_obj->aref );
}

TEST_F(cpptri_test, cpptri1) {
    EXPECT_NEAR(0.0, cpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpptri_test, cpptri2) {
    EXPECT_NEAR(0.0, cpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpptri_test, cpptri3) {
    EXPECT_NEAR(0.0, cpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cpptri_test, cpptri4) {
    EXPECT_NEAR(0.0, cpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin pptri_dcomplex_parameters  class definition */
class pptri_dcomplex_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
      
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      lapack_int n; // The order of A; the number of rows in B
      char uplo; //  Must be 'U' or 'L'.

      /* Input/ Output parameters */
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      pptri_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i);
              
      ~pptri_dcomplex_parameters (); 
};  /* end of pptri_dcomplex_parameters  class definition */


/* Constructor pptri_dcomplex_parameters definition */
pptri_dcomplex_parameters:: pptri_dcomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n pptri lapack_complex_double:  n: %d, Uplo: %c  \n", n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       pptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

pptri_dcomplex_parameters:: ~pptri_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" pptri_dcomplex_parameters object: destructor invoked. \n");
#endif
   pptri_free();
}

//  Test fixture class definition
class zpptri_test  : public  ::testing::Test {
public:
   pptri_dcomplex_parameters  *zpptri_obj;
   void SetUp();  
   void TearDown () { delete zpptri_obj; }
};


void zpptri_test::SetUp(){

    /* LAPACKE ZPPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_zpptri) (int matrix_layout, char uplo,
                        lapack_int n,  lapack_complex_double * ap );

    Fptr_NL_LAPACKE_zpptri ZPPTRI;

    typedef int (*Fptr_NL_LAPACKE_zpptrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *ap );
    Fptr_NL_LAPACKE_zpptrf ZPPTRF;

    zpptri_obj = new pptri_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);

    zpptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zpptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zpptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zpptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZPPTRI = (Fptr_NL_LAPACKE_zpptri)dlsym(zpptri_obj->hModule, "LAPACKE_zpptri");
    ASSERT_TRUE(ZPPTRI != NULL) << "failed to get the Netlib LAPACKE_zpptri symbol";
    
    ZPPTRF = (Fptr_NL_LAPACKE_zpptrf)dlsym(zpptri_obj->hModule,"LAPACKE_zpptrf");
    ASSERT_TRUE(ZPPTRF != NULL) << "failed to get the Netlib LAPACKE_zpptrf symbol";

    /* Pre condition: need to call pptrf - before calling pptri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zpptri_obj->inforef = ZPPTRF(   zpptri_obj->matrix_layout,
                                    zpptri_obj->uplo,
                                    zpptri_obj->n,
                                    zpptri_obj->aref);

    zpptri_obj->inforef = ZPPTRI(   zpptri_obj->matrix_layout,
                                    zpptri_obj->uplo, 
                                    zpptri_obj->n,
                                    zpptri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    zpptri_obj->info  = LAPACKE_zpptrf( zpptri_obj->matrix_layout,
                                        zpptri_obj->uplo,
                                        zpptri_obj->n,
                                        zpptri_obj->a);

    zpptri_obj->info = LAPACKE_zpptri(  zpptri_obj->matrix_layout, 
                                        zpptri_obj->uplo,
                                        zpptri_obj->n, 
                                        zpptri_obj->a);

    if( zpptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zpptri is wrong\n", zpptri_obj->info );
    }
    if( zpptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zpptri is wrong\n", 
        zpptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zpptri_obj->diff =  computeDiff_z( (zpptri_obj->n)*(zpptri_obj->n+1)/2,
                           zpptri_obj->a, zpptri_obj->aref );
}

TEST_F(zpptri_test, zpptri1) {
    EXPECT_NEAR(0.0, zpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpptri_test, zpptri2) {
    EXPECT_NEAR(0.0, zpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpptri_test, zpptri3) {
    EXPECT_NEAR(0.0, zpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zpptri_test, zpptri4) {
    EXPECT_NEAR(0.0, zpptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
