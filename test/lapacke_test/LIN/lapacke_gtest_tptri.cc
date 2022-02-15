#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define tptri_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin tptri_float_parameters  class definition */
class tptri_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B

      /* Input/ Output parameters */
      float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tptri_float_parameters ( int matrix_layout_i, char uplo_i, 
                        char diag_i, lapack_int n_i);
              
      ~tptri_float_parameters (); 
};  /* end of tptri_float_parameters  class definition */


/* Constructor tptri_float_parameters definition */
tptri_float_parameters:: tptri_float_parameters ( int matrix_layout_i, 
            char uplo_i, char diag_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	diag = diag_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tptri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       tptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

tptri_float_parameters:: ~tptri_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tptri_float_parameters object: destructor invoked. \n");
#endif
   tptri_free();
}

//  Test fixture class definition
class stptri_test  : public  ::testing::Test {
public:
   tptri_float_parameters  *stptri_obj;
   void SetUp();  
   void TearDown () { delete stptri_obj; }
};


void stptri_test::SetUp(){

    /* LAPACKE STPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_stptri) (int matrix_layout, char uplo,
                   char diag, lapack_int n,  float * a );

    Fptr_NL_LAPACKE_stptri STPTRI;

    idx = Circular_Increment_Index(idx);
    stptri_obj = new tptri_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n);


    stptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    STPTRI = (Fptr_NL_LAPACKE_stptri)dlsym(stptri_obj->hModule, "LAPACKE_stptri");
    ASSERT_TRUE(STPTRI != NULL) << "failed to get the Netlib LAPACKE_stptri symbol";
    
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    stptri_obj->inforef = STPTRI( 	stptri_obj->matrix_layout,
									stptri_obj->uplo, 
									stptri_obj->diag, 
									stptri_obj->n,
									stptri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    stptri_obj->info = LAPACKE_stptri( 	stptri_obj->matrix_layout, 
										stptri_obj->uplo,
										stptri_obj->diag, 
										stptri_obj->n, 
										stptri_obj->a);

    if( stptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stptri is wrong\n", stptri_obj->info );
    }
    if( stptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stptri is wrong\n", 
        stptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    stptri_obj->diff =  computeDiff_s( (stptri_obj->n)*(stptri_obj->n),
						   stptri_obj->a, stptri_obj->aref );
}

TEST_F(stptri_test, stptri1) {
    EXPECT_NEAR(0.0, stptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stptri_test, stptri2) {
    EXPECT_NEAR(0.0, stptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stptri_test, stptri3) {
    EXPECT_NEAR(0.0, stptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stptri_test, stptri4) {
    EXPECT_NEAR(0.0, stptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin tptri_double_parameters  class definition */
class tptri_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B

      /* Input/ Output parameters */
      double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tptri_double_parameters ( int matrix_layout_i, char uplo_i, 
                        char diag_i, lapack_int n_i);
              
      ~tptri_double_parameters (); 
};  /* end of tptri_double_parameters  class definition */


/* Constructor tptri_double_parameters definition */
tptri_double_parameters:: tptri_double_parameters ( int matrix_layout_i, 
            char uplo_i, char diag_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	diag = diag_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tptri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       tptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

tptri_double_parameters:: ~tptri_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tptri_double_parameters object: destructor invoked. \n");
#endif
   tptri_free();
}

//  Test fixture class definition
class dtptri_test  : public  ::testing::Test {
public:
   tptri_double_parameters  *dtptri_obj;
   void SetUp();  
   void TearDown () { delete dtptri_obj; }
};


void dtptri_test::SetUp(){

    /* LAPACKE DTPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_dtptri) (int matrix_layout, char uplo,
                   char diag, lapack_int n,  double * a );

    Fptr_NL_LAPACKE_dtptri DTPTRI;

    dtptri_obj = new tptri_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    dtptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DTPTRI = (Fptr_NL_LAPACKE_dtptri)dlsym(dtptri_obj->hModule, "LAPACKE_dtptri");
    ASSERT_TRUE(DTPTRI != NULL) << "failed to get the Netlib LAPACKE_dtptri symbol";

    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dtptri_obj->inforef = DTPTRI( 	dtptri_obj->matrix_layout,
									dtptri_obj->uplo, 
									dtptri_obj->diag, 
									dtptri_obj->n,
									dtptri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    dtptri_obj->info = LAPACKE_dtptri( 	dtptri_obj->matrix_layout, 
										dtptri_obj->uplo,
										dtptri_obj->diag, 
										dtptri_obj->n, 
										dtptri_obj->a);

    if( dtptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtptri is wrong\n", dtptri_obj->info );
    }
    if( dtptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtptri is wrong\n", 
        dtptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dtptri_obj->diff =  computeDiff_d( (dtptri_obj->n)*(dtptri_obj->n),
						   dtptri_obj->a, dtptri_obj->aref );
}

TEST_F(dtptri_test, dtptri1) {
    EXPECT_NEAR(0.0, dtptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtptri_test, dtptri2) {
    EXPECT_NEAR(0.0, dtptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtptri_test, dtptri3) {
    EXPECT_NEAR(0.0, dtptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtptri_test, dtptri4) {
    EXPECT_NEAR(0.0, dtptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin tptri_lapack_complex_float_parameters  class definition */
class tptri_lapack_complex_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tptri_lapack_complex_float_parameters ( int matrix_layout_i, char uplo_i, 
                        char diag_i, lapack_int n_i);
              
      ~tptri_lapack_complex_float_parameters (); 
};  /* end of tptri_lapack_complex_float_parameters  class definition */


/* Constructor tptri_lapack_complex_float_parameters definition */
tptri_lapack_complex_float_parameters:: tptri_lapack_complex_float_parameters ( int matrix_layout_i, 
            char uplo_i, char diag_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	diag = diag_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tptri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_lapack_complex_float_parameters object: malloc error.";
       tptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

tptri_lapack_complex_float_parameters:: ~tptri_lapack_complex_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tptri_lapack_complex_float_parameters object: destructor invoked. \n");
#endif
   tptri_free();
}

//  Test fixture class definition
class ctptri_test  : public  ::testing::Test {
public:
   tptri_lapack_complex_float_parameters  *ctptri_obj;
   void SetUp();  
   void TearDown () { delete ctptri_obj; }
};


void ctptri_test::SetUp(){

    /* LAPACKE CTPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_ctptri) (int matrix_layout, char uplo,
                   char diag, lapack_int n,  lapack_complex_float * a );

    Fptr_NL_LAPACKE_ctptri CTPTRI;

    ctptri_obj = new tptri_lapack_complex_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    ctptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CTPTRI = (Fptr_NL_LAPACKE_ctptri)dlsym(ctptri_obj->hModule, "LAPACKE_ctptri");
    ASSERT_TRUE(CTPTRI != NULL) << "failed to get the Netlib LAPACKE_ctptri symbol";
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    ctptri_obj->inforef = CTPTRI( 	ctptri_obj->matrix_layout,
									ctptri_obj->uplo, 
									ctptri_obj->diag, 
									ctptri_obj->n,
									ctptri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    ctptri_obj->info = LAPACKE_ctptri( 	ctptri_obj->matrix_layout, 
										ctptri_obj->uplo,
										ctptri_obj->diag, 
										ctptri_obj->n, 
										ctptri_obj->a);

    if( ctptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctptri is wrong\n", ctptri_obj->info );
    }
    if( ctptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctptri is wrong\n", 
        ctptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ctptri_obj->diff =  computeDiff_c( (ctptri_obj->n)*(ctptri_obj->n),
						   ctptri_obj->a, ctptri_obj->aref );
}

TEST_F(ctptri_test, ctptri1) {
    EXPECT_NEAR(0.0, ctptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctptri_test, ctptri2) {
    EXPECT_NEAR(0.0, ctptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctptri_test, ctptri3) {
    EXPECT_NEAR(0.0, ctptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctptri_test, ctptri4) {
    EXPECT_NEAR(0.0, ctptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin tptri_lapack_complex_double_parameters  class definition */
class tptri_lapack_complex_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B

      /* Input/ Output parameters */
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tptri_lapack_complex_double_parameters ( int matrix_layout_i, char uplo_i, 
                        char diag_i, lapack_int n_i);
              
      ~tptri_lapack_complex_double_parameters (); 
};  /* end of tptri_lapack_complex_double_parameters  class definition */


/* Constructor tptri_lapack_complex_double_parameters definition */
tptri_lapack_complex_double_parameters:: tptri_lapack_complex_double_parameters ( int matrix_layout_i, 
            char uplo_i, char diag_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	diag = diag_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tptri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_lapack_complex_double_parameters object: malloc error.";
       tptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

tptri_lapack_complex_double_parameters:: ~tptri_lapack_complex_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tptri_lapack_complex_double_parameters object: destructor invoked. \n");
#endif
   tptri_free();
}

//  Test fixture class definition
class ztptri_test  : public  ::testing::Test {
public:
   tptri_lapack_complex_double_parameters  *ztptri_obj;
   void SetUp();  
   void TearDown () { delete ztptri_obj; }
};


void ztptri_test::SetUp(){

    /* LAPACKE ZTPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_ztptri) (int matrix_layout, char uplo,
                   char diag, lapack_int n,  lapack_complex_double * a );

    Fptr_NL_LAPACKE_ztptri ZTPTRI;

    ztptri_obj = new tptri_lapack_complex_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    ztptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZTPTRI = (Fptr_NL_LAPACKE_ztptri)dlsym(ztptri_obj->hModule, "LAPACKE_ztptri");
    ASSERT_TRUE(ZTPTRI != NULL) << "failed to get the Netlib LAPACKE_ztptri symbol";
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    ztptri_obj->inforef = ZTPTRI( 	ztptri_obj->matrix_layout,
									ztptri_obj->uplo, 
									ztptri_obj->diag, 
									ztptri_obj->n,
									ztptri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    ztptri_obj->info = LAPACKE_ztptri( 	ztptri_obj->matrix_layout, 
										ztptri_obj->uplo,
										ztptri_obj->diag, 
										ztptri_obj->n, 
										ztptri_obj->a);

    if( ztptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztptri is wrong\n", ztptri_obj->info );
    }
    if( ztptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztptri is wrong\n", 
        ztptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ztptri_obj->diff =  computeDiff_z( (ztptri_obj->n)*(ztptri_obj->n),
						   ztptri_obj->a, ztptri_obj->aref );
}

TEST_F(ztptri_test, ztptri1) {
    EXPECT_NEAR(0.0, ztptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztptri_test, ztptri2) {
    EXPECT_NEAR(0.0, ztptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztptri_test, ztptri3) {
    EXPECT_NEAR(0.0, ztptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztptri_test, ztptri4) {
    EXPECT_NEAR(0.0, ztptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
