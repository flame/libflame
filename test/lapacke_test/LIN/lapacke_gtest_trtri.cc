#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define trtri_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin trtri_float_parameters  class definition */
class trtri_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      trtri_float_parameters ( int matrix_layout_i, char uplo_i, 
                        char diag_i, lapack_int n_i, lapack_int lda_i);
              
      ~trtri_float_parameters (); 
};  /* end of trtri_float_parameters  class definition */


/* Constructor trtri_float_parameters definition */
trtri_float_parameters:: trtri_float_parameters ( int matrix_layout_i, 
            char uplo_i, char diag_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	diag = diag_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n trtri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       trtri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

trtri_float_parameters:: ~trtri_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trtri_float_parameters object: destructor invoked. \n");
#endif
   trtri_free();
}

//  Test fixture class definition
class strtri_test  : public  ::testing::Test {
public:
   trtri_float_parameters  *strtri_obj;
   void SetUp();  
   void TearDown () { delete strtri_obj; }
};


void strtri_test::SetUp(){

    /* LAPACKE STRTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_strtri) (int matrix_layout, char uplo,
                   char diag, lapack_int n,  float * a, lapack_int lda  );

    Fptr_NL_LAPACKE_strtri STRTRI;

    idx = Circular_Increment_Index(idx);
    strtri_obj = new trtri_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );


    strtri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    strtri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(strtri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(strtri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    STRTRI = (Fptr_NL_LAPACKE_strtri)dlsym(strtri_obj->hModule, "LAPACKE_strtri");
    ASSERT_TRUE(STRTRI != NULL) << "failed to get the Netlib LAPACKE_strtri symbol";
    
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    strtri_obj->inforef = STRTRI( 	strtri_obj->matrix_layout,
									strtri_obj->uplo, 
									strtri_obj->diag, 
									strtri_obj->n,
									strtri_obj->aref,
									strtri_obj->lda);
                          

    /* Compute libflame's Lapacke o/p  */
    strtri_obj->info = LAPACKE_strtri( 	strtri_obj->matrix_layout, 
										strtri_obj->uplo,
										strtri_obj->diag, 
										strtri_obj->n, 
										strtri_obj->a, 
										strtri_obj->lda);

    if( strtri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_strtri is wrong\n", strtri_obj->info );
    }
    if( strtri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_strtri is wrong\n", 
        strtri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    strtri_obj->diff =  computeDiff_s( (strtri_obj->n)*(strtri_obj->lda),
						   strtri_obj->a, strtri_obj->aref );
}

TEST_F(strtri_test, strtri1) {
    EXPECT_NEAR(0.0, strtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strtri_test, strtri2) {
    EXPECT_NEAR(0.0, strtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strtri_test, strtri3) {
    EXPECT_NEAR(0.0, strtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(strtri_test, strtri4) {
    EXPECT_NEAR(0.0, strtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin trtri_double_parameters  class definition */
class trtri_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      trtri_double_parameters ( int matrix_layout_i, char uplo_i, 
                        char diag_i, lapack_int n_i, lapack_int lda_i);
              
      ~trtri_double_parameters (); 
};  /* end of trtri_double_parameters  class definition */


/* Constructor trtri_double_parameters definition */
trtri_double_parameters:: trtri_double_parameters ( int matrix_layout_i, 
            char uplo_i, char diag_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	diag = diag_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n trtri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       trtri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

trtri_double_parameters:: ~trtri_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trtri_double_parameters object: destructor invoked. \n");
#endif
   trtri_free();
}

//  Test fixture class definition
class dtrtri_test  : public  ::testing::Test {
public:
   trtri_double_parameters  *dtrtri_obj;
   void SetUp();  
   void TearDown () { delete dtrtri_obj; }
};


void dtrtri_test::SetUp(){

    /* LAPACKE DTRTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_dtrtri) (int matrix_layout, char uplo,
                   char diag, lapack_int n,  double * a, lapack_int lda  );

    Fptr_NL_LAPACKE_dtrtri DTRTRI;

    dtrtri_obj = new trtri_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    dtrtri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtrtri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtrtri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtrtri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DTRTRI = (Fptr_NL_LAPACKE_dtrtri)dlsym(dtrtri_obj->hModule, "LAPACKE_dtrtri");
    ASSERT_TRUE(DTRTRI != NULL) << "failed to get the Netlib LAPACKE_dtrtri symbol";

    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dtrtri_obj->inforef = DTRTRI( 	dtrtri_obj->matrix_layout,
									dtrtri_obj->uplo, 
									dtrtri_obj->diag, 
									dtrtri_obj->n,
									dtrtri_obj->aref,
									dtrtri_obj->lda);
                          

    /* Compute libflame's Lapacke o/p  */
    dtrtri_obj->info = LAPACKE_dtrtri( 	dtrtri_obj->matrix_layout, 
										dtrtri_obj->uplo,
										dtrtri_obj->diag, 
										dtrtri_obj->n, 
										dtrtri_obj->a, 
										dtrtri_obj->lda);

    if( dtrtri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtrtri is wrong\n", dtrtri_obj->info );
    }
    if( dtrtri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtrtri is wrong\n", 
        dtrtri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dtrtri_obj->diff =  computeDiff_d( (dtrtri_obj->n)*(dtrtri_obj->lda),
						   dtrtri_obj->a, dtrtri_obj->aref );
}

TEST_F(dtrtri_test, dtrtri1) {
    EXPECT_NEAR(0.0, dtrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrtri_test, dtrtri2) {
    EXPECT_NEAR(0.0, dtrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrtri_test, dtrtri3) {
    EXPECT_NEAR(0.0, dtrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtrtri_test, dtrtri4) {
    EXPECT_NEAR(0.0, dtrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin trtri_lapack_complex_float_parameters  class definition */
class trtri_lapack_complex_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      trtri_lapack_complex_float_parameters ( int matrix_layout_i, char uplo_i, 
                        char diag_i, lapack_int n_i, lapack_int lda_i);
              
      ~trtri_lapack_complex_float_parameters (); 
};  /* end of trtri_lapack_complex_float_parameters  class definition */


/* Constructor trtri_lapack_complex_float_parameters definition */
trtri_lapack_complex_float_parameters:: trtri_lapack_complex_float_parameters ( int matrix_layout_i, 
            char uplo_i, char diag_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	diag = diag_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n trtri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_lapack_complex_float_parameters object: malloc error.";
       trtri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

trtri_lapack_complex_float_parameters:: ~trtri_lapack_complex_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trtri_lapack_complex_float_parameters object: destructor invoked. \n");
#endif
   trtri_free();
}

//  Test fixture class definition
class ctrtri_test  : public  ::testing::Test {
public:
   trtri_lapack_complex_float_parameters  *ctrtri_obj;
   void SetUp();  
   void TearDown () { delete ctrtri_obj; }
};


void ctrtri_test::SetUp(){

    /* LAPACKE CTRTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_ctrtri) (int matrix_layout, char uplo,
                   char diag, lapack_int n,  lapack_complex_float * a, lapack_int lda  );

    Fptr_NL_LAPACKE_ctrtri CTRTRI;

    ctrtri_obj = new trtri_lapack_complex_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    ctrtri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctrtri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctrtri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctrtri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CTRTRI = (Fptr_NL_LAPACKE_ctrtri)dlsym(ctrtri_obj->hModule, "LAPACKE_ctrtri");
    ASSERT_TRUE(CTRTRI != NULL) << "failed to get the Netlib LAPACKE_ctrtri symbol";
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    ctrtri_obj->inforef = CTRTRI( 	ctrtri_obj->matrix_layout,
									ctrtri_obj->uplo, 
									ctrtri_obj->diag, 
									ctrtri_obj->n,
									ctrtri_obj->aref,
									ctrtri_obj->lda);
                          

    /* Compute libflame's Lapacke o/p  */
    ctrtri_obj->info = LAPACKE_ctrtri( 	ctrtri_obj->matrix_layout, 
										ctrtri_obj->uplo,
										ctrtri_obj->diag, 
										ctrtri_obj->n, 
										ctrtri_obj->a, 
										ctrtri_obj->lda);

    if( ctrtri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctrtri is wrong\n", ctrtri_obj->info );
    }
    if( ctrtri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctrtri is wrong\n", 
        ctrtri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ctrtri_obj->diff =  computeDiff_c( (ctrtri_obj->n)*(ctrtri_obj->lda),
						   ctrtri_obj->a, ctrtri_obj->aref );
}

TEST_F(ctrtri_test, ctrtri1) {
    EXPECT_NEAR(0.0, ctrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrtri_test, ctrtri2) {
    EXPECT_NEAR(0.0, ctrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrtri_test, ctrtri3) {
    EXPECT_NEAR(0.0, ctrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctrtri_test, ctrtri4) {
    EXPECT_NEAR(0.0, ctrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin trtri_lapack_complex_double_parameters  class definition */
class trtri_lapack_complex_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'

      /* Input/ Output parameters */
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      trtri_lapack_complex_double_parameters ( int matrix_layout_i, char uplo_i, 
                        char diag_i, lapack_int n_i, lapack_int lda_i);
              
      ~trtri_lapack_complex_double_parameters (); 
};  /* end of trtri_lapack_complex_double_parameters  class definition */


/* Constructor trtri_lapack_complex_double_parameters definition */
trtri_lapack_complex_double_parameters:: trtri_lapack_complex_double_parameters ( int matrix_layout_i, 
            char uplo_i, char diag_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	diag = diag_i;
    lda = lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n trtri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_lapack_complex_double_parameters object: malloc error.";
       trtri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

trtri_lapack_complex_double_parameters:: ~trtri_lapack_complex_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" trtri_lapack_complex_double_parameters object: destructor invoked. \n");
#endif
   trtri_free();
}

//  Test fixture class definition
class ztrtri_test  : public  ::testing::Test {
public:
   trtri_lapack_complex_double_parameters  *ztrtri_obj;
   void SetUp();  
   void TearDown () { delete ztrtri_obj; }
};


void ztrtri_test::SetUp(){

    /* LAPACKE ZTRTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_ztrtri) (int matrix_layout, char uplo,
                   char diag, lapack_int n,  lapack_complex_double * a, lapack_int lda  );

    Fptr_NL_LAPACKE_ztrtri ZTRTRI;

    ztrtri_obj = new trtri_lapack_complex_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    ztrtri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztrtri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztrtri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztrtri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZTRTRI = (Fptr_NL_LAPACKE_ztrtri)dlsym(ztrtri_obj->hModule, "LAPACKE_ztrtri");
    ASSERT_TRUE(ZTRTRI != NULL) << "failed to get the Netlib LAPACKE_ztrtri symbol";
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    ztrtri_obj->inforef = ZTRTRI( 	ztrtri_obj->matrix_layout,
									ztrtri_obj->uplo, 
									ztrtri_obj->diag, 
									ztrtri_obj->n,
									ztrtri_obj->aref,
									ztrtri_obj->lda);
                          

    /* Compute libflame's Lapacke o/p  */
    ztrtri_obj->info = LAPACKE_ztrtri( 	ztrtri_obj->matrix_layout, 
										ztrtri_obj->uplo,
										ztrtri_obj->diag, 
										ztrtri_obj->n, 
										ztrtri_obj->a, 
										ztrtri_obj->lda);

    if( ztrtri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztrtri is wrong\n", ztrtri_obj->info );
    }
    if( ztrtri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztrtri is wrong\n", 
        ztrtri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ztrtri_obj->diff =  computeDiff_z( (ztrtri_obj->n)*(ztrtri_obj->lda),
						   ztrtri_obj->a, ztrtri_obj->aref );
}

TEST_F(ztrtri_test, ztrtri1) {
    EXPECT_NEAR(0.0, ztrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrtri_test, ztrtri2) {
    EXPECT_NEAR(0.0, ztrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrtri_test, ztrtri3) {
    EXPECT_NEAR(0.0, ztrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztrtri_test, ztrtri4) {
    EXPECT_NEAR(0.0, ztrtri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
