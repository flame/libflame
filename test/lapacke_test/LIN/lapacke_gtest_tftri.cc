#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define tftri_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin tftri_float_parameters  class definition */
class tftri_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	  char transr; //Must be 'N', 'T' (for real data) or 'C' (for complex data).
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B

      /* Input/ Output parameters */
      float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tftri_float_parameters ( int matrix_layout_i, char transr_i, 
						char uplo_i, char diag_i, lapack_int n_i);
              
      ~tftri_float_parameters (); 
};  /* end of tftri_float_parameters  class definition */


/* Constructor tftri_float_parameters definition */
tftri_float_parameters:: tftri_float_parameters ( int matrix_layout_i, 
            char transr_i, char uplo_i, char diag_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	transr = transr_i;
	diag = diag_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tftri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_float_parameters object: malloc error.";
       tftri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_float_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

tftri_float_parameters:: ~tftri_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tftri_float_parameters object: destructor invoked. \n");
#endif
   tftri_free();
}

//  Test fixture class definition
class stftri_test  : public  ::testing::Test {
public:
   tftri_float_parameters  *stftri_obj;
   void SetUp();  
   void TearDown () { delete stftri_obj; }
};


void stftri_test::SetUp(){

    /* LAPACKE STFTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_stftri) (int matrix_layout, char transr,
					char uplo, char diag, lapack_int n,  float * a );

    Fptr_NL_LAPACKE_stftri STFTRI;

    idx = Circular_Increment_Index(idx);
    stftri_obj = new tftri_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n);


    stftri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stftri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stftri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stftri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    STFTRI = (Fptr_NL_LAPACKE_stftri)dlsym(stftri_obj->hModule, "LAPACKE_stftri");
    ASSERT_TRUE(STFTRI != NULL) << "failed to get the Netlib LAPACKE_stftri symbol";
    
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    stftri_obj->inforef = STFTRI( 	stftri_obj->matrix_layout,
									stftri_obj->transr, 
									stftri_obj->uplo, 
									stftri_obj->diag, 
									stftri_obj->n,
									stftri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    stftri_obj->info = LAPACKE_stftri( 	stftri_obj->matrix_layout, 
										stftri_obj->transr, 
										stftri_obj->uplo,
										stftri_obj->diag, 
										stftri_obj->n, 
										stftri_obj->a);

    if( stftri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stftri is wrong\n", stftri_obj->info );
    }
    if( stftri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stftri is wrong\n", 
        stftri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    stftri_obj->diff =  computeDiff_s( (stftri_obj->n)*(stftri_obj->n),
						   stftri_obj->a, stftri_obj->aref );
}

TEST_F(stftri_test, stftri1) {
    EXPECT_NEAR(0.0, stftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stftri_test, stftri2) {
    EXPECT_NEAR(0.0, stftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stftri_test, stftri3) {
    EXPECT_NEAR(0.0, stftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stftri_test, stftri4) {
    EXPECT_NEAR(0.0, stftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin tftri_double_parameters  class definition */
class tftri_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	  char transr; //Must be 'N', 'T' (for real data) or 'C' (for complex data).
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B

      /* Input/ Output parameters */
      double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tftri_double_parameters ( int matrix_layout_i, char transr_i, 
							char uplo_i, char diag_i, lapack_int n_i);
              
      ~tftri_double_parameters (); 
};  /* end of tftri_double_parameters  class definition */


/* Constructor tftri_double_parameters definition */
tftri_double_parameters:: tftri_double_parameters ( int matrix_layout_i, 
            char transr_i, char uplo_i, char diag_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	transr = transr_i;
	diag = diag_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tftri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_double_parameters object: malloc error.";
       tftri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_double_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

tftri_double_parameters:: ~tftri_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tftri_double_parameters object: destructor invoked. \n");
#endif
   tftri_free();
}

//  Test fixture class definition
class dtftri_test  : public  ::testing::Test {
public:
   tftri_double_parameters  *dtftri_obj;
   void SetUp();  
   void TearDown () { delete dtftri_obj; }
};


void dtftri_test::SetUp(){

    /* LAPACKE DTFTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_dtftri) (int matrix_layout, char transr,
					char uplo, char diag, lapack_int n,  double * a );

    Fptr_NL_LAPACKE_dtftri DTFTRI;

    dtftri_obj = new tftri_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    dtftri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtftri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtftri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtftri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DTFTRI = (Fptr_NL_LAPACKE_dtftri)dlsym(dtftri_obj->hModule, "LAPACKE_dtftri");
    ASSERT_TRUE(DTFTRI != NULL) << "failed to get the Netlib LAPACKE_dtftri symbol";

    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    dtftri_obj->inforef = DTFTRI( 	dtftri_obj->matrix_layout,
									dtftri_obj->transr, 
									dtftri_obj->uplo, 
									dtftri_obj->diag, 
									dtftri_obj->n,
									dtftri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    dtftri_obj->info = LAPACKE_dtftri( 	dtftri_obj->matrix_layout, 
										dtftri_obj->transr, 
										dtftri_obj->uplo,
										dtftri_obj->diag, 
										dtftri_obj->n, 
										dtftri_obj->a);

    if( dtftri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtftri is wrong\n", dtftri_obj->info );
    }
    if( dtftri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtftri is wrong\n", 
        dtftri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    dtftri_obj->diff =  computeDiff_d( (dtftri_obj->n)*(dtftri_obj->n),
						   dtftri_obj->a, dtftri_obj->aref );
}

TEST_F(dtftri_test, dtftri1) {
    EXPECT_NEAR(0.0, dtftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtftri_test, dtftri2) {
    EXPECT_NEAR(0.0, dtftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtftri_test, dtftri3) {
    EXPECT_NEAR(0.0, dtftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtftri_test, dtftri4) {
    EXPECT_NEAR(0.0, dtftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin tftri_lapack_complex_float_parameters  class definition */
class tftri_lapack_complex_float_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	  char transr; //Must be 'N', 'T' (for real data) or 'C' (for complex data).
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tftri_lapack_complex_float_parameters ( int matrix_layout_i, char transr_i, 
									char uplo_i, char diag_i, lapack_int n_i);
              
      ~tftri_lapack_complex_float_parameters (); 
};  /* end of tftri_lapack_complex_float_parameters  class definition */


/* Constructor tftri_lapack_complex_float_parameters definition */
tftri_lapack_complex_float_parameters:: tftri_lapack_complex_float_parameters ( int matrix_layout_i, 
            char transr_i, char uplo_i, char diag_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	transr = transr_i;
	diag = diag_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tftri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_lapack_complex_float_parameters object: malloc error.";
       tftri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

tftri_lapack_complex_float_parameters:: ~tftri_lapack_complex_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tftri_lapack_complex_float_parameters object: destructor invoked. \n");
#endif
   tftri_free();
}

//  Test fixture class definition
class ctftri_test  : public  ::testing::Test {
public:
   tftri_lapack_complex_float_parameters  *ctftri_obj;
   void SetUp();  
   void TearDown () { delete ctftri_obj; }
};


void ctftri_test::SetUp(){

    /* LAPACKE CTFTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_ctftri) (int matrix_layout, char transr,
					char uplo, char diag, lapack_int n,  lapack_complex_float *a );

    Fptr_NL_LAPACKE_ctftri CTFTRI;

    ctftri_obj = new tftri_lapack_complex_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    ctftri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctftri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctftri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctftri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CTFTRI = (Fptr_NL_LAPACKE_ctftri)dlsym(ctftri_obj->hModule, "LAPACKE_ctftri");
    ASSERT_TRUE(CTFTRI != NULL) << "failed to get the Netlib LAPACKE_ctftri symbol";
    /* Compute the reference o/p by invoking Netlib-LapackE's API  */
    ctftri_obj->inforef = CTFTRI( 	ctftri_obj->matrix_layout,
									ctftri_obj->transr, 
									ctftri_obj->uplo, 
									ctftri_obj->diag, 
									ctftri_obj->n,
									ctftri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    ctftri_obj->info = LAPACKE_ctftri( 	ctftri_obj->matrix_layout, 
										ctftri_obj->transr, 
										ctftri_obj->uplo,
										ctftri_obj->diag, 
										ctftri_obj->n, 
										ctftri_obj->a);

    if( ctftri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctftri is wrong\n", ctftri_obj->info );
    }
    if( ctftri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctftri is wrong\n", 
        ctftri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ctftri_obj->diff =  computeDiff_c( (ctftri_obj->n)*(ctftri_obj->n),
						   ctftri_obj->a, ctftri_obj->aref );
}

TEST_F(ctftri_test, ctftri1) {
    EXPECT_NEAR(0.0, ctftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctftri_test, ctftri2) {
    EXPECT_NEAR(0.0, ctftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctftri_test, ctftri3) {
    EXPECT_NEAR(0.0, ctftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctftri_test, ctftri4) {
    EXPECT_NEAR(0.0, ctftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin tftri_lapack_complex_double_parameters  class definition */
class tftri_lapack_complex_double_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
	 
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
	  char transr; //Must be 'N', 'T' (for real data) or 'C' (for complex data).
      char uplo; //  Must be 'U' or 'L'.
      char diag; //  Must be 'U' or 'N'.
      lapack_int n; // The order of A; the number of rows in B

      /* Input/ Output parameters */
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      tftri_lapack_complex_double_parameters ( int matrix_layout_i, char transr_i, 
								char uplo_i, char diag_i, lapack_int n_i);
              
      ~tftri_lapack_complex_double_parameters (); 
};  /* end of tftri_lapack_complex_double_parameters  class definition */


/* Constructor tftri_lapack_complex_double_parameters definition */
tftri_lapack_complex_double_parameters:: tftri_lapack_complex_double_parameters ( int matrix_layout_i, 
            char transr_i, char uplo_i, char diag_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
	transr = transr_i;
	diag = diag_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n tftri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_lapack_complex_double_parameters object: malloc error.";
       tftri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

tftri_lapack_complex_double_parameters:: ~tftri_lapack_complex_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" tftri_lapack_complex_double_parameters object: destructor invoked. \n");
#endif
   tftri_free();
}

//  Test fixture class definition
class ztftri_test  : public  ::testing::Test {
public:
   tftri_lapack_complex_double_parameters  *ztftri_obj;
   void SetUp();  
   void TearDown () { delete ztftri_obj; }
};


void ztftri_test::SetUp(){

    /* LAPACKE ZTFTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_ztftri) (int matrix_layout, char transr,
					char uplo, char diag, lapack_int n,  lapack_complex_double * a );

    Fptr_NL_LAPACKE_ztftri ZTFTRI;

    ztftri_obj = new tftri_lapack_complex_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].diag,
                           lin_solver_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    ztftri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztftri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztftri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztftri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZTFTRI = (Fptr_NL_LAPACKE_ztftri)dlsym(ztftri_obj->hModule, "LAPACKE_ztftri");
    ASSERT_TRUE(ZTFTRI != NULL) << "failed to get the Netlib LAPACKE_ztftri symbol";
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    ztftri_obj->inforef = ZTFTRI( 	ztftri_obj->matrix_layout,
									ztftri_obj->transr, 
									ztftri_obj->uplo,
									ztftri_obj->diag, 
									ztftri_obj->n,
									ztftri_obj->aref);
                          

    /* Compute libflame's Lapacke o/p  */
    ztftri_obj->info = LAPACKE_ztftri( 	ztftri_obj->matrix_layout, 
										ztftri_obj->transr, 
										ztftri_obj->uplo,
										ztftri_obj->diag, 
										ztftri_obj->n, 
										ztftri_obj->a);

    if( ztftri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztftri is wrong\n", ztftri_obj->info );
    }
    if( ztftri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztftri is wrong\n", 
        ztftri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    ztftri_obj->diff =  computeDiff_z( (ztftri_obj->n)*(ztftri_obj->n),
						   ztftri_obj->a, ztftri_obj->aref );
}

TEST_F(ztftri_test, ztftri1) {
    EXPECT_NEAR(0.0, ztftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztftri_test, ztftri2) {
    EXPECT_NEAR(0.0, ztftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztftri_test, ztftri3) {
    EXPECT_NEAR(0.0, ztftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztftri_test, ztftri4) {
    EXPECT_NEAR(0.0, ztftri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
