#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define hetri2_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (ipiv != NULL) free (ipiv); \
  if (ipivref != NULL)free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin hetri2_scomplex_parameters  class definition */
class hetri2_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      hetri2_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~hetri2_scomplex_parameters (); 
};  /* end of hetri2_scomplex_parameters  class definition */


/* Constructor hetri2_scomplex_parameters definition */
hetri2_scomplex_parameters:: hetri2_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n hetri2 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       hetri2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    //lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(a, aref, n, n, 'H');
    
   } /* end of Constructor  */

hetri2_scomplex_parameters:: ~hetri2_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetri2_scomplex_parameters object: destructor invoked. \n");
#endif
   hetri2_free();
}

//  Test fixture class definition
class chetri2_test  : public  ::testing::Test {
public:
   hetri2_scomplex_parameters  *chetri2_obj;
   void SetUp();  
   void TearDown () { delete chetri2_obj; }
};


void chetri2_test::SetUp(){

    /* LAPACKE CHETRI2 prototype */
    typedef int (*Fptr_NL_LAPACKE_chetri2) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_float * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_chetri2 CHETRI2;

    typedef int (*Fptr_NL_LAPACKE_chetrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_chetrf CHETRF;

    chetri2_obj = new hetri2_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    chetri2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetri2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetri2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetri2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CHETRI2 = (Fptr_NL_LAPACKE_chetri2)dlsym(chetri2_obj->hModule, "LAPACKE_chetri2");
    ASSERT_TRUE(CHETRI2 != NULL) << "failed to get the Netlib LAPACKE_chetri2 symbol";
    
    CHETRF = (Fptr_NL_LAPACKE_chetrf)dlsym(chetri2_obj->hModule,"LAPACKE_chetrf");
    ASSERT_TRUE(CHETRF != NULL) << "failed to get the Netlib LAPACKE_chetrf symbol";

    /* Pre condition: need to call hetrf - before calling hetri2 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    chetri2_obj->inforef = CHETRF(   chetri2_obj->matrix_layout,
                                    chetri2_obj->uplo,
                                    chetri2_obj->n,
                                    chetri2_obj->aref,
                                    chetri2_obj->lda,
                                    chetri2_obj->ipivref);

    chetri2_obj->inforef = CHETRI2(   chetri2_obj->matrix_layout,
                                    chetri2_obj->uplo, 
                                    chetri2_obj->n,
                                    chetri2_obj->aref,
                                    chetri2_obj->lda,
                                    (const lapack_int *)chetri2_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    chetri2_obj->info  = LAPACKE_chetrf( chetri2_obj->matrix_layout,
                                        chetri2_obj->uplo,
                                        chetri2_obj->n,
                                        chetri2_obj->a,
                                        chetri2_obj->lda,
                                        chetri2_obj->ipiv);

    chetri2_obj->info = LAPACKE_chetri2(  chetri2_obj->matrix_layout, 
                                        chetri2_obj->uplo,
                                        chetri2_obj->n, 
                                        chetri2_obj->a, 
                                        chetri2_obj->lda,
                                        (const lapack_int *)chetri2_obj->ipiv);

    if( chetri2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chetri2 is wrong\n", chetri2_obj->info );
    }
    if( chetri2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetri2 is wrong\n", 
        chetri2_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    chetri2_obj->diff =  computeDiff_c( (chetri2_obj->n)*(chetri2_obj->lda),
                           chetri2_obj->a, chetri2_obj->aref );
}

TEST_F(chetri2_test, chetri21) {
    EXPECT_NEAR(0.0, chetri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri2_test, chetri22) {
    EXPECT_NEAR(0.0, chetri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri2_test, chetri23) {
    EXPECT_NEAR(0.0, chetri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri2_test, chetri24) {
    EXPECT_NEAR(0.0, chetri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin hetri2_dcomplex_parameters  class definition */
class hetri2_dcomplex_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Input/ Output parameters */
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      hetri2_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~hetri2_dcomplex_parameters (); 
};  /* end of hetri2_dcomplex_parameters  class definition */


/* Constructor hetri2_dcomplex_parameters definition */
hetri2_dcomplex_parameters:: hetri2_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n hetri2 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       hetri2_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    //lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(a, aref, n, n, 'H');
    
   } /* end of Constructor  */

hetri2_dcomplex_parameters:: ~hetri2_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetri2_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetri2_free();
}

//  Test fixture class definition
class zhetri2_test  : public  ::testing::Test {
public:
   hetri2_dcomplex_parameters  *zhetri2_obj;
   void SetUp();  
   void TearDown () { delete zhetri2_obj; }
};


void zhetri2_test::SetUp(){

    /* LAPACKE ZHETRI2 prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetri2) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_double * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_zhetri2 ZHETRI2;

    typedef int (*Fptr_NL_LAPACKE_zhetrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_zhetrf ZHETRF;

    zhetri2_obj = new hetri2_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    zhetri2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetri2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetri2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetri2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZHETRI2 = (Fptr_NL_LAPACKE_zhetri2)dlsym(zhetri2_obj->hModule, "LAPACKE_zhetri2");
    ASSERT_TRUE(ZHETRI2 != NULL) << "failed to get the Netlib LAPACKE_zhetri2 symbol";
    
    ZHETRF = (Fptr_NL_LAPACKE_zhetrf)dlsym(zhetri2_obj->hModule,"LAPACKE_zhetrf");
    ASSERT_TRUE(ZHETRF != NULL) << "failed to get the Netlib LAPACKE_zhetrf symbol";

    /* Pre condition: need to call hetrf - before calling hetri2 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zhetri2_obj->inforef = ZHETRF(   zhetri2_obj->matrix_layout,
                                    zhetri2_obj->uplo,
                                    zhetri2_obj->n,
                                    zhetri2_obj->aref,
                                    zhetri2_obj->lda,
                                    zhetri2_obj->ipivref);

    zhetri2_obj->inforef = ZHETRI2(   zhetri2_obj->matrix_layout,
                                    zhetri2_obj->uplo, 
                                    zhetri2_obj->n,
                                    zhetri2_obj->aref,
                                    zhetri2_obj->lda,
                                    zhetri2_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    zhetri2_obj->info  = LAPACKE_zhetrf( zhetri2_obj->matrix_layout,
                                        zhetri2_obj->uplo,
                                        zhetri2_obj->n,
                                        zhetri2_obj->a,
                                        zhetri2_obj->lda,
                                        zhetri2_obj->ipiv);

    zhetri2_obj->info = LAPACKE_zhetri2(  zhetri2_obj->matrix_layout, 
                                        zhetri2_obj->uplo,
                                        zhetri2_obj->n, 
                                        zhetri2_obj->a, 
                                        zhetri2_obj->lda,
                                        zhetri2_obj->ipiv);

    if( zhetri2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhetri2 is wrong\n", zhetri2_obj->info );
    }
    if( zhetri2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetri2 is wrong\n", 
        zhetri2_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zhetri2_obj->diff =  computeDiff_z( (zhetri2_obj->n)*(zhetri2_obj->lda),
                           zhetri2_obj->a, zhetri2_obj->aref );
}

TEST_F(zhetri2_test, zhetri21) {
    EXPECT_NEAR(0.0, zhetri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri2_test, zhetri22) {
    EXPECT_NEAR(0.0, zhetri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri2_test, zhetri23) {
    EXPECT_NEAR(0.0, zhetri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri2_test, zhetri24) {
    EXPECT_NEAR(0.0, zhetri2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

