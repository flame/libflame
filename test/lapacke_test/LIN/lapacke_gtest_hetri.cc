#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define hetri_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (ipiv != NULL) free (ipiv); \
  if (ipivref != NULL)free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin hetri_scomplex_parameters  class definition */
class hetri_scomplex_parameters{
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
      hetri_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~hetri_scomplex_parameters (); 
};  /* end of hetri_scomplex_parameters  class definition */


/* Constructor hetri_scomplex_parameters definition */
hetri_scomplex_parameters:: hetri_scomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = n;//lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n hetri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       hetri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    //lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, n,n, 'H');
    
   } /* end of Constructor  */

hetri_scomplex_parameters:: ~hetri_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetri_scomplex_parameters object: destructor invoked. \n");
#endif
   hetri_free();
}

//  Test fixture class definition
class chetri_test  : public  ::testing::Test {
public:
   hetri_scomplex_parameters  *chetri_obj;
   void SetUp();  
   void TearDown () { delete chetri_obj; }
};


void chetri_test::SetUp(){

    /* LAPACKE CHETRI prototype */
    typedef int (*Fptr_NL_LAPACKE_chetri) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_float * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_chetri CHETRI;

    typedef int (*Fptr_NL_LAPACKE_chetrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_chetrf CHETRF;

    chetri_obj = new hetri_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    chetri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CHETRI = (Fptr_NL_LAPACKE_chetri)dlsym(chetri_obj->hModule, "LAPACKE_chetri");
    ASSERT_TRUE(CHETRI != NULL) << "failed to get the Netlib LAPACKE_chetri symbol";
    
    CHETRF = (Fptr_NL_LAPACKE_chetrf)dlsym(chetri_obj->hModule,"LAPACKE_chetrf");
    ASSERT_TRUE(CHETRF != NULL) << "failed to get the Netlib LAPACKE_chetrf symbol";

    /* Pre condition: need to call hetrf - before calling hetri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    chetri_obj->inforef = CHETRF(   chetri_obj->matrix_layout,
                                    chetri_obj->uplo,
                                    chetri_obj->n,
                                    chetri_obj->aref,
                                    chetri_obj->lda,
                                    chetri_obj->ipivref);

    chetri_obj->inforef = CHETRI(   chetri_obj->matrix_layout,
                                    chetri_obj->uplo, 
                                    chetri_obj->n,
                                    chetri_obj->aref,
                                    chetri_obj->lda,
                                    (const lapack_int *)chetri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    chetri_obj->info  = LAPACKE_chetrf( chetri_obj->matrix_layout,
                                        chetri_obj->uplo,
                                        chetri_obj->n,
                                        chetri_obj->a,
                                        chetri_obj->lda,
                                        chetri_obj->ipiv);

    chetri_obj->info = LAPACKE_chetri(  chetri_obj->matrix_layout, 
                                        chetri_obj->uplo,
                                        chetri_obj->n, 
                                        chetri_obj->a, 
                                        chetri_obj->lda,
                                        (const lapack_int *)chetri_obj->ipiv);

    if( chetri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chetri is wrong\n", chetri_obj->info );
    }
    if( chetri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetri is wrong\n", 
        chetri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    chetri_obj->diff =  computeDiff_c( (chetri_obj->n)*(chetri_obj->lda),
                           chetri_obj->a, chetri_obj->aref );
}

TEST_F(chetri_test, chetri1) {
    EXPECT_NEAR(0.0, chetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri_test, chetri2) {
    EXPECT_NEAR(0.0, chetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri_test, chetri3) {
    EXPECT_NEAR(0.0, chetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri_test, chetri4) {
    EXPECT_NEAR(0.0, chetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin hetri_dcomplex_parameters  class definition */
class hetri_dcomplex_parameters{
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
      hetri_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~hetri_dcomplex_parameters (); 
};  /* end of hetri_dcomplex_parameters  class definition */


/* Constructor hetri_dcomplex_parameters definition */
hetri_dcomplex_parameters:: hetri_dcomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = n;//lda_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n hetri Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       hetri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

hetri_dcomplex_parameters:: ~hetri_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetri_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetri_free();
}

//  Test fixture class definition
class zhetri_test  : public  ::testing::Test {
public:
   hetri_dcomplex_parameters  *zhetri_obj;
   void SetUp();  
   void TearDown () { delete zhetri_obj; }
};


void zhetri_test::SetUp(){

    /* LAPACKE ZHETRI prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetri) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_double * a, lapack_int lda,
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_zhetri ZHETRI;

    typedef int (*Fptr_NL_LAPACKE_zhetrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_zhetrf ZHETRF;

    zhetri_obj = new hetri_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    zhetri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZHETRI = (Fptr_NL_LAPACKE_zhetri)dlsym(zhetri_obj->hModule, "LAPACKE_zhetri");
    ASSERT_TRUE(ZHETRI != NULL) << "failed to get the Netlib LAPACKE_zhetri symbol";
    
    ZHETRF = (Fptr_NL_LAPACKE_zhetrf)dlsym(zhetri_obj->hModule,"LAPACKE_zhetrf");
    ASSERT_TRUE(ZHETRF != NULL) << "failed to get the Netlib LAPACKE_zhetrf symbol";

    /* Pre condition: need to call hetrf - before calling hetri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zhetri_obj->inforef = ZHETRF(   zhetri_obj->matrix_layout,
                                    zhetri_obj->uplo,
                                    zhetri_obj->n,
                                    zhetri_obj->aref,
                                    zhetri_obj->lda,
                                    zhetri_obj->ipivref);

    zhetri_obj->inforef = ZHETRI(   zhetri_obj->matrix_layout,
                                    zhetri_obj->uplo, 
                                    zhetri_obj->n,
                                    zhetri_obj->aref,
                                    zhetri_obj->lda,
                                    zhetri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    zhetri_obj->info  = LAPACKE_zhetrf( zhetri_obj->matrix_layout,
                                        zhetri_obj->uplo,
                                        zhetri_obj->n,
                                        zhetri_obj->a,
                                        zhetri_obj->lda,
                                        zhetri_obj->ipiv);

    zhetri_obj->info = LAPACKE_zhetri(  zhetri_obj->matrix_layout, 
                                        zhetri_obj->uplo,
                                        zhetri_obj->n, 
                                        zhetri_obj->a, 
                                        zhetri_obj->lda,
                                        zhetri_obj->ipiv);

    if( zhetri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhetri is wrong\n", zhetri_obj->info );
    }
    if( zhetri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetri is wrong\n", 
        zhetri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zhetri_obj->diff =  computeDiff_z( (zhetri_obj->n)*(zhetri_obj->lda),
                           zhetri_obj->a, zhetri_obj->aref );
}

TEST_F(zhetri_test, zhetri1) {
    EXPECT_NEAR(0.0, zhetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri_test, zhetri2) {
    EXPECT_NEAR(0.0, zhetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri_test, zhetri3) {
    EXPECT_NEAR(0.0, zhetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri_test, zhetri4) {
    EXPECT_NEAR(0.0, zhetri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

