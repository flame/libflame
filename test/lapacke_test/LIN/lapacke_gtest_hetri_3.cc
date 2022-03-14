#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define hetri_3_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (e != NULL)    free (e   ); \
  if (eref != NULL) free (eref); \
  if (ipiv != NULL) free (ipiv); \
  if (ipivref != NULL)free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin hetri_3_scomplex_parameters  class definition */
class hetri_3_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_complex_float *e, *eref; // superdiagonal (or subdiagonal) elements of
      //  the symmetric block diagonal matrix D


      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      hetri_3_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~hetri_3_scomplex_parameters (); 
};  /* end of hetri_3_scomplex_parameters  class definition */


/* Constructor hetri_3_scomplex_parameters definition */
hetri_3_scomplex_parameters:: hetri_3_scomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n hetri_3 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       hetri_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    //lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(a, aref, n, n, 'H');
    
   } /* end of Constructor  */

hetri_3_scomplex_parameters:: ~hetri_3_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetri_3_scomplex_parameters object: destructor invoked. \n");
#endif
   hetri_3_free();
}

//  Test fixture class definition
class chetri_3_test  : public  ::testing::Test {
public:
   hetri_3_scomplex_parameters  *chetri_3_obj;
   void SetUp();  
   void TearDown () { delete chetri_3_obj; }
};


void chetri_3_test::SetUp(){

    /* LAPACKE CHETRI_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_chetri_3) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_float * a, lapack_int lda,
                                const lapack_complex_float * e, const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_chetri_3 CHETRI_3;

    typedef int (*Fptr_NL_LAPACKE_chetrf_rk) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda,
                               lapack_complex_float * e, lapack_int *ipiv );

    Fptr_NL_LAPACKE_chetrf_rk CHETRF_RK;

    chetri_3_obj = new hetri_3_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    chetri_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetri_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetri_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetri_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CHETRI_3 = (Fptr_NL_LAPACKE_chetri_3)dlsym(chetri_3_obj->hModule, "LAPACKE_chetri_3");
    ASSERT_TRUE(CHETRI_3 != NULL) << "failed to get the Netlib LAPACKE_chetri_3 symbol";
    
    CHETRF_RK = (Fptr_NL_LAPACKE_chetrf_rk)dlsym(chetri_3_obj->hModule,"LAPACKE_chetrf_rk");
    ASSERT_TRUE(CHETRF_RK != NULL) << "failed to get the Netlib LAPACKE_chetrf_rk symbol";

    /* Pre condition: need to call hetrf_rk - before calling hetri_3 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    chetri_3_obj->inforef = CHETRF_RK(   chetri_3_obj->matrix_layout,
                                    chetri_3_obj->uplo,
                                    chetri_3_obj->n,
                                    chetri_3_obj->aref,
                                    chetri_3_obj->lda,
                                    chetri_3_obj->eref,
                                    chetri_3_obj->ipivref);

    chetri_3_obj->inforef = CHETRI_3(   chetri_3_obj->matrix_layout,
                                    chetri_3_obj->uplo, 
                                    chetri_3_obj->n,
                                    chetri_3_obj->aref,
                                    chetri_3_obj->lda,
                             (const lapack_complex_float *)chetri_3_obj->eref,
                                    (const lapack_int *)chetri_3_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    chetri_3_obj->info  = LAPACKE_chetrf_rk( chetri_3_obj->matrix_layout,
                                        chetri_3_obj->uplo,
                                        chetri_3_obj->n,
                                        chetri_3_obj->a,
                                        chetri_3_obj->lda,
										chetri_3_obj->e,
                                        chetri_3_obj->ipiv);

    chetri_3_obj->info = LAPACKE_chetri_3(  chetri_3_obj->matrix_layout, 
                                        chetri_3_obj->uplo,
                                        chetri_3_obj->n, 
                                        chetri_3_obj->a, 
                                        chetri_3_obj->lda,
                             (const lapack_complex_float *)chetri_3_obj->e,
                                        (const lapack_int *)chetri_3_obj->ipiv);

    if( chetri_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chetri_3 is wrong\n", chetri_3_obj->info );
    }
    if( chetri_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetri_3 is wrong\n", 
        chetri_3_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    chetri_3_obj->diff =  computeDiff_c( (chetri_3_obj->n)*(chetri_3_obj->lda),
                           chetri_3_obj->a, chetri_3_obj->aref );
}

TEST_F(chetri_3_test, chetri_31) {
    EXPECT_NEAR(0.0, chetri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri_3_test, chetri_32) {
    EXPECT_NEAR(0.0, chetri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri_3_test, chetri_33) {
    EXPECT_NEAR(0.0, chetri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri_3_test, chetri_34) {
    EXPECT_NEAR(0.0, chetri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin hetri_3_dcomplex_parameters  class definition */
class hetri_3_dcomplex_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_complex_double *e, *eref; // superdiagonal (or subdiagonal) elements of
      //  the symmetric block diagonal matrix D


      /* Input/ Output parameters */
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      hetri_3_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i);
              
      ~hetri_3_dcomplex_parameters (); 
};  /* end of hetri_3_dcomplex_parameters  class definition */


/* Constructor hetri_3_dcomplex_parameters definition */
hetri_3_dcomplex_parameters:: hetri_3_dcomplex_parameters ( int matrix_layout_i, 
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
   printf(" \n hetri_3 Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &e, &eref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       hetri_3_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    //lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(a, aref, n, n, 'H');
    
   } /* end of Constructor  */

hetri_3_dcomplex_parameters:: ~hetri_3_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetri_3_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetri_3_free();
}

//  Test fixture class definition
class zhetri_3_test  : public  ::testing::Test {
public:
   hetri_3_dcomplex_parameters  *zhetri_3_obj;
   void SetUp();  
   void TearDown () { delete zhetri_3_obj; }
};


void zhetri_3_test::SetUp(){

    /* LAPACKE ZHETRI_3 prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetri_3) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_double * a, lapack_int lda,
                                const lapack_complex_double * e, const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_zhetri_3 ZHETRI_3;

    typedef int (*Fptr_NL_LAPACKE_zhetrf_rk) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda,
                               lapack_complex_double * e, lapack_int *ipiv );

    Fptr_NL_LAPACKE_zhetrf_rk ZHETRF_RK;

    zhetri_3_obj = new hetri_3_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda );

    idx = Circular_Increment_Index(idx);

    zhetri_3_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetri_3_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetri_3_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetri_3_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZHETRI_3 = (Fptr_NL_LAPACKE_zhetri_3)dlsym(zhetri_3_obj->hModule, "LAPACKE_zhetri_3");
    ASSERT_TRUE(ZHETRI_3 != NULL) << "failed to get the Netlib LAPACKE_zhetri_3 symbol";
    
    ZHETRF_RK = (Fptr_NL_LAPACKE_zhetrf_rk)dlsym(zhetri_3_obj->hModule,"LAPACKE_zhetrf_rk");
    ASSERT_TRUE(ZHETRF_RK != NULL) << "failed to get the Netlib LAPACKE_zhetrf_rk symbol";

    /* Pre condition: need to call hetrf_rk - before calling hetri_3 function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zhetri_3_obj->inforef = ZHETRF_RK(   zhetri_3_obj->matrix_layout,
                                    zhetri_3_obj->uplo,
                                    zhetri_3_obj->n,
                                    zhetri_3_obj->aref,
                                    zhetri_3_obj->lda,
                                    zhetri_3_obj->eref,
                                    zhetri_3_obj->ipivref);

    zhetri_3_obj->inforef = ZHETRI_3(   zhetri_3_obj->matrix_layout,
                                    zhetri_3_obj->uplo, 
                                    zhetri_3_obj->n,
                                    zhetri_3_obj->aref,
                                    zhetri_3_obj->lda,
                             (const lapack_complex_double *)zhetri_3_obj->eref,
                                    (const lapack_int *)zhetri_3_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    zhetri_3_obj->info  = LAPACKE_zhetrf_rk( zhetri_3_obj->matrix_layout,
                                        zhetri_3_obj->uplo,
                                        zhetri_3_obj->n,
                                        zhetri_3_obj->a,
                                        zhetri_3_obj->lda,
										zhetri_3_obj->e,
                                        zhetri_3_obj->ipiv);

    zhetri_3_obj->info = LAPACKE_zhetri_3(  zhetri_3_obj->matrix_layout, 
                                        zhetri_3_obj->uplo,
                                        zhetri_3_obj->n, 
                                        zhetri_3_obj->a, 
                                        zhetri_3_obj->lda,
                             (const lapack_complex_double *)zhetri_3_obj->e,
                                        (const lapack_int *)zhetri_3_obj->ipiv);

    if( zhetri_3_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhetri_3 is wrong\n", zhetri_3_obj->info );
    }
    if( zhetri_3_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetri_3 is wrong\n", 
        zhetri_3_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zhetri_3_obj->diff =  computeDiff_z( (zhetri_3_obj->n)*(zhetri_3_obj->lda),
                           zhetri_3_obj->a, zhetri_3_obj->aref );
}

TEST_F(zhetri_3_test, zhetri_31) {
    EXPECT_NEAR(0.0, zhetri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri_3_test, zhetri_32) {
    EXPECT_NEAR(0.0, zhetri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri_3_test, zhetri_33) {
    EXPECT_NEAR(0.0, zhetri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri_3_test, zhetri_34) {
    EXPECT_NEAR(0.0, zhetri_3_obj->diff, LAPACKE_GTEST_THRESHOLD);
}