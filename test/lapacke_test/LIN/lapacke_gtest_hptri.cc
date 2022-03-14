#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define hptri_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (ipiv != NULL) free (ipiv); \
  if (ipivref != NULL)free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;


/* Begin hptri_scomplex_parameters  class definition */
class hptri_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int *ipiv, *ipivref; // The pivot indices

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      hptri_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i);
              
      ~hptri_scomplex_parameters (); 
};  /* end of hptri_scomplex_parameters  class definition */


/* Constructor hptri_scomplex_parameters definition */
hptri_scomplex_parameters:: hptri_scomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n hptri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       hptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

hptri_scomplex_parameters:: ~hptri_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hptri_scomplex_parameters object: destructor invoked. \n");
#endif
   hptri_free();
}

//  Test fixture class definition
class chptri_test  : public  ::testing::Test {
public:
   hptri_scomplex_parameters  *chptri_obj;
   void SetUp();  
   void TearDown () { delete chptri_obj; }
};


void chptri_test::SetUp(){

    /* LAPACKE CHPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_chptri) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_float * a, 
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_chptri CHPTRI;

    typedef int (*Fptr_NL_LAPACKE_chptrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, 
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_chptrf CHPTRF;

    chptri_obj = new hptri_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);

    chptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CHPTRI = (Fptr_NL_LAPACKE_chptri)dlsym(chptri_obj->hModule, "LAPACKE_chptri");
    ASSERT_TRUE(CHPTRI != NULL) << "failed to get the Netlib LAPACKE_chptri symbol";
    
    CHPTRF = (Fptr_NL_LAPACKE_chptrf)dlsym(chptri_obj->hModule,"LAPACKE_chptrf");
    ASSERT_TRUE(CHPTRF != NULL) << "failed to get the Netlib LAPACKE_chptrf symbol";

    /* Pre condition: need to call hptrf - before calling hptri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    chptri_obj->inforef = CHPTRF(   chptri_obj->matrix_layout,
                                    chptri_obj->uplo,
                                    chptri_obj->n,
                                    chptri_obj->aref,
                                    chptri_obj->ipivref);

    chptri_obj->inforef = CHPTRI(   chptri_obj->matrix_layout,
                                    chptri_obj->uplo, 
                                    chptri_obj->n,
                                    chptri_obj->aref,
                                    chptri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    chptri_obj->info  = LAPACKE_chptrf( chptri_obj->matrix_layout,
                                        chptri_obj->uplo,
                                        chptri_obj->n,
                                        chptri_obj->a,
                                        chptri_obj->ipiv);

    chptri_obj->info = LAPACKE_chptri(  chptri_obj->matrix_layout, 
                                        chptri_obj->uplo,
                                        chptri_obj->n, 
                                        chptri_obj->a, 
                                        chptri_obj->ipiv);

    if( chptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chptri is wrong\n", chptri_obj->info );
    }
    if( chptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chptri is wrong\n", 
        chptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    chptri_obj->diff =  computeDiff_c( (chptri_obj->n)*(chptri_obj->n),
                           chptri_obj->a, chptri_obj->aref );
}

TEST_F(chptri_test, chptri1) {
    EXPECT_NEAR(0.0, chptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chptri_test, chptri2) {
    EXPECT_NEAR(0.0, chptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chptri_test, chptri3) {
    EXPECT_NEAR(0.0, chptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chptri_test, chptri4) {
    EXPECT_NEAR(0.0, chptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin hptri_dcomplex_parameters  class definition */
class hptri_dcomplex_parameters{
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
      hptri_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i);
              
      ~hptri_dcomplex_parameters (); 
};  /* end of hptri_dcomplex_parameters  class definition */


/* Constructor hptri_dcomplex_parameters definition */
hptri_dcomplex_parameters:: hptri_dcomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n hptri Double:  n: %d, Uplo: %c  \n",
             n, uplo);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       hptri_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
    
   } /* end of Constructor  */

hptri_dcomplex_parameters:: ~hptri_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hptri_dcomplex_parameters object: destructor invoked. \n");
#endif
   hptri_free();
}

//  Test fixture class definition
class zhptri_test  : public  ::testing::Test {
public:
   hptri_dcomplex_parameters  *zhptri_obj;
   void SetUp();  
   void TearDown () { delete zhptri_obj; }
};


void zhptri_test::SetUp(){

    /* LAPACKE ZHPTRI prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptri) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_double * a, 
                                const lapack_int *ipiv  );

    Fptr_NL_LAPACKE_zhptri ZHPTRI;

    typedef int (*Fptr_NL_LAPACKE_zhptrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, 
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_zhptrf ZHPTRF;

    zhptri_obj = new hptri_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n );

    idx = Circular_Increment_Index(idx);

    zhptri_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhptri_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhptri_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhptri_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZHPTRI = (Fptr_NL_LAPACKE_zhptri)dlsym(zhptri_obj->hModule, "LAPACKE_zhptri");
    ASSERT_TRUE(ZHPTRI != NULL) << "failed to get the Netlib LAPACKE_zhptri symbol";
    
    ZHPTRF = (Fptr_NL_LAPACKE_zhptrf)dlsym(zhptri_obj->hModule,"LAPACKE_zhptrf");
    ASSERT_TRUE(ZHPTRF != NULL) << "failed to get the Netlib LAPACKE_zhptrf symbol";

    /* Pre condition: need to call hptrf - before calling hptri function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zhptri_obj->inforef = ZHPTRF(   zhptri_obj->matrix_layout,
                                    zhptri_obj->uplo,
                                    zhptri_obj->n,
                                    zhptri_obj->aref,
                                    zhptri_obj->ipivref);

    zhptri_obj->inforef = ZHPTRI(   zhptri_obj->matrix_layout,
                                    zhptri_obj->uplo, 
                                    zhptri_obj->n,
                                    zhptri_obj->aref,
                                    zhptri_obj->ipivref);
                        

    /* Compute libflame's Lapacke o/p  */
    zhptri_obj->info  = LAPACKE_zhptrf( zhptri_obj->matrix_layout,
                                        zhptri_obj->uplo,
                                        zhptri_obj->n,
                                        zhptri_obj->a,
                                        zhptri_obj->ipiv);

    zhptri_obj->info = LAPACKE_zhptri(  zhptri_obj->matrix_layout, 
                                        zhptri_obj->uplo,
                                        zhptri_obj->n, 
                                        zhptri_obj->a, 
                                        zhptri_obj->ipiv);

    if( zhptri_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhptri is wrong\n", zhptri_obj->info );
    }
    if( zhptri_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhptri is wrong\n", 
        zhptri_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zhptri_obj->diff =  computeDiff_z( (zhptri_obj->n)*(zhptri_obj->n),
                           zhptri_obj->a, zhptri_obj->aref );
}

TEST_F(zhptri_test, zhptri1) {
    EXPECT_NEAR(0.0, zhptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhptri_test, zhptri2) {
    EXPECT_NEAR(0.0, zhptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhptri_test, zhptri3) {
    EXPECT_NEAR(0.0, zhptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhptri_test, zhptri4) {
    EXPECT_NEAR(0.0, zhptri_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

