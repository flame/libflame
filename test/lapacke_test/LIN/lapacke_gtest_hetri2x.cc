#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define hetri2x_free() \
  if (a != NULL)    free (a  ); \
  if (aref != NULL) free (aref); \
  if (ipiv != NULL) free (ipiv); \
  if (ipivref != NULL)free (ipivref); \
  if( hModule != NULL) dlclose(hModule); \
  if( dModule != NULL) dlclose(dModule)


// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin hetri2x_scomplex_parameters  class definition */
class hetri2x_scomplex_parameters{
   public:
      float diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices
	  lapack_int nb; // block size

      /* Input/ Output parameters */
      lapack_complex_float *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      hetri2x_scomplex_parameters ( int matrix_layout_i, char uplo_i, 
                                lapack_int n_i, lapack_int lda_i,
								lapack_int nb_i);
              
      ~hetri2x_scomplex_parameters (); 
};  /* end of hetri2x_scomplex_parameters  class definition */


/* Constructor hetri2x_scomplex_parameters definition */
hetri2x_scomplex_parameters:: hetri2x_scomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i,
						lapack_int nb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
	nb = nb_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n hetri2x Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_scomplex_parameters object: malloc error.";
       hetri2x_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    //lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, n*n);
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix(a, aref, n, n, 'H');
    
   } /* end of Constructor  */

hetri2x_scomplex_parameters:: ~hetri2x_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetri2x_scomplex_parameters object: destructor invoked. \n");
#endif
   hetri2x_free();
}

//  Test fixture class definition
class chetri2x_test  : public  ::testing::Test {
public:
   hetri2x_scomplex_parameters  *chetri2x_obj;
   void SetUp();  
   void TearDown () { delete chetri2x_obj; }
};


void chetri2x_test::SetUp(){

    /* LAPACKE CHETRI2X prototype */
    typedef int (*Fptr_NL_LAPACKE_chetri2x) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_float * a, lapack_int lda,
                                const lapack_int *ipiv, lapack_int nb );

    Fptr_NL_LAPACKE_chetri2x CHETRI2X;

    typedef int (*Fptr_NL_LAPACKE_chetrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_float *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_chetrf CHETRF;

    chetri2x_obj = new hetri2x_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
						   lin_solver_paramslist[idx].lda/2);

    idx = Circular_Increment_Index(idx);

    chetri2x_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetri2x_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetri2x_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetri2x_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CHETRI2X = (Fptr_NL_LAPACKE_chetri2x)dlsym(chetri2x_obj->hModule, "LAPACKE_chetri2x");
    ASSERT_TRUE(CHETRI2X != NULL) << "failed to get the Netlib LAPACKE_chetri2x symbol";
    
    CHETRF = (Fptr_NL_LAPACKE_chetrf)dlsym(chetri2x_obj->hModule,"LAPACKE_chetrf");
    ASSERT_TRUE(CHETRF != NULL) << "failed to get the Netlib LAPACKE_chetrf symbol";

    /* Pre condition: need to call hetrf - before calling hetri2x function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    chetri2x_obj->inforef = CHETRF(   chetri2x_obj->matrix_layout,
                                    chetri2x_obj->uplo,
                                    chetri2x_obj->n,
                                    chetri2x_obj->aref,
                                    chetri2x_obj->lda,
                                    chetri2x_obj->ipivref);

    chetri2x_obj->inforef = CHETRI2X(   chetri2x_obj->matrix_layout,
                                    chetri2x_obj->uplo, 
                                    chetri2x_obj->n,
                                    chetri2x_obj->aref,
                                    chetri2x_obj->lda,
                                    (const lapack_int *)chetri2x_obj->ipivref,
									chetri2x_obj->nb);
                        

    /* Compute libflame's Lapacke o/p  */
    chetri2x_obj->info  = LAPACKE_chetrf( chetri2x_obj->matrix_layout,
                                        chetri2x_obj->uplo,
                                        chetri2x_obj->n,
                                        chetri2x_obj->a,
                                        chetri2x_obj->lda,
                                        chetri2x_obj->ipiv);

    chetri2x_obj->info = LAPACKE_chetri2x(  chetri2x_obj->matrix_layout, 
                                        chetri2x_obj->uplo,
                                        chetri2x_obj->n, 
                                        chetri2x_obj->a, 
                                        chetri2x_obj->lda,
                                        (const lapack_int *)chetri2x_obj->ipiv,
										chetri2x_obj->nb);

    if( chetri2x_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chetri2x is wrong\n", chetri2x_obj->info );
    }
    if( chetri2x_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetri2x is wrong\n", 
        chetri2x_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    chetri2x_obj->diff =  computeDiff_c( (chetri2x_obj->n)*(chetri2x_obj->lda),
                           chetri2x_obj->a, chetri2x_obj->aref );
}

TEST_F(chetri2x_test, chetri2x1) {
    EXPECT_NEAR(0.0, chetri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri2x_test, chetri2x2) {
    EXPECT_NEAR(0.0, chetri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri2x_test, chetri2x3) {
    EXPECT_NEAR(0.0, chetri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetri2x_test, chetri2x4) {
    EXPECT_NEAR(0.0, chetri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin hetri2x_dcomplex_parameters  class definition */
class hetri2x_dcomplex_parameters{
   public:
      double diff; //to capture reference and libflame o/ps' difference
      void *hModule, *dModule;
     
      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; //  Must be 'U' or 'L'.
      lapack_int n; // The order of A; the number of rows in B
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int *ipiv, *ipivref; // The pivot indices
	  lapack_int nb; // block size

      /* Input/ Output parameters */
      lapack_complex_double *a, *aref; //The array ab contains the matrix A

      /* Return Values */
      lapack_int info, inforef;

   public: 
      hetri2x_dcomplex_parameters ( int matrix_layout_i, char uplo_i, 
									lapack_int n_i, lapack_int lda_i,
									  lapack_int nb_i);
              
      ~hetri2x_dcomplex_parameters (); 
};  /* end of hetri2x_dcomplex_parameters  class definition */


/* Constructor hetri2x_dcomplex_parameters definition */
hetri2x_dcomplex_parameters:: hetri2x_dcomplex_parameters ( int matrix_layout_i, 
                       char uplo_i, lapack_int n_i, lapack_int lda_i,
								lapack_int nb_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = lda_i;
	nb = nb_i;
    hModule = NULL;
    dModule = NULL;
    info = 0;
    inforef = 0;
#if LAPACKE_TEST_VERBOSE
   printf(" \n hetri2x Double:  n: %d, Uplo: %c lda: %d \n",
             n, uplo, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &a, &aref, (n*n));
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);

    if( (a==NULL) || (aref==NULL)  ){
       EXPECT_FALSE( true) << "pbtrs_dcomplex_parameters object: malloc error.";
       hetri2x_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
    //lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, n*n);
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix(a, aref, n, n, 'H');
    
   } /* end of Constructor  */

hetri2x_dcomplex_parameters:: ~hetri2x_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" hetri2x_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetri2x_free();
}

//  Test fixture class definition
class zhetri2x_test  : public  ::testing::Test {
public:
   hetri2x_dcomplex_parameters  *zhetri2x_obj;
   void SetUp();  
   void TearDown () { delete zhetri2x_obj; }
};


void zhetri2x_test::SetUp(){

    /* LAPACKE ZHETRI2X prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetri2x) ( int matrix_layout, char uplo,
                                lapack_int n,  lapack_complex_double * a,
								lapack_int lda,
                                const lapack_int *ipiv, lapack_int nb  );

    Fptr_NL_LAPACKE_zhetri2x ZHETRI2X;

    typedef int (*Fptr_NL_LAPACKE_zhetrf) ( int matrix_layout ,char uplo,
                               lapack_int n, lapack_complex_double *a, lapack_int lda,
                               lapack_int *ipiv );
    Fptr_NL_LAPACKE_zhetrf ZHETRF;

    zhetri2x_obj = new hetri2x_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
                           lin_solver_paramslist[idx].n,
                           lin_solver_paramslist[idx].lda,
						   lin_solver_paramslist[idx].lda/2);

    idx = Circular_Increment_Index(idx);

    zhetri2x_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetri2x_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetri2x_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetri2x_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ZHETRI2X = (Fptr_NL_LAPACKE_zhetri2x)dlsym(zhetri2x_obj->hModule, "LAPACKE_zhetri2x");
    ASSERT_TRUE(ZHETRI2X != NULL) << "failed to get the Netlib LAPACKE_zhetri2x symbol";
    
    ZHETRF = (Fptr_NL_LAPACKE_zhetrf)dlsym(zhetri2x_obj->hModule,"LAPACKE_zhetrf");
    ASSERT_TRUE(ZHETRF != NULL) << "failed to get the Netlib LAPACKE_zhetrf symbol";

    /* Pre condition: need to call hetrf - before calling hetri2x function */
    /* Compute the reference o/p by invoking Netlib-LapackE's API */ 
    zhetri2x_obj->inforef = ZHETRF(   zhetri2x_obj->matrix_layout,
                                    zhetri2x_obj->uplo,
                                    zhetri2x_obj->n,
                                    zhetri2x_obj->aref,
                                    zhetri2x_obj->lda,
                                    zhetri2x_obj->ipivref);

    zhetri2x_obj->inforef = ZHETRI2X(   zhetri2x_obj->matrix_layout,
                                    zhetri2x_obj->uplo, 
                                    zhetri2x_obj->n,
                                    zhetri2x_obj->aref,
                                    zhetri2x_obj->lda,
                                    zhetri2x_obj->ipivref,
									zhetri2x_obj->nb);
                        

    /* Compute libflame's Lapacke o/p  */
    zhetri2x_obj->info  = LAPACKE_zhetrf( zhetri2x_obj->matrix_layout,
                                        zhetri2x_obj->uplo,
                                        zhetri2x_obj->n,
                                        zhetri2x_obj->a,
                                        zhetri2x_obj->lda,
                                        zhetri2x_obj->ipiv);

    zhetri2x_obj->info = LAPACKE_zhetri2x(  zhetri2x_obj->matrix_layout, 
                                        zhetri2x_obj->uplo,
                                        zhetri2x_obj->n, 
                                        zhetri2x_obj->a, 
                                        zhetri2x_obj->lda,
                                        zhetri2x_obj->ipiv,
										zhetri2x_obj->nb);

    if( zhetri2x_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhetri2x is wrong\n", zhetri2x_obj->info );
    }
    if( zhetri2x_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetri2x is wrong\n", 
        zhetri2x_obj->inforef );
    }
    /* Compute Difference between libflame and Netlib o/ps  */
    zhetri2x_obj->diff =  computeDiff_z( (zhetri2x_obj->n)*(zhetri2x_obj->lda),
                           zhetri2x_obj->a, zhetri2x_obj->aref );
}

TEST_F(zhetri2x_test, zhetri2x1) {
    EXPECT_NEAR(0.0, zhetri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri2x_test, zhetri2x2) {
    EXPECT_NEAR(0.0, zhetri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri2x_test, zhetri2x3) {
    EXPECT_NEAR(0.0, zhetri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetri2x_test, zhetri2x4) {
    EXPECT_NEAR(0.0, zhetri2x_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

