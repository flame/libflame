#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define sytrf_aa_2stage_free() \
       if (tbref != NULL) free (tbref); \
       if (tb != NULL)    free (tb   ); \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (ipiv != NULL) free (ipiv); \
       if (ipivref != NULL)free (ipivref); \
       if (ipiv2 != NULL) free (ipiv2); \
       if (ipiv2ref != NULL)free (ipiv2ref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ssytrf_aa_2stage_float_parameters  class definition */
class ssytrf_aa_2stage_float_parameters{

   public:
      void *hModule,*dModule;
      float diff_tb;
      int diff_ipiv, diff_ipiv2;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ltb;  //  leading dimension of 'b'
      float *a, *aref; //The array 'a' contains the matrix A

      /* Output parameters */
      float *tb, *tbref; //right-hand sides for the systems of equations.
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int *ipiv2, *ipiv2ref; // The pivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      ssytrf_aa_2stage_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i);
      ~ssytrf_aa_2stage_float_parameters ();
};  /* end of ssytrf_aa_2stage_float_parameters  class definition */


/* Constructor ssytrf_aa_2stage_float_parameters definition */
ssytrf_aa_2stage_float_parameters:: ssytrf_aa_2stage_float_parameters (
					int matrix_layout_i, char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = n;
    ltb = 4*n;

    hModule = NULL;
    dModule = NULL;
	diff_tb = 0;
    diff_ipiv = 0;
	diff_ipiv2 = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ssytrf_aa_2stage float:  n: %d, uplo: %c  ltb: %d lda: %d \n",
             n, uplo, ltb, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair(  &a, &aref, (n*lda));
    lapacke_gtest_alloc_float_buffer_pair( &tb, &tbref, ltb*ltb);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv2, &ipiv2ref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (tb==NULL) || (tbref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ||\
        (ipiv2==NULL) || (ipiv2ref==NULL) ){
       EXPECT_FALSE( true) << "ssytrf_aa_2stage_double_parameters object: malloc error.";
       sytrf_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

ssytrf_aa_2stage_float_parameters:: ~ssytrf_aa_2stage_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ssytrf_aa_2stage_float_parameters object: destructor invoked. \n");
#endif
   sytrf_aa_2stage_free();
}


//  Test fixture class definition
class ssytrf_aa_2stage_test  : public  ::testing::Test {
public:
   ssytrf_aa_2stage_float_parameters  *ssytrf_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete ssytrf_aa_2stage_obj; }
};


void ssytrf_aa_2stage_test::SetUp(){

    /* LAPACKE_ssytrf_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,float *a ,lapack_int lda,
								 float *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_ssytrf_aa_2stage SSYTRF_AA_2STAGE;

    ssytrf_aa_2stage_obj = new  ssytrf_aa_2stage_float_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n);

    ssytrf_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytrf_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytrf_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytrf_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_ssytrf_aa_2stage)dlsym(
							ssytrf_aa_2stage_obj->hModule, "LAPACKE_ssytrf_aa_2stage");
    ASSERT_TRUE(SSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_ssytrf_aa_2stage symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    //ssytrf_aa_2stage_obj->inforef = SSYTRF_AA_2STAGE(
    ssytrf_aa_2stage_obj->inforef = LAPACKE_ssytrf_aa_2stage(
										ssytrf_aa_2stage_obj->matrix_layout,
										ssytrf_aa_2stage_obj->uplo,
										ssytrf_aa_2stage_obj->n,
										ssytrf_aa_2stage_obj->aref,
										ssytrf_aa_2stage_obj->lda,
										ssytrf_aa_2stage_obj->tbref,
										ssytrf_aa_2stage_obj->ltb,
										ssytrf_aa_2stage_obj->ipivref,
										ssytrf_aa_2stage_obj->ipiv2ref
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssytrf_aa_2stage_obj->info = LAPACKE_ssytrf_aa_2stage(
										ssytrf_aa_2stage_obj->matrix_layout,
										ssytrf_aa_2stage_obj->uplo,
										ssytrf_aa_2stage_obj->n,
										ssytrf_aa_2stage_obj->a,
										ssytrf_aa_2stage_obj->lda,
										ssytrf_aa_2stage_obj->tb,
										ssytrf_aa_2stage_obj->ltb,
										ssytrf_aa_2stage_obj->ipiv,
										ssytrf_aa_2stage_obj->ipiv2
										);
								 
								 
    if( ssytrf_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrf_aa_2stage is wrong\n",
                    ssytrf_aa_2stage_obj->info );
    }
    if( ssytrf_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrf_aa_2stage is wrong\n",
        ssytrf_aa_2stage_obj->inforef );
    }
	
	ssytrf_aa_2stage_obj->diff_tb =  computeDiff_s( ssytrf_aa_2stage_obj->ltb,
                           ssytrf_aa_2stage_obj->tb,
						   ssytrf_aa_2stage_obj->tbref );

    ssytrf_aa_2stage_obj->diff_ipiv = computeDiff_i( ssytrf_aa_2stage_obj->n,
							ssytrf_aa_2stage_obj->ipiv,
							ssytrf_aa_2stage_obj->ipivref);

    ssytrf_aa_2stage_obj->diff_ipiv2 = computeDiff_i( ssytrf_aa_2stage_obj->n,
							ssytrf_aa_2stage_obj->ipiv2,
							ssytrf_aa_2stage_obj->ipiv2ref);

}

TEST_F(ssytrf_aa_2stage_test, ssytrf_aa_2stage1) {
    EXPECT_NEAR(0.0, ssytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssytrf_aa_2stage_test, ssytrf_aa_2stage2) {
    EXPECT_NEAR(0.0, ssytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssytrf_aa_2stage_test, ssytrf_aa_2stage3) {
    EXPECT_NEAR(0.0, ssytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssytrf_aa_2stage_test, ssytrf_aa_2stage4) {
    EXPECT_NEAR(0.0, ssytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

/* Begin dsytrf_aa_2stage_double_parameters  class definition */
class dsytrf_aa_2stage_double_parameters{

   public:
      void *hModule,*dModule;
      double diff_tb;
      int diff_ipiv, diff_ipiv2;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ltb;  //  leading dimension of 'b'
      double *a, *aref; //The array 'a' contains the matrix A

      /* Output parameters */
      double *tb, *tbref; //right-hand sides for the systems of equations.
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int *ipiv2, *ipiv2ref; // The pivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      dsytrf_aa_2stage_double_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i);
      ~dsytrf_aa_2stage_double_parameters ();
};  /* end of dsytrf_aa_2stage_double_parameters  class definition */


/* Constructor dsytrf_aa_2stage_double_parameters definition */
dsytrf_aa_2stage_double_parameters:: dsytrf_aa_2stage_double_parameters (
					int matrix_layout_i, char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = n;
    ltb = 4*n;

    hModule = NULL;
    dModule = NULL;
	diff_tb = 0;
    diff_ipiv = 0;
	diff_ipiv2 = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n dsytrf_aa_2stage double:  n: %d, uplo: %c  ltb: %d lda: %d \n",
             n, uplo, ltb, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair(  &a, &aref, (n*lda));
    lapacke_gtest_alloc_double_buffer_pair( &tb, &tbref, ltb*ltb);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv2, &ipiv2ref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (tb==NULL) || (tbref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ||\
        (ipiv2==NULL) || (ipiv2ref==NULL) ){
       EXPECT_FALSE( true) << "dsytrf_aa_2stage_double_parameters object: malloc error.";
       sytrf_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

dsytrf_aa_2stage_double_parameters:: ~dsytrf_aa_2stage_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" dsytrf_aa_2stage_double_parameters object: destructor invoked. \n");
#endif
   sytrf_aa_2stage_free();
}


//  Test fixture class definition
class dsytrf_aa_2stage_test  : public  ::testing::Test {
public:
   dsytrf_aa_2stage_double_parameters  *dsytrf_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete dsytrf_aa_2stage_obj; }
};


void dsytrf_aa_2stage_test::SetUp(){

    /* LAPACKE_dsytrf_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,double *a ,lapack_int lda,
								 double *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_dsytrf_aa_2stage DSYTRF_AA_2STAGE;

    dsytrf_aa_2stage_obj = new  dsytrf_aa_2stage_double_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n);

    dsytrf_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytrf_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytrf_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytrf_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_dsytrf_aa_2stage)dlsym(
							dsytrf_aa_2stage_obj->hModule, "LAPACKE_dsytrf_aa_2stage");
    ASSERT_TRUE(DSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_dsytrf_aa_2stage symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    //dsytrf_aa_2stage_obj->inforef = DSYTRF_AA_2STAGE(
    dsytrf_aa_2stage_obj->inforef = LAPACKE_dsytrf_aa_2stage(
										dsytrf_aa_2stage_obj->matrix_layout,
										dsytrf_aa_2stage_obj->uplo,
										dsytrf_aa_2stage_obj->n,
										dsytrf_aa_2stage_obj->aref,
										dsytrf_aa_2stage_obj->lda,
										dsytrf_aa_2stage_obj->tbref,
										dsytrf_aa_2stage_obj->ltb,
										dsytrf_aa_2stage_obj->ipivref,
										dsytrf_aa_2stage_obj->ipiv2ref
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    dsytrf_aa_2stage_obj->info = LAPACKE_dsytrf_aa_2stage(
										dsytrf_aa_2stage_obj->matrix_layout,
										dsytrf_aa_2stage_obj->uplo,
										dsytrf_aa_2stage_obj->n,
										dsytrf_aa_2stage_obj->a,
										dsytrf_aa_2stage_obj->lda,
										dsytrf_aa_2stage_obj->tb,
										dsytrf_aa_2stage_obj->ltb,
										dsytrf_aa_2stage_obj->ipiv,
										dsytrf_aa_2stage_obj->ipiv2
										);
								 
								 
    if( dsytrf_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrf_aa_2stage is wrong\n",
                    dsytrf_aa_2stage_obj->info );
    }
    if( dsytrf_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrf_aa_2stage is wrong\n",
        dsytrf_aa_2stage_obj->inforef );
    }
	
	dsytrf_aa_2stage_obj->diff_tb =  computeDiff_d( dsytrf_aa_2stage_obj->ltb,
                           dsytrf_aa_2stage_obj->tb,
						   dsytrf_aa_2stage_obj->tbref );

    dsytrf_aa_2stage_obj->diff_ipiv = computeDiff_i( dsytrf_aa_2stage_obj->n,
							dsytrf_aa_2stage_obj->ipiv,
							dsytrf_aa_2stage_obj->ipivref);

    dsytrf_aa_2stage_obj->diff_ipiv2 = computeDiff_i( dsytrf_aa_2stage_obj->n,
							dsytrf_aa_2stage_obj->ipiv2,
							dsytrf_aa_2stage_obj->ipiv2ref);

}

TEST_F(dsytrf_aa_2stage_test, dsytrf_aa_2stage1) {
    EXPECT_NEAR(0.0, dsytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsytrf_aa_2stage_test, dsytrf_aa_2stage2) {
    EXPECT_NEAR(0.0, dsytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsytrf_aa_2stage_test, dsytrf_aa_2stage3) {
    EXPECT_NEAR(0.0, dsytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsytrf_aa_2stage_test, dsytrf_aa_2stage4) {
    EXPECT_NEAR(0.0, dsytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

/* Begin csytrf_aa_2stage_scomplex_parameters  class definition */
class csytrf_aa_2stage_scomplex_parameters{

   public:
      void *hModule,*dModule;
      float diff_tb;
      int diff_ipiv, diff_ipiv2;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ltb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A

      /* Output parameters */
      lapack_complex_float *tb, *tbref; //right-hand sides for the systems of equations.
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int *ipiv2, *ipiv2ref; // The pivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      csytrf_aa_2stage_scomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i);
      ~csytrf_aa_2stage_scomplex_parameters ();
};  /* end of csytrf_aa_2stage_scomplex_parameters  class definition */


/* Constructor csytrf_aa_2stage_scomplex_parameters definition */
csytrf_aa_2stage_scomplex_parameters:: csytrf_aa_2stage_scomplex_parameters (
					int matrix_layout_i, char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = n;
    ltb = 4*n;

    hModule = NULL;
    dModule = NULL;
	diff_tb = 0;
    diff_ipiv = 0;
	diff_ipiv2 = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n csytrf_aa_2stage lapack_complex_float:  n: %d, uplo: %c  ltb: %d lda: %d \n",
             n, uplo, ltb, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair(  &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &tb, &tbref, ltb*ltb);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv2, &ipiv2ref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (tb==NULL) || (tbref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ||\
        (ipiv2==NULL) || (ipiv2ref==NULL) ){
       EXPECT_FALSE( true) << "csytrf_aa_2stage_double_parameters object: malloc error.";
       sytrf_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

csytrf_aa_2stage_scomplex_parameters:: ~csytrf_aa_2stage_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" csytrf_aa_2stage_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrf_aa_2stage_free();
}


//  Test fixture class definition
class csytrf_aa_2stage_test  : public  ::testing::Test {
public:
   csytrf_aa_2stage_scomplex_parameters  *csytrf_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete csytrf_aa_2stage_obj; }
};


void csytrf_aa_2stage_test::SetUp(){

    /* LAPACKE_csytrf_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,lapack_complex_float *a ,lapack_int lda,
								 lapack_complex_float *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_csytrf_aa_2stage CSYTRF_AA_2STAGE;

    csytrf_aa_2stage_obj = new  csytrf_aa_2stage_scomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n);

    csytrf_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytrf_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytrf_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytrf_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_csytrf_aa_2stage)dlsym(
							csytrf_aa_2stage_obj->hModule, "LAPACKE_csytrf_aa_2stage");
    ASSERT_TRUE(CSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_csytrf_aa_2stage symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    //csytrf_aa_2stage_obj->inforef = CSYTRF_AA_2STAGE(
    csytrf_aa_2stage_obj->inforef = LAPACKE_csytrf_aa_2stage(
										csytrf_aa_2stage_obj->matrix_layout,
										csytrf_aa_2stage_obj->uplo,
										csytrf_aa_2stage_obj->n,
										csytrf_aa_2stage_obj->aref,
										csytrf_aa_2stage_obj->lda,
										csytrf_aa_2stage_obj->tbref,
										csytrf_aa_2stage_obj->ltb,
										csytrf_aa_2stage_obj->ipivref,
										csytrf_aa_2stage_obj->ipiv2ref
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    csytrf_aa_2stage_obj->info = LAPACKE_csytrf_aa_2stage(
										csytrf_aa_2stage_obj->matrix_layout,
										csytrf_aa_2stage_obj->uplo,
										csytrf_aa_2stage_obj->n,
										csytrf_aa_2stage_obj->a,
										csytrf_aa_2stage_obj->lda,
										csytrf_aa_2stage_obj->tb,
										csytrf_aa_2stage_obj->ltb,
										csytrf_aa_2stage_obj->ipiv,
										csytrf_aa_2stage_obj->ipiv2
										);
								 
								 
    if( csytrf_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrf_aa_2stage is wrong\n",
                    csytrf_aa_2stage_obj->info );
    }
    if( csytrf_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrf_aa_2stage is wrong\n",
        csytrf_aa_2stage_obj->inforef );
    }
	
	csytrf_aa_2stage_obj->diff_tb =  computeDiff_c( csytrf_aa_2stage_obj->ltb,
                           csytrf_aa_2stage_obj->tb,
						   csytrf_aa_2stage_obj->tbref );

    csytrf_aa_2stage_obj->diff_ipiv = computeDiff_i( csytrf_aa_2stage_obj->n,
							csytrf_aa_2stage_obj->ipiv,
							csytrf_aa_2stage_obj->ipivref);

    csytrf_aa_2stage_obj->diff_ipiv2 = computeDiff_i( csytrf_aa_2stage_obj->n,
							csytrf_aa_2stage_obj->ipiv2,
							csytrf_aa_2stage_obj->ipiv2ref);

}

TEST_F(csytrf_aa_2stage_test, csytrf_aa_2stage1) {
    EXPECT_NEAR(0.0, csytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(csytrf_aa_2stage_test, csytrf_aa_2stage2) {
    EXPECT_NEAR(0.0, csytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(csytrf_aa_2stage_test, csytrf_aa_2stage3) {
    EXPECT_NEAR(0.0, csytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(csytrf_aa_2stage_test, csytrf_aa_2stage4) {
    EXPECT_NEAR(0.0, csytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

/* Begin zsytrf_aa_2stage_dcomplex_parameters  class definition */
class zsytrf_aa_2stage_dcomplex_parameters{

   public:
      void *hModule,*dModule;
      double diff_tb;
      int diff_ipiv, diff_ipiv2;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ltb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A

      /* Output parameters */
      lapack_complex_double *tb, *tbref; //right-hand sides for the systems of equations.
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int *ipiv2, *ipiv2ref; // The pivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      zsytrf_aa_2stage_dcomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i);
      ~zsytrf_aa_2stage_dcomplex_parameters ();
};  /* end of zsytrf_aa_2stage_dcomplex_parameters  class definition */


/* Constructor zsytrf_aa_2stage_dcomplex_parameters definition */
zsytrf_aa_2stage_dcomplex_parameters:: zsytrf_aa_2stage_dcomplex_parameters (
					int matrix_layout_i, char uplo_i, lapack_int n_i ) {

    matrix_layout = matrix_layout_i;
    n = n_i;
    uplo = uplo_i;
    lda = n;
    ltb = 4*n;

    hModule = NULL;
    dModule = NULL;
	diff_tb = 0;
    diff_ipiv = 0;
	diff_ipiv2 = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n zsytrf_aa_2stage lapack_complex_double:  n: %d, uplo: %c  ltb: %d lda: %d \n",
             n, uplo, ltb, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(  &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &tb, &tbref, ltb*ltb);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv2, &ipiv2ref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (tb==NULL) || (tbref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ||\
        (ipiv2==NULL) || (ipiv2ref==NULL) ){
       EXPECT_FALSE( true) << "zsytrf_aa_2stage_double_parameters object: malloc error.";
       sytrf_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, uplo);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

zsytrf_aa_2stage_dcomplex_parameters:: ~zsytrf_aa_2stage_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" zsytrf_aa_2stage_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrf_aa_2stage_free();
}


//  Test fixture class definition
class zsytrf_aa_2stage_test  : public  ::testing::Test {
public:
   zsytrf_aa_2stage_dcomplex_parameters  *zsytrf_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete zsytrf_aa_2stage_obj; }
};


void zsytrf_aa_2stage_test::SetUp(){

    /* LAPACKE_zsytrf_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,lapack_complex_double *a ,lapack_int lda,
								 lapack_complex_double *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_zsytrf_aa_2stage ZSYTRF_AA_2STAGE;

    zsytrf_aa_2stage_obj = new  zsytrf_aa_2stage_dcomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n);

    zsytrf_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytrf_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytrf_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytrf_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_zsytrf_aa_2stage)dlsym(
							zsytrf_aa_2stage_obj->hModule, "LAPACKE_zsytrf_aa_2stage");
    ASSERT_TRUE(ZSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_zsytrf_aa_2stage symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zsytrf_aa_2stage_obj->inforef = ZSYTRF_AA_2STAGE(
										zsytrf_aa_2stage_obj->matrix_layout,
										zsytrf_aa_2stage_obj->uplo,
										zsytrf_aa_2stage_obj->n,
										zsytrf_aa_2stage_obj->aref,
										zsytrf_aa_2stage_obj->lda,
										zsytrf_aa_2stage_obj->tbref,
										zsytrf_aa_2stage_obj->ltb,
										zsytrf_aa_2stage_obj->ipivref,
										zsytrf_aa_2stage_obj->ipiv2ref
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    zsytrf_aa_2stage_obj->info = LAPACKE_zsytrf_aa_2stage(
										zsytrf_aa_2stage_obj->matrix_layout,
										zsytrf_aa_2stage_obj->uplo,
										zsytrf_aa_2stage_obj->n,
										zsytrf_aa_2stage_obj->a,
										zsytrf_aa_2stage_obj->lda,
										zsytrf_aa_2stage_obj->tb,
										zsytrf_aa_2stage_obj->ltb,
										zsytrf_aa_2stage_obj->ipiv,
										zsytrf_aa_2stage_obj->ipiv2
										);
								 
								 
    if( zsytrf_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrf_aa_2stage is wrong\n",
                    zsytrf_aa_2stage_obj->info );
    }
    if( zsytrf_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrf_aa_2stage is wrong\n",
        zsytrf_aa_2stage_obj->inforef );
    }
	
	zsytrf_aa_2stage_obj->diff_tb =  computeDiff_z( zsytrf_aa_2stage_obj->ltb,
                           zsytrf_aa_2stage_obj->tb,
						   zsytrf_aa_2stage_obj->tbref );

    zsytrf_aa_2stage_obj->diff_ipiv = computeDiff_i( zsytrf_aa_2stage_obj->n,
							zsytrf_aa_2stage_obj->ipiv,
							zsytrf_aa_2stage_obj->ipivref);

    zsytrf_aa_2stage_obj->diff_ipiv2 = computeDiff_i( zsytrf_aa_2stage_obj->n,
							zsytrf_aa_2stage_obj->ipiv2,
							zsytrf_aa_2stage_obj->ipiv2ref);

}

TEST_F(zsytrf_aa_2stage_test, zsytrf_aa_2stage1) {
    EXPECT_NEAR(0.0, zsytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsytrf_aa_2stage_test, zsytrf_aa_2stage2) {
    EXPECT_NEAR(0.0, zsytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsytrf_aa_2stage_test, zsytrf_aa_2stage3) {
    EXPECT_NEAR(0.0, zsytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsytrf_aa_2stage_test, zsytrf_aa_2stage4) {
    EXPECT_NEAR(0.0, zsytrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsytrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsytrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

