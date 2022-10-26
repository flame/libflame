#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"
#include "../../lapacke_gtest_helper.h"

#define sytrs_aa_2stage_free() \
       if (tbref != NULL) free (tbref); \
       if (tb != NULL)    free (tb   ); \
       if (a != NULL)    free (a   ); \
       if (aref != NULL) free (aref); \
       if (b != NULL)    free (b   ); \
       if (bref != NULL) free (bref); \
       if (ipiv != NULL) free (ipiv); \
       if (ipivref != NULL)free (ipivref); \
       if (ipiv2 != NULL) free (ipiv2); \
       if (ipiv2ref != NULL)free (ipiv2ref); \
       if( hModule != NULL) dlclose(hModule); \
       if(dModule != NULL) dlclose(dModule)

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin ssytrs_aa_2stage_float_parameters  class definition */
class ssytrs_aa_2stage_float_parameters{

   public:
      void *hModule,*dModule;
      float diff_tb, diff_b;
      int diff_ipiv, diff_ipiv2;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ltb;  //  dimension of 'tb'
      lapack_int ldb;  //  leading dimension of 'b'
      float *a, *aref; //The array 'a' contains the matrix A
      float *b, *bref; //The array 'a' contains the matrix A

      /* Output parameters */
      float *tb, *tbref; //right-hand sides for the systems of equations.
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int *ipiv2, *ipiv2ref; // The pivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      ssytrs_aa_2stage_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~ssytrs_aa_2stage_float_parameters ();
};  /* end of ssytrs_aa_2stage_float_parameters  class definition */


/* Constructor ssytrs_aa_2stage_float_parameters definition */
ssytrs_aa_2stage_float_parameters:: ssytrs_aa_2stage_float_parameters (
					int matrix_layout_i, char uplo_i, lapack_int n_i,
								lapack_int nrhs_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
	nrhs = nrhs_i;
    uplo = uplo_i;
    lda = n;
    ltb = 4*n;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

    hModule = NULL;
    dModule = NULL;
	diff_tb = 0;
	diff_b = 0;
    diff_ipiv = 0;
	diff_ipiv2 = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n ssytrs_aa_2stage float:  n: %d, uplo: %c  ltb: %d lda: %d \n",
             n, uplo, ltb, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_float_buffer_pair(  &a, &aref, (n*lda));
    lapacke_gtest_alloc_float_buffer_pair(  &b, &bref, (n*nrhs));
    lapacke_gtest_alloc_float_buffer_pair( &tb, &tbref, ltb*ltb);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv2, &ipiv2ref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (tb==NULL) || (tbref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ||\
        (ipiv2==NULL) || (ipiv2ref==NULL) ){
       EXPECT_FALSE( true) << "ssytrs_aa_2stage_double_parameters object: malloc error.";
       sytrs_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

ssytrs_aa_2stage_float_parameters:: ~ssytrs_aa_2stage_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ssytrs_aa_2stage_float_parameters object: destructor invoked. \n");
#endif
   sytrs_aa_2stage_free();
}


//  Test fixture class definition
class ssytrs_aa_2stage_test  : public  ::testing::Test {
public:
   ssytrs_aa_2stage_float_parameters  *ssytrs_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete ssytrs_aa_2stage_obj; }
};


void ssytrs_aa_2stage_test::SetUp(){

    /* LAPACKE_ssytrs_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,float *a ,lapack_int lda,
								 float *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_ssytrf_aa_2stage SSYTRF_AA_2STAGE;

    typedef int (*Fptr_NL_LAPACKE_ssytrs_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, float *a,
			lapack_int lda, float *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, float *b, lapack_int ldb );

    Fptr_NL_LAPACKE_ssytrs_aa_2stage SSYTRS_AA_2STAGE;

    ssytrs_aa_2stage_obj = new  ssytrs_aa_2stage_float_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    ssytrs_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssytrs_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssytrs_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssytrs_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_ssytrf_aa_2stage)dlsym(
							ssytrs_aa_2stage_obj->hModule, "LAPACKE_ssytrf_aa_2stage");
    ASSERT_TRUE(SSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_ssytrs_aa_2stage symbol";

    SSYTRS_AA_2STAGE = (Fptr_NL_LAPACKE_ssytrs_aa_2stage)dlsym(
							ssytrs_aa_2stage_obj->hModule, "LAPACKE_ssytrs_aa_2stage");
    ASSERT_TRUE(SSYTRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_ssytrs_aa_2stage symbol";

    ssytrs_aa_2stage_obj->inforef = SSYTRF_AA_2STAGE(
										ssytrs_aa_2stage_obj->matrix_layout,
										ssytrs_aa_2stage_obj->uplo,
										ssytrs_aa_2stage_obj->n,
										ssytrs_aa_2stage_obj->aref,
										ssytrs_aa_2stage_obj->lda,
										ssytrs_aa_2stage_obj->tbref,
										ssytrs_aa_2stage_obj->ltb,
										ssytrs_aa_2stage_obj->ipivref,
										ssytrs_aa_2stage_obj->ipiv2ref
										);

    ssytrs_aa_2stage_obj->info = LAPACKE_ssytrf_aa_2stage(
										ssytrs_aa_2stage_obj->matrix_layout,
										ssytrs_aa_2stage_obj->uplo,
										ssytrs_aa_2stage_obj->n,
										ssytrs_aa_2stage_obj->a,
										ssytrs_aa_2stage_obj->lda,
										ssytrs_aa_2stage_obj->tb,
										ssytrs_aa_2stage_obj->ltb,
										ssytrs_aa_2stage_obj->ipiv,
										ssytrs_aa_2stage_obj->ipiv2
										);

    /* Compute the Netlib-Lapacke's reference o/p */

    ssytrs_aa_2stage_obj->inforef = SSYTRS_AA_2STAGE(
										ssytrs_aa_2stage_obj->matrix_layout,
										ssytrs_aa_2stage_obj->uplo,
										ssytrs_aa_2stage_obj->n,
										ssytrs_aa_2stage_obj->nrhs,
										ssytrs_aa_2stage_obj->aref,
										ssytrs_aa_2stage_obj->lda,
										ssytrs_aa_2stage_obj->tbref,
										ssytrs_aa_2stage_obj->ltb,
										ssytrs_aa_2stage_obj->ipivref,
										ssytrs_aa_2stage_obj->ipiv2ref,
										ssytrs_aa_2stage_obj->bref,
										ssytrs_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssytrs_aa_2stage_obj->info = LAPACKE_ssytrs_aa_2stage(
										ssytrs_aa_2stage_obj->matrix_layout,
										ssytrs_aa_2stage_obj->uplo,
										ssytrs_aa_2stage_obj->n,
										ssytrs_aa_2stage_obj->nrhs,
										ssytrs_aa_2stage_obj->a,
										ssytrs_aa_2stage_obj->lda,
										ssytrs_aa_2stage_obj->tb,
										ssytrs_aa_2stage_obj->ltb,
										ssytrs_aa_2stage_obj->ipiv,
										ssytrs_aa_2stage_obj->ipiv2,
										ssytrs_aa_2stage_obj->b,
										ssytrs_aa_2stage_obj->ldb
										);
								 
								 
    if( ssytrs_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssytrs_aa_2stage is wrong\n",
                    ssytrs_aa_2stage_obj->info );
    }
    if( ssytrs_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssytrs_aa_2stage is wrong\n",
        ssytrs_aa_2stage_obj->inforef );
    }
	
	ssytrs_aa_2stage_obj->diff_tb =  computeDiff_s( ssytrs_aa_2stage_obj->ltb,
                           ssytrs_aa_2stage_obj->tb,
						   ssytrs_aa_2stage_obj->tbref );

	ssytrs_aa_2stage_obj->diff_b =  computeDiff_s( ssytrs_aa_2stage_obj->n*ssytrs_aa_2stage_obj->nrhs,
                           ssytrs_aa_2stage_obj->b,
						   ssytrs_aa_2stage_obj->bref );

    ssytrs_aa_2stage_obj->diff_ipiv = computeDiff_i( ssytrs_aa_2stage_obj->n,
							ssytrs_aa_2stage_obj->ipiv,
							ssytrs_aa_2stage_obj->ipivref);

    ssytrs_aa_2stage_obj->diff_ipiv2 = computeDiff_i( ssytrs_aa_2stage_obj->n,
							ssytrs_aa_2stage_obj->ipiv2,
							ssytrs_aa_2stage_obj->ipiv2ref);

}

TEST_F(ssytrs_aa_2stage_test, ssytrs_aa_2stage1) {
    EXPECT_NEAR(0.0, ssytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssytrs_aa_2stage_test, ssytrs_aa_2stage2) {
    EXPECT_NEAR(0.0, ssytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssytrs_aa_2stage_test, ssytrs_aa_2stage3) {
    EXPECT_NEAR(0.0, ssytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssytrs_aa_2stage_test, ssytrs_aa_2stage4) {
    EXPECT_NEAR(0.0, ssytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}


/* Begin dsytrs_aa_2stage_double_parameters  class definition */
class dsytrs_aa_2stage_double_parameters{

   public:
      void *hModule,*dModule;
      double diff_tb, diff_b;
      int diff_ipiv, diff_ipiv2;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ltb;  //  dimension of 'tb'
      lapack_int ldb;  //  leading dimension of 'b'
      double *a, *aref; //The array 'a' contains the matrix A
      double *b, *bref; //The array 'a' contains the matrix A

      /* Output parameters */
      double *tb, *tbref; //right-hand sides for the systems of equations.
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int *ipiv2, *ipiv2ref; // The pivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      dsytrs_aa_2stage_double_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~dsytrs_aa_2stage_double_parameters ();
};  /* end of dsytrs_aa_2stage_double_parameters  class definition */


/* Constructor dsytrs_aa_2stage_double_parameters definition */
dsytrs_aa_2stage_double_parameters:: dsytrs_aa_2stage_double_parameters (
					int matrix_layout_i, char uplo_i, lapack_int n_i,
								lapack_int nrhs_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
	nrhs = nrhs_i;
    uplo = uplo_i;
    lda = n;
    ltb = 4*n;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

    hModule = NULL;
    dModule = NULL;
	diff_tb = 0;
	diff_b = 0;
    diff_ipiv = 0;
	diff_ipiv2 = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n dsytrs_aa_2stage double:  n: %d, uplo: %c  ltb: %d lda: %d \n",
             n, uplo, ltb, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_double_buffer_pair(  &a, &aref, (n*lda));
    lapacke_gtest_alloc_double_buffer_pair(  &b, &bref, (n*nrhs));
    lapacke_gtest_alloc_double_buffer_pair( &tb, &tbref, ltb*ltb);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv2, &ipiv2ref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (tb==NULL) || (tbref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ||\
        (ipiv2==NULL) || (ipiv2ref==NULL) ){
       EXPECT_FALSE( true) << "dsytrs_aa_2stage_double_parameters object: malloc error.";
       sytrs_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

dsytrs_aa_2stage_double_parameters:: ~dsytrs_aa_2stage_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" dsytrs_aa_2stage_double_parameters object: destructor invoked. \n");
#endif
   sytrs_aa_2stage_free();
}


//  Test fixture class definition
class dsytrs_aa_2stage_test  : public  ::testing::Test {
public:
   dsytrs_aa_2stage_double_parameters  *dsytrs_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete dsytrs_aa_2stage_obj; }
};


void dsytrs_aa_2stage_test::SetUp(){

    /* LAPACKE_dsytrs_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,double *a ,lapack_int lda,
								 double *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_dsytrf_aa_2stage DSYTRF_AA_2STAGE;

    typedef int (*Fptr_NL_LAPACKE_dsytrs_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, double *a,
			lapack_int lda, double *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, double *b, lapack_int ldb );

    Fptr_NL_LAPACKE_dsytrs_aa_2stage DSYTRS_AA_2STAGE;

    dsytrs_aa_2stage_obj = new  dsytrs_aa_2stage_double_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    dsytrs_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsytrs_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsytrs_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsytrs_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    DSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_dsytrf_aa_2stage)dlsym(
							dsytrs_aa_2stage_obj->hModule, "LAPACKE_dsytrf_aa_2stage");
    ASSERT_TRUE(DSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_dsytrs_aa_2stage symbol";

    DSYTRS_AA_2STAGE = (Fptr_NL_LAPACKE_dsytrs_aa_2stage)dlsym(
							dsytrs_aa_2stage_obj->hModule, "LAPACKE_dsytrs_aa_2stage");
    ASSERT_TRUE(DSYTRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_dsytrs_aa_2stage symbol";

    dsytrs_aa_2stage_obj->inforef = DSYTRF_AA_2STAGE(
										dsytrs_aa_2stage_obj->matrix_layout,
										dsytrs_aa_2stage_obj->uplo,
										dsytrs_aa_2stage_obj->n,
										dsytrs_aa_2stage_obj->aref,
										dsytrs_aa_2stage_obj->lda,
										dsytrs_aa_2stage_obj->tbref,
										dsytrs_aa_2stage_obj->ltb,
										dsytrs_aa_2stage_obj->ipivref,
										dsytrs_aa_2stage_obj->ipiv2ref
										);

    dsytrs_aa_2stage_obj->info = LAPACKE_dsytrf_aa_2stage(
										dsytrs_aa_2stage_obj->matrix_layout,
										dsytrs_aa_2stage_obj->uplo,
										dsytrs_aa_2stage_obj->n,
										dsytrs_aa_2stage_obj->a,
										dsytrs_aa_2stage_obj->lda,
										dsytrs_aa_2stage_obj->tb,
										dsytrs_aa_2stage_obj->ltb,
										dsytrs_aa_2stage_obj->ipiv,
										dsytrs_aa_2stage_obj->ipiv2
										);

    /* Compute the Netlib-Lapacke's reference o/p */

    dsytrs_aa_2stage_obj->inforef = DSYTRS_AA_2STAGE(
										dsytrs_aa_2stage_obj->matrix_layout,
										dsytrs_aa_2stage_obj->uplo,
										dsytrs_aa_2stage_obj->n,
										dsytrs_aa_2stage_obj->nrhs,
										dsytrs_aa_2stage_obj->aref,
										dsytrs_aa_2stage_obj->lda,
										dsytrs_aa_2stage_obj->tbref,
										dsytrs_aa_2stage_obj->ltb,
										dsytrs_aa_2stage_obj->ipivref,
										dsytrs_aa_2stage_obj->ipiv2ref,
										dsytrs_aa_2stage_obj->bref,
										dsytrs_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    dsytrs_aa_2stage_obj->info = LAPACKE_dsytrs_aa_2stage(
										dsytrs_aa_2stage_obj->matrix_layout,
										dsytrs_aa_2stage_obj->uplo,
										dsytrs_aa_2stage_obj->n,
										dsytrs_aa_2stage_obj->nrhs,
										dsytrs_aa_2stage_obj->a,
										dsytrs_aa_2stage_obj->lda,
										dsytrs_aa_2stage_obj->tb,
										dsytrs_aa_2stage_obj->ltb,
										dsytrs_aa_2stage_obj->ipiv,
										dsytrs_aa_2stage_obj->ipiv2,
										dsytrs_aa_2stage_obj->b,
										dsytrs_aa_2stage_obj->ldb
										);
								 
								 
    if( dsytrs_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsytrs_aa_2stage is wrong\n",
                    dsytrs_aa_2stage_obj->info );
    }
    if( dsytrs_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsytrs_aa_2stage is wrong\n",
        dsytrs_aa_2stage_obj->inforef );
    }
	
	dsytrs_aa_2stage_obj->diff_tb =  computeDiff_d( dsytrs_aa_2stage_obj->ltb,
                           dsytrs_aa_2stage_obj->tb,
						   dsytrs_aa_2stage_obj->tbref );

	dsytrs_aa_2stage_obj->diff_b =  computeDiff_d( dsytrs_aa_2stage_obj->n*dsytrs_aa_2stage_obj->nrhs,
                           dsytrs_aa_2stage_obj->b,
						   dsytrs_aa_2stage_obj->bref );

    dsytrs_aa_2stage_obj->diff_ipiv = computeDiff_i( dsytrs_aa_2stage_obj->n,
							dsytrs_aa_2stage_obj->ipiv,
							dsytrs_aa_2stage_obj->ipivref);

    dsytrs_aa_2stage_obj->diff_ipiv2 = computeDiff_i( dsytrs_aa_2stage_obj->n,
							dsytrs_aa_2stage_obj->ipiv2,
							dsytrs_aa_2stage_obj->ipiv2ref);

}

TEST_F(dsytrs_aa_2stage_test, dsytrs_aa_2stage1) {
    EXPECT_NEAR(0.0, dsytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsytrs_aa_2stage_test, dsytrs_aa_2stage2) {
    EXPECT_NEAR(0.0, dsytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsytrs_aa_2stage_test, dsytrs_aa_2stage3) {
    EXPECT_NEAR(0.0, dsytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsytrs_aa_2stage_test, dsytrs_aa_2stage4) {
    EXPECT_NEAR(0.0, dsytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

/* Begin zsytrs_aa_2stage_dcomplex_parameters  class definition */
class zsytrs_aa_2stage_dcomplex_parameters{

   public:
      void *hModule,*dModule;
      double diff_tb, diff_b;
      int diff_ipiv, diff_ipiv2;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ltb;  //  dimension of 'tb'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_double *a, *aref; //The array 'a' contains the matrix A
      lapack_complex_double *b, *bref; //The array 'a' contains the matrix A

      /* Output parameters */
      lapack_complex_double *tb, *tbref; //right-hand sides for the systems of equations.
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int *ipiv2, *ipiv2ref; // The pivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      zsytrs_aa_2stage_dcomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~zsytrs_aa_2stage_dcomplex_parameters ();
};  /* end of zsytrs_aa_2stage_dcomplex_parameters  class definition */


/* Constructor zsytrs_aa_2stage_dcomplex_parameters definition */
zsytrs_aa_2stage_dcomplex_parameters:: zsytrs_aa_2stage_dcomplex_parameters (
					int matrix_layout_i, char uplo_i, lapack_int n_i,
								lapack_int nrhs_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
	nrhs = nrhs_i;
    uplo = uplo_i;
    lda = n;
    ltb = 4*n;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

    hModule = NULL;
    dModule = NULL;
	diff_tb = 0;
	diff_b = 0;
    diff_ipiv = 0;
	diff_ipiv2 = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n zsytrs_aa_2stage lapack_complex_double:  n: %d, uplo: %c  ltb: %d lda: %d \n",
             n, uplo, ltb, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(  &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(  &b, &bref, (n*nrhs));
    lapacke_gtest_alloc_lapack_dcomplex_buffer_pair( &tb, &tbref, ltb*ltb);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv2, &ipiv2ref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (tb==NULL) || (tbref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ||\
        (ipiv2==NULL) || (ipiv2ref==NULL) ){
       EXPECT_FALSE( true) << "zsytrs_aa_2stage_double_parameters object: malloc error.";
       sytrs_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

zsytrs_aa_2stage_dcomplex_parameters:: ~zsytrs_aa_2stage_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" zsytrs_aa_2stage_dcomplex_parameters object: destructor invoked. \n");
#endif
   sytrs_aa_2stage_free();
}


//  Test fixture class definition
class zsytrs_aa_2stage_test  : public  ::testing::Test {
public:
   zsytrs_aa_2stage_dcomplex_parameters  *zsytrs_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete zsytrs_aa_2stage_obj; }
};


void zsytrs_aa_2stage_test::SetUp(){

    /* LAPACKE_zsytrs_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,lapack_complex_double *a ,lapack_int lda,
								 lapack_complex_double *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_zsytrf_aa_2stage ZSYTRF_AA_2STAGE;

    typedef int (*Fptr_NL_LAPACKE_zsytrs_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, lapack_complex_double *a,
			lapack_int lda, lapack_complex_double *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, lapack_complex_double *b, lapack_int ldb );

    Fptr_NL_LAPACKE_zsytrs_aa_2stage ZSYTRS_AA_2STAGE;

    zsytrs_aa_2stage_obj = new  zsytrs_aa_2stage_dcomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    zsytrs_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsytrs_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsytrs_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsytrs_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_zsytrf_aa_2stage)dlsym(
							zsytrs_aa_2stage_obj->hModule, "LAPACKE_zsytrf_aa_2stage");
    ASSERT_TRUE(ZSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_zsytrs_aa_2stage symbol";

    ZSYTRS_AA_2STAGE = (Fptr_NL_LAPACKE_zsytrs_aa_2stage)dlsym(
							zsytrs_aa_2stage_obj->hModule, "LAPACKE_zsytrs_aa_2stage");
    ASSERT_TRUE(ZSYTRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_zsytrs_aa_2stage symbol";

    zsytrs_aa_2stage_obj->inforef = ZSYTRF_AA_2STAGE(
										zsytrs_aa_2stage_obj->matrix_layout,
										zsytrs_aa_2stage_obj->uplo,
										zsytrs_aa_2stage_obj->n,
										zsytrs_aa_2stage_obj->aref,
										zsytrs_aa_2stage_obj->lda,
										zsytrs_aa_2stage_obj->tbref,
										zsytrs_aa_2stage_obj->ltb,
										zsytrs_aa_2stage_obj->ipivref,
										zsytrs_aa_2stage_obj->ipiv2ref
										);

    zsytrs_aa_2stage_obj->info = LAPACKE_zsytrf_aa_2stage(
										zsytrs_aa_2stage_obj->matrix_layout,
										zsytrs_aa_2stage_obj->uplo,
										zsytrs_aa_2stage_obj->n,
										zsytrs_aa_2stage_obj->a,
										zsytrs_aa_2stage_obj->lda,
										zsytrs_aa_2stage_obj->tb,
										zsytrs_aa_2stage_obj->ltb,
										zsytrs_aa_2stage_obj->ipiv,
										zsytrs_aa_2stage_obj->ipiv2
										);

    /* Compute the Netlib-Lapacke's reference o/p */

    zsytrs_aa_2stage_obj->inforef = ZSYTRS_AA_2STAGE(
										zsytrs_aa_2stage_obj->matrix_layout,
										zsytrs_aa_2stage_obj->uplo,
										zsytrs_aa_2stage_obj->n,
										zsytrs_aa_2stage_obj->nrhs,
										zsytrs_aa_2stage_obj->aref,
										zsytrs_aa_2stage_obj->lda,
										zsytrs_aa_2stage_obj->tbref,
										zsytrs_aa_2stage_obj->ltb,
										zsytrs_aa_2stage_obj->ipivref,
										zsytrs_aa_2stage_obj->ipiv2ref,
										zsytrs_aa_2stage_obj->bref,
										zsytrs_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    zsytrs_aa_2stage_obj->info = LAPACKE_zsytrs_aa_2stage(
										zsytrs_aa_2stage_obj->matrix_layout,
										zsytrs_aa_2stage_obj->uplo,
										zsytrs_aa_2stage_obj->n,
										zsytrs_aa_2stage_obj->nrhs,
										zsytrs_aa_2stage_obj->a,
										zsytrs_aa_2stage_obj->lda,
										zsytrs_aa_2stage_obj->tb,
										zsytrs_aa_2stage_obj->ltb,
										zsytrs_aa_2stage_obj->ipiv,
										zsytrs_aa_2stage_obj->ipiv2,
										zsytrs_aa_2stage_obj->b,
										zsytrs_aa_2stage_obj->ldb
										);
								 
								 
    if( zsytrs_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsytrs_aa_2stage is wrong\n",
                    zsytrs_aa_2stage_obj->info );
    }
    if( zsytrs_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsytrs_aa_2stage is wrong\n",
        zsytrs_aa_2stage_obj->inforef );
    }
	
	zsytrs_aa_2stage_obj->diff_tb =  computeDiff_z( zsytrs_aa_2stage_obj->ltb,
                           zsytrs_aa_2stage_obj->tb,
						   zsytrs_aa_2stage_obj->tbref );

	zsytrs_aa_2stage_obj->diff_b =  computeDiff_z( zsytrs_aa_2stage_obj->n*zsytrs_aa_2stage_obj->nrhs,
                           zsytrs_aa_2stage_obj->b,
						   zsytrs_aa_2stage_obj->bref );

    zsytrs_aa_2stage_obj->diff_ipiv = computeDiff_i( zsytrs_aa_2stage_obj->n,
							zsytrs_aa_2stage_obj->ipiv,
							zsytrs_aa_2stage_obj->ipivref);

    zsytrs_aa_2stage_obj->diff_ipiv2 = computeDiff_i( zsytrs_aa_2stage_obj->n,
							zsytrs_aa_2stage_obj->ipiv2,
							zsytrs_aa_2stage_obj->ipiv2ref);

}

TEST_F(zsytrs_aa_2stage_test, zsytrs_aa_2stage1) {
    EXPECT_NEAR(0.0, zsytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsytrs_aa_2stage_test, zsytrs_aa_2stage2) {
    EXPECT_NEAR(0.0, zsytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsytrs_aa_2stage_test, zsytrs_aa_2stage3) {
    EXPECT_NEAR(0.0, zsytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsytrs_aa_2stage_test, zsytrs_aa_2stage4) {
    EXPECT_NEAR(0.0, zsytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}


/* Begin csytrs_aa_2stage_scomplex_parameters  class definition */
class csytrs_aa_2stage_scomplex_parameters{

   public:
      void *hModule,*dModule;
      float diff_tb, diff_b;
      int diff_ipiv, diff_ipiv2;

      /* Input parameters */
      int matrix_layout; //  storage layout LAPACK_ROW_MAJOR or LAPACK_COL_MAJOR
      char uplo; // Must be 'U' or 'L'
      lapack_int n; // The order of A; the number of rows in B
      lapack_int nrhs; // The number of right-hand sides
      lapack_int lda;  //  leading dimension of 'a'
      lapack_int ltb;  //  dimension of 'tb'
      lapack_int ldb;  //  leading dimension of 'b'
      lapack_complex_float *a, *aref; //The array 'a' contains the matrix A
      lapack_complex_float *b, *bref; //The array 'a' contains the matrix A

      /* Output parameters */
      lapack_complex_float *tb, *tbref; //right-hand sides for the systems of equations.
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int *ipiv2, *ipiv2ref; // The pivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      csytrs_aa_2stage_scomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~csytrs_aa_2stage_scomplex_parameters ();
};  /* end of csytrs_aa_2stage_scomplex_parameters  class definition */


/* Constructor csytrs_aa_2stage_scomplex_parameters definition */
csytrs_aa_2stage_scomplex_parameters:: csytrs_aa_2stage_scomplex_parameters (
					int matrix_layout_i, char uplo_i, lapack_int n_i,
								lapack_int nrhs_i) {

    matrix_layout = matrix_layout_i;
    n = n_i;
	nrhs = nrhs_i;
    uplo = uplo_i;
    lda = n;
    ltb = 4*n;

    if(matrix_layout==LAPACK_COL_MAJOR){
        ldb = n;
    }
    else if(matrix_layout==LAPACK_ROW_MAJOR){
        ldb = nrhs;
    }
    else
    {
        EXPECT_TRUE(false) << "matrix_layout invalid";
    }

    hModule = NULL;
    dModule = NULL;
	diff_tb = 0;
	diff_b = 0;
    diff_ipiv = 0;
	diff_ipiv2 = 0;

#if LAPACKE_TEST_VERBOSE
   printf(" \n csytrs_aa_2stage lapack_complex_float:  n: %d, uplo: %c  ltb: %d lda: %d \n",
             n, uplo, ltb, lda);
#endif

    /* Memory allocation of the buffers */
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair(  &a, &aref, (n*lda));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair(  &b, &bref, (n*nrhs));
    lapacke_gtest_alloc_lapack_scomplex_buffer_pair( &tb, &tbref, ltb*ltb);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv, &ipivref, n);
    lapacke_gtest_alloc_int_buffer_pair ( &ipiv2, &ipiv2ref, n);

    if( (a==NULL) || (aref==NULL) ||  \
        (tb==NULL) || (tbref==NULL) || \
        (ipiv==NULL) || (ipivref==NULL) ||\
        (ipiv2==NULL) || (ipiv2ref==NULL) ){
       EXPECT_FALSE( true) << "csytrs_aa_2stage_double_parameters object: malloc error.";
       sytrs_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

csytrs_aa_2stage_scomplex_parameters:: ~csytrs_aa_2stage_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" csytrs_aa_2stage_scomplex_parameters object: destructor invoked. \n");
#endif
   sytrs_aa_2stage_free();
}


//  Test fixture class definition
class csytrs_aa_2stage_test  : public  ::testing::Test {
public:
   csytrs_aa_2stage_scomplex_parameters  *csytrs_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete csytrs_aa_2stage_obj; }
};


void csytrs_aa_2stage_test::SetUp(){

    /* LAPACKE_csytrs_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,lapack_complex_float *a ,lapack_int lda,
								 lapack_complex_float *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_csytrf_aa_2stage CSYTRF_AA_2STAGE;

    typedef int (*Fptr_NL_LAPACKE_csytrs_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, lapack_complex_float *a,
			lapack_int lda, lapack_complex_float *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, lapack_complex_float *b, lapack_int ldb );

    Fptr_NL_LAPACKE_csytrs_aa_2stage CSYTRS_AA_2STAGE;

    csytrs_aa_2stage_obj = new  csytrs_aa_2stage_scomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    csytrs_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csytrs_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csytrs_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csytrs_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_csytrf_aa_2stage)dlsym(
							csytrs_aa_2stage_obj->hModule, "LAPACKE_csytrf_aa_2stage");
    ASSERT_TRUE(CSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_csytrs_aa_2stage symbol";

    CSYTRS_AA_2STAGE = (Fptr_NL_LAPACKE_csytrs_aa_2stage)dlsym(
							csytrs_aa_2stage_obj->hModule, "LAPACKE_csytrs_aa_2stage");
    ASSERT_TRUE(CSYTRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_csytrs_aa_2stage symbol";

    csytrs_aa_2stage_obj->inforef = CSYTRF_AA_2STAGE(
										csytrs_aa_2stage_obj->matrix_layout,
										csytrs_aa_2stage_obj->uplo,
										csytrs_aa_2stage_obj->n,
										csytrs_aa_2stage_obj->aref,
										csytrs_aa_2stage_obj->lda,
										csytrs_aa_2stage_obj->tbref,
										csytrs_aa_2stage_obj->ltb,
										csytrs_aa_2stage_obj->ipivref,
										csytrs_aa_2stage_obj->ipiv2ref
										);

    csytrs_aa_2stage_obj->info = LAPACKE_csytrf_aa_2stage(
										csytrs_aa_2stage_obj->matrix_layout,
										csytrs_aa_2stage_obj->uplo,
										csytrs_aa_2stage_obj->n,
										csytrs_aa_2stage_obj->a,
										csytrs_aa_2stage_obj->lda,
										csytrs_aa_2stage_obj->tb,
										csytrs_aa_2stage_obj->ltb,
										csytrs_aa_2stage_obj->ipiv,
										csytrs_aa_2stage_obj->ipiv2
										);

    /* Compute the Netlib-Lapacke's reference o/p */

    csytrs_aa_2stage_obj->inforef = LAPACKE_csytrs_aa_2stage(
										csytrs_aa_2stage_obj->matrix_layout,
										csytrs_aa_2stage_obj->uplo,
										csytrs_aa_2stage_obj->n,
										csytrs_aa_2stage_obj->nrhs,
										csytrs_aa_2stage_obj->aref,
										csytrs_aa_2stage_obj->lda,
										csytrs_aa_2stage_obj->tbref,
										csytrs_aa_2stage_obj->ltb,
										csytrs_aa_2stage_obj->ipivref,
										csytrs_aa_2stage_obj->ipiv2ref,
										csytrs_aa_2stage_obj->bref,
										csytrs_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    csytrs_aa_2stage_obj->info = LAPACKE_csytrs_aa_2stage(
										csytrs_aa_2stage_obj->matrix_layout,
										csytrs_aa_2stage_obj->uplo,
										csytrs_aa_2stage_obj->n,
										csytrs_aa_2stage_obj->nrhs,
										csytrs_aa_2stage_obj->a,
										csytrs_aa_2stage_obj->lda,
										csytrs_aa_2stage_obj->tb,
										csytrs_aa_2stage_obj->ltb,
										csytrs_aa_2stage_obj->ipiv,
										csytrs_aa_2stage_obj->ipiv2,
										csytrs_aa_2stage_obj->b,
										csytrs_aa_2stage_obj->ldb
										);
								 
								 
    if( csytrs_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csytrs_aa_2stage is wrong\n",
                    csytrs_aa_2stage_obj->info );
    }
    if( csytrs_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csytrs_aa_2stage is wrong\n",
        csytrs_aa_2stage_obj->inforef );
    }
	
	csytrs_aa_2stage_obj->diff_tb =  computeDiff_c( csytrs_aa_2stage_obj->ltb,
                           csytrs_aa_2stage_obj->tb,
						   csytrs_aa_2stage_obj->tbref );

	csytrs_aa_2stage_obj->diff_b =  computeDiff_c( csytrs_aa_2stage_obj->n*csytrs_aa_2stage_obj->nrhs,
                           csytrs_aa_2stage_obj->b,
						   csytrs_aa_2stage_obj->bref );

    csytrs_aa_2stage_obj->diff_ipiv = computeDiff_i( csytrs_aa_2stage_obj->n,
							csytrs_aa_2stage_obj->ipiv,
							csytrs_aa_2stage_obj->ipivref);

    csytrs_aa_2stage_obj->diff_ipiv2 = computeDiff_i( csytrs_aa_2stage_obj->n,
							csytrs_aa_2stage_obj->ipiv2,
							csytrs_aa_2stage_obj->ipiv2ref);

}

TEST_F(csytrs_aa_2stage_test, csytrs_aa_2stage1) {
    EXPECT_NEAR(0.0, csytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(csytrs_aa_2stage_test, csytrs_aa_2stage2) {
    EXPECT_NEAR(0.0, csytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(csytrs_aa_2stage_test, csytrs_aa_2stage3) {
    EXPECT_NEAR(0.0, csytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(csytrs_aa_2stage_test, csytrs_aa_2stage4) {
    EXPECT_NEAR(0.0, csytrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csytrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csytrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csytrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}
