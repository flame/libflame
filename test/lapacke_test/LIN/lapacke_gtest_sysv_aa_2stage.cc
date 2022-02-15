#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define sysv_aa_2stage_free() \
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

/* Begin ssysv_aa_2stage_float_parameters  class definition */
class ssysv_aa_2stage_float_parameters{

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
      ssysv_aa_2stage_float_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~ssysv_aa_2stage_float_parameters ();
};  /* end of ssysv_aa_2stage_float_parameters  class definition */


/* Constructor ssysv_aa_2stage_float_parameters definition */
ssysv_aa_2stage_float_parameters:: ssysv_aa_2stage_float_parameters (
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
   printf(" \n ssysv_aa_2stage float:  n: %d, uplo: %c  ltb: %d lda: %d \n",
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
       EXPECT_FALSE( true) << "ssysv_aa_2stage_double_parameters object: malloc error.";
       sysv_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_float_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

ssysv_aa_2stage_float_parameters:: ~ssysv_aa_2stage_float_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" ssysv_aa_2stage_float_parameters object: destructor invoked. \n");
#endif
   sysv_aa_2stage_free();
}


//  Test fixture class definition
class ssysv_aa_2stage_test  : public  ::testing::Test {
public:
   ssysv_aa_2stage_float_parameters  *ssysv_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete ssysv_aa_2stage_obj; }
};


void ssysv_aa_2stage_test::SetUp(){

    /* LAPACKE_ssysv_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_ssytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,float *a ,lapack_int lda,
								 float *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_ssytrf_aa_2stage SSYTRF_AA_2STAGE;

    typedef int (*Fptr_NL_LAPACKE_ssysv_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, float *a,
			lapack_int lda, float *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, float *b, lapack_int ldb );

    Fptr_NL_LAPACKE_ssysv_aa_2stage SSYTRS_AA_2STAGE;

    ssysv_aa_2stage_obj = new  ssysv_aa_2stage_float_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    ssysv_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ssysv_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ssysv_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ssysv_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    SSYTRS_AA_2STAGE = (Fptr_NL_LAPACKE_ssysv_aa_2stage)dlsym(
							ssysv_aa_2stage_obj->hModule, "LAPACKE_ssysv_aa_2stage");
    ASSERT_TRUE(SSYTRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_ssysv_aa_2stage symbol";

	/*
    SSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_ssytrf_aa_2stage)dlsym(
							ssysv_aa_2stage_obj->hModule, "LAPACKE_ssytrf_aa_2stage");
    ASSERT_TRUE(SSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_ssysv_aa_2stage symbol";

    ssysv_aa_2stage_obj->inforef = SSYTRF_AA_2STAGE(
										ssysv_aa_2stage_obj->matrix_layout,
										ssysv_aa_2stage_obj->uplo,
										ssysv_aa_2stage_obj->n,
										ssysv_aa_2stage_obj->aref,
										ssysv_aa_2stage_obj->lda,
										ssysv_aa_2stage_obj->tbref,
										ssysv_aa_2stage_obj->ltb,
										ssysv_aa_2stage_obj->ipivref,
										ssysv_aa_2stage_obj->ipiv2ref
										);  

    ssysv_aa_2stage_obj->info = LAPACKE_ssytrf_aa_2stage(
										ssysv_aa_2stage_obj->matrix_layout,
										ssysv_aa_2stage_obj->uplo,
										ssysv_aa_2stage_obj->n,
										ssysv_aa_2stage_obj->a,
										ssysv_aa_2stage_obj->lda,
										ssysv_aa_2stage_obj->tb,
										ssysv_aa_2stage_obj->ltb,
										ssysv_aa_2stage_obj->ipiv,
										ssysv_aa_2stage_obj->ipiv2
										);*/

    /* Compute the Netlib-Lapacke's reference o/p */

    ssysv_aa_2stage_obj->inforef = SSYTRS_AA_2STAGE(
										ssysv_aa_2stage_obj->matrix_layout,
										ssysv_aa_2stage_obj->uplo,
										ssysv_aa_2stage_obj->n,
										ssysv_aa_2stage_obj->nrhs,
										ssysv_aa_2stage_obj->aref,
										ssysv_aa_2stage_obj->lda,
										ssysv_aa_2stage_obj->tbref,
										ssysv_aa_2stage_obj->ltb,
										ssysv_aa_2stage_obj->ipivref,
										ssysv_aa_2stage_obj->ipiv2ref,
										ssysv_aa_2stage_obj->bref,
										ssysv_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    ssysv_aa_2stage_obj->info = LAPACKE_ssysv_aa_2stage(
										ssysv_aa_2stage_obj->matrix_layout,
										ssysv_aa_2stage_obj->uplo,
										ssysv_aa_2stage_obj->n,
										ssysv_aa_2stage_obj->nrhs,
										ssysv_aa_2stage_obj->a,
										ssysv_aa_2stage_obj->lda,
										ssysv_aa_2stage_obj->tb,
										ssysv_aa_2stage_obj->ltb,
										ssysv_aa_2stage_obj->ipiv,
										ssysv_aa_2stage_obj->ipiv2,
										ssysv_aa_2stage_obj->b,
										ssysv_aa_2stage_obj->ldb
										);
								 
								 
    if( ssysv_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_ssysv_aa_2stage is wrong\n",
                    ssysv_aa_2stage_obj->info );
    }
    if( ssysv_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ssysv_aa_2stage is wrong\n",
        ssysv_aa_2stage_obj->inforef );
    }
	
	ssysv_aa_2stage_obj->diff_tb =  computeDiff_s( ssysv_aa_2stage_obj->ltb,
                           ssysv_aa_2stage_obj->tb,
						   ssysv_aa_2stage_obj->tbref );

	ssysv_aa_2stage_obj->diff_b =  computeDiff_s( ssysv_aa_2stage_obj->n*ssysv_aa_2stage_obj->nrhs,
                           ssysv_aa_2stage_obj->b,
						   ssysv_aa_2stage_obj->bref );

    ssysv_aa_2stage_obj->diff_ipiv = computeDiff_i( ssysv_aa_2stage_obj->n,
							ssysv_aa_2stage_obj->ipiv,
							ssysv_aa_2stage_obj->ipivref);

    ssysv_aa_2stage_obj->diff_ipiv2 = computeDiff_i( ssysv_aa_2stage_obj->n,
							ssysv_aa_2stage_obj->ipiv2,
							ssysv_aa_2stage_obj->ipiv2ref);

}

TEST_F(ssysv_aa_2stage_test, ssysv_aa_2stage1) {
    EXPECT_NEAR(0.0, ssysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssysv_aa_2stage_test, ssysv_aa_2stage2) {
    EXPECT_NEAR(0.0, ssysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssysv_aa_2stage_test, ssysv_aa_2stage3) {
    EXPECT_NEAR(0.0, ssysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(ssysv_aa_2stage_test, ssysv_aa_2stage4) {
    EXPECT_NEAR(0.0, ssysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, ssysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, ssysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, ssysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}


/* Begin dsysv_aa_2stage_double_parameters  class definition */
class dsysv_aa_2stage_double_parameters{

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
      dsysv_aa_2stage_double_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~dsysv_aa_2stage_double_parameters ();
};  /* end of dsysv_aa_2stage_double_parameters  class definition */


/* Constructor dsysv_aa_2stage_double_parameters definition */
dsysv_aa_2stage_double_parameters:: dsysv_aa_2stage_double_parameters (
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
   printf(" \n dsysv_aa_2stage double:  n: %d, uplo: %c  ltb: %d lda: %d \n",
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
       EXPECT_FALSE( true) << "dsysv_aa_2stage_double_parameters object: malloc error.";
       sysv_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_double_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

dsysv_aa_2stage_double_parameters:: ~dsysv_aa_2stage_double_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" dsysv_aa_2stage_double_parameters object: destructor invoked. \n");
#endif
   sysv_aa_2stage_free();
}


//  Test fixture class definition
class dsysv_aa_2stage_test  : public  ::testing::Test {
public:
   dsysv_aa_2stage_double_parameters  *dsysv_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete dsysv_aa_2stage_obj; }
};


void dsysv_aa_2stage_test::SetUp(){

    /* LAPACKE_dsysv_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_dsytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,double *a ,lapack_int lda,
								 double *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_dsytrf_aa_2stage DSYTRF_AA_2STAGE;

    typedef int (*Fptr_NL_LAPACKE_dsysv_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, double *a,
			lapack_int lda, double *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, double *b, lapack_int ldb );

    Fptr_NL_LAPACKE_dsysv_aa_2stage DSYTRS_AA_2STAGE;

    dsysv_aa_2stage_obj = new  dsysv_aa_2stage_double_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    dsysv_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dsysv_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dsysv_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dsysv_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    DSYTRS_AA_2STAGE = (Fptr_NL_LAPACKE_dsysv_aa_2stage)dlsym(
							dsysv_aa_2stage_obj->hModule, "LAPACKE_dsysv_aa_2stage");
    ASSERT_TRUE(DSYTRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_dsysv_aa_2stage symbol";
/*
    DSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_dsytrf_aa_2stage)dlsym(
							dsysv_aa_2stage_obj->hModule, "LAPACKE_dsytrf_aa_2stage");
    ASSERT_TRUE(DSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_dsysv_aa_2stage symbol";

    dsysv_aa_2stage_obj->inforef = DSYTRF_AA_2STAGE(
										dsysv_aa_2stage_obj->matrix_layout,
										dsysv_aa_2stage_obj->uplo,
										dsysv_aa_2stage_obj->n,
										dsysv_aa_2stage_obj->aref,
										dsysv_aa_2stage_obj->lda,
										dsysv_aa_2stage_obj->tbref,
										dsysv_aa_2stage_obj->ltb,
										dsysv_aa_2stage_obj->ipivref,
										dsysv_aa_2stage_obj->ipiv2ref
										);

    dsysv_aa_2stage_obj->info = LAPACKE_dsytrf_aa_2stage(
										dsysv_aa_2stage_obj->matrix_layout,
										dsysv_aa_2stage_obj->uplo,
										dsysv_aa_2stage_obj->n,
										dsysv_aa_2stage_obj->a,
										dsysv_aa_2stage_obj->lda,
										dsysv_aa_2stage_obj->tb,
										dsysv_aa_2stage_obj->ltb,
										dsysv_aa_2stage_obj->ipiv,
										dsysv_aa_2stage_obj->ipiv2
										);

     Compute the Netlib-Lapacke's reference o/p */

    dsysv_aa_2stage_obj->inforef = DSYTRS_AA_2STAGE(
										dsysv_aa_2stage_obj->matrix_layout,
										dsysv_aa_2stage_obj->uplo,
										dsysv_aa_2stage_obj->n,
										dsysv_aa_2stage_obj->nrhs,
										dsysv_aa_2stage_obj->aref,
										dsysv_aa_2stage_obj->lda,
										dsysv_aa_2stage_obj->tbref,
										dsysv_aa_2stage_obj->ltb,
										dsysv_aa_2stage_obj->ipivref,
										dsysv_aa_2stage_obj->ipiv2ref,
										dsysv_aa_2stage_obj->bref,
										dsysv_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    dsysv_aa_2stage_obj->info = LAPACKE_dsysv_aa_2stage(
										dsysv_aa_2stage_obj->matrix_layout,
										dsysv_aa_2stage_obj->uplo,
										dsysv_aa_2stage_obj->n,
										dsysv_aa_2stage_obj->nrhs,
										dsysv_aa_2stage_obj->a,
										dsysv_aa_2stage_obj->lda,
										dsysv_aa_2stage_obj->tb,
										dsysv_aa_2stage_obj->ltb,
										dsysv_aa_2stage_obj->ipiv,
										dsysv_aa_2stage_obj->ipiv2,
										dsysv_aa_2stage_obj->b,
										dsysv_aa_2stage_obj->ldb
										);
								 
								 
    if( dsysv_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_dsysv_aa_2stage is wrong\n",
                    dsysv_aa_2stage_obj->info );
    }
    if( dsysv_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dsysv_aa_2stage is wrong\n",
        dsysv_aa_2stage_obj->inforef );
    }
	
	dsysv_aa_2stage_obj->diff_tb =  computeDiff_d( dsysv_aa_2stage_obj->ltb,
                           dsysv_aa_2stage_obj->tb,
						   dsysv_aa_2stage_obj->tbref );

	dsysv_aa_2stage_obj->diff_b =  computeDiff_d( dsysv_aa_2stage_obj->n*dsysv_aa_2stage_obj->nrhs,
                           dsysv_aa_2stage_obj->b,
						   dsysv_aa_2stage_obj->bref );

    dsysv_aa_2stage_obj->diff_ipiv = computeDiff_i( dsysv_aa_2stage_obj->n,
							dsysv_aa_2stage_obj->ipiv,
							dsysv_aa_2stage_obj->ipivref);

    dsysv_aa_2stage_obj->diff_ipiv2 = computeDiff_i( dsysv_aa_2stage_obj->n,
							dsysv_aa_2stage_obj->ipiv2,
							dsysv_aa_2stage_obj->ipiv2ref);

}

TEST_F(dsysv_aa_2stage_test, dsysv_aa_2stage1) {
    EXPECT_NEAR(0.0, dsysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsysv_aa_2stage_test, dsysv_aa_2stage2) {
    EXPECT_NEAR(0.0, dsysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsysv_aa_2stage_test, dsysv_aa_2stage3) {
    EXPECT_NEAR(0.0, dsysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(dsysv_aa_2stage_test, dsysv_aa_2stage4) {
    EXPECT_NEAR(0.0, dsysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, dsysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, dsysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, dsysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

/* Begin zsysv_aa_2stage_dcomplex_parameters  class definition */
class zsysv_aa_2stage_dcomplex_parameters{

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
      zsysv_aa_2stage_dcomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~zsysv_aa_2stage_dcomplex_parameters ();
};  /* end of zsysv_aa_2stage_dcomplex_parameters  class definition */


/* Constructor zsysv_aa_2stage_dcomplex_parameters definition */
zsysv_aa_2stage_dcomplex_parameters:: zsysv_aa_2stage_dcomplex_parameters (
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
   printf(" \n zsysv_aa_2stage lapack_complex_double:  n: %d, uplo: %c  ltb: %d lda: %d \n",
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
       EXPECT_FALSE( true) << "zsysv_aa_2stage_double_parameters object: malloc error.";
       sysv_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

zsysv_aa_2stage_dcomplex_parameters:: ~zsysv_aa_2stage_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" zsysv_aa_2stage_dcomplex_parameters object: destructor invoked. \n");
#endif
   sysv_aa_2stage_free();
}


//  Test fixture class definition
class zsysv_aa_2stage_test  : public  ::testing::Test {
public:
   zsysv_aa_2stage_dcomplex_parameters  *zsysv_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete zsysv_aa_2stage_obj; }
};


void zsysv_aa_2stage_test::SetUp(){

    /* LAPACKE_zsysv_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zsytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,lapack_complex_double *a ,lapack_int lda,
								 lapack_complex_double *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_zsytrf_aa_2stage ZSYTRF_AA_2STAGE;

    typedef int (*Fptr_NL_LAPACKE_zsysv_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, lapack_complex_double *a,
			lapack_int lda, lapack_complex_double *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, lapack_complex_double *b, lapack_int ldb );

    Fptr_NL_LAPACKE_zsysv_aa_2stage ZSYTRS_AA_2STAGE;

    zsysv_aa_2stage_obj = new  zsysv_aa_2stage_dcomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    zsysv_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsysv_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsysv_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsysv_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZSYTRS_AA_2STAGE = (Fptr_NL_LAPACKE_zsysv_aa_2stage)dlsym(
							zsysv_aa_2stage_obj->hModule, "LAPACKE_zsysv_aa_2stage");
    ASSERT_TRUE(ZSYTRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_zsysv_aa_2stage symbol";
    /*
    ZSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_zsytrf_aa_2stage)dlsym(
							zsysv_aa_2stage_obj->hModule, "LAPACKE_zsytrf_aa_2stage");
    ASSERT_TRUE(ZSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_zsysv_aa_2stage symbol";

    zsysv_aa_2stage_obj->inforef = ZSYTRF_AA_2STAGE(
										zsysv_aa_2stage_obj->matrix_layout,
										zsysv_aa_2stage_obj->uplo,
										zsysv_aa_2stage_obj->n,
										zsysv_aa_2stage_obj->aref,
										zsysv_aa_2stage_obj->lda,
										zsysv_aa_2stage_obj->tbref,
										zsysv_aa_2stage_obj->ltb,
										zsysv_aa_2stage_obj->ipivref,
										zsysv_aa_2stage_obj->ipiv2ref
										);

    zsysv_aa_2stage_obj->info = LAPACKE_zsytrf_aa_2stage(
										zsysv_aa_2stage_obj->matrix_layout,
										zsysv_aa_2stage_obj->uplo,
										zsysv_aa_2stage_obj->n,
										zsysv_aa_2stage_obj->a,
										zsysv_aa_2stage_obj->lda,
										zsysv_aa_2stage_obj->tb,
										zsysv_aa_2stage_obj->ltb,
										zsysv_aa_2stage_obj->ipiv,
										zsysv_aa_2stage_obj->ipiv2
										);

 Compute the Netlib-Lapacke's reference o/p */

    zsysv_aa_2stage_obj->inforef = ZSYTRS_AA_2STAGE(
										zsysv_aa_2stage_obj->matrix_layout,
										zsysv_aa_2stage_obj->uplo,
										zsysv_aa_2stage_obj->n,
										zsysv_aa_2stage_obj->nrhs,
										zsysv_aa_2stage_obj->aref,
										zsysv_aa_2stage_obj->lda,
										zsysv_aa_2stage_obj->tbref,
										zsysv_aa_2stage_obj->ltb,
										zsysv_aa_2stage_obj->ipivref,
										zsysv_aa_2stage_obj->ipiv2ref,
										zsysv_aa_2stage_obj->bref,
										zsysv_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    zsysv_aa_2stage_obj->info = LAPACKE_zsysv_aa_2stage(
										zsysv_aa_2stage_obj->matrix_layout,
										zsysv_aa_2stage_obj->uplo,
										zsysv_aa_2stage_obj->n,
										zsysv_aa_2stage_obj->nrhs,
										zsysv_aa_2stage_obj->a,
										zsysv_aa_2stage_obj->lda,
										zsysv_aa_2stage_obj->tb,
										zsysv_aa_2stage_obj->ltb,
										zsysv_aa_2stage_obj->ipiv,
										zsysv_aa_2stage_obj->ipiv2,
										zsysv_aa_2stage_obj->b,
										zsysv_aa_2stage_obj->ldb
										);
								 
								 
    if( zsysv_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zsysv_aa_2stage is wrong\n",
                    zsysv_aa_2stage_obj->info );
    }
    if( zsysv_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsysv_aa_2stage is wrong\n",
        zsysv_aa_2stage_obj->inforef );
    }
	
	zsysv_aa_2stage_obj->diff_tb =  computeDiff_z( zsysv_aa_2stage_obj->ltb,
                           zsysv_aa_2stage_obj->tb,
						   zsysv_aa_2stage_obj->tbref );

	zsysv_aa_2stage_obj->diff_b =  computeDiff_z( zsysv_aa_2stage_obj->n*zsysv_aa_2stage_obj->nrhs,
                           zsysv_aa_2stage_obj->b,
						   zsysv_aa_2stage_obj->bref );

    zsysv_aa_2stage_obj->diff_ipiv = computeDiff_i( zsysv_aa_2stage_obj->n,
							zsysv_aa_2stage_obj->ipiv,
							zsysv_aa_2stage_obj->ipivref);

    zsysv_aa_2stage_obj->diff_ipiv2 = computeDiff_i( zsysv_aa_2stage_obj->n,
							zsysv_aa_2stage_obj->ipiv2,
							zsysv_aa_2stage_obj->ipiv2ref);

}

TEST_F(zsysv_aa_2stage_test, zsysv_aa_2stage1) {
    EXPECT_NEAR(0.0, zsysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsysv_aa_2stage_test, zsysv_aa_2stage2) {
    EXPECT_NEAR(0.0, zsysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsysv_aa_2stage_test, zsysv_aa_2stage3) {
    EXPECT_NEAR(0.0, zsysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zsysv_aa_2stage_test, zsysv_aa_2stage4) {
    EXPECT_NEAR(0.0, zsysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zsysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zsysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zsysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}


/* Begin csysv_aa_2stage_scomplex_parameters  class definition */
class csysv_aa_2stage_scomplex_parameters{

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
      csysv_aa_2stage_scomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~csysv_aa_2stage_scomplex_parameters ();
};  /* end of csysv_aa_2stage_scomplex_parameters  class definition */


/* Constructor csysv_aa_2stage_scomplex_parameters definition */
csysv_aa_2stage_scomplex_parameters:: csysv_aa_2stage_scomplex_parameters (
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
   printf(" \n csysv_aa_2stage lapack_complex_float:  n: %d, uplo: %c  ltb: %d lda: %d \n",
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
       EXPECT_FALSE( true) << "csysv_aa_2stage_double_parameters object: malloc error.";
       sysv_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

csysv_aa_2stage_scomplex_parameters:: ~csysv_aa_2stage_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" csysv_aa_2stage_scomplex_parameters object: destructor invoked. \n");
#endif
   sysv_aa_2stage_free();
}


//  Test fixture class definition
class csysv_aa_2stage_test  : public  ::testing::Test {
public:
   csysv_aa_2stage_scomplex_parameters  *csysv_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete csysv_aa_2stage_obj; }
};


void csysv_aa_2stage_test::SetUp(){

    /* LAPACKE_csysv_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_csytrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,lapack_complex_float *a ,lapack_int lda,
								 lapack_complex_float *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_csytrf_aa_2stage CSYTRF_AA_2STAGE;

    typedef int (*Fptr_NL_LAPACKE_csysv_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, lapack_complex_float *a,
			lapack_int lda, lapack_complex_float *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, lapack_complex_float *b, lapack_int ldb );

    Fptr_NL_LAPACKE_csysv_aa_2stage CSYTRS_AA_2STAGE;

    csysv_aa_2stage_obj = new  csysv_aa_2stage_scomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    csysv_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csysv_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csysv_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csysv_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    CSYTRS_AA_2STAGE = (Fptr_NL_LAPACKE_csysv_aa_2stage)dlsym(
							csysv_aa_2stage_obj->hModule, "LAPACKE_csysv_aa_2stage");
    ASSERT_TRUE(CSYTRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_csysv_aa_2stage symbol";
    /*
     CSYTRF_AA_2STAGE = (Fptr_NL_LAPACKE_csytrf_aa_2stage)dlsym(
							csysv_aa_2stage_obj->hModule, "LAPACKE_csytrf_aa_2stage");
    ASSERT_TRUE(CSYTRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_csysv_aa_2stage symbol";
   csysv_aa_2stage_obj->inforef = CSYTRF_AA_2STAGE(
										csysv_aa_2stage_obj->matrix_layout,
										csysv_aa_2stage_obj->uplo,
										csysv_aa_2stage_obj->n,
										csysv_aa_2stage_obj->aref,
										csysv_aa_2stage_obj->lda,
										csysv_aa_2stage_obj->tbref,
										csysv_aa_2stage_obj->ltb,
										csysv_aa_2stage_obj->ipivref,
										csysv_aa_2stage_obj->ipiv2ref
										);

    csysv_aa_2stage_obj->info = LAPACKE_csytrf_aa_2stage(
										csysv_aa_2stage_obj->matrix_layout,
										csysv_aa_2stage_obj->uplo,
										csysv_aa_2stage_obj->n,
										csysv_aa_2stage_obj->a,
										csysv_aa_2stage_obj->lda,
										csysv_aa_2stage_obj->tb,
										csysv_aa_2stage_obj->ltb,
										csysv_aa_2stage_obj->ipiv,
										csysv_aa_2stage_obj->ipiv2
										);

 Compute the Netlib-Lapacke's reference o/p */

    csysv_aa_2stage_obj->inforef = LAPACKE_csysv_aa_2stage(
										csysv_aa_2stage_obj->matrix_layout,
										csysv_aa_2stage_obj->uplo,
										csysv_aa_2stage_obj->n,
										csysv_aa_2stage_obj->nrhs,
										csysv_aa_2stage_obj->aref,
										csysv_aa_2stage_obj->lda,
										csysv_aa_2stage_obj->tbref,
										csysv_aa_2stage_obj->ltb,
										csysv_aa_2stage_obj->ipivref,
										csysv_aa_2stage_obj->ipiv2ref,
										csysv_aa_2stage_obj->bref,
										csysv_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    csysv_aa_2stage_obj->info = LAPACKE_csysv_aa_2stage(
										csysv_aa_2stage_obj->matrix_layout,
										csysv_aa_2stage_obj->uplo,
										csysv_aa_2stage_obj->n,
										csysv_aa_2stage_obj->nrhs,
										csysv_aa_2stage_obj->a,
										csysv_aa_2stage_obj->lda,
										csysv_aa_2stage_obj->tb,
										csysv_aa_2stage_obj->ltb,
										csysv_aa_2stage_obj->ipiv,
										csysv_aa_2stage_obj->ipiv2,
										csysv_aa_2stage_obj->b,
										csysv_aa_2stage_obj->ldb
										);
								 
								 
    if( csysv_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_csysv_aa_2stage is wrong\n",
                    csysv_aa_2stage_obj->info );
    }
    if( csysv_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csysv_aa_2stage is wrong\n",
        csysv_aa_2stage_obj->inforef );
    }
	
	csysv_aa_2stage_obj->diff_tb =  computeDiff_c( csysv_aa_2stage_obj->ltb,
                           csysv_aa_2stage_obj->tb,
						   csysv_aa_2stage_obj->tbref );

	csysv_aa_2stage_obj->diff_b =  computeDiff_c( csysv_aa_2stage_obj->n*csysv_aa_2stage_obj->nrhs,
                           csysv_aa_2stage_obj->b,
						   csysv_aa_2stage_obj->bref );

    csysv_aa_2stage_obj->diff_ipiv = computeDiff_i( csysv_aa_2stage_obj->n,
							csysv_aa_2stage_obj->ipiv,
							csysv_aa_2stage_obj->ipivref);

    csysv_aa_2stage_obj->diff_ipiv2 = computeDiff_i( csysv_aa_2stage_obj->n,
							csysv_aa_2stage_obj->ipiv2,
							csysv_aa_2stage_obj->ipiv2ref);

}

TEST_F(csysv_aa_2stage_test, csysv_aa_2stage1) {
    EXPECT_NEAR(0.0, csysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(csysv_aa_2stage_test, csysv_aa_2stage2) {
    EXPECT_NEAR(0.0, csysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(csysv_aa_2stage_test, csysv_aa_2stage3) {
    EXPECT_NEAR(0.0, csysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(csysv_aa_2stage_test, csysv_aa_2stage4) {
    EXPECT_NEAR(0.0, csysv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, csysv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, csysv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, csysv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}
