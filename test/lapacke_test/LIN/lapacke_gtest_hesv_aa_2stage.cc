#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define hesv_aa_2stage_free() \
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

/* Begin zhesv_aa_2stage_dcomplex_parameters  class definition */
class zhesv_aa_2stage_dcomplex_parameters{

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
      zhesv_aa_2stage_dcomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~zhesv_aa_2stage_dcomplex_parameters ();
};  /* end of zhesv_aa_2stage_dcomplex_parameters  class definition */


/* Constructor zhesv_aa_2stage_dcomplex_parameters definition */
zhesv_aa_2stage_dcomplex_parameters:: zhesv_aa_2stage_dcomplex_parameters (
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
   printf(" \n zhesv_aa_2stage lapack_complex_double:  n: %d, uplo: %c  ltb: %d lda: %d \n",
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
       EXPECT_FALSE( true) << "zhesv_aa_2stage_double_parameters object: malloc error.";
       hesv_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

zhesv_aa_2stage_dcomplex_parameters:: ~zhesv_aa_2stage_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" zhesv_aa_2stage_dcomplex_parameters object: destructor invoked. \n");
#endif
   hesv_aa_2stage_free();
}


//  Test fixture class definition
class zhesv_aa_2stage_test  : public  ::testing::Test {
public:
   zhesv_aa_2stage_dcomplex_parameters  *zhesv_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete zhesv_aa_2stage_obj; }
};


void zhesv_aa_2stage_test::SetUp(){

    /* LAPACKE_zhesv_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zhesv_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, lapack_complex_double *a,
			lapack_int lda, lapack_complex_double *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, lapack_complex_double *b, lapack_int ldb );

    Fptr_NL_LAPACKE_zhesv_aa_2stage ZHETRS_AA_2STAGE;

    zhesv_aa_2stage_obj = new  zhesv_aa_2stage_dcomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    zhesv_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhesv_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhesv_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhesv_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHETRS_AA_2STAGE = (Fptr_NL_LAPACKE_zhesv_aa_2stage)dlsym(
							zhesv_aa_2stage_obj->hModule, "LAPACKE_zhesv_aa_2stage");
    ASSERT_TRUE(ZHETRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_zhesv_aa_2stage symbol";

    /* Compute the Netlib-Lapacke's reference o/p */

    zhesv_aa_2stage_obj->inforef = ZHETRS_AA_2STAGE(
										zhesv_aa_2stage_obj->matrix_layout,
										zhesv_aa_2stage_obj->uplo,
										zhesv_aa_2stage_obj->n,
										zhesv_aa_2stage_obj->nrhs,
										zhesv_aa_2stage_obj->aref,
										zhesv_aa_2stage_obj->lda,
										zhesv_aa_2stage_obj->tbref,
										zhesv_aa_2stage_obj->ltb,
										zhesv_aa_2stage_obj->ipivref,
										zhesv_aa_2stage_obj->ipiv2ref,
										zhesv_aa_2stage_obj->bref,
										zhesv_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    zhesv_aa_2stage_obj->info = LAPACKE_zhesv_aa_2stage(
										zhesv_aa_2stage_obj->matrix_layout,
										zhesv_aa_2stage_obj->uplo,
										zhesv_aa_2stage_obj->n,
										zhesv_aa_2stage_obj->nrhs,
										zhesv_aa_2stage_obj->a,
										zhesv_aa_2stage_obj->lda,
										zhesv_aa_2stage_obj->tb,
										zhesv_aa_2stage_obj->ltb,
										zhesv_aa_2stage_obj->ipiv,
										zhesv_aa_2stage_obj->ipiv2,
										zhesv_aa_2stage_obj->b,
										zhesv_aa_2stage_obj->ldb
										);
								 
								 
    if( zhesv_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhesv_aa_2stage is wrong\n",
                    zhesv_aa_2stage_obj->info );
    }
    if( zhesv_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhesv_aa_2stage is wrong\n",
        zhesv_aa_2stage_obj->inforef );
    }
	
	zhesv_aa_2stage_obj->diff_tb =  computeDiff_z( zhesv_aa_2stage_obj->ltb,
                           zhesv_aa_2stage_obj->tb,
						   zhesv_aa_2stage_obj->tbref );

	zhesv_aa_2stage_obj->diff_b =  computeDiff_z( zhesv_aa_2stage_obj->n*zhesv_aa_2stage_obj->nrhs,
                           zhesv_aa_2stage_obj->b,
						   zhesv_aa_2stage_obj->bref );

    zhesv_aa_2stage_obj->diff_ipiv = computeDiff_i( zhesv_aa_2stage_obj->n,
							zhesv_aa_2stage_obj->ipiv,
							zhesv_aa_2stage_obj->ipivref);

    zhesv_aa_2stage_obj->diff_ipiv2 = computeDiff_i( zhesv_aa_2stage_obj->n,
							zhesv_aa_2stage_obj->ipiv2,
							zhesv_aa_2stage_obj->ipiv2ref);

}

TEST_F(zhesv_aa_2stage_test, zhesv_aa_2stage1) {
    EXPECT_NEAR(0.0, zhesv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zhesv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zhesv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zhesv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zhesv_aa_2stage_test, zhesv_aa_2stage2) {
    EXPECT_NEAR(0.0, zhesv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zhesv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zhesv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zhesv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zhesv_aa_2stage_test, zhesv_aa_2stage3) {
    EXPECT_NEAR(0.0, zhesv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zhesv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zhesv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zhesv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zhesv_aa_2stage_test, zhesv_aa_2stage4) {
    EXPECT_NEAR(0.0, zhesv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zhesv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zhesv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zhesv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

/* Begin chesv_aa_2stage_scomplex_parameters  class definition */
class chesv_aa_2stage_scomplex_parameters{

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
      lapack_complex_float *b, *bref; //The array 'b' contains the matrix B

      /* Output parameters */
      lapack_complex_float *tb, *tbref; //right-hand sides for the systems of equations.
      lapack_int *ipiv, *ipivref; // The pivot indices
      lapack_int *ipiv2, *ipiv2ref; // The pivot indices

      /* Return Values */
      lapack_int info, inforef;

   public:
      chesv_aa_2stage_scomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~chesv_aa_2stage_scomplex_parameters ();
};  /* end of chesv_aa_2stage_scomplex_parameters  class definition */


/* Constructor chesv_aa_2stage_scomplex_parameters definition */
chesv_aa_2stage_scomplex_parameters:: chesv_aa_2stage_scomplex_parameters (
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
   printf(" \n chesv_aa_2stage lapack_complex_float:  n: %d, uplo: %c  ltb: %d lda: %d \n",
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
       EXPECT_FALSE( true) << "chesv_aa_2stage_double_parameters object: malloc error.";
       hesv_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

chesv_aa_2stage_scomplex_parameters:: ~chesv_aa_2stage_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" chesv_aa_2stage_scomplex_parameters object: destructor invoked. \n");
#endif
   hesv_aa_2stage_free();
}


//  Test fixture class definition
class chesv_aa_2stage_test  : public  ::testing::Test {
public:
   chesv_aa_2stage_scomplex_parameters  *chesv_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete chesv_aa_2stage_obj; }
};


void chesv_aa_2stage_test::SetUp(){

    /* LAPACKE_chesv_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_chesv_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, lapack_complex_float *a,
			lapack_int lda, lapack_complex_float *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, lapack_complex_float *b, lapack_int ldb );

    Fptr_NL_LAPACKE_chesv_aa_2stage CHETRS_AA_2STAGE;

    chesv_aa_2stage_obj = new  chesv_aa_2stage_scomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    chesv_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chesv_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chesv_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chesv_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHETRS_AA_2STAGE = (Fptr_NL_LAPACKE_chesv_aa_2stage)dlsym(
							chesv_aa_2stage_obj->hModule, "LAPACKE_chesv_aa_2stage");
    ASSERT_TRUE(CHETRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_chesv_aa_2stage symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    chesv_aa_2stage_obj->inforef = CHETRS_AA_2STAGE(
										chesv_aa_2stage_obj->matrix_layout,
										chesv_aa_2stage_obj->uplo,
										chesv_aa_2stage_obj->n,
										chesv_aa_2stage_obj->nrhs,
										chesv_aa_2stage_obj->aref,
										chesv_aa_2stage_obj->lda,
										chesv_aa_2stage_obj->tbref,
										chesv_aa_2stage_obj->ltb,
										chesv_aa_2stage_obj->ipivref,
										chesv_aa_2stage_obj->ipiv2ref,
										chesv_aa_2stage_obj->bref,
										chesv_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    chesv_aa_2stage_obj->info = LAPACKE_chesv_aa_2stage(
										chesv_aa_2stage_obj->matrix_layout,
										chesv_aa_2stage_obj->uplo,
										chesv_aa_2stage_obj->n,
										chesv_aa_2stage_obj->nrhs,
										chesv_aa_2stage_obj->a,
										chesv_aa_2stage_obj->lda,
										chesv_aa_2stage_obj->tb,
										chesv_aa_2stage_obj->ltb,
										chesv_aa_2stage_obj->ipiv,
										chesv_aa_2stage_obj->ipiv2,
										chesv_aa_2stage_obj->b,
										chesv_aa_2stage_obj->ldb
										);
								 
								 
    if( chesv_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chesv_aa_2stage is wrong\n",
                    chesv_aa_2stage_obj->info );
    }
    if( chesv_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chesv_aa_2stage is wrong\n",
        chesv_aa_2stage_obj->inforef );
    }
	
	chesv_aa_2stage_obj->diff_tb =  computeDiff_c( chesv_aa_2stage_obj->ltb,
                           chesv_aa_2stage_obj->tb,
						   chesv_aa_2stage_obj->tbref );

	chesv_aa_2stage_obj->diff_b =  computeDiff_c( chesv_aa_2stage_obj->n*chesv_aa_2stage_obj->nrhs,
                           chesv_aa_2stage_obj->b,
						   chesv_aa_2stage_obj->bref );

    chesv_aa_2stage_obj->diff_ipiv = computeDiff_i( chesv_aa_2stage_obj->n,
							chesv_aa_2stage_obj->ipiv,
							chesv_aa_2stage_obj->ipivref);

    chesv_aa_2stage_obj->diff_ipiv2 = computeDiff_i( chesv_aa_2stage_obj->n,
							chesv_aa_2stage_obj->ipiv2,
							chesv_aa_2stage_obj->ipiv2ref);

}

TEST_F(chesv_aa_2stage_test, chesv_aa_2stage1) {
    EXPECT_NEAR(0.0, chesv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, chesv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chesv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chesv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(chesv_aa_2stage_test, chesv_aa_2stage2) {
    EXPECT_NEAR(0.0, chesv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, chesv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chesv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chesv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(chesv_aa_2stage_test, chesv_aa_2stage3) {
    EXPECT_NEAR(0.0, chesv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, chesv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chesv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chesv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(chesv_aa_2stage_test, chesv_aa_2stage4) {
    EXPECT_NEAR(0.0, chesv_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, chesv_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chesv_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chesv_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}
