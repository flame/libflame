#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define hetrs_aa_2stage_free() \
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

/* Begin zhetrs_aa_2stage_dcomplex_parameters  class definition */
class zhetrs_aa_2stage_dcomplex_parameters{

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
      zhetrs_aa_2stage_dcomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~zhetrs_aa_2stage_dcomplex_parameters ();
};  /* end of zhetrs_aa_2stage_dcomplex_parameters  class definition */


/* Constructor zhetrs_aa_2stage_dcomplex_parameters definition */
zhetrs_aa_2stage_dcomplex_parameters:: zhetrs_aa_2stage_dcomplex_parameters (
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
   printf(" \n zhetrs_aa_2stage lapack_complex_double:  n: %d, uplo: %c  ltb: %d lda: %d \n",
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
       EXPECT_FALSE( true) << "zhetrs_aa_2stage_double_parameters object: malloc error.";
       hetrs_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

zhetrs_aa_2stage_dcomplex_parameters:: ~zhetrs_aa_2stage_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" zhetrs_aa_2stage_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetrs_aa_2stage_free();
}


//  Test fixture class definition
class zhetrs_aa_2stage_test  : public  ::testing::Test {
public:
   zhetrs_aa_2stage_dcomplex_parameters  *zhetrs_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete zhetrs_aa_2stage_obj; }
};


void zhetrs_aa_2stage_test::SetUp(){

    /* LAPACKE_zhetrs_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrs_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, lapack_complex_double *a,
			lapack_int lda, lapack_complex_double *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, lapack_complex_double *b, lapack_int ldb );

    Fptr_NL_LAPACKE_zhetrs_aa_2stage ZHETRS_AA_2STAGE;

    typedef int (*Fptr_NL_LAPACKE_zhetrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,lapack_complex_double *a ,lapack_int lda,
								 lapack_complex_double *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_zhetrf_aa_2stage ZHETRF_AA_2STAGE;

    zhetrs_aa_2stage_obj = new  zhetrs_aa_2stage_dcomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    zhetrs_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetrs_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetrs_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetrs_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHETRS_AA_2STAGE = (Fptr_NL_LAPACKE_zhetrs_aa_2stage)dlsym(
							zhetrs_aa_2stage_obj->hModule, "LAPACKE_zhetrs_aa_2stage");
    ASSERT_TRUE(ZHETRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_zhetrs_aa_2stage symbol";

    ZHETRF_AA_2STAGE = (Fptr_NL_LAPACKE_zhetrf_aa_2stage)dlsym(
							zhetrs_aa_2stage_obj->hModule, "LAPACKE_zhetrf_aa_2stage");
    ASSERT_TRUE(ZHETRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_zhetrf_aa_2stage symbol";

    /* Compute the Netlib-Lapacke's reference o/p */

    zhetrs_aa_2stage_obj->inforef = ZHETRF_AA_2STAGE(
										zhetrs_aa_2stage_obj->matrix_layout,
										zhetrs_aa_2stage_obj->uplo,
										zhetrs_aa_2stage_obj->n,
										zhetrs_aa_2stage_obj->aref,
										zhetrs_aa_2stage_obj->lda,
										zhetrs_aa_2stage_obj->tbref,
										zhetrs_aa_2stage_obj->ltb,
										zhetrs_aa_2stage_obj->ipivref,
										zhetrs_aa_2stage_obj->ipiv2ref
										);

    zhetrs_aa_2stage_obj->inforef = ZHETRS_AA_2STAGE(
										zhetrs_aa_2stage_obj->matrix_layout,
										zhetrs_aa_2stage_obj->uplo,
										zhetrs_aa_2stage_obj->n,
										zhetrs_aa_2stage_obj->nrhs,
										zhetrs_aa_2stage_obj->aref,
										zhetrs_aa_2stage_obj->lda,
										zhetrs_aa_2stage_obj->tbref,
										zhetrs_aa_2stage_obj->ltb,
										zhetrs_aa_2stage_obj->ipivref,
										zhetrs_aa_2stage_obj->ipiv2ref,
										zhetrs_aa_2stage_obj->bref,
										zhetrs_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
	
	    zhetrs_aa_2stage_obj->info = LAPACKE_zhetrf_aa_2stage(
										zhetrs_aa_2stage_obj->matrix_layout,
										zhetrs_aa_2stage_obj->uplo,
										zhetrs_aa_2stage_obj->n,
										zhetrs_aa_2stage_obj->a,
										zhetrs_aa_2stage_obj->lda,
										zhetrs_aa_2stage_obj->tb,
										zhetrs_aa_2stage_obj->ltb,
										zhetrs_aa_2stage_obj->ipiv,
										zhetrs_aa_2stage_obj->ipiv2
										);

    zhetrs_aa_2stage_obj->info = LAPACKE_zhetrs_aa_2stage(
										zhetrs_aa_2stage_obj->matrix_layout,
										zhetrs_aa_2stage_obj->uplo,
										zhetrs_aa_2stage_obj->n,
										zhetrs_aa_2stage_obj->nrhs,
										zhetrs_aa_2stage_obj->a,
										zhetrs_aa_2stage_obj->lda,
										zhetrs_aa_2stage_obj->tb,
										zhetrs_aa_2stage_obj->ltb,
										zhetrs_aa_2stage_obj->ipiv,
										zhetrs_aa_2stage_obj->ipiv2,
										zhetrs_aa_2stage_obj->b,
										zhetrs_aa_2stage_obj->ldb
										);
								 
								 
    if( zhetrs_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrs_aa_2stage is wrong\n",
                    zhetrs_aa_2stage_obj->info );
    }
    if( zhetrs_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrs_aa_2stage is wrong\n",
        zhetrs_aa_2stage_obj->inforef );
    }
	
	zhetrs_aa_2stage_obj->diff_tb =  computeDiff_z( zhetrs_aa_2stage_obj->ltb,
                           zhetrs_aa_2stage_obj->tb,
						   zhetrs_aa_2stage_obj->tbref );

	zhetrs_aa_2stage_obj->diff_b =  computeDiff_z( zhetrs_aa_2stage_obj->n*zhetrs_aa_2stage_obj->nrhs,
                           zhetrs_aa_2stage_obj->b,
						   zhetrs_aa_2stage_obj->bref );

    zhetrs_aa_2stage_obj->diff_ipiv = computeDiff_i( zhetrs_aa_2stage_obj->n,
							zhetrs_aa_2stage_obj->ipiv,
							zhetrs_aa_2stage_obj->ipivref);

    zhetrs_aa_2stage_obj->diff_ipiv2 = computeDiff_i( zhetrs_aa_2stage_obj->n,
							zhetrs_aa_2stage_obj->ipiv2,
							zhetrs_aa_2stage_obj->ipiv2ref);

}

TEST_F(zhetrs_aa_2stage_test, zhetrs_aa_2stage1) {
    EXPECT_NEAR(0.0, zhetrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zhetrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zhetrs_aa_2stage_test, zhetrs_aa_2stage2) {
    EXPECT_NEAR(0.0, zhetrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zhetrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zhetrs_aa_2stage_test, zhetrs_aa_2stage3) {
    EXPECT_NEAR(0.0, zhetrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zhetrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zhetrs_aa_2stage_test, zhetrs_aa_2stage4) {
    EXPECT_NEAR(0.0, zhetrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, zhetrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    idx = Circular_Increment_Index(idx);
}

/* Begin chetrs_aa_2stage_scomplex_parameters  class definition */
class chetrs_aa_2stage_scomplex_parameters{

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
      chetrs_aa_2stage_scomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i, lapack_int nrhs_i);
      ~chetrs_aa_2stage_scomplex_parameters ();
};  /* end of chetrs_aa_2stage_scomplex_parameters  class definition */


/* Constructor chetrs_aa_2stage_scomplex_parameters definition */
chetrs_aa_2stage_scomplex_parameters:: chetrs_aa_2stage_scomplex_parameters (
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
   printf(" \n chetrs_aa_2stage lapack_complex_float:  n: %d, uplo: %c  ltb: %d lda: %d \n",
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
       EXPECT_FALSE( true) << "chetrs_aa_2stage_double_parameters object: malloc error.";
       hetrs_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, nrhs*n);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

chetrs_aa_2stage_scomplex_parameters:: ~chetrs_aa_2stage_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" chetrs_aa_2stage_scomplex_parameters object: destructor invoked. \n");
#endif
   hetrs_aa_2stage_free();
}


//  Test fixture class definition
class chetrs_aa_2stage_test  : public  ::testing::Test {
public:
   chetrs_aa_2stage_scomplex_parameters  *chetrs_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete chetrs_aa_2stage_obj; }
};


void chetrs_aa_2stage_test::SetUp(){

    /* LAPACKE_chetrs_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrs_aa_2stage) ( int matrix_layout,
			char uplo, lapack_int n, lapack_int nrhs, lapack_complex_float *a,
			lapack_int lda, lapack_complex_float *tb, lapack_int ltb, lapack_int *ipiv,
			lapack_int *ipiv2, lapack_complex_float *b, lapack_int ldb );

    Fptr_NL_LAPACKE_chetrs_aa_2stage CHETRS_AA_2STAGE;

    typedef int (*Fptr_NL_LAPACKE_chetrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,lapack_complex_float *a ,lapack_int lda,
								 lapack_complex_float *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_chetrf_aa_2stage CHETRF_AA_2STAGE;

    chetrs_aa_2stage_obj = new  chetrs_aa_2stage_scomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n,
										 lin_solver_paramslist[idx].nrhs);

    chetrs_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetrs_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetrs_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetrs_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHETRS_AA_2STAGE = (Fptr_NL_LAPACKE_chetrs_aa_2stage)dlsym(
							chetrs_aa_2stage_obj->hModule, "LAPACKE_chetrs_aa_2stage");
    ASSERT_TRUE(CHETRS_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_chetrs_aa_2stage symbol";

    CHETRF_AA_2STAGE = (Fptr_NL_LAPACKE_chetrf_aa_2stage)dlsym(
							chetrs_aa_2stage_obj->hModule, "LAPACKE_chetrf_aa_2stage");
    ASSERT_TRUE(CHETRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_chetrf_aa_2stage symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    chetrs_aa_2stage_obj->inforef = CHETRF_AA_2STAGE(
										chetrs_aa_2stage_obj->matrix_layout,
										chetrs_aa_2stage_obj->uplo,
										chetrs_aa_2stage_obj->n,
										chetrs_aa_2stage_obj->aref,
										chetrs_aa_2stage_obj->lda,
										chetrs_aa_2stage_obj->tbref,
										chetrs_aa_2stage_obj->ltb,
										chetrs_aa_2stage_obj->ipivref,
										chetrs_aa_2stage_obj->ipiv2ref
										);

    chetrs_aa_2stage_obj->inforef = CHETRS_AA_2STAGE(
										chetrs_aa_2stage_obj->matrix_layout,
										chetrs_aa_2stage_obj->uplo,
										chetrs_aa_2stage_obj->n,
										chetrs_aa_2stage_obj->nrhs,
										chetrs_aa_2stage_obj->aref,
										chetrs_aa_2stage_obj->lda,
										chetrs_aa_2stage_obj->tbref,
										chetrs_aa_2stage_obj->ltb,
										chetrs_aa_2stage_obj->ipivref,
										chetrs_aa_2stage_obj->ipiv2ref,
										chetrs_aa_2stage_obj->bref,
										chetrs_aa_2stage_obj->ldb
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    chetrs_aa_2stage_obj->info = LAPACKE_chetrf_aa_2stage(
										chetrs_aa_2stage_obj->matrix_layout,
										chetrs_aa_2stage_obj->uplo,
										chetrs_aa_2stage_obj->n,
										chetrs_aa_2stage_obj->a,
										chetrs_aa_2stage_obj->lda,
										chetrs_aa_2stage_obj->tb,
										chetrs_aa_2stage_obj->ltb,
										chetrs_aa_2stage_obj->ipiv,
										chetrs_aa_2stage_obj->ipiv2
										);

    chetrs_aa_2stage_obj->info = LAPACKE_chetrs_aa_2stage(
										chetrs_aa_2stage_obj->matrix_layout,
										chetrs_aa_2stage_obj->uplo,
										chetrs_aa_2stage_obj->n,
										chetrs_aa_2stage_obj->nrhs,
										chetrs_aa_2stage_obj->a,
										chetrs_aa_2stage_obj->lda,
										chetrs_aa_2stage_obj->tb,
										chetrs_aa_2stage_obj->ltb,
										chetrs_aa_2stage_obj->ipiv,
										chetrs_aa_2stage_obj->ipiv2,
										chetrs_aa_2stage_obj->b,
										chetrs_aa_2stage_obj->ldb
										);
								 
								 
    if( chetrs_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrs_aa_2stage is wrong\n",
                    chetrs_aa_2stage_obj->info );
    }
    if( chetrs_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrs_aa_2stage is wrong\n",
        chetrs_aa_2stage_obj->inforef );
    }
	
	chetrs_aa_2stage_obj->diff_tb =  computeDiff_c( chetrs_aa_2stage_obj->ltb,
                           chetrs_aa_2stage_obj->tb,
						   chetrs_aa_2stage_obj->tbref );

	chetrs_aa_2stage_obj->diff_b =  computeDiff_c( chetrs_aa_2stage_obj->n*chetrs_aa_2stage_obj->nrhs,
                           chetrs_aa_2stage_obj->b,
						   chetrs_aa_2stage_obj->bref );

    chetrs_aa_2stage_obj->diff_ipiv = computeDiff_i( chetrs_aa_2stage_obj->n,
							chetrs_aa_2stage_obj->ipiv,
							chetrs_aa_2stage_obj->ipivref);

    chetrs_aa_2stage_obj->diff_ipiv2 = computeDiff_i( chetrs_aa_2stage_obj->n,
							chetrs_aa_2stage_obj->ipiv2,
							chetrs_aa_2stage_obj->ipiv2ref);

}

TEST_F(chetrs_aa_2stage_test, chetrs_aa_2stage1) {
    EXPECT_NEAR(0.0, chetrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, chetrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chetrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chetrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(chetrs_aa_2stage_test, chetrs_aa_2stage2) {
    EXPECT_NEAR(0.0, chetrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, chetrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chetrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chetrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(chetrs_aa_2stage_test, chetrs_aa_2stage3) {
    EXPECT_NEAR(0.0, chetrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, chetrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chetrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chetrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(chetrs_aa_2stage_test, chetrs_aa_2stage4) {
    EXPECT_NEAR(0.0, chetrs_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_NEAR(0.0, chetrs_aa_2stage_obj->diff_b, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chetrs_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chetrs_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}
