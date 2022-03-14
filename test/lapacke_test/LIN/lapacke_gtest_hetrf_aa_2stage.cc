#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define hetrf_aa_2stage_free() \
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

/* Begin chetrf_aa_2stage_scomplex_parameters  class definition */
class chetrf_aa_2stage_scomplex_parameters{

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
      chetrf_aa_2stage_scomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i);
      ~chetrf_aa_2stage_scomplex_parameters ();
};  /* end of chetrf_aa_2stage_scomplex_parameters  class definition */


/* Constructor chetrf_aa_2stage_scomplex_parameters definition */
chetrf_aa_2stage_scomplex_parameters:: chetrf_aa_2stage_scomplex_parameters (
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
   printf(" \n chetrf_aa_2stage lapack_complex_float:  n: %d, uplo: %c  ltb: %d lda: %d \n",
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
       EXPECT_FALSE( true) << "chetrf_aa_2stage_double_parameters object: malloc error.";
       hetrf_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_scomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

chetrf_aa_2stage_scomplex_parameters:: ~chetrf_aa_2stage_scomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" chetrf_aa_2stage_scomplex_parameters object: destructor invoked. \n");
#endif
   hetrf_aa_2stage_free();
}


//  Test fixture class definition
class chetrf_aa_2stage_test  : public  ::testing::Test {
public:
   chetrf_aa_2stage_scomplex_parameters  *chetrf_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete chetrf_aa_2stage_obj; }
};


void chetrf_aa_2stage_test::SetUp(){

    /* LAPACKE_chetrf_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,lapack_complex_float *a ,lapack_int lda,
								 lapack_complex_float *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_chetrf_aa_2stage CHETRF_AA_2STAGE;

    chetrf_aa_2stage_obj = new  chetrf_aa_2stage_scomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n);

    chetrf_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetrf_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetrf_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetrf_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    CHETRF_AA_2STAGE = (Fptr_NL_LAPACKE_chetrf_aa_2stage)dlsym(
							chetrf_aa_2stage_obj->hModule, "LAPACKE_chetrf_aa_2stage");
    ASSERT_TRUE(CHETRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_chetrf_aa_2stage symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    chetrf_aa_2stage_obj->inforef = CHETRF_AA_2STAGE(
										chetrf_aa_2stage_obj->matrix_layout,
										chetrf_aa_2stage_obj->uplo,
										chetrf_aa_2stage_obj->n,
										chetrf_aa_2stage_obj->aref,
										chetrf_aa_2stage_obj->lda,
										chetrf_aa_2stage_obj->tbref,
										chetrf_aa_2stage_obj->ltb,
										chetrf_aa_2stage_obj->ipivref,
										chetrf_aa_2stage_obj->ipiv2ref
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    chetrf_aa_2stage_obj->info = LAPACKE_chetrf_aa_2stage(
										chetrf_aa_2stage_obj->matrix_layout,
										chetrf_aa_2stage_obj->uplo,
										chetrf_aa_2stage_obj->n,
										chetrf_aa_2stage_obj->a,
										chetrf_aa_2stage_obj->lda,
										chetrf_aa_2stage_obj->tb,
										chetrf_aa_2stage_obj->ltb,
										chetrf_aa_2stage_obj->ipiv,
										chetrf_aa_2stage_obj->ipiv2
										);
								 
								 
    if( chetrf_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_chetrf_aa_2stage is wrong\n",
                    chetrf_aa_2stage_obj->info );
    }
    if( chetrf_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrf_aa_2stage is wrong\n",
        chetrf_aa_2stage_obj->inforef );
    }
	
	chetrf_aa_2stage_obj->diff_tb =  computeDiff_c( chetrf_aa_2stage_obj->ltb,
                           chetrf_aa_2stage_obj->tb,
						   chetrf_aa_2stage_obj->tbref );

    chetrf_aa_2stage_obj->diff_ipiv = computeDiff_i( chetrf_aa_2stage_obj->n,
							chetrf_aa_2stage_obj->ipiv,
							chetrf_aa_2stage_obj->ipivref);

    chetrf_aa_2stage_obj->diff_ipiv2 = computeDiff_i( chetrf_aa_2stage_obj->n,
							chetrf_aa_2stage_obj->ipiv2,
							chetrf_aa_2stage_obj->ipiv2ref);

}

TEST_F(chetrf_aa_2stage_test, chetrf_aa_2stage1) {
    EXPECT_NEAR(0.0, chetrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chetrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chetrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(chetrf_aa_2stage_test, chetrf_aa_2stage2) {
    EXPECT_NEAR(0.0, chetrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chetrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chetrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(chetrf_aa_2stage_test, chetrf_aa_2stage3) {
    EXPECT_NEAR(0.0, chetrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chetrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chetrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(chetrf_aa_2stage_test, chetrf_aa_2stage4) {
    EXPECT_NEAR(0.0, chetrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, chetrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, chetrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

/* Begin zhetrf_aa_2stage_dcomplex_parameters  class definition */
class zhetrf_aa_2stage_dcomplex_parameters{

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
      zhetrf_aa_2stage_dcomplex_parameters (int matrix_layout_i, char uplo_i,
                               lapack_int n_i);
      ~zhetrf_aa_2stage_dcomplex_parameters ();
};  /* end of zhetrf_aa_2stage_dcomplex_parameters  class definition */


/* Constructor zhetrf_aa_2stage_dcomplex_parameters definition */
zhetrf_aa_2stage_dcomplex_parameters:: zhetrf_aa_2stage_dcomplex_parameters (
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
   printf(" \n zhetrf_aa_2stage lapack_complex_double:  n: %d, uplo: %c  ltb: %d lda: %d \n",
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
       EXPECT_FALSE( true) << "zhetrf_aa_2stage_double_parameters object: malloc error.";
       hetrf_aa_2stage_free();
       exit(-1);
    }

    /* Initialization of input Buffers */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( a, aref, lda, n, 'S');
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv, ipivref, n, 0);
    lapacke_gtest_init_int_buffer_pair_with_constant(ipiv2, ipiv2ref, n, 0);
    lapacke_gtest_init_dcomplex_buffer_pair_with_constant(tb, tbref, ltb*ltb, 0.0);

   } /* end of Constructor  */

zhetrf_aa_2stage_dcomplex_parameters:: ~zhetrf_aa_2stage_dcomplex_parameters ()
{
#if LAPACKE_TEST_VERBOSE
   printf(" zhetrf_aa_2stage_dcomplex_parameters object: destructor invoked. \n");
#endif
   hetrf_aa_2stage_free();
}


//  Test fixture class definition
class zhetrf_aa_2stage_test  : public  ::testing::Test {
public:
   zhetrf_aa_2stage_dcomplex_parameters  *zhetrf_aa_2stage_obj;
   void SetUp();  
   void TearDown () { delete zhetrf_aa_2stage_obj; }
};


void zhetrf_aa_2stage_test::SetUp(){

    /* LAPACKE_zhetrf_aa_2stage prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrf_aa_2stage) ( int matrix_layout ,char uplo,
                                 lapack_int n ,lapack_complex_double *a ,lapack_int lda,
								 lapack_complex_double *tb ,lapack_int ltb,
                                 lapack_int *ipiv, lapack_int *ipiv2 );

    Fptr_NL_LAPACKE_zhetrf_aa_2stage ZHETRF_AA_2STAGE;

    zhetrf_aa_2stage_obj = new  zhetrf_aa_2stage_dcomplex_parameters(
								lin_solver_paramslist[idx].matrix_layout,
                                         lin_solver_paramslist[idx].Uplo,
                                         lin_solver_paramslist[idx].n);

    zhetrf_aa_2stage_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetrf_aa_2stage_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetrf_aa_2stage_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetrf_aa_2stage_obj->hModule != NULL) << "Netlib lapacke handle NULL";

    ZHETRF_AA_2STAGE = (Fptr_NL_LAPACKE_zhetrf_aa_2stage)dlsym(
							zhetrf_aa_2stage_obj->hModule, "LAPACKE_zhetrf_aa_2stage");
    ASSERT_TRUE(ZHETRF_AA_2STAGE != NULL) << "failed to get the Netlib LAPACKE_zhetrf_aa_2stage symbol";

    /* Compute the Netlib-Lapacke's reference o/p */
    zhetrf_aa_2stage_obj->inforef = ZHETRF_AA_2STAGE(
										zhetrf_aa_2stage_obj->matrix_layout,
										zhetrf_aa_2stage_obj->uplo,
										zhetrf_aa_2stage_obj->n,
										zhetrf_aa_2stage_obj->aref,
										zhetrf_aa_2stage_obj->lda,
										zhetrf_aa_2stage_obj->tbref,
										zhetrf_aa_2stage_obj->ltb,
										zhetrf_aa_2stage_obj->ipivref,
										zhetrf_aa_2stage_obj->ipiv2ref
										);

    /* Compute the Libflme lapacke o/p by invoking Netlib-Lapack's API */
    zhetrf_aa_2stage_obj->info = LAPACKE_zhetrf_aa_2stage(
										zhetrf_aa_2stage_obj->matrix_layout,
										zhetrf_aa_2stage_obj->uplo,
										zhetrf_aa_2stage_obj->n,
										zhetrf_aa_2stage_obj->a,
										zhetrf_aa_2stage_obj->lda,
										zhetrf_aa_2stage_obj->tb,
										zhetrf_aa_2stage_obj->ltb,
										zhetrf_aa_2stage_obj->ipiv,
										zhetrf_aa_2stage_obj->ipiv2
										);
								 
								 
    if( zhetrf_aa_2stage_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame LAPACKE_zhetrf_aa_2stage is wrong\n",
                    zhetrf_aa_2stage_obj->info );
    }
    if( zhetrf_aa_2stage_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrf_aa_2stage is wrong\n",
        zhetrf_aa_2stage_obj->inforef );
    }
	
	zhetrf_aa_2stage_obj->diff_tb =  computeDiff_z( zhetrf_aa_2stage_obj->ltb,
                           zhetrf_aa_2stage_obj->tb,
						   zhetrf_aa_2stage_obj->tbref );

    zhetrf_aa_2stage_obj->diff_ipiv = computeDiff_i( zhetrf_aa_2stage_obj->n,
							zhetrf_aa_2stage_obj->ipiv,
							zhetrf_aa_2stage_obj->ipivref);

    zhetrf_aa_2stage_obj->diff_ipiv2 = computeDiff_i( zhetrf_aa_2stage_obj->n,
							zhetrf_aa_2stage_obj->ipiv2,
							zhetrf_aa_2stage_obj->ipiv2ref);

}

TEST_F(zhetrf_aa_2stage_test, zhetrf_aa_2stage1) {
    EXPECT_NEAR(0.0, zhetrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zhetrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zhetrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zhetrf_aa_2stage_test, zhetrf_aa_2stage2) {
    EXPECT_NEAR(0.0, zhetrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zhetrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zhetrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zhetrf_aa_2stage_test, zhetrf_aa_2stage3) {
    EXPECT_NEAR(0.0, zhetrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zhetrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zhetrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

TEST_F(zhetrf_aa_2stage_test, zhetrf_aa_2stage4) {
    EXPECT_NEAR(0.0, zhetrf_aa_2stage_obj->diff_tb, lin_solver_paramslist[idx].solver_threhold);
    EXPECT_EQ(0, zhetrf_aa_2stage_obj->diff_ipiv);
    EXPECT_EQ(0, zhetrf_aa_2stage_obj->diff_ipiv2);
    idx = Circular_Increment_Index(idx);
}

