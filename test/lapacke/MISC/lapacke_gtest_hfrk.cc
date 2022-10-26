#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define hfrk_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (c != NULL) free(c); \
if (cref != NULL) free(cref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hfrk_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_c;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int ka;
	lapack_complex_float* A;
	lapack_int lda;
	char uplo;
	char trans;
	char transr;
	float alpha, beta;
	/*Output Parameter*/
	lapack_complex_float *c;
	lapack_complex_float *Aref, *cref;
	/*Return Values*/
	lapack_int info, inforef;
   public:
      hfrk_scomplex_parameters (int matrix_layout,  char trans, char uplo, lapack_int m, lapack_int n, float alpha, float beta);
      ~hfrk_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hfrk_scomplex_parameters:: hfrk_scomplex_parameters (int matrix_layout_i, char trans_i, char uplo_i, lapack_int m_i, lapack_int n_i, float alpha_i, float beta_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	trans = trans_i;
	transr = trans;
	alpha = alpha_i;
	beta = beta_i;
	m = m_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hfrk scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	if (trans == 'T')
	{	trans = 'C';
		transr = trans;
	}	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
			if (trans == 'N')
			{
				k = m;
				lda = n;
				ka = k;
			}else {
				k = n;
				lda = k;
				ka = n;
			}
			
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
			if (trans == 'N')
			{
				k = m;
				lda = k;
				ka = n;
			}else { 
				 k = n;
				 lda = n;
				 ka = k;
			}
		
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_c = (n*(n +1)/2);
	bufsize_a = lda*ka;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_c);
	if ((A==NULL) || (Aref==NULL) || \
		(c == NULL) || (cref == NULL)){
		EXPECT_FALSE( true) << "hfrk_float_parameters object: malloc error.";
		hfrk_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hfrk_scomplex_parameters :: ~hfrk_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hfrk_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hfrk_free();

}
/*  Test fixture class definition */
class chfrk_test  : public  ::testing::Test {
public:
   hfrk_scomplex_parameters  *chfrk_obj;
   void SetUp();
   void TearDown () { delete chfrk_obj;}
};

void chfrk_test::SetUp(){

	 /* LAPACKE chfrk prototype */
    typedef int (*Fptr_NL_LAPACKE_chfrk) (int matrix_layout, char transr, char uplo, char trans, lapack_int n,  
											lapack_int k, float alpha, const lapack_complex_float *A, lapack_int lda,
											float beta, lapack_complex_float* c);

    Fptr_NL_LAPACKE_chfrk chfrk;

    chfrk_obj = new hfrk_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].m,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);
						   

    idx = Circular_Increment_Index(idx);

    chfrk_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chfrk_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chfrk_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chfrk_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*chfrk library call */
    chfrk = (Fptr_NL_LAPACKE_chfrk)dlsym(chfrk_obj->hModule, "LAPACKE_chfrk");
    ASSERT_TRUE(chfrk != NULL) << "failed to get the Netlib LAPACKE_chfrk symbol";
    
/*Compute chfrk's  o/p */
    chfrk_obj->inforef = chfrk( chfrk_obj->matrix_layout, chfrk_obj->transr, chfrk_obj->uplo, chfrk_obj->trans, chfrk_obj->n,
								chfrk_obj->k, chfrk_obj->alpha, (const lapack_complex_float *)chfrk_obj->Aref, chfrk_obj->lda, 
								chfrk_obj->beta, chfrk_obj->cref);

    /* Compute libflame's Lapacke o/p  */
    chfrk_obj->info = LAPACKE_chfrk( chfrk_obj->matrix_layout, chfrk_obj->transr, chfrk_obj->uplo, chfrk_obj->trans, chfrk_obj->n, 
										chfrk_obj->k,chfrk_obj->alpha, (const lapack_complex_float *)chfrk_obj->A,  chfrk_obj->lda, chfrk_obj->beta,
										chfrk_obj->c);

    if( chfrk_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chfrk is wrong\n", chfrk_obj->info );
    }
    if( chfrk_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chfrk is wrong\n", 
        chfrk_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chfrk_obj->diff =  computeDiff_c( chfrk_obj->bufsize_c, 
                chfrk_obj->c, chfrk_obj->cref );

}

TEST_F(chfrk_test, chfrk1) {
    EXPECT_NEAR(0.0, chfrk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chfrk_test, chfrk2) {
    EXPECT_NEAR(0.0, chfrk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chfrk_test, chfrk3) {
    EXPECT_NEAR(0.0, chfrk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chfrk_test, chfrk4) {
    EXPECT_NEAR(0.0, chfrk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class hfrk_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_c;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_int ka;
	lapack_complex_double* A;
	lapack_int lda;
	char uplo;
	char trans;
	char transr;
	double alpha, beta;
	/*Output Parameter*/
	lapack_complex_double *c;
	lapack_complex_double *Aref, *cref;
	/*Return Values*/
	lapack_int info, inforef;
   public:
      hfrk_dcomplex_parameters (int matrix_layout,  char trans, char uplo, lapack_int m, lapack_int n, double alpha, double beta);
      ~hfrk_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
hfrk_dcomplex_parameters:: hfrk_dcomplex_parameters (int matrix_layout_i, char trans_i, char uplo_i, lapack_int m_i, lapack_int n_i, double alpha_i, double beta_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	trans = trans_i;
	transr = trans;
	alpha = alpha_i;
	beta = beta_i;
	m = m_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hfrk dcomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	if (trans == 'T')
	{	trans = 'C';
		transr = trans;
	}	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
			if (trans == 'N')
			{
				k = m;
				lda = n;
				ka = k;
			}else {
				k = n;
				lda = k;
				ka = n;
			}
			
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
			if (trans == 'N')
			{
				k = m;
				lda = k;
				ka = n;
			}else { 
				 k = n;
				 lda = n;
				 ka = k;
			}
		
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_c = (n*(n +1)/2);
	bufsize_a = lda*ka;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	if ((A==NULL) || (Aref==NULL) || \
		(c == NULL) || (cref == NULL)){
		EXPECT_FALSE( true) << "hfrk_double_parameters object: malloc error.";
		hfrk_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
hfrk_dcomplex_parameters :: ~hfrk_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hfrk_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hfrk_free();

}
/*  Test fixture class definition */
class zhfrk_test  : public  ::testing::Test {
public:
   hfrk_dcomplex_parameters  *zhfrk_obj;
   void SetUp();
   void TearDown () { delete zhfrk_obj;}
};

void zhfrk_test::SetUp(){

	 /* LAPACKE zhfrk prototype */
    typedef int (*Fptr_NL_LAPACKE_zhfrk) (int matrix_layout, char transr, char uplo, char trans, lapack_int n,  
											lapack_int k, double alpha, const lapack_complex_double *A, lapack_int lda,
											double beta, lapack_complex_double* c);

    Fptr_NL_LAPACKE_zhfrk zhfrk;

    zhfrk_obj = new hfrk_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].transr,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].m,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);
						   

    idx = Circular_Increment_Index(idx);

    zhfrk_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhfrk_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhfrk_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhfrk_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zhfrk library call */
    zhfrk = (Fptr_NL_LAPACKE_zhfrk)dlsym(zhfrk_obj->hModule, "LAPACKE_zhfrk");
    ASSERT_TRUE(zhfrk != NULL) << "failed to get the Netlib LAPACKE_zhfrk symbol";
    
/*Compute zhfrk's  o/p */
    zhfrk_obj->inforef = zhfrk( zhfrk_obj->matrix_layout, zhfrk_obj->transr, zhfrk_obj->uplo, zhfrk_obj->trans, zhfrk_obj->n,
								zhfrk_obj->k, zhfrk_obj->alpha, (const lapack_complex_double *)zhfrk_obj->Aref, zhfrk_obj->lda, 
								zhfrk_obj->beta, zhfrk_obj->cref);

    /* Compute libflame's Lapacke o/p  */
    zhfrk_obj->info = LAPACKE_zhfrk( zhfrk_obj->matrix_layout, zhfrk_obj->transr, zhfrk_obj->uplo, zhfrk_obj->trans, zhfrk_obj->n, 
										zhfrk_obj->k,zhfrk_obj->alpha, (const lapack_complex_double *)zhfrk_obj->A,  zhfrk_obj->lda, zhfrk_obj->beta,
										zhfrk_obj->c);

    if( zhfrk_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhfrk is wrong\n", zhfrk_obj->info );
    }
    if( zhfrk_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhfrk is wrong\n", 
        zhfrk_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhfrk_obj->diff =  computeDiff_z( zhfrk_obj->bufsize_c, 
                zhfrk_obj->c, zhfrk_obj->cref );

}

TEST_F(zhfrk_test, zhfrk1) {
    EXPECT_NEAR(0.0, zhfrk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhfrk_test, zhfrk2) {
    EXPECT_NEAR(0.0, zhfrk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhfrk_test, zhfrk3) {
    EXPECT_NEAR(0.0, zhfrk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhfrk_test, zhfrk4) {
    EXPECT_NEAR(0.0, zhfrk_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
