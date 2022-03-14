#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define unglq_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class unglq_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_complex_float* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float* tau;
	lapack_complex_float *Aref, *tauref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      unglq_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int k, lapack_int lda);
      ~unglq_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
unglq_scomplex_parameters:: unglq_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	k = k_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n unglq scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(lapack_complex_float);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(lapack_complex_float);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (k*sizeof(lapack_complex_float)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "unglq_float_parameters object: malloc error.";
		unglq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(k);i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
unglq_scomplex_parameters :: ~unglq_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unglq_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unglq_free();

}
/*  Test fixture class definition */
class cunglq_test  : public  ::testing::Test {
public:
   unglq_scomplex_parameters  *cunglq_obj;
   void SetUp();
   void TearDown () { delete cunglq_obj; }
};

void cunglq_test::SetUp(){

    /* LAPACKE cunglq prototype */
    typedef int (*Fptr_NL_LAPACKE_cunglq) (int matrix_layout, lapack_int m,lapack_int n, lapack_int k,
											lapack_complex_float *A, lapack_int lda, const lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cunglq cunglq;

    cunglq_obj = new unglq_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
			   eig_paramslist[idx].k,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cunglq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cunglq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cunglq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cunglq_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cunglq = (Fptr_NL_LAPACKE_cunglq)dlsym(cunglq_obj->hModule, "LAPACKE_cunglq");
    ASSERT_TRUE(cunglq != NULL) << "failed to get the Netlib LAPACKE_cunglq symbol";
    

    cunglq_obj->inforef = cunglq( cunglq_obj->matrix_layout, cunglq_obj->m,
								cunglq_obj->n, cunglq_obj->k, cunglq_obj->Aref,
								cunglq_obj->lda, cunglq_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cunglq_obj->info = LAPACKE_cunglq( cunglq_obj->matrix_layout, cunglq_obj->m,
										cunglq_obj->n,cunglq_obj->k, cunglq_obj->A, 
										cunglq_obj->lda, cunglq_obj->tau);

    if( cunglq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cunglq is wrong\n", cunglq_obj->info );
    }
    if( cunglq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cunglq is wrong\n", 
        cunglq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cunglq_obj->diff =  computeDiff_c( cunglq_obj->bufsize, 
                cunglq_obj->A, cunglq_obj->Aref );

}

TEST_F(cunglq_test, cunglq1) {
    EXPECT_NEAR(0.0, cunglq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunglq_test, cunglq2) {
    EXPECT_NEAR(0.0, cunglq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunglq_test, cunglq3) {
    EXPECT_NEAR(0.0, cunglq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cunglq_test, cunglq4) {
    EXPECT_NEAR(0.0, cunglq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class unglq_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int k;
	lapack_complex_double* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double* tau;
	lapack_complex_double *Aref, *tauref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      unglq_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int k, lapack_int lda);
      ~unglq_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
unglq_dcomplex_parameters:: unglq_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	k = k_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n unglq dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(lapack_complex_double);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(lapack_complex_double);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (k*sizeof(lapack_complex_double)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "unglq_float_parameters object: malloc error.";
		unglq_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(k);i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
unglq_dcomplex_parameters :: ~unglq_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" unglq_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   unglq_free();

}
/*  Test fixture class definition */
class zunglq_test  : public  ::testing::Test {
public:
   unglq_dcomplex_parameters  *zunglq_obj;
   void SetUp();
   void TearDown () { delete zunglq_obj; }
};

void zunglq_test::SetUp(){

    /* LAPACKE zunglq prototype */
    typedef int (*Fptr_NL_LAPACKE_zunglq) (int matrix_layout, lapack_int m,lapack_int n, lapack_int k,
											lapack_complex_double *A, lapack_int lda, const lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zunglq zunglq;

    zunglq_obj = new unglq_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
			   eig_paramslist[idx].k,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zunglq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zunglq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zunglq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zunglq_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zunglq = (Fptr_NL_LAPACKE_zunglq)dlsym(zunglq_obj->hModule, "LAPACKE_zunglq");
    ASSERT_TRUE(zunglq != NULL) << "failed to get the Netlib LAPACKE_zunglq symbol";
    

    zunglq_obj->inforef = zunglq( zunglq_obj->matrix_layout, zunglq_obj->m,
								zunglq_obj->n, zunglq_obj->k, zunglq_obj->Aref,
								zunglq_obj->lda, zunglq_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zunglq_obj->info = LAPACKE_zunglq( zunglq_obj->matrix_layout, zunglq_obj->m,
										zunglq_obj->n,zunglq_obj->k, zunglq_obj->A, 
										zunglq_obj->lda, zunglq_obj->tau);

    if( zunglq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zunglq is wrong\n", zunglq_obj->info );
    }
    if( zunglq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zunglq is wrong\n", 
        zunglq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zunglq_obj->diff =  computeDiff_z( zunglq_obj->bufsize, 
                zunglq_obj->A, zunglq_obj->Aref );

}

TEST_F(zunglq_test, zunglq1) {
    EXPECT_NEAR(0.0, zunglq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunglq_test, zunglq2) {
    EXPECT_NEAR(0.0, zunglq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunglq_test, zunglq3) {
    EXPECT_NEAR(0.0, zunglq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zunglq_test, zunglq4) {
    EXPECT_NEAR(0.0, zunglq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


