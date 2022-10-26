#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define ungrq_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class ungrq_scomplex_parameters{

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
      ungrq_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int k, lapack_int lda);
      ~ungrq_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
ungrq_scomplex_parameters:: ungrq_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	k = k_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ungrq scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "ungrq_float_parameters object: malloc error.";
		ungrq_free();
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
ungrq_scomplex_parameters :: ~ungrq_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ungrq_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ungrq_free();

}
/*  Test fixture class definition */
class cungrq_test  : public  ::testing::Test {
public:
   ungrq_scomplex_parameters  *cungrq_obj;
   void SetUp();
   void TearDown () { delete cungrq_obj; }
};

void cungrq_test::SetUp(){

    /* LAPACKE cungrq prototype */
    typedef int (*Fptr_NL_LAPACKE_cungrq) (int matrix_layout, lapack_int m,lapack_int n, lapack_int k,
											lapack_complex_float *A, lapack_int lda, const lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cungrq cungrq;

    cungrq_obj = new ungrq_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
			   eig_paramslist[idx].k,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cungrq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cungrq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cungrq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cungrq_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cungrq = (Fptr_NL_LAPACKE_cungrq)dlsym(cungrq_obj->hModule, "LAPACKE_cungrq");
    ASSERT_TRUE(cungrq != NULL) << "failed to get the Netlib LAPACKE_cungrq symbol";
    

    cungrq_obj->inforef = cungrq( cungrq_obj->matrix_layout, cungrq_obj->m,
								cungrq_obj->n, cungrq_obj->k, cungrq_obj->Aref,
								cungrq_obj->lda, cungrq_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cungrq_obj->info = LAPACKE_cungrq( cungrq_obj->matrix_layout, cungrq_obj->m,
										cungrq_obj->n,cungrq_obj->k, cungrq_obj->A, 
										cungrq_obj->lda, cungrq_obj->tau);

    if( cungrq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cungrq is wrong\n", cungrq_obj->info );
    }
    if( cungrq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cungrq is wrong\n", 
        cungrq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cungrq_obj->diff =  computeDiff_c( cungrq_obj->bufsize, 
                cungrq_obj->A, cungrq_obj->Aref );

}

TEST_F(cungrq_test, cungrq1) {
    EXPECT_NEAR(0.0, cungrq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungrq_test, cungrq2) {
    EXPECT_NEAR(0.0, cungrq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungrq_test, cungrq3) {
    EXPECT_NEAR(0.0, cungrq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungrq_test, cungrq4) {
    EXPECT_NEAR(0.0, cungrq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class ungrq_dcomplex_parameters{

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
      ungrq_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int k, lapack_int lda);
      ~ungrq_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
ungrq_dcomplex_parameters:: ungrq_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	k = k_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ungrq dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "ungrq_float_parameters object: malloc error.";
		ungrq_free();
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
ungrq_dcomplex_parameters :: ~ungrq_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ungrq_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ungrq_free();

}
/*  Test fixture class definition */
class zungrq_test  : public  ::testing::Test {
public:
   ungrq_dcomplex_parameters  *zungrq_obj;
   void SetUp();
   void TearDown () { delete zungrq_obj; }
};

void zungrq_test::SetUp(){

    /* LAPACKE zungrq prototype */
    typedef int (*Fptr_NL_LAPACKE_zungrq) (int matrix_layout, lapack_int m,lapack_int n, lapack_int k,
											lapack_complex_double *A, lapack_int lda, const lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zungrq zungrq;

    zungrq_obj = new ungrq_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
			   eig_paramslist[idx].k,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zungrq_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zungrq_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zungrq_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zungrq_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zungrq = (Fptr_NL_LAPACKE_zungrq)dlsym(zungrq_obj->hModule, "LAPACKE_zungrq");
    ASSERT_TRUE(zungrq != NULL) << "failed to get the Netlib LAPACKE_zungrq symbol";
    

    zungrq_obj->inforef = zungrq( zungrq_obj->matrix_layout, zungrq_obj->m,
								zungrq_obj->n, zungrq_obj->k, zungrq_obj->Aref,
								zungrq_obj->lda, zungrq_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zungrq_obj->info = LAPACKE_zungrq( zungrq_obj->matrix_layout, zungrq_obj->m,
										zungrq_obj->n,zungrq_obj->k, zungrq_obj->A, 
										zungrq_obj->lda, zungrq_obj->tau);

    if( zungrq_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zungrq is wrong\n", zungrq_obj->info );
    }
    if( zungrq_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zungrq is wrong\n", 
        zungrq_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zungrq_obj->diff =  computeDiff_z( zungrq_obj->bufsize, 
                zungrq_obj->A, zungrq_obj->Aref );

}

TEST_F(zungrq_test, zungrq1) {
    EXPECT_NEAR(0.0, zungrq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungrq_test, zungrq2) {
    EXPECT_NEAR(0.0, zungrq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungrq_test, zungrq3) {
    EXPECT_NEAR(0.0, zungrq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungrq_test, zungrq4) {
    EXPECT_NEAR(0.0, zungrq_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


