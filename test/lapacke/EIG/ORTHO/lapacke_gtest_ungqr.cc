#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define ungqr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class ungqr_scomplex_parameters{

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
      ungqr_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int k, lapack_int lda);
      ~ungqr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
ungqr_scomplex_parameters:: ungqr_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	k = k_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ungqr scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "ungqr_float_parameters object: malloc error.";
		ungqr_free();
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
ungqr_scomplex_parameters :: ~ungqr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ungqr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ungqr_free();

}
/*  Test fixture class definition */
class cungqr_test  : public  ::testing::Test {
public:
   ungqr_scomplex_parameters  *cungqr_obj;
   void SetUp();
   void TearDown () { delete cungqr_obj; }
};

void cungqr_test::SetUp(){

    /* LAPACKE cungqr prototype */
    typedef int (*Fptr_NL_LAPACKE_cungqr) (int matrix_layout, lapack_int m,lapack_int n, lapack_int k,
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cungqr cungqr;

    cungqr_obj = new ungqr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].k,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cungqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cungqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cungqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cungqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cungqr = (Fptr_NL_LAPACKE_cungqr)dlsym(cungqr_obj->hModule, "LAPACKE_cungqr");
    ASSERT_TRUE(cungqr != NULL) << "failed to get the Netlib LAPACKE_cungqr symbol";
    

    cungqr_obj->inforef = cungqr( cungqr_obj->matrix_layout, cungqr_obj->m,
								cungqr_obj->n, cungqr_obj->k, cungqr_obj->Aref,
								cungqr_obj->lda, cungqr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cungqr_obj->info = LAPACKE_cungqr( cungqr_obj->matrix_layout, cungqr_obj->m,
										cungqr_obj->n,cungqr_obj->k, cungqr_obj->A, 
										cungqr_obj->lda, cungqr_obj->tau);

    if( cungqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cungqr is wrong\n", cungqr_obj->info );
    }
    if( cungqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cungqr is wrong\n", 
        cungqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cungqr_obj->diff =  computeDiff_c( cungqr_obj->bufsize, 
                cungqr_obj->A, cungqr_obj->Aref );

}

TEST_F(cungqr_test, cungqr1) {
    EXPECT_NEAR(0.0, cungqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungqr_test, cungqr2) {
    EXPECT_NEAR(0.0, cungqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungqr_test, cungqr3) {
    EXPECT_NEAR(0.0, cungqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungqr_test, cungqr4) {
    EXPECT_NEAR(0.0, cungqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class ungqr_dcomplex_parameters{

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
      ungqr_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int k, lapack_int lda);
      ~ungqr_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
ungqr_dcomplex_parameters:: ungqr_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	k = k_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ungqr dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "ungqr_float_parameters object: malloc error.";
		ungqr_free();
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
ungqr_dcomplex_parameters :: ~ungqr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ungqr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ungqr_free();

}
/*  Test fixture class definition */
class zungqr_test  : public  ::testing::Test {
public:
   ungqr_dcomplex_parameters  *zungqr_obj;
   void SetUp();
   void TearDown () { delete zungqr_obj; }
};

void zungqr_test::SetUp(){

    /* LAPACKE zungqr prototype */
    typedef int (*Fptr_NL_LAPACKE_zungqr) (int matrix_layout, lapack_int m,lapack_int n, lapack_int k,
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zungqr zungqr;

    zungqr_obj = new ungqr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].k,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zungqr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zungqr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zungqr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zungqr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zungqr = (Fptr_NL_LAPACKE_zungqr)dlsym(zungqr_obj->hModule, "LAPACKE_zungqr");
    ASSERT_TRUE(zungqr != NULL) << "failed to get the Netlib LAPACKE_zungqr symbol";
    

    zungqr_obj->inforef = zungqr( zungqr_obj->matrix_layout, zungqr_obj->m,
								zungqr_obj->n, zungqr_obj->k, zungqr_obj->Aref,
								zungqr_obj->lda, zungqr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zungqr_obj->info = LAPACKE_zungqr( zungqr_obj->matrix_layout, zungqr_obj->m,
										zungqr_obj->n,zungqr_obj->k, zungqr_obj->A, 
										zungqr_obj->lda, zungqr_obj->tau);

    if( zungqr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zungqr is wrong\n", zungqr_obj->info );
    }
    if( zungqr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zungqr is wrong\n", 
        zungqr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zungqr_obj->diff =  computeDiff_z( zungqr_obj->bufsize, 
                zungqr_obj->A, zungqr_obj->Aref );

}

TEST_F(zungqr_test, zungqr1) {
    EXPECT_NEAR(0.0, zungqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungqr_test, zungqr2) {
    EXPECT_NEAR(0.0, zungqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungqr_test, zungqr3) {
    EXPECT_NEAR(0.0, zungqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungqr_test, zungqr4) {
    EXPECT_NEAR(0.0, zungqr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


