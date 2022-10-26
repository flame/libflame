#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define ungql_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class ungql_scomplex_parameters{

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
      ungql_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int k, lapack_int lda);
      ~ungql_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
ungql_scomplex_parameters:: ungql_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	k = k_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ungql scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "ungql_float_parameters object: malloc error.";
		ungql_free();
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
ungql_scomplex_parameters :: ~ungql_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ungql_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ungql_free();

}
/*  Test fixture class definition */
class cungql_test  : public  ::testing::Test {
public:
   ungql_scomplex_parameters  *cungql_obj;
   void SetUp();
   void TearDown () { delete cungql_obj; }
};

void cungql_test::SetUp(){

    /* LAPACKE cungql prototype */
    typedef int (*Fptr_NL_LAPACKE_cungql) (int matrix_layout, lapack_int m,lapack_int n, lapack_int k,
											lapack_complex_float *A, lapack_int lda, const lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cungql cungql;

    cungql_obj = new ungql_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].k,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cungql_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cungql_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cungql_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cungql_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cungql = (Fptr_NL_LAPACKE_cungql)dlsym(cungql_obj->hModule, "LAPACKE_cungql");
    ASSERT_TRUE(cungql != NULL) << "failed to get the Netlib LAPACKE_cungql symbol";
    

    cungql_obj->inforef = cungql( cungql_obj->matrix_layout, cungql_obj->m,
								cungql_obj->n, cungql_obj->k, cungql_obj->Aref,
								cungql_obj->lda, cungql_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cungql_obj->info = LAPACKE_cungql( cungql_obj->matrix_layout, cungql_obj->m,
										cungql_obj->n,cungql_obj->k, cungql_obj->A, 
										cungql_obj->lda, cungql_obj->tau);

    if( cungql_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cungql is wrong\n", cungql_obj->info );
    }
    if( cungql_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cungql is wrong\n", 
        cungql_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cungql_obj->diff =  computeDiff_c( cungql_obj->bufsize, 
                cungql_obj->A, cungql_obj->Aref );

}

TEST_F(cungql_test, cungql1) {
    EXPECT_NEAR(0.0, cungql_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungql_test, cungql2) {
    EXPECT_NEAR(0.0, cungql_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungql_test, cungql3) {
    EXPECT_NEAR(0.0, cungql_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungql_test, cungql4) {
    EXPECT_NEAR(0.0, cungql_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class ungql_dcomplex_parameters{

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
      ungql_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int k, lapack_int lda);
      ~ungql_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
ungql_dcomplex_parameters:: ungql_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int k_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	k = k_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ungql dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "ungql_float_parameters object: malloc error.";
		ungql_free();
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
ungql_dcomplex_parameters :: ~ungql_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ungql_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ungql_free();

}
/*  Test fixture class definition */
class zungql_test  : public  ::testing::Test {
public:
   ungql_dcomplex_parameters  *zungql_obj;
   void SetUp();
   void TearDown () { delete zungql_obj; }
};

void zungql_test::SetUp(){

    /* LAPACKE zungql prototype */
    typedef int (*Fptr_NL_LAPACKE_zungql) (int matrix_layout, lapack_int m,lapack_int n, lapack_int k,
											lapack_complex_double *A, lapack_int lda, const lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zungql zungql;

    zungql_obj = new ungql_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].k,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zungql_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zungql_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zungql_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zungql_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zungql = (Fptr_NL_LAPACKE_zungql)dlsym(zungql_obj->hModule, "LAPACKE_zungql");
    ASSERT_TRUE(zungql != NULL) << "failed to get the Netlib LAPACKE_zungql symbol";
    

    zungql_obj->inforef = zungql( zungql_obj->matrix_layout, zungql_obj->m,
								zungql_obj->n, zungql_obj->k, zungql_obj->Aref,
								zungql_obj->lda, zungql_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zungql_obj->info = LAPACKE_zungql( zungql_obj->matrix_layout, zungql_obj->m,
										zungql_obj->n,zungql_obj->k, zungql_obj->A, 
										zungql_obj->lda, zungql_obj->tau);

    if( zungql_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zungql is wrong\n", zungql_obj->info );
    }
    if( zungql_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zungql is wrong\n", 
        zungql_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zungql_obj->diff =  computeDiff_z( zungql_obj->bufsize, 
                zungql_obj->A, zungql_obj->Aref );

}

TEST_F(zungql_test, zungql1) {
    EXPECT_NEAR(0.0, zungql_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungql_test, zungql2) {
    EXPECT_NEAR(0.0, zungql_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungql_test, zungql3) {
    EXPECT_NEAR(0.0, zungql_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungql_test, zungql4) {
    EXPECT_NEAR(0.0, zungql_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


