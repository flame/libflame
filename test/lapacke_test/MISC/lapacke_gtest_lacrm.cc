#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define lacrm_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (c!=NULL)  free(c);\
if (cref!=NULL) free(cref); \
if (b!=NULL)  free(b);\
if (bref!=NULL) free(bref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class lacrm_scomplex_parameters{

   public:
	int bufsize;
	int bufsize_b;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_complex_float* A;
	float *b;
	lapack_int lda, ldb, ldc;
	/*Output Parameter*/
	lapack_complex_float* c;
	float *bref;
	lapack_complex_float *Aref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacrm_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb, lapack_int ldc);
      ~lacrm_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lacrm_scomplex_parameters:: lacrm_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i, lapack_int ldc_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldb = ldb_i;
	ldc = ldc_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lacrm scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	
	bufsize = lda*n;
	bufsize_b = ldb*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) || \
		(c==NULL) || (cref==NULL) || \
		(b == NULL) || (bref == NULL)){
		EXPECT_FALSE( true) << "lacrm_float_parameters object: malloc error.";
		lacrm_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lacrm_scomplex_parameters :: ~lacrm_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacrm_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacrm_free();

}
/*  Test fixture class definition */
class clacrm_test  : public  ::testing::Test {
public:
   lacrm_scomplex_parameters  *clacrm_obj;
   void SetUp();
   void TearDown () { delete clacrm_obj; }
};

void clacrm_test::SetUp(){

    /* LAPACKE clacrm prototype */
    typedef int (*Fptr_NL_LAPACKE_clacrm) (int matrix_layout, lapack_int m,lapack_int n, const lapack_complex_float *A, 
											lapack_int lda, const float *b, lapack_int ldb, lapack_complex_float* c, lapack_int ldc);

    Fptr_NL_LAPACKE_clacrm clacrm;

    clacrm_obj = new lacrm_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].n,
						eig_paramslist[idx].m,
						eig_paramslist[idx].lda,
						eig_paramslist[idx].lda,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    clacrm_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clacrm_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clacrm_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clacrm_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clacrm = (Fptr_NL_LAPACKE_clacrm)dlsym(clacrm_obj->hModule, "LAPACKE_clacrm");
    ASSERT_TRUE(clacrm != NULL) << "failed to get the Netlib LAPACKE_clacrm symbol";
    

    clacrm_obj->inforef = clacrm( clacrm_obj->matrix_layout, clacrm_obj->m,
								clacrm_obj->n, (const lapack_complex_float*)clacrm_obj->Aref,
								clacrm_obj->lda, (const float *)clacrm_obj->bref, clacrm_obj->ldb, clacrm_obj->cref, clacrm_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    clacrm_obj->info = LAPACKE_clacrm( clacrm_obj->matrix_layout, clacrm_obj->m,
										clacrm_obj->n, (const lapack_complex_float*)clacrm_obj->A, 
										clacrm_obj->lda, (const float *)clacrm_obj->b, clacrm_obj->ldb, clacrm_obj->c, clacrm_obj->ldc);

    if( clacrm_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clacrm is wrong\n", clacrm_obj->info );
    }
    if( clacrm_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clacrm is wrong\n", 
        clacrm_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clacrm_obj->diff =  computeDiff_c( clacrm_obj->bufsize_b, 
                clacrm_obj->c, clacrm_obj->cref );

}

TEST_F(clacrm_test, clacrm1) {
    EXPECT_NEAR(0.0, clacrm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(clacrm_test, clacrm2) {
    EXPECT_NEAR(0.0, clacrm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(clacrm_test, clacrm3) {
    EXPECT_NEAR(0.0, clacrm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(clacrm_test, clacrm4) {
    EXPECT_NEAR(0.0, clacrm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lacrm_dcomplex_parameters{

   public:
	int bufsize;
	int bufsize_b;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_complex_double* A;
	double *b;
	lapack_int lda, ldb, ldc;
	/*Output Parameter*/
	lapack_complex_double* c;
	double *bref;
	lapack_complex_double *Aref, *cref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacrm_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb, lapack_int ldc);
      ~lacrm_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lacrm_dcomplex_parameters:: lacrm_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i, lapack_int ldc_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldb = ldb_i;
	ldc = ldc_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lacrm dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	bufsize = lda*n;
	bufsize_b = ldb*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) || \
		(c==NULL) || (cref==NULL) || \
		(b == NULL) || (bref == NULL)){
		EXPECT_FALSE( true) << "lacrm_double_parameters object: malloc error.";
		lacrm_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lacrm_dcomplex_parameters :: ~lacrm_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacrm_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacrm_free();

}
/*  Test fixture class definition */
class zlacrm_test  : public  ::testing::Test {
public:
   lacrm_dcomplex_parameters  *zlacrm_obj;
   void SetUp();
   void TearDown () { delete zlacrm_obj; }
};

void zlacrm_test::SetUp(){

    /* LAPACKE zlacrm prototype */
    typedef int (*Fptr_NL_LAPACKE_zlacrm) (int matrix_layout, lapack_int m,lapack_int n, const lapack_complex_double *A, 
											lapack_int lda, const double *b, lapack_int ldb, lapack_complex_double* c, lapack_int ldc);

    Fptr_NL_LAPACKE_zlacrm zlacrm;

    zlacrm_obj = new lacrm_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].n,
						eig_paramslist[idx].m,
						eig_paramslist[idx].lda,
						eig_paramslist[idx].lda,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zlacrm_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlacrm_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlacrm_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlacrm_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlacrm = (Fptr_NL_LAPACKE_zlacrm)dlsym(zlacrm_obj->hModule, "LAPACKE_zlacrm");
    ASSERT_TRUE(zlacrm != NULL) << "failed to get the Netlib LAPACKE_zlacrm symbol";
    

    zlacrm_obj->inforef = zlacrm( zlacrm_obj->matrix_layout, zlacrm_obj->m,
								zlacrm_obj->n, (const lapack_complex_double*)zlacrm_obj->Aref,
								zlacrm_obj->lda, (const double *)zlacrm_obj->bref, zlacrm_obj->ldb, zlacrm_obj->cref, zlacrm_obj->ldc);

    /* Compute libflame's Lapacke o/p  */
    zlacrm_obj->info = LAPACKE_zlacrm( zlacrm_obj->matrix_layout, zlacrm_obj->m,
										zlacrm_obj->n, (const lapack_complex_double*)zlacrm_obj->A, 
										zlacrm_obj->lda, (const double *)zlacrm_obj->b, zlacrm_obj->ldb, zlacrm_obj->c, zlacrm_obj->ldc);

    if( zlacrm_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlacrm is wrong\n", zlacrm_obj->info );
    }
    if( zlacrm_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlacrm is wrong\n", 
        zlacrm_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlacrm_obj->diff =  computeDiff_z( zlacrm_obj->bufsize_b, 
                zlacrm_obj->c, zlacrm_obj->cref );

}

TEST_F(zlacrm_test, zlacrm1) {
    EXPECT_NEAR(0.0, zlacrm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zlacrm_test, zlacrm2) {
    EXPECT_NEAR(0.0, zlacrm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zlacrm_test, zlacrm3) {
    EXPECT_NEAR(0.0, zlacrm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zlacrm_test, zlacrm4) {
    EXPECT_NEAR(0.0, zlacrm_obj->diff, LAPACKE_EIG_THRESHOLD);
}
