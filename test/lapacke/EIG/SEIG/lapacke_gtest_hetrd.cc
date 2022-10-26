#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"




#define hetrd_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if (d != NULL) free(d); \
if (e != NULL) free(e); \
if (dref!=NULL) free(dref);\
if (eref!=NULL) free(eref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class hetrd_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A;
	lapack_int lda;
	char uplo;
	/*Output Parameter*/
	lapack_complex_float* tau;
	float* d;
	float* e;
	lapack_complex_float *Aref, *tauref;
	float* dref;
	float* eref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hetrd_scomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~hetrd_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hetrd_scomplex_parameters:: hetrd_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hetrd scomplex: matrix_layout = %d, uplo:%c, n: %d lda: %d \n", matrix_layout, uplo, n, lda);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (n-1));
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, n);
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, (n-1));
	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL)){
		EXPECT_FALSE( true) << "hetrd_float_parameters object: malloc error.";
		hetrd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(n-1);i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}
	/*initialize output matrix by 0 */
	for(i=0;i<(n-1);i++) {
		e[i] = 0;
		eref[i] = e[i];
	}
	/*initialize output matrix by 0 */
	for(i=0;i<n;i++) {
		d[i] = 0;
		dref[i] = d[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hetrd_scomplex_parameters :: ~hetrd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hetrd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hetrd_free();

}
/*  Test fixture class definition */
class chetrd_test  : public  ::testing::Test {
public:
   hetrd_scomplex_parameters  *chetrd_obj;
   void SetUp();
   void TearDown () { delete chetrd_obj; }
};

void chetrd_test::SetUp(){

    /* LAPACKE chetrd prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_float *A, 
											lapack_int lda, float* d, float* e, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_chetrd chetrd;

    chetrd_obj = new hetrd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chetrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chetrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chetrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chetrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chetrd = (Fptr_NL_LAPACKE_chetrd)dlsym(chetrd_obj->hModule, "LAPACKE_chetrd");
    ASSERT_TRUE(chetrd != NULL) << "failed to get the Netlib LAPACKE_chetrd symbol";
    

    chetrd_obj->inforef = chetrd( chetrd_obj->matrix_layout, chetrd_obj->uplo,
								chetrd_obj->n, chetrd_obj->Aref, chetrd_obj->lda, 
								chetrd_obj->dref, chetrd_obj->eref, chetrd_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    chetrd_obj->info = LAPACKE_chetrd( chetrd_obj->matrix_layout, chetrd_obj->uplo,
										chetrd_obj->n,chetrd_obj->A, chetrd_obj->lda, 
										chetrd_obj->d, chetrd_obj->e,chetrd_obj->tau);

    if( chetrd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chetrd is wrong\n", chetrd_obj->info );
    }
    if( chetrd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrd is wrong\n", 
        chetrd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chetrd_obj->diff =  computeDiff_c( chetrd_obj->bufsize, 
                chetrd_obj->A, chetrd_obj->Aref );

}

TEST_F(chetrd_test, chetrd1) {
    EXPECT_NEAR(0.0, chetrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetrd_test, chetrd2) {
    EXPECT_NEAR(0.0, chetrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetrd_test, chetrd3) {
    EXPECT_NEAR(0.0, chetrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chetrd_test, chetrd4) {
    EXPECT_NEAR(0.0, chetrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class hetrd_dcomplex_parameters{
	
	 public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* A;
	lapack_int lda;
	char uplo;
	/*Output Parameter*/
	lapack_complex_double* tau;
	double* d;
	double* e;
	lapack_complex_double *Aref, *tauref;
	double* dref;
	double* eref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      hetrd_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~hetrd_dcomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hetrd_dcomplex_parameters:: hetrd_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;	
	lda = lda_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hetrd scomplex:  n: %d lda: %d \n", n, lda);
	#endif

	bufsize = lda*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (n-1));
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, n);
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, (n-1));
	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL)){
		EXPECT_FALSE( true) << "hetrd_float_parameters object: malloc error.";
		hetrd_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(n-1);i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}
	/*initialize output matrix by 0 */
	for(i=0;i<(n-1);i++) {
		e[i] = 0;
		eref[i] = e[i];
	}
	/*initialize output matrix by 0 */
	for(i=0;i<n;i++) {
		d[i] = 0;
		dref[i] = d[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
hetrd_dcomplex_parameters :: ~hetrd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hetrd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hetrd_free();

}
/*  Test fixture class definition */
class zhetrd_test  : public  ::testing::Test {
public:
   hetrd_dcomplex_parameters  *zhetrd_obj;
   void SetUp();
   void TearDown () { delete zhetrd_obj; }
};

void zhetrd_test::SetUp(){

    /* LAPACKE zhetrd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_double *A, 
											lapack_int lda, double* d, double* e, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zhetrd zhetrd;

    zhetrd_obj = new hetrd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhetrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhetrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhetrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhetrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zhetrd = (Fptr_NL_LAPACKE_zhetrd)dlsym(zhetrd_obj->hModule, "LAPACKE_zhetrd");
    ASSERT_TRUE(zhetrd != NULL) << "failed to get the Netlib LAPACKE_zhetrd symbol";
    

    zhetrd_obj->inforef = zhetrd( zhetrd_obj->matrix_layout, zhetrd_obj->uplo,
								zhetrd_obj->n, zhetrd_obj->Aref, zhetrd_obj->lda, 
								zhetrd_obj->dref, zhetrd_obj->eref, zhetrd_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zhetrd_obj->info = LAPACKE_zhetrd( zhetrd_obj->matrix_layout, zhetrd_obj->uplo,
										zhetrd_obj->n,zhetrd_obj->A, zhetrd_obj->lda, 
										zhetrd_obj->d, zhetrd_obj->e,zhetrd_obj->tau);

    if( zhetrd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhetrd is wrong\n", zhetrd_obj->info );
    }
    if( zhetrd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrd is wrong\n", 
        zhetrd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhetrd_obj->diff =  computeDiff_z( zhetrd_obj->bufsize, 
                zhetrd_obj->A, zhetrd_obj->Aref );

}

TEST_F(zhetrd_test, zhetrd1) {
    EXPECT_NEAR(0.0, zhetrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetrd_test, zhetrd2) {
    EXPECT_NEAR(0.0, zhetrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetrd_test, zhetrd3) {
    EXPECT_NEAR(0.0, zhetrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhetrd_test, zhetrd4) {
    EXPECT_NEAR(0.0, zhetrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

	