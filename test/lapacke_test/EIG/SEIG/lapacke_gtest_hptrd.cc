#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"




#define hptrd_free() \
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
class hptrd_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A;
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
      hptrd_scomplex_parameters (int matrix_layout, char uplo, lapack_int n);
      ~hptrd_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hptrd_scomplex_parameters:: hptrd_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hptrd scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif

	bufsize = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (n-1));
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, n);
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, (n-1));
	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL)){
		EXPECT_FALSE( true) << "hptrd_float_parameters object: malloc error.";
		hptrd_free();
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
hptrd_scomplex_parameters :: ~hptrd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hptrd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hptrd_free();

}
/*  Test fixture class definition */
class chptrd_test  : public  ::testing::Test {
public:
   hptrd_scomplex_parameters  *chptrd_obj;
   void SetUp();
   void TearDown () { delete chptrd_obj; }
};

void chptrd_test::SetUp(){

    /* LAPACKE chptrd prototype */
    typedef int (*Fptr_NL_LAPACKE_chptrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_float *A, 
											float* d, float* e, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_chptrd chptrd;

    chptrd_obj = new hptrd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    chptrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    chptrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(chptrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(chptrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    chptrd = (Fptr_NL_LAPACKE_chptrd)dlsym(chptrd_obj->hModule, "LAPACKE_chptrd");
    ASSERT_TRUE(chptrd != NULL) << "failed to get the Netlib LAPACKE_chptrd symbol";
    

    chptrd_obj->inforef = chptrd( chptrd_obj->matrix_layout, chptrd_obj->uplo,
								chptrd_obj->n, chptrd_obj->Aref,chptrd_obj->dref, 
								chptrd_obj->eref, chptrd_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    chptrd_obj->info = LAPACKE_chptrd( chptrd_obj->matrix_layout, chptrd_obj->uplo,
										chptrd_obj->n,chptrd_obj->A, chptrd_obj->d,
										chptrd_obj->e,chptrd_obj->tau);

    if( chptrd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chptrd is wrong\n", chptrd_obj->info );
    }
    if( chptrd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chptrd is wrong\n", 
        chptrd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    chptrd_obj->diff =  computeDiff_c( chptrd_obj->bufsize, 
                chptrd_obj->A, chptrd_obj->Aref );

}

TEST_F(chptrd_test, chptrd1) {
    EXPECT_NEAR(0.0, chptrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chptrd_test, chptrd2) {
    EXPECT_NEAR(0.0, chptrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chptrd_test, chptrd3) {
    EXPECT_NEAR(0.0, chptrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(chptrd_test, chptrd4) {
    EXPECT_NEAR(0.0, chptrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class hptrd_dcomplex_parameters{
	
	 public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_complex_double* A;
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
      hptrd_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n);
      ~hptrd_dcomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
hptrd_dcomplex_parameters:: hptrd_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n hptrd scomplex:  n: %d lda: %d \n", n);
	#endif

	bufsize = (n *(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (n-1));
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, n);
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, (n-1));
	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL)){
		EXPECT_FALSE( true) << "hptrd_float_parameters object: malloc error.";
		hptrd_free();
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
hptrd_dcomplex_parameters :: ~hptrd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" hptrd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   hptrd_free();

}
/*  Test fixture class definition */
class zhptrd_test  : public  ::testing::Test {
public:
   hptrd_dcomplex_parameters  *zhptrd_obj;
   void SetUp();
   void TearDown () { delete zhptrd_obj; }
};

void zhptrd_test::SetUp(){

    /* LAPACKE zhptrd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_double *A, 
											double* d, double* e, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zhptrd zhptrd;

    zhptrd_obj = new hptrd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    zhptrd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zhptrd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zhptrd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zhptrd_obj->hModule != NULL) << "Netlib lapacke handle NULL";	


    zhptrd = (Fptr_NL_LAPACKE_zhptrd)dlsym(zhptrd_obj->hModule, "LAPACKE_zhptrd");
    ASSERT_TRUE(zhptrd != NULL) << "failed to get the Netlib LAPACKE_zhptrd symbol";
    

    zhptrd_obj->inforef = zhptrd( zhptrd_obj->matrix_layout, zhptrd_obj->uplo,
								zhptrd_obj->n, zhptrd_obj->Aref, zhptrd_obj->dref,
								zhptrd_obj->eref, zhptrd_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zhptrd_obj->info = LAPACKE_zhptrd( zhptrd_obj->matrix_layout, zhptrd_obj->uplo,
										zhptrd_obj->n,zhptrd_obj->A, zhptrd_obj->d,
										zhptrd_obj->e,zhptrd_obj->tau);

    if( zhptrd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhptrd is wrong\n", zhptrd_obj->info );
    }
    if( zhptrd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhptrd is wrong\n", 
        zhptrd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zhptrd_obj->diff =  computeDiff_z( zhptrd_obj->bufsize, 
                zhptrd_obj->A, zhptrd_obj->Aref );

}

TEST_F(zhptrd_test, zhptrd1) {
    EXPECT_NEAR(0.0, zhptrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhptrd_test, zhptrd2) {
    EXPECT_NEAR(0.0, zhptrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhptrd_test, zhptrd3) {
    EXPECT_NEAR(0.0, zhptrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zhptrd_test, zhptrd4) {
    EXPECT_NEAR(0.0, zhptrd_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

	