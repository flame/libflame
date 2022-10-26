#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"




#define upgtr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if (d != NULL) free(d); \
if (e != NULL) free(e); \
if (dref!=NULL) free(dref);\
if (eref!=NULL) free(eref);\
if (q != NULL) free(q); \
if (qref != NULL) free(qref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;

/* Begin scomplex_common_parameters  class definition */
class upgtr_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	void *bModule, *lModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_float* A;
	char uplo;
	lapack_int ldq;
	/*Output Parameter*/
	lapack_complex_float* tau, *q;
	float* d;
	float* e;
	lapack_complex_float *Aref, *tauref, *qref;
	float* dref;
	float* eref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_hptrd, inforef_hptrd;

   public:
      upgtr_scomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int ldq);
      ~upgtr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
upgtr_scomplex_parameters:: upgtr_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int ldq_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	ldq = ldq_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n upgtr scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif

	bufsize = (n*(n+1)/2);
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (n-1));
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, n);
	lapacke_gtest_alloc_float_buffer_pair(&e, &eref, (n-1));
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&q, &qref, (ldq*n));
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL) ||
		(q == NULL) || (qref == NULL)){
		EXPECT_FALSE( true) << "upgtr_float_parameters object: malloc error.";
		upgtr_free();
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
	for(i=0;i<(ldq*n);i++) {
		q[i] = 0;
		qref[i] = q[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
upgtr_scomplex_parameters :: ~upgtr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" upgtr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   upgtr_free();

}
/*  Test fixture class definition */
class cupgtr_test  : public  ::testing::Test {
public:
   upgtr_scomplex_parameters  *cupgtr_obj;
   void SetUp();
   void TearDown () { delete cupgtr_obj;}
};

void cupgtr_test::SetUp(){

	/* LAPACKE chptrd prototype */
    typedef int (*Fptr_NL_LAPACKE_chptrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_float *A, 
											float* d, float* e, lapack_complex_float* tau);
	 Fptr_NL_LAPACKE_chptrd chptrd;
	 /* LAPACKE cupgtr prototype */
    typedef int (*Fptr_NL_LAPACKE_cupgtr) (int matrix_layout, char uplo, lapack_int n, const lapack_complex_float *A, 
											const lapack_complex_float* tau, lapack_complex_float* q, lapack_int ldq);

    Fptr_NL_LAPACKE_cupgtr cupgtr;

    cupgtr_obj = new upgtr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    cupgtr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cupgtr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cupgtr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cupgtr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cupgtr library call */
    cupgtr = (Fptr_NL_LAPACKE_cupgtr)dlsym(cupgtr_obj->hModule, "LAPACKE_cupgtr");
    ASSERT_TRUE(cupgtr != NULL) << "failed to get the Netlib LAPACKE_cupgtr symbol";

	/*chptrd library call*/
	cupgtr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cupgtr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cupgtr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cupgtr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    chptrd = (Fptr_NL_LAPACKE_chptrd)dlsym(cupgtr_obj->lModule, "LAPACKE_chptrd");
    ASSERT_TRUE(chptrd != NULL) << "failed to get the Netlib LAPACKE_chptrd symbol";    
    

    cupgtr_obj->inforef_hptrd = chptrd( cupgtr_obj->matrix_layout, cupgtr_obj->uplo,
								cupgtr_obj->n, cupgtr_obj->Aref,cupgtr_obj->dref, 
								cupgtr_obj->eref, cupgtr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cupgtr_obj->info_hptrd = LAPACKE_chptrd( cupgtr_obj->matrix_layout, cupgtr_obj->uplo,
										cupgtr_obj->n,cupgtr_obj->A, cupgtr_obj->d,
										cupgtr_obj->e,cupgtr_obj->tau);

    if( cupgtr_obj->info_hptrd < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chptrd is wrong\n", cupgtr_obj->info_hptrd );
    }
    if( cupgtr_obj->inforef_hptrd < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chptrd is wrong\n", 
        cupgtr_obj->inforef_hptrd );
    }  
/*Compute cupgtr's  o/p */
    cupgtr_obj->inforef = cupgtr( cupgtr_obj->matrix_layout, cupgtr_obj->uplo,
								cupgtr_obj->n, (const lapack_complex_float *)cupgtr_obj->Aref,(const lapack_complex_float *)cupgtr_obj->tauref, 
								cupgtr_obj->qref, cupgtr_obj->ldq);

    /* Compute libflame's Lapacke o/p  */
    cupgtr_obj->info = LAPACKE_cupgtr( cupgtr_obj->matrix_layout, cupgtr_obj->uplo,
										cupgtr_obj->n,(const lapack_complex_float *)cupgtr_obj->A, (const lapack_complex_float*)cupgtr_obj->tau,
										cupgtr_obj->q,cupgtr_obj->ldq);

    if( cupgtr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cupgtr is wrong\n", cupgtr_obj->info );
    }
    if( cupgtr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cupgtr is wrong\n", 
        cupgtr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cupgtr_obj->diff =  computeDiff_c( cupgtr_obj->bufsize, 
                cupgtr_obj->A, cupgtr_obj->Aref );

}

TEST_F(cupgtr_test, cupgtr1) {
    EXPECT_NEAR(0.0, cupgtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cupgtr_test, cupgtr2) {
    EXPECT_NEAR(0.0, cupgtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cupgtr_test, cupgtr3) {
    EXPECT_NEAR(0.0, cupgtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cupgtr_test, cupgtr4) {
    EXPECT_NEAR(0.0, cupgtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class upgtr_dcomplex_parameters{
	
	 public:
	int bufsize;
	void *hModule, *dModule;
	void *bModule, *lModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;	
	lapack_complex_double* A;
	char uplo;
	lapack_int ldq;
	/*Output Parameter*/
	lapack_complex_double* tau, *q;
	double* d;
	double* e;
	lapack_complex_double *Aref, *tauref, *qref;
	double* dref;
	double* eref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_hptrd, inforef_hptrd;

   public:
      upgtr_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int ldq);
      ~upgtr_dcomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
upgtr_dcomplex_parameters:: upgtr_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int ldq_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	ldq = ldq_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n upgtr scomplex:  n: %d lda: %d \n", n);
	#endif

	bufsize = (n *(n+1)/2);;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (n-1));
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, n);
	lapacke_gtest_alloc_double_buffer_pair(&e, &eref, (n-1));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&q, &qref, ldq*n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL) ||
		(d == NULL) || (dref == NULL) ||
		(e == NULL) || (eref == NULL) ||
		(q == NULL) || (qref == NULL)){
		EXPECT_FALSE( true) << "upgtr_float_parameters object: malloc error.";
		upgtr_free();
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
	/*initialize output matrix by 0 */
	for(i=0;i<(ldq*n);i++) {
		q[i] = 0;
		qref[i] = q[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
upgtr_dcomplex_parameters :: ~upgtr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" upgtr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   upgtr_free();

}
/*  Test fixture class definition */
class zupgtr_test  : public  ::testing::Test {
public:
   upgtr_dcomplex_parameters  *zupgtr_obj;
   void SetUp();
   void TearDown () { delete zupgtr_obj;}
};

void zupgtr_test::SetUp(){
	/* LAPACKE zhptrd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhptrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_double *A, 
											double* d, double* e, lapack_complex_double* tau);
	Fptr_NL_LAPACKE_zhptrd zhptrd;

    /* LAPACKE zupgtr prototype */
    typedef int (*Fptr_NL_LAPACKE_zupgtr) (int matrix_layout, char uplo, lapack_int n, const lapack_complex_double *A, 
											const lapack_complex_double* tau, lapack_complex_double* q, lapack_int ldq );

    Fptr_NL_LAPACKE_zupgtr zupgtr;

    zupgtr_obj = new upgtr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);
	
	zupgtr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zupgtr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
	
	ASSERT_TRUE(zupgtr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zupgtr_obj->hModule != NULL) << "Netlib lapacke handle NULL";
	
	zhptrd = (Fptr_NL_LAPACKE_zhptrd)dlsym(zupgtr_obj->hModule, "LAPACKE_zhptrd");
    ASSERT_TRUE(zhptrd != NULL) << "failed to get the Netlib LAPACKE_zhptrd symbol";
	
    zupgtr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zupgtr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zupgtr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zupgtr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zupgtr = (Fptr_NL_LAPACKE_zupgtr)dlsym(zupgtr_obj->hModule, "LAPACKE_zupgtr");
    ASSERT_TRUE(zupgtr != NULL) << "failed to get the Netlib LAPACKE_zupgtr symbol";
	/*compute Zhptrd o/p"*/
	zupgtr_obj->inforef_hptrd = zhptrd( zupgtr_obj->matrix_layout, zupgtr_obj->uplo,
								zupgtr_obj->n, zupgtr_obj->Aref, zupgtr_obj->dref,
								zupgtr_obj->eref, zupgtr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zupgtr_obj->info_hptrd = LAPACKE_zhptrd( zupgtr_obj->matrix_layout, zupgtr_obj->uplo,
										zupgtr_obj->n,zupgtr_obj->A, zupgtr_obj->d,
										zupgtr_obj->e,zupgtr_obj->tau);

    if( zupgtr_obj->info_hptrd < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhptrd is wrong\n", zupgtr_obj->info_hptrd );
    }
    if( zupgtr_obj->inforef_hptrd < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhptrd is wrong\n", 
        zupgtr_obj->inforef_hptrd );
    }
    
	/*coumpute zupgtr o/p */
    zupgtr_obj->inforef = zupgtr( zupgtr_obj->matrix_layout, zupgtr_obj->uplo,
								zupgtr_obj->n, (const lapack_complex_double*)zupgtr_obj->Aref,\
								(const lapack_complex_double*)zupgtr_obj->tauref, zupgtr_obj->qref, zupgtr_obj->ldq);

    /* Compute libflame's Lapacke o/p  */
    zupgtr_obj->info = LAPACKE_zupgtr(  zupgtr_obj->matrix_layout, zupgtr_obj->uplo,
								zupgtr_obj->n, (const lapack_complex_double* )zupgtr_obj->A,\
								(const lapack_complex_double*)zupgtr_obj->tau, zupgtr_obj->q, zupgtr_obj->ldq);

    if( zupgtr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zupgtr is wrong\n", zupgtr_obj->info );
    }
    if( zupgtr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zupgtr is wrong\n", 
        zupgtr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zupgtr_obj->diff =  computeDiff_z( zupgtr_obj->bufsize, 
                zupgtr_obj->A, zupgtr_obj->Aref );

}

TEST_F(zupgtr_test, zupgtr1) {
    EXPECT_NEAR(0.0, zupgtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zupgtr_test, zupgtr2) {
    EXPECT_NEAR(0.0, zupgtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zupgtr_test, zupgtr3) {
    EXPECT_NEAR(0.0, zupgtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zupgtr_test, zupgtr4) {
    EXPECT_NEAR(0.0, zupgtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

	
