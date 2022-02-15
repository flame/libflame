#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"




#define ungtr_free() \
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
class ungtr_scomplex_parameters{

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
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float* tau;
	float* d;
	float* e;
	lapack_complex_float *Aref, *tauref;
	float* dref;
	float* eref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_hetrd, inforef_hetrd;

   public:
      ungtr_scomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~ungtr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
ungtr_scomplex_parameters:: ungtr_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ungtr scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
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
		EXPECT_FALSE( true) << "ungtr_float_parameters object: malloc error.";
		ungtr_free();
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
ungtr_scomplex_parameters :: ~ungtr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ungtr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ungtr_free();

}
/*  Test fixture class definition */
class cungtr_test  : public  ::testing::Test {
public:
   ungtr_scomplex_parameters  *cungtr_obj;
   void SetUp();
   void TearDown () { delete cungtr_obj;}
};

void cungtr_test::SetUp(){

	/* LAPACKE chetrd prototype */
    typedef int (*Fptr_NL_LAPACKE_chetrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_float *A, lapack_int lda, \
											float* d, float* e, lapack_complex_float* tau);
	 Fptr_NL_LAPACKE_chetrd chetrd;
	 /* LAPACKE cungtr prototype */
    typedef int (*Fptr_NL_LAPACKE_cungtr) (int matrix_layout, char uplo, lapack_int n, lapack_complex_float *A, 
											lapack_int lda, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cungtr cungtr;

    cungtr_obj = new ungtr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);

    cungtr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cungtr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cungtr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cungtr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*cungtr library call */
    cungtr = (Fptr_NL_LAPACKE_cungtr)dlsym(cungtr_obj->hModule, "LAPACKE_cungtr");
    ASSERT_TRUE(cungtr != NULL) << "failed to get the Netlib LAPACKE_cungtr symbol";

	/*chetrd library call*/
	cungtr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cungtr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cungtr_obj->bModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cungtr_obj->lModule != NULL) << "Netlib lapacke handle NULL";

    chetrd = (Fptr_NL_LAPACKE_chetrd)dlsym(cungtr_obj->lModule, "LAPACKE_chetrd");
    ASSERT_TRUE(chetrd != NULL) << "failed to get the Netlib LAPACKE_chetrd symbol";    
    

    cungtr_obj->inforef_hetrd = chetrd( cungtr_obj->matrix_layout, cungtr_obj->uplo,
								cungtr_obj->n, cungtr_obj->Aref, cungtr_obj->lda, cungtr_obj->dref, 
								cungtr_obj->eref, cungtr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cungtr_obj->info_hetrd = LAPACKE_chetrd( cungtr_obj->matrix_layout, cungtr_obj->uplo,
										cungtr_obj->n,cungtr_obj->A, cungtr_obj->lda, cungtr_obj->d,
										cungtr_obj->e,cungtr_obj->tau);

    if( cungtr_obj->info_hetrd < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_chetrd is wrong\n", cungtr_obj->info_hetrd );
    }
    if( cungtr_obj->inforef_hetrd < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_chetrd is wrong\n", 
        cungtr_obj->inforef_hetrd );
    }  
/*Compute cungtr's  o/p */
    cungtr_obj->inforef = cungtr( cungtr_obj->matrix_layout, cungtr_obj->uplo,cungtr_obj->n,
								cungtr_obj->Aref,cungtr_obj->lda, cungtr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cungtr_obj->info = LAPACKE_cungtr( cungtr_obj->matrix_layout, cungtr_obj->uplo,
										cungtr_obj->n,cungtr_obj->A, cungtr_obj->lda, cungtr_obj->tau);

    if( cungtr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cungtr is wrong\n", cungtr_obj->info );
    }
    if( cungtr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cungtr is wrong\n", 
        cungtr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cungtr_obj->diff =  computeDiff_c( cungtr_obj->bufsize, 
                cungtr_obj->A, cungtr_obj->Aref );

}

TEST_F(cungtr_test, cungtr1) {
    EXPECT_NEAR(0.0, cungtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungtr_test, cungtr2) {
    EXPECT_NEAR(0.0, cungtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungtr_test, cungtr3) {
    EXPECT_NEAR(0.0, cungtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cungtr_test, cungtr4) {
    EXPECT_NEAR(0.0, cungtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class ungtr_dcomplex_parameters{
	
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
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double* tau;
	double* d;
	double* e;
	lapack_complex_double *Aref, *tauref;
	double* dref;
	double* eref;
	/*Return Values*/	
	lapack_int info, inforef;
	lapack_int info_hetrd, inforef_hetrd;

   public:
      ungtr_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_int lda);
      ~ungtr_dcomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
ungtr_dcomplex_parameters:: ungtr_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ungtr scomplex:  n: %d lda: %d \n", n);
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
		EXPECT_FALSE( true) << "ungtr_float_parameters object: malloc error.";
		ungtr_free();
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
ungtr_dcomplex_parameters :: ~ungtr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ungtr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ungtr_free();

}
/*  Test fixture class definition */
class zungtr_test  : public  ::testing::Test {
public:
   ungtr_dcomplex_parameters  *zungtr_obj;
   void SetUp();
   void TearDown () { delete zungtr_obj;}
};

void zungtr_test::SetUp(){
	/* LAPACKE zhetrd prototype */
    typedef int (*Fptr_NL_LAPACKE_zhetrd) (int matrix_layout, char uplo, lapack_int n,lapack_complex_double *A, 
											lapack_int lda, double* d, double* e, lapack_complex_double* tau);
	Fptr_NL_LAPACKE_zhetrd zhetrd;

    /* LAPACKE zungtr prototype */
    typedef int (*Fptr_NL_LAPACKE_zungtr) (int matrix_layout, char uplo, lapack_int n, lapack_complex_double *A, 
											lapack_int lda, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zungtr zungtr;

    zungtr_obj = new ungtr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].uplo,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].n);

    idx = Circular_Increment_Index(idx);
	
	zungtr_obj->bModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zungtr_obj->lModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);
	
	ASSERT_TRUE(zungtr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zungtr_obj->hModule != NULL) << "Netlib lapacke handle NULL";
	
	zhetrd = (Fptr_NL_LAPACKE_zhetrd)dlsym(zungtr_obj->hModule, "LAPACKE_zhetrd");
    ASSERT_TRUE(zhetrd != NULL) << "failed to get the Netlib LAPACKE_zhetrd symbol";
	
    zungtr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zungtr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zungtr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zungtr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zungtr = (Fptr_NL_LAPACKE_zungtr)dlsym(zungtr_obj->hModule, "LAPACKE_zungtr");
    ASSERT_TRUE(zungtr != NULL) << "failed to get the Netlib LAPACKE_zungtr symbol";
	/*compute Zhetrd o/p"*/
	zungtr_obj->inforef_hetrd = zhetrd( zungtr_obj->matrix_layout, zungtr_obj->uplo,
								zungtr_obj->n, zungtr_obj->Aref, zungtr_obj->lda, zungtr_obj->dref,
								zungtr_obj->eref, zungtr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zungtr_obj->info_hetrd = LAPACKE_zhetrd( zungtr_obj->matrix_layout, zungtr_obj->uplo,
										zungtr_obj->n,zungtr_obj->A, zungtr_obj->lda, zungtr_obj->d,
										zungtr_obj->e,zungtr_obj->tau);

    if( zungtr_obj->info_hetrd < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zhetrd is wrong\n", zungtr_obj->info_hetrd );
    }
    if( zungtr_obj->inforef_hetrd < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zhetrd is wrong\n", 
        zungtr_obj->inforef_hetrd );
    }
    
	/*coumpute zungtr o/p */
    zungtr_obj->inforef = zungtr( zungtr_obj->matrix_layout, zungtr_obj->uplo,zungtr_obj->n, \
								zungtr_obj->Aref, zungtr_obj->lda, zungtr_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zungtr_obj->info = LAPACKE_zungtr(  zungtr_obj->matrix_layout, zungtr_obj->uplo,
								zungtr_obj->n,zungtr_obj->A, zungtr_obj->lda, \
								zungtr_obj->tau);

    if( zungtr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zungtr is wrong\n", zungtr_obj->info );
    }
    if( zungtr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zungtr is wrong\n", 
        zungtr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zungtr_obj->diff =  computeDiff_z( zungtr_obj->bufsize, 
                zungtr_obj->A, zungtr_obj->Aref );

}

TEST_F(zungtr_test, zungtr1) {
    EXPECT_NEAR(0.0, zungtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungtr_test, zungtr2) {
    EXPECT_NEAR(0.0, zungtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungtr_test, zungtr3) {
    EXPECT_NEAR(0.0, zungtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zungtr_test, zungtr4) {
    EXPECT_NEAR(0.0, zungtr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

	
