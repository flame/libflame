#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define geqrf_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class geqrf_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	float* A;
	lapack_int lda;
	/*Output Parameter*/
	float* tau;
	float *Aref, *tauref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqrf_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrf_float_parameters ();

};

/* Constructor definition  float_common_parameters */
geqrf_float_parameters:: geqrf_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrf float:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(float);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(float);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&tau, &tauref, (min(m,n)*sizeof(float)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqrf_float_parameters object: malloc error.";
		geqrf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(min(m,n));i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
geqrf_float_parameters :: ~geqrf_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrf_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrf_free();

}
/*  Test fixture class definition */
class sgeqrf_test  : public  ::testing::Test {
public:
   geqrf_float_parameters  *sgeqrf_obj;
   void SetUp();
   void TearDown () { delete sgeqrf_obj; }
};

void sgeqrf_test::SetUp(){

    /* LAPACKE sgeqrf prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeqrf) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, float* tau);

    Fptr_NL_LAPACKE_sgeqrf sgeqrf;

    sgeqrf_obj = new geqrf_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgeqrf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeqrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeqrf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeqrf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgeqrf = (Fptr_NL_LAPACKE_sgeqrf)dlsym(sgeqrf_obj->hModule, "LAPACKE_sgeqrf");
    ASSERT_TRUE(sgeqrf != NULL) << "failed to get the Netlib LAPACKE_sgeqrf symbol";
    

    sgeqrf_obj->inforef = sgeqrf( sgeqrf_obj->matrix_layout, sgeqrf_obj->m,
								sgeqrf_obj->n,sgeqrf_obj->Aref,
								sgeqrf_obj->lda, sgeqrf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    sgeqrf_obj->info = LAPACKE_sgeqrf( sgeqrf_obj->matrix_layout, sgeqrf_obj->m,
										sgeqrf_obj->n,sgeqrf_obj->A, 
										sgeqrf_obj->lda, sgeqrf_obj->tau);

    if( sgeqrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgeqrf is wrong\n", sgeqrf_obj->info );
    }
    if( sgeqrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeqrf is wrong\n", 
        sgeqrf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgeqrf_obj->diff =  computeDiff_s( sgeqrf_obj->bufsize, 
                sgeqrf_obj->A, sgeqrf_obj->Aref );

}

TEST_F(sgeqrf_test, sgeqrf1) {
    EXPECT_NEAR(0.0, sgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqrf_test, sgeqrf2) {
    EXPECT_NEAR(0.0, sgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqrf_test, sgeqrf3) {
    EXPECT_NEAR(0.0, sgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqrf_test, sgeqrf4) {
    EXPECT_NEAR(0.0, sgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class geqrf_double_parameters{

   public:
	int bufsize;
	double diff;
	void *hModule, *dModule;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	double* A;
	lapack_int lda;
	/*Output Parameter*/
	double* tau;
	double *Aref, *tauref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqrf_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrf_double_parameters ();

};

/* Constructor definition  double_common_parameters */
geqrf_double_parameters:: geqrf_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrf double:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(double);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(double);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&tau, &tauref, (min(m,n)*sizeof(double)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqrf_double_parameters object: malloc error.";
		geqrf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(min(m,n));i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
geqrf_double_parameters :: ~geqrf_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrf_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrf_free();

}
/*  Test fixture class definition */
class dgeqrf_test  : public  ::testing::Test {
public:
   geqrf_double_parameters  *dgeqrf_obj;
   void SetUp();
   void TearDown () { delete dgeqrf_obj; }
};

void dgeqrf_test::SetUp(){

    /* LAPACKE dgeqrf prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeqrf) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, double* tau);

    Fptr_NL_LAPACKE_dgeqrf dgeqrf;

    dgeqrf_obj = new geqrf_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);
    idx = Circular_Increment_Index(idx);
    dgeqrf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeqrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeqrf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeqrf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgeqrf = (Fptr_NL_LAPACKE_dgeqrf)dlsym(dgeqrf_obj->hModule, "LAPACKE_dgeqrf");
    ASSERT_TRUE(dgeqrf != NULL) << "failed to get the Netlib LAPACKE_dgeqrf symbol";
    

    dgeqrf_obj->inforef = dgeqrf( dgeqrf_obj->matrix_layout, dgeqrf_obj->m,
								dgeqrf_obj->n,dgeqrf_obj->Aref,
								dgeqrf_obj->lda, dgeqrf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dgeqrf_obj->info = LAPACKE_dgeqrf( dgeqrf_obj->matrix_layout, dgeqrf_obj->m,
										dgeqrf_obj->n,dgeqrf_obj->A, 
										dgeqrf_obj->lda, dgeqrf_obj->tau);

    if( dgeqrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgeqrf is wrong\n", dgeqrf_obj->info );
    }
    if( dgeqrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeqrf is wrong\n", 
        dgeqrf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgeqrf_obj->diff =  computeDiff_d( dgeqrf_obj->bufsize, 
                dgeqrf_obj->A, dgeqrf_obj->Aref );

}

TEST_F(dgeqrf_test, dgeqrf1) {
    EXPECT_NEAR(0.0, dgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqrf_test, dgeqrf2) {
    EXPECT_NEAR(0.0, dgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqrf_test, dgeqrf3) {
    EXPECT_NEAR(0.0, dgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqrf_test, dgeqrf4) {
    EXPECT_NEAR(0.0, dgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class geqrf_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_complex_float* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float* tau;
	lapack_complex_float *Aref, *tauref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqrf_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrf_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
geqrf_scomplex_parameters:: geqrf_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrf scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(lapack_complex_float);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(lapack_complex_float);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (min(m,n)*sizeof(lapack_complex_float)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqrf_float_parameters object: malloc error.";
		geqrf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(min(m,n));i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
geqrf_scomplex_parameters :: ~geqrf_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrf_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrf_free();

}
/*  Test fixture class definition */
class cgeqrf_test  : public  ::testing::Test {
public:
   geqrf_scomplex_parameters  *cgeqrf_obj;
   void SetUp();
   void TearDown () { delete cgeqrf_obj; }
};

void cgeqrf_test::SetUp(){

    /* LAPACKE cgeqrf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeqrf) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cgeqrf cgeqrf;

    cgeqrf_obj = new geqrf_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cgeqrf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeqrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeqrf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeqrf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgeqrf = (Fptr_NL_LAPACKE_cgeqrf)dlsym(cgeqrf_obj->hModule, "LAPACKE_cgeqrf");
    ASSERT_TRUE(cgeqrf != NULL) << "failed to get the Netlib LAPACKE_cgeqrf symbol";
    

    cgeqrf_obj->inforef = cgeqrf( cgeqrf_obj->matrix_layout, cgeqrf_obj->m,
								cgeqrf_obj->n,cgeqrf_obj->Aref,
								cgeqrf_obj->lda, cgeqrf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cgeqrf_obj->info = LAPACKE_cgeqrf( cgeqrf_obj->matrix_layout, cgeqrf_obj->m,
										cgeqrf_obj->n,cgeqrf_obj->A, 
										cgeqrf_obj->lda, cgeqrf_obj->tau);

    if( cgeqrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgeqrf is wrong\n", cgeqrf_obj->info );
    }
    if( cgeqrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeqrf is wrong\n", 
        cgeqrf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgeqrf_obj->diff =  computeDiff_c( cgeqrf_obj->bufsize, 
                cgeqrf_obj->A, cgeqrf_obj->Aref );

}

TEST_F(cgeqrf_test, cgeqrf1) {
    EXPECT_NEAR(0.0, cgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqrf_test, cgeqrf2) {
    EXPECT_NEAR(0.0, cgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqrf_test, cgeqrf3) {
    EXPECT_NEAR(0.0, cgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqrf_test, cgeqrf4) {
    EXPECT_NEAR(0.0, cgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class geqrf_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_complex_double* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double* tau;
	lapack_complex_double *Aref, *tauref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqrf_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrf_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
geqrf_dcomplex_parameters:: geqrf_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrf dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(lapack_complex_double);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(lapack_complex_double);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (min(m,n)*sizeof(lapack_complex_double)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqrf_float_parameters object: malloc error.";
		geqrf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(min(m,n));i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
geqrf_dcomplex_parameters :: ~geqrf_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrf_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrf_free();

}
/*  Test fixture class definition */
class zgeqrf_test  : public  ::testing::Test {
public:
   geqrf_dcomplex_parameters  *zgeqrf_obj;
   void SetUp();
   void TearDown () { delete zgeqrf_obj; }
};

void zgeqrf_test::SetUp(){

    /* LAPACKE zgeqrf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeqrf) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zgeqrf zgeqrf;

    zgeqrf_obj = new geqrf_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zgeqrf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeqrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeqrf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeqrf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgeqrf = (Fptr_NL_LAPACKE_zgeqrf)dlsym(zgeqrf_obj->hModule, "LAPACKE_zgeqrf");
    ASSERT_TRUE(zgeqrf != NULL) << "failed to get the Netlib LAPACKE_zgeqrf symbol";
    

    zgeqrf_obj->inforef = zgeqrf( zgeqrf_obj->matrix_layout, zgeqrf_obj->m,
								zgeqrf_obj->n,zgeqrf_obj->Aref,
								zgeqrf_obj->lda, zgeqrf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zgeqrf_obj->info = LAPACKE_zgeqrf( zgeqrf_obj->matrix_layout, zgeqrf_obj->m,
										zgeqrf_obj->n,zgeqrf_obj->A, 
										zgeqrf_obj->lda, zgeqrf_obj->tau);

    if( zgeqrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgeqrf is wrong\n", zgeqrf_obj->info );
    }
    if( zgeqrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeqrf is wrong\n", 
        zgeqrf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgeqrf_obj->diff =  computeDiff_z( zgeqrf_obj->bufsize, 
                zgeqrf_obj->A, zgeqrf_obj->Aref );

}

TEST_F(zgeqrf_test, zgeqrf1) {
    EXPECT_NEAR(0.0, zgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqrf_test, zgeqrf2) {
    EXPECT_NEAR(0.0, zgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqrf_test, zgeqrf3) {
    EXPECT_NEAR(0.0, zgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqrf_test, zgeqrf4) {
    EXPECT_NEAR(0.0, zgeqrf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


