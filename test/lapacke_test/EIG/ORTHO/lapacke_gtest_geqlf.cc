#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define geqlf_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class geqlf_float_parameters{

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
      geqlf_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqlf_float_parameters ();

};

/* Constructor definition  float_common_parameters */
geqlf_float_parameters:: geqlf_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqlf float:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "geqlf_float_parameters object: malloc error.";
		geqlf_free();
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
geqlf_float_parameters :: ~geqlf_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqlf_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqlf_free();

}
/*  Test fixture class definition */
class sgeqlf_test  : public  ::testing::Test {
public:
   geqlf_float_parameters  *sgeqlf_obj;
   void SetUp();
   void TearDown () { delete sgeqlf_obj; }
};

void sgeqlf_test::SetUp(){

    /* LAPACKE sgeqlf prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeqlf) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, float* tau);

    Fptr_NL_LAPACKE_sgeqlf sgeqlf;

    sgeqlf_obj = new geqlf_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgeqlf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeqlf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeqlf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeqlf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgeqlf = (Fptr_NL_LAPACKE_sgeqlf)dlsym(sgeqlf_obj->hModule, "LAPACKE_sgeqlf");
    ASSERT_TRUE(sgeqlf != NULL) << "failed to get the Netlib LAPACKE_sgeqlf symbol";
    

    sgeqlf_obj->inforef = sgeqlf( sgeqlf_obj->matrix_layout, sgeqlf_obj->m,
								sgeqlf_obj->n,sgeqlf_obj->Aref,
								sgeqlf_obj->lda, sgeqlf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    sgeqlf_obj->info = LAPACKE_sgeqlf( sgeqlf_obj->matrix_layout, sgeqlf_obj->m,
										sgeqlf_obj->n,sgeqlf_obj->A, 
										sgeqlf_obj->lda, sgeqlf_obj->tau);

    if( sgeqlf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgeqlf is wrong\n", sgeqlf_obj->info );
    }
    if( sgeqlf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeqlf is wrong\n", 
        sgeqlf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgeqlf_obj->diff =  computeDiff_s( sgeqlf_obj->bufsize, 
                sgeqlf_obj->A, sgeqlf_obj->Aref );

}

TEST_F(sgeqlf_test, sgeqlf1) {
    EXPECT_NEAR(0.0, sgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqlf_test, sgeqlf2) {
    EXPECT_NEAR(0.0, sgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqlf_test, sgeqlf3) {
    EXPECT_NEAR(0.0, sgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqlf_test, sgeqlf4) {
    EXPECT_NEAR(0.0, sgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class geqlf_double_parameters{

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
      geqlf_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqlf_double_parameters ();

};

/* Constructor definition  double_common_parameters */
geqlf_double_parameters:: geqlf_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqlf double:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "geqlf_double_parameters object: malloc error.";
		geqlf_free();
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
geqlf_double_parameters :: ~geqlf_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqlf_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqlf_free();

}
/*  Test fixture class definition */
class dgeqlf_test  : public  ::testing::Test {
public:
   geqlf_double_parameters  *dgeqlf_obj;
   void SetUp();
   void TearDown () { delete dgeqlf_obj; }
};

void dgeqlf_test::SetUp(){

    /* LAPACKE dgeqlf prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeqlf) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, double* tau);

    Fptr_NL_LAPACKE_dgeqlf dgeqlf;

    dgeqlf_obj = new geqlf_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);
    idx = Circular_Increment_Index(idx);
    dgeqlf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeqlf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeqlf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeqlf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgeqlf = (Fptr_NL_LAPACKE_dgeqlf)dlsym(dgeqlf_obj->hModule, "LAPACKE_dgeqlf");
    ASSERT_TRUE(dgeqlf != NULL) << "failed to get the Netlib LAPACKE_dgeqlf symbol";
    

    dgeqlf_obj->inforef = dgeqlf( dgeqlf_obj->matrix_layout, dgeqlf_obj->m,
								dgeqlf_obj->n,dgeqlf_obj->Aref,
								dgeqlf_obj->lda, dgeqlf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dgeqlf_obj->info = LAPACKE_dgeqlf( dgeqlf_obj->matrix_layout, dgeqlf_obj->m,
										dgeqlf_obj->n,dgeqlf_obj->A, 
										dgeqlf_obj->lda, dgeqlf_obj->tau);

    if( dgeqlf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgeqlf is wrong\n", dgeqlf_obj->info );
    }
    if( dgeqlf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeqlf is wrong\n", 
        dgeqlf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgeqlf_obj->diff =  computeDiff_d( dgeqlf_obj->bufsize, 
                dgeqlf_obj->A, dgeqlf_obj->Aref );

}

TEST_F(dgeqlf_test, dgeqlf1) {
    EXPECT_NEAR(0.0, dgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqlf_test, dgeqlf2) {
    EXPECT_NEAR(0.0, dgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqlf_test, dgeqlf3) {
    EXPECT_NEAR(0.0, dgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqlf_test, dgeqlf4) {
    EXPECT_NEAR(0.0, dgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class geqlf_scomplex_parameters{

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
      geqlf_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqlf_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
geqlf_scomplex_parameters:: geqlf_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqlf scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "geqlf_float_parameters object: malloc error.";
		geqlf_free();
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
geqlf_scomplex_parameters :: ~geqlf_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqlf_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqlf_free();

}
/*  Test fixture class definition */
class cgeqlf_test  : public  ::testing::Test {
public:
   geqlf_scomplex_parameters  *cgeqlf_obj;
   void SetUp();
   void TearDown () { delete cgeqlf_obj; }
};

void cgeqlf_test::SetUp(){

    /* LAPACKE cgeqlf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeqlf) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cgeqlf cgeqlf;

    cgeqlf_obj = new geqlf_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cgeqlf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeqlf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeqlf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeqlf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgeqlf = (Fptr_NL_LAPACKE_cgeqlf)dlsym(cgeqlf_obj->hModule, "LAPACKE_cgeqlf");
    ASSERT_TRUE(cgeqlf != NULL) << "failed to get the Netlib LAPACKE_cgeqlf symbol";
    

    cgeqlf_obj->inforef = cgeqlf( cgeqlf_obj->matrix_layout, cgeqlf_obj->m,
								cgeqlf_obj->n,cgeqlf_obj->Aref,
								cgeqlf_obj->lda, cgeqlf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cgeqlf_obj->info = LAPACKE_cgeqlf( cgeqlf_obj->matrix_layout, cgeqlf_obj->m,
										cgeqlf_obj->n,cgeqlf_obj->A, 
										cgeqlf_obj->lda, cgeqlf_obj->tau);

    if( cgeqlf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgeqlf is wrong\n", cgeqlf_obj->info );
    }
    if( cgeqlf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeqlf is wrong\n", 
        cgeqlf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgeqlf_obj->diff =  computeDiff_c( cgeqlf_obj->bufsize, 
                cgeqlf_obj->A, cgeqlf_obj->Aref );

}

TEST_F(cgeqlf_test, cgeqlf1) {
    EXPECT_NEAR(0.0, cgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqlf_test, cgeqlf2) {
    EXPECT_NEAR(0.0, cgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqlf_test, cgeqlf3) {
    EXPECT_NEAR(0.0, cgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqlf_test, cgeqlf4) {
    EXPECT_NEAR(0.0, cgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class geqlf_dcomplex_parameters{

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
      geqlf_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqlf_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
geqlf_dcomplex_parameters:: geqlf_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqlf dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "geqlf_float_parameters object: malloc error.";
		geqlf_free();
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
geqlf_dcomplex_parameters :: ~geqlf_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqlf_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqlf_free();

}
/*  Test fixture class definition */
class zgeqlf_test  : public  ::testing::Test {
public:
   geqlf_dcomplex_parameters  *zgeqlf_obj;
   void SetUp();
   void TearDown () { delete zgeqlf_obj; }
};

void zgeqlf_test::SetUp(){

    /* LAPACKE zgeqlf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeqlf) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zgeqlf zgeqlf;

    zgeqlf_obj = new geqlf_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zgeqlf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeqlf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeqlf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeqlf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgeqlf = (Fptr_NL_LAPACKE_zgeqlf)dlsym(zgeqlf_obj->hModule, "LAPACKE_zgeqlf");
    ASSERT_TRUE(zgeqlf != NULL) << "failed to get the Netlib LAPACKE_zgeqlf symbol";
    

    zgeqlf_obj->inforef = zgeqlf( zgeqlf_obj->matrix_layout, zgeqlf_obj->m,
								zgeqlf_obj->n,zgeqlf_obj->Aref,
								zgeqlf_obj->lda, zgeqlf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zgeqlf_obj->info = LAPACKE_zgeqlf( zgeqlf_obj->matrix_layout, zgeqlf_obj->m,
										zgeqlf_obj->n,zgeqlf_obj->A, 
										zgeqlf_obj->lda, zgeqlf_obj->tau);

    if( zgeqlf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgeqlf is wrong\n", zgeqlf_obj->info );
    }
    if( zgeqlf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeqlf is wrong\n", 
        zgeqlf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgeqlf_obj->diff =  computeDiff_z( zgeqlf_obj->bufsize, 
                zgeqlf_obj->A, zgeqlf_obj->Aref );

}

TEST_F(zgeqlf_test, zgeqlf1) {
    EXPECT_NEAR(0.0, zgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqlf_test, zgeqlf2) {
    EXPECT_NEAR(0.0, zgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqlf_test, zgeqlf3) {
    EXPECT_NEAR(0.0, zgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqlf_test, zgeqlf4) {
    EXPECT_NEAR(0.0, zgeqlf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


