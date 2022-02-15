#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"


#define geqrfp_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class geqrfp_float_parameters{

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
      geqrfp_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrfp_float_parameters ();

};

/* Constructor definition  float_common_parameters */
geqrfp_float_parameters:: geqrfp_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrfp float:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&tau, &tauref, (min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqrfp_float_parameters object: malloc error.";
		geqrfp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
geqrfp_float_parameters :: ~geqrfp_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrfp_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrfp_free();

}
/*  Test fixture class definition */
class sgeqrfp_test  : public  ::testing::Test {
public:
   geqrfp_float_parameters  *sgeqrfp_obj;
   void SetUp();
   void TearDown () { delete sgeqrfp_obj; }
};

void sgeqrfp_test::SetUp(){

    /* LAPACKE sgeqrfp prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeqrfp) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, float* tau);

    Fptr_NL_LAPACKE_sgeqrfp sgeqrfp;

    sgeqrfp_obj = new geqrfp_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgeqrfp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeqrfp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeqrfp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeqrfp_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgeqrfp = (Fptr_NL_LAPACKE_sgeqrfp)dlsym(sgeqrfp_obj->hModule, "LAPACKE_sgeqrfp");
    ASSERT_TRUE(sgeqrfp != NULL) << "failed to get the Netlib LAPACKE_sgeqrfp symbol";
    

    sgeqrfp_obj->inforef = sgeqrfp( sgeqrfp_obj->matrix_layout, sgeqrfp_obj->m,
								sgeqrfp_obj->n,sgeqrfp_obj->Aref,
								sgeqrfp_obj->lda, sgeqrfp_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    sgeqrfp_obj->info = LAPACKE_sgeqrfp( sgeqrfp_obj->matrix_layout, sgeqrfp_obj->m,
										sgeqrfp_obj->n,sgeqrfp_obj->A, 
										sgeqrfp_obj->lda, sgeqrfp_obj->tau);

    if( sgeqrfp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgeqrfp is wrong\n", sgeqrfp_obj->info );
    }
    if( sgeqrfp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeqrfp is wrong\n", 
        sgeqrfp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgeqrfp_obj->diff =  computeDiff_s( sgeqrfp_obj->bufsize, 
                sgeqrfp_obj->A, sgeqrfp_obj->Aref );

}

TEST_F(sgeqrfp_test, sgeqrfp1) {
    EXPECT_NEAR(0.0, sgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgeqrfp_test, sgeqrfp2) {
    EXPECT_NEAR(0.0, sgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgeqrfp_test, sgeqrfp3) {
    EXPECT_NEAR(0.0, sgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgeqrfp_test, sgeqrfp4) {
    EXPECT_NEAR(0.0, sgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class geqrfp_double_parameters{

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
      geqrfp_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrfp_double_parameters ();

};

/* Constructor definition  double_common_parameters */
geqrfp_double_parameters:: geqrfp_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrfp double:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&tau, &tauref, (min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqrfp_double_parameters object: malloc error.";
		geqrfp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
geqrfp_double_parameters :: ~geqrfp_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrfp_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrfp_free();

}
/*  Test fixture class definition */
class dgeqrfp_test  : public  ::testing::Test {
public:
   geqrfp_double_parameters  *dgeqrfp_obj;
   void SetUp();
   void TearDown () { delete dgeqrfp_obj; }
};

void dgeqrfp_test::SetUp(){

    /* LAPACKE dgeqrfp prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeqrfp) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, double* tau);

    Fptr_NL_LAPACKE_dgeqrfp dgeqrfp;

    dgeqrfp_obj = new geqrfp_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);
    idx = Circular_Increment_Index(idx);
    dgeqrfp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeqrfp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeqrfp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeqrfp_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgeqrfp = (Fptr_NL_LAPACKE_dgeqrfp)dlsym(dgeqrfp_obj->hModule, "LAPACKE_dgeqrfp");
    ASSERT_TRUE(dgeqrfp != NULL) << "failed to get the Netlib LAPACKE_dgeqrfp symbol";
    

    dgeqrfp_obj->inforef = dgeqrfp( dgeqrfp_obj->matrix_layout, dgeqrfp_obj->m,
								dgeqrfp_obj->n,dgeqrfp_obj->Aref,
								dgeqrfp_obj->lda, dgeqrfp_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dgeqrfp_obj->info = LAPACKE_dgeqrfp( dgeqrfp_obj->matrix_layout, dgeqrfp_obj->m,
										dgeqrfp_obj->n,dgeqrfp_obj->A, 
										dgeqrfp_obj->lda, dgeqrfp_obj->tau);

    if( dgeqrfp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgeqrfp is wrong\n", dgeqrfp_obj->info );
    }
    if( dgeqrfp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeqrfp is wrong\n", 
        dgeqrfp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgeqrfp_obj->diff =  computeDiff_d( dgeqrfp_obj->bufsize, 
                dgeqrfp_obj->A, dgeqrfp_obj->Aref );

}

TEST_F(dgeqrfp_test, dgeqrfp1) {
    EXPECT_NEAR(0.0, dgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgeqrfp_test, dgeqrfp2) {
    EXPECT_NEAR(0.0, dgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgeqrfp_test, dgeqrfp3) {
    EXPECT_NEAR(0.0, dgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgeqrfp_test, dgeqrfp4) {
    EXPECT_NEAR(0.0, dgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class geqrfp_scomplex_parameters{

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
      geqrfp_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrfp_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
geqrfp_scomplex_parameters:: geqrfp_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrfp scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqrfp_float_parameters object: malloc error.";
		geqrfp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
geqrfp_scomplex_parameters :: ~geqrfp_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrfp_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrfp_free();

}
/*  Test fixture class definition */
class cgeqrfp_test  : public  ::testing::Test {
public:
   geqrfp_scomplex_parameters  *cgeqrfp_obj;
   void SetUp();
   void TearDown () { delete cgeqrfp_obj; }
};

void cgeqrfp_test::SetUp(){

    /* LAPACKE cgeqrfp prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeqrfp) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cgeqrfp cgeqrfp;

    cgeqrfp_obj = new geqrfp_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cgeqrfp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeqrfp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeqrfp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeqrfp_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgeqrfp = (Fptr_NL_LAPACKE_cgeqrfp)dlsym(cgeqrfp_obj->hModule, "LAPACKE_cgeqrfp");
    ASSERT_TRUE(cgeqrfp != NULL) << "failed to get the Netlib LAPACKE_cgeqrfp symbol";
    

    cgeqrfp_obj->inforef = cgeqrfp( cgeqrfp_obj->matrix_layout, cgeqrfp_obj->m,
								cgeqrfp_obj->n,cgeqrfp_obj->Aref,
								cgeqrfp_obj->lda, cgeqrfp_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cgeqrfp_obj->info = LAPACKE_cgeqrfp( cgeqrfp_obj->matrix_layout, cgeqrfp_obj->m,
										cgeqrfp_obj->n,cgeqrfp_obj->A, 
										cgeqrfp_obj->lda, cgeqrfp_obj->tau);

    if( cgeqrfp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgeqrfp is wrong\n", cgeqrfp_obj->info );
    }
    if( cgeqrfp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeqrfp is wrong\n", 
        cgeqrfp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgeqrfp_obj->diff =  computeDiff_c( cgeqrfp_obj->bufsize, 
                cgeqrfp_obj->A, cgeqrfp_obj->Aref );

}

TEST_F(cgeqrfp_test, cgeqrfp1) {
    EXPECT_NEAR(0.0, cgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgeqrfp_test, cgeqrfp2) {
    EXPECT_NEAR(0.0, cgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgeqrfp_test, cgeqrfp3) {
    EXPECT_NEAR(0.0, cgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgeqrfp_test, cgeqrfp4) {
    EXPECT_NEAR(0.0, cgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class geqrfp_dcomplex_parameters{

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
      geqrfp_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrfp_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
geqrfp_dcomplex_parameters:: geqrfp_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrfp dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqrfp_float_parameters object: malloc error.";
		geqrfp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	
} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
geqrfp_dcomplex_parameters :: ~geqrfp_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrfp_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrfp_free();

}
/*  Test fixture class definition */
class zgeqrfp_test  : public  ::testing::Test {
public:
   geqrfp_dcomplex_parameters  *zgeqrfp_obj;
   void SetUp();
   void TearDown () { delete zgeqrfp_obj; }
};

void zgeqrfp_test::SetUp(){

    /* LAPACKE zgeqrfp prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeqrfp) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zgeqrfp zgeqrfp;

    zgeqrfp_obj = new geqrfp_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zgeqrfp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeqrfp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeqrfp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeqrfp_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgeqrfp = (Fptr_NL_LAPACKE_zgeqrfp)dlsym(zgeqrfp_obj->hModule, "LAPACKE_zgeqrfp");
    ASSERT_TRUE(zgeqrfp != NULL) << "failed to get the Netlib LAPACKE_zgeqrfp symbol";
    

    zgeqrfp_obj->inforef = zgeqrfp( zgeqrfp_obj->matrix_layout, zgeqrfp_obj->m,
								zgeqrfp_obj->n,zgeqrfp_obj->Aref,
								zgeqrfp_obj->lda, zgeqrfp_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zgeqrfp_obj->info = LAPACKE_zgeqrfp( zgeqrfp_obj->matrix_layout, zgeqrfp_obj->m,
										zgeqrfp_obj->n,zgeqrfp_obj->A, 
										zgeqrfp_obj->lda, zgeqrfp_obj->tau);

    if( zgeqrfp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgeqrfp is wrong\n", zgeqrfp_obj->info );
    }
    if( zgeqrfp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeqrfp is wrong\n", 
        zgeqrfp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgeqrfp_obj->diff =  computeDiff_z( zgeqrfp_obj->bufsize, 
                zgeqrfp_obj->A, zgeqrfp_obj->Aref );

}

TEST_F(zgeqrfp_test, zgeqrfp1) {
    EXPECT_NEAR(0.0, zgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgeqrfp_test, zgeqrfp2) {
    EXPECT_NEAR(0.0, zgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgeqrfp_test, zgeqrfp3) {
    EXPECT_NEAR(0.0, zgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgeqrfp_test, zgeqrfp4) {
    EXPECT_NEAR(0.0, zgeqrfp_obj->diff, LAPACKE_EIG_THRESHOLD);
}


