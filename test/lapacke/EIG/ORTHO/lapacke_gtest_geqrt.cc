#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define geqrt_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class geqrt_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt;
	float* A;
	lapack_int lda;
	/*Output Parameter*/
	float* tau;
	float *Aref, *tauref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqrt_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nb, lapack_int lda, lapack_int ldt);
      ~geqrt_float_parameters ();

};

/* Constructor definition  float_common_parameters */
geqrt_float_parameters:: geqrt_float_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nb_i, lapack_int lda_i, lapack_int ldt_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	nb = nb_i;
	ldt = ldt_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrt float:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "geqrt_float_parameters object: malloc error.";
		geqrt_free();
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
geqrt_float_parameters :: ~geqrt_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrt_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrt_free();

}
/*  Test fixture class definition */
class sgeqrt_test  : public  ::testing::Test {
public:
   geqrt_float_parameters  *sgeqrt_obj;
   void SetUp();
   void TearDown () { delete sgeqrt_obj; }
};

void sgeqrt_test::SetUp(){

    /* LAPACKE sgeqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeqrt) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, float* tau);

    Fptr_NL_LAPACKE_sgeqrt sgeqrt;

    sgeqrt_obj = new geqrt_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgeqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgeqrt = (Fptr_NL_LAPACKE_sgeqrt)dlsym(sgeqrt_obj->hModule, "LAPACKE_sgeqrt");
    ASSERT_TRUE(sgeqrt != NULL) << "failed to get the Netlib LAPACKE_sgeqrt symbol";
    

    sgeqrt_obj->inforef = sgeqrt( sgeqrt_obj->matrix_layout, sgeqrt_obj->m,
								sgeqrt_obj->n,sgeqrt_obj->Aref,
								sgeqrt_obj->lda, sgeqrt_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    sgeqrt_obj->info = LAPACKE_sgeqrt( sgeqrt_obj->matrix_layout, sgeqrt_obj->m,
										sgeqrt_obj->n,sgeqrt_obj->A, 
										sgeqrt_obj->lda, sgeqrt_obj->tau);

    if( sgeqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgeqrt is wrong\n", sgeqrt_obj->info );
    }
    if( sgeqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeqrt is wrong\n", 
        sgeqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgeqrt_obj->diff =  computeDiff_s( sgeqrt_obj->bufsize, 
                sgeqrt_obj->A, sgeqrt_obj->Aref );

}

TEST_F(sgeqrt_test, sgeqrt1) {
    EXPECT_NEAR(0.0, sgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqrt_test, sgeqrt2) {
    EXPECT_NEAR(0.0, sgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqrt_test, sgeqrt3) {
    EXPECT_NEAR(0.0, sgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgeqrt_test, sgeqrt4) {
    EXPECT_NEAR(0.0, sgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class geqrt_double_parameters{

   public:
	int bufsize;
	double diff;
	void *hModule, *dModule;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt
	double* A;
	lapack_int lda;
	/*Output Parameter*/
	double* tau;
	double *Aref, *tauref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqrt_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrt_double_parameters ();

};

/* Constructor definition  double_common_parameters */
geqrt_double_parameters:: geqrt_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrt double:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "geqrt_double_parameters object: malloc error.";
		geqrt_free();
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
geqrt_double_parameters :: ~geqrt_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrt_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrt_free();

}
/*  Test fixture class definition */
class dgeqrt_test  : public  ::testing::Test {
public:
   geqrt_double_parameters  *dgeqrt_obj;
   void SetUp();
   void TearDown () { delete dgeqrt_obj; }
};

void dgeqrt_test::SetUp(){

    /* LAPACKE dgeqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeqrt) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, double* tau);

    Fptr_NL_LAPACKE_dgeqrt dgeqrt;

    dgeqrt_obj = new geqrt_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);
    idx = Circular_Increment_Index(idx);
    dgeqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgeqrt = (Fptr_NL_LAPACKE_dgeqrt)dlsym(dgeqrt_obj->hModule, "LAPACKE_dgeqrt");
    ASSERT_TRUE(dgeqrt != NULL) << "failed to get the Netlib LAPACKE_dgeqrt symbol";
    

    dgeqrt_obj->inforef = dgeqrt( dgeqrt_obj->matrix_layout, dgeqrt_obj->m,
								dgeqrt_obj->n,dgeqrt_obj->Aref,
								dgeqrt_obj->lda, dgeqrt_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dgeqrt_obj->info = LAPACKE_dgeqrt( dgeqrt_obj->matrix_layout, dgeqrt_obj->m,
										dgeqrt_obj->n,dgeqrt_obj->A, 
										dgeqrt_obj->lda, dgeqrt_obj->tau);

    if( dgeqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgeqrt is wrong\n", dgeqrt_obj->info );
    }
    if( dgeqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeqrt is wrong\n", 
        dgeqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgeqrt_obj->diff =  computeDiff_d( dgeqrt_obj->bufsize, 
                dgeqrt_obj->A, dgeqrt_obj->Aref );

}

TEST_F(dgeqrt_test, dgeqrt1) {
    EXPECT_NEAR(0.0, dgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqrt_test, dgeqrt2) {
    EXPECT_NEAR(0.0, dgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqrt_test, dgeqrt3) {
    EXPECT_NEAR(0.0, dgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgeqrt_test, dgeqrt4) {
    EXPECT_NEAR(0.0, dgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class geqrt_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt;
	lapack_complex_float* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float* tau;
	lapack_complex_float *Aref, *tauref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqrt_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrt_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
geqrt_scomplex_parameters:: geqrt_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrt scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "geqrt_float_parameters object: malloc error.";
		geqrt_free();
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
geqrt_scomplex_parameters :: ~geqrt_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrt_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrt_free();

}
/*  Test fixture class definition */
class cgeqrt_test  : public  ::testing::Test {
public:
   geqrt_scomplex_parameters  *cgeqrt_obj;
   void SetUp();
   void TearDown () { delete cgeqrt_obj; }
};

void cgeqrt_test::SetUp(){

    /* LAPACKE cgeqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeqrt) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cgeqrt cgeqrt;

    cgeqrt_obj = new geqrt_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cgeqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgeqrt = (Fptr_NL_LAPACKE_cgeqrt)dlsym(cgeqrt_obj->hModule, "LAPACKE_cgeqrt");
    ASSERT_TRUE(cgeqrt != NULL) << "failed to get the Netlib LAPACKE_cgeqrt symbol";
    

    cgeqrt_obj->inforef = cgeqrt( cgeqrt_obj->matrix_layout, cgeqrt_obj->m,
								cgeqrt_obj->n,cgeqrt_obj->Aref,
								cgeqrt_obj->lda, cgeqrt_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cgeqrt_obj->info = LAPACKE_cgeqrt( cgeqrt_obj->matrix_layout, cgeqrt_obj->m,
										cgeqrt_obj->n,cgeqrt_obj->A, 
										cgeqrt_obj->lda, cgeqrt_obj->tau);

    if( cgeqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgeqrt is wrong\n", cgeqrt_obj->info );
    }
    if( cgeqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeqrt is wrong\n", 
        cgeqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgeqrt_obj->diff =  computeDiff_c( cgeqrt_obj->bufsize, 
                cgeqrt_obj->A, cgeqrt_obj->Aref );

}

TEST_F(cgeqrt_test, cgeqrt1) {
    EXPECT_NEAR(0.0, cgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqrt_test, cgeqrt2) {
    EXPECT_NEAR(0.0, cgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqrt_test, cgeqrt3) {
    EXPECT_NEAR(0.0, cgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgeqrt_test, cgeqrt4) {
    EXPECT_NEAR(0.0, cgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class geqrt_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt;
	lapack_complex_double* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double* tau;
	lapack_complex_double *Aref, *tauref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      geqrt_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqrt_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
geqrt_dcomplex_parameters:: geqrt_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqrt dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "geqrt_float_parameters object: malloc error.";
		geqrt_free();
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
geqrt_dcomplex_parameters :: ~geqrt_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqrt_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqrt_free();

}
/*  Test fixture class definition */
class zgeqrt_test  : public  ::testing::Test {
public:
   geqrt_dcomplex_parameters  *zgeqrt_obj;
   void SetUp();
   void TearDown () { delete zgeqrt_obj; }
};

void zgeqrt_test::SetUp(){

    /* LAPACKE zgeqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeqrt) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zgeqrt zgeqrt;

    zgeqrt_obj = new geqrt_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zgeqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgeqrt = (Fptr_NL_LAPACKE_zgeqrt)dlsym(zgeqrt_obj->hModule, "LAPACKE_zgeqrt");
    ASSERT_TRUE(zgeqrt != NULL) << "failed to get the Netlib LAPACKE_zgeqrt symbol";
    

    zgeqrt_obj->inforef = zgeqrt( zgeqrt_obj->matrix_layout, zgeqrt_obj->m,
								zgeqrt_obj->n,zgeqrt_obj->Aref,
								zgeqrt_obj->lda, zgeqrt_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zgeqrt_obj->info = LAPACKE_zgeqrt( zgeqrt_obj->matrix_layout, zgeqrt_obj->m,
										zgeqrt_obj->n,zgeqrt_obj->A, 
										zgeqrt_obj->lda, zgeqrt_obj->tau);

    if( zgeqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgeqrt is wrong\n", zgeqrt_obj->info );
    }
    if( zgeqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeqrt is wrong\n", 
        zgeqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgeqrt_obj->diff =  computeDiff_z( zgeqrt_obj->bufsize, 
                zgeqrt_obj->A, zgeqrt_obj->Aref );

}

TEST_F(zgeqrt_test, zgeqrt1) {
    EXPECT_NEAR(0.0, zgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqrt_test, zgeqrt2) {
    EXPECT_NEAR(0.0, zgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqrt_test, zgeqrt3) {
    EXPECT_NEAR(0.0, zgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgeqrt_test, zgeqrt4) {
    EXPECT_NEAR(0.0, zgeqrt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


