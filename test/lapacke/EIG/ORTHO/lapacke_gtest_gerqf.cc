#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define gerqf_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class gerqf_float_parameters{

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
      gerqf_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gerqf_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gerqf_float_parameters:: gerqf_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gerqf float:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(float);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(float);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&tau, &tauref, (fla_min(m,n)*sizeof(float)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "gerqf_float_parameters object: malloc error.";
		gerqf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(fla_min(m,n));i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gerqf_float_parameters :: ~gerqf_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gerqf_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gerqf_free();

}
/*  Test fixture class definition */
class sgerqf_test  : public  ::testing::Test {
public:
   gerqf_float_parameters  *sgerqf_obj;
   void SetUp();
   void TearDown () { delete sgerqf_obj; }
};

void sgerqf_test::SetUp(){

    /* LAPACKE sgerqf prototype */
    typedef int (*Fptr_NL_LAPACKE_sgerqf) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, float* tau);

    Fptr_NL_LAPACKE_sgerqf sgerqf;

    sgerqf_obj = new gerqf_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgerqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgerqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgerqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgerqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgerqf = (Fptr_NL_LAPACKE_sgerqf)dlsym(sgerqf_obj->hModule, "LAPACKE_sgerqf");
    ASSERT_TRUE(sgerqf != NULL) << "failed to get the Netlib LAPACKE_sgerqf symbol";
    

    sgerqf_obj->inforef = sgerqf( sgerqf_obj->matrix_layout, sgerqf_obj->m,
								sgerqf_obj->n,sgerqf_obj->Aref,
								sgerqf_obj->lda, sgerqf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    sgerqf_obj->info = LAPACKE_sgerqf( sgerqf_obj->matrix_layout, sgerqf_obj->m,
										sgerqf_obj->n,sgerqf_obj->A, 
										sgerqf_obj->lda, sgerqf_obj->tau);

    if( sgerqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgerqf is wrong\n", sgerqf_obj->info );
    }
    if( sgerqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgerqf is wrong\n", 
        sgerqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgerqf_obj->diff =  computeDiff_s( sgerqf_obj->bufsize, 
                sgerqf_obj->A, sgerqf_obj->Aref );

}

TEST_F(sgerqf_test, sgerqf1) {
    EXPECT_NEAR(0.0, sgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgerqf_test, sgerqf2) {
    EXPECT_NEAR(0.0, sgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgerqf_test, sgerqf3) {
    EXPECT_NEAR(0.0, sgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgerqf_test, sgerqf4) {
    EXPECT_NEAR(0.0, sgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class gerqf_double_parameters{

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
      gerqf_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gerqf_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gerqf_double_parameters:: gerqf_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gerqf double:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(double);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(double);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&tau, &tauref, (fla_min(m,n)*sizeof(double)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "gerqf_double_parameters object: malloc error.";
		gerqf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(fla_min(m,n));i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gerqf_double_parameters :: ~gerqf_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gerqf_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gerqf_free();

}
/*  Test fixture class definition */
class dgerqf_test  : public  ::testing::Test {
public:
   gerqf_double_parameters  *dgerqf_obj;
   void SetUp();
   void TearDown () { delete dgerqf_obj; }
};

void dgerqf_test::SetUp(){

    /* LAPACKE dgerqf prototype */
    typedef int (*Fptr_NL_LAPACKE_dgerqf) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, double* tau);

    Fptr_NL_LAPACKE_dgerqf dgerqf;

    dgerqf_obj = new gerqf_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);
    idx = Circular_Increment_Index(idx);
    dgerqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgerqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgerqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgerqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgerqf = (Fptr_NL_LAPACKE_dgerqf)dlsym(dgerqf_obj->hModule, "LAPACKE_dgerqf");
    ASSERT_TRUE(dgerqf != NULL) << "failed to get the Netlib LAPACKE_dgerqf symbol";
    

    dgerqf_obj->inforef = dgerqf( dgerqf_obj->matrix_layout, dgerqf_obj->m,
								dgerqf_obj->n,dgerqf_obj->Aref,
								dgerqf_obj->lda, dgerqf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dgerqf_obj->info = LAPACKE_dgerqf( dgerqf_obj->matrix_layout, dgerqf_obj->m,
										dgerqf_obj->n,dgerqf_obj->A, 
										dgerqf_obj->lda, dgerqf_obj->tau);

    if( dgerqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgerqf is wrong\n", dgerqf_obj->info );
    }
    if( dgerqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgerqf is wrong\n", 
        dgerqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgerqf_obj->diff =  computeDiff_d( dgerqf_obj->bufsize, 
                dgerqf_obj->A, dgerqf_obj->Aref );

}

TEST_F(dgerqf_test, dgerqf1) {
    EXPECT_NEAR(0.0, dgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgerqf_test, dgerqf2) {
    EXPECT_NEAR(0.0, dgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgerqf_test, dgerqf3) {
    EXPECT_NEAR(0.0, dgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgerqf_test, dgerqf4) {
    EXPECT_NEAR(0.0, dgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class gerqf_scomplex_parameters{

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
      gerqf_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gerqf_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
gerqf_scomplex_parameters:: gerqf_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gerqf scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(lapack_complex_float);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(lapack_complex_float);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (fla_min(m,n)*sizeof(lapack_complex_float)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "gerqf_float_parameters object: malloc error.";
		gerqf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(fla_min(m,n));i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gerqf_scomplex_parameters :: ~gerqf_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gerqf_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gerqf_free();

}
/*  Test fixture class definition */
class cgerqf_test  : public  ::testing::Test {
public:
   gerqf_scomplex_parameters  *cgerqf_obj;
   void SetUp();
   void TearDown () { delete cgerqf_obj; }
};

void cgerqf_test::SetUp(){

    /* LAPACKE cgerqf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgerqf) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cgerqf cgerqf;

    cgerqf_obj = new gerqf_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cgerqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgerqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgerqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgerqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgerqf = (Fptr_NL_LAPACKE_cgerqf)dlsym(cgerqf_obj->hModule, "LAPACKE_cgerqf");
    ASSERT_TRUE(cgerqf != NULL) << "failed to get the Netlib LAPACKE_cgerqf symbol";
    

    cgerqf_obj->inforef = cgerqf( cgerqf_obj->matrix_layout, cgerqf_obj->m,
								cgerqf_obj->n,cgerqf_obj->Aref,
								cgerqf_obj->lda, cgerqf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cgerqf_obj->info = LAPACKE_cgerqf( cgerqf_obj->matrix_layout, cgerqf_obj->m,
										cgerqf_obj->n,cgerqf_obj->A, 
										cgerqf_obj->lda, cgerqf_obj->tau);

    if( cgerqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgerqf is wrong\n", cgerqf_obj->info );
    }
    if( cgerqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgerqf is wrong\n", 
        cgerqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgerqf_obj->diff =  computeDiff_c( cgerqf_obj->bufsize, 
                cgerqf_obj->A, cgerqf_obj->Aref );

}

TEST_F(cgerqf_test, cgerqf1) {
    EXPECT_NEAR(0.0, cgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgerqf_test, cgerqf2) {
    EXPECT_NEAR(0.0, cgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgerqf_test, cgerqf3) {
    EXPECT_NEAR(0.0, cgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgerqf_test, cgerqf4) {
    EXPECT_NEAR(0.0, cgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class gerqf_dcomplex_parameters{

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
      gerqf_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gerqf_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
gerqf_dcomplex_parameters:: gerqf_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gerqf dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(lapack_complex_double);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(lapack_complex_double);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (fla_min(m,n)*sizeof(lapack_complex_double)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "gerqf_float_parameters object: malloc error.";
		gerqf_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(fla_min(m,n));i++) {
		tau[i] = 0;
		tauref[i] = tau[i];
	}

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
gerqf_dcomplex_parameters :: ~gerqf_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gerqf_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gerqf_free();

}
/*  Test fixture class definition */
class zgerqf_test  : public  ::testing::Test {
public:
   gerqf_dcomplex_parameters  *zgerqf_obj;
   void SetUp();
   void TearDown () { delete zgerqf_obj; }
};

void zgerqf_test::SetUp(){

    /* LAPACKE zgerqf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgerqf) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zgerqf zgerqf;

    zgerqf_obj = new gerqf_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zgerqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgerqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgerqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgerqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgerqf = (Fptr_NL_LAPACKE_zgerqf)dlsym(zgerqf_obj->hModule, "LAPACKE_zgerqf");
    ASSERT_TRUE(zgerqf != NULL) << "failed to get the Netlib LAPACKE_zgerqf symbol";
    

    zgerqf_obj->inforef = zgerqf( zgerqf_obj->matrix_layout, zgerqf_obj->m,
								zgerqf_obj->n,zgerqf_obj->Aref,
								zgerqf_obj->lda, zgerqf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zgerqf_obj->info = LAPACKE_zgerqf( zgerqf_obj->matrix_layout, zgerqf_obj->m,
										zgerqf_obj->n,zgerqf_obj->A, 
										zgerqf_obj->lda, zgerqf_obj->tau);

    if( zgerqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgerqf is wrong\n", zgerqf_obj->info );
    }
    if( zgerqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgerqf is wrong\n", 
        zgerqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgerqf_obj->diff =  computeDiff_z( zgerqf_obj->bufsize, 
                zgerqf_obj->A, zgerqf_obj->Aref );

}

TEST_F(zgerqf_test, zgerqf1) {
    EXPECT_NEAR(0.0, zgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgerqf_test, zgerqf2) {
    EXPECT_NEAR(0.0, zgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgerqf_test, zgerqf3) {
    EXPECT_NEAR(0.0, zgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgerqf_test, zgerqf4) {
    EXPECT_NEAR(0.0, zgerqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


