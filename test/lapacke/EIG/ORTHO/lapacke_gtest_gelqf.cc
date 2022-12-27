#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define gelqf_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class gelqf_float_parameters{

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
      gelqf_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gelqf_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gelqf_float_parameters:: gelqf_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelqf float:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "gelqf_float_parameters object: malloc error.";
		gelqf_free();
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
gelqf_float_parameters :: ~gelqf_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelqf_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelqf_free();

}
/*  Test fixture class definition */
class sgelqf_test  : public  ::testing::Test {
public:
   gelqf_float_parameters  *sgelqf_obj;
   void SetUp();
   void TearDown () { delete sgelqf_obj; }
};

void sgelqf_test::SetUp(){

    /* LAPACKE sgelqf prototype */
    typedef int (*Fptr_NL_LAPACKE_sgelqf) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, float* tau);

    Fptr_NL_LAPACKE_sgelqf sgelqf;

    sgelqf_obj = new gelqf_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgelqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgelqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgelqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgelqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgelqf = (Fptr_NL_LAPACKE_sgelqf)dlsym(sgelqf_obj->hModule, "LAPACKE_sgelqf");
    ASSERT_TRUE(sgelqf != NULL) << "failed to get the Netlib LAPACKE_sgelqf symbol";
    

    sgelqf_obj->inforef = sgelqf( sgelqf_obj->matrix_layout, sgelqf_obj->m,
								sgelqf_obj->n,sgelqf_obj->Aref,
								sgelqf_obj->lda, sgelqf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    sgelqf_obj->info = LAPACKE_sgelqf( sgelqf_obj->matrix_layout, sgelqf_obj->m,
										sgelqf_obj->n,sgelqf_obj->A, 
										sgelqf_obj->lda, sgelqf_obj->tau);

    if( sgelqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgelqf is wrong\n", sgelqf_obj->info );
    }
    if( sgelqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgelqf is wrong\n", 
        sgelqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgelqf_obj->diff =  computeDiff_s( sgelqf_obj->bufsize, 
                sgelqf_obj->A, sgelqf_obj->Aref );

}

TEST_F(sgelqf_test, sgelqf1) {
    EXPECT_NEAR(0.0, sgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelqf_test, sgelqf2) {
    EXPECT_NEAR(0.0, sgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelqf_test, sgelqf3) {
    EXPECT_NEAR(0.0, sgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelqf_test, sgelqf4) {
    EXPECT_NEAR(0.0, sgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class gelqf_double_parameters{

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
      gelqf_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gelqf_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gelqf_double_parameters:: gelqf_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelqf double:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "gelqf_double_parameters object: malloc error.";
		gelqf_free();
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
gelqf_double_parameters :: ~gelqf_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelqf_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelqf_free();

}
/*  Test fixture class definition */
class dgelqf_test  : public  ::testing::Test {
public:
   gelqf_double_parameters  *dgelqf_obj;
   void SetUp();
   void TearDown () { delete dgelqf_obj; }
};

void dgelqf_test::SetUp(){

    /* LAPACKE dgelqf prototype */
    typedef int (*Fptr_NL_LAPACKE_dgelqf) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, double* tau);

    Fptr_NL_LAPACKE_dgelqf dgelqf;

    dgelqf_obj = new gelqf_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);
    idx = Circular_Increment_Index(idx);
    dgelqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgelqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgelqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgelqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgelqf = (Fptr_NL_LAPACKE_dgelqf)dlsym(dgelqf_obj->hModule, "LAPACKE_dgelqf");
    ASSERT_TRUE(dgelqf != NULL) << "failed to get the Netlib LAPACKE_dgelqf symbol";
    

    dgelqf_obj->inforef = dgelqf( dgelqf_obj->matrix_layout, dgelqf_obj->m,
								dgelqf_obj->n,dgelqf_obj->Aref,
								dgelqf_obj->lda, dgelqf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dgelqf_obj->info = LAPACKE_dgelqf( dgelqf_obj->matrix_layout, dgelqf_obj->m,
										dgelqf_obj->n,dgelqf_obj->A, 
										dgelqf_obj->lda, dgelqf_obj->tau);

    if( dgelqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgelqf is wrong\n", dgelqf_obj->info );
    }
    if( dgelqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgelqf is wrong\n", 
        dgelqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgelqf_obj->diff =  computeDiff_d( dgelqf_obj->bufsize, 
                dgelqf_obj->A, dgelqf_obj->Aref );

}

TEST_F(dgelqf_test, dgelqf1) {
    EXPECT_NEAR(0.0, dgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelqf_test, dgelqf2) {
    EXPECT_NEAR(0.0, dgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelqf_test, dgelqf3) {
    EXPECT_NEAR(0.0, dgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelqf_test, dgelqf4) {
    EXPECT_NEAR(0.0, dgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class gelqf_scomplex_parameters{

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
      gelqf_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gelqf_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
gelqf_scomplex_parameters:: gelqf_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelqf scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "gelqf_float_parameters object: malloc error.";
		gelqf_free();
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
gelqf_scomplex_parameters :: ~gelqf_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelqf_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelqf_free();

}
/*  Test fixture class definition */
class cgelqf_test  : public  ::testing::Test {
public:
   gelqf_scomplex_parameters  *cgelqf_obj;
   void SetUp();
   void TearDown () { delete cgelqf_obj; }
};

void cgelqf_test::SetUp(){

    /* LAPACKE cgelqf prototype */
    typedef int (*Fptr_NL_LAPACKE_cgelqf) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cgelqf cgelqf;

    cgelqf_obj = new gelqf_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cgelqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgelqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgelqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgelqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgelqf = (Fptr_NL_LAPACKE_cgelqf)dlsym(cgelqf_obj->hModule, "LAPACKE_cgelqf");
    ASSERT_TRUE(cgelqf != NULL) << "failed to get the Netlib LAPACKE_cgelqf symbol";
    

    cgelqf_obj->inforef = cgelqf( cgelqf_obj->matrix_layout, cgelqf_obj->m,
								cgelqf_obj->n,cgelqf_obj->Aref,
								cgelqf_obj->lda, cgelqf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cgelqf_obj->info = LAPACKE_cgelqf( cgelqf_obj->matrix_layout, cgelqf_obj->m,
										cgelqf_obj->n,cgelqf_obj->A, 
										cgelqf_obj->lda, cgelqf_obj->tau);

    if( cgelqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgelqf is wrong\n", cgelqf_obj->info );
    }
    if( cgelqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgelqf is wrong\n", 
        cgelqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgelqf_obj->diff =  computeDiff_c( cgelqf_obj->bufsize, 
                cgelqf_obj->A, cgelqf_obj->Aref );

}

TEST_F(cgelqf_test, cgelqf1) {
    EXPECT_NEAR(0.0, cgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelqf_test, cgelqf2) {
    EXPECT_NEAR(0.0, cgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelqf_test, cgelqf3) {
    EXPECT_NEAR(0.0, cgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelqf_test, cgelqf4) {
    EXPECT_NEAR(0.0, cgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class gelqf_dcomplex_parameters{

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
      gelqf_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gelqf_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
gelqf_dcomplex_parameters:: gelqf_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelqf dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
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
		EXPECT_FALSE( true) << "gelqf_float_parameters object: malloc error.";
		gelqf_free();
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
gelqf_dcomplex_parameters :: ~gelqf_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelqf_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelqf_free();

}
/*  Test fixture class definition */
class zgelqf_test  : public  ::testing::Test {
public:
   gelqf_dcomplex_parameters  *zgelqf_obj;
   void SetUp();
   void TearDown () { delete zgelqf_obj; }
};

void zgelqf_test::SetUp(){

    /* LAPACKE zgelqf prototype */
    typedef int (*Fptr_NL_LAPACKE_zgelqf) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zgelqf zgelqf;

    zgelqf_obj = new gelqf_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zgelqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgelqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgelqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgelqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgelqf = (Fptr_NL_LAPACKE_zgelqf)dlsym(zgelqf_obj->hModule, "LAPACKE_zgelqf");
    ASSERT_TRUE(zgelqf != NULL) << "failed to get the Netlib LAPACKE_zgelqf symbol";
    

    zgelqf_obj->inforef = zgelqf( zgelqf_obj->matrix_layout, zgelqf_obj->m,
								zgelqf_obj->n,zgelqf_obj->Aref,
								zgelqf_obj->lda, zgelqf_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zgelqf_obj->info = LAPACKE_zgelqf( zgelqf_obj->matrix_layout, zgelqf_obj->m,
										zgelqf_obj->n,zgelqf_obj->A, 
										zgelqf_obj->lda, zgelqf_obj->tau);

    if( zgelqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgelqf is wrong\n", zgelqf_obj->info );
    }
    if( zgelqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgelqf is wrong\n", 
        zgelqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgelqf_obj->diff =  computeDiff_z( zgelqf_obj->bufsize, 
                zgelqf_obj->A, zgelqf_obj->Aref );

}

TEST_F(zgelqf_test, zgelqf1) {
    EXPECT_NEAR(0.0, zgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelqf_test, zgelqf2) {
    EXPECT_NEAR(0.0, zgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelqf_test, zgelqf3) {
    EXPECT_NEAR(0.0, zgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelqf_test, zgelqf4) {
    EXPECT_NEAR(0.0, zgelqf_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


