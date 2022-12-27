#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define geqr2_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class geqr2_float_parameters{

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
      geqr2_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqr2_float_parameters ();

};

/* Constructor definition  float_common_parameters */
geqr2_float_parameters:: geqr2_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqr2 float:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&tau, &tauref, (fla_min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqr2_float_parameters object: malloc error.";
		geqr2_free();
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
geqr2_float_parameters :: ~geqr2_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqr2_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqr2_free();

}
/*  Test fixture class definition */
class sgeqr2_test  : public  ::testing::Test {
public:
   geqr2_float_parameters  *sgeqr2_obj;
   void SetUp();
   void TearDown () { delete sgeqr2_obj; }
};

void sgeqr2_test::SetUp(){

    /* LAPACKE sgeqr2 prototype */
    typedef int (*Fptr_NL_LAPACKE_sgeqr2) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, float* tau);

    Fptr_NL_LAPACKE_sgeqr2 sgeqr2;

    sgeqr2_obj = new geqr2_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgeqr2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgeqr2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgeqr2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgeqr2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgeqr2 = (Fptr_NL_LAPACKE_sgeqr2)dlsym(sgeqr2_obj->hModule, "LAPACKE_sgeqr2");
    ASSERT_TRUE(sgeqr2 != NULL) << "failed to get the Netlib LAPACKE_sgeqr2 symbol";
    

    sgeqr2_obj->inforef = sgeqr2( sgeqr2_obj->matrix_layout, sgeqr2_obj->m,
								sgeqr2_obj->n,sgeqr2_obj->Aref,
								sgeqr2_obj->lda, sgeqr2_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    sgeqr2_obj->info = LAPACKE_sgeqr2( sgeqr2_obj->matrix_layout, sgeqr2_obj->m,
										sgeqr2_obj->n,sgeqr2_obj->A, 
										sgeqr2_obj->lda, sgeqr2_obj->tau);

    if( sgeqr2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgeqr2 is wrong\n", sgeqr2_obj->info );
    }
    if( sgeqr2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgeqr2 is wrong\n", 
        sgeqr2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgeqr2_obj->diff =  computeDiff_s( sgeqr2_obj->bufsize, 
                sgeqr2_obj->A, sgeqr2_obj->Aref );

}

TEST_F(sgeqr2_test, sgeqr21) {
    EXPECT_NEAR(0.0, sgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgeqr2_test, sgeqr22) {
    EXPECT_NEAR(0.0, sgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgeqr2_test, sgeqr23) {
    EXPECT_NEAR(0.0, sgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgeqr2_test, sgeqr24) {
    EXPECT_NEAR(0.0, sgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class geqr2_double_parameters{

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
      geqr2_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqr2_double_parameters ();

};

/* Constructor definition  double_common_parameters */
geqr2_double_parameters:: geqr2_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqr2 double:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&tau, &tauref, (fla_min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqr2_double_parameters object: malloc error.";
		geqr2_free();
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
geqr2_double_parameters :: ~geqr2_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqr2_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqr2_free();

}
/*  Test fixture class definition */
class dgeqr2_test  : public  ::testing::Test {
public:
   geqr2_double_parameters  *dgeqr2_obj;
   void SetUp();
   void TearDown () { delete dgeqr2_obj; }
};

void dgeqr2_test::SetUp(){

    /* LAPACKE dgeqr2 prototype */
    typedef int (*Fptr_NL_LAPACKE_dgeqr2) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, double* tau);

    Fptr_NL_LAPACKE_dgeqr2 dgeqr2;

    dgeqr2_obj = new geqr2_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);
    idx = Circular_Increment_Index(idx);
    dgeqr2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgeqr2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgeqr2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgeqr2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgeqr2 = (Fptr_NL_LAPACKE_dgeqr2)dlsym(dgeqr2_obj->hModule, "LAPACKE_dgeqr2");
    ASSERT_TRUE(dgeqr2 != NULL) << "failed to get the Netlib LAPACKE_dgeqr2 symbol";
    

    dgeqr2_obj->inforef = dgeqr2( dgeqr2_obj->matrix_layout, dgeqr2_obj->m,
								dgeqr2_obj->n,dgeqr2_obj->Aref,
								dgeqr2_obj->lda, dgeqr2_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dgeqr2_obj->info = LAPACKE_dgeqr2( dgeqr2_obj->matrix_layout, dgeqr2_obj->m,
										dgeqr2_obj->n,dgeqr2_obj->A, 
										dgeqr2_obj->lda, dgeqr2_obj->tau);

    if( dgeqr2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgeqr2 is wrong\n", dgeqr2_obj->info );
    }
    if( dgeqr2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgeqr2 is wrong\n", 
        dgeqr2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgeqr2_obj->diff =  computeDiff_d( dgeqr2_obj->bufsize, 
                dgeqr2_obj->A, dgeqr2_obj->Aref );

}

TEST_F(dgeqr2_test, dgeqr21) {
    EXPECT_NEAR(0.0, dgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgeqr2_test, dgeqr22) {
    EXPECT_NEAR(0.0, dgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgeqr2_test, dgeqr23) {
    EXPECT_NEAR(0.0, dgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgeqr2_test, dgeqr24) {
    EXPECT_NEAR(0.0, dgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class geqr2_scomplex_parameters{

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
      geqr2_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqr2_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
geqr2_scomplex_parameters:: geqr2_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqr2 scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (fla_min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqr2_float_parameters object: malloc error.";
		geqr2_free();
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
geqr2_scomplex_parameters :: ~geqr2_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqr2_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqr2_free();

}
/*  Test fixture class definition */
class cgeqr2_test  : public  ::testing::Test {
public:
   geqr2_scomplex_parameters  *cgeqr2_obj;
   void SetUp();
   void TearDown () { delete cgeqr2_obj; }
};

void cgeqr2_test::SetUp(){

    /* LAPACKE cgeqr2 prototype */
    typedef int (*Fptr_NL_LAPACKE_cgeqr2) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cgeqr2 cgeqr2;

    cgeqr2_obj = new geqr2_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cgeqr2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgeqr2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgeqr2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgeqr2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgeqr2 = (Fptr_NL_LAPACKE_cgeqr2)dlsym(cgeqr2_obj->hModule, "LAPACKE_cgeqr2");
    ASSERT_TRUE(cgeqr2 != NULL) << "failed to get the Netlib LAPACKE_cgeqr2 symbol";
    

    cgeqr2_obj->inforef = cgeqr2( cgeqr2_obj->matrix_layout, cgeqr2_obj->m,
								cgeqr2_obj->n,cgeqr2_obj->Aref,
								cgeqr2_obj->lda, cgeqr2_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cgeqr2_obj->info = LAPACKE_cgeqr2( cgeqr2_obj->matrix_layout, cgeqr2_obj->m,
										cgeqr2_obj->n,cgeqr2_obj->A, 
										cgeqr2_obj->lda, cgeqr2_obj->tau);

    if( cgeqr2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgeqr2 is wrong\n", cgeqr2_obj->info );
    }
    if( cgeqr2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgeqr2 is wrong\n", 
        cgeqr2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgeqr2_obj->diff =  computeDiff_c( cgeqr2_obj->bufsize, 
                cgeqr2_obj->A, cgeqr2_obj->Aref );

}

TEST_F(cgeqr2_test, cgeqr21) {
    EXPECT_NEAR(0.0, cgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgeqr2_test, cgeqr22) {
    EXPECT_NEAR(0.0, cgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgeqr2_test, cgeqr23) {
    EXPECT_NEAR(0.0, cgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgeqr2_test, cgeqr24) {
    EXPECT_NEAR(0.0, cgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class geqr2_dcomplex_parameters{

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
      geqr2_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~geqr2_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
geqr2_dcomplex_parameters:: geqr2_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n geqr2 dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (fla_min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "geqr2_float_parameters object: malloc error.";
		geqr2_free();
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
geqr2_dcomplex_parameters :: ~geqr2_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" geqr2_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   geqr2_free();

}
/*  Test fixture class definition */
class zgeqr2_test  : public  ::testing::Test {
public:
   geqr2_dcomplex_parameters  *zgeqr2_obj;
   void SetUp();
   void TearDown () { delete zgeqr2_obj; }
};

void zgeqr2_test::SetUp(){

    /* LAPACKE zgeqr2 prototype */
    typedef int (*Fptr_NL_LAPACKE_zgeqr2) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zgeqr2 zgeqr2;

    zgeqr2_obj = new geqr2_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zgeqr2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgeqr2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgeqr2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgeqr2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgeqr2 = (Fptr_NL_LAPACKE_zgeqr2)dlsym(zgeqr2_obj->hModule, "LAPACKE_zgeqr2");
    ASSERT_TRUE(zgeqr2 != NULL) << "failed to get the Netlib LAPACKE_zgeqr2 symbol";
    

    zgeqr2_obj->inforef = zgeqr2( zgeqr2_obj->matrix_layout, zgeqr2_obj->m,
								zgeqr2_obj->n,zgeqr2_obj->Aref,
								zgeqr2_obj->lda, zgeqr2_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zgeqr2_obj->info = LAPACKE_zgeqr2( zgeqr2_obj->matrix_layout, zgeqr2_obj->m,
										zgeqr2_obj->n,zgeqr2_obj->A, 
										zgeqr2_obj->lda, zgeqr2_obj->tau);

    if( zgeqr2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgeqr2 is wrong\n", zgeqr2_obj->info );
    }
    if( zgeqr2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgeqr2 is wrong\n", 
        zgeqr2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgeqr2_obj->diff =  computeDiff_z( zgeqr2_obj->bufsize, 
                zgeqr2_obj->A, zgeqr2_obj->Aref );

}

TEST_F(zgeqr2_test, zgeqr21) {
    EXPECT_NEAR(0.0, zgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgeqr2_test, zgeqr22) {
    EXPECT_NEAR(0.0, zgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgeqr2_test, zgeqr23) {
    EXPECT_NEAR(0.0, zgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgeqr2_test, zgeqr24) {
    EXPECT_NEAR(0.0, zgeqr2_obj->diff, LAPACKE_EIG_THRESHOLD);
}


