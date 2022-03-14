#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define gelq2_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (tau!=NULL)  free(tau);\
if (tauref!=NULL) free(tauref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class gelq2_float_parameters{

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
      gelq2_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gelq2_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gelq2_float_parameters:: gelq2_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelq2 float:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{ 	lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&tau, &tauref, (min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "gelq2_float_parameters object: malloc error.";
		gelq2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gelq2_float_parameters :: ~gelq2_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelq2_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelq2_free();

}
/*  Test fixture class definition */
class sgelq2_test  : public  ::testing::Test {
public:
   gelq2_float_parameters  *sgelq2_obj;
   void SetUp();
   void TearDown () { delete sgelq2_obj; }
};

void sgelq2_test::SetUp(){

    /* LAPACKE sgelq2 prototype */
    typedef int (*Fptr_NL_LAPACKE_sgelq2) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, float* tau);

    Fptr_NL_LAPACKE_sgelq2 sgelq2;

    sgelq2_obj = new gelq2_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgelq2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgelq2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgelq2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgelq2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgelq2 = (Fptr_NL_LAPACKE_sgelq2)dlsym(sgelq2_obj->hModule, "LAPACKE_sgelq2");
    ASSERT_TRUE(sgelq2 != NULL) << "failed to get the Netlib LAPACKE_sgelq2 symbol";
    

    sgelq2_obj->inforef = sgelq2( sgelq2_obj->matrix_layout, sgelq2_obj->m,
								sgelq2_obj->n,sgelq2_obj->Aref,
								sgelq2_obj->lda, sgelq2_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    sgelq2_obj->info = LAPACKE_sgelq2( sgelq2_obj->matrix_layout, sgelq2_obj->m,
										sgelq2_obj->n,sgelq2_obj->A, 
										sgelq2_obj->lda, sgelq2_obj->tau);

    if( sgelq2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgelq2 is wrong\n", sgelq2_obj->info );
    }
    if( sgelq2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgelq2 is wrong\n", 
        sgelq2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgelq2_obj->diff =  computeDiff_s( sgelq2_obj->bufsize, 
                sgelq2_obj->A, sgelq2_obj->Aref );

}

TEST_F(sgelq2_test, sgelq21) {
    EXPECT_NEAR(0.0, sgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelq2_test, sgelq22) {
    EXPECT_NEAR(0.0, sgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelq2_test, sgelq23) {
    EXPECT_NEAR(0.0, sgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelq2_test, sgelq24) {
    EXPECT_NEAR(0.0, sgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class gelq2_double_parameters{

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
      gelq2_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gelq2_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gelq2_double_parameters:: gelq2_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelq2 double:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) 
	{	lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&tau, &tauref, (min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "gelq2_double_parameters object: malloc error.";
		gelq2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gelq2_double_parameters :: ~gelq2_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelq2_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelq2_free();

}
/*  Test fixture class definition */
class dgelq2_test  : public  ::testing::Test {
public:
   gelq2_double_parameters  *dgelq2_obj;
   void SetUp();
   void TearDown () { delete dgelq2_obj; }
};

void dgelq2_test::SetUp(){

    /* LAPACKE dgelq2 prototype */
    typedef int (*Fptr_NL_LAPACKE_dgelq2) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, double* tau);

    Fptr_NL_LAPACKE_dgelq2 dgelq2;

    dgelq2_obj = new gelq2_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].lda);
	
    idx = Circular_Increment_Index(idx);
    dgelq2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgelq2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgelq2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgelq2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgelq2 = (Fptr_NL_LAPACKE_dgelq2)dlsym(dgelq2_obj->hModule, "LAPACKE_dgelq2");
    ASSERT_TRUE(dgelq2 != NULL) << "failed to get the Netlib LAPACKE_dgelq2 symbol";
    

    dgelq2_obj->inforef = dgelq2( dgelq2_obj->matrix_layout, dgelq2_obj->m,
								dgelq2_obj->n,dgelq2_obj->Aref,
								dgelq2_obj->lda, dgelq2_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    dgelq2_obj->info = LAPACKE_dgelq2( dgelq2_obj->matrix_layout, dgelq2_obj->m,
										dgelq2_obj->n,dgelq2_obj->A, 
										dgelq2_obj->lda, dgelq2_obj->tau);

    if( dgelq2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgelq2 is wrong\n", dgelq2_obj->info );
    }
    if( dgelq2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgelq2 is wrong\n", 
        dgelq2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgelq2_obj->diff =  computeDiff_d( dgelq2_obj->bufsize, 
                dgelq2_obj->A, dgelq2_obj->Aref );

}

TEST_F(dgelq2_test, dgelq21) {
    EXPECT_NEAR(0.0, dgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelq2_test, dgelq22) {
    EXPECT_NEAR(0.0, dgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelq2_test, dgelq23) {
    EXPECT_NEAR(0.0, dgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelq2_test, dgelq24) {
    EXPECT_NEAR(0.0, dgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class gelq2_scomplex_parameters{

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
      gelq2_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gelq2_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
gelq2_scomplex_parameters:: gelq2_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelq2 scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = m;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&tau, &tauref, (min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "gelq2_float_parameters object: malloc error.";
		gelq2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gelq2_scomplex_parameters :: ~gelq2_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelq2_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelq2_free();

}
/*  Test fixture class definition */
class cgelq2_test  : public  ::testing::Test {
public:
   gelq2_scomplex_parameters  *cgelq2_obj;
   void SetUp();
   void TearDown () { delete cgelq2_obj; }
};

void cgelq2_test::SetUp(){

    /* LAPACKE cgelq2 prototype */
    typedef int (*Fptr_NL_LAPACKE_cgelq2) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_complex_float* tau);

    Fptr_NL_LAPACKE_cgelq2 cgelq2;

    cgelq2_obj = new gelq2_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cgelq2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgelq2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgelq2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgelq2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgelq2 = (Fptr_NL_LAPACKE_cgelq2)dlsym(cgelq2_obj->hModule, "LAPACKE_cgelq2");
    ASSERT_TRUE(cgelq2 != NULL) << "failed to get the Netlib LAPACKE_cgelq2 symbol";
    

    cgelq2_obj->inforef = cgelq2( cgelq2_obj->matrix_layout, cgelq2_obj->m,
								cgelq2_obj->n,cgelq2_obj->Aref,
								cgelq2_obj->lda, cgelq2_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    cgelq2_obj->info = LAPACKE_cgelq2( cgelq2_obj->matrix_layout, cgelq2_obj->m,
										cgelq2_obj->n,cgelq2_obj->A, 
										cgelq2_obj->lda, cgelq2_obj->tau);

    if( cgelq2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgelq2 is wrong\n", cgelq2_obj->info );
    }
    if( cgelq2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgelq2 is wrong\n", 
        cgelq2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgelq2_obj->diff =  computeDiff_c( cgelq2_obj->bufsize, 
                cgelq2_obj->A, cgelq2_obj->Aref );

}

TEST_F(cgelq2_test, cgelq21) {
    EXPECT_NEAR(0.0, cgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelq2_test, cgelq22) {
    EXPECT_NEAR(0.0, cgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelq2_test, cgelq23) {
    EXPECT_NEAR(0.0, cgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelq2_test, cgelq24) {
    EXPECT_NEAR(0.0, cgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class gelq2_dcomplex_parameters{

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
      gelq2_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~gelq2_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
gelq2_dcomplex_parameters:: gelq2_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelq2 dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR) {
		lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&tau, &tauref, (min(m,n)));	
	if ((A==NULL) || (Aref==NULL) || \
		(tau==NULL) || (tauref==NULL)){
		EXPECT_FALSE( true) << "gelq2_float_parameters object: malloc error.";
		gelq2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
gelq2_dcomplex_parameters :: ~gelq2_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelq2_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelq2_free();

}
/*  Test fixture class definition */
class zgelq2_test  : public  ::testing::Test {
public:
   gelq2_dcomplex_parameters  *zgelq2_obj;
   void SetUp();
   void TearDown () { delete zgelq2_obj; }
};

void zgelq2_test::SetUp(){

    /* LAPACKE zgelq2 prototype */
    typedef int (*Fptr_NL_LAPACKE_zgelq2) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_complex_double* tau);

    Fptr_NL_LAPACKE_zgelq2 zgelq2;

    zgelq2_obj = new gelq2_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zgelq2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgelq2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgelq2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgelq2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgelq2 = (Fptr_NL_LAPACKE_zgelq2)dlsym(zgelq2_obj->hModule, "LAPACKE_zgelq2");
    ASSERT_TRUE(zgelq2 != NULL) << "failed to get the Netlib LAPACKE_zgelq2 symbol";
    

    zgelq2_obj->inforef = zgelq2( zgelq2_obj->matrix_layout, zgelq2_obj->m,
								zgelq2_obj->n,zgelq2_obj->Aref,
								zgelq2_obj->lda, zgelq2_obj->tauref);

    /* Compute libflame's Lapacke o/p  */
    zgelq2_obj->info = LAPACKE_zgelq2( zgelq2_obj->matrix_layout, zgelq2_obj->m,
										zgelq2_obj->n,zgelq2_obj->A, 
										zgelq2_obj->lda, zgelq2_obj->tau);

    if( zgelq2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgelq2 is wrong\n", zgelq2_obj->info );
    }
    if( zgelq2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgelq2 is wrong\n", 
        zgelq2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgelq2_obj->diff =  computeDiff_z( zgelq2_obj->bufsize, 
                zgelq2_obj->A, zgelq2_obj->Aref );

}

TEST_F(zgelq2_test, zgelq21) {
    EXPECT_NEAR(0.0, zgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelq2_test, zgelq22) {
    EXPECT_NEAR(0.0, zgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelq2_test, zgelq23) {
    EXPECT_NEAR(0.0, zgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelq2_test, zgelq24) {
    EXPECT_NEAR(0.0, zgelq2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


