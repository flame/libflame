#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define getf2_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (ipiv!=NULL)  free(ipiv);\
if (ipivref!=NULL) free(ipivref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class getf2_float_parameters{

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
	lapack_int* ipiv, *ipivref;
	float *Aref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      getf2_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~getf2_float_parameters ();

};

/* Constructor definition  float_common_parameters */
getf2_float_parameters:: getf2_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n getf2 float:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(float);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(float);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, (min(m,n)*sizeof(lapack_int)));	
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv==NULL) || (ipivref==NULL)){
		EXPECT_FALSE( true) << "getf2_float_parameters object: malloc error.";
		getf2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(min(m,n));i++) {
		ipiv[i] = 0;
		ipivref[i] = ipiv[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
getf2_float_parameters :: ~getf2_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" getf2_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   getf2_free();

}
/*  Test fixture class definition */
class sgetf2_test  : public  ::testing::Test {
public:
   getf2_float_parameters  *sgetf2_obj;
   void SetUp();
   void TearDown () { delete sgetf2_obj; }
};

void sgetf2_test::SetUp(){

    /* LAPACKE sgetf2 prototype */
    typedef int (*Fptr_NL_LAPACKE_sgetf2) (int matrix_layout, lapack_int m,lapack_int n, 
											float *A, lapack_int lda, lapack_int* ipiv);

    Fptr_NL_LAPACKE_sgetf2 sgetf2;

    sgetf2_obj = new getf2_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    sgetf2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgetf2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgetf2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgetf2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgetf2 = (Fptr_NL_LAPACKE_sgetf2)dlsym(sgetf2_obj->hModule, "LAPACKE_sgetf2");
    ASSERT_TRUE(sgetf2 != NULL) << "failed to get the Netlib LAPACKE_sgetf2 symbol";
    

    sgetf2_obj->inforef = sgetf2( sgetf2_obj->matrix_layout, sgetf2_obj->m,
								sgetf2_obj->n,sgetf2_obj->Aref,
								sgetf2_obj->lda, sgetf2_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    sgetf2_obj->info = LAPACKE_sgetf2( sgetf2_obj->matrix_layout, sgetf2_obj->m,
										sgetf2_obj->n,sgetf2_obj->A, 
										sgetf2_obj->lda, sgetf2_obj->ipiv);

    if( sgetf2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgetf2 is wrong\n", sgetf2_obj->info );
    }
    if( sgetf2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgetf2 is wrong\n", 
        sgetf2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgetf2_obj->diff =  computeDiff_s( sgetf2_obj->bufsize, 
                sgetf2_obj->A, sgetf2_obj->Aref );

}

TEST_F(sgetf2_test, sgetf21) {
    EXPECT_NEAR(0.0, sgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgetf2_test, sgetf22) {
    EXPECT_NEAR(0.0, sgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgetf2_test, sgetf23) {
    EXPECT_NEAR(0.0, sgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgetf2_test, sgetf24) {
    EXPECT_NEAR(0.0, sgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class getf2_double_parameters{

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
	lapack_int* ipiv, *ipivref;
	double *Aref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      getf2_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~getf2_double_parameters ();

};

/* Constructor definition  double_common_parameters */
getf2_double_parameters:: getf2_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n getf2 double:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(double);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(double);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, (min(m,n)*sizeof(lapack_int)));	
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv==NULL) || (ipivref==NULL)){
		EXPECT_FALSE( true) << "getf2_double_parameters object: malloc error.";
		getf2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(min(m,n));i++) {
		ipiv[i] = 0;
		ipivref[i] = ipiv[i];
	}

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
getf2_double_parameters :: ~getf2_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" getf2_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   getf2_free();

}
/*  Test fixture class definition */
class dgetf2_test  : public  ::testing::Test {
public:
   getf2_double_parameters  *dgetf2_obj;
   void SetUp();
   void TearDown () { delete dgetf2_obj; }
};

void dgetf2_test::SetUp(){

    /* LAPACKE dgetf2 prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetf2) (int matrix_layout, lapack_int m,lapack_int n, 
											double *A, lapack_int lda, lapack_int* ipiv);

    Fptr_NL_LAPACKE_dgetf2 dgetf2;

    dgetf2_obj = new getf2_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);
    idx = Circular_Increment_Index(idx);
    dgetf2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgetf2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgetf2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgetf2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgetf2 = (Fptr_NL_LAPACKE_dgetf2)dlsym(dgetf2_obj->hModule, "LAPACKE_dgetf2");
    ASSERT_TRUE(dgetf2 != NULL) << "failed to get the Netlib LAPACKE_dgetf2 symbol";
    

    dgetf2_obj->inforef = dgetf2( dgetf2_obj->matrix_layout, dgetf2_obj->m,
								dgetf2_obj->n,dgetf2_obj->Aref,
								dgetf2_obj->lda, dgetf2_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    dgetf2_obj->info = LAPACKE_dgetf2( dgetf2_obj->matrix_layout, dgetf2_obj->m,
										dgetf2_obj->n,dgetf2_obj->A, 
										dgetf2_obj->lda, dgetf2_obj->ipiv);

    if( dgetf2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgetf2 is wrong\n", dgetf2_obj->info );
    }
    if( dgetf2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgetf2 is wrong\n", 
        dgetf2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgetf2_obj->diff =  computeDiff_d( dgetf2_obj->bufsize, 
                dgetf2_obj->A, dgetf2_obj->Aref );

}

TEST_F(dgetf2_test, dgetf21) {
    EXPECT_NEAR(0.0, dgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetf2_test, dgetf22) {
    EXPECT_NEAR(0.0, dgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetf2_test, dgetf23) {
    EXPECT_NEAR(0.0, dgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetf2_test, dgetf24) {
    EXPECT_NEAR(0.0, dgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class getf2_scomplex_parameters{

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
	lapack_int* ipiv, *ipivref;
	lapack_complex_float *Aref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      getf2_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~getf2_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
getf2_scomplex_parameters:: getf2_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n getf2 scomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(lapack_complex_float);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(lapack_complex_float);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, (min(m,n)*sizeof(lapack_int)));	
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv==NULL) || (ipivref==NULL)){
		EXPECT_FALSE( true) << "getf2_float_parameters object: malloc error.";
		getf2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(min(m,n));i++) {
		ipiv[i] = 0;
		ipivref[i] = ipiv[i];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
getf2_scomplex_parameters :: ~getf2_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" getf2_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   getf2_free();

}
/*  Test fixture class definition */
class cgetf2_test  : public  ::testing::Test {
public:
   getf2_scomplex_parameters  *cgetf2_obj;
   void SetUp();
   void TearDown () { delete cgetf2_obj; }
};

void cgetf2_test::SetUp(){

    /* LAPACKE cgetf2 prototype */
    typedef int (*Fptr_NL_LAPACKE_cgetf2) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_float *A, lapack_int lda, lapack_int* ipiv);

    Fptr_NL_LAPACKE_cgetf2 cgetf2;

    cgetf2_obj = new getf2_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    cgetf2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgetf2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgetf2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgetf2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgetf2 = (Fptr_NL_LAPACKE_cgetf2)dlsym(cgetf2_obj->hModule, "LAPACKE_cgetf2");
    ASSERT_TRUE(cgetf2 != NULL) << "failed to get the Netlib LAPACKE_cgetf2 symbol";
    

    cgetf2_obj->inforef = cgetf2( cgetf2_obj->matrix_layout, cgetf2_obj->m,
								cgetf2_obj->n,cgetf2_obj->Aref,
								cgetf2_obj->lda, cgetf2_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    cgetf2_obj->info = LAPACKE_cgetf2( cgetf2_obj->matrix_layout, cgetf2_obj->m,
										cgetf2_obj->n,cgetf2_obj->A, 
										cgetf2_obj->lda, cgetf2_obj->ipiv);

    if( cgetf2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgetf2 is wrong\n", cgetf2_obj->info );
    }
    if( cgetf2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgetf2 is wrong\n", 
        cgetf2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgetf2_obj->diff =  computeDiff_c( cgetf2_obj->bufsize, 
                cgetf2_obj->A, cgetf2_obj->Aref );

}

TEST_F(cgetf2_test, cgetf21) {
    EXPECT_NEAR(0.0, cgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgetf2_test, cgetf22) {
    EXPECT_NEAR(0.0, cgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgetf2_test, cgetf23) {
    EXPECT_NEAR(0.0, cgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgetf2_test, cgetf24) {
    EXPECT_NEAR(0.0, cgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class getf2_dcomplex_parameters{

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
	lapack_int* ipiv, *ipivref;
	lapack_complex_double *Aref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      getf2_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int lda);
      ~getf2_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
getf2_dcomplex_parameters:: getf2_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n getf2 dcomplex:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = lda*n*sizeof(lapack_complex_double);
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = lda*m*sizeof(lapack_complex_double);
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, (min(m,n)*sizeof(lapack_int)));	
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv==NULL) || (ipivref==NULL)){
		EXPECT_FALSE( true) << "getf2_float_parameters object: malloc error.";
		getf2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	/*initialize output matrix by 0 */
	for(i=0;i<(min(m,n));i++) {
		ipiv[i] = 0;
		ipivref[i] = ipiv[i];
	}

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
getf2_dcomplex_parameters :: ~getf2_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" getf2_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   getf2_free();

}
/*  Test fixture class definition */
class zgetf2_test  : public  ::testing::Test {
public:
   getf2_dcomplex_parameters  *zgetf2_obj;
   void SetUp();
   void TearDown () { delete zgetf2_obj; }
};

void zgetf2_test::SetUp(){

    /* LAPACKE zgetf2 prototype */
    typedef int (*Fptr_NL_LAPACKE_zgetf2) (int matrix_layout, lapack_int m,lapack_int n, 
											lapack_complex_double *A, lapack_int lda, lapack_int* ipiv);

    Fptr_NL_LAPACKE_zgetf2 zgetf2;

    zgetf2_obj = new getf2_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zgetf2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgetf2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgetf2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgetf2_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgetf2 = (Fptr_NL_LAPACKE_zgetf2)dlsym(zgetf2_obj->hModule, "LAPACKE_zgetf2");
    ASSERT_TRUE(zgetf2 != NULL) << "failed to get the Netlib LAPACKE_zgetf2 symbol";
    

    zgetf2_obj->inforef = zgetf2( zgetf2_obj->matrix_layout, zgetf2_obj->m,
								zgetf2_obj->n,zgetf2_obj->Aref,
								zgetf2_obj->lda, zgetf2_obj->ipivref);

    /* Compute libflame's Lapacke o/p  */
    zgetf2_obj->info = LAPACKE_zgetf2( zgetf2_obj->matrix_layout, zgetf2_obj->m,
										zgetf2_obj->n,zgetf2_obj->A, 
										zgetf2_obj->lda, zgetf2_obj->ipiv);

    if( zgetf2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgetf2 is wrong\n", zgetf2_obj->info );
    }
    if( zgetf2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgetf2 is wrong\n", 
        zgetf2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgetf2_obj->diff =  computeDiff_z( zgetf2_obj->bufsize, 
                zgetf2_obj->A, zgetf2_obj->Aref );

}

TEST_F(zgetf2_test, zgetf21) {
    EXPECT_NEAR(0.0, zgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgetf2_test, zgetf22) {
    EXPECT_NEAR(0.0, zgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgetf2_test, zgetf23) {
    EXPECT_NEAR(0.0, zgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgetf2_test, zgetf24) {
    EXPECT_NEAR(0.0, zgetf2_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


