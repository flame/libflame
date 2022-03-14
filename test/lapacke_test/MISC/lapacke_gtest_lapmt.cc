#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define lapmt_free() \
if (x!=NULL)    free(x); \
if (xref!=NULL) free(xref);\
if (k!=NULL)  free(k);\
if (kref!=NULL) free(kref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class lapmt_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	float* x;
	lapack_int ldx;
	lapack_logical forwrd;
	/*Output Parameter*/
	lapack_int* k;
	float* xref; 
	lapack_int* kref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      lapmt_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int ldx, lapack_logical forwrd);
      ~lapmt_float_parameters ();

};

/* Constructor definition  float_common_parameters */
lapmt_float_parameters:: lapmt_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int ldx_i, lapack_logical forwrd_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldx = ldx_i;
	forwrd = forwrd_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lapmt float:  m: %d, n: %d ldx: %d, forwrd:%d: \n",  m, n, ldx, forwrd);
	#endif
	
	bufsize = ldx *n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&k, &kref, n);	
	if ((x==NULL) || (xref==NULL) || (k==NULL) || (kref==NULL)){
		EXPECT_FALSE( true) << "lapmt_float_parameters object: malloc error.";
		lapmt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( x, xref, bufsize);
	
	lapacke_gtest_init_int_buffer_pair_with_constant(k, kref, n, forwrd);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lapmt_float_parameters :: ~lapmt_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lapmt_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lapmt_free();

}
/*  Test fixture class definition */
class slapmt_test  : public  ::testing::Test {
public:
   lapmt_float_parameters  *slapmt_obj;
   void SetUp();
   void TearDown () { delete slapmt_obj; }
};

void slapmt_test::SetUp(){

    /* LAPACKE slapmt prototype */
    typedef int (*Fptr_NL_LAPACKE_slapmt) (int matrix_layout, lapack_logical forwrd, lapack_int m,lapack_int n, 
											float* x, lapack_int ldx, lapack_int* k);

    Fptr_NL_LAPACKE_slapmt slapmt;

    slapmt_obj = new lapmt_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda,
						   eig_paramslist[idx].isgn);

    idx = Circular_Increment_Index(idx);
	

    slapmt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slapmt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slapmt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slapmt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slapmt = (Fptr_NL_LAPACKE_slapmt)dlsym(slapmt_obj->hModule, "LAPACKE_slapmt");
    ASSERT_TRUE(slapmt != NULL) << "failed to get the Netlib LAPACKE_slapmt symbol";
    

    slapmt_obj->inforef = slapmt( slapmt_obj->matrix_layout, slapmt_obj->forwrd, slapmt_obj->m, \
								slapmt_obj->n,slapmt_obj->xref, slapmt_obj->ldx, slapmt_obj->kref);

    /* Compute libflame's Lapacke o/p  */
    slapmt_obj->info = LAPACKE_slapmt( slapmt_obj->matrix_layout, slapmt_obj->forwrd, slapmt_obj->m,
										slapmt_obj->n,slapmt_obj->x, slapmt_obj->ldx, slapmt_obj->k); 

    if( slapmt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slapmt is wrong\n", slapmt_obj->info );
    }
    if( slapmt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slapmt is wrong\n", 
        slapmt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slapmt_obj->diff =  computeDiff_s( slapmt_obj->bufsize, slapmt_obj->x, slapmt_obj->xref );

}

TEST_F(slapmt_test, slapmt1) {
    EXPECT_NEAR(0.0, slapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slapmt_test, slapmt2) {
    EXPECT_NEAR(0.0, slapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slapmt_test, slapmt3) {
    EXPECT_NEAR(0.0, slapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slapmt_test, slapmt4) {
    EXPECT_NEAR(0.0, slapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class lapmt_double_parameters{

   public:
	int bufsize;
	double diff;
	void *hModule, *dModule;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	double* x;
	lapack_int ldx;
	lapack_logical forwrd;
	/*Output Parameter*/
	lapack_int* k;
	double *xref;
	lapack_int* kref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      lapmt_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int ldx, lapack_logical forwrd);
      ~lapmt_double_parameters ();

};

/* Constructor definition  double_common_parameters */
lapmt_double_parameters:: lapmt_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int ldx_i, lapack_logical forwrd_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldx = ldx_i;
	forwrd = forwrd_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lapmt double:  m: %d, n: %d ldx: %d \n",  m, n, ldx);
	#endif
	
	bufsize = ldx*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&k, &kref, n);	
	if ((x==NULL) || (xref==NULL) || \
		(k==NULL) || (kref==NULL)){
		EXPECT_FALSE( true) << "lapmt_double_parameters object: malloc error.";
		lapmt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( x, xref, bufsize);
	lapacke_gtest_init_int_buffer_pair_with_constant(k, kref, n, forwrd);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lapmt_double_parameters :: ~lapmt_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lapmt_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lapmt_free();

}
/*  Test fixture class definition */
class dlapmt_test  : public  ::testing::Test {
public:
   lapmt_double_parameters  *dlapmt_obj;
   void SetUp();
   void TearDown () { delete dlapmt_obj; }
};

void dlapmt_test::SetUp(){

    /* LAPACKE dlapmt prototype */
    typedef int (*Fptr_NL_LAPACKE_dlapmt) (int matrix_layout, lapack_logical forwrd, lapack_int m,lapack_int n, 
											double *x, lapack_int ldx, lapack_int* k);

    Fptr_NL_LAPACKE_dlapmt dlapmt;

    dlapmt_obj = new lapmt_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda, 
						   eig_paramslist[idx].isgn);

    idx = Circular_Increment_Index(idx);
    dlapmt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlapmt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlapmt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlapmt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlapmt = (Fptr_NL_LAPACKE_dlapmt)dlsym(dlapmt_obj->hModule, "LAPACKE_dlapmt");
    ASSERT_TRUE(dlapmt != NULL) << "failed to get the Netlib LAPACKE_dlapmt symbol";
    

    dlapmt_obj->inforef = dlapmt( dlapmt_obj->matrix_layout, dlapmt_obj->forwrd, dlapmt_obj->m,
								dlapmt_obj->n,dlapmt_obj->xref,
								dlapmt_obj->ldx, dlapmt_obj->kref);

    /* Compute libflame's Lapacke o/p  */
    dlapmt_obj->info = LAPACKE_dlapmt( dlapmt_obj->matrix_layout, dlapmt_obj->forwrd, dlapmt_obj->m,
										dlapmt_obj->n,dlapmt_obj->x, 
										dlapmt_obj->ldx, dlapmt_obj->k);

    if( dlapmt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlapmt is wrong\n", dlapmt_obj->info );
    }
    if( dlapmt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlapmt is wrong\n", 
        dlapmt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlapmt_obj->diff =  computeDiff_d( dlapmt_obj->bufsize, 
                dlapmt_obj->x, dlapmt_obj->xref );

}

TEST_F(dlapmt_test, dlapmt1) {
    EXPECT_NEAR(0.0, dlapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlapmt_test, dlapmt2) {
    EXPECT_NEAR(0.0, dlapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlapmt_test, dlapmt3) {
    EXPECT_NEAR(0.0, dlapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlapmt_test, dlapmt4) {
    EXPECT_NEAR(0.0, dlapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class lapmt_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_complex_float* x;
	lapack_int ldx;
	lapack_logical forwrd;
	/*Output Parameter*/
	lapack_int* k;
	lapack_complex_float *xref; 
	lapack_int *kref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      lapmt_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int ldx, lapack_logical forwrd);
      ~lapmt_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lapmt_scomplex_parameters:: lapmt_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int ldx_i, lapack_logical forwrd_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldx = ldx_i;
	forwrd = forwrd_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lapmt scomplex:  m: %d, n: %d ldx: %d \n",  m, n, ldx);
	#endif	
	
	bufsize = ldx*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&k, &kref, n);	
	if ((x==NULL) || (xref==NULL) || \
		(k==NULL) || (kref==NULL)){
		EXPECT_FALSE( true) << "lapmt_float_parameters object: malloc error.";
		lapmt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( x, xref, bufsize);
	
	lapacke_gtest_init_int_buffer_pair_with_constant(k, kref, n, forwrd);


} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lapmt_scomplex_parameters :: ~lapmt_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lapmt_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lapmt_free();

}
/*  Test fixture class definition */
class clapmt_test  : public  ::testing::Test {
public:
   lapmt_scomplex_parameters  *clapmt_obj;
   void SetUp();
   void TearDown () { delete clapmt_obj; }
};

void clapmt_test::SetUp(){

    /* LAPACKE clapmt prototype */
    typedef int (*Fptr_NL_LAPACKE_clapmt) (int matrix_layout, lapack_logical forwrd, lapack_int m,lapack_int n, 
											lapack_complex_float *x, lapack_int ldx, lapack_int* k);

    Fptr_NL_LAPACKE_clapmt clapmt;

    clapmt_obj = new lapmt_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda, 
						   eig_paramslist[idx].isgn);

    idx = Circular_Increment_Index(idx);

    clapmt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clapmt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clapmt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clapmt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clapmt = (Fptr_NL_LAPACKE_clapmt)dlsym(clapmt_obj->hModule, "LAPACKE_clapmt");
    ASSERT_TRUE(clapmt != NULL) << "failed to get the Netlib LAPACKE_clapmt symbol";
    

    clapmt_obj->inforef = clapmt( clapmt_obj->matrix_layout, clapmt_obj->forwrd, clapmt_obj->m,
								clapmt_obj->n,clapmt_obj->xref,
								clapmt_obj->ldx, clapmt_obj->kref);

    /* Compute libflame's Lapacke o/p  */
    clapmt_obj->info = LAPACKE_clapmt( clapmt_obj->matrix_layout, clapmt_obj->forwrd,clapmt_obj->m,
										clapmt_obj->n,clapmt_obj->x, 
										clapmt_obj->ldx, clapmt_obj->k);

    if( clapmt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clapmt is wrong\n", clapmt_obj->info );
    }
    if( clapmt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clapmt is wrong\n", 
        clapmt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clapmt_obj->diff =  computeDiff_c( clapmt_obj->bufsize, 
                clapmt_obj->x, clapmt_obj->xref );

}

TEST_F(clapmt_test, clapmt1) {
    EXPECT_NEAR(0.0, clapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clapmt_test, clapmt2) {
    EXPECT_NEAR(0.0, clapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clapmt_test, clapmt3) {
    EXPECT_NEAR(0.0, clapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clapmt_test, clapmt4) {
    EXPECT_NEAR(0.0, clapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lapmt_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_complex_double* x;
	lapack_int ldx;
	lapack_logical forwrd;
	/*Output Parameter*/
	lapack_int* k;
	lapack_complex_double *xref; 
	lapack_int *kref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      lapmt_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int ldx, lapack_logical forwrd);
      ~lapmt_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
lapmt_dcomplex_parameters:: lapmt_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int ldx_i, lapack_logical forwrd_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldx = ldx_i;
	forwrd = forwrd_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lapmt dcomplex:  m: %d, n: %d ldx: %d \n",  m, n, ldx);
	#endif
	
	bufsize = ldx*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&k, &kref, n);	
	if ((x==NULL) || (xref==NULL) || \
		(k==NULL) || (kref==NULL)){
		EXPECT_FALSE( true) << "lapmt_float_parameters object: malloc error.";
		lapmt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x, xref, bufsize);
	lapacke_gtest_init_int_buffer_pair_with_constant(k, kref, n, forwrd);

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
lapmt_dcomplex_parameters :: ~lapmt_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lapmt_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lapmt_free();

}
/*  Test fixture class definition */
class zlapmt_test  : public  ::testing::Test {
public:
   lapmt_dcomplex_parameters  *zlapmt_obj;
   void SetUp();
   void TearDown () { delete zlapmt_obj; }
};

void zlapmt_test::SetUp(){

    /* LAPACKE zlapmt prototype */
    typedef int (*Fptr_NL_LAPACKE_zlapmt) (int matrix_layout, lapack_logical forwrd, lapack_int m,lapack_int n, 
											lapack_complex_double *x, lapack_int ldx, lapack_int* k);

    Fptr_NL_LAPACKE_zlapmt zlapmt;

    zlapmt_obj = new lapmt_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda, 
						   eig_paramslist[idx].isgn);

    idx = Circular_Increment_Index(idx);

    zlapmt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlapmt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlapmt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlapmt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlapmt = (Fptr_NL_LAPACKE_zlapmt)dlsym(zlapmt_obj->hModule, "LAPACKE_zlapmt");
    ASSERT_TRUE(zlapmt != NULL) << "failed to get the Netlib LAPACKE_zlapmt symbol";
    

    zlapmt_obj->inforef = zlapmt( zlapmt_obj->matrix_layout, zlapmt_obj->forwrd, zlapmt_obj->m,
								zlapmt_obj->n,zlapmt_obj->xref,
								zlapmt_obj->ldx, zlapmt_obj->kref);

    /* Compute libflame's Lapacke o/p  */
    zlapmt_obj->info = LAPACKE_zlapmt( zlapmt_obj->matrix_layout, zlapmt_obj->forwrd, zlapmt_obj->m,
										zlapmt_obj->n,zlapmt_obj->x, 
										zlapmt_obj->ldx, zlapmt_obj->k);

    if( zlapmt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlapmt is wrong\n", zlapmt_obj->info );
    }
    if( zlapmt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlapmt is wrong\n", 
        zlapmt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlapmt_obj->diff =  computeDiff_z( zlapmt_obj->bufsize, 
                zlapmt_obj->x, zlapmt_obj->xref );

}

TEST_F(zlapmt_test, zlapmt1) {
    EXPECT_NEAR(0.0, zlapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlapmt_test, zlapmt2) {
    EXPECT_NEAR(0.0, zlapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlapmt_test, zlapmt3) {
    EXPECT_NEAR(0.0, zlapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlapmt_test, zlapmt4) {
    EXPECT_NEAR(0.0, zlapmt_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


