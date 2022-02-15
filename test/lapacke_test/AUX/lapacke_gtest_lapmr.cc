#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define lapmr_free() \
if (x!=NULL)    free(x); \
if (xref!=NULL) free(xref);\
if (k!=NULL)  free(k);\
if (kref!=NULL) free(kref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class lapmr_float_parameters{

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
      lapmr_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int ldx, lapack_logical forwrd);
      ~lapmr_float_parameters ();

};

/* Constructor definition  float_common_parameters */
lapmr_float_parameters:: lapmr_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int ldx_i, lapack_logical forwrd_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldx = ldx_i;
	forwrd = forwrd_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lapmr float:  m: %d, n: %d ldx: %d, forwrd:%d: \n",  m, n, ldx, forwrd);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	ldx = m;
		bufsize = ldx*n;
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	ldx = n;
		bufsize = ldx*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&k, &kref, m);	
	if ((x==NULL) || (xref==NULL) || (k==NULL) || (kref==NULL)){
		EXPECT_FALSE( true) << "lapmr_float_parameters object: malloc error.";
		lapmr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( x, xref, bufsize);
	lapacke_gtest_init_int_buffer_pair_with_constant(k, kref, m, forwrd);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lapmr_float_parameters :: ~lapmr_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lapmr_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lapmr_free();

}
/*  Test fixture class definition */
class slapmr_test  : public  ::testing::Test {
public:
   lapmr_float_parameters  *slapmr_obj;
   void SetUp();
   void TearDown () { delete slapmr_obj; }
};

void slapmr_test::SetUp(){

    /* LAPACKE slapmr prototype */
    typedef int (*Fptr_NL_LAPACKE_slapmr) (int matrix_layout, lapack_logical forwrd, lapack_int m,lapack_int n, 
											float* x, lapack_int ldx, lapack_int* k);

    Fptr_NL_LAPACKE_slapmr slapmr;

    slapmr_obj = new lapmr_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda,
						   eig_paramslist[idx].isgn);

    idx = Circular_Increment_Index(idx);
	

    slapmr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slapmr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slapmr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slapmr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slapmr = (Fptr_NL_LAPACKE_slapmr)dlsym(slapmr_obj->hModule, "LAPACKE_slapmr");
    ASSERT_TRUE(slapmr != NULL) << "failed to get the Netlib LAPACKE_slapmr symbol";
    

    slapmr_obj->inforef = slapmr( slapmr_obj->matrix_layout, slapmr_obj->forwrd, slapmr_obj->m, \
								slapmr_obj->n,slapmr_obj->xref, slapmr_obj->ldx, slapmr_obj->kref);

    /* Compute libflame's Lapacke o/p  */
    slapmr_obj->info = LAPACKE_slapmr( slapmr_obj->matrix_layout, slapmr_obj->forwrd, slapmr_obj->m,
										slapmr_obj->n,slapmr_obj->x, slapmr_obj->ldx, slapmr_obj->k); 

    if( slapmr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slapmr is wrong\n", slapmr_obj->info );
    }
    if( slapmr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slapmr is wrong\n", 
        slapmr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slapmr_obj->diff =  computeDiff_s( slapmr_obj->bufsize, slapmr_obj->x, slapmr_obj->xref );

}

TEST_F(slapmr_test, slapmr1) {
    EXPECT_NEAR(0.0, slapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slapmr_test, slapmr2) {
    EXPECT_NEAR(0.0, slapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slapmr_test, slapmr3) {
    EXPECT_NEAR(0.0, slapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slapmr_test, slapmr4) {
    EXPECT_NEAR(0.0, slapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class lapmr_double_parameters{

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
      lapmr_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int ldx, lapack_logical forwrd);
      ~lapmr_double_parameters ();

};

/* Constructor definition  double_common_parameters */
lapmr_double_parameters:: lapmr_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int ldx_i, lapack_logical forwrd_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldx = ldx_i;
	forwrd = forwrd_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lapmr double:  m: %d, n: %d ldx: %d \n",  m, n, ldx);
	#endif
	

	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldx*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldx*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&k, &kref, m);	
	if ((x==NULL) || (xref==NULL) || \
		(k==NULL) || (kref==NULL)){
		EXPECT_FALSE( true) << "lapmr_double_parameters object: malloc error.";
		lapmr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( x, xref, bufsize);
	lapacke_gtest_init_int_buffer_pair_with_constant(k, kref, m, forwrd);


} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lapmr_double_parameters :: ~lapmr_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lapmr_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lapmr_free();

}
/*  Test fixture class definition */
class dlapmr_test  : public  ::testing::Test {
public:
   lapmr_double_parameters  *dlapmr_obj;
   void SetUp();
   void TearDown () { delete dlapmr_obj; }
};

void dlapmr_test::SetUp(){

    /* LAPACKE dlapmr prototype */
    typedef int (*Fptr_NL_LAPACKE_dlapmr) (int matrix_layout, lapack_logical forwrd, lapack_int m,lapack_int n, 
											double *x, lapack_int ldx, lapack_int* k);

    Fptr_NL_LAPACKE_dlapmr dlapmr;

    dlapmr_obj = new lapmr_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda, 
						   eig_paramslist[idx].isgn);

    idx = Circular_Increment_Index(idx);
    dlapmr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlapmr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlapmr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlapmr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlapmr = (Fptr_NL_LAPACKE_dlapmr)dlsym(dlapmr_obj->hModule, "LAPACKE_dlapmr");
    ASSERT_TRUE(dlapmr != NULL) << "failed to get the Netlib LAPACKE_dlapmr symbol";
    

    dlapmr_obj->inforef = dlapmr( dlapmr_obj->matrix_layout, dlapmr_obj->forwrd, dlapmr_obj->m,
								dlapmr_obj->n,dlapmr_obj->xref,
								dlapmr_obj->ldx, dlapmr_obj->kref);

    /* Compute libflame's Lapacke o/p  */
    dlapmr_obj->info = LAPACKE_dlapmr( dlapmr_obj->matrix_layout, dlapmr_obj->forwrd, dlapmr_obj->m,
										dlapmr_obj->n,dlapmr_obj->x, 
										dlapmr_obj->ldx, dlapmr_obj->k);

    if( dlapmr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlapmr is wrong\n", dlapmr_obj->info );
    }
    if( dlapmr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlapmr is wrong\n", 
        dlapmr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlapmr_obj->diff =  computeDiff_d( dlapmr_obj->bufsize, 
                dlapmr_obj->x, dlapmr_obj->xref );

}

TEST_F(dlapmr_test, dlapmr1) {
    EXPECT_NEAR(0.0, dlapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlapmr_test, dlapmr2) {
    EXPECT_NEAR(0.0, dlapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlapmr_test, dlapmr3) {
    EXPECT_NEAR(0.0, dlapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlapmr_test, dlapmr4) {
    EXPECT_NEAR(0.0, dlapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class lapmr_scomplex_parameters{

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
      lapmr_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int ldx, lapack_logical forwrd);
      ~lapmr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lapmr_scomplex_parameters:: lapmr_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int ldx_i, lapack_logical forwrd_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldx = ldx_i;
	forwrd = forwrd_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lapmr scomplex:  m: %d, n: %d ldx: %d \n",  m, n, ldx);
	#endif
	
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldx*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldx*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&k, &kref, m);	
	if ((x==NULL) || (xref==NULL) || \
		(k==NULL) || (kref==NULL)){
		EXPECT_FALSE( true) << "lapmr_float_parameters object: malloc error.";
		lapmr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( x, xref, bufsize);
	lapacke_gtest_init_int_buffer_pair_with_constant(k, kref, m, forwrd);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lapmr_scomplex_parameters :: ~lapmr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lapmr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lapmr_free();

}
/*  Test fixture class definition */
class clapmr_test  : public  ::testing::Test {
public:
   lapmr_scomplex_parameters  *clapmr_obj;
   void SetUp();
   void TearDown () { delete clapmr_obj; }
};

void clapmr_test::SetUp(){

    /* LAPACKE clapmr prototype */
    typedef int (*Fptr_NL_LAPACKE_clapmr) (int matrix_layout, lapack_logical forwrd, lapack_int m,lapack_int n, 
											lapack_complex_float *x, lapack_int ldx, lapack_int* k);

    Fptr_NL_LAPACKE_clapmr clapmr;

    clapmr_obj = new lapmr_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda, 
						   eig_paramslist[idx].isgn);

    idx = Circular_Increment_Index(idx);

    clapmr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clapmr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clapmr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clapmr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clapmr = (Fptr_NL_LAPACKE_clapmr)dlsym(clapmr_obj->hModule, "LAPACKE_clapmr");
    ASSERT_TRUE(clapmr != NULL) << "failed to get the Netlib LAPACKE_clapmr symbol";
    

    clapmr_obj->inforef = clapmr( clapmr_obj->matrix_layout, clapmr_obj->forwrd, clapmr_obj->m,
								clapmr_obj->n,clapmr_obj->xref,
								clapmr_obj->ldx, clapmr_obj->kref);

    /* Compute libflame's Lapacke o/p  */
    clapmr_obj->info = LAPACKE_clapmr( clapmr_obj->matrix_layout, clapmr_obj->forwrd,clapmr_obj->m,
										clapmr_obj->n,clapmr_obj->x, 
										clapmr_obj->ldx, clapmr_obj->k);

    if( clapmr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clapmr is wrong\n", clapmr_obj->info );
    }
    if( clapmr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clapmr is wrong\n", 
        clapmr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clapmr_obj->diff =  computeDiff_c( clapmr_obj->bufsize, 
                clapmr_obj->x, clapmr_obj->xref );

}

TEST_F(clapmr_test, clapmr1) {
    EXPECT_NEAR(0.0, clapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clapmr_test, clapmr2) {
    EXPECT_NEAR(0.0, clapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clapmr_test, clapmr3) {
    EXPECT_NEAR(0.0, clapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clapmr_test, clapmr4) {
    EXPECT_NEAR(0.0, clapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lapmr_dcomplex_parameters{

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
      lapmr_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int ldx, lapack_logical forwrd);
      ~lapmr_dcomplex_parameters ();

};

/* Constructor definition  domplex_common_parameters */
lapmr_dcomplex_parameters:: lapmr_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int ldx_i, lapack_logical forwrd_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	ldx = ldx_i;
	forwrd = forwrd_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lapmr dcomplex:  m: %d, n: %d ldx: %d \n",  m, n, ldx);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldx*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldx*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&k, &kref, m);	
	if ((x==NULL) || (xref==NULL) || \
		(k==NULL) || (kref==NULL)){
		EXPECT_FALSE( true) << "lapmr_float_parameters object: malloc error.";
		lapmr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( x, xref, bufsize);
	lapacke_gtest_init_int_buffer_pair_with_constant(k, kref, m, forwrd);

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
lapmr_dcomplex_parameters :: ~lapmr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lapmr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lapmr_free();

}
/*  Test fixture class definition */
class zlapmr_test  : public  ::testing::Test {
public:
   lapmr_dcomplex_parameters  *zlapmr_obj;
   void SetUp();
   void TearDown () { delete zlapmr_obj; }
};

void zlapmr_test::SetUp(){

    /* LAPACKE zlapmr prototype */
    typedef int (*Fptr_NL_LAPACKE_zlapmr) (int matrix_layout, lapack_logical forwrd, lapack_int m,lapack_int n, 
											lapack_complex_double *x, lapack_int ldx, lapack_int* k);

    Fptr_NL_LAPACKE_zlapmr zlapmr;

    zlapmr_obj = new lapmr_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].lda, 
						   eig_paramslist[idx].isgn);

    idx = Circular_Increment_Index(idx);

    zlapmr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlapmr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlapmr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlapmr_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlapmr = (Fptr_NL_LAPACKE_zlapmr)dlsym(zlapmr_obj->hModule, "LAPACKE_zlapmr");
    ASSERT_TRUE(zlapmr != NULL) << "failed to get the Netlib LAPACKE_zlapmr symbol";
    

    zlapmr_obj->inforef = zlapmr( zlapmr_obj->matrix_layout, zlapmr_obj->forwrd, zlapmr_obj->m,
								zlapmr_obj->n,zlapmr_obj->xref,
								zlapmr_obj->ldx, zlapmr_obj->kref);

    /* Compute libflame's Lapacke o/p  */
    zlapmr_obj->info = LAPACKE_zlapmr( zlapmr_obj->matrix_layout, zlapmr_obj->forwrd, zlapmr_obj->m,
										zlapmr_obj->n,zlapmr_obj->x, 
										zlapmr_obj->ldx, zlapmr_obj->k);

    if( zlapmr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlapmr is wrong\n", zlapmr_obj->info );
    }
    if( zlapmr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlapmr is wrong\n", 
        zlapmr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlapmr_obj->diff =  computeDiff_z( zlapmr_obj->bufsize, 
                zlapmr_obj->x, zlapmr_obj->xref );

}

TEST_F(zlapmr_test, zlapmr1) {
    EXPECT_NEAR(0.0, zlapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlapmr_test, zlapmr2) {
    EXPECT_NEAR(0.0, zlapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlapmr_test, zlapmr3) {
    EXPECT_NEAR(0.0, zlapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlapmr_test, zlapmr4) {
    EXPECT_NEAR(0.0, zlapmr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


