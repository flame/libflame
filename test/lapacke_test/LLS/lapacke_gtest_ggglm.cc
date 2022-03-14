#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define ggglm_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)  free(b);\
if (bref!=NULL) free(bref); \
if (d!=NULL)    free(d); \
if (dref!=NULL) free(dref);\
if (x!=NULL)    free(x); \
if (xref!=NULL)    free(xref); \
if (y!=NULL)    free(y); \
if (yref!=NULL)    free(yref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class ggglm_float_parameters{

   public:
	int bufsize, bufsize_b, bufsize_x, bufsize_y;
	void *hModule, *dModule;
	float diff, diff_b, diff_x, diff_y;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int p;
	float* A, *b, *d;
	lapack_int lda, ldb;
	/*Output Parameter*/
	float* x, *y;
	float *Aref, *bref, *dref, *xref, *yref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      ggglm_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggglm_float_parameters ();

};

/* Constructor definition  float_common_parameters */
ggglm_float_parameters:: ggglm_float_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int p_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ggglm float: m: %d, n: %d , p: %d, rcond:%d \n",  m, n, p;
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_x = m;
	bufsize_y = p;
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&x, &xref, bufsize_x);
	lapacke_gtest_alloc_float_buffer_pair(&y, &yref, bufsize_y);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, n);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(d==NULL) || (dref==NULL) ||\
		(x==NULL) || (xref==NULL) ||\
		(y==NULL) || (yref==NULL)){
		EXPECT_FALSE( true) << "ggglm_float_parameters object: malloc error.";
		ggglm_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_float_buffer_pair_rand( d, dref, n);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
ggglm_float_parameters :: ~ggglm_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggglm_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggglm_free();

}
/*  Test fixture class definition */
class sggglmest  : public  ::testing::Test {
public:
   ggglm_float_parameters  *sggglm_obj;
   void SetUp();
   void TearDown () { delete sggglm_obj;}
};

void sggglmest::SetUp(){

    /* LAPACKE sggglm prototype */
    typedef int (*Fptr_NL_LAPACKE_sggglm) (int matrix_layout, lapack_int n, lapack_int m,\
	lapack_int p, float* a, lapack_int lda, float* b, lapack_int ldb, float* d, float* x, float* y);

    Fptr_NL_LAPACKE_sggglm sggglm;

    sggglm_obj = new ggglm_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    sggglm_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sggglm_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggglm_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggglm_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sggglm = (Fptr_NL_LAPACKE_sggglm)dlsym(sggglm_obj->hModule, "LAPACKE_sggglm");
    ASSERT_TRUE(sggglm != NULL) << "failed to get the Netlib LAPACKE_sggglm symbol";
    

    sggglm_obj->inforef = sggglm( sggglm_obj->matrix_layout, sggglm_obj->n, sggglm_obj->m,sggglm_obj->p, 
	sggglm_obj->Aref, sggglm_obj->lda, sggglm_obj->bref, sggglm_obj->ldb, sggglm_obj->dref, sggglm_obj->xref, sggglm_obj->yref);

    /* Compute libflame'd Lapacke o/p  */
    sggglm_obj->info = LAPACKE_sggglm( sggglm_obj->matrix_layout, sggglm_obj->n,sggglm_obj->m,sggglm_obj->p, 
	sggglm_obj->A, sggglm_obj->lda, sggglm_obj->b, sggglm_obj->ldb, sggglm_obj->d, sggglm_obj->x, sggglm_obj->y);
	
	if( sggglm_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sggglm is wrong\n", sggglm_obj->info );
    }
    if( sggglm_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sggglm is wrong\n", 
        sggglm_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sggglm_obj->diff =  computeDiff_s( sggglm_obj->bufsize, sggglm_obj->A, sggglm_obj->Aref );
	sggglm_obj->diff_b =  computeDiff_s( sggglm_obj->bufsize_b, sggglm_obj->b, sggglm_obj->bref );
	sggglm_obj->diff_x =  computeDiff_s( sggglm_obj->bufsize_x, sggglm_obj->x, sggglm_obj->xref );
	sggglm_obj->diff_y =  computeDiff_s( sggglm_obj->bufsize_y, sggglm_obj->y, sggglm_obj->yref );

}

TEST_F(sggglmest, sggglm1) {
    EXPECT_NEAR(0.0, sggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sggglmest, sggglm2) {
    EXPECT_NEAR(0.0, sggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sggglmest, sggglm3) {
    EXPECT_NEAR(0.0, sggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sggglmest, sggglm4) {
    EXPECT_NEAR(0.0, sggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}
/* Begin double_common_parameters  class definition */
class ggglm_double_parameters{

   public:
	int bufsize, bufsize_b, bufsize_x, bufsize_y;
	void *hModule, *dModule;
	double diff, diff_b, diff_x, diff_y;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int p;
	double* A, *b, *d;
	lapack_int lda, ldb;
	/*Output Parameter*/
	double* x, *y;
	double *Aref, *bref, *dref, *xref, *yref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      ggglm_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggglm_double_parameters ();

};

/* Constructor definition  double_common_parameters */
ggglm_double_parameters:: ggglm_double_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int p_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ggglm double: m: %d, n: %d , p: %d, rcond:%d \n",  m, n, p;
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_x = m;
	bufsize_y = p;
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&x, &xref, bufsize_x);
	lapacke_gtest_alloc_double_buffer_pair(&y, &yref, bufsize_y);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, n);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(d==NULL) || (dref==NULL) ||\
		(x==NULL) || (xref==NULL) ||\
		(y==NULL) || (yref==NULL)){
		EXPECT_FALSE( true) << "ggglm_double_parameters object: malloc error.";
		ggglm_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_double_buffer_pair_rand( d, dref, n);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
ggglm_double_parameters :: ~ggglm_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggglm_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggglm_free();

}
/*  Test fixture class definition */
class dggglmest  : public  ::testing::Test {
public:
   ggglm_double_parameters  *dggglm_obj;
   void SetUp();
   void TearDown () { delete dggglm_obj;}
};

void dggglmest::SetUp(){

    /* LAPACKE dggglm prototype */
    typedef int (*Fptr_NL_LAPACKE_dggglm) (int matrix_layout, lapack_int n, lapack_int m,\
	lapack_int p, double* a, lapack_int lda, double* b, lapack_int ldb, double* d, double* x, double* y);

    Fptr_NL_LAPACKE_dggglm dggglm;

    dggglm_obj = new ggglm_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    dggglm_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dggglm_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dggglm_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dggglm_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dggglm = (Fptr_NL_LAPACKE_dggglm)dlsym(dggglm_obj->hModule, "LAPACKE_dggglm");
    ASSERT_TRUE(dggglm != NULL) << "failed to get the Netlib LAPACKE_dggglm symbol";
    

    dggglm_obj->inforef = dggglm( dggglm_obj->matrix_layout, dggglm_obj->n, dggglm_obj->m,dggglm_obj->p, 
	dggglm_obj->Aref, dggglm_obj->lda, dggglm_obj->bref, dggglm_obj->ldb, dggglm_obj->dref, dggglm_obj->xref, dggglm_obj->yref);

    /* Compute libflame'd Lapacke o/p  */
    dggglm_obj->info = LAPACKE_dggglm( dggglm_obj->matrix_layout, dggglm_obj->n,dggglm_obj->m,dggglm_obj->p, 
	dggglm_obj->A, dggglm_obj->lda, dggglm_obj->b, dggglm_obj->ldb, dggglm_obj->d, dggglm_obj->x, dggglm_obj->y);
	
	if( dggglm_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dggglm is wrong\n", dggglm_obj->info );
    }
    if( dggglm_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dggglm is wrong\n", 
        dggglm_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dggglm_obj->diff =  computeDiff_d( dggglm_obj->bufsize, dggglm_obj->A, dggglm_obj->Aref );
	dggglm_obj->diff_b =  computeDiff_d( dggglm_obj->bufsize_b, dggglm_obj->b, dggglm_obj->bref );
	dggglm_obj->diff_x =  computeDiff_d( dggglm_obj->bufsize_x, dggglm_obj->x, dggglm_obj->xref );
	dggglm_obj->diff_y =  computeDiff_d( dggglm_obj->bufsize_y, dggglm_obj->y, dggglm_obj->yref );

}

TEST_F(dggglmest, dggglm1) {
    EXPECT_NEAR(0.0, dggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dggglmest, dggglm2) {
    EXPECT_NEAR(0.0, dggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dggglmest, dggglm3) {
    EXPECT_NEAR(0.0, dggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dggglmest, dggglm4) {
    EXPECT_NEAR(0.0, dggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class ggglm_scomplex_parameters{

   public:
	int bufsize, bufsize_b, bufsize_x, bufsize_y;
	void *hModule, *dModule;
	float diff, diff_b, diff_x, diff_y;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int p;
	lapack_complex_float* A, *b, *d;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_float* x, *y;
	lapack_complex_float *Aref, *bref, *dref, *xref, *yref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      ggglm_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggglm_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
ggglm_scomplex_parameters:: ggglm_scomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int p_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ggglm scomplex: m: %d, n: %d , p: %d, rcond:%d \n",  m, n, p;
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_x = m;
	bufsize_y = p;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, bufsize_x);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&y, &yref, bufsize_y);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&d, &dref, n);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(d==NULL) || (dref==NULL) ||\
		(x==NULL) || (xref==NULL) ||\
		(y==NULL) || (yref==NULL)){
		EXPECT_FALSE( true) << "ggglm_scomplex_parameters object: malloc error.";
		ggglm_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_scomplex_buffer_pair_rand( d, dref, n);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
ggglm_scomplex_parameters :: ~ggglm_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggglm_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggglm_free();

}
/*  Test fixture class definition */
class cggglmest  : public  ::testing::Test {
public:
   ggglm_scomplex_parameters  *cggglm_obj;
   void SetUp();
   void TearDown () { delete cggglm_obj;}
};

void cggglmest::SetUp(){

    /* LAPACKE cggglm prototype */
    typedef int (*Fptr_NL_LAPACKE_cggglm) (int matrix_layout, lapack_int n, lapack_int m,\
	lapack_int p, lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, lapack_complex_float* d, lapack_complex_float* x, lapack_complex_float* y);

    Fptr_NL_LAPACKE_cggglm cggglm;

    cggglm_obj = new ggglm_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    cggglm_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cggglm_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cggglm_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cggglm_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cggglm = (Fptr_NL_LAPACKE_cggglm)dlsym(cggglm_obj->hModule, "LAPACKE_cggglm");
    ASSERT_TRUE(cggglm != NULL) << "failed to get the Netlib LAPACKE_cggglm symbol";
    

    cggglm_obj->inforef = cggglm( cggglm_obj->matrix_layout, cggglm_obj->n, cggglm_obj->m,cggglm_obj->p, 
	cggglm_obj->Aref, cggglm_obj->lda, cggglm_obj->bref, cggglm_obj->ldb, cggglm_obj->dref, cggglm_obj->xref, cggglm_obj->yref);

    /* Compute libflame'd Lapacke o/p  */
    cggglm_obj->info = LAPACKE_cggglm( cggglm_obj->matrix_layout, cggglm_obj->n,cggglm_obj->m,cggglm_obj->p, 
	cggglm_obj->A, cggglm_obj->lda, cggglm_obj->b, cggglm_obj->ldb, cggglm_obj->d, cggglm_obj->x, cggglm_obj->y);
	
	if( cggglm_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cggglm is wrong\n", cggglm_obj->info );
    }
    if( cggglm_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cggglm is wrong\n", 
        cggglm_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cggglm_obj->diff =  computeDiff_c( cggglm_obj->bufsize, cggglm_obj->A, cggglm_obj->Aref );
	cggglm_obj->diff_b =  computeDiff_c( cggglm_obj->bufsize_b, cggglm_obj->b, cggglm_obj->bref );
	cggglm_obj->diff_x =  computeDiff_c( cggglm_obj->bufsize_x, cggglm_obj->x, cggglm_obj->xref );
	cggglm_obj->diff_y =  computeDiff_c( cggglm_obj->bufsize_y, cggglm_obj->y, cggglm_obj->yref );

}

TEST_F(cggglmest, cggglm1) {
    EXPECT_NEAR(0.0, cggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cggglmest, cggglm2) {
    EXPECT_NEAR(0.0, cggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cggglmest, cggglm3) {
    EXPECT_NEAR(0.0, cggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cggglmest, cggglm4) {
    EXPECT_NEAR(0.0, cggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class ggglm_dcomplex_parameters{

   public:
	int bufsize, bufsize_b, bufsize_x, bufsize_y;
	void *hModule, *dModule;
	double diff, diff_b, diff_x, diff_y;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int p;
	lapack_complex_double* A, *b, *d;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_double* x, *y;
	lapack_complex_double *Aref, *bref, *dref, *xref, *yref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      ggglm_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggglm_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
ggglm_dcomplex_parameters:: ggglm_dcomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int p_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n ggglm dcomplex: m: %d, n: %d , p: %d, rcond:%d \n",  m, n, p;
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_x = m;
	bufsize_y = p;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, bufsize_x);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&y, &yref, bufsize_y);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&d, &dref, n);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(d==NULL) || (dref==NULL) ||\
		(x==NULL) || (xref==NULL) ||\
		(y==NULL) || (yref==NULL)){
		EXPECT_FALSE( true) << "ggglm_dcomplex_parameters object: malloc error.";
		ggglm_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( d, dref, n);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
ggglm_dcomplex_parameters :: ~ggglm_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggglm_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggglm_free();

}
/*  Test fixture class definition */
class zggglmest  : public  ::testing::Test {
public:
   ggglm_dcomplex_parameters  *zggglm_obj;
   void SetUp();
   void TearDown () { delete zggglm_obj;}
};

void zggglmest::SetUp(){

    /* LAPACKE zggglm prototype */
    typedef int (*Fptr_NL_LAPACKE_zggglm) (int matrix_layout, lapack_int n, lapack_int m,\
	lapack_int p, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, lapack_complex_double* d, lapack_complex_double* x, lapack_complex_double* y);

    Fptr_NL_LAPACKE_zggglm zggglm;

    zggglm_obj = new ggglm_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    zggglm_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zggglm_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zggglm_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zggglm_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zggglm = (Fptr_NL_LAPACKE_zggglm)dlsym(zggglm_obj->hModule, "LAPACKE_zggglm");
    ASSERT_TRUE(zggglm != NULL) << "failed to get the Netlib LAPACKE_zggglm symbol";
    

    zggglm_obj->inforef = zggglm( zggglm_obj->matrix_layout, zggglm_obj->n, zggglm_obj->m,zggglm_obj->p, 
	zggglm_obj->Aref, zggglm_obj->lda, zggglm_obj->bref, zggglm_obj->ldb, zggglm_obj->dref, zggglm_obj->xref, zggglm_obj->yref);

    /* Compute libflame'd Lapacke o/p  */
    zggglm_obj->info = LAPACKE_zggglm( zggglm_obj->matrix_layout, zggglm_obj->n,zggglm_obj->m,zggglm_obj->p, 
	zggglm_obj->A, zggglm_obj->lda, zggglm_obj->b, zggglm_obj->ldb, zggglm_obj->d, zggglm_obj->x, zggglm_obj->y);
	
	if( zggglm_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zggglm is wrong\n", zggglm_obj->info );
    }
    if( zggglm_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zggglm is wrong\n", 
        zggglm_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zggglm_obj->diff =  computeDiff_z( zggglm_obj->bufsize, zggglm_obj->A, zggglm_obj->Aref );
	zggglm_obj->diff_b =  computeDiff_z( zggglm_obj->bufsize_b, zggglm_obj->b, zggglm_obj->bref );
	zggglm_obj->diff_x =  computeDiff_z( zggglm_obj->bufsize_x, zggglm_obj->x, zggglm_obj->xref );
	zggglm_obj->diff_y =  computeDiff_z( zggglm_obj->bufsize_y, zggglm_obj->y, zggglm_obj->yref );

}

TEST_F(zggglmest, zggglm1) {
    EXPECT_NEAR(0.0, zggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zggglmest, zggglm2) {
    EXPECT_NEAR(0.0, zggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zggglmest, zggglm3) {
    EXPECT_NEAR(0.0, zggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zggglmest, zggglm4) {
    EXPECT_NEAR(0.0, zggglm_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_x, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zggglm_obj->diff_y, LAPACKE_GTEST_THRESHOLD);
}