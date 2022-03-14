#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define gglse_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)  free(b);\
if (bref!=NULL) free(bref); \
if (d!=NULL)    free(d); \
if (dref!=NULL) free(dref);\
if (c!=NULL)    free(c); \
if (cref!=NULL)    free(cref); \
if (x!=NULL)    free(x); \
if (xref!=NULL)    free(xref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);\

	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class gglse_float_parameters{

   public:
	int bufsize, bufsize_b, bufsize_x, bufsize_c;
	void *hModule, *dModule;
	float diff, diff_b, diff_x, diff_c;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int p;
	float* A, *b, *d;
	lapack_int lda, ldb;
	/*Output Parameter*/
	float* c, *x;
	float *Aref, *bref, *dref, *cref, *xref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gglse_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~gglse_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gglse_float_parameters:: gglse_float_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gglse float: m: %d, n: %d , p: %d, rcond:%d \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_x = n;
	bufsize_c = m;
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&c, &cref, bufsize_c);
	lapacke_gtest_alloc_float_buffer_pair(&x, &xref, bufsize_x);
	lapacke_gtest_alloc_float_buffer_pair(&d, &dref, p);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(d==NULL) || (dref==NULL) ||\
		(c==NULL) || (cref==NULL) ||\
		(x==NULL) || (xref==NULL)){
		EXPECT_FALSE( true) << "gglse_float_parameters object: malloc error.";
		gglse_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_float_buffer_pair_rand( d, dref, p);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gglse_float_parameters :: ~gglse_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gglse_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gglse_free();

}
/*  Test fixture class definition */
class sgglseest  : public  ::testing::Test {
public:
   gglse_float_parameters  *sgglse_obj;
   void SetUp();
   void TearDown () { delete sgglse_obj;}
};

void sgglseest::SetUp(){

    /* LAPACKE sgglse prototype */
    typedef int (*Fptr_NL_LAPACKE_sgglse) (int matrix_layout, lapack_int m, lapack_int n, lapack_int p,\
	float* a, lapack_int lda, float* b, lapack_int ldb, float* c, float* d, float* x);

    Fptr_NL_LAPACKE_sgglse sgglse;
	
	/* LAPACKE sggrqf prototype */
	//typedef int (*Fptr_NL_LAPACKE_sggrqf)

    sgglse_obj = new gglse_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    sgglse_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgglse_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgglse_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgglse_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgglse = (Fptr_NL_LAPACKE_sgglse)dlsym(sgglse_obj->hModule, "LAPACKE_sgglse");
    ASSERT_TRUE(sgglse != NULL) << "failed to get the Netlib LAPACKE_sgglse symbol";
    

    sgglse_obj->inforef = sgglse( sgglse_obj->matrix_layout, sgglse_obj->m, sgglse_obj->n,sgglse_obj->p, 
	sgglse_obj->Aref, sgglse_obj->lda, sgglse_obj->bref, sgglse_obj->ldb, sgglse_obj->cref, sgglse_obj->dref, sgglse_obj->xref);

    /* Compute libflame'd Lapacke o/p  */
    sgglse_obj->info = LAPACKE_sgglse( sgglse_obj->matrix_layout, sgglse_obj->m,sgglse_obj->n,sgglse_obj->p, 
	sgglse_obj->A, sgglse_obj->lda, sgglse_obj->b, sgglse_obj->ldb, sgglse_obj->c, sgglse_obj->d, sgglse_obj->x);
	
	if( sgglse_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgglse is wrong\n", sgglse_obj->info );
    }
    if( sgglse_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgglse is wrong\n", 
        sgglse_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgglse_obj->diff =  computeDiff_s( sgglse_obj->bufsize, sgglse_obj->A, sgglse_obj->Aref );
	sgglse_obj->diff_b =  computeDiff_s( sgglse_obj->bufsize_b, sgglse_obj->b, sgglse_obj->bref );
	sgglse_obj->diff_x =  computeDiff_s( sgglse_obj->bufsize_x, sgglse_obj->c, sgglse_obj->cref );
	sgglse_obj->diff_c =  computeDiff_s( sgglse_obj->bufsize_c, sgglse_obj->x, sgglse_obj->xref );

}

TEST_F(sgglseest, sgglse1) {
    EXPECT_NEAR(0.0, sgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgglseest, sgglse2) {
    EXPECT_NEAR(0.0, sgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgglseest, sgglse3) {
    EXPECT_NEAR(0.0, sgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(sgglseest, sgglse4) {
    EXPECT_NEAR(0.0, sgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, sgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class gglse_double_parameters{

   public:
	int bufsize, bufsize_b, bufsize_x, bufsize_c;
	void *hModule, *dModule;
	double diff, diff_b, diff_x, diff_c;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int p;
	double* A, *b, *d;
	lapack_int lda, ldb;
	/*Output Parameter*/
	double* c, *x;
	double *Aref, *bref, *dref, *cref, *xref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gglse_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~gglse_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gglse_double_parameters:: gglse_double_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gglse double: m: %d, n: %d , p: %d, rcond:%d \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_x = n;
	bufsize_c = m;
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&c, &cref, bufsize_c);
	lapacke_gtest_alloc_double_buffer_pair(&x, &xref, bufsize_x);
	lapacke_gtest_alloc_double_buffer_pair(&d, &dref, p);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(d==NULL) || (dref==NULL) ||\
		(c==NULL) || (cref==NULL) ||\
		(x==NULL) || (xref==NULL)){
		EXPECT_FALSE( true) << "gglse_double_parameters object: malloc error.";
		gglse_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_double_buffer_pair_rand( d, dref, p);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gglse_double_parameters :: ~gglse_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gglse_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gglse_free();

}
/*  Test fixture class definition */
class dgglseest  : public  ::testing::Test {
public:
   gglse_double_parameters  *dgglse_obj;
   void SetUp();
   void TearDown () { delete dgglse_obj;}
};

void dgglseest::SetUp(){

    /* LAPACKE dgglse prototype */
    typedef int (*Fptr_NL_LAPACKE_dgglse) (int matrix_layout, lapack_int m, lapack_int n, lapack_int p,\
	double* a, lapack_int lda, double* b, lapack_int ldb, double* c, double* d, double* x);

    Fptr_NL_LAPACKE_dgglse dgglse;

    dgglse_obj = new gglse_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n); 
	
	idx = Circular_Increment_Index(idx);

    dgglse_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgglse_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgglse_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgglse_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgglse = (Fptr_NL_LAPACKE_dgglse)dlsym(dgglse_obj->hModule, "LAPACKE_dgglse");
    ASSERT_TRUE(dgglse != NULL) << "failed to get the Netlib LAPACKE_dgglse symbol";
    

    dgglse_obj->inforef = dgglse( dgglse_obj->matrix_layout, dgglse_obj->m, dgglse_obj->n,dgglse_obj->p, 
	dgglse_obj->Aref, dgglse_obj->lda, dgglse_obj->bref, dgglse_obj->ldb, dgglse_obj->cref, dgglse_obj->dref, dgglse_obj->xref);

    /* Compute libflame'd Lapacke o/p  */
    dgglse_obj->info = LAPACKE_dgglse( dgglse_obj->matrix_layout, dgglse_obj->m,dgglse_obj->n,dgglse_obj->p, 
	dgglse_obj->A, dgglse_obj->lda, dgglse_obj->b, dgglse_obj->ldb, dgglse_obj->c, dgglse_obj->d, dgglse_obj->x);
	
	if( dgglse_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgglse is wrong\n", dgglse_obj->info );
    }
    if( dgglse_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgglse is wrong\n", 
        dgglse_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgglse_obj->diff =  computeDiff_d( dgglse_obj->bufsize, dgglse_obj->A, dgglse_obj->Aref );
	dgglse_obj->diff_b =  computeDiff_d( dgglse_obj->bufsize_b, dgglse_obj->b, dgglse_obj->bref );
	dgglse_obj->diff_x =  computeDiff_d( dgglse_obj->bufsize_x, dgglse_obj->c, dgglse_obj->cref );
	dgglse_obj->diff_c =  computeDiff_d( dgglse_obj->bufsize_c, dgglse_obj->x, dgglse_obj->xref );

}

TEST_F(dgglseest, dgglse1) {
    EXPECT_NEAR(0.0, dgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgglseest, dgglse2) {
    EXPECT_NEAR(0.0, dgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgglseest, dgglse3) {
    EXPECT_NEAR(0.0, dgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(dgglseest, dgglse4) {
    EXPECT_NEAR(0.0, dgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, dgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class gglse_scomplex_parameters{

   public:
	int bufsize, bufsize_b, bufsize_x, bufsize_c;
	void *hModule, *dModule;
	float diff, diff_b, diff_x, diff_c;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int p;
	lapack_complex_float* A, *b, *d;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_float* c, *x;
	lapack_complex_float *Aref, *bref, *dref, *cref, *xref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gglse_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~gglse_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
gglse_scomplex_parameters:: gglse_scomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gglse scomplex: m: %d, n: %d , p: %d, rcond:%d \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_x = n;
	bufsize_c = m;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize_c);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, bufsize_x);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&d, &dref, p);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(d==NULL) || (dref==NULL) ||\
		(c==NULL) || (cref==NULL) ||\
		(x==NULL) || (xref==NULL)){
		EXPECT_FALSE( true) << "gglse_scomplex_parameters object: malloc error.";
		gglse_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_scomplex_buffer_pair_rand( d, dref, p);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
gglse_scomplex_parameters :: ~gglse_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gglse_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gglse_free();

}
/*  Test fixture class definition */
class cgglseest  : public  ::testing::Test {
public:
   gglse_scomplex_parameters  *cgglse_obj;
   void SetUp();
   void TearDown () { delete cgglse_obj;}
};

void cgglseest::SetUp(){

    /* LAPACKE cgglse prototype */
    typedef int (*Fptr_NL_LAPACKE_cgglse) (int matrix_layout, lapack_int m, lapack_int n, lapack_int p,\
	lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, lapack_complex_float* c, lapack_complex_float* d, lapack_complex_float* x);

    Fptr_NL_LAPACKE_cgglse cgglse;

    cgglse_obj = new gglse_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n); 
	
	idx = Circular_Increment_Index(idx);

    cgglse_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgglse_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgglse_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgglse_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgglse = (Fptr_NL_LAPACKE_cgglse)dlsym(cgglse_obj->hModule, "LAPACKE_cgglse");
    ASSERT_TRUE(cgglse != NULL) << "failed to get the Netlib LAPACKE_cgglse symbol";
    

    cgglse_obj->inforef = cgglse( cgglse_obj->matrix_layout, cgglse_obj->m, cgglse_obj->n,cgglse_obj->p, 
	cgglse_obj->Aref, cgglse_obj->lda, cgglse_obj->bref, cgglse_obj->ldb, cgglse_obj->cref, cgglse_obj->dref, cgglse_obj->xref);

    /* Compute libflame'd Lapacke o/p  */
    cgglse_obj->info = LAPACKE_cgglse( cgglse_obj->matrix_layout, cgglse_obj->m,cgglse_obj->n,cgglse_obj->p, 
	cgglse_obj->A, cgglse_obj->lda, cgglse_obj->b, cgglse_obj->ldb, cgglse_obj->c, cgglse_obj->d, cgglse_obj->x);
	
	if( cgglse_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgglse is wrong\n", cgglse_obj->info );
    }
    if( cgglse_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgglse is wrong\n", 
        cgglse_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgglse_obj->diff =  computeDiff_c( cgglse_obj->bufsize, cgglse_obj->A, cgglse_obj->Aref );
	cgglse_obj->diff_b =  computeDiff_c( cgglse_obj->bufsize_b, cgglse_obj->b, cgglse_obj->bref );
	cgglse_obj->diff_x =  computeDiff_c( cgglse_obj->bufsize_x, cgglse_obj->c, cgglse_obj->cref );
	cgglse_obj->diff_c =  computeDiff_c( cgglse_obj->bufsize_c, cgglse_obj->x, cgglse_obj->xref );

}

TEST_F(cgglseest, cgglse1) {
    EXPECT_NEAR(0.0, cgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgglseest, cgglse2) {
    EXPECT_NEAR(0.0, cgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgglseest, cgglse3) {
    EXPECT_NEAR(0.0, cgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(cgglseest, cgglse4) {
    EXPECT_NEAR(0.0, cgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, cgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class gglse_dcomplex_parameters{

   public:
	int bufsize, bufsize_b, bufsize_x, bufsize_c;
	void *hModule, *dModule;
	double diff, diff_b, diff_x, diff_c;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int p;
	lapack_complex_double* A, *b, *d;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_double* c, *x;
	lapack_complex_double *Aref, *bref, *dref, *cref, *xref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gglse_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~gglse_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
gglse_dcomplex_parameters:: gglse_dcomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gglse dcomplex: m: %d, n: %d , p: %d, rcond:%d \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_x = n;
	bufsize_c = m;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize_c);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, bufsize_x);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&d, &dref, p);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(d==NULL) || (dref==NULL) ||\
		(c==NULL) || (cref==NULL) ||\
		(x==NULL) || (xref==NULL)){
		EXPECT_FALSE( true) << "gglse_dcomplex_parameters object: malloc error.";
		gglse_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( d, dref, p);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
gglse_dcomplex_parameters :: ~gglse_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gglse_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gglse_free();

}
/*  Test fixture class definition */
class zgglseest  : public  ::testing::Test {
public:
   gglse_dcomplex_parameters  *zgglse_obj;
   void SetUp();
   void TearDown () { delete zgglse_obj;}
};

void zgglseest::SetUp(){

    /* LAPACKE zgglse prototype */
    typedef int (*Fptr_NL_LAPACKE_zgglse) (int matrix_layout, lapack_int m, lapack_int n, lapack_int p,\
	lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, lapack_complex_double* c, lapack_complex_double* d, lapack_complex_double* x);

    Fptr_NL_LAPACKE_zgglse zgglse;

    zgglse_obj = new gglse_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n); 
	
	idx = Circular_Increment_Index(idx);

    zgglse_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgglse_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgglse_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgglse_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgglse = (Fptr_NL_LAPACKE_zgglse)dlsym(zgglse_obj->hModule, "LAPACKE_zgglse");
    ASSERT_TRUE(zgglse != NULL) << "failed to get the Netlib LAPACKE_zgglse symbol";
    

    zgglse_obj->inforef = zgglse( zgglse_obj->matrix_layout, zgglse_obj->m, zgglse_obj->n,zgglse_obj->p, 
	zgglse_obj->Aref, zgglse_obj->lda, zgglse_obj->bref, zgglse_obj->ldb, zgglse_obj->cref, zgglse_obj->dref, zgglse_obj->xref);

    /* Compute libflame'd Lapacke o/p  */
    zgglse_obj->info = LAPACKE_zgglse( zgglse_obj->matrix_layout, zgglse_obj->m,zgglse_obj->n,zgglse_obj->p, 
	zgglse_obj->A, zgglse_obj->lda, zgglse_obj->b, zgglse_obj->ldb, zgglse_obj->c, zgglse_obj->d, zgglse_obj->x);
	
	if( zgglse_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgglse is wrong\n", zgglse_obj->info );
    }
    if( zgglse_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgglse is wrong\n", 
        zgglse_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgglse_obj->diff =  computeDiff_z( zgglse_obj->bufsize, zgglse_obj->A, zgglse_obj->Aref );
	zgglse_obj->diff_b =  computeDiff_z( zgglse_obj->bufsize_b, zgglse_obj->b, zgglse_obj->bref );
	zgglse_obj->diff_x =  computeDiff_z( zgglse_obj->bufsize_x, zgglse_obj->c, zgglse_obj->cref );
	zgglse_obj->diff_c =  computeDiff_z( zgglse_obj->bufsize_c, zgglse_obj->x, zgglse_obj->xref );

}

TEST_F(zgglseest, zgglse1) {
    EXPECT_NEAR(0.0, zgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, zgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, zgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgglseest, zgglse2) {
    EXPECT_NEAR(0.0, zgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgglseest, zgglse3) {
    EXPECT_NEAR(0.0, zgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, zgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}

TEST_F(zgglseest, zgglse4) {
    EXPECT_NEAR(0.0, zgglse_obj->diff,LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zgglse_obj->diff_b,LAPACKE_EIG_THRESHOLD);
//	EXPECT_NEAR(0.0, zgglse_obj->diff_x,LAPACKE_EIG_THRESHOLD);
	//EXPECT_NEAR(0.0, zgglse_obj->diff_c,LAPACKE_EIG_THRESHOLD);
}