#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define getsls_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)  free(b);\
if (bref!=NULL) free(bref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class getsls_float_parameters{

   public:
	int bufsize, bufsize_b;
	void *hModule, *dModule;
	float diff, diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char trans;
	lapack_int nrhs;
	float* A;
	lapack_int lda, ldb;
	/*Output Parameter*/
	float* b;
	float *Aref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      getsls_float_parameters (int matrix_layout, char trans, lapack_int m , lapack_int n, lapack_int nrhs);
      ~getsls_float_parameters ();

};

/* Constructor definition  float_common_parameters */
getsls_float_parameters:: getsls_float_parameters (int matrix_layout_i,char trans_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	trans = trans_i;
	nrhs = nrhs_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n getsls float: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = max(m,n);
		bufsize = lda*n;
		bufsize_b = ldb*nrhs;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
		ldb =nrhs;
		bufsize = lda*m;
		bufsize_b = ldb*max(m,n);
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	if (trans == 'C')
		trans = 'N';
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "getsls_float_parameters object: malloc error.";
		getsls_free();
		exit(0);
	}	

	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
getsls_float_parameters :: ~getsls_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" getsls_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   getsls_free();

}
/*  Test fixture class definition */
class sgetslsest  : public  ::testing::Test {
public:
   getsls_float_parameters  *sgetsls_obj;
   void SetUp();
   void TearDown () { delete sgetsls_obj;}
};

void sgetslsest::SetUp(){

    /* LAPACKE sgetsls prototype */
    typedef int (*Fptr_NL_LAPACKE_sgetsls) (int matrix_layout, char trans, lapack_int m,\
	lapack_int n, lapack_int nrhs, float* a, lapack_int lda, float* b, lapack_int ldb);

    Fptr_NL_LAPACKE_sgetsls sgetsls;

    sgetsls_obj = new getsls_float_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].trans,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs); 
	
	idx = Circular_Increment_Index(idx);

    sgetsls_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgetsls_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgetsls_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgetsls_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgetsls = (Fptr_NL_LAPACKE_sgetsls)dlsym(sgetsls_obj->hModule, "LAPACKE_sgetsls");
    ASSERT_TRUE(sgetsls != NULL) << "failed to get the Netlib LAPACKE_sgetsls symbol";
    

    sgetsls_obj->inforef = sgetsls( sgetsls_obj->matrix_layout, sgetsls_obj->trans, sgetsls_obj->m, sgetsls_obj->n,sgetsls_obj->nrhs, 
	sgetsls_obj->Aref, sgetsls_obj->lda, sgetsls_obj->bref, sgetsls_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    sgetsls_obj->info = LAPACKE_sgetsls( sgetsls_obj->matrix_layout,  sgetsls_obj->trans, sgetsls_obj->m,sgetsls_obj->n,sgetsls_obj->nrhs, 
	sgetsls_obj->A, sgetsls_obj->lda, sgetsls_obj->b, sgetsls_obj->ldb);
	
	if( sgetsls_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgetsls is wrong\n", sgetsls_obj->info );
    }
    if( sgetsls_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgetsls is wrong\n", 
        sgetsls_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgetsls_obj->diff =  computeDiff_s( sgetsls_obj->bufsize, sgetsls_obj->A, sgetsls_obj->Aref );
	sgetsls_obj->diff_b =  computeDiff_s( sgetsls_obj->bufsize_b, sgetsls_obj->b, sgetsls_obj->bref );


}

TEST_F(sgetslsest, sgetsls1) {
    EXPECT_NEAR(0.0, sgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(sgetslsest, sgetsls2) {
    EXPECT_NEAR(0.0, sgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgetslsest, sgetsls3) {
    EXPECT_NEAR(0.0, sgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgetslsest, sgetsls4) {
    EXPECT_NEAR(0.0, sgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);

}

/* Begin double_common_parameters  class definition */
class getsls_double_parameters{

   public:
	int bufsize, bufsize_b;
	void *hModule, *dModule;
	double diff, diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char trans;
	lapack_int nrhs;
	double* A;
	lapack_int lda, ldb;
	/*Output Parameter*/
	double* b;
	double *Aref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      getsls_double_parameters (int matrix_layout, char trans, lapack_int m , lapack_int n, lapack_int nrhs);
      ~getsls_double_parameters ();

};

/* Constructor definition  double_common_parameters */
getsls_double_parameters:: getsls_double_parameters (int matrix_layout_i,char trans_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	trans = trans_i;
	nrhs = nrhs_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n getsls double: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = max(m,n);
		bufsize = lda*n;
		bufsize_b = ldb*nrhs;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
		ldb =nrhs;
		bufsize = lda*m;
		bufsize_b = ldb*max(m,n);
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}	
	if (trans == 'C')
		trans = 'N';
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "getsls_double_parameters object: malloc error.";
		getsls_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
getsls_double_parameters :: ~getsls_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" getsls_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   getsls_free();

}
/*  Test fixture class definition */
class dgetslsest  : public  ::testing::Test {
public:
   getsls_double_parameters  *dgetsls_obj;
   void SetUp();
   void TearDown () { delete dgetsls_obj;}
};

void dgetslsest::SetUp(){

    /* LAPACKE dgetsls prototype */
    typedef int (*Fptr_NL_LAPACKE_dgetsls) (int matrix_layout, char trans, lapack_int m,\
	lapack_int n, lapack_int nrhs, double* a, lapack_int lda, double* b, lapack_int ldb);

    Fptr_NL_LAPACKE_dgetsls dgetsls;

    dgetsls_obj = new getsls_double_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].trans,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs); 
	
	idx = Circular_Increment_Index(idx);

    dgetsls_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgetsls_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgetsls_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgetsls_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgetsls = (Fptr_NL_LAPACKE_dgetsls)dlsym(dgetsls_obj->hModule, "LAPACKE_dgetsls");
    ASSERT_TRUE(dgetsls != NULL) << "failed to get the Netlib LAPACKE_dgetsls symbol";
    

    dgetsls_obj->inforef = dgetsls( dgetsls_obj->matrix_layout, dgetsls_obj->trans, dgetsls_obj->m, dgetsls_obj->n,dgetsls_obj->nrhs, 
	dgetsls_obj->Aref, dgetsls_obj->lda, dgetsls_obj->bref, dgetsls_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    dgetsls_obj->info = LAPACKE_dgetsls( dgetsls_obj->matrix_layout,  dgetsls_obj->trans, dgetsls_obj->m,dgetsls_obj->n,dgetsls_obj->nrhs, 
	dgetsls_obj->A, dgetsls_obj->lda, dgetsls_obj->b, dgetsls_obj->ldb);
	
	if( dgetsls_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgetsls is wrong\n", dgetsls_obj->info );
    }
    if( dgetsls_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgetsls is wrong\n", 
        dgetsls_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgetsls_obj->diff =  computeDiff_d( dgetsls_obj->bufsize, dgetsls_obj->A, dgetsls_obj->Aref );
	dgetsls_obj->diff_b =  computeDiff_d( dgetsls_obj->bufsize_b, dgetsls_obj->b, dgetsls_obj->bref );


}

TEST_F(dgetslsest, dgetsls1) {
    EXPECT_NEAR(0.0, dgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(dgetslsest, dgetsls2) {
    EXPECT_NEAR(0.0, dgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetslsest, dgetsls3) {
    EXPECT_NEAR(0.0, dgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgetslsest, dgetsls4) {
    EXPECT_NEAR(0.0, dgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);

}




/* Begin scomplex_common_parameters  class definition */
class getsls_scomplex_parameters{

   public:
	int bufsize, bufsize_b;
	void *hModule, *dModule;
	float diff, diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nrhs;
	lapack_complex_float* A;
	char trans;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_float* b;
	lapack_complex_float *Aref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      getsls_scomplex_parameters (int matrix_layout, char trans, lapack_int m , lapack_int n, lapack_int nrhs);
      ~getsls_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
getsls_scomplex_parameters:: getsls_scomplex_parameters (int matrix_layout_i, char trans_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	trans = trans_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n getsls scomplex: m: %d, n: %d , nrhs: %d\n",  m, n, nrhs);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = max(m,n);
		bufsize = lda*n;
		bufsize_b = ldb*nrhs;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
		ldb =nrhs;
		bufsize = lda*m;
		bufsize_b = ldb*max(m,n);
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}	
	trans = 'N'; // 'T' is not supported
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "getsls_scomplex_parameters object: malloc error.";
		getsls_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
getsls_scomplex_parameters :: ~getsls_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" getsls_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   getsls_free();

}
/*  Test fixture class definition */
class cgetslsest  : public  ::testing::Test {
public:
   getsls_scomplex_parameters  *cgetsls_obj;
   void SetUp();
   void TearDown () { delete cgetsls_obj;}
};

void cgetslsest::SetUp(){

    /* LAPACKE cgetsls prototype */
    typedef int (*Fptr_NL_LAPACKE_cgetsls) (int matrix_layout, char trans, lapack_int m,\
                            lapack_int n, lapack_int nrhs,\
                            lapack_complex_float* a, lapack_int lda,\
                            lapack_complex_float* b, lapack_int ldb);

    Fptr_NL_LAPACKE_cgetsls cgetsls;

    cgetsls_obj = new getsls_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].trans,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs); 
	
	idx = Circular_Increment_Index(idx);

    cgetsls_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgetsls_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgetsls_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgetsls_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgetsls = (Fptr_NL_LAPACKE_cgetsls)dlsym(cgetsls_obj->hModule, "LAPACKE_cgetsls");
    ASSERT_TRUE(cgetsls != NULL) << "failed to get the Netlib LAPACKE_cgetsls symbol";
    

    cgetsls_obj->inforef = cgetsls( cgetsls_obj->matrix_layout, cgetsls_obj->trans,  cgetsls_obj->m, cgetsls_obj->n,cgetsls_obj->nrhs, 
	cgetsls_obj->Aref, cgetsls_obj->lda, cgetsls_obj->bref, cgetsls_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    cgetsls_obj->info = LAPACKE_cgetsls( cgetsls_obj->matrix_layout,  cgetsls_obj->trans, cgetsls_obj->m,cgetsls_obj->n,cgetsls_obj->nrhs, 
	cgetsls_obj->A, cgetsls_obj->lda, cgetsls_obj->b, cgetsls_obj->ldb);
	
	if( cgetsls_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgetsls is wrong\n", cgetsls_obj->info );
    }
    if( cgetsls_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgetsls is wrong\n", 
        cgetsls_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgetsls_obj->diff =  computeDiff_c( cgetsls_obj->bufsize, cgetsls_obj->A, cgetsls_obj->Aref );
	cgetsls_obj->diff_b =  computeDiff_c( cgetsls_obj->bufsize_b, cgetsls_obj->b, cgetsls_obj->bref );


}

TEST_F(cgetslsest, cgetsls1) {
    EXPECT_NEAR(0.0, cgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgetslsest, cgetsls2) {
    EXPECT_NEAR(0.0, cgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgetslsest, cgetsls3) {
    EXPECT_NEAR(0.0, cgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgetslsest, cgetsls4) {
    EXPECT_NEAR(0.0, cgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	
}

/* Begin dcomplex_common_parameters  class definition */
class getsls_dcomplex_parameters{

   public:
	int bufsize, bufsize_b;
	void *hModule, *dModule;
	double diff, diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nrhs;
	lapack_complex_double* A;
	char trans;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_double* b;
	lapack_complex_double *Aref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      getsls_dcomplex_parameters (int matrix_layout, char trans, lapack_int m , lapack_int n, lapack_int nrhs);
      ~getsls_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
getsls_dcomplex_parameters:: getsls_dcomplex_parameters (int matrix_layout_i, char trans_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	trans = trans_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n getsls dcomplex: m: %d, n: %d , nrhs: %d\n",  m, n, nrhs);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = max(m,n);
		bufsize = lda*n;
		bufsize_b = ldb*nrhs;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		lda = n;
		ldb =nrhs;
		bufsize = lda*m;
		bufsize_b = ldb*max(m,n);
	}else {
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	trans = 'N'; // 'T' is not supported
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "getsls_dcomplex_parameters object: malloc error.";
		getsls_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
getsls_dcomplex_parameters :: ~getsls_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" getsls_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   getsls_free();

}
/*  Test fixture class definition */
class zgetslsest  : public  ::testing::Test {
public:
   getsls_dcomplex_parameters  *zgetsls_obj;
   void SetUp();
   void TearDown () { delete zgetsls_obj;}
};

void zgetslsest::SetUp(){

    /* LAPACKE zgetsls prototype */
    typedef int (*Fptr_NL_LAPACKE_zgetsls) (int matrix_layout, char trans, lapack_int m,\
                            lapack_int n, lapack_int nrhs,\
                            lapack_complex_double* a, lapack_int lda,\
                            lapack_complex_double* b, lapack_int ldb);

    Fptr_NL_LAPACKE_zgetsls zgetsls;

    zgetsls_obj = new getsls_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
							eig_paramslist[idx].trans,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs); 
	
	idx = Circular_Increment_Index(idx);

    zgetsls_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgetsls_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgetsls_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgetsls_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgetsls = (Fptr_NL_LAPACKE_zgetsls)dlsym(zgetsls_obj->hModule, "LAPACKE_zgetsls");
    ASSERT_TRUE(zgetsls != NULL) << "failed to get the Netlib LAPACKE_zgetsls symbol";
    

    zgetsls_obj->inforef = zgetsls( zgetsls_obj->matrix_layout, zgetsls_obj->trans,  zgetsls_obj->m, zgetsls_obj->n,zgetsls_obj->nrhs, 
	zgetsls_obj->Aref, zgetsls_obj->lda, zgetsls_obj->bref, zgetsls_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    zgetsls_obj->info = LAPACKE_zgetsls( zgetsls_obj->matrix_layout,  zgetsls_obj->trans, zgetsls_obj->m,zgetsls_obj->n,zgetsls_obj->nrhs, 
	zgetsls_obj->A, zgetsls_obj->lda, zgetsls_obj->b, zgetsls_obj->ldb);
	
	if( zgetsls_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgetsls is wrong\n", zgetsls_obj->info );
    }
    if( zgetsls_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgetsls is wrong\n", 
        zgetsls_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgetsls_obj->diff =  computeDiff_z( zgetsls_obj->bufsize, zgetsls_obj->A, zgetsls_obj->Aref );
	zgetsls_obj->diff_b =  computeDiff_z( zgetsls_obj->bufsize_b, zgetsls_obj->b, zgetsls_obj->bref );


}

TEST_F(zgetslsest, zgetsls1) {
    EXPECT_NEAR(0.0, zgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgetslsest, zgetsls2) {
    EXPECT_NEAR(0.0, zgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgetslsest, zgetsls3) {
    EXPECT_NEAR(0.0, zgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgetslsest, zgetsls4) {
    EXPECT_NEAR(0.0, zgetsls_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgetsls_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	
}