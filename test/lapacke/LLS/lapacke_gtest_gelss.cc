#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define gelss_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)  free(b);\
if (bref!=NULL) free(bref); \
if (s!=NULL)    free(s); \
if (sref!=NULL) free(sref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class gelss_float_parameters{

   public:
	int bufsize, bufsize_b, bufsize_s;
	void *hModule, *dModule;
	float diff, diff_b, diff_s;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nrhs, rank;
	float* A, rcond;
	lapack_int lda, ldb;
	/*Output Parameter*/
	float* b, *s;
	float *Aref, *bref, *sref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gelss_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, float rcond);
      ~gelss_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gelss_float_parameters:: gelss_float_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, float rcond_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelss float: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = fla_max(m,n);
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
	
	bufsize_s = fla_min(m,n);
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&s, &sref, bufsize_s);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(s==NULL) || (sref==NULL)){
		EXPECT_FALSE( true) << "gelss_float_parameters object: malloc error.";
		gelss_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gelss_float_parameters :: ~gelss_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelss_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelss_free();

}
/*  Test fixture class definition */
class sgelssest  : public  ::testing::Test {
public:
   gelss_float_parameters  *sgelss_obj;
   void SetUp();
   void TearDown () { delete sgelss_obj;}
};

void sgelssest::SetUp(){

    /* LAPACKE sgelss prototype */
    typedef int (*Fptr_NL_LAPACKE_sgelss) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, float* a, lapack_int lda, float* b, lapack_int ldb, float* s, float rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_sgelss sgelss;

    sgelss_obj = new gelss_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].isgn); 
	
	idx = Circular_Increment_Index(idx);

    sgelss_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgelss_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgelss_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgelss_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgelss = (Fptr_NL_LAPACKE_sgelss)dlsym(sgelss_obj->hModule, "LAPACKE_sgelss");
    ASSERT_TRUE(sgelss != NULL) << "failed to get the Netlib LAPACKE_sgelss symbol";
    

    sgelss_obj->inforef = sgelss( sgelss_obj->matrix_layout, sgelss_obj->m, sgelss_obj->n,sgelss_obj->nrhs, 
	sgelss_obj->Aref, sgelss_obj->lda, sgelss_obj->bref, sgelss_obj->ldb, sgelss_obj->sref, sgelss_obj->rcond, &sgelss_obj->rank);

    /* Compute libflame's Lapacke o/p  */
    sgelss_obj->info = LAPACKE_sgelss( sgelss_obj->matrix_layout, sgelss_obj->m,sgelss_obj->n,sgelss_obj->nrhs, 
	sgelss_obj->A, sgelss_obj->lda, sgelss_obj->b, sgelss_obj->ldb, sgelss_obj->s, sgelss_obj->rcond, &sgelss_obj->rank);
	
	if( sgelss_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgelss is wrong\n", sgelss_obj->info );
    }
    if( sgelss_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgelss is wrong\n", 
        sgelss_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    //sgelss_obj->diff =  computeDiff_s( sgelss_obj->bufsize, sgelss_obj->A, sgelss_obj->Aref );
	//sgelss_obj->diff_b =  computeDiff_s( sgelss_obj->bufsize_b, sgelss_obj->b, sgelss_obj->bref );
	sgelss_obj->diff_s =  computeDiff_s( sgelss_obj->bufsize_s, sgelss_obj->s, sgelss_obj->sref );

}

TEST_F(sgelssest, sgelss1) {
    //EXPECT_NEAR(0.0, sgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, sgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelssest, sgelss2) {
    //EXPECT_NEAR(0.0, sgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, sgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelssest, sgelss3) {
    //EXPECT_NEAR(0.0, sgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, sgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelssest, sgelss4) {
    //EXPECT_NEAR(0.0, sgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, sgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class gelss_double_parameters{

   public:
	int bufsize, bufsize_b, bufsize_s;
	void *hModule, *dModule;
	double diff, diff_b, diff_s;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nrhs, rank;
	double* A, rcond;
	lapack_int lda, ldb;
	/*Output Parameter*/
	double* b, *s;
	double *Aref, *bref, *sref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gelss_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, double rcond);
      ~gelss_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gelss_double_parameters:: gelss_double_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, double rcond_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelss double: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = fla_max(m,n);
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
	
	bufsize_s = fla_min(m,n);
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&s, &sref, bufsize_s);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(s==NULL) || (sref==NULL)){
		EXPECT_FALSE( true) << "gelss_double_parameters object: malloc error.";
		gelss_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gelss_double_parameters :: ~gelss_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelss_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelss_free();

}
/*  Test fixture class definition */
class dgelssest  : public  ::testing::Test {
public:
   gelss_double_parameters  *dgelss_obj;
   void SetUp();
   void TearDown () { delete dgelss_obj;}
};

void dgelssest::SetUp(){

    /* LAPACKE dgelss prototype */
    typedef int (*Fptr_NL_LAPACKE_dgelss) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, double* a, lapack_int lda, double* b, lapack_int ldb, double* s, double rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_dgelss dgelss;

    dgelss_obj = new gelss_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].isgn); 
	
	idx = Circular_Increment_Index(idx);

    dgelss_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgelss_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgelss_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgelss_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgelss = (Fptr_NL_LAPACKE_dgelss)dlsym(dgelss_obj->hModule, "LAPACKE_dgelss");
    ASSERT_TRUE(dgelss != NULL) << "failed to get the Netlib LAPACKE_dgelss symbol";
    

    dgelss_obj->inforef = dgelss( dgelss_obj->matrix_layout, dgelss_obj->m, dgelss_obj->n,dgelss_obj->nrhs, 
	dgelss_obj->Aref, dgelss_obj->lda, dgelss_obj->bref, dgelss_obj->ldb, dgelss_obj->sref, dgelss_obj->rcond, &dgelss_obj->rank);

    /* Compute libflame's Lapacke o/p  */
    dgelss_obj->info = LAPACKE_dgelss( dgelss_obj->matrix_layout, dgelss_obj->m,dgelss_obj->n,dgelss_obj->nrhs, 
	dgelss_obj->A, dgelss_obj->lda, dgelss_obj->b, dgelss_obj->ldb, dgelss_obj->s, dgelss_obj->rcond, &dgelss_obj->rank);
	
	if( dgelss_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgelss is wrong\n", dgelss_obj->info );
    }
    if( dgelss_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgelss is wrong\n", 
        dgelss_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    //dgelss_obj->diff =  computeDiff_d( dgelss_obj->bufsize, dgelss_obj->A, dgelss_obj->Aref );
	//dgelss_obj->diff_b =  computeDiff_d( dgelss_obj->bufsize_b, dgelss_obj->b, dgelss_obj->bref );
	dgelss_obj->diff_s =  computeDiff_d( dgelss_obj->bufsize_s, dgelss_obj->s, dgelss_obj->sref );

}

TEST_F(dgelssest, dgelss1) {
    //EXPECT_NEAR(0.0, dgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, dgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelssest, dgelss2) {
    //EXPECT_NEAR(0.0, dgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, dgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelssest, dgelss3) {
    //EXPECT_NEAR(0.0, dgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, dgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelssest, dgelss4) {
    //EXPECT_NEAR(0.0, dgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, dgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class gelss_scomplex_parameters{

   public:
	int bufsize, bufsize_b, bufsize_s;
	void *hModule, *dModule;
	float diff, diff_b, diff_s;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nrhs, rank;
	lapack_complex_float* A;
	float rcond;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_float* b;
	float *s, *sref;
	lapack_complex_float *Aref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gelss_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, float rcond);
      ~gelss_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
gelss_scomplex_parameters:: gelss_scomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, float rcond_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelss scomplex: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = fla_max(m,n);
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
	
	bufsize_s = fla_min(m,n);
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&s, &sref, bufsize_s);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(s==NULL) || (sref==NULL)){
		EXPECT_FALSE( true) << "gelss_scomplex_parameters object: malloc error.";
		gelss_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
gelss_scomplex_parameters :: ~gelss_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelss_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelss_free();

}
/*  Test fixture class definition */
class cgelssest  : public  ::testing::Test {
public:
   gelss_scomplex_parameters  *cgelss_obj;
   void SetUp();
   void TearDown () { delete cgelss_obj;}
};

void cgelssest::SetUp(){

    /* LAPACKE cgelss prototype */
    typedef int (*Fptr_NL_LAPACKE_cgelss) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, float* s, float rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_cgelss cgelss;

    cgelss_obj = new gelss_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].isgn); 
	
	idx = Circular_Increment_Index(idx);

    cgelss_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgelss_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgelss_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgelss_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgelss = (Fptr_NL_LAPACKE_cgelss)dlsym(cgelss_obj->hModule, "LAPACKE_cgelss");
    ASSERT_TRUE(cgelss != NULL) << "failed to get the Netlib LAPACKE_cgelss symbol";
    

    cgelss_obj->inforef = cgelss( cgelss_obj->matrix_layout, cgelss_obj->m, cgelss_obj->n,cgelss_obj->nrhs, 
	cgelss_obj->Aref, cgelss_obj->lda, cgelss_obj->bref, cgelss_obj->ldb, cgelss_obj->sref, cgelss_obj->rcond, &cgelss_obj->rank);

    /* Compute libflame's Lapacke o/p  */
    cgelss_obj->info = LAPACKE_cgelss( cgelss_obj->matrix_layout, cgelss_obj->m,cgelss_obj->n,cgelss_obj->nrhs, 
	cgelss_obj->A, cgelss_obj->lda, cgelss_obj->b, cgelss_obj->ldb, cgelss_obj->s, cgelss_obj->rcond, &cgelss_obj->rank);
	
	if( cgelss_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgelss is wrong\n", cgelss_obj->info );
    }
    if( cgelss_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgelss is wrong\n", 
        cgelss_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    //cgelss_obj->diff =  computeDiff_c( cgelss_obj->bufsize, cgelss_obj->A, cgelss_obj->Aref );
	//cgelss_obj->diff_b =  computeDiff_c( cgelss_obj->bufsize_b, cgelss_obj->b, cgelss_obj->bref );
	cgelss_obj->diff_s =  computeDiff_s( cgelss_obj->bufsize_s, cgelss_obj->s, cgelss_obj->sref );

}

TEST_F(cgelssest, cgelss1) {
    //EXPECT_NEAR(0.0, cgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, cgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelssest, cgelss2) {
    //EXPECT_NEAR(0.0, cgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, cgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelssest, cgelss3) {
    //EXPECT_NEAR(0.0, cgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, cgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelssest, cgelss4) {
    //EXPECT_NEAR(0.0, cgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, cgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class gelss_dcomplex_parameters{

   public:
	int bufsize, bufsize_b, bufsize_s;
	void *hModule, *dModule;
	double diff, diff_b, diff_s;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nrhs, rank;
	lapack_complex_double* A;
	double rcond;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_double* b;
	double *s, *sref;
	lapack_complex_double *Aref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gelss_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, double rcond);
      ~gelss_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
gelss_dcomplex_parameters:: gelss_dcomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, double rcond_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelss dcomplex: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = fla_max(m,n);
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
	
	bufsize_s = fla_min(m,n);
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&s, &sref, bufsize_s);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(s==NULL) || (sref==NULL)){
		EXPECT_FALSE( true) << "gelss_dcomplex_parameters object: malloc error.";
		gelss_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
gelss_dcomplex_parameters :: ~gelss_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelss_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelss_free();

}
/*  Test fixture class definition */
class zgelssest  : public  ::testing::Test {
public:
   gelss_dcomplex_parameters  *zgelss_obj;
   void SetUp();
   void TearDown () { delete zgelss_obj;}
};

void zgelssest::SetUp(){

    /* LAPACKE zgelss prototype */
    typedef int (*Fptr_NL_LAPACKE_zgelss) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, double* s, double rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_zgelss zgelss;

    zgelss_obj = new gelss_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].isgn); 
	
	idx = Circular_Increment_Index(idx);

    zgelss_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgelss_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgelss_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgelss_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgelss = (Fptr_NL_LAPACKE_zgelss)dlsym(zgelss_obj->hModule, "LAPACKE_zgelss");
    ASSERT_TRUE(zgelss != NULL) << "failed to get the Netlib LAPACKE_zgelss symbol";
    

    zgelss_obj->inforef = zgelss( zgelss_obj->matrix_layout, zgelss_obj->m, zgelss_obj->n,zgelss_obj->nrhs, 
	zgelss_obj->Aref, zgelss_obj->lda, zgelss_obj->bref, zgelss_obj->ldb, zgelss_obj->sref, zgelss_obj->rcond, &zgelss_obj->rank);

    /* Compute libflame's Lapacke o/p  */
    zgelss_obj->info = LAPACKE_zgelss( zgelss_obj->matrix_layout, zgelss_obj->m,zgelss_obj->n,zgelss_obj->nrhs, 
	zgelss_obj->A, zgelss_obj->lda, zgelss_obj->b, zgelss_obj->ldb, zgelss_obj->s, zgelss_obj->rcond, &zgelss_obj->rank);
	
	if( zgelss_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgelss is wrong\n", zgelss_obj->info );
    }
    if( zgelss_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgelss is wrong\n", 
        zgelss_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    //zgelss_obj->diff =  computeDiff_z( zgelss_obj->bufsize, zgelss_obj->A, zgelss_obj->Aref );
	//zgelss_obj->diff_b =  computeDiff_z( zgelss_obj->bufsize_b, zgelss_obj->b, zgelss_obj->bref );
	zgelss_obj->diff_s =  computeDiff_d( zgelss_obj->bufsize_s, zgelss_obj->s, zgelss_obj->sref );

}

TEST_F(zgelssest, zgelss1) {
    //EXPECT_NEAR(0.0, zgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, zgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelssest, zgelss2) {
    //EXPECT_NEAR(0.0, zgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, zgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelssest, zgelss3) {
    //EXPECT_NEAR(0.0, zgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, zgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelssest, zgelss4) {
    //EXPECT_NEAR(0.0, zgelss_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, zgelss_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgelss_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}