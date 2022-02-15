#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define gelsd_free() \
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
class gelsd_float_parameters{

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
      gelsd_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, float rcond);
      ~gelsd_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gelsd_float_parameters:: gelsd_float_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, float rcond_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelsd float: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
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
	
	bufsize_s = min(m,n);
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&s, &sref, bufsize_s);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(s==NULL) || (sref==NULL)){
		EXPECT_FALSE( true) << "gelsd_float_parameters object: malloc error.";
		gelsd_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gelsd_float_parameters :: ~gelsd_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelsd_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelsd_free();

}
/*  Test fixture class definition */
class sgelsdest  : public  ::testing::Test {
public:
   gelsd_float_parameters  *sgelsd_obj;
   void SetUp();
   void TearDown () { delete sgelsd_obj;}
};

void sgelsdest::SetUp(){

    /* LAPACKE sgelsd prototype */
    typedef int (*Fptr_NL_LAPACKE_sgelsd) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, float* a, lapack_int lda, float* b, lapack_int ldb, float* s, float rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_sgelsd sgelsd;

    sgelsd_obj = new gelsd_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].isgn); 
	
	idx = Circular_Increment_Index(idx);

    sgelsd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgelsd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgelsd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgelsd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgelsd = (Fptr_NL_LAPACKE_sgelsd)dlsym(sgelsd_obj->hModule, "LAPACKE_sgelsd");
    ASSERT_TRUE(sgelsd != NULL) << "failed to get the Netlib LAPACKE_sgelsd symbol";
    

    sgelsd_obj->inforef = sgelsd( sgelsd_obj->matrix_layout, sgelsd_obj->m, sgelsd_obj->n,sgelsd_obj->nrhs, 
	sgelsd_obj->Aref, sgelsd_obj->lda, sgelsd_obj->bref, sgelsd_obj->ldb, sgelsd_obj->sref, sgelsd_obj->rcond, &sgelsd_obj->rank);

    /* Compute libflame's Lapacke o/p  */
    sgelsd_obj->info = LAPACKE_sgelsd( sgelsd_obj->matrix_layout, sgelsd_obj->m,sgelsd_obj->n,sgelsd_obj->nrhs, 
	sgelsd_obj->A, sgelsd_obj->lda, sgelsd_obj->b, sgelsd_obj->ldb, sgelsd_obj->s, sgelsd_obj->rcond, &sgelsd_obj->rank);
	
	if( sgelsd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgelsd is wrong\n", sgelsd_obj->info );
    }
    if( sgelsd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgelsd is wrong\n", 
        sgelsd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    //sgelsd_obj->diff =  computeDiff_s( sgelsd_obj->bufsize, sgelsd_obj->A, sgelsd_obj->Aref );
	//sgelsd_obj->diff_b =  computeDiff_s( sgelsd_obj->bufsize_b, sgelsd_obj->b, sgelsd_obj->bref );
	sgelsd_obj->diff_s =  computeDiff_s( sgelsd_obj->bufsize_s, sgelsd_obj->s, sgelsd_obj->sref );

}

TEST_F(sgelsdest, sgelsd1) {
    //EXPECT_NEAR(0.0, sgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, sgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelsdest, sgelsd2) {
    //EXPECT_NEAR(0.0, sgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, sgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelsdest, sgelsd3) {
    //EXPECT_NEAR(0.0, sgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, sgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(sgelsdest, sgelsd4) {
    //EXPECT_NEAR(0.0, sgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, sgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, sgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class gelsd_double_parameters{

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
      gelsd_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, double rcond);
      ~gelsd_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gelsd_double_parameters:: gelsd_double_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, double rcond_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelsd double: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
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
	
	bufsize_s = min(m,n);
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&s, &sref, bufsize_s);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(s==NULL) || (sref==NULL)){
		EXPECT_FALSE( true) << "gelsd_double_parameters object: malloc error.";
		gelsd_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gelsd_double_parameters :: ~gelsd_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelsd_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelsd_free();

}
/*  Test fixture class definition */
class dgelsdest  : public  ::testing::Test {
public:
   gelsd_double_parameters  *dgelsd_obj;
   void SetUp();
   void TearDown () { delete dgelsd_obj;}
};

void dgelsdest::SetUp(){

    /* LAPACKE dgelsd prototype */
    typedef int (*Fptr_NL_LAPACKE_dgelsd) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, double* a, lapack_int lda, double* b, lapack_int ldb, double* s, double rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_dgelsd dgelsd;

    dgelsd_obj = new gelsd_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].isgn); 
	
	idx = Circular_Increment_Index(idx);

    dgelsd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgelsd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgelsd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgelsd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgelsd = (Fptr_NL_LAPACKE_dgelsd)dlsym(dgelsd_obj->hModule, "LAPACKE_dgelsd");
    ASSERT_TRUE(dgelsd != NULL) << "failed to get the Netlib LAPACKE_dgelsd symbol";
    

    dgelsd_obj->inforef = dgelsd( dgelsd_obj->matrix_layout, dgelsd_obj->m, dgelsd_obj->n,dgelsd_obj->nrhs, 
	dgelsd_obj->Aref, dgelsd_obj->lda, dgelsd_obj->bref, dgelsd_obj->ldb, dgelsd_obj->sref, dgelsd_obj->rcond, &dgelsd_obj->rank);

    /* Compute libflame's Lapacke o/p  */
    dgelsd_obj->info = LAPACKE_dgelsd( dgelsd_obj->matrix_layout, dgelsd_obj->m,dgelsd_obj->n,dgelsd_obj->nrhs, 
	dgelsd_obj->A, dgelsd_obj->lda, dgelsd_obj->b, dgelsd_obj->ldb, dgelsd_obj->s, dgelsd_obj->rcond, &dgelsd_obj->rank);
	
	if( dgelsd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgelsd is wrong\n", dgelsd_obj->info );
    }
    if( dgelsd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgelsd is wrong\n", 
        dgelsd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    //dgelsd_obj->diff =  computeDiff_d( dgelsd_obj->bufsize, dgelsd_obj->A, dgelsd_obj->Aref );
	//dgelsd_obj->diff_b =  computeDiff_d( dgelsd_obj->bufsize_b, dgelsd_obj->b, dgelsd_obj->bref );
	dgelsd_obj->diff_s =  computeDiff_d( dgelsd_obj->bufsize_s, dgelsd_obj->s, dgelsd_obj->sref );

}

TEST_F(dgelsdest, dgelsd1) {
    //EXPECT_NEAR(0.0, dgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, dgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelsdest, dgelsd2) {
    //EXPECT_NEAR(0.0, dgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, dgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelsdest, dgelsd3) {
    //EXPECT_NEAR(0.0, dgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, dgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dgelsdest, dgelsd4) {
    //EXPECT_NEAR(0.0, dgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, dgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class gelsd_scomplex_parameters{

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
      gelsd_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, float rcond);
      ~gelsd_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
gelsd_scomplex_parameters:: gelsd_scomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, float rcond_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelsd scomplex: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
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
	
	bufsize_s = min(m,n);
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&s, &sref, bufsize_s);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(s==NULL) || (sref==NULL)){
		EXPECT_FALSE( true) << "gelsd_scomplex_parameters object: malloc error.";
		gelsd_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
gelsd_scomplex_parameters :: ~gelsd_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelsd_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelsd_free();

}
/*  Test fixture class definition */
class cgelsdest  : public  ::testing::Test {
public:
   gelsd_scomplex_parameters  *cgelsd_obj;
   void SetUp();
   void TearDown () { delete cgelsd_obj;}
};

void cgelsdest::SetUp(){

    /* LAPACKE cgelsd prototype */
    typedef int (*Fptr_NL_LAPACKE_cgelsd) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, float* s, float rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_cgelsd cgelsd;

    cgelsd_obj = new gelsd_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].isgn); 
	
	idx = Circular_Increment_Index(idx);

    cgelsd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgelsd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgelsd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgelsd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgelsd = (Fptr_NL_LAPACKE_cgelsd)dlsym(cgelsd_obj->hModule, "LAPACKE_cgelsd");
    ASSERT_TRUE(cgelsd != NULL) << "failed to get the Netlib LAPACKE_cgelsd symbol";
    

    cgelsd_obj->inforef = cgelsd( cgelsd_obj->matrix_layout, cgelsd_obj->m, cgelsd_obj->n,cgelsd_obj->nrhs, 
	cgelsd_obj->Aref, cgelsd_obj->lda, cgelsd_obj->bref, cgelsd_obj->ldb, cgelsd_obj->sref, cgelsd_obj->rcond, &cgelsd_obj->rank);

    /* Compute libflame's Lapacke o/p  */
    cgelsd_obj->info = LAPACKE_cgelsd( cgelsd_obj->matrix_layout, cgelsd_obj->m,cgelsd_obj->n,cgelsd_obj->nrhs, 
	cgelsd_obj->A, cgelsd_obj->lda, cgelsd_obj->b, cgelsd_obj->ldb, cgelsd_obj->s, cgelsd_obj->rcond, &cgelsd_obj->rank);
	
	if( cgelsd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgelsd is wrong\n", cgelsd_obj->info );
    }
    if( cgelsd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgelsd is wrong\n", 
        cgelsd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    //cgelsd_obj->diff =  computeDiff_c( cgelsd_obj->bufsize, cgelsd_obj->A, cgelsd_obj->Aref );
	//cgelsd_obj->diff_b =  computeDiff_c( cgelsd_obj->bufsize_b, cgelsd_obj->b, cgelsd_obj->bref );
	cgelsd_obj->diff_s =  computeDiff_s( cgelsd_obj->bufsize_s, cgelsd_obj->s, cgelsd_obj->sref );

}

TEST_F(cgelsdest, cgelsd1) {
    //EXPECT_NEAR(0.0, cgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, cgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelsdest, cgelsd2) {
    //EXPECT_NEAR(0.0, cgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, cgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelsdest, cgelsd3) {
    //EXPECT_NEAR(0.0, cgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, cgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(cgelsdest, cgelsd4) {
    //EXPECT_NEAR(0.0, cgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, cgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, cgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class gelsd_dcomplex_parameters{

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
      gelsd_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, double rcond);
      ~gelsd_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
gelsd_dcomplex_parameters:: gelsd_dcomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, double rcond_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelsd dcomplex: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
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
	
	bufsize_s = min(m,n);
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&s, &sref, bufsize_s);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(s==NULL) || (sref==NULL)){
		EXPECT_FALSE( true) << "gelsd_dcomplex_parameters object: malloc error.";
		gelsd_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
gelsd_dcomplex_parameters :: ~gelsd_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelsd_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelsd_free();

}
/*  Test fixture class definition */
class zgelsdest  : public  ::testing::Test {
public:
   gelsd_dcomplex_parameters  *zgelsd_obj;
   void SetUp();
   void TearDown () { delete zgelsd_obj;}
};

void zgelsdest::SetUp(){

    /* LAPACKE zgelsd prototype */
    typedef int (*Fptr_NL_LAPACKE_zgelsd) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, double* s, double rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_zgelsd zgelsd;

    zgelsd_obj = new gelsd_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].nrhs,
						   eig_paramslist[idx].isgn); 
	
	idx = Circular_Increment_Index(idx);

    zgelsd_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgelsd_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgelsd_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgelsd_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgelsd = (Fptr_NL_LAPACKE_zgelsd)dlsym(zgelsd_obj->hModule, "LAPACKE_zgelsd");
    ASSERT_TRUE(zgelsd != NULL) << "failed to get the Netlib LAPACKE_zgelsd symbol";
    

    zgelsd_obj->inforef = zgelsd( zgelsd_obj->matrix_layout, zgelsd_obj->m, zgelsd_obj->n,zgelsd_obj->nrhs, 
	zgelsd_obj->Aref, zgelsd_obj->lda, zgelsd_obj->bref, zgelsd_obj->ldb, zgelsd_obj->sref, zgelsd_obj->rcond, &zgelsd_obj->rank);

    /* Compute libflame's Lapacke o/p  */
    zgelsd_obj->info = LAPACKE_zgelsd( zgelsd_obj->matrix_layout, zgelsd_obj->m,zgelsd_obj->n,zgelsd_obj->nrhs, 
	zgelsd_obj->A, zgelsd_obj->lda, zgelsd_obj->b, zgelsd_obj->ldb, zgelsd_obj->s, zgelsd_obj->rcond, &zgelsd_obj->rank);
	
	if( zgelsd_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgelsd is wrong\n", zgelsd_obj->info );
    }
    if( zgelsd_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgelsd is wrong\n", 
        zgelsd_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    //zgelsd_obj->diff =  computeDiff_z( zgelsd_obj->bufsize, zgelsd_obj->A, zgelsd_obj->Aref );
	//zgelsd_obj->diff_b =  computeDiff_z( zgelsd_obj->bufsize_b, zgelsd_obj->b, zgelsd_obj->bref );
	zgelsd_obj->diff_s =  computeDiff_d( zgelsd_obj->bufsize_s, zgelsd_obj->s, zgelsd_obj->sref );

}

TEST_F(zgelsdest, zgelsd1) {
    //EXPECT_NEAR(0.0, zgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, zgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelsdest, zgelsd2) {
    //EXPECT_NEAR(0.0, zgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, zgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelsdest, zgelsd3) {
    //EXPECT_NEAR(0.0, zgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
//	EXPECT_NEAR(0.0, zgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zgelsdest, zgelsd4) {
    //EXPECT_NEAR(0.0, zgelsd_obj->diff, LAPACKE_GTEST_THRESHOLD);
	//EXPECT_NEAR(0.0, zgelsd_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, zgelsd_obj->diff_s, LAPACKE_GTEST_THRESHOLD);
}
