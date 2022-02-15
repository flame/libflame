#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"


#define gelsy_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)  free(b);\
if (bref!=NULL) free(bref); \
if (jpvt!=NULL)    free(jpvt); \
if (jpvtref!=NULL) free(jpvtref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class gelsy_float_parameters{

   public:
	int bufsize, bufsize_b, bufsize_jpvt;
	void *hModule, *dModule;
	float diff, diff_b, diff_s;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nrhs, rank;
	float* A, rcond, eps;
	lapack_int lda, ldb;
	lapack_int gelsy_threshold;
	/*Output Parameter*/
	float* b;
	lapack_int *jpvt;
	float *Aref, *bref;
	lapack_int *jpvtref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gelsy_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, float rcond, lapack_int gelsy_threshold);
      ~gelsy_float_parameters ();

};

/* Constructor definition  float_common_parameters */
gelsy_float_parameters:: gelsy_float_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, float rcond_i, lapack_int gelsy_threshold_i)
{

	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	gelsy_threshold = gelsy_threshold_i;
	
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelsy float: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
	#endif
	
	/*getting rcond variable from eps */
	eps =  LAPACKE_slamch( 'E' );
	rcond = sqrt( eps ) - ( sqrt( eps )-eps ) / 2;
	
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
	
	bufsize_jpvt = n;
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_int_buffer_pair(&jpvt, &jpvtref, bufsize_jpvt);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(jpvt==NULL) || (jpvtref==NULL)){
		EXPECT_FALSE( true) << "gelsy_float_parameters object: malloc error.";
		gelsy_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_int_buffer_pair_rand( jpvt, jpvtref, bufsize_jpvt);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
gelsy_float_parameters :: ~gelsy_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelsy_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelsy_free();

}
/*  Test fixture class definition */
class sgelsyest  : public  ::testing::Test {
public:
   gelsy_float_parameters  *sgelsy_obj;
   void SetUp();
   void TearDown () { delete sgelsy_obj;}
};

void sgelsyest::SetUp(){

    /* LAPACKE sgelsy prototype */
    typedef int (*Fptr_NL_LAPACKE_sgelsy) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, float* a, lapack_int lda, float* b, lapack_int ldb, lapack_int* jpvt, float rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_sgelsy sgelsy;

    sgelsy_obj = new gelsy_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].isgn,
						   eig_paramslist[idx].threshold_value); 
	
	idx = Circular_Increment_Index(idx);

    sgelsy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sgelsy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sgelsy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sgelsy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sgelsy = (Fptr_NL_LAPACKE_sgelsy)dlsym(sgelsy_obj->hModule, "LAPACKE_sgelsy");
    ASSERT_TRUE(sgelsy != NULL) << "failed to get the Netlib LAPACKE_sgelsy symbol";
    

    sgelsy_obj->inforef = sgelsy( sgelsy_obj->matrix_layout, sgelsy_obj->m, sgelsy_obj->n,sgelsy_obj->nrhs, 
	sgelsy_obj->Aref, sgelsy_obj->lda, sgelsy_obj->bref, sgelsy_obj->ldb, sgelsy_obj->jpvtref, sgelsy_obj->rcond, &sgelsy_obj->rank);

    /* Compute libflame'jpvt Lapacke o/p  */
    sgelsy_obj->info = LAPACKE_sgelsy( sgelsy_obj->matrix_layout, sgelsy_obj->m,sgelsy_obj->n,sgelsy_obj->nrhs, 
	sgelsy_obj->A, sgelsy_obj->lda, sgelsy_obj->b, sgelsy_obj->ldb, sgelsy_obj->jpvt, sgelsy_obj->rcond, &sgelsy_obj->rank);
	
	if( sgelsy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sgelsy is wrong\n", sgelsy_obj->info );
    }
    if( sgelsy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sgelsy is wrong\n", 
        sgelsy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sgelsy_obj->diff =  computeDiff_s( sgelsy_obj->bufsize, sgelsy_obj->A, sgelsy_obj->Aref );
	sgelsy_obj->diff_b =  computeDiff_s( sgelsy_obj->bufsize_b, sgelsy_obj->b, sgelsy_obj->bref );
	sgelsy_obj->diff_s =  computeDiff_i( sgelsy_obj->bufsize_jpvt, sgelsy_obj->jpvt, sgelsy_obj->jpvtref );

}

TEST_F(sgelsyest, sgelsy1) {
    EXPECT_NEAR(0.0, sgelsy_obj->diff, sgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, sgelsy_obj->diff_b, sgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, sgelsy_obj->diff_s, sgelsy_obj->gelsy_threshold);
}

TEST_F(sgelsyest, sgelsy2) {
    EXPECT_NEAR(0.0, sgelsy_obj->diff, sgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, sgelsy_obj->diff_b, sgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, sgelsy_obj->diff_s, sgelsy_obj->gelsy_threshold);
}

TEST_F(sgelsyest, sgelsy3) {
    EXPECT_NEAR(0.0, sgelsy_obj->diff, sgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, sgelsy_obj->diff_b, sgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, sgelsy_obj->diff_s, sgelsy_obj->gelsy_threshold);
}

TEST_F(sgelsyest, sgelsy4) {
    EXPECT_NEAR(0.0, sgelsy_obj->diff, sgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, sgelsy_obj->diff_b, sgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, sgelsy_obj->diff_s, sgelsy_obj->gelsy_threshold);
}


/* Begin double_common_parameters  class definition */
class gelsy_double_parameters{

   public:
	int bufsize, bufsize_b, bufsize_jpvt;
	void *hModule, *dModule;
	double diff, diff_b, diff_s;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nrhs, rank;
	double* A, rcond;
	lapack_int lda, ldb;
	double eps;
	lapack_int gelsy_threshold;
	/*Output Parameter*/
	double* b;
	lapack_int *jpvt;
	double *Aref, *bref;
	lapack_int *jpvtref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gelsy_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, double rcond, lapack_int gelsy_threshold);
      ~gelsy_double_parameters ();

};

/* Constructor definition  double_common_parameters */
gelsy_double_parameters:: gelsy_double_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, double rcond_i, lapack_int gelsy_threshold_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	gelsy_threshold = gelsy_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelsy double: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
	#endif
	
	/*getting rcond variable from eps */
	eps =  LAPACKE_dlamch( 'E' );
	rcond = sqrt( eps ) - ( sqrt( eps )-eps ) / 2;
	
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
	
	bufsize_jpvt = n;
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_int_buffer_pair(&jpvt, &jpvtref, bufsize_jpvt);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(jpvt==NULL) || (jpvtref==NULL)){
		EXPECT_FALSE( true) << "gelsy_double_parameters object: malloc error.";
		gelsy_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_int_buffer_pair_rand( jpvt, jpvtref, bufsize_jpvt);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
gelsy_double_parameters :: ~gelsy_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelsy_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelsy_free();

}
/*  Test fixture class definition */
class dgelsyest  : public  ::testing::Test {
public:
   gelsy_double_parameters  *dgelsy_obj;
   void SetUp();
   void TearDown () { delete dgelsy_obj;}
};

void dgelsyest::SetUp(){

    /* LAPACKE dgelsy prototype */
    typedef int (*Fptr_NL_LAPACKE_dgelsy) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, double* a, lapack_int lda, double* b, lapack_int ldb, lapack_int* jpvt, double rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_dgelsy dgelsy;

    dgelsy_obj = new gelsy_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].isgn,
						   eig_paramslist[idx].threshold_value); 
	
	idx = Circular_Increment_Index(idx);

    dgelsy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dgelsy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dgelsy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dgelsy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dgelsy = (Fptr_NL_LAPACKE_dgelsy)dlsym(dgelsy_obj->hModule, "LAPACKE_dgelsy");
    ASSERT_TRUE(dgelsy != NULL) << "failed to get the Netlib LAPACKE_dgelsy symbol";
    

    dgelsy_obj->inforef = dgelsy( dgelsy_obj->matrix_layout, dgelsy_obj->m, dgelsy_obj->n,dgelsy_obj->nrhs, 
	dgelsy_obj->Aref, dgelsy_obj->lda, dgelsy_obj->bref, dgelsy_obj->ldb, dgelsy_obj->jpvtref, dgelsy_obj->rcond, &dgelsy_obj->rank);

    /* Compute libflame'jpvt Lapacke o/p  */
    dgelsy_obj->info = LAPACKE_dgelsy( dgelsy_obj->matrix_layout, dgelsy_obj->m,dgelsy_obj->n,dgelsy_obj->nrhs, 
	dgelsy_obj->A, dgelsy_obj->lda, dgelsy_obj->b, dgelsy_obj->ldb, dgelsy_obj->jpvt, dgelsy_obj->rcond, &dgelsy_obj->rank);
	
	if( dgelsy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dgelsy is wrong\n", dgelsy_obj->info );
    }
    if( dgelsy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dgelsy is wrong\n", 
        dgelsy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dgelsy_obj->diff =  computeDiff_d( dgelsy_obj->bufsize, dgelsy_obj->A, dgelsy_obj->Aref );
	dgelsy_obj->diff_b =  computeDiff_d( dgelsy_obj->bufsize_b, dgelsy_obj->b, dgelsy_obj->bref );
	dgelsy_obj->diff_s =  computeDiff_i( dgelsy_obj->bufsize_jpvt, dgelsy_obj->jpvt, dgelsy_obj->jpvtref );

}

TEST_F(dgelsyest, dgelsy1) {
    EXPECT_NEAR(0.0, dgelsy_obj->diff, dgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, dgelsy_obj->diff_b, dgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, dgelsy_obj->diff_s, dgelsy_obj->gelsy_threshold);
}

TEST_F(dgelsyest, dgelsy2) {
    EXPECT_NEAR(0.0, dgelsy_obj->diff, dgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, dgelsy_obj->diff_b, dgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, dgelsy_obj->diff_s, dgelsy_obj->gelsy_threshold);
}

TEST_F(dgelsyest, dgelsy3) {
    EXPECT_NEAR(0.0, dgelsy_obj->diff, dgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, dgelsy_obj->diff_b, dgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, dgelsy_obj->diff_s, dgelsy_obj->gelsy_threshold);
}

TEST_F(dgelsyest, dgelsy4) {
    EXPECT_NEAR(0.0, dgelsy_obj->diff, dgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, dgelsy_obj->diff_b, dgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, dgelsy_obj->diff_s, dgelsy_obj->gelsy_threshold);
}


/* Begin scomplex_common_parameters  class definition */
class gelsy_scomplex_parameters{

   public:
	int bufsize, bufsize_b, bufsize_jpvt;
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
	lapack_int gelsy_threshold;
	/*Output Parameter*/
	lapack_complex_float* b;
	lapack_int *jpvt, *jpvtref;
	lapack_complex_float *Aref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gelsy_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, float rcond, lapack_int gelsy_threshold);
      ~gelsy_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
gelsy_scomplex_parameters:: gelsy_scomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, float rcond_i, lapack_int gelsy_threshold_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	gelsy_threshold = gelsy_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelsy scomplex: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
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
	
	bufsize_jpvt = n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_int_buffer_pair(&jpvt, &jpvtref, bufsize_jpvt);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(jpvt==NULL) || (jpvtref==NULL)){
		EXPECT_FALSE( true) << "gelsy_scomplex_parameters object: malloc error.";
		gelsy_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_int_buffer_pair_rand( jpvt, jpvtref, bufsize_jpvt);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
gelsy_scomplex_parameters :: ~gelsy_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelsy_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelsy_free();

}
/*  Test fixture class definition */
class cgelsyest  : public  ::testing::Test {
public:
   gelsy_scomplex_parameters  *cgelsy_obj;
   void SetUp();
   void TearDown () { delete cgelsy_obj;}
};

void cgelsyest::SetUp(){

    /* LAPACKE cgelsy prototype */
    typedef int (*Fptr_NL_LAPACKE_cgelsy) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, lapack_int* jpvt, float rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_cgelsy cgelsy;

    cgelsy_obj = new gelsy_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].isgn,
						   eig_paramslist[idx].threshold_value); 
	
	idx = Circular_Increment_Index(idx);

    cgelsy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cgelsy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cgelsy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cgelsy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cgelsy = (Fptr_NL_LAPACKE_cgelsy)dlsym(cgelsy_obj->hModule, "LAPACKE_cgelsy");
    ASSERT_TRUE(cgelsy != NULL) << "failed to get the Netlib LAPACKE_cgelsy symbol";
    

    cgelsy_obj->inforef = cgelsy( cgelsy_obj->matrix_layout, cgelsy_obj->m, cgelsy_obj->n,cgelsy_obj->nrhs, 
	cgelsy_obj->Aref, cgelsy_obj->lda, cgelsy_obj->bref, cgelsy_obj->ldb, cgelsy_obj->jpvtref, cgelsy_obj->rcond, &cgelsy_obj->rank);

    /* Compute libflame'jpvt Lapacke o/p  */
    cgelsy_obj->info = LAPACKE_cgelsy( cgelsy_obj->matrix_layout, cgelsy_obj->m,cgelsy_obj->n,cgelsy_obj->nrhs, 
	cgelsy_obj->A, cgelsy_obj->lda, cgelsy_obj->b, cgelsy_obj->ldb, cgelsy_obj->jpvt, cgelsy_obj->rcond, &cgelsy_obj->rank);
	
	if( cgelsy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cgelsy is wrong\n", cgelsy_obj->info );
    }
    if( cgelsy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cgelsy is wrong\n", 
        cgelsy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cgelsy_obj->diff =  computeDiff_c( cgelsy_obj->bufsize, cgelsy_obj->A, cgelsy_obj->Aref );
	cgelsy_obj->diff_b =  computeDiff_c( cgelsy_obj->bufsize_b, cgelsy_obj->b, cgelsy_obj->bref );
	cgelsy_obj->diff_s =  computeDiff_i( cgelsy_obj->bufsize_jpvt, cgelsy_obj->jpvt, cgelsy_obj->jpvtref );

}

TEST_F(cgelsyest, cgelsy1) {
    EXPECT_NEAR(0.0, cgelsy_obj->diff, cgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, cgelsy_obj->diff_b, cgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, cgelsy_obj->diff_s, cgelsy_obj->gelsy_threshold);
}

TEST_F(cgelsyest, cgelsy2) {
    EXPECT_NEAR(0.0, cgelsy_obj->diff, cgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, cgelsy_obj->diff_b, cgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, cgelsy_obj->diff_s, cgelsy_obj->gelsy_threshold);
}

TEST_F(cgelsyest, cgelsy3) {
    EXPECT_NEAR(0.0, cgelsy_obj->diff, cgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, cgelsy_obj->diff_b, cgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, cgelsy_obj->diff_s, cgelsy_obj->gelsy_threshold);
}

TEST_F(cgelsyest, cgelsy4) {
    EXPECT_NEAR(0.0, cgelsy_obj->diff, cgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, cgelsy_obj->diff_b, cgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, cgelsy_obj->diff_s, cgelsy_obj->gelsy_threshold);
}

/* Begin dcomplex_common_parameters  class definition */
class gelsy_dcomplex_parameters{

   public:
	int bufsize, bufsize_b, bufsize_jpvt;
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
	lapack_int gelsy_threshold;
	/*Output Parameter*/
	lapack_complex_double* b;
	lapack_int *jpvt, *jpvtref;
	lapack_complex_double *Aref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      gelsy_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nrhs, double rcond, lapack_int gelsy_threshold);
      ~gelsy_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
gelsy_dcomplex_parameters:: gelsy_dcomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nrhs_i, double rcond_i, lapack_int gelsy_threshold_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nrhs = nrhs_i;
	rcond = rcond_i;
	gelsy_threshold = gelsy_threshold_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n gelsy dcomplex: m: %d, n: %d , nrhs: %d, rcond:%d \n",  m, n, nrhs, rcond);
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
	
	bufsize_jpvt = n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_int_buffer_pair(&jpvt, &jpvtref, bufsize_jpvt);
	
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(jpvt==NULL) || (jpvtref==NULL)){
		EXPECT_FALSE( true) << "gelsy_dcomplex_parameters object: malloc error.";
		gelsy_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	lapacke_gtest_init_int_buffer_pair_rand( jpvt, jpvtref, bufsize_jpvt);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
gelsy_dcomplex_parameters :: ~gelsy_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" gelsy_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   gelsy_free();

}
/*  Test fixture class definition */
class zgelsyest  : public  ::testing::Test {
public:
   gelsy_dcomplex_parameters  *zgelsy_obj;
   void SetUp();
   void TearDown () { delete zgelsy_obj;}
};

void zgelsyest::SetUp(){

    /* LAPACKE zgelsy prototype */
    typedef int (*Fptr_NL_LAPACKE_zgelsy) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int nrhs, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, lapack_int* jpvt, double rcond, lapack_int* rank );

    Fptr_NL_LAPACKE_zgelsy zgelsy;

    zgelsy_obj = new gelsy_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].isgn,
						   eig_paramslist[idx].threshold_value); 
	
	idx = Circular_Increment_Index(idx);

    zgelsy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zgelsy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zgelsy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zgelsy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zgelsy = (Fptr_NL_LAPACKE_zgelsy)dlsym(zgelsy_obj->hModule, "LAPACKE_zgelsy");
    ASSERT_TRUE(zgelsy != NULL) << "failed to get the Netlib LAPACKE_zgelsy symbol";
    

    zgelsy_obj->inforef = zgelsy( zgelsy_obj->matrix_layout, zgelsy_obj->m, zgelsy_obj->n,zgelsy_obj->nrhs, 
	zgelsy_obj->Aref, zgelsy_obj->lda, zgelsy_obj->bref, zgelsy_obj->ldb, zgelsy_obj->jpvtref, zgelsy_obj->rcond, &zgelsy_obj->rank);

    /* Compute libflame'jpvt Lapacke o/p  */
    zgelsy_obj->info = LAPACKE_zgelsy( zgelsy_obj->matrix_layout, zgelsy_obj->m,zgelsy_obj->n,zgelsy_obj->nrhs, 
	zgelsy_obj->A, zgelsy_obj->lda, zgelsy_obj->b, zgelsy_obj->ldb, zgelsy_obj->jpvt, zgelsy_obj->rcond, &zgelsy_obj->rank);
	
	if( zgelsy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zgelsy is wrong\n", zgelsy_obj->info );
    }
    if( zgelsy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zgelsy is wrong\n", 
        zgelsy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zgelsy_obj->diff =  computeDiff_z( zgelsy_obj->bufsize, zgelsy_obj->A, zgelsy_obj->Aref );
	zgelsy_obj->diff_b =  computeDiff_z( zgelsy_obj->bufsize_b, zgelsy_obj->b, zgelsy_obj->bref );
	zgelsy_obj->diff_s =  computeDiff_i( zgelsy_obj->bufsize_jpvt, zgelsy_obj->jpvt, zgelsy_obj->jpvtref );

}

TEST_F(zgelsyest, zgelsy1) {
    EXPECT_NEAR(0.0, zgelsy_obj->diff, zgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, zgelsy_obj->diff_b, zgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, zgelsy_obj->diff_s, zgelsy_obj->gelsy_threshold);
}

TEST_F(zgelsyest, zgelsy2) {
    EXPECT_NEAR(0.0, zgelsy_obj->diff, zgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, zgelsy_obj->diff_b, zgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, zgelsy_obj->diff_s, zgelsy_obj->gelsy_threshold);
}

TEST_F(zgelsyest, zgelsy3) {
    EXPECT_NEAR(0.0, zgelsy_obj->diff, zgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, zgelsy_obj->diff_b, zgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, zgelsy_obj->diff_s, zgelsy_obj->gelsy_threshold);
}

TEST_F(zgelsyest, zgelsy4) {
    EXPECT_NEAR(0.0, zgelsy_obj->diff, zgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, zgelsy_obj->diff_b, zgelsy_obj->gelsy_threshold);
	EXPECT_NEAR(0.0, zgelsy_obj->diff_s, zgelsy_obj->gelsy_threshold);
}
