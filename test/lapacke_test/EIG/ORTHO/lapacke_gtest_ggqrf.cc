#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"


#define ggqrf_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)    free(b); \
if (bref!=NULL) free(bref); \
if (taua!=NULL)  free(taua);\
if (tauaref!=NULL) free(tauaref); \
if (taub!=NULL)  free(taub);\
if (taubref!=NULL) free(taubref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class ggqrf_float_parameters{

   public:
	int bufsize, bufsize_b;
	void *hModule, *dModule;
	float diff, diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m, p;
	float* A, *b;
	lapack_int lda, ldb;	
	/*Output Parameter*/
	float* taua, *taub;
	float *Aref, *bref;
	float *tauaref, *taubref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      ggqrf_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggqrf_float_parameters ();

};

/* Constructor definition  float_common_parameters */
ggqrf_float_parameters:: ggqrf_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n ggqrf float:  m: %d, n: %d, p: %d, \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else if (matrix_layout == LAPACK_ROW_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&taua, &tauaref, min(n,m));
	lapacke_gtest_alloc_float_buffer_pair(&taub, &taubref, min(n,p));
		if ((A==NULL) || (Aref==NULL) ||\
		(b == NULL) || (bref == NULL) ||
		(taua == NULL) || (tauaref == NULL) ||
		(taub == NULL) || (taubref == NULL)){
		EXPECT_FALSE( true) << "ggqrf_float_parameters object: malloc error.";
		ggqrf_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
ggqrf_float_parameters :: ~ggqrf_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggqrf_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggqrf_free();

}
/*  Test fixture class definition */
class sggqrfest  : public  ::testing::Test {
public:
   ggqrf_float_parameters  *sggqrf_obj;
   void SetUp();
   void TearDown () { delete sggqrf_obj;}
};

void sggqrfest::SetUp(){

    /* LAPACKE sggqrf prototype */
    typedef int (*Fptr_NL_LAPACKE_sggqrf) (int matrix_layout, lapack_int n, lapack_int m, lapack_int p,\
											float* a, lapack_int lda, float* taua, float* b, lapack_int ldb, float* taub);

    Fptr_NL_LAPACKE_sggqrf sggqrf;

    sggqrf_obj = new ggqrf_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n); 
	
	idx = Circular_Increment_Index(idx);

    sggqrf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sggqrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggqrf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggqrf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sggqrf = (Fptr_NL_LAPACKE_sggqrf)dlsym(sggqrf_obj->hModule, "LAPACKE_sggqrf");
    ASSERT_TRUE(sggqrf != NULL) << "failed to get the Netlib LAPACKE_sggqrf symbol";
    

    sggqrf_obj->inforef = sggqrf( sggqrf_obj->matrix_layout, sggqrf_obj->n, sggqrf_obj->m,
								sggqrf_obj->p, sggqrf_obj->Aref, sggqrf_obj->lda, sggqrf_obj->tauaref,
								sggqrf_obj->bref, sggqrf_obj->ldb, sggqrf_obj->taubref);

    /* Compute libflame's Lapacke o/p  */
    sggqrf_obj->info = LAPACKE_sggqrf( sggqrf_obj->matrix_layout, sggqrf_obj->n, sggqrf_obj->m,
								sggqrf_obj->p, sggqrf_obj->A, sggqrf_obj->lda, sggqrf_obj->taua,
								sggqrf_obj->b, sggqrf_obj->ldb, sggqrf_obj->taub);
	
	if( sggqrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sggqrf is wrong\n", sggqrf_obj->info );
    }
    if( sggqrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sggqrf is wrong\n", 
        sggqrf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sggqrf_obj->diff =  computeDiff_s( sggqrf_obj->bufsize, 
                sggqrf_obj->A, sggqrf_obj->Aref );
	sggqrf_obj->diff_b =  computeDiff_s( sggqrf_obj->bufsize_b, 
                sggqrf_obj->b, sggqrf_obj->bref );

}

TEST_F(sggqrfest, sggqrf1) {
    EXPECT_NEAR(0.0, sggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sggqrfest, sggqrf2) {
    EXPECT_NEAR(0.0, sggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sggqrfest, sggqrf3) {
    EXPECT_NEAR(0.0, sggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sggqrfest, sggqrf4) {
    EXPECT_NEAR(0.0, sggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class ggqrf_double_parameters{

   public:
	int bufsize, bufsize_b;
	void *hModule, *dModule;
	double diff, diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m, p;
	double* A, *b;
	lapack_int lda, ldb;	
	/*Output Parameter*/
	double* taua, *taub;
	double *Aref, *bref;
	double *tauaref, *taubref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      ggqrf_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggqrf_double_parameters ();

};

/* Constructor definition  double_common_parameters */
ggqrf_double_parameters:: ggqrf_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n ggqrf double:  m: %d, n: %d, p: %d, \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else if (matrix_layout == LAPACK_ROW_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&taua, &tauaref, min(n,m));
	lapacke_gtest_alloc_double_buffer_pair(&taub, &taubref, min(n,p));
		if ((A==NULL) || (Aref==NULL) ||\
		(b == NULL) || (bref == NULL) ||
		(taua == NULL) || (tauaref == NULL) ||
		(taub == NULL) || (taubref == NULL)){
		EXPECT_FALSE( true) << "ggqrf_double_parameters object: malloc error.";
		ggqrf_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
ggqrf_double_parameters :: ~ggqrf_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggqrf_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggqrf_free();

}
/*  Test fixture class definition */
class dggqrfest  : public  ::testing::Test {
public:
   ggqrf_double_parameters  *dggqrf_obj;
   void SetUp();
   void TearDown () { delete dggqrf_obj;}
};

void dggqrfest::SetUp(){

    /* LAPACKE dggqrf prototype */
    typedef int (*Fptr_NL_LAPACKE_dggqrf) (int matrix_layout, lapack_int n, lapack_int m, lapack_int p,\
											double* a, lapack_int lda, double* taua, double* b, lapack_int ldb, double* taub);

    Fptr_NL_LAPACKE_dggqrf dggqrf;

    dggqrf_obj = new ggqrf_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    dggqrf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dggqrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dggqrf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dggqrf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dggqrf = (Fptr_NL_LAPACKE_dggqrf)dlsym(dggqrf_obj->hModule, "LAPACKE_dggqrf");
    ASSERT_TRUE(dggqrf != NULL) << "failed to get the Netlib LAPACKE_dggqrf symbol";
    

    dggqrf_obj->inforef = dggqrf( dggqrf_obj->matrix_layout, dggqrf_obj->n, dggqrf_obj->m,
								dggqrf_obj->p, dggqrf_obj->Aref, dggqrf_obj->lda, dggqrf_obj->tauaref,
								dggqrf_obj->bref, dggqrf_obj->ldb, dggqrf_obj->taubref);

    /* Compute libflame's Lapacke o/p  */
    dggqrf_obj->info = LAPACKE_dggqrf( dggqrf_obj->matrix_layout, dggqrf_obj->n, dggqrf_obj->m,
								dggqrf_obj->p, dggqrf_obj->A, dggqrf_obj->lda, dggqrf_obj->taua,
								dggqrf_obj->b, dggqrf_obj->ldb, dggqrf_obj->taub);
	
	if( dggqrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dggqrf is wrong\n", dggqrf_obj->info );
    }
    if( dggqrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dggqrf is wrong\n", 
        dggqrf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dggqrf_obj->diff =  computeDiff_d( dggqrf_obj->bufsize, 
                dggqrf_obj->A, dggqrf_obj->Aref );
	dggqrf_obj->diff_b =  computeDiff_d( dggqrf_obj->bufsize_b, 
                dggqrf_obj->b, dggqrf_obj->bref );

}

TEST_F(dggqrfest, dggqrf1) {
    EXPECT_NEAR(0.0, dggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dggqrfest, dggqrf2) {
    EXPECT_NEAR(0.0, dggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dggqrfest, dggqrf3) {
    EXPECT_NEAR(0.0, dggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dggqrfest, dggqrf4) {
    EXPECT_NEAR(0.0, dggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class ggqrf_scomplex_parameters{

   public:
	int bufsize, bufsize_b;
	void *hModule, *dModule;
	float diff, diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m, p;
	lapack_complex_float* A, *b;
	lapack_int lda, ldb;	
	/*Output Parameter*/
	lapack_complex_float* taua, *taub;
	lapack_complex_float *Aref, *bref;
	lapack_complex_float *tauaref, *taubref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      ggqrf_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggqrf_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
ggqrf_scomplex_parameters:: ggqrf_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n ggqrf scomplex:  m: %d, n: %d, p: %d, \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else if (matrix_layout == LAPACK_ROW_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&taua, &tauaref, min(n,m));
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&taub, &taubref, min(n,p));
		if ((A==NULL) || (Aref==NULL) ||\
		(b == NULL) || (bref == NULL) ||
		(taua == NULL) || (tauaref == NULL) ||
		(taub == NULL) || (taubref == NULL)){
		EXPECT_FALSE( true) << "ggqrf_scomplex_parameters object: malloc error.";
		ggqrf_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
ggqrf_scomplex_parameters :: ~ggqrf_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggqrf_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggqrf_free();

}
/*  Test fixture class definition */
class cggqrfest  : public  ::testing::Test {
public:
   ggqrf_scomplex_parameters  *cggqrf_obj;
   void SetUp();
   void TearDown () { delete cggqrf_obj;}
};

void cggqrfest::SetUp(){

    /* LAPACKE cggqrf prototype */
    typedef int (*Fptr_NL_LAPACKE_cggqrf) (int matrix_layout, lapack_int n, lapack_int m, lapack_int p, lapack_complex_float* a, lapack_int lda,\
												lapack_complex_float* taua, lapack_complex_float* b, lapack_int ldb, lapack_complex_float* taub);

    Fptr_NL_LAPACKE_cggqrf cggqrf;

    cggqrf_obj = new ggqrf_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n); 
	
	idx = Circular_Increment_Index(idx);

    cggqrf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cggqrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cggqrf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cggqrf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cggqrf = (Fptr_NL_LAPACKE_cggqrf)dlsym(cggqrf_obj->hModule, "LAPACKE_cggqrf");
    ASSERT_TRUE(cggqrf != NULL) << "failed to get the Netlib LAPACKE_cggqrf symbol";
    

    cggqrf_obj->inforef = cggqrf( cggqrf_obj->matrix_layout, cggqrf_obj->n, cggqrf_obj->m,
								cggqrf_obj->p, cggqrf_obj->Aref, cggqrf_obj->lda, cggqrf_obj->tauaref,
								cggqrf_obj->bref, cggqrf_obj->ldb, cggqrf_obj->taubref);

    /* Compute libflame's Lapacke o/p  */
    cggqrf_obj->info = LAPACKE_cggqrf( cggqrf_obj->matrix_layout, cggqrf_obj->n, cggqrf_obj->m,
								cggqrf_obj->p, cggqrf_obj->A, cggqrf_obj->lda, cggqrf_obj->taua,
								cggqrf_obj->b, cggqrf_obj->ldb, cggqrf_obj->taub);
	
	if( cggqrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cggqrf is wrong\n", cggqrf_obj->info );
    }
    if( cggqrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cggqrf is wrong\n", 
        cggqrf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cggqrf_obj->diff =  computeDiff_c( cggqrf_obj->bufsize, 
                cggqrf_obj->A, cggqrf_obj->Aref );
	cggqrf_obj->diff_b =  computeDiff_c( cggqrf_obj->bufsize_b, 
                cggqrf_obj->b, cggqrf_obj->bref );

}

TEST_F(cggqrfest, cggqrf1) {
    EXPECT_NEAR(0.0, cggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cggqrfest, cggqrf2) {
    EXPECT_NEAR(0.0, cggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cggqrfest, cggqrf3) {
    EXPECT_NEAR(0.0, cggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cggqrfest, cggqrf4) {
    EXPECT_NEAR(0.0, cggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class ggqrf_dcomplex_parameters{

   public:
	int bufsize, bufsize_b;
	void *hModule, *dModule;
	double diff, diff_b;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m, p;
	lapack_complex_double* A, *b;
	lapack_int lda, ldb;	
	/*Output Parameter*/
	lapack_complex_double* taua, *taub;
	lapack_complex_double *Aref, *bref;
	lapack_complex_double *tauaref, *taubref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      ggqrf_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggqrf_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
ggqrf_dcomplex_parameters:: ggqrf_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n ggqrf dcomplex:  m: %d, n: %d, p: %d, \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else if (matrix_layout == LAPACK_ROW_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&taua, &tauaref, min(n,m));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&taub, &taubref, min(n,p));
		if ((A==NULL) || (Aref==NULL) ||\
		(b == NULL) || (bref == NULL) ||
		(taua == NULL) || (tauaref == NULL) ||
		(taub == NULL) || (taubref == NULL)){
		EXPECT_FALSE( true) << "ggqrf_dcomplex_parameters object: malloc error.";
		ggqrf_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
ggqrf_dcomplex_parameters :: ~ggqrf_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggqrf_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggqrf_free();

}
/*  Test fixture class definition */
class zggqrfest  : public  ::testing::Test {
public:
   ggqrf_dcomplex_parameters  *zggqrf_obj;
   void SetUp();
   void TearDown () { delete zggqrf_obj;}
};

void zggqrfest::SetUp(){

    /* LAPACKE zggqrf prototype */
    typedef int (*Fptr_NL_LAPACKE_zggqrf) (int matrix_layout, lapack_int n, lapack_int m, lapack_int p, lapack_complex_double* a,\
											lapack_int lda, lapack_complex_double* taua, lapack_complex_double* b, lapack_int ldb, lapack_complex_double* taub);

    Fptr_NL_LAPACKE_zggqrf zggqrf;

    zggqrf_obj = new ggqrf_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].n); 
	
	idx = Circular_Increment_Index(idx);

    zggqrf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zggqrf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zggqrf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zggqrf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zggqrf = (Fptr_NL_LAPACKE_zggqrf)dlsym(zggqrf_obj->hModule, "LAPACKE_zggqrf");
    ASSERT_TRUE(zggqrf != NULL) << "failed to get the Netlib LAPACKE_zggqrf symbol";
    

    zggqrf_obj->inforef = zggqrf( zggqrf_obj->matrix_layout, zggqrf_obj->n, zggqrf_obj->m,
								zggqrf_obj->p, zggqrf_obj->Aref, zggqrf_obj->lda, zggqrf_obj->tauaref,
								zggqrf_obj->bref, zggqrf_obj->ldb, zggqrf_obj->taubref);

    /* Compute libflame's Lapacke o/p  */
    zggqrf_obj->info = LAPACKE_zggqrf( zggqrf_obj->matrix_layout, zggqrf_obj->n, zggqrf_obj->m,
								zggqrf_obj->p, zggqrf_obj->A, zggqrf_obj->lda, zggqrf_obj->taua,
								zggqrf_obj->b, zggqrf_obj->ldb, zggqrf_obj->taub);
	
	if( zggqrf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zggqrf is wrong\n", zggqrf_obj->info );
    }
    if( zggqrf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zggqrf is wrong\n", 
        zggqrf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zggqrf_obj->diff =  computeDiff_z( zggqrf_obj->bufsize, 
                zggqrf_obj->A, zggqrf_obj->Aref );
	zggqrf_obj->diff_b =  computeDiff_z( zggqrf_obj->bufsize_b, 
                zggqrf_obj->b, zggqrf_obj->bref );

}

TEST_F(zggqrfest, zggqrf1) {
    EXPECT_NEAR(0.0, zggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zggqrfest, zggqrf2) {
    EXPECT_NEAR(0.0, zggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zggqrfest, zggqrf3) {
    EXPECT_NEAR(0.0, zggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zggqrfest, zggqrf4) {
    EXPECT_NEAR(0.0, zggqrf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zggqrf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}