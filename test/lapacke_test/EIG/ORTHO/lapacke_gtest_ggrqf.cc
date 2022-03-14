#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define ggrqf_free() \
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
class ggrqf_float_parameters{

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
      ggrqf_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggrqf_float_parameters ();

};

/* Constructor definition  float_common_parameters */
ggrqf_float_parameters:: ggrqf_float_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n ggrqf float:  m: %d, n: %d, p: %d, \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&taua, &tauaref, min(m,n));
	lapacke_gtest_alloc_float_buffer_pair(&taub, &taubref, min(p,n));
		if ((A==NULL) || (Aref==NULL) ||\
		(b == NULL) || (bref == NULL) ||
		(taua == NULL) || (tauaref == NULL) ||
		(taub == NULL) || (taubref == NULL)){
		EXPECT_FALSE( true) << "ggrqf_float_parameters object: malloc error.";
		ggrqf_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
ggrqf_float_parameters :: ~ggrqf_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggrqf_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggrqf_free();

}
/*  Test fixture class definition */
class sggrqfest  : public  ::testing::Test {
public:
   ggrqf_float_parameters  *sggrqf_obj;
   void SetUp();
   void TearDown () { delete sggrqf_obj;}
};

void sggrqfest::SetUp(){

    /* LAPACKE sggrqf prototype */
    typedef int (*Fptr_NL_LAPACKE_sggrqf) (int matrix_layout, lapack_int m, lapack_int p, lapack_int n, float* a,\
											lapack_int lda, float* taua, float* b, lapack_int ldb, float* taub);

    Fptr_NL_LAPACKE_sggrqf sggrqf;

    sggrqf_obj = new ggrqf_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    sggrqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    sggrqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(sggrqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(sggrqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    sggrqf = (Fptr_NL_LAPACKE_sggrqf)dlsym(sggrqf_obj->hModule, "LAPACKE_sggrqf");
    ASSERT_TRUE(sggrqf != NULL) << "failed to get the Netlib LAPACKE_sggrqf symbol";
    

    sggrqf_obj->inforef = sggrqf( sggrqf_obj->matrix_layout, sggrqf_obj->m, sggrqf_obj->p,
								sggrqf_obj->n, sggrqf_obj->Aref, sggrqf_obj->lda, sggrqf_obj->tauaref,
								sggrqf_obj->bref, sggrqf_obj->ldb, sggrqf_obj->taubref);

    /* Compute libflame's Lapacke o/p  */
    sggrqf_obj->info = LAPACKE_sggrqf( sggrqf_obj->matrix_layout, sggrqf_obj->m, sggrqf_obj->p,
								sggrqf_obj->n, sggrqf_obj->A, sggrqf_obj->lda, sggrqf_obj->taua,
								sggrqf_obj->b, sggrqf_obj->ldb, sggrqf_obj->taub);
	
	if( sggrqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_sggrqf is wrong\n", sggrqf_obj->info );
    }
    if( sggrqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_sggrqf is wrong\n", 
        sggrqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    sggrqf_obj->diff =  computeDiff_s( sggrqf_obj->bufsize, 
                sggrqf_obj->A, sggrqf_obj->Aref );
	sggrqf_obj->diff_b =  computeDiff_s( sggrqf_obj->bufsize_b, 
                sggrqf_obj->b, sggrqf_obj->bref );

}

TEST_F(sggrqfest, sggrqf1) {
    EXPECT_NEAR(0.0, sggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sggrqfest, sggrqf2) {
    EXPECT_NEAR(0.0, sggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sggrqfest, sggrqf3) {
    EXPECT_NEAR(0.0, sggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(sggrqfest, sggrqf4) {
    EXPECT_NEAR(0.0, sggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, sggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class ggrqf_double_parameters{

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
      ggrqf_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggrqf_double_parameters ();

};

/* Constructor definition  double_common_parameters */
ggrqf_double_parameters:: ggrqf_double_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n ggrqf double:  m: %d, n: %d, p: %d, \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&taua, &tauaref, min(m,n));
	lapacke_gtest_alloc_double_buffer_pair(&taub, &taubref, min(p,n));
		if ((A==NULL) || (Aref==NULL) ||\
		(b == NULL) || (bref == NULL) ||
		(taua == NULL) || (tauaref == NULL) ||
		(taub == NULL) || (taubref == NULL)){
		EXPECT_FALSE( true) << "ggrqf_double_parameters object: malloc error.";
		ggrqf_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
ggrqf_double_parameters :: ~ggrqf_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggrqf_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggrqf_free();

}
/*  Test fixture class definition */
class dggrqfest  : public  ::testing::Test {
public:
   ggrqf_double_parameters  *dggrqf_obj;
   void SetUp();
   void TearDown () { delete dggrqf_obj;}
};

void dggrqfest::SetUp(){

    /* LAPACKE dggrqf prototype */
    typedef int (*Fptr_NL_LAPACKE_dggrqf) (int matrix_layout, lapack_int m, lapack_int p, lapack_int n, double* a,\
											lapack_int lda, double* taua, double* b, lapack_int ldb, double* taub);

    Fptr_NL_LAPACKE_dggrqf dggrqf;

    dggrqf_obj = new ggrqf_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    dggrqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dggrqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dggrqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dggrqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dggrqf = (Fptr_NL_LAPACKE_dggrqf)dlsym(dggrqf_obj->hModule, "LAPACKE_dggrqf");
    ASSERT_TRUE(dggrqf != NULL) << "failed to get the Netlib LAPACKE_dggrqf symbol";
    

    dggrqf_obj->inforef = dggrqf( dggrqf_obj->matrix_layout, dggrqf_obj->m, dggrqf_obj->p,
								dggrqf_obj->n, dggrqf_obj->Aref, dggrqf_obj->lda, dggrqf_obj->tauaref,
								dggrqf_obj->bref, dggrqf_obj->ldb, dggrqf_obj->taubref);

    /* Compute libflame's Lapacke o/p  */
    dggrqf_obj->info = LAPACKE_dggrqf( dggrqf_obj->matrix_layout, dggrqf_obj->m, dggrqf_obj->p,
								dggrqf_obj->n, dggrqf_obj->A, dggrqf_obj->lda, dggrqf_obj->taua,
								dggrqf_obj->b, dggrqf_obj->ldb, dggrqf_obj->taub);
	
	if( dggrqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dggrqf is wrong\n", dggrqf_obj->info );
    }
    if( dggrqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dggrqf is wrong\n", 
        dggrqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dggrqf_obj->diff =  computeDiff_d( dggrqf_obj->bufsize, 
                dggrqf_obj->A, dggrqf_obj->Aref );
	dggrqf_obj->diff_b =  computeDiff_d( dggrqf_obj->bufsize_b, 
                dggrqf_obj->b, dggrqf_obj->bref );

}

TEST_F(dggrqfest, dggrqf1) {
    EXPECT_NEAR(0.0, dggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dggrqfest, dggrqf2) {
    EXPECT_NEAR(0.0, dggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dggrqfest, dggrqf3) {
    EXPECT_NEAR(0.0, dggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dggrqfest, dggrqf4) {
    EXPECT_NEAR(0.0, dggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class ggrqf_scomplex_parameters{

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
      ggrqf_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggrqf_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
ggrqf_scomplex_parameters:: ggrqf_scomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n ggrqf scomplex:  m: %d, n: %d, p: %d, \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&taua, &tauaref, min(m,n));
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&taub, &taubref, min(p,n));
		if ((A==NULL) || (Aref==NULL) ||\
		(b == NULL) || (bref == NULL) ||
		(taua == NULL) || (tauaref == NULL) ||
		(taub == NULL) || (taubref == NULL)){
		EXPECT_FALSE( true) << "ggrqf_scomplex_parameters object: malloc error.";
		ggrqf_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
ggrqf_scomplex_parameters :: ~ggrqf_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggrqf_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggrqf_free();

}
/*  Test fixture class definition */
class cggrqfest  : public  ::testing::Test {
public:
   ggrqf_scomplex_parameters  *cggrqf_obj;
   void SetUp();
   void TearDown () { delete cggrqf_obj;}
};

void cggrqfest::SetUp(){

    /* LAPACKE cggrqf prototype */
    typedef int (*Fptr_NL_LAPACKE_cggrqf) (int matrix_layout, lapack_int m, lapack_int p, lapack_int n, lapack_complex_float* a,\
											lapack_int lda, lapack_complex_float* taua, lapack_complex_float* b, lapack_int ldb, lapack_complex_float* taub);

    Fptr_NL_LAPACKE_cggrqf cggrqf;

    cggrqf_obj = new ggrqf_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    cggrqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    cggrqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(cggrqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(cggrqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    cggrqf = (Fptr_NL_LAPACKE_cggrqf)dlsym(cggrqf_obj->hModule, "LAPACKE_cggrqf");
    ASSERT_TRUE(cggrqf != NULL) << "failed to get the Netlib LAPACKE_cggrqf symbol";
    

    cggrqf_obj->inforef = cggrqf( cggrqf_obj->matrix_layout, cggrqf_obj->m, cggrqf_obj->p,
								cggrqf_obj->n, cggrqf_obj->Aref, cggrqf_obj->lda, cggrqf_obj->tauaref,
								cggrqf_obj->bref, cggrqf_obj->ldb, cggrqf_obj->taubref);

    /* Compute libflame's Lapacke o/p  */
    cggrqf_obj->info = LAPACKE_cggrqf( cggrqf_obj->matrix_layout, cggrqf_obj->m, cggrqf_obj->p,
								cggrqf_obj->n, cggrqf_obj->A, cggrqf_obj->lda, cggrqf_obj->taua,
								cggrqf_obj->b, cggrqf_obj->ldb, cggrqf_obj->taub);
	
	if( cggrqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_cggrqf is wrong\n", cggrqf_obj->info );
    }
    if( cggrqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_cggrqf is wrong\n", 
        cggrqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    cggrqf_obj->diff =  computeDiff_c( cggrqf_obj->bufsize, 
                cggrqf_obj->A, cggrqf_obj->Aref );
	cggrqf_obj->diff_b =  computeDiff_c( cggrqf_obj->bufsize_b, 
                cggrqf_obj->b, cggrqf_obj->bref );

}

TEST_F(cggrqfest, cggrqf1) {
    EXPECT_NEAR(0.0, cggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cggrqfest, cggrqf2) {
    EXPECT_NEAR(0.0, cggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cggrqfest, cggrqf3) {
    EXPECT_NEAR(0.0, cggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(cggrqfest, cggrqf4) {
    EXPECT_NEAR(0.0, cggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, cggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class ggrqf_dcomplex_parameters{

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
      ggrqf_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int p);
      ~ggrqf_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
ggrqf_dcomplex_parameters:: ggrqf_dcomplex_parameters (int matrix_layout_i, lapack_int m_i, lapack_int n_i, lapack_int p_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	p = p_i;
	
	#if LAPACKE_TEST_VERBOSE
		printf(" \n ggrqf dcomplex:  m: %d, n: %d, p: %d, \n",  m, n, p);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR){
		lda = m;
		ldb = p;
		bufsize = lda*n;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR){
		lda = n;
		ldb = n;
		bufsize = lda*m;
		bufsize_b = ldb*p;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&taua, &tauaref, min(m,n));
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&taub, &taubref, min(p,n));
		if ((A==NULL) || (Aref==NULL) ||\
		(b == NULL) || (bref == NULL) ||
		(taua == NULL) || (tauaref == NULL) ||
		(taub == NULL) || (taubref == NULL)){
		EXPECT_FALSE( true) << "ggrqf_dcomplex_parameters object: malloc error.";
		ggrqf_free();
		exit(0);
	}	
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
ggrqf_dcomplex_parameters :: ~ggrqf_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" ggrqf_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   ggrqf_free();

}
/*  Test fixture class definition */
class zggrqfest  : public  ::testing::Test {
public:
   ggrqf_dcomplex_parameters  *zggrqf_obj;
   void SetUp();
   void TearDown () { delete zggrqf_obj;}
};

void zggrqfest::SetUp(){

    /* LAPACKE zggrqf prototype */
    typedef int (*Fptr_NL_LAPACKE_zggrqf) (int matrix_layout, lapack_int m, lapack_int p, lapack_int n, lapack_complex_double* a,\
											lapack_int lda, lapack_complex_double* taua, lapack_complex_double* b, lapack_int ldb, lapack_complex_double* taub);

    Fptr_NL_LAPACKE_zggrqf zggrqf;

    zggrqf_obj = new ggrqf_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
                           eig_paramslist[idx].m); 
	
	idx = Circular_Increment_Index(idx);

    zggrqf_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zggrqf_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zggrqf_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zggrqf_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zggrqf = (Fptr_NL_LAPACKE_zggrqf)dlsym(zggrqf_obj->hModule, "LAPACKE_zggrqf");
    ASSERT_TRUE(zggrqf != NULL) << "failed to get the Netlib LAPACKE_zggrqf symbol";
    

    zggrqf_obj->inforef = zggrqf( zggrqf_obj->matrix_layout, zggrqf_obj->m, zggrqf_obj->p,
								zggrqf_obj->n, zggrqf_obj->Aref, zggrqf_obj->lda, zggrqf_obj->tauaref,
								zggrqf_obj->bref, zggrqf_obj->ldb, zggrqf_obj->taubref);

    /* Compute libflame's Lapacke o/p  */
    zggrqf_obj->info = LAPACKE_zggrqf( zggrqf_obj->matrix_layout, zggrqf_obj->m, zggrqf_obj->p,
								zggrqf_obj->n, zggrqf_obj->A, zggrqf_obj->lda, zggrqf_obj->taua,
								zggrqf_obj->b, zggrqf_obj->ldb, zggrqf_obj->taub);
	
	if( zggrqf_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zggrqf is wrong\n", zggrqf_obj->info );
    }
    if( zggrqf_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zggrqf is wrong\n", 
        zggrqf_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zggrqf_obj->diff =  computeDiff_z( zggrqf_obj->bufsize, 
                zggrqf_obj->A, zggrqf_obj->Aref );
	zggrqf_obj->diff_b =  computeDiff_z( zggrqf_obj->bufsize_b, 
                zggrqf_obj->b, zggrqf_obj->bref );

}

TEST_F(zggrqfest, zggrqf1) {
    EXPECT_NEAR(0.0, zggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zggrqfest, zggrqf2) {
    EXPECT_NEAR(0.0, zggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zggrqfest, zggrqf3) {
    EXPECT_NEAR(0.0, zggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}

TEST_F(zggrqfest, zggrqf4) {
    EXPECT_NEAR(0.0, zggrqf_obj->diff, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, zggrqf_obj->diff_b, LAPACKE_EIG_THRESHOLD);
}