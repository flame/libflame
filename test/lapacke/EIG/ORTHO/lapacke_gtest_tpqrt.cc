#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_helper.h"


#define tpqrt_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (t!=NULL)  free(t);\
if (tref!=NULL) free(tref); \
if (b!=NULL)    free(b); \
if (bref!=NULL)    free(bref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class tpqrt_float_parameters{

   public:
	int bufsize;
	int bufsize_b;
	int bufsize_t;
	void *hModule, *dModule;
	float diff_a;
	float diff_b;
	float diff_t;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt, lda, ldb;
	lapack_int l;
	float* A, *b;
	/*Output Parameter*/
	float* t;
	float *Aref, *tref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      tpqrt_float_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nb);
      ~tpqrt_float_parameters ();

};

/* Constructor definition  float_common_parameters */
tpqrt_float_parameters:: tpqrt_float_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nb = nb_i;
	l = min(m,n);
	lda =n;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpqrt float:  m: %d, n: %d  \n",  m, n);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR) {
		ldb = m;
		ldt = nb;
		bufsize_b = ldb*n;
		bufsize_t = ldt*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		ldb = n;
		ldt = n;
		bufsize_b = ldb*m;
		bufsize_t = ldt*nb;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize = lda*n;
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, bufsize_t);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "tpqrt_float_parameters object: malloc error.";
		tpqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand_custom_matrix( A, Aref, lda, n, 'U');
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tpqrt_float_parameters :: ~tpqrt_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpqrt_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpqrt_free();

}
/*  Test fixture class definition */
class stpqrt_test  : public  ::testing::Test {
public:
   tpqrt_float_parameters  *stpqrt_obj;
   void SetUp();
   void TearDown () { delete stpqrt_obj; }
};

void stpqrt_test::SetUp(){

    /* LAPACKE stpqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_stpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int l, lapack_int nb, float* a, lapack_int lda, float* b, lapack_int ldb, float* t, lapack_int ldt);
											

    Fptr_NL_LAPACKE_stpqrt stpqrt;

    stpqrt_obj = new tpqrt_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);

    idx = Circular_Increment_Index(idx);
	

    stpqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stpqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stpqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stpqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    stpqrt = (Fptr_NL_LAPACKE_stpqrt)dlsym(stpqrt_obj->hModule, "LAPACKE_stpqrt");
    ASSERT_TRUE(stpqrt != NULL) << "failed to get the Netlib LAPACKE_stpqrt symbol";
    

    stpqrt_obj->inforef = stpqrt( stpqrt_obj->matrix_layout, stpqrt_obj->m,
								stpqrt_obj->n, stpqrt_obj->l, stpqrt_obj->nb, stpqrt_obj->Aref,
								stpqrt_obj->lda, stpqrt_obj->bref, stpqrt_obj->ldb, stpqrt_obj->tref, stpqrt_obj->ldt);

    /* Compute libflame's Lapacke o/p  */
    stpqrt_obj->info = LAPACKE_stpqrt( stpqrt_obj->matrix_layout, stpqrt_obj->m, stpqrt_obj->n,
	stpqrt_obj->l, stpqrt_obj->nb, stpqrt_obj->A, stpqrt_obj->lda,stpqrt_obj->b,
	stpqrt_obj->ldb, stpqrt_obj->t, stpqrt_obj->ldt);

    if( stpqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stpqrt is wrong\n", stpqrt_obj->info );
    }
    if( stpqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stpqrt is wrong\n", 
        stpqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stpqrt_obj->diff_a =  computeDiff_s( stpqrt_obj->bufsize, 
                stpqrt_obj->A, stpqrt_obj->Aref );
	stpqrt_obj->diff_b =  computeDiff_s( stpqrt_obj->bufsize_b, stpqrt_obj->b, stpqrt_obj->bref );
	stpqrt_obj->diff_t =  computeDiff_s( stpqrt_obj->bufsize_t, stpqrt_obj->t, stpqrt_obj->tref );
								

}

TEST_F(stpqrt_test, stpqrt1) {
    EXPECT_NEAR(0.0, stpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stpqrt_test, stpqrt2) {
    EXPECT_NEAR(0.0, stpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stpqrt_test, stpqrt3) {
    EXPECT_NEAR(0.0, stpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stpqrt_test, stpqrt4) {
    EXPECT_NEAR(0.0, stpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class tpqrt_double_parameters{

   public:
	int bufsize;
	int bufsize_b;
	int bufsize_t;
	void *hModule, *dModule;
	double diff_a;
	double diff_b;
	double diff_t;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt, lda, ldb;
	lapack_int l;
	double* A, *b;
	/*Output Parameter*/
	double* t;
	double *Aref, *tref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      tpqrt_double_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nb);
      ~tpqrt_double_parameters ();

};

/* Constructor definition  double_common_parameters */
tpqrt_double_parameters:: tpqrt_double_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nb = nb_i;
	l = min(m,n);
	lda =n;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpqrt double:  m: %d, n: %d  \n",  m, n);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR) {
		ldb = m;
		ldt = nb;
		bufsize_b = ldb*n;
		bufsize_t = ldt*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		ldb = n;
		ldt = n;
		bufsize_b = ldb*m;
		bufsize_t = ldt*nb;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize = lda*n;
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&t, &tref, bufsize_t);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "tpqrt_double_parameters object: malloc error.";
		tpqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand_custom_matrix( A, Aref, lda, n, 'U');
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
tpqrt_double_parameters :: ~tpqrt_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpqrt_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpqrt_free();

}
/*  Test fixture class definition */
class dtpqrt_test  : public  ::testing::Test {
public:
   tpqrt_double_parameters  *dtpqrt_obj;
   void SetUp();
   void TearDown () { delete dtpqrt_obj; }
};

void dtpqrt_test::SetUp(){

    /* LAPACKE dtpqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_dtpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int l, lapack_int nb, double* a, lapack_int lda, double* b, lapack_int ldb, double* t, lapack_int ldt);
											

    Fptr_NL_LAPACKE_dtpqrt dtpqrt;

    dtpqrt_obj = new tpqrt_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);

    idx = Circular_Increment_Index(idx);
	

    dtpqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtpqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtpqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtpqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dtpqrt = (Fptr_NL_LAPACKE_dtpqrt)dlsym(dtpqrt_obj->hModule, "LAPACKE_dtpqrt");
    ASSERT_TRUE(dtpqrt != NULL) << "failed to get the Netlib LAPACKE_dtpqrt symbol";
    

    dtpqrt_obj->inforef = dtpqrt( dtpqrt_obj->matrix_layout, dtpqrt_obj->m,
								dtpqrt_obj->n, dtpqrt_obj->l, dtpqrt_obj->nb, dtpqrt_obj->Aref,
								dtpqrt_obj->lda, dtpqrt_obj->bref, dtpqrt_obj->ldb, dtpqrt_obj->tref, dtpqrt_obj->ldt);

    /* Compute libflame's Lapacke o/p  */
    dtpqrt_obj->info = LAPACKE_dtpqrt( dtpqrt_obj->matrix_layout, dtpqrt_obj->m, dtpqrt_obj->n,
	dtpqrt_obj->l, dtpqrt_obj->nb, dtpqrt_obj->A, dtpqrt_obj->lda,dtpqrt_obj->b,
	dtpqrt_obj->ldb, dtpqrt_obj->t, dtpqrt_obj->ldt);

    if( dtpqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtpqrt is wrong\n", dtpqrt_obj->info );
    }
    if( dtpqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtpqrt is wrong\n", 
        dtpqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtpqrt_obj->diff_a =  computeDiff_d( dtpqrt_obj->bufsize, 
                dtpqrt_obj->A, dtpqrt_obj->Aref );
	dtpqrt_obj->diff_b =  computeDiff_d( dtpqrt_obj->bufsize_b, dtpqrt_obj->b, dtpqrt_obj->bref );
	dtpqrt_obj->diff_t =  computeDiff_d( dtpqrt_obj->bufsize_t, dtpqrt_obj->t, dtpqrt_obj->tref );
								

}

TEST_F(dtpqrt_test, dtpqrt1) {
    EXPECT_NEAR(0.0, dtpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtpqrt_test, dtpqrt2) {
    EXPECT_NEAR(0.0, dtpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtpqrt_test, dtpqrt3) {
    EXPECT_NEAR(0.0, dtpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtpqrt_test, dtpqrt4) {
    EXPECT_NEAR(0.0, dtpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class tpqrt_scomplex_parameters{

   public:
	int bufsize;
	int bufsize_b;
	int bufsize_t;
	void *hModule, *dModule;
	float diff_a;
	float diff_b;
	float diff_t;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt, lda, ldb;
	lapack_int l;
	lapack_complex_float* A, *b;
	/*Output Parameter*/
	lapack_complex_float* t;
	lapack_complex_float *Aref, *tref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      tpqrt_scomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nb);
      ~tpqrt_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
tpqrt_scomplex_parameters:: tpqrt_scomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nb = nb_i;
	l = min(m,n);
	lda =n;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpqrt scomplex:  m: %d, n: %d  \n",  m, n);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR) {
		ldb = m;
		ldt = nb;
		bufsize_b = ldb*n;
		bufsize_t = ldt*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		ldb = n;
		ldt = n;
		bufsize_b = ldb*m;
		bufsize_t = ldt*nb;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize = lda*n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, bufsize_t);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "tpqrt_scomplex_parameters object: malloc error.";
		tpqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( A, Aref, lda, n, 'U');
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
tpqrt_scomplex_parameters :: ~tpqrt_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpqrt_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpqrt_free();

}
/*  Test fixture class definition */
class ctpqrt_test  : public  ::testing::Test {
public:
   tpqrt_scomplex_parameters  *ctpqrt_obj;
   void SetUp();
   void TearDown () { delete ctpqrt_obj; }
};

void ctpqrt_test::SetUp(){

    /* LAPACKE ctpqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_ctpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int l, lapack_int nb, lapack_complex_float* a, lapack_int lda, lapack_complex_float* b, lapack_int ldb, lapack_complex_float* t, lapack_int ldt);
											

    Fptr_NL_LAPACKE_ctpqrt ctpqrt;

    ctpqrt_obj = new tpqrt_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);

    idx = Circular_Increment_Index(idx);
	

    ctpqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctpqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctpqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctpqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ctpqrt = (Fptr_NL_LAPACKE_ctpqrt)dlsym(ctpqrt_obj->hModule, "LAPACKE_ctpqrt");
    ASSERT_TRUE(ctpqrt != NULL) << "failed to get the Netlib LAPACKE_ctpqrt symbol";
    

    ctpqrt_obj->inforef = ctpqrt( ctpqrt_obj->matrix_layout, ctpqrt_obj->m,
								ctpqrt_obj->n, ctpqrt_obj->l, ctpqrt_obj->nb, ctpqrt_obj->Aref,
								ctpqrt_obj->lda, ctpqrt_obj->bref, ctpqrt_obj->ldb, ctpqrt_obj->tref, ctpqrt_obj->ldt);

    /* Compute libflame's Lapacke o/p  */
    ctpqrt_obj->info = LAPACKE_ctpqrt( ctpqrt_obj->matrix_layout, ctpqrt_obj->m, ctpqrt_obj->n,
	ctpqrt_obj->l, ctpqrt_obj->nb, ctpqrt_obj->A, ctpqrt_obj->lda,ctpqrt_obj->b,
	ctpqrt_obj->ldb, ctpqrt_obj->t, ctpqrt_obj->ldt);

    if( ctpqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctpqrt is wrong\n", ctpqrt_obj->info );
    }
    if( ctpqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctpqrt is wrong\n", 
        ctpqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctpqrt_obj->diff_a =  computeDiff_c( ctpqrt_obj->bufsize, 
                ctpqrt_obj->A, ctpqrt_obj->Aref );
	ctpqrt_obj->diff_b =  computeDiff_c( ctpqrt_obj->bufsize_b, ctpqrt_obj->b, ctpqrt_obj->bref );
	ctpqrt_obj->diff_t =  computeDiff_c( ctpqrt_obj->bufsize_t, ctpqrt_obj->t, ctpqrt_obj->tref );
								

}

TEST_F(ctpqrt_test, ctpqrt1) {
    EXPECT_NEAR(0.0, ctpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctpqrt_test, ctpqrt2) {
    EXPECT_NEAR(0.0, ctpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctpqrt_test, ctpqrt3) {
    EXPECT_NEAR(0.0, ctpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctpqrt_test, ctpqrt4) {
    EXPECT_NEAR(0.0, ctpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}
/* Begin dcomplex_common_parameters  class definition */
class tpqrt_dcomplex_parameters{

   public:
	int bufsize;
	int bufsize_b;
	int bufsize_t;
	void *hModule, *dModule;
	double diff_a;
	double diff_b;
	double diff_t;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int nb;
	lapack_int ldt, lda, ldb;
	lapack_int l;
	lapack_complex_double* A, *b;
	/*Output Parameter*/
	lapack_complex_double* t;
	lapack_complex_double *Aref, *tref, *bref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      tpqrt_dcomplex_parameters (int matrix_layout, lapack_int m , lapack_int n, lapack_int nb);
      ~tpqrt_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
tpqrt_dcomplex_parameters:: tpqrt_dcomplex_parameters (int matrix_layout_i, lapack_int m_i , lapack_int n_i, lapack_int nb_i)
{
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	nb = nb_i;
	l = min(m,n);
	lda =n;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpqrt dcomplex:  m: %d, n: %d  \n",  m, n);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR) {
		ldb = m;
		ldt = nb;
		bufsize_b = ldb*n;
		bufsize_t = ldt*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR) {
		ldb = n;
		ldt = n;
		bufsize_b = ldb*m;
		bufsize_t = ldt*nb;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize = lda*n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, bufsize_t);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) ||\
		(b==NULL) || (bref==NULL) ||\
		(t==NULL) || (tref==NULL)){
		EXPECT_FALSE( true) << "tpqrt_dcomplex_parameters object: malloc error.";
		tpqrt_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( A, Aref, lda, n, 'U');
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
tpqrt_dcomplex_parameters :: ~tpqrt_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpqrt_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpqrt_free();

}
/*  Test fixture class definition */
class ztpqrt_test  : public  ::testing::Test {
public:
   tpqrt_dcomplex_parameters  *ztpqrt_obj;
   void SetUp();
   void TearDown () { delete ztpqrt_obj; }
};

void ztpqrt_test::SetUp(){

    /* LAPACKE ztpqrt prototype */
    typedef int (*Fptr_NL_LAPACKE_ztpqrt) (int matrix_layout, lapack_int m, lapack_int n,\
	lapack_int l, lapack_int nb, lapack_complex_double* a, lapack_int lda, lapack_complex_double* b, lapack_int ldb, lapack_complex_double* t, lapack_int ldt);
											

    Fptr_NL_LAPACKE_ztpqrt ztpqrt;

    ztpqrt_obj = new tpqrt_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].m,
                           eig_paramslist[idx].n,
						   eig_paramslist[idx].nb);

    idx = Circular_Increment_Index(idx);
	

    ztpqrt_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztpqrt_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztpqrt_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztpqrt_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ztpqrt = (Fptr_NL_LAPACKE_ztpqrt)dlsym(ztpqrt_obj->hModule, "LAPACKE_ztpqrt");
    ASSERT_TRUE(ztpqrt != NULL) << "failed to get the Netlib LAPACKE_ztpqrt symbol";
    

    ztpqrt_obj->inforef = ztpqrt( ztpqrt_obj->matrix_layout, ztpqrt_obj->m,
								ztpqrt_obj->n, ztpqrt_obj->l, ztpqrt_obj->nb, ztpqrt_obj->Aref,
								ztpqrt_obj->lda, ztpqrt_obj->bref, ztpqrt_obj->ldb, ztpqrt_obj->tref, ztpqrt_obj->ldt);

    /* Compute libflame's Lapacke o/p  */
    ztpqrt_obj->info = LAPACKE_ztpqrt( ztpqrt_obj->matrix_layout, ztpqrt_obj->m, ztpqrt_obj->n,
	ztpqrt_obj->l, ztpqrt_obj->nb, ztpqrt_obj->A, ztpqrt_obj->lda,ztpqrt_obj->b,
	ztpqrt_obj->ldb, ztpqrt_obj->t, ztpqrt_obj->ldt);

    if( ztpqrt_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztpqrt is wrong\n", ztpqrt_obj->info );
    }
    if( ztpqrt_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztpqrt is wrong\n", 
        ztpqrt_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztpqrt_obj->diff_a =  computeDiff_z( ztpqrt_obj->bufsize, 
                ztpqrt_obj->A, ztpqrt_obj->Aref );
	ztpqrt_obj->diff_b =  computeDiff_z( ztpqrt_obj->bufsize_b, ztpqrt_obj->b, ztpqrt_obj->bref );
	ztpqrt_obj->diff_t =  computeDiff_z( ztpqrt_obj->bufsize_t, ztpqrt_obj->t, ztpqrt_obj->tref );
								

}

TEST_F(ztpqrt_test, ztpqrt1) {
    EXPECT_NEAR(0.0, ztpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztpqrt_test, ztpqrt2) {
    EXPECT_NEAR(0.0, ztpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztpqrt_test, ztpqrt3) {
    EXPECT_NEAR(0.0, ztpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztpqrt_test, ztpqrt4) {
    EXPECT_NEAR(0.0, ztpqrt_obj->diff_a, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt_obj->diff_b, LAPACKE_EIG_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt_obj->diff_t, LAPACKE_EIG_THRESHOLD);
}
