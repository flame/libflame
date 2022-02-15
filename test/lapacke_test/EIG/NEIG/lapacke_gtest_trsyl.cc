#include "gtest/gtest.h"
#include "../../lapacke_gtest_main.h"
#include "../../lapacke_gtest_aux.h"


#define trsyl_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)    free(b); \
if (bref!=NULL) free(bref);\
if (c!=NULL)    free(c); \
if (cref!=NULL) free(cref);\
if (scale!=NULL)  free(scale);\
if (scaleref!=NULL) free(scaleref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule);

// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class trsyl_float_parameters{

   public:
	int bufsize;
	int bufsize_a, bufsize_b;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char transa;
	char transb;
	float* A;
	float* b;
	float* c;	
	lapack_int lda;
	lapack_int ldb;
	lapack_int ldc;
	lapack_int isgn;	
	/*Output Parameter*/
	float* scale;
	float *Aref, *scaleref;
	float *bref, *cref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      trsyl_float_parameters (int matrix_layout, char transa, char transb, lapack_int isgn,\
							lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb, lapack_int ldc);
      ~trsyl_float_parameters ();

};

/* Constructor definition  float_common_parameters */
trsyl_float_parameters:: trsyl_float_parameters (int matrix_layout_i, char transa_i, char transb_i, lapack_int isgn_i,\
													lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i, lapack_int ldc_i)
{
	
	matrix_layout = matrix_layout_i;
	transa = transa_i;
	transb = transb_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldb = ldb_i;
	ldc = ldc_i;
	isgn = isgn_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trsyl float:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldc*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldc*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_a = lda*m;
	bufsize_b = ldb*n;
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_float_buffer_pair(&c, &cref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&scale, &scaleref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL) || \
		(c == NULL) || (cref == NULL) ||
		(scale == NULL) || (scaleref == NULL)){
		EXPECT_FALSE( true) << "trsyl_float_parameters object: malloc error.";
		trsyl_free();
		exit(0);
	}
	/* Initialization of input matrices as quasi-triangular matrix */
	lapacke_gtest_init_float_buffer_pair_with_constant(A, Aref, bufsize_a, 0);
	lapacke_gtest_init_float_buffer_pair_with_constant(b, bref, bufsize_b, 0);
	lapacke_gtest_init_float_buffer_pair_rand( c, cref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
trsyl_float_parameters :: ~trsyl_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trsyl_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trsyl_free();

}
/*  Test fixture class definition */
class strsyl_test  : public  ::testing::Test {
public:
   trsyl_float_parameters  *strsyl_obj;
   void SetUp();
   void TearDown () { delete strsyl_obj; }
};

void strsyl_test::SetUp(){

    /* LAPACKE strsyl prototype */
    typedef int (*Fptr_NL_LAPACKE_strsyl) (int matrix_layout, char transa, char transb, lapack_int isgn,\
											lapack_int m, lapack_int n, const float *A, lapack_int lda, \ 
											const float* b, lapack_int ldb, const float* c, lapack_int ldc, float* scale);

    Fptr_NL_LAPACKE_strsyl strsyl;

    strsyl_obj = new trsyl_float_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].isgn,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
						   eig_paramslist[idx].lda,
						   eig_paramslist[idx].ldb,
						   eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    strsyl_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    strsyl_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(strsyl_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(strsyl_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    strsyl = (Fptr_NL_LAPACKE_strsyl)dlsym(strsyl_obj->hModule, "LAPACKE_strsyl");
    ASSERT_TRUE(strsyl != NULL) << "failed to get the Netlib LAPACKE_strsyl symbol";
    

    strsyl_obj->inforef = strsyl( strsyl_obj->matrix_layout,strsyl_obj->transa, strsyl_obj->transb,\
								strsyl_obj->isgn, strsyl_obj->m, strsyl_obj->n, strsyl_obj->Aref, \
								strsyl_obj->lda, strsyl_obj->bref, strsyl_obj->ldb, strsyl_obj->cref,\
								strsyl_obj->ldc, strsyl_obj->scaleref);

    /* Compute libflame's Lapacke o/p  */
    strsyl_obj->info = LAPACKE_strsyl( strsyl_obj->matrix_layout,strsyl_obj->transa, strsyl_obj->transb,\
								strsyl_obj->isgn, strsyl_obj->m, strsyl_obj->n, strsyl_obj->A, \
								strsyl_obj->lda, strsyl_obj->b, strsyl_obj->ldb, strsyl_obj->c,\
								strsyl_obj->ldc, strsyl_obj->scale);

    if( strsyl_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_strsyl is wrong\n", strsyl_obj->info );
    }
    if( strsyl_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_strsyl is wrong\n", 
        strsyl_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    strsyl_obj->diff =  computeDiff_s( strsyl_obj->bufsize, 
                strsyl_obj->c, strsyl_obj->cref );

}

TEST_F(strsyl_test, strsyl1) {
    EXPECT_NEAR(0.0, strsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(strsyl_test, strsyl2) {
    EXPECT_NEAR(0.0, strsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(strsyl_test, strsyl3) {
    EXPECT_NEAR(0.0, strsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(strsyl_test, strsyl4) {
    EXPECT_NEAR(0.0, strsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class trsyl_double_parameters{

   public:
	int bufsize, bufsize_a, bufsize_b;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char transa;
	char transb;
	double* A;
	double* b;
	double* c;	
	lapack_int lda;
	lapack_int ldb;
	lapack_int ldc;
	lapack_int isgn;	
	/*Output Parameter*/
	double* scale;
	double *Aref, *scaleref;
	double *bref, *cref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      trsyl_double_parameters (int matrix_layout, char transa, char transb, lapack_int isgn,\
							lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb, lapack_int ldc);
      ~trsyl_double_parameters ();

};

/* Constructor definition  double_common_parameters */
trsyl_double_parameters:: trsyl_double_parameters (int matrix_layout_i, char transa_i, char transb_i, lapack_int isgn_i,\
													lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i, lapack_int ldc_i)
{
	
	matrix_layout = matrix_layout_i;
	transa = transa_i;
	transb = transb_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldb = ldb_i;
	ldc = ldc_i;
	isgn = isgn_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trsyl double:  m: %d, n: %d lda: %d \n",  m, n, lda);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldc*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldc*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_a = lda*m;
	bufsize_b = ldb*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_double_buffer_pair(&c, &cref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&scale, &scaleref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL) || \
		(c == NULL) || (cref == NULL) ||
		(scale == NULL) ||(scaleref == NULL)){
		EXPECT_FALSE( true) << "trsyl_double_parameters object: malloc error.";
		trsyl_free();
		exit(0);
	}
	/* Initialization of input matrices as quasi-triangular matrix */
	lapacke_gtest_init_double_buffer_pair_with_constant(A, Aref, bufsize_a, 0);
	lapacke_gtest_init_double_buffer_pair_with_constant(b, bref, bufsize_b, 0);
	lapacke_gtest_init_double_buffer_pair_rand( c, cref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
trsyl_double_parameters :: ~trsyl_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trsyl_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trsyl_free();

}
/*  Test fixture class definition */
class dtrsyl_test  : public  ::testing::Test {
public:
   trsyl_double_parameters  *dtrsyl_obj;
   void SetUp();
   void TearDown () { delete dtrsyl_obj; }
};

void dtrsyl_test::SetUp(){

    /* LAPACKE dtrsyl prototype */
    typedef int (*Fptr_NL_LAPACKE_dtrsyl) (int matrix_layout, char transa, char transb, lapack_int isgn,\
											lapack_int m, lapack_int n, const double *A, lapack_int lda, \ 
											const double* b, lapack_int ldb, const double* c, lapack_int ldc, double* scale);

    Fptr_NL_LAPACKE_dtrsyl dtrsyl;

    dtrsyl_obj = new trsyl_double_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].isgn,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
						   eig_paramslist[idx].lda,
						   eig_paramslist[idx].ldb,
						   eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    dtrsyl_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtrsyl_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtrsyl_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtrsyl_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dtrsyl = (Fptr_NL_LAPACKE_dtrsyl)dlsym(dtrsyl_obj->hModule, "LAPACKE_dtrsyl");
    ASSERT_TRUE(dtrsyl != NULL) << "failed to get the Netlib LAPACKE_dtrsyl symbol";
    

    dtrsyl_obj->inforef = dtrsyl( dtrsyl_obj->matrix_layout,dtrsyl_obj->transa, dtrsyl_obj->transb,\
								dtrsyl_obj->isgn, dtrsyl_obj->m, dtrsyl_obj->n, dtrsyl_obj->Aref, \
								dtrsyl_obj->lda, dtrsyl_obj->bref, dtrsyl_obj->ldb, dtrsyl_obj->cref,\
								dtrsyl_obj->ldc, dtrsyl_obj->scaleref);

    /* Compute libflame's Lapacke o/p  */
    dtrsyl_obj->info = LAPACKE_dtrsyl( dtrsyl_obj->matrix_layout,dtrsyl_obj->transa, dtrsyl_obj->transb,\
								dtrsyl_obj->isgn, dtrsyl_obj->m, dtrsyl_obj->n, dtrsyl_obj->A, \
								dtrsyl_obj->lda, dtrsyl_obj->b, dtrsyl_obj->ldb, dtrsyl_obj->c,\
								dtrsyl_obj->ldc, dtrsyl_obj->scale);

    if( dtrsyl_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtrsyl is wrong\n", dtrsyl_obj->info );
    }
    if( dtrsyl_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtrsyl is wrong\n", 
        dtrsyl_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtrsyl_obj->diff =  computeDiff_d( dtrsyl_obj->bufsize, 
                dtrsyl_obj->c, dtrsyl_obj->cref );

}

TEST_F(dtrsyl_test, dtrsyl1) {
    EXPECT_NEAR(0.0, dtrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtrsyl_test, dtrsyl2) {
    EXPECT_NEAR(0.0, dtrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtrsyl_test, dtrsyl3) {
    EXPECT_NEAR(0.0, dtrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtrsyl_test, dtrsyl4) {
    EXPECT_NEAR(0.0, dtrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class trsyl_scomplex_parameters{

   public:
	int bufsize, bufsize_a, bufsize_b;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char transa;
	char transb;
	lapack_complex_float* A;
	lapack_complex_float* b;
	lapack_complex_float* c;	
	lapack_int lda;
	lapack_int ldb;
	lapack_int ldc;
	lapack_int isgn;	
	/*Output Parameter*/
	float* scale, *scaleref;
	lapack_complex_float *Aref;
	lapack_complex_float *bref, *cref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      trsyl_scomplex_parameters (int matrix_layout, char transa, char transb, lapack_int isgn,\
							lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb, lapack_int ldc);
      ~trsyl_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
trsyl_scomplex_parameters:: trsyl_scomplex_parameters (int matrix_layout_i, char transa_i, char transb_i, lapack_int isgn_i,\
													lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i, lapack_int ldc_i)
{
	
	matrix_layout = matrix_layout_i;
	transa = transa_i;
	transb = transb_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldb = ldb_i;
	ldc = ldc_i;
	isgn = isgn_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trsyl lapack_complex_float:  m: %d, n: %d lda: %d, ldb:%d, ldc:%d, transa:%c \n",  m, n, lda, ldb,ldc, transa);
	#endif
	
	if ((transa == 'T') ||(transb == 'T')) // trans = 'T' supports  real flavours only
	{
		transa = transb = 'C';
	}
	
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldc*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldc*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_a = lda*m;
	bufsize_b = ldb*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&c, &cref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&scale, &scaleref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL) || \
		(c == NULL) || (cref == NULL) ||
		(scale == NULL) || (scaleref == NULL)){
		EXPECT_FALSE( true) << "trsyl_lapack_complex_float_parameters object: malloc error.";
		trsyl_free();
		exit(0);
	}
	/* Initialization of input matrices as quasi-triangular matrix */
	lapacke_gtest_init_scomplex_buffer_pair_with_constant(A, Aref, bufsize_a, 0);
	lapacke_gtest_init_scomplex_buffer_pair_with_constant(b, bref, bufsize_b, 0);
	lapacke_gtest_init_scomplex_buffer_pair_rand( c, cref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
trsyl_scomplex_parameters :: ~trsyl_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trsyl_lapack_complex_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trsyl_free();

}
/*  Test fixture class definition */
class ctrsyl_test  : public  ::testing::Test {
public:
   trsyl_scomplex_parameters  *ctrsyl_obj;
   void SetUp();
   void TearDown () { delete ctrsyl_obj; }
};

void ctrsyl_test::SetUp(){

    /* LAPACKE ctrsyl prototype */
    typedef int (*Fptr_NL_LAPACKE_ctrsyl) (int matrix_layout, char transa, char transb, lapack_int isgn,\
											lapack_int m, lapack_int n, const lapack_complex_float *A, lapack_int lda, \ 
											const lapack_complex_float* b, lapack_int ldb, const lapack_complex_float* c, lapack_int ldc, float* scale);

    Fptr_NL_LAPACKE_ctrsyl ctrsyl;

    ctrsyl_obj = new trsyl_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].isgn,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
						   eig_paramslist[idx].lda,
						   eig_paramslist[idx].ldb,
						   eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    ctrsyl_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctrsyl_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctrsyl_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctrsyl_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ctrsyl = (Fptr_NL_LAPACKE_ctrsyl)dlsym(ctrsyl_obj->hModule, "LAPACKE_ctrsyl");
    ASSERT_TRUE(ctrsyl != NULL) << "failed to get the Netlib LAPACKE_ctrsyl symbol";
    

    ctrsyl_obj->inforef = ctrsyl( ctrsyl_obj->matrix_layout,ctrsyl_obj->transa, ctrsyl_obj->transb,\
								ctrsyl_obj->isgn, ctrsyl_obj->m, ctrsyl_obj->n, ctrsyl_obj->Aref, \
								ctrsyl_obj->lda, ctrsyl_obj->bref, ctrsyl_obj->ldb, ctrsyl_obj->cref,\
								ctrsyl_obj->ldc, ctrsyl_obj->scaleref);

    /* Compute libflame's Lapacke o/p  */
    ctrsyl_obj->info = LAPACKE_ctrsyl( ctrsyl_obj->matrix_layout,ctrsyl_obj->transa, ctrsyl_obj->transb,\
								ctrsyl_obj->isgn, ctrsyl_obj->m, ctrsyl_obj->n, ctrsyl_obj->A, \
								ctrsyl_obj->lda, ctrsyl_obj->b, ctrsyl_obj->ldb, ctrsyl_obj->c,\
								ctrsyl_obj->ldc, ctrsyl_obj->scale);

    if( ctrsyl_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctrsyl is wrong\n", ctrsyl_obj->info );
    }
    if( ctrsyl_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctrsyl is wrong\n", 
        ctrsyl_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctrsyl_obj->diff =  computeDiff_c( ctrsyl_obj->bufsize, 
                ctrsyl_obj->c, ctrsyl_obj->cref );

}

TEST_F(ctrsyl_test, ctrsyl1) {
    EXPECT_NEAR(0.0, ctrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctrsyl_test, ctrsyl2) {
    EXPECT_NEAR(0.0, ctrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctrsyl_test, ctrsyl3) {
    EXPECT_NEAR(0.0, ctrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctrsyl_test, ctrsyl4) {
    EXPECT_NEAR(0.0, ctrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class trsyl_dcomplex_parameters{

   public:
	int bufsize, bufsize_a, bufsize_b;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char transa;
	char transb;
	lapack_complex_double* A;
	lapack_complex_double* b;
	lapack_complex_double* c;	
	lapack_int lda;
	lapack_int ldb;
	lapack_int ldc;
	lapack_int isgn;	
	/*Output Parameter*/
	double* scale, *scaleref;
	lapack_complex_double *Aref;
	lapack_complex_double *bref, *cref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      trsyl_dcomplex_parameters (int matrix_layout, char transa, char transb, lapack_int isgn,\
							lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb, lapack_int ldc);
      ~trsyl_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
trsyl_dcomplex_parameters:: trsyl_dcomplex_parameters (int matrix_layout_i, char transa_i, char transb_i, lapack_int isgn_i,\
													lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i, lapack_int ldc_i)
{
	
	matrix_layout = matrix_layout_i;
	transa = transa_i;
	transb = transb_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	ldb = ldb_i;
	ldc = ldc_i;
	isgn = isgn_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n trsyl lapack_complex_double:  m: %d, n: %d lda: %d, ldb:%d, ldc:%d, transa:%c \n",  m, n, lda, ldb,ldc, transa);
	#endif
	
	if ((transa == 'T') ||(transb == 'T')) // trans = 'T' supports  real flavours only
	{
		transa = transb = 'C';
	}
	
	
	if (matrix_layout == LAPACK_COL_MAJOR)
		bufsize = ldc*n;
	else if (matrix_layout == LAPACK_ROW_MAJOR)
		bufsize = ldc*m;
	else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	bufsize_a = lda*m;
	bufsize_b = ldb*n;
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&c, &cref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&scale, &scaleref, n);
	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL) || \
		(c == NULL) || (cref == NULL) ||
		(scale == NULL) || (scaleref == NULL)){
		EXPECT_FALSE( true) << "trsyl_lapack_complex_double_parameters object: malloc error.";
		trsyl_free();
		exit(0);
	}
	/* Initialization of input matrices as quasi-triangular matrix */
	lapacke_gtest_init_dcomplex_buffer_pair_with_constant(A, Aref, bufsize_a, 0);
	lapacke_gtest_init_dcomplex_buffer_pair_with_constant(b, bref, bufsize_b, 0);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( c, cref, bufsize);
	

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
trsyl_dcomplex_parameters :: ~trsyl_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" trsyl_lapack_complex_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   trsyl_free();

}
/*  Test fixture class definition */
class ztrsyl_test  : public  ::testing::Test {
public:
   trsyl_dcomplex_parameters  *ztrsyl_obj;
   void SetUp();
   void TearDown () { delete ztrsyl_obj; }
};

void ztrsyl_test::SetUp(){

    /* LAPACKE ztrsyl prototype */
    typedef int (*Fptr_NL_LAPACKE_ztrsyl) (int matrix_layout, char transa, char transb, lapack_int isgn,\
											lapack_int m, lapack_int n, const lapack_complex_double *A, lapack_int lda, \ 
											const lapack_complex_double* b, lapack_int ldb, const lapack_complex_double* c, lapack_int ldc, double* scale);

    Fptr_NL_LAPACKE_ztrsyl ztrsyl;

    ztrsyl_obj = new trsyl_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].trans,
                           eig_paramslist[idx].isgn,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
						   eig_paramslist[idx].lda,
						   eig_paramslist[idx].ldb,
						   eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);
	

    ztrsyl_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztrsyl_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztrsyl_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztrsyl_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ztrsyl = (Fptr_NL_LAPACKE_ztrsyl)dlsym(ztrsyl_obj->hModule, "LAPACKE_ztrsyl");
    ASSERT_TRUE(ztrsyl != NULL) << "failed to get the Netlib LAPACKE_ztrsyl symbol";
    

    ztrsyl_obj->inforef = ztrsyl( ztrsyl_obj->matrix_layout,ztrsyl_obj->transa, ztrsyl_obj->transb,\
								ztrsyl_obj->isgn, ztrsyl_obj->m, ztrsyl_obj->n, ztrsyl_obj->Aref, \
								ztrsyl_obj->lda, ztrsyl_obj->bref, ztrsyl_obj->ldb, ztrsyl_obj->cref,\
								ztrsyl_obj->ldc, ztrsyl_obj->scaleref);

    /* Compute libflame's Lapacke o/p  */
    ztrsyl_obj->info = LAPACKE_ztrsyl( ztrsyl_obj->matrix_layout,ztrsyl_obj->transa, ztrsyl_obj->transb,\
								ztrsyl_obj->isgn, ztrsyl_obj->m, ztrsyl_obj->n, ztrsyl_obj->A, \
								ztrsyl_obj->lda, ztrsyl_obj->b, ztrsyl_obj->ldb, ztrsyl_obj->c,\
								ztrsyl_obj->ldc, ztrsyl_obj->scale);

    if( ztrsyl_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztrsyl is wrong\n", ztrsyl_obj->info );
    }
    if( ztrsyl_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztrsyl is wrong\n", 
        ztrsyl_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztrsyl_obj->diff =  computeDiff_z( ztrsyl_obj->bufsize, 
                ztrsyl_obj->c, ztrsyl_obj->cref );

}

TEST_F(ztrsyl_test, ztrsyl1) {
    EXPECT_NEAR(0.0, ztrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztrsyl_test, ztrsyl2) {
    EXPECT_NEAR(0.0, ztrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztrsyl_test, ztrsyl3) {
    EXPECT_NEAR(0.0, ztrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztrsyl_test, ztrsyl4) {
    EXPECT_NEAR(0.0, ztrsyl_obj->diff, LAPACKE_EIG_THRESHOLD);
}
