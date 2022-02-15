#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define tpqrt2_free() \
if (a!=NULL)    free(a); \
if (aref!=NULL) free(aref);\
if (t!=NULL)  free(t);\
if (tref!=NULL) free(tref); \
if (b != NULL) free(b); \
if (bref != NULL) free(bref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class tpqrt2_float_parameters{

   public:
	int bufsize_t;
	int bufsize_b;
	void *hModule, *dModule;
	float diff;
	float diff_b;
	float diff_t;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int l;
	lapack_int lda, ldb, ldt;
	float* a;
	/*Output Parameter*/
	float* t, *b;
	float *aref, *tref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      tpqrt2_float_parameters (int matrix_layout,  lapack_int m, lapack_int n, lapack_int l);
      ~tpqrt2_float_parameters ();

};

/* Constructor definition  float_common_parameters */
tpqrt2_float_parameters:: tpqrt2_float_parameters (int matrix_layout_i , lapack_int m_i, lapack_int n_i, lapack_int l_i)
{
	
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	l = l_i;
	ldt = lda = n;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpqrt2 float: matrix_layout = %d,  m:%d, n: %d, l%d \n", matrix_layout,  m, n, l);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	
		ldb = m;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		ldb = n;
		bufsize_b = ldb*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*n;
	l = min(m,n);

	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&a, &aref, bufsize_t);
	lapacke_gtest_alloc_float_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	
	if ((a==NULL) || (aref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(b==NULL) || (bref==NULL))
	{
		EXPECT_FALSE( true) << "tpqrt2_float_parameters object: malloc error.";
		tpqrt2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( a, aref, bufsize_t);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tpqrt2_float_parameters :: ~tpqrt2_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpqrt2_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpqrt2_free();

}
/*  Test fixture class definition */
class stpqrt2_test  : public  ::testing::Test {
public:
   tpqrt2_float_parameters  *stpqrt2_obj;
   void SetUp();
   void TearDown () { delete stpqrt2_obj;}
};

void stpqrt2_test::SetUp(){
	 /* LAPACKE stpqrt2 prototype */
    typedef int (*Fptr_NL_LAPACKE_stpqrt2) (int matrix_layout, lapack_int m, lapack_int n, lapack_int l,\
	float * a, lapack_int lda, float * b, lapack_int ldb, float * t, lapack_int ldt);

    Fptr_NL_LAPACKE_stpqrt2 stpqrt2;

    stpqrt2_obj = new tpqrt2_float_parameters ( eig_paramslist[idx].matrix_layout,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    stpqrt2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stpqrt2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stpqrt2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stpqrt2_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*stpqrt2 library call */
    stpqrt2 = (Fptr_NL_LAPACKE_stpqrt2)dlsym(stpqrt2_obj->hModule, "LAPACKE_stpqrt2");
    ASSERT_TRUE(stpqrt2 != NULL) << "failed to get the Netlib LAPACKE_stpqrt2 symbol";
  
/*Compute stpqrt2's  o/p */
    stpqrt2_obj->inforef = stpqrt2( stpqrt2_obj->matrix_layout,	stpqrt2_obj->m, stpqrt2_obj->n,stpqrt2_obj->l, stpqrt2_obj->aref,\
	stpqrt2_obj->lda, stpqrt2_obj->bref, stpqrt2_obj->ldb,	stpqrt2_obj->tref, stpqrt2_obj->ldt);

    /* Compute libflame's Lapacke o/p  */
    stpqrt2_obj->info = LAPACKE_stpqrt2( stpqrt2_obj->matrix_layout, stpqrt2_obj->m, stpqrt2_obj->n,stpqrt2_obj->l, stpqrt2_obj->a,\
	stpqrt2_obj->lda, stpqrt2_obj->b, stpqrt2_obj->ldb,	stpqrt2_obj->t, stpqrt2_obj->ldt);
	
    if( stpqrt2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_stpqrt2 is wrong\n", stpqrt2_obj->info );
    }
    if( stpqrt2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stpqrt2 is wrong\n", 
        stpqrt2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stpqrt2_obj->diff =  computeDiff_s( stpqrt2_obj->bufsize_t, 
                stpqrt2_obj->a, stpqrt2_obj->aref );
	stpqrt2_obj->diff_b =  computeDiff_s( stpqrt2_obj->bufsize_b, 
                stpqrt2_obj->b, stpqrt2_obj->bref );
	stpqrt2_obj->diff_t =  computeDiff_s( stpqrt2_obj->bufsize_t, 
                stpqrt2_obj->t, stpqrt2_obj->tref );

}

TEST_F(stpqrt2_test, stpqrt21) {
    EXPECT_NEAR(0.0, stpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(stpqrt2_test, stpqrt22) {
    EXPECT_NEAR(0.0, stpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(stpqrt2_test, stpqrt23) {
    EXPECT_NEAR(0.0, stpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(stpqrt2_test, stpqrt24) {
    EXPECT_NEAR(0.0, stpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, stpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class tpqrt2_double_parameters{

   public:
	int bufsize_t;
	int bufsize_b;
	void *hModule, *dModule;
	double diff;
	double diff_b;
	double diff_t;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int l;
	lapack_int lda, ldb, ldt;
	double* a;
	/*Output Parameter*/
	double* t, *b;
	double *aref, *tref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      tpqrt2_double_parameters (int matrix_layout,  lapack_int m, lapack_int n, lapack_int l);
      ~tpqrt2_double_parameters ();

};

/* Constructor definition  double_common_parameters */
tpqrt2_double_parameters:: tpqrt2_double_parameters (int matrix_layout_i , lapack_int m_i, lapack_int n_i, lapack_int l_i)
{
	
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	l = l_i;
	ldt = lda = n;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpqrt2 double: matrix_layout = %d,  m:%d, n: %d, l%d \n", matrix_layout,  m, n, l);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	
		ldb = m;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		ldb = n;
		bufsize_b = ldb*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*n;
	l = min(m,n);

	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&a, &aref, bufsize_t);
	lapacke_gtest_alloc_double_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	
	if ((a==NULL) || (aref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(b==NULL) || (bref==NULL))
	{
		EXPECT_FALSE( true) << "tpqrt2_double_parameters object: malloc error.";
		tpqrt2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( a, aref, bufsize_t);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
tpqrt2_double_parameters :: ~tpqrt2_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpqrt2_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpqrt2_free();

}
/*  Test fixture class definition */
class dtpqrt2_test  : public  ::testing::Test {
public:
   tpqrt2_double_parameters  *dtpqrt2_obj;
   void SetUp();
   void TearDown () { delete dtpqrt2_obj;}
};

void dtpqrt2_test::SetUp(){
	 /* LAPACKE dtpqrt2 prototype */
    typedef int (*Fptr_NL_LAPACKE_dtpqrt2) (int matrix_layout, lapack_int m, lapack_int n, lapack_int l,\
	double * a, lapack_int lda, double * b, lapack_int ldb, double * t, lapack_int ldt);

    Fptr_NL_LAPACKE_dtpqrt2 dtpqrt2;

    dtpqrt2_obj = new tpqrt2_double_parameters ( eig_paramslist[idx].matrix_layout,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    dtpqrt2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtpqrt2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtpqrt2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtpqrt2_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dtpqrt2 library call */
    dtpqrt2 = (Fptr_NL_LAPACKE_dtpqrt2)dlsym(dtpqrt2_obj->hModule, "LAPACKE_dtpqrt2");
    ASSERT_TRUE(dtpqrt2 != NULL) << "failed to get the Netlib LAPACKE_dtpqrt2 symbol";
  
/*Compute dtpqrt2's  o/p */
    dtpqrt2_obj->inforef = dtpqrt2( dtpqrt2_obj->matrix_layout,	dtpqrt2_obj->m, dtpqrt2_obj->n,dtpqrt2_obj->l, dtpqrt2_obj->aref,\
	dtpqrt2_obj->lda, dtpqrt2_obj->bref, dtpqrt2_obj->ldb,	dtpqrt2_obj->tref, dtpqrt2_obj->ldt);

    /* Compute libflame's Lapacke o/p  */
    dtpqrt2_obj->info = LAPACKE_dtpqrt2( dtpqrt2_obj->matrix_layout, dtpqrt2_obj->m, dtpqrt2_obj->n,dtpqrt2_obj->l, dtpqrt2_obj->a,\
	dtpqrt2_obj->lda, dtpqrt2_obj->b, dtpqrt2_obj->ldb,	dtpqrt2_obj->t, dtpqrt2_obj->ldt);
	
    if( dtpqrt2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_dtpqrt2 is wrong\n", dtpqrt2_obj->info );
    }
    if( dtpqrt2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtpqrt2 is wrong\n", 
        dtpqrt2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtpqrt2_obj->diff =  computeDiff_d( dtpqrt2_obj->bufsize_t, 
                dtpqrt2_obj->a, dtpqrt2_obj->aref );
	dtpqrt2_obj->diff_b =  computeDiff_d( dtpqrt2_obj->bufsize_b, 
                dtpqrt2_obj->b, dtpqrt2_obj->bref );
	dtpqrt2_obj->diff_t =  computeDiff_d( dtpqrt2_obj->bufsize_t, 
                dtpqrt2_obj->t, dtpqrt2_obj->tref );

}

TEST_F(dtpqrt2_test, dtpqrt21) {
    EXPECT_NEAR(0.0, dtpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(dtpqrt2_test, dtpqrt22) {
    EXPECT_NEAR(0.0, dtpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dtpqrt2_test, dtpqrt23) {
    EXPECT_NEAR(0.0, dtpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}
TEST_F(dtpqrt2_test, dtpqrt24) {
    EXPECT_NEAR(0.0, dtpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, dtpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

/*Beging scomplex common parameter class definition*/
class tpqrt2_scomplex_parameters{

   public:
	int bufsize_t;
	int bufsize_b;
	void *hModule, *dModule;
	float diff;
	float diff_b;
	float diff_t;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int l;
	lapack_int lda, ldb, ldt;
	lapack_complex_float* a;
	/*Output Parameter*/
	lapack_complex_float* t, *b;
	lapack_complex_float *aref, *tref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      tpqrt2_scomplex_parameters (int matrix_layout,  lapack_int m, lapack_int n, lapack_int l);
      ~tpqrt2_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
tpqrt2_scomplex_parameters:: tpqrt2_scomplex_parameters (int matrix_layout_i , lapack_int m_i, lapack_int n_i, lapack_int l_i)
{
	
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	l = l_i;
	ldt = lda = n;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpqrt2 scomplex: matrix_layout = %d,  m:%d, n: %d, l%d \n", matrix_layout,  m, n, l);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	
		ldb = m;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		ldb = n;
		bufsize_b = ldb*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*n;
	l = min(m,n);

	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&a, &aref, bufsize_t);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	
	if ((a==NULL) || (aref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(b==NULL) || (bref==NULL))
	{
		EXPECT_FALSE( true) << "tpqrt2_scomplex_parameters object: malloc error.";
		tpqrt2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( a, aref, bufsize_t);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);

} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
tpqrt2_scomplex_parameters :: ~tpqrt2_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpqrt2_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpqrt2_free();

}
/*  Test fixture class definition */
class ctpqrt2_test  : public  ::testing::Test {
public:
   tpqrt2_scomplex_parameters  *ctpqrt2_obj;
   void SetUp();
   void TearDown () { delete ctpqrt2_obj;}
};

void ctpqrt2_test::SetUp(){
	 /* LAPACKE ctpqrt2 prototype */
    typedef int (*Fptr_NL_LAPACKE_ctpqrt2) (int matrix_layout, lapack_int m, lapack_int n, lapack_int l,\
	lapack_complex_float * a, lapack_int lda, lapack_complex_float * b, lapack_int ldb, lapack_complex_float * t, lapack_int ldt);

    Fptr_NL_LAPACKE_ctpqrt2 ctpqrt2;

    ctpqrt2_obj = new tpqrt2_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    ctpqrt2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctpqrt2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctpqrt2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctpqrt2_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ctpqrt2 library call */
    ctpqrt2 = (Fptr_NL_LAPACKE_ctpqrt2)dlsym(ctpqrt2_obj->hModule, "LAPACKE_ctpqrt2");
    ASSERT_TRUE(ctpqrt2 != NULL) << "failed to get the Netlib LAPACKE_ctpqrt2 symbol";
  
/*Compute ctpqrt2's  o/p */
    ctpqrt2_obj->inforef = ctpqrt2( ctpqrt2_obj->matrix_layout,	ctpqrt2_obj->m, ctpqrt2_obj->n,ctpqrt2_obj->l, ctpqrt2_obj->aref,\
	ctpqrt2_obj->lda, ctpqrt2_obj->bref, ctpqrt2_obj->ldb,	ctpqrt2_obj->tref, ctpqrt2_obj->ldt);

    /* Compute libflame's Lapacke o/p  */
    ctpqrt2_obj->info = LAPACKE_ctpqrt2( ctpqrt2_obj->matrix_layout, ctpqrt2_obj->m, ctpqrt2_obj->n,ctpqrt2_obj->l, ctpqrt2_obj->a,\
	ctpqrt2_obj->lda, ctpqrt2_obj->b, ctpqrt2_obj->ldb,	ctpqrt2_obj->t, ctpqrt2_obj->ldt);
	
    if( ctpqrt2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_ctpqrt2 is wrong\n", ctpqrt2_obj->info );
    }
    if( ctpqrt2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctpqrt2 is wrong\n", 
        ctpqrt2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctpqrt2_obj->diff =  computeDiff_c( ctpqrt2_obj->bufsize_t, 
                ctpqrt2_obj->a, ctpqrt2_obj->aref );
	ctpqrt2_obj->diff_b =  computeDiff_c( ctpqrt2_obj->bufsize_b, 
                ctpqrt2_obj->b, ctpqrt2_obj->bref );
	ctpqrt2_obj->diff_t =  computeDiff_c( ctpqrt2_obj->bufsize_t, 
                ctpqrt2_obj->t, ctpqrt2_obj->tref );

}

TEST_F(ctpqrt2_test, ctpqrt21) {
    EXPECT_NEAR(0.0, ctpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(ctpqrt2_test, ctpqrt22) {
    EXPECT_NEAR(0.0, ctpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctpqrt2_test, ctpqrt23) {
    EXPECT_NEAR(0.0, ctpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ctpqrt2_test, ctpqrt24) {
    EXPECT_NEAR(0.0, ctpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ctpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}
/*Beging scomplex common parameter class definition*/
class tpqrt2_dcomplex_parameters{

   public:
	int bufsize_t;
	int bufsize_b;
	void *hModule, *dModule;
	double diff;
	double diff_b;
	double diff_t;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	lapack_int l;
	lapack_int lda, ldb, ldt;
	lapack_complex_double* a;
	/*Output Parameter*/
	lapack_complex_double* t, *b;
	lapack_complex_double *aref, *tref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      tpqrt2_dcomplex_parameters (int matrix_layout,  lapack_int m, lapack_int n, lapack_int l);
      ~tpqrt2_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
tpqrt2_dcomplex_parameters:: tpqrt2_dcomplex_parameters (int matrix_layout_i , lapack_int m_i, lapack_int n_i, lapack_int l_i)
{
	
	matrix_layout = matrix_layout_i;
	m = m_i;
	n = n_i;
	l = l_i;
	ldt = lda = n;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tpqrt2 dcomplex: matrix_layout = %d,  m:%d, n: %d, l%d \n", matrix_layout,  m, n, l);
	#endif

	if (matrix_layout == LAPACK_COL_MAJOR)
	{	
		ldb = m;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		ldb = n;
		bufsize_b = ldb*m;
	}else
	{
		EXPECT_TRUE(false) << "matrix_layout invalid";
	}
	
	bufsize_t = ldt*n;
	l = min(m,n);
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&a, &aref, bufsize_t);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&t, &tref, bufsize_t);	
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	
	if ((a==NULL) || (aref==NULL) ||
		(t==NULL) || (tref==NULL) ||
		(b==NULL) || (bref==NULL))
	{
		EXPECT_FALSE( true) << "tpqrt2_dcomplex_parameters object: malloc error.";
		tpqrt2_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( a, aref, bufsize_t);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);

} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
tpqrt2_dcomplex_parameters :: ~tpqrt2_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tpqrt2_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tpqrt2_free();

}
/*  Test fixture class definition */
class ztpqrt2_test  : public  ::testing::Test {
public:
   tpqrt2_dcomplex_parameters  *ztpqrt2_obj;
   void SetUp();
   void TearDown () { delete ztpqrt2_obj;}
};

void ztpqrt2_test::SetUp(){
	 /* LAPACKE ztpqrt2 prototype */
    typedef int (*Fptr_NL_LAPACKE_ztpqrt2) (int matrix_layout, lapack_int m, lapack_int n, lapack_int l,\
	lapack_complex_double * a, lapack_int lda, lapack_complex_double * b, lapack_int ldb, lapack_complex_double * t, lapack_int ldt);

    Fptr_NL_LAPACKE_ztpqrt2 ztpqrt2;

    ztpqrt2_obj = new tpqrt2_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						   eig_paramslist[idx].m,
						   eig_paramslist[idx].n,
                           eig_paramslist[idx].n);
						   

    idx = Circular_Increment_Index(idx);

    ztpqrt2_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztpqrt2_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztpqrt2_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztpqrt2_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*ztpqrt2 library call */
    ztpqrt2 = (Fptr_NL_LAPACKE_ztpqrt2)dlsym(ztpqrt2_obj->hModule, "LAPACKE_ztpqrt2");
    ASSERT_TRUE(ztpqrt2 != NULL) << "failed to get the Netlib LAPACKE_ztpqrt2 symbol";
  
/*Compute ztpqrt2's  o/p */
    ztpqrt2_obj->inforef = ztpqrt2( ztpqrt2_obj->matrix_layout,	ztpqrt2_obj->m, ztpqrt2_obj->n,ztpqrt2_obj->l, ztpqrt2_obj->aref,\
	ztpqrt2_obj->lda, ztpqrt2_obj->bref, ztpqrt2_obj->ldb,	ztpqrt2_obj->tref, ztpqrt2_obj->ldt);

    /* Compute libflame's Lapacke o/p  */
    ztpqrt2_obj->info = LAPACKE_ztpqrt2( ztpqrt2_obj->matrix_layout, ztpqrt2_obj->m, ztpqrt2_obj->n,ztpqrt2_obj->l, ztpqrt2_obj->a,\
	ztpqrt2_obj->lda, ztpqrt2_obj->b, ztpqrt2_obj->ldb,	ztpqrt2_obj->t, ztpqrt2_obj->ldt);
	
    if( ztpqrt2_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame\
        LAPACKE_ztpqrt2 is wrong\n", ztpqrt2_obj->info );
    }
    if( ztpqrt2_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztpqrt2 is wrong\n", 
        ztpqrt2_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztpqrt2_obj->diff =  computeDiff_z( ztpqrt2_obj->bufsize_t, 
                ztpqrt2_obj->a, ztpqrt2_obj->aref );
	ztpqrt2_obj->diff_b =  computeDiff_z( ztpqrt2_obj->bufsize_b, 
                ztpqrt2_obj->b, ztpqrt2_obj->bref );
	ztpqrt2_obj->diff_t =  computeDiff_z( ztpqrt2_obj->bufsize_t, 
                ztpqrt2_obj->t, ztpqrt2_obj->tref );

}

TEST_F(ztpqrt2_test, ztpqrt21) {
    EXPECT_NEAR(0.0, ztpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
	
}

TEST_F(ztpqrt2_test, ztpqrt22) {
    EXPECT_NEAR(0.0, ztpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztpqrt2_test, ztpqrt23) {
    EXPECT_NEAR(0.0, ztpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(ztpqrt2_test, ztpqrt24) {
    EXPECT_NEAR(0.0, ztpqrt2_obj->diff, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt2_obj->diff_b, LAPACKE_GTEST_THRESHOLD);
	EXPECT_NEAR(0.0, ztpqrt2_obj->diff_t, LAPACKE_GTEST_THRESHOLD);
}