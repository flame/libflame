#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define lacpy_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)  free(b);\
if (bref!=NULL) free(bref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
class lacpy_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	float* A;
	lapack_int lda, ldb;
	/*Output Parameter*/
	float* b;
	float *Aref;
	float *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacpy_float_parameters (int matrix_layout, char uplo, lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb);
      ~lacpy_float_parameters ();

};

/* Constructor definition  float_common_parameters */
lacpy_float_parameters::lacpy_float_parameters (int matrix_layout_i,char uplo_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	uplo = uplo_i;
	lda = lda_i;
	ldb = ldb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lacpy float:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, uplo, ldb);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		ldb = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		ldb = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize);	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "lacpy_float_parameters object: malloc error.";
		lacpy_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lacpy_float_parameters :: ~lacpy_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacpy_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacpy_free();

}
/*  Test fixture class definition */
class slacpy_test  : public  ::testing::Test {
public:
   lacpy_float_parameters  *slacpy_obj;
   void SetUp();
   void TearDown () { delete slacpy_obj; }
};

void slacpy_test::SetUp(){

    /* LAPACKE slacpy prototype */
    typedef int (*Fptr_NL_LAPACKE_slacpy) (int matrix_layout, char uplo, lapack_int m,lapack_int n, const float *A, lapack_int lda,  float* b, lapack_int ldb);

    Fptr_NL_LAPACKE_slacpy slacpy;

    slacpy_obj = new lacpy_float_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda,
						eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    slacpy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slacpy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slacpy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slacpy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slacpy = (Fptr_NL_LAPACKE_slacpy)dlsym(slacpy_obj->hModule, "LAPACKE_slacpy");
    ASSERT_TRUE(slacpy != NULL) << "failed to get the Netlib LAPACKE_slacpy symbol";
    

    slacpy_obj->inforef = slacpy( slacpy_obj->matrix_layout, slacpy_obj->uplo,
								slacpy_obj->m, slacpy_obj->n, (const float *)slacpy_obj->Aref,
								slacpy_obj->lda, slacpy_obj->bref, slacpy_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    slacpy_obj->info = LAPACKE_slacpy( slacpy_obj->matrix_layout, slacpy_obj->uplo,
										slacpy_obj->m, slacpy_obj->n, (const float*)slacpy_obj->A, 
										slacpy_obj->lda, slacpy_obj->b, slacpy_obj->ldb);

    if( slacpy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slacpy is wrong\n", slacpy_obj->info );
    }
    if( slacpy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slacpy is wrong\n", 
        slacpy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slacpy_obj->diff =  computeDiff_s( slacpy_obj->bufsize, 
                slacpy_obj->b, slacpy_obj->bref );

}

TEST_F(slacpy_test, slacpy1) {
    EXPECT_NEAR(0.0, slacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slacpy_test, slacpy2) {
    EXPECT_NEAR(0.0, slacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slacpy_test, slacpy3) {
    EXPECT_NEAR(0.0, slacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slacpy_test, slacpy4) {
    EXPECT_NEAR(0.0, slacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class lacpy_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	double* A;
	lapack_int lda, ldb;
	/*Output Parameter*/
	double* b;
	double *Aref;
	double *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacpy_double_parameters (int matrix_layout, char uplo, lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb);
      ~lacpy_double_parameters ();

};

/* Constructor definition  double_common_parameters */
lacpy_double_parameters::lacpy_double_parameters (int matrix_layout_i,char uplo_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	uplo = uplo_i;
	lda = lda_i;
	ldb = ldb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lacpy double:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, uplo, ldb);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		ldb = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		ldb = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize);	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "lacpy_double_parameters object: malloc error.";
		lacpy_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lacpy_double_parameters :: ~lacpy_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacpy_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacpy_free();

}
/*  Test fixture class definition */
class dlacpy_test  : public  ::testing::Test {
public:
   lacpy_double_parameters  *dlacpy_obj;
   void SetUp();
   void TearDown () { delete dlacpy_obj; }
};

void dlacpy_test::SetUp(){

    /* LAPACKE dlacpy prototype */
    typedef int (*Fptr_NL_LAPACKE_dlacpy) (int matrix_layout, char uplo, lapack_int m,lapack_int n, const double *A, lapack_int lda,  double* b, lapack_int ldb);

    Fptr_NL_LAPACKE_dlacpy dlacpy;

    dlacpy_obj = new lacpy_double_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda,
						eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    dlacpy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlacpy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlacpy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlacpy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlacpy = (Fptr_NL_LAPACKE_dlacpy)dlsym(dlacpy_obj->hModule, "LAPACKE_dlacpy");
    ASSERT_TRUE(dlacpy != NULL) << "failed to get the Netlib LAPACKE_dlacpy symbol";
    

    dlacpy_obj->inforef = dlacpy( dlacpy_obj->matrix_layout, dlacpy_obj->uplo,
								dlacpy_obj->m, dlacpy_obj->n, (const double *)dlacpy_obj->Aref,
								dlacpy_obj->lda, dlacpy_obj->bref, dlacpy_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    dlacpy_obj->info = LAPACKE_dlacpy( dlacpy_obj->matrix_layout, dlacpy_obj->uplo,
										dlacpy_obj->m, dlacpy_obj->n, (const double*)dlacpy_obj->A, 
										dlacpy_obj->lda, dlacpy_obj->b, dlacpy_obj->ldb);

    if( dlacpy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlacpy is wrong\n", dlacpy_obj->info );
    }
    if( dlacpy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlacpy is wrong\n", 
        dlacpy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlacpy_obj->diff =  computeDiff_d( dlacpy_obj->bufsize, 
                dlacpy_obj->b, dlacpy_obj->bref );

}

TEST_F(dlacpy_test, dlacpy1) {
    EXPECT_NEAR(0.0, dlacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlacpy_test, dlacpy2) {
    EXPECT_NEAR(0.0, dlacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlacpy_test, dlacpy3) {
    EXPECT_NEAR(0.0, dlacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlacpy_test, dlacpy4) {
    EXPECT_NEAR(0.0, dlacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class lacpy_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	lapack_complex_float* A;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_float* b;
	lapack_complex_float *Aref;
	lapack_complex_float *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacpy_scomplex_parameters (int matrix_layout, char uplo, lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb);
      ~lacpy_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lacpy_scomplex_parameters::lacpy_scomplex_parameters (int matrix_layout_i,char uplo_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	uplo = uplo_i;
	lda = lda_i;
	ldb = ldb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lacpy scomplex:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, uplo, ldb);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		ldb = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		ldb = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize);	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "lacpy_float_parameters object: malloc error.";
		lacpy_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lacpy_scomplex_parameters :: ~lacpy_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacpy_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacpy_free();

}
/*  Test fixture class definition */
class clacpy_test  : public  ::testing::Test {
public:
   lacpy_scomplex_parameters  *clacpy_obj;
   void SetUp();
   void TearDown () { delete clacpy_obj; }
};

void clacpy_test::SetUp(){

    /* LAPACKE clacpy prototype */
    typedef int (*Fptr_NL_LAPACKE_clacpy) (int matrix_layout, char uplo, lapack_int m,lapack_int n, const lapack_complex_float *A, lapack_int lda,  lapack_complex_float* b, lapack_int ldb);

    Fptr_NL_LAPACKE_clacpy clacpy;

    clacpy_obj = new lacpy_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda,
						eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    clacpy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clacpy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clacpy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clacpy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clacpy = (Fptr_NL_LAPACKE_clacpy)dlsym(clacpy_obj->hModule, "LAPACKE_clacpy");
    ASSERT_TRUE(clacpy != NULL) << "failed to get the Netlib LAPACKE_clacpy symbol";
    

    clacpy_obj->inforef = clacpy( clacpy_obj->matrix_layout, clacpy_obj->uplo,
								clacpy_obj->m, clacpy_obj->n, (const lapack_complex_float *)clacpy_obj->Aref,
								clacpy_obj->lda, clacpy_obj->bref, clacpy_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    clacpy_obj->info = LAPACKE_clacpy( clacpy_obj->matrix_layout, clacpy_obj->uplo,
										clacpy_obj->m, clacpy_obj->n, (const lapack_complex_float*)clacpy_obj->A, 
										clacpy_obj->lda, clacpy_obj->b, clacpy_obj->ldb);

    if( clacpy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clacpy is wrong\n", clacpy_obj->info );
    }
    if( clacpy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clacpy is wrong\n", 
        clacpy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clacpy_obj->diff =  computeDiff_c( clacpy_obj->bufsize, 
                clacpy_obj->b, clacpy_obj->bref );

}

TEST_F(clacpy_test, clacpy1) {
    EXPECT_NEAR(0.0, clacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacpy_test, clacpy2) {
    EXPECT_NEAR(0.0, clacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacpy_test, clacpy3) {
    EXPECT_NEAR(0.0, clacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clacpy_test, clacpy4) {
    EXPECT_NEAR(0.0, clacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lacpy_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	lapack_complex_double* A;
	lapack_int lda, ldb;
	/*Output Parameter*/
	lapack_complex_double* b;
	lapack_complex_double *Aref;
	lapack_complex_double *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lacpy_dcomplex_parameters (int matrix_layout, char uplo, lapack_int m , lapack_int n, lapack_int lda, lapack_int ldb);
      ~lacpy_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lacpy_dcomplex_parameters::lacpy_dcomplex_parameters (int matrix_layout_i,char uplo_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i, lapack_int ldb_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	uplo = uplo_i;
	lda = lda_i;
	ldb = ldb_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lacpy dcomplex:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, uplo, ldb);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		ldb = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		ldb = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize);	
	if ((A==NULL) || (Aref==NULL) || \
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "lacpy_double_parameters object: malloc error.";
		lacpy_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lacpy_dcomplex_parameters :: ~lacpy_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lacpy_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lacpy_free();

}
/*  Test fixture class definition */
class zlacpy_test  : public  ::testing::Test {
public:
   lacpy_dcomplex_parameters  *zlacpy_obj;
   void SetUp();
   void TearDown () { delete zlacpy_obj; }
};

void zlacpy_test::SetUp(){

    /* LAPACKE zlacpy prototype */
    typedef int (*Fptr_NL_LAPACKE_zlacpy) (int matrix_layout, char uplo, lapack_int m,lapack_int n, const lapack_complex_double *A, lapack_int lda,  lapack_complex_double* b, lapack_int ldb);

    Fptr_NL_LAPACKE_zlacpy zlacpy;

    zlacpy_obj = new lacpy_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda,
						eig_paramslist[idx].ldb);

    idx = Circular_Increment_Index(idx);

    zlacpy_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlacpy_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlacpy_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlacpy_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlacpy = (Fptr_NL_LAPACKE_zlacpy)dlsym(zlacpy_obj->hModule, "LAPACKE_zlacpy");
    ASSERT_TRUE(zlacpy != NULL) << "failed to get the Netlib LAPACKE_zlacpy symbol";
    

    zlacpy_obj->inforef = zlacpy( zlacpy_obj->matrix_layout, zlacpy_obj->uplo,
								zlacpy_obj->m, zlacpy_obj->n, (const lapack_complex_double *)zlacpy_obj->Aref,
								zlacpy_obj->lda, zlacpy_obj->bref, zlacpy_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    zlacpy_obj->info = LAPACKE_zlacpy( zlacpy_obj->matrix_layout, zlacpy_obj->uplo,
										zlacpy_obj->m, zlacpy_obj->n, (const lapack_complex_double*)zlacpy_obj->A, 
										zlacpy_obj->lda, zlacpy_obj->b, zlacpy_obj->ldb);

    if( zlacpy_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlacpy is wrong\n", zlacpy_obj->info );
    }
    if( zlacpy_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlacpy is wrong\n", 
        zlacpy_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlacpy_obj->diff =  computeDiff_z( zlacpy_obj->bufsize, 
                zlacpy_obj->b, zlacpy_obj->bref );

}

TEST_F(zlacpy_test, zlacpy1) {
    EXPECT_NEAR(0.0, zlacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacpy_test, zlacpy2) {
    EXPECT_NEAR(0.0, zlacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacpy_test, zlacpy3) {
    EXPECT_NEAR(0.0, zlacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlacpy_test, zlacpy4) {
    EXPECT_NEAR(0.0, zlacpy_obj->diff, LAPACKE_GTEST_THRESHOLD);
}