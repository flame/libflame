#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define lange_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
class lange_float_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char norm;
	float* A;
	lapack_int lda;
	/*Output Parameter*/
	float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lange_float_parameters (int matrix_layout, char norm, lapack_int m , lapack_int n, lapack_int lda);
      ~lange_float_parameters ();

};

/* Constructor definition  float_common_parameters */
lange_float_parameters::lange_float_parameters (int matrix_layout_i,char norm_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	norm = norm_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lange float:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, norm);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lange_float_parameters object: malloc error.";
		lange_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lange_float_parameters :: ~lange_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lange_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lange_free();

}
/*  Test fixture class definition */
class slange_test  : public  ::testing::Test {
public:
   lange_float_parameters  *slange_obj;
   void SetUp();
   void TearDown () { delete slange_obj; }
};

void slange_test::SetUp(){

    /* LAPACKE slange prototype */
    typedef int (*Fptr_NL_LAPACKE_slange) (int matrix_layout, char norm, lapack_int m,lapack_int n, const float *A, lapack_int lda);

    Fptr_NL_LAPACKE_slange slange;

    slange_obj = new lange_float_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    slange_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slange_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slange_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slange_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    slange = (Fptr_NL_LAPACKE_slange)dlsym(slange_obj->hModule, "LAPACKE_slange");
    ASSERT_TRUE(slange != NULL) << "failed to get the Netlib LAPACKE_slange symbol";
    

    slange_obj->inforef = slange( slange_obj->matrix_layout, slange_obj->norm,
								slange_obj->m, slange_obj->n, (const float *)slange_obj->Aref,
								slange_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    slange_obj->info = LAPACKE_slange( slange_obj->matrix_layout, slange_obj->norm,
										slange_obj->m, slange_obj->n, (const float*)slange_obj->A, 
										slange_obj->lda);

    if( slange_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slange is wrong\n", slange_obj->info );
    }
    if( slange_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slange is wrong\n", 
        slange_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slange_obj->diff =  computeDiff_s( slange_obj->bufsize, 
                slange_obj->A, slange_obj->Aref );

}

TEST_F(slange_test, slange1) {
    EXPECT_NEAR(0.0, slange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slange_test, slange2) {
    EXPECT_NEAR(0.0, slange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slange_test, slange3) {
    EXPECT_NEAR(0.0, slange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slange_test, slange4) {
    EXPECT_NEAR(0.0, slange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class lange_double_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char norm;
	double* A;
	lapack_int lda;
	/*Output Parameter*/
	double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lange_double_parameters (int matrix_layout, char norm, lapack_int m , lapack_int n, lapack_int lda);
      ~lange_double_parameters ();

};

/* Constructor definition  double_common_parameters */
lange_double_parameters::lange_double_parameters (int matrix_layout_i,char norm_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	norm = norm_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lange double:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, norm);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lange_double_parameters object: malloc error.";
		lange_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lange_double_parameters :: ~lange_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lange_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lange_free();

}
/*  Test fixture class definition */
class dlange_test  : public  ::testing::Test {
public:
   lange_double_parameters  *dlange_obj;
   void SetUp();
   void TearDown () { delete dlange_obj; }
};

void dlange_test::SetUp(){

    /* LAPACKE dlange prototype */
    typedef int (*Fptr_NL_LAPACKE_dlange) (int matrix_layout, char norm, lapack_int m,lapack_int n, const double *A, lapack_int lda);

    Fptr_NL_LAPACKE_dlange dlange;

    dlange_obj = new lange_double_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    dlange_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlange_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlange_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlange_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dlange = (Fptr_NL_LAPACKE_dlange)dlsym(dlange_obj->hModule, "LAPACKE_dlange");
    ASSERT_TRUE(dlange != NULL) << "failed to get the Netlib LAPACKE_dlange symbol";
    

    dlange_obj->inforef = dlange( dlange_obj->matrix_layout, dlange_obj->norm,
								dlange_obj->m, dlange_obj->n, (const double *)dlange_obj->Aref,
								dlange_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    dlange_obj->info = LAPACKE_dlange( dlange_obj->matrix_layout, dlange_obj->norm,
										dlange_obj->m, dlange_obj->n, (const double*)dlange_obj->A, 
										dlange_obj->lda);

    if( dlange_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlange is wrong\n", dlange_obj->info );
    }
    if( dlange_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlange is wrong\n", 
        dlange_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlange_obj->diff =  computeDiff_d( dlange_obj->bufsize, 
                dlange_obj->A, dlange_obj->Aref );

}

TEST_F(dlange_test, dlange1) {
    EXPECT_NEAR(0.0, dlange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlange_test, dlange2) {
    EXPECT_NEAR(0.0, dlange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlange_test, dlange3) {
    EXPECT_NEAR(0.0, dlange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlange_test, dlange4) {
    EXPECT_NEAR(0.0, dlange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class lange_scomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char norm;
	lapack_complex_float* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lange_scomplex_parameters (int matrix_layout, char norm, lapack_int m , lapack_int n, lapack_int lda);
      ~lange_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
lange_scomplex_parameters::lange_scomplex_parameters (int matrix_layout_i,char norm_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	norm = norm_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lange scomplex:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, norm );
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lange_float_parameters object: malloc error.";
		lange_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
lange_scomplex_parameters :: ~lange_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lange_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lange_free();

}
/*  Test fixture class definition */
class clange_test  : public  ::testing::Test {
public:
   lange_scomplex_parameters  *clange_obj;
   void SetUp();
   void TearDown () { delete clange_obj; }
};

void clange_test::SetUp(){

    /* LAPACKE clange prototype */
    typedef int (*Fptr_NL_LAPACKE_clange) (int matrix_layout, char norm, lapack_int m,lapack_int n, const lapack_complex_float *A, lapack_int lda);

    Fptr_NL_LAPACKE_clange clange;

    clange_obj = new lange_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    clange_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    clange_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(clange_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(clange_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    clange = (Fptr_NL_LAPACKE_clange)dlsym(clange_obj->hModule, "LAPACKE_clange");
    ASSERT_TRUE(clange != NULL) << "failed to get the Netlib LAPACKE_clange symbol";
    

    clange_obj->inforef = clange( clange_obj->matrix_layout, clange_obj->norm,
								clange_obj->m, clange_obj->n, (const lapack_complex_float *)clange_obj->Aref,
								clange_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    clange_obj->info = LAPACKE_clange( clange_obj->matrix_layout, clange_obj->norm,
										clange_obj->m, clange_obj->n, (const lapack_complex_float*)clange_obj->A, 
										clange_obj->lda);

    if( clange_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_clange is wrong\n", clange_obj->info );
    }
    if( clange_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_clange is wrong\n", 
        clange_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    clange_obj->diff =  computeDiff_c( clange_obj->bufsize, 
                clange_obj->A, clange_obj->Aref );

}

TEST_F(clange_test, clange1) {
    EXPECT_NEAR(0.0, clange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clange_test, clange2) {
    EXPECT_NEAR(0.0, clange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clange_test, clange3) {
    EXPECT_NEAR(0.0, clange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(clange_test, clange4) {
    EXPECT_NEAR(0.0, clange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class lange_dcomplex_parameters{

   public:
	int bufsize;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char norm;
	lapack_complex_double* A;
	lapack_int lda;
	/*Output Parameter*/
	lapack_complex_double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      lange_dcomplex_parameters (int matrix_layout, char norm, lapack_int m , lapack_int n, lapack_int lda);
      ~lange_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
lange_dcomplex_parameters::lange_dcomplex_parameters (int matrix_layout_i,char norm_i, lapack_int m_i, lapack_int n_i, lapack_int lda_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	norm = norm_i;
	lda = lda_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n lange dcomplex:  m: %d, n: %d lda: %d, uplo:%c, ldb:%d\n",  m, n, lda, norm );
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	lda = m;
		bufsize = lda*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	lda = n;
		bufsize = lda*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "lange_double_parameters object: malloc error.";
		lange_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
lange_dcomplex_parameters :: ~lange_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" lange_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   lange_free();

}
/*  Test fixture class definition */
class zlange_test  : public  ::testing::Test {
public:
   lange_dcomplex_parameters  *zlange_obj;
   void SetUp();
   void TearDown () { delete zlange_obj; }
};

void zlange_test::SetUp(){

    /* LAPACKE zlange prototype */
    typedef int (*Fptr_NL_LAPACKE_zlange) (int matrix_layout, char norm, lapack_int m,lapack_int n, const lapack_complex_double *A, lapack_int lda);

    Fptr_NL_LAPACKE_zlange zlange;

    zlange_obj = new lange_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						eig_paramslist[idx].norm,
						eig_paramslist[idx].m,
						eig_paramslist[idx].n,
						eig_paramslist[idx].lda);

    idx = Circular_Increment_Index(idx);

    zlange_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlange_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlange_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlange_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    zlange = (Fptr_NL_LAPACKE_zlange)dlsym(zlange_obj->hModule, "LAPACKE_zlange");
    ASSERT_TRUE(zlange != NULL) << "failed to get the Netlib LAPACKE_zlange symbol";
    

    zlange_obj->inforef = zlange( zlange_obj->matrix_layout, zlange_obj->norm,
								zlange_obj->m, zlange_obj->n, (const lapack_complex_double *)zlange_obj->Aref,
								zlange_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    zlange_obj->info = LAPACKE_zlange( zlange_obj->matrix_layout, zlange_obj->norm,
										zlange_obj->m, zlange_obj->n, (const lapack_complex_double*)zlange_obj->A, 
										zlange_obj->lda);

    if( zlange_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlange is wrong\n", zlange_obj->info );
    }
    if( zlange_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlange is wrong\n", 
        zlange_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlange_obj->diff =  computeDiff_z( zlange_obj->bufsize, 
                zlange_obj->A, zlange_obj->Aref );

}

TEST_F(zlange_test, zlange1) {
    EXPECT_NEAR(0.0, zlange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlange_test, zlange2) {
    EXPECT_NEAR(0.0, zlange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlange_test, zlange3) {
    EXPECT_NEAR(0.0, zlange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlange_test, zlange4) {
    EXPECT_NEAR(0.0, zlange_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
