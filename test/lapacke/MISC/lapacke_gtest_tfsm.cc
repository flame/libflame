#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"


#define tfsm_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (b!=NULL)    free(b); \
if (bref!=NULL)    free(bref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
class tfsm_float_parameters{

   public:
	int bufsize;
	int bufsize_b;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	char diag;
	char trans, transr;
	char side;
	float* A;
	float alpha; 
	lapack_int  ldb=0;
	/*Output Parameter*/
	float *b;
	float *Aref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      tfsm_float_parameters (int matrix_layout, char transr, char side, char uplo,  char diag, lapack_int m ,lapack_int n, float alpha);
      ~tfsm_float_parameters ();

};

/* Constructor definition  float_common_parameters */
tfsm_float_parameters::tfsm_float_parameters (int matrix_layout_i, char transr_i, char side_i, char uplo_i, \ 
char diag_i, lapack_int m_i ,lapack_int n_i, float alpha_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	transr = transr_i;
	trans = transr;
	uplo = uplo_i;
	diag = diag_i; 
	alpha = 1;
	side = side_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfsm float:   n: %d, uplo:%c, trans:%c\n",  n, uplo, trans);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	ldb = m;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	ldb = n;
		bufsize_b = ldb*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	
	bufsize = (n*(n+1)/2);
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_float_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) ||
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "tfsm_float_parameters object: malloc error.";
		tfsm_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_float_buffer_pair_rand( b, bref, bufsize_b);
} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
tfsm_float_parameters :: ~tfsm_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfsm_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfsm_free();

}
/*  Test fixture class definition */
class stfsm_test  : public  ::testing::Test {
public:
   tfsm_float_parameters  *stfsm_obj;
   void SetUp();
   void TearDown () { delete stfsm_obj; }
};

void stfsm_test::SetUp(){

    /* LAPACKE stfsm prototype */
    typedef int (*Fptr_NL_LAPACKE_stfsm) (int matrix_layout , char transr , char side , char uplo , char trans ,\
	char diag , lapack_int m , lapack_int n , float alpha , const float * a , float * b , lapack_int ldb);

    Fptr_NL_LAPACKE_stfsm stfsm;

    stfsm_obj = new tfsm_float_parameters ( eig_paramslist[idx].matrix_layout,
						lin_solver_paramslist[idx].transr,
						eig_paramslist[idx].side,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].diag,
						eig_paramslist[idx].n,
						eig_paramslist[idx].m, 
						eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    stfsm_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    stfsm_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(stfsm_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(stfsm_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    stfsm = (Fptr_NL_LAPACKE_stfsm)dlsym(stfsm_obj->hModule, "LAPACKE_stfsm");
    ASSERT_TRUE(stfsm != NULL) << "failed to get the Netlib LAPACKE_stfsm symbol";	
	
    stfsm_obj->inforef = stfsm( stfsm_obj->matrix_layout, stfsm_obj->transr, stfsm_obj->side, stfsm_obj->uplo, 
								stfsm_obj->trans, stfsm_obj->diag, stfsm_obj->m, stfsm_obj->n, stfsm_obj->alpha, (const float *)stfsm_obj->Aref,
								stfsm_obj->bref, stfsm_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    stfsm_obj->info = LAPACKE_stfsm( stfsm_obj->matrix_layout, stfsm_obj->transr, stfsm_obj->side,stfsm_obj->uplo,
										stfsm_obj->trans, stfsm_obj->diag, stfsm_obj->m, stfsm_obj->n, stfsm_obj->alpha, (const float*)stfsm_obj->A, 
										stfsm_obj->b, stfsm_obj->ldb);

    if( stfsm_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_stfsm is wrong\n", stfsm_obj->info );
    }
    if( stfsm_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_stfsm is wrong\n", 
        stfsm_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    stfsm_obj->diff =  computeDiff_s( stfsm_obj->bufsize_b, 
                stfsm_obj->b, stfsm_obj->bref );

}

TEST_F(stfsm_test, stfsm1) {
    EXPECT_NEAR(0.0, stfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stfsm_test, stfsm2) {
    EXPECT_NEAR(0.0, stfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stfsm_test, stfsm3) {
    EXPECT_NEAR(0.0, stfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(stfsm_test, stfsm4) {
    EXPECT_NEAR(0.0, stfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class tfsm_double_parameters{

   public:
	int bufsize;
	int bufsize_b;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	char diag;
	char trans, transr;
	char side;
	double* A;
	double alpha; 
	lapack_int  ldb=0;
	/*Output Parameter*/
	double *b;
	double *Aref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      tfsm_double_parameters (int matrix_layout, char transr, char side, char uplo,  char diag, lapack_int m ,lapack_int n, double alpha);
      ~tfsm_double_parameters ();

};

/* Constructor definition  double_common_parameters */
tfsm_double_parameters::tfsm_double_parameters (int matrix_layout_i, char transr_i, char side_i, char uplo_i, \ 
char diag_i, lapack_int m_i ,lapack_int n_i, double alpha_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	transr = transr_i;
	trans = transr;
	uplo = uplo_i;
	diag = diag_i;
	alpha = 1; 
	side = side_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfsm double:   n: %d, uplo:%c, trans:%c\n",  n, uplo, trans);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	ldb = m;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	ldb = n;
		bufsize_b = ldb*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	
	bufsize = (n*(n+1)/2);
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_double_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) ||
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "tfsm_double_parameters object: malloc error.";
		tfsm_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_double_buffer_pair_rand( b, bref, bufsize_b);
} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
tfsm_double_parameters :: ~tfsm_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfsm_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfsm_free();

}
/*  Test fixture class definition */
class dtfsm_test  : public  ::testing::Test {
public:
   tfsm_double_parameters  *dtfsm_obj;
   void SetUp();
   void TearDown () { delete dtfsm_obj; }
};

void dtfsm_test::SetUp(){

    /* LAPACKE dtfsm prototype */
    typedef int (*Fptr_NL_LAPACKE_dtfsm) (int matrix_layout , char transr , char side , char uplo , char trans ,\
	char diag , lapack_int m , lapack_int n , double alpha , const double * a , double * b , lapack_int ldb);

    Fptr_NL_LAPACKE_dtfsm dtfsm;

    dtfsm_obj = new tfsm_double_parameters ( eig_paramslist[idx].matrix_layout,
						lin_solver_paramslist[idx].transr,
						eig_paramslist[idx].side,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].diag,
						eig_paramslist[idx].n,
						eig_paramslist[idx].m, 
						eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    dtfsm_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dtfsm_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dtfsm_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dtfsm_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    dtfsm = (Fptr_NL_LAPACKE_dtfsm)dlsym(dtfsm_obj->hModule, "LAPACKE_dtfsm");
    ASSERT_TRUE(dtfsm != NULL) << "failed to get the Netlib LAPACKE_dtfsm symbol";	
	
    dtfsm_obj->inforef = dtfsm( dtfsm_obj->matrix_layout, dtfsm_obj->transr, dtfsm_obj->side, dtfsm_obj->uplo, 
								dtfsm_obj->trans, dtfsm_obj->diag, dtfsm_obj->m, dtfsm_obj->n, dtfsm_obj->alpha, (const double *)dtfsm_obj->Aref,
								dtfsm_obj->bref, dtfsm_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    dtfsm_obj->info = LAPACKE_dtfsm( dtfsm_obj->matrix_layout, dtfsm_obj->transr, dtfsm_obj->side,dtfsm_obj->uplo,
										dtfsm_obj->trans, dtfsm_obj->diag, dtfsm_obj->m, dtfsm_obj->n, dtfsm_obj->alpha, (const double*)dtfsm_obj->A, 
										dtfsm_obj->b, dtfsm_obj->ldb);

    if( dtfsm_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dtfsm is wrong\n", dtfsm_obj->info );
    }
    if( dtfsm_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dtfsm is wrong\n", 
        dtfsm_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dtfsm_obj->diff =  computeDiff_d( dtfsm_obj->bufsize_b, 
                dtfsm_obj->b, dtfsm_obj->bref );

}

TEST_F(dtfsm_test, dtfsm1) {
    EXPECT_NEAR(0.0, dtfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtfsm_test, dtfsm2) {
    EXPECT_NEAR(0.0, dtfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtfsm_test, dtfsm3) {
    EXPECT_NEAR(0.0, dtfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(dtfsm_test, dtfsm4) {
    EXPECT_NEAR(0.0, dtfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class tfsm_scomplex_parameters{

   public:
	int bufsize;
	int bufsize_b;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	char diag;
	char trans, transr;
	char side;
	lapack_complex_float* A;
	lapack_complex_float alpha; 
	lapack_int  ldb=0;
	/*Output Parameter*/
	lapack_complex_float *b;
	lapack_complex_float *Aref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      tfsm_scomplex_parameters (int matrix_layout, char transr, char side, char uplo,  char diag, lapack_int m ,lapack_int n, lapack_complex_float alpha);
      ~tfsm_scomplex_parameters ();

};

/* Constructor definition  scomplex_common_parameters */
tfsm_scomplex_parameters::tfsm_scomplex_parameters (int matrix_layout_i, char transr_i, char side_i, char uplo_i, \ 
char diag_i, lapack_int m_i ,lapack_int n_i, lapack_complex_float alpha_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	transr = transr_i;
	trans = transr;
	uplo = uplo_i;
	diag = diag_i;
	alpha = 1;
	side = side_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfsm lapack_complex_float:   n: %d, uplo:%c, trans:%c\n",  n, uplo, trans);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	ldb = m;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	ldb = n;
		bufsize_b = ldb*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	
	if (transr == 'T')
	{	transr = 'C';
		trans = transr;
	} 
	bufsize = (n*(n+1)/2);
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) ||
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "tfsm_scomplex_parameters object: malloc error.";
		tfsm_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_scomplex_buffer_pair_rand( b, bref, bufsize_b);
} /* end of Constructor  */

/* Destructor definition  'scomplex_common_parameters' class  */
tfsm_scomplex_parameters :: ~tfsm_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfsm_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfsm_free();

}
/*  Test fixture class definition */
class ctfsm_test  : public  ::testing::Test {
public:
   tfsm_scomplex_parameters  *ctfsm_obj;
   void SetUp();
   void TearDown () { delete ctfsm_obj; }
};

void ctfsm_test::SetUp(){

    /* LAPACKE ctfsm prototype */
    typedef int (*Fptr_NL_LAPACKE_ctfsm) (int matrix_layout , char transr , char side , char uplo , char trans ,\
	char diag , lapack_int m , lapack_int n , lapack_complex_float alpha , const lapack_complex_float * a , lapack_complex_float * b , lapack_int ldb);

    Fptr_NL_LAPACKE_ctfsm ctfsm;

    ctfsm_obj = new tfsm_scomplex_parameters ( eig_paramslist[idx].matrix_layout,
						lin_solver_paramslist[idx].transr,
						eig_paramslist[idx].side,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].diag,
						eig_paramslist[idx].n,
						eig_paramslist[idx].m, 
						eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    ctfsm_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ctfsm_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ctfsm_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ctfsm_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ctfsm = (Fptr_NL_LAPACKE_ctfsm)dlsym(ctfsm_obj->hModule, "LAPACKE_ctfsm");
    ASSERT_TRUE(ctfsm != NULL) << "failed to get the Netlib LAPACKE_ctfsm symbol";	
	
    ctfsm_obj->inforef = ctfsm( ctfsm_obj->matrix_layout, ctfsm_obj->transr, ctfsm_obj->side, ctfsm_obj->uplo, 
								ctfsm_obj->trans, ctfsm_obj->diag, ctfsm_obj->m, ctfsm_obj->n, ctfsm_obj->alpha, (const lapack_complex_float *)ctfsm_obj->Aref,
								ctfsm_obj->bref, ctfsm_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    ctfsm_obj->info = LAPACKE_ctfsm( ctfsm_obj->matrix_layout, ctfsm_obj->transr, ctfsm_obj->side,ctfsm_obj->uplo,
										ctfsm_obj->trans, ctfsm_obj->diag, ctfsm_obj->m, ctfsm_obj->n, ctfsm_obj->alpha, (const lapack_complex_float*)ctfsm_obj->A, 
										ctfsm_obj->b, ctfsm_obj->ldb);

    if( ctfsm_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ctfsm is wrong\n", ctfsm_obj->info );
    }
    if( ctfsm_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ctfsm is wrong\n", 
        ctfsm_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ctfsm_obj->diff =  computeDiff_c( ctfsm_obj->bufsize_b, 
                ctfsm_obj->b, ctfsm_obj->bref );

}

TEST_F(ctfsm_test, ctfsm1) {
    EXPECT_NEAR(0.0, ctfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctfsm_test, ctfsm2) {
    EXPECT_NEAR(0.0, ctfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctfsm_test, ctfsm3) {
    EXPECT_NEAR(0.0, ctfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ctfsm_test, ctfsm4) {
    EXPECT_NEAR(0.0, ctfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class tfsm_dcomplex_parameters{

   public:
	int bufsize;
	int bufsize_b;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;
	char uplo;
	char diag;
	char trans, transr;
	char side;
	lapack_complex_double* A;
	lapack_complex_double alpha; 
	lapack_int  ldb=0;
	/*Output Parameter*/
	lapack_complex_double *b;
	lapack_complex_double *Aref, *bref;
	/*Return Values*/	
	lapack_int info, inforef;

   public:
      tfsm_dcomplex_parameters (int matrix_layout, char transr, char side, char uplo,  char diag, lapack_int m ,lapack_int n, lapack_complex_double alpha);
      ~tfsm_dcomplex_parameters ();

};

/* Constructor definition  dcomplex_common_parameters */
tfsm_dcomplex_parameters::tfsm_dcomplex_parameters (int matrix_layout_i, char transr_i, char side_i, char uplo_i, \ 
char diag_i, lapack_int m_i ,lapack_int n_i, lapack_complex_double alpha_i)
{
	
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	transr = transr_i;
	trans = transr;
	uplo = uplo_i;
	diag = diag_i;
	alpha = 1;
	side = side_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n tfsm dcomplex:   n: %d, uplo:%c, trans:%c\n",  n, uplo, trans);
	#endif
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{	ldb = m;
		bufsize_b = ldb*n;
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	ldb = n;
		bufsize_b = ldb*m;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	
	if (transr == 'T')
	{	transr = 'C';
		trans = transr;
	}
	bufsize = (n*(n+1)/2);
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&b, &bref, bufsize_b);
	if ((A==NULL) || (Aref==NULL) ||
		(b==NULL) || (bref==NULL)){
		EXPECT_FALSE( true) << "tfsm_dcomplex_parameters object: malloc error.";
		tfsm_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_dcomplex_buffer_pair_rand( b, bref, bufsize_b);
} /* end of Constructor  */

/* Destructor definition  'dcomplex_common_parameters' class  */
tfsm_dcomplex_parameters :: ~tfsm_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" tfsm_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   tfsm_free();

}
/*  Test fixture class definition */
class ztfsm_test  : public  ::testing::Test {
public:
   tfsm_dcomplex_parameters  *ztfsm_obj;
   void SetUp();
   void TearDown () { delete ztfsm_obj; }
};

void ztfsm_test::SetUp(){

    /* LAPACKE ztfsm prototype */
    typedef int (*Fptr_NL_LAPACKE_ztfsm) (int matrix_layout , char transr , char side , char uplo , char trans ,\
	char diag , lapack_int m , lapack_int n , lapack_complex_double alpha , const lapack_complex_double * a , lapack_complex_double * b , lapack_int ldb);

    Fptr_NL_LAPACKE_ztfsm ztfsm;

    ztfsm_obj = new tfsm_dcomplex_parameters ( eig_paramslist[idx].matrix_layout,
						lin_solver_paramslist[idx].transr,
						eig_paramslist[idx].side,
						eig_paramslist[idx].uplo,
						eig_paramslist[idx].diag,
						eig_paramslist[idx].n,
						eig_paramslist[idx].m, 
						eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    ztfsm_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    ztfsm_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(ztfsm_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(ztfsm_obj->hModule != NULL) << "Netlib lapacke handle NULL";


    ztfsm = (Fptr_NL_LAPACKE_ztfsm)dlsym(ztfsm_obj->hModule, "LAPACKE_ztfsm");
    ASSERT_TRUE(ztfsm != NULL) << "failed to get the Netlib LAPACKE_ztfsm symbol";	
	
    ztfsm_obj->inforef = ztfsm( ztfsm_obj->matrix_layout, ztfsm_obj->transr, ztfsm_obj->side, ztfsm_obj->uplo, 
								ztfsm_obj->trans, ztfsm_obj->diag, ztfsm_obj->m, ztfsm_obj->n, ztfsm_obj->alpha, (const lapack_complex_double *)ztfsm_obj->Aref,
								ztfsm_obj->bref, ztfsm_obj->ldb);

    /* Compute libflame's Lapacke o/p  */
    ztfsm_obj->info = LAPACKE_ztfsm( ztfsm_obj->matrix_layout, ztfsm_obj->transr, ztfsm_obj->side,ztfsm_obj->uplo,
										ztfsm_obj->trans, ztfsm_obj->diag, ztfsm_obj->m, ztfsm_obj->n, ztfsm_obj->alpha, (const lapack_complex_double*)ztfsm_obj->A, 
										ztfsm_obj->b, ztfsm_obj->ldb);

    if( ztfsm_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_ztfsm is wrong\n", ztfsm_obj->info );
    }
    if( ztfsm_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_ztfsm is wrong\n", 
        ztfsm_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    ztfsm_obj->diff =  computeDiff_z( ztfsm_obj->bufsize_b, 
                ztfsm_obj->b, ztfsm_obj->bref );

}

TEST_F(ztfsm_test, ztfsm1) {
    EXPECT_NEAR(0.0, ztfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztfsm_test, ztfsm2) {
    EXPECT_NEAR(0.0, ztfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztfsm_test, ztfsm3) {
    EXPECT_NEAR(0.0, ztfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}

TEST_F(ztfsm_test, ztfsm4) {
    EXPECT_NEAR(0.0, ztfsm_obj->diff, LAPACKE_EIG_THRESHOLD);
}
