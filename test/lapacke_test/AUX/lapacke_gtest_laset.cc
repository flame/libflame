#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define laset_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
class laset_float_parameters{

   public:
	int bufsize_a;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;	
	float* A;
	lapack_int lda;
	char uplo;	
	float alpha, beta;
	/*Output Parameter*/
	float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      laset_float_parameters (int matrix_layout, char uplo, lapack_int m, lapack_int n, float alpha, float beta);
      ~laset_float_parameters ();

};

/* Constructor definition  float_common_parameters */
laset_float_parameters:: laset_float_parameters (int matrix_layout_i, char uplo_i, lapack_int m_i, lapack_int n_i, float alpha_i, float beta_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	alpha = alpha_i;
	beta = beta_i;
	m = m_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n laset float: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		lda = m;
		bufsize_a = lda*n;
			
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		lda = n;
		bufsize_a = lda*m;
		
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize_a);	
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "laset_float_parameters object: malloc error.";
		laset_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
laset_float_parameters :: ~laset_float_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" laset_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   laset_free();

}
/*  Test fixture class definition */
class slaset_test  : public  ::testing::Test {
public:
   laset_float_parameters  *slaset_obj;
   void SetUp();
   void TearDown () { delete slaset_obj;}
};

void slaset_test::SetUp(){

	 /* LAPACKE slaset prototype */
    typedef int (*Fptr_NL_LAPACKE_slaset) (int matrix_layout,  char uplo, lapack_int m, lapack_int n, float alpha,  float beta , float *A, lapack_int lda);

    Fptr_NL_LAPACKE_slaset slaset;

    slaset_obj = new laset_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].m,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);
						   

    idx = Circular_Increment_Index(idx);

    slaset_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slaset_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slaset_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slaset_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*slaset library call */
    slaset = (Fptr_NL_LAPACKE_slaset)dlsym(slaset_obj->hModule, "LAPACKE_slaset");
    ASSERT_TRUE(slaset != NULL) << "failed to get the Netlib LAPACKE_slaset symbol";
    
/*Compute slaset's  o/p */
    slaset_obj->inforef = slaset( slaset_obj->matrix_layout,  slaset_obj->uplo, slaset_obj->m, slaset_obj->n,
								slaset_obj->alpha, slaset_obj->beta, slaset_obj->Aref, slaset_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    slaset_obj->info = LAPACKE_slaset( slaset_obj->matrix_layout,slaset_obj->uplo,slaset_obj->m,slaset_obj->n, 
										slaset_obj->alpha, slaset_obj->beta, slaset_obj->A,  slaset_obj->lda);

    if( slaset_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slaset is wrong\n", slaset_obj->info );
    }
    if( slaset_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slaset is wrong\n", 
        slaset_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slaset_obj->diff =  computeDiff_s( slaset_obj->bufsize_a, 
                slaset_obj->A, slaset_obj->Aref );

}

TEST_F(slaset_test, slaset1) {
    EXPECT_NEAR(0.0, slaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slaset_test, slaset2) {
    EXPECT_NEAR(0.0, slaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slaset_test, slaset3) {
    EXPECT_NEAR(0.0, slaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slaset_test, slaset4) {
    EXPECT_NEAR(0.0, slaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin double_common_parameters  class definition */
class laset_double_parameters{

   public:
	int bufsize_a;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;	
	double* A;
	lapack_int lda;
	char uplo;	
	double alpha, beta;
	/*Output Parameter*/
	double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      laset_double_parameters (int matrix_layout, char uplo, lapack_int m, lapack_int n, double alpha, double beta);
      ~laset_double_parameters ();

};

/* Constructor definition  double_common_parameters */
laset_double_parameters:: laset_double_parameters (int matrix_layout_i, char uplo_i, lapack_int m_i, lapack_int n_i, double alpha_i, double beta_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	alpha = alpha_i;
	beta = beta_i;
	m = m_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n laset double: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		lda = m;
		bufsize_a = lda*n;
			
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		lda = n;
		bufsize_a = lda*m;
		
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize_a);	
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "laset_double_parameters object: malloc error.";
		laset_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
laset_double_parameters :: ~laset_double_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" laset_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   laset_free();

}
/*  Test fixture class definition */
class dlaset_test  : public  ::testing::Test {
public:
   laset_double_parameters  *dlaset_obj;
   void SetUp();
   void TearDown () { delete dlaset_obj;}
};

void dlaset_test::SetUp(){

	 /* LAPACKE dlaset prototype */
    typedef int (*Fptr_NL_LAPACKE_dlaset) (int matrix_layout,  char uplo, lapack_int m, lapack_int n, double alpha,  double beta , double *A, lapack_int lda);

    Fptr_NL_LAPACKE_dlaset dlaset;

    dlaset_obj = new laset_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].m,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);
						   

    idx = Circular_Increment_Index(idx);

    dlaset_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlaset_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlaset_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlaset_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dlaset library call */
    dlaset = (Fptr_NL_LAPACKE_dlaset)dlsym(dlaset_obj->hModule, "LAPACKE_dlaset");
    ASSERT_TRUE(dlaset != NULL) << "failed to get the Netlib LAPACKE_dlaset symbol";
    
/*Compute dlaset's  o/p */
    dlaset_obj->inforef = dlaset( dlaset_obj->matrix_layout,  dlaset_obj->uplo, dlaset_obj->m, dlaset_obj->n,
								dlaset_obj->alpha, dlaset_obj->beta, dlaset_obj->Aref, dlaset_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    dlaset_obj->info = LAPACKE_dlaset( dlaset_obj->matrix_layout,dlaset_obj->uplo,dlaset_obj->m,dlaset_obj->n, 
										dlaset_obj->alpha, dlaset_obj->beta, dlaset_obj->A,  dlaset_obj->lda);

    if( dlaset_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlaset is wrong\n", dlaset_obj->info );
    }
    if( dlaset_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlaset is wrong\n", 
        dlaset_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlaset_obj->diff =  computeDiff_d( dlaset_obj->bufsize_a, 
                dlaset_obj->A, dlaset_obj->Aref );

}

TEST_F(dlaset_test, dlaset1) {
    EXPECT_NEAR(0.0, dlaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlaset_test, dlaset2) {
    EXPECT_NEAR(0.0, dlaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlaset_test, dlaset3) {
    EXPECT_NEAR(0.0, dlaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlaset_test, dlaset4) {
    EXPECT_NEAR(0.0, dlaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin scomplex_common_parameters  class definition */
class laset_scomplex_parameters{

   public:
	int bufsize_a;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;	
	lapack_complex_float* A;
	lapack_int lda;
	char uplo;	
	lapack_complex_float alpha, beta;
	/*Output Parameter*/
	lapack_complex_float *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      laset_scomplex_parameters (int matrix_layout, char uplo, lapack_int m, lapack_int n, lapack_complex_float alpha, lapack_complex_float beta);
      ~laset_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
laset_scomplex_parameters:: laset_scomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int m_i, lapack_int n_i, lapack_complex_float alpha_i,lapack_complex_float beta_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	alpha = alpha_i;
	beta = beta_i;
	m = m_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n laset scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		lda = m;
		bufsize_a = lda*n;
			
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		lda = n;
		bufsize_a = lda*m;
		
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);	
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "laset_float_parameters object: malloc error.";
		laset_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
laset_scomplex_parameters :: ~laset_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" laset_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   laset_free();

}
/*  Test fixture class definition */
class claset_test  : public  ::testing::Test {
public:
   laset_scomplex_parameters  *claset_obj;
   void SetUp();
   void TearDown () { delete claset_obj;}
};

void claset_test::SetUp(){

	 /* LAPACKE claset prototype */
    typedef int (*Fptr_NL_LAPACKE_claset) (int matrix_layout,  char uplo, lapack_int m, lapack_int n, lapack_complex_float alpha,  lapack_complex_float beta , lapack_complex_float *A, lapack_int lda);

    Fptr_NL_LAPACKE_claset claset;

    claset_obj = new laset_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].m,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);
						   

    idx = Circular_Increment_Index(idx);

    claset_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    claset_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(claset_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(claset_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*claset library call */
    claset = (Fptr_NL_LAPACKE_claset)dlsym(claset_obj->hModule, "LAPACKE_claset");
    ASSERT_TRUE(claset != NULL) << "failed to get the Netlib LAPACKE_claset symbol";
    
/*Compute claset's  o/p */
    claset_obj->inforef = claset( claset_obj->matrix_layout,  claset_obj->uplo, claset_obj->m, claset_obj->n,
								claset_obj->alpha, claset_obj->beta, claset_obj->Aref, claset_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    claset_obj->info = LAPACKE_claset( claset_obj->matrix_layout,claset_obj->uplo,claset_obj->m,claset_obj->n, 
										claset_obj->alpha, claset_obj->beta, claset_obj->A,  claset_obj->lda);

    if( claset_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_claset is wrong\n", claset_obj->info );
    }
    if( claset_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_claset is wrong\n", 
        claset_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    claset_obj->diff =  computeDiff_c( claset_obj->bufsize_a, 
                claset_obj->A, claset_obj->Aref );

}

TEST_F(claset_test, claset1) {
    EXPECT_NEAR(0.0, claset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(claset_test, claset2) {
    EXPECT_NEAR(0.0, claset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(claset_test, claset3) {
    EXPECT_NEAR(0.0, claset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(claset_test, claset4) {
    EXPECT_NEAR(0.0, claset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class laset_dcomplex_parameters{

   public:
	int bufsize_a;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n;
	lapack_int m;	
	lapack_complex_double* A;
	lapack_int lda;
	char uplo;	
	lapack_complex_double alpha, beta;
	/*Output Parameter*/
	lapack_complex_double *Aref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      laset_dcomplex_parameters (int matrix_layout, char uplo, lapack_int m, lapack_int n, lapack_complex_double alpha, lapack_complex_double beta);
      ~laset_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
laset_dcomplex_parameters:: laset_dcomplex_parameters (int matrix_layout_i, char uplo_i, lapack_int m_i, lapack_int n_i, lapack_complex_double alpha_i, lapack_complex_double beta_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	alpha = alpha_i;
	beta = beta_i;
	m = m_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n laset dcomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		lda = m;
		bufsize_a = lda*n;
			
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
		lda = n;
		bufsize_a = lda*m;
		
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);	
	if ((A==NULL) || (Aref==NULL)){
		EXPECT_FALSE( true) << "laset_double_parameters object: malloc error.";
		laset_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize_a);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
laset_dcomplex_parameters :: ~laset_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" laset_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   laset_free();

}
/*  Test fixture class definition */
class zlaset_test  : public  ::testing::Test {
public:
   laset_dcomplex_parameters  *zlaset_obj;
   void SetUp();
   void TearDown () { delete zlaset_obj;}
};

void zlaset_test::SetUp(){

	 /* LAPACKE zlaset prototype */
    typedef int (*Fptr_NL_LAPACKE_zlaset) (int matrix_layout,  char uplo, lapack_int m, lapack_int n, lapack_complex_double alpha,  lapack_complex_double beta , lapack_complex_double *A, lapack_int lda);

    Fptr_NL_LAPACKE_zlaset zlaset;

    zlaset_obj = new laset_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].m,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);
						   

    idx = Circular_Increment_Index(idx);

    zlaset_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlaset_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlaset_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlaset_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zlaset library call */
    zlaset = (Fptr_NL_LAPACKE_zlaset)dlsym(zlaset_obj->hModule, "LAPACKE_zlaset");
    ASSERT_TRUE(zlaset != NULL) << "failed to get the Netlib LAPACKE_zlaset symbol";
    
/*Compute zlaset's  o/p */
    zlaset_obj->inforef = zlaset( zlaset_obj->matrix_layout,  zlaset_obj->uplo, zlaset_obj->m, zlaset_obj->n,
								zlaset_obj->alpha, zlaset_obj->beta, zlaset_obj->Aref, zlaset_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    zlaset_obj->info = LAPACKE_zlaset( zlaset_obj->matrix_layout,zlaset_obj->uplo,zlaset_obj->m,zlaset_obj->n, 
										zlaset_obj->alpha, zlaset_obj->beta, zlaset_obj->A,  zlaset_obj->lda);

    if( zlaset_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlaset is wrong\n", zlaset_obj->info );
    }
    if( zlaset_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlaset is wrong\n", 
        zlaset_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlaset_obj->diff =  computeDiff_z( zlaset_obj->bufsize_a, 
                zlaset_obj->A, zlaset_obj->Aref );

}

TEST_F(zlaset_test, zlaset1) {
    EXPECT_NEAR(0.0, zlaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlaset_test, zlaset2) {
    EXPECT_NEAR(0.0, zlaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlaset_test, zlaset3) {
    EXPECT_NEAR(0.0, zlaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlaset_test, zlaset4) {
    EXPECT_NEAR(0.0, zlaset_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
