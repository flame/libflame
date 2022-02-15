#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_aux.h"

#define laswp_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (ipiv!=NULL)    free(ipiv); \
if (ipivref!=NULL)    free(ipivref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule); \
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin float_common_parameters  class definition */
class laswp_float_parameters{

   public:
	int bufsize;
	int bufsize_ipiv;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n,m;	
	float* A;
	char uplo;
	lapack_int k1, k2, incx;
	lapack_int lda, *ipiv;
	float* Aref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      laswp_float_parameters (int matrix_layout, lapack_int n, lapack_int m, lapack_int lda, lapack_int k1, lapack_int k2, lapack_int incx);
      ~laswp_float_parameters ();

};	

/* Constructor definition  float_common_parameters */
laswp_float_parameters:: laswp_float_parameters (int matrix_layout_i, lapack_int n_i, lapack_int m_i, lapack_int lda_i, lapack_int k1_i, lapack_int k2_i, lapack_int incx_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	k1 = k1_i;
	k2 = k2_i;
	incx = incx_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n laswp float: matrix_layout = %d,  n: %d \n", matrix_layout,  n);
	#endif
	
	if (incx == 3)
		incx = 1;
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		bufsize = lda*n;
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
//		mm = ipiv[k1-(1+(k2-k1)*abs(incx))];
		bufsize = lda;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	 
	bufsize_ipiv = (k1 +((k2-k1)*(incx)));
	
	/*Memory allocation */
	lapacke_gtest_alloc_float_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, bufsize_ipiv);
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "laswp_float_parameters object: malloc error.";
		laswp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_float_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_pivot_buffer_pair_rand(ipiv, ipivref, bufsize_ipiv);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
laswp_float_parameters :: ~laswp_float_parameters ()
{ 
	#if LAPACKE_TEST_VERBOSE
	printf(" laswp_float_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   laswp_free();

}
/*  Test fixture class definition */
class slaswp_test  : public  ::testing::Test {
public:
   laswp_float_parameters  *slaswp_obj;
   void SetUp();
   void TearDown () { delete slaswp_obj;}
};

void slaswp_test::SetUp(){

 /* LAPACKE slaswp prototype */
    typedef int (*Fptr_NL_LAPACKE_slaswp) (int matrix_layout , lapack_int n , float *a , lapack_int lda,\
	lapack_int k1 , lapack_int k2 , const lapack_int * ipiv , lapack_int incx );

    Fptr_NL_LAPACKE_slaswp slaswp;

    slaswp_obj = new laswp_float_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].m,
						   lin_solver_paramslist[idx].ldb,
						   lin_solver_paramslist[idx].ku,
						   lin_solver_paramslist[idx].kl,
						   eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    slaswp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    slaswp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(slaswp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(slaswp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*slaswp library call */
    slaswp = (Fptr_NL_LAPACKE_slaswp)dlsym(slaswp_obj->hModule, "LAPACKE_slaswp");
    ASSERT_TRUE(slaswp != NULL) << "failed to get the Netlib LAPACKE_slaswp symbol";

/*Compute slaswp's  o/p */
    slaswp_obj->inforef = slaswp( slaswp_obj->matrix_layout, slaswp_obj->n, slaswp_obj->Aref, slaswp_obj->lda,\
	slaswp_obj->k1, slaswp_obj->k2, (const lapack_int*)slaswp_obj->ipivref, slaswp_obj->incx );

    /* Compute libflame's Lapacke o/p  */
	
    slaswp_obj->info = LAPACKE_slaswp( slaswp_obj->matrix_layout, slaswp_obj->n, slaswp_obj->A, slaswp_obj->lda,\
	slaswp_obj->k1, slaswp_obj->k2, (const lapack_int*)slaswp_obj->ipiv, slaswp_obj->incx );

    if( slaswp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_slaswp is wrong\n", slaswp_obj->info );
    }
    if( slaswp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_slaswp is wrong\n", 
        slaswp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    slaswp_obj->diff =  computeDiff_s( slaswp_obj->bufsize, 
                slaswp_obj->A, slaswp_obj->Aref );

}

TEST_F(slaswp_test, slaswp1) {
    EXPECT_NEAR(0.0, slaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slaswp_test, slaswp2) {
    EXPECT_NEAR(0.0, slaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slaswp_test, slaswp3) {
    EXPECT_NEAR(0.0, slaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(slaswp_test, slaswp4) {
    EXPECT_NEAR(0.0, slaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}


/* Begin double_common_parameters  class definition */
class laswp_double_parameters{

   public:
	int bufsize;
	int bufsize_ipiv;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n,m;	
	double* A;
	char uplo;
	lapack_int k1, k2, incx;
	lapack_int lda, *ipiv;
	double* Aref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      laswp_double_parameters (int matrix_layout, lapack_int n, lapack_int m, lapack_int lda, lapack_int k1, lapack_int k2, lapack_int incx);
      ~laswp_double_parameters ();

};	

/* Constructor definition  double_common_parameters */
laswp_double_parameters:: laswp_double_parameters (int matrix_layout_i, lapack_int n_i, lapack_int m_i, lapack_int lda_i, lapack_int k1_i, lapack_int k2_i, lapack_int incx_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	k1 = k1_i;
	k2 = k2_i;
	incx = incx_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n laswp double: matrix_layout = %d,  n: %d \n", matrix_layout,  n);
	#endif
	
	if (incx == 3)
		incx = 1;
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		bufsize = lda*n;
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
//		mm = ipiv[k1-(1+(k2-k1)*abs(incx))];
		bufsize = lda;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	 
	bufsize_ipiv = (k1 +((k2-k1)*(incx)));
	
	/*Memory allocation */
	lapacke_gtest_alloc_double_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, bufsize_ipiv);
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "laswp_double_parameters object: malloc error.";
		laswp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_double_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_pivot_buffer_pair_rand(ipiv, ipivref, bufsize_ipiv);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
laswp_double_parameters :: ~laswp_double_parameters ()
{ 
	#if LAPACKE_TEST_VERBOSE
	printf(" laswp_double_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   laswp_free();

}
/*  Test fixture class definition */
class dlaswp_test  : public  ::testing::Test {
public:
   laswp_double_parameters  *dlaswp_obj;
   void SetUp();
   void TearDown () { delete dlaswp_obj;}
};

void dlaswp_test::SetUp(){

 /* LAPACKE dlaswp prototype */
    typedef int (*Fptr_NL_LAPACKE_dlaswp) (int matrix_layout , lapack_int n , double *a , lapack_int lda,\
	lapack_int k1 , lapack_int k2 , const lapack_int * ipiv , lapack_int incx );

    Fptr_NL_LAPACKE_dlaswp dlaswp;

    dlaswp_obj = new laswp_double_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].m,
						   lin_solver_paramslist[idx].ldb,
						   lin_solver_paramslist[idx].ku,
						   lin_solver_paramslist[idx].kl,
						   eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    dlaswp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    dlaswp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(dlaswp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(dlaswp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*dlaswp library call */
    dlaswp = (Fptr_NL_LAPACKE_dlaswp)dlsym(dlaswp_obj->hModule, "LAPACKE_dlaswp");
    ASSERT_TRUE(dlaswp != NULL) << "failed to get the Netlib LAPACKE_dlaswp symbol";

/*Compute dlaswp's  o/p */
    dlaswp_obj->inforef = dlaswp( dlaswp_obj->matrix_layout, dlaswp_obj->n, dlaswp_obj->Aref, dlaswp_obj->lda,\
	dlaswp_obj->k1, dlaswp_obj->k2, (const lapack_int*)dlaswp_obj->ipivref, dlaswp_obj->incx );

    /* Compute libflame's Lapacke o/p  */
	
    dlaswp_obj->info = LAPACKE_dlaswp( dlaswp_obj->matrix_layout, dlaswp_obj->n, dlaswp_obj->A, dlaswp_obj->lda,\
	dlaswp_obj->k1, dlaswp_obj->k2, (const lapack_int*)dlaswp_obj->ipiv, dlaswp_obj->incx );

    if( dlaswp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_dlaswp is wrong\n", dlaswp_obj->info );
    }
    if( dlaswp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_dlaswp is wrong\n", 
        dlaswp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    dlaswp_obj->diff =  computeDiff_d( dlaswp_obj->bufsize, 
                dlaswp_obj->A, dlaswp_obj->Aref );

}

TEST_F(dlaswp_test, dlaswp1) {
    EXPECT_NEAR(0.0, dlaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlaswp_test, dlaswp2) {
    EXPECT_NEAR(0.0, dlaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlaswp_test, dlaswp3) {
    EXPECT_NEAR(0.0, dlaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(dlaswp_test, dlaswp4) {
    EXPECT_NEAR(0.0, dlaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin scomplex_common_parameters  class definition */
class laswp_scomplex_parameters{

   public:
	int bufsize;
	int bufsize_ipiv;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n,m;	
	lapack_complex_float* A;
	char uplo;
	lapack_int k1, k2, incx;
	lapack_int lda, *ipiv;
	lapack_complex_float* Aref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      laswp_scomplex_parameters (int matrix_layout, lapack_int n, lapack_int m, lapack_int lda, lapack_int k1, lapack_int k2, lapack_int incx);
      ~laswp_scomplex_parameters ();

};	

/* Constructor definition  float_common_parameters */
laswp_scomplex_parameters:: laswp_scomplex_parameters (int matrix_layout_i, lapack_int n_i, lapack_int m_i, lapack_int lda_i, lapack_int k1_i, lapack_int k2_i, lapack_int incx_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	k1 = k1_i;
	k2 = k2_i;
	incx = incx_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n laswp scomplex: matrix_layout = %d,  n: %d \n", matrix_layout,  n);
	#endif
	
	if (incx == 3)
		incx = 1;
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		bufsize = lda*n;
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
//		mm = ipiv[k1-(1+(k2-k1)*abs(incx))];
		bufsize = lda;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	 
	bufsize_ipiv = (k1 +((k2-k1)*(incx)));
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, bufsize_ipiv);
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "laswp_float_parameters object: malloc error.";
		laswp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_pivot_buffer_pair_rand(ipiv, ipivref, bufsize_ipiv);

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
laswp_scomplex_parameters :: ~laswp_scomplex_parameters ()
{ 
	#if LAPACKE_TEST_VERBOSE
	printf(" laswp_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   laswp_free();

}
/*  Test fixture class definition */
class claswp_test  : public  ::testing::Test {
public:
   laswp_scomplex_parameters  *claswp_obj;
   void SetUp();
   void TearDown () { delete claswp_obj;}
};

void claswp_test::SetUp(){

 /* LAPACKE claswp prototype */
    typedef int (*Fptr_NL_LAPACKE_claswp) (int matrix_layout , lapack_int n , lapack_complex_float *a , lapack_int lda,\
	lapack_int k1 , lapack_int k2 , const lapack_int * ipiv , lapack_int incx );

    Fptr_NL_LAPACKE_claswp claswp;

    claswp_obj = new laswp_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].m,
						   lin_solver_paramslist[idx].ldb,
						   lin_solver_paramslist[idx].ku,
						   lin_solver_paramslist[idx].kl,
						   eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    claswp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    claswp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(claswp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(claswp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*claswp library call */
    claswp = (Fptr_NL_LAPACKE_claswp)dlsym(claswp_obj->hModule, "LAPACKE_claswp");
    ASSERT_TRUE(claswp != NULL) << "failed to get the Netlib LAPACKE_claswp symbol";

/*Compute claswp's  o/p */
    claswp_obj->inforef = claswp( claswp_obj->matrix_layout, claswp_obj->n, claswp_obj->Aref, claswp_obj->lda,\
	claswp_obj->k1, claswp_obj->k2, (const lapack_int*)claswp_obj->ipivref, claswp_obj->incx );

    /* Compute libflame's Lapacke o/p  */
	
    claswp_obj->info = LAPACKE_claswp( claswp_obj->matrix_layout, claswp_obj->n, claswp_obj->A, claswp_obj->lda,\
	claswp_obj->k1, claswp_obj->k2, (const lapack_int*)claswp_obj->ipiv, claswp_obj->incx );

    if( claswp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_claswp is wrong\n", claswp_obj->info );
    }
    if( claswp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_claswp is wrong\n", 
        claswp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    claswp_obj->diff =  computeDiff_c( claswp_obj->bufsize, 
                claswp_obj->A, claswp_obj->Aref );

}

TEST_F(claswp_test, claswp1) {
    EXPECT_NEAR(0.0, claswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(claswp_test, claswp2) {
    EXPECT_NEAR(0.0, claswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(claswp_test, claswp3) {
    EXPECT_NEAR(0.0, claswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(claswp_test, claswp4) {
    EXPECT_NEAR(0.0, claswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class laswp_dcomplex_parameters{

   public:
	int bufsize;
	int bufsize_ipiv;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n,m;	
	lapack_complex_double* A;
	char uplo;
	lapack_int k1, k2, incx;
	lapack_int lda, *ipiv;
	lapack_complex_double* Aref;
	lapack_int *ipivref;
	/*Return Values*/
	lapack_int info, inforef;

   public:
      laswp_dcomplex_parameters (int matrix_layout, lapack_int n, lapack_int m, lapack_int lda, lapack_int k1, lapack_int k2, lapack_int incx);
      ~laswp_dcomplex_parameters ();

};	

/* Constructor definition  double_common_parameters */
laswp_dcomplex_parameters:: laswp_dcomplex_parameters (int matrix_layout_i, lapack_int n_i, lapack_int m_i, lapack_int lda_i, lapack_int k1_i, lapack_int k2_i, lapack_int incx_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	m = m_i;
	lda = lda_i;
	k1 = k1_i;
	k2 = k2_i;
	incx = incx_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n laswp dcomplex: matrix_layout = %d,  n: %d \n", matrix_layout,  n);
	#endif
	
	if (incx == 3)
		incx = 1;
	
	if (matrix_layout == LAPACK_COL_MAJOR)
	{
		bufsize = lda*n;
		
	}else if (matrix_layout == LAPACK_ROW_MAJOR)
	{	
//		mm = ipiv[k1-(1+(k2-k1)*abs(incx))];
		bufsize = lda;
	}else
		EXPECT_TRUE(false) << "matrix_layout invalid";
	 
	bufsize_ipiv = (k1 +((k2-k1)*(incx)));
	
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize);
	lapacke_gtest_alloc_int_buffer_pair(&ipiv, &ipivref, bufsize_ipiv);
	if ((A==NULL) || (Aref==NULL) || \
		(ipiv == NULL) || (ipivref == NULL)) {	
		EXPECT_FALSE( true) << "laswp_double_parameters object: malloc error.";
		laswp_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand( A, Aref, bufsize);
	lapacke_gtest_init_pivot_buffer_pair_rand(ipiv, ipivref, bufsize_ipiv);

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
laswp_dcomplex_parameters :: ~laswp_dcomplex_parameters ()
{ 
	#if LAPACKE_TEST_VERBOSE
	printf(" laswp_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   laswp_free();

}
/*  Test fixture class definition */
class zlaswp_test  : public  ::testing::Test {
public:
   laswp_dcomplex_parameters  *zlaswp_obj;
   void SetUp();
   void TearDown () { delete zlaswp_obj;}
};

void zlaswp_test::SetUp(){

 /* LAPACKE zlaswp prototype */
    typedef int (*Fptr_NL_LAPACKE_zlaswp) (int matrix_layout , lapack_int n , lapack_complex_double *a , lapack_int lda,\
	lapack_int k1 , lapack_int k2 , const lapack_int * ipiv , lapack_int incx );

    Fptr_NL_LAPACKE_zlaswp zlaswp;

    zlaswp_obj = new laswp_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].m,
						   lin_solver_paramslist[idx].ldb,
						   lin_solver_paramslist[idx].ku,
						   lin_solver_paramslist[idx].kl,
						   eig_paramslist[idx].itype);

    idx = Circular_Increment_Index(idx);

    zlaswp_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zlaswp_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zlaswp_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zlaswp_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zlaswp library call */
    zlaswp = (Fptr_NL_LAPACKE_zlaswp)dlsym(zlaswp_obj->hModule, "LAPACKE_zlaswp");
    ASSERT_TRUE(zlaswp != NULL) << "failed to get the Netlib LAPACKE_zlaswp symbol";

/*Compute zlaswp's  o/p */
    zlaswp_obj->inforef = zlaswp( zlaswp_obj->matrix_layout, zlaswp_obj->n, zlaswp_obj->Aref, zlaswp_obj->lda,\
	zlaswp_obj->k1, zlaswp_obj->k2, (const lapack_int*)zlaswp_obj->ipivref, zlaswp_obj->incx );

    /* Compute libflame's Lapacke o/p  */
	
    zlaswp_obj->info = LAPACKE_zlaswp( zlaswp_obj->matrix_layout, zlaswp_obj->n, zlaswp_obj->A, zlaswp_obj->lda,\
	zlaswp_obj->k1, zlaswp_obj->k2, (const lapack_int*)zlaswp_obj->ipiv, zlaswp_obj->incx );

    if( zlaswp_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zlaswp is wrong\n", zlaswp_obj->info );
    }
    if( zlaswp_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zlaswp is wrong\n", 
        zlaswp_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zlaswp_obj->diff =  computeDiff_z( zlaswp_obj->bufsize, 
                zlaswp_obj->A, zlaswp_obj->Aref );

}

TEST_F(zlaswp_test, zlaswp1) {
    EXPECT_NEAR(0.0, zlaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlaswp_test, zlaswp2) {
    EXPECT_NEAR(0.0, zlaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlaswp_test, zlaswp3) {
    EXPECT_NEAR(0.0, zlaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zlaswp_test, zlaswp4) {
    EXPECT_NEAR(0.0, zlaswp_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
