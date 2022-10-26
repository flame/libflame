#include "gtest/gtest.h"
#include "../lapacke_gtest_main.h"
#include "../lapacke_gtest_helper.h"

#define syr_free() \
if (A!=NULL)    free(A); \
if (Aref!=NULL) free(Aref);\
if (x!=NULL)    free(x); \
if (xref!=NULL)    free(xref); \
if( hModule != NULL) dlclose(hModule); \
if(dModule != NULL) dlclose(dModule)
	
// index of the 'config parameter struct' array to enable multiple sub-tests.
static int idx = 0;
/* Begin scomplex_common_parameters  class definition */
class syr_scomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_x;
	void *hModule, *dModule;
	float diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n, incx;
	lapack_complex_float* A;
	lapack_complex_float* x;
	lapack_int lda;
	char uplo;	
	lapack_complex_float alpha;
	/*Output Parameter*/
	lapack_complex_float *Aref, *xref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      syr_scomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_complex_float alpha, lapack_int incx);
      ~syr_scomplex_parameters ();

};

/* Constructor definition  float_common_parameters */
syr_scomplex_parameters:: syr_scomplex_parameters (int matrix_layout_i, char uplo_i,  lapack_int n_i, lapack_complex_float alpha_i,lapack_int incx_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	alpha = alpha_i;
	incx = incx_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n syr scomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif
	
	lda = n;
	/*Matrix sizes*/
	bufsize_x = (1+(n-1)*abs(incx));
	bufsize_a = lda*n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_scomplex_buffer_pair(&x, &xref, bufsize_x);	
	if ((A==NULL) || (Aref==NULL)||\
		(x == NULL) || (xref == NULL)){
		EXPECT_FALSE( true) << "syr_float_parameters object: malloc error.";
		syr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_scomplex_buffer_pair_rand_custom_matrix( A, Aref,lda, n, uplo);
	
	for (idx =0; idx <bufsize_x; idx++)
	{
		x[idx] = x[idx]+1;
		xref[idx] = x[idx];
	}

} /* end of Constructor  */

/* Destructor definition  'float_common_parameters' class  */
syr_scomplex_parameters :: ~syr_scomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" syr_scomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   syr_free();

}
/*  Test fixture class definition */
class csyr_test  : public  ::testing::Test {
public:
   syr_scomplex_parameters  *csyr_obj;
   void SetUp();
   void TearDown () { delete csyr_obj;}
};

void csyr_test::SetUp(){

	 /* LAPACKE csyr prototype */
    typedef int (*Fptr_NL_LAPACKE_csyr) (int matrix_layout, char uplo, lapack_int n, lapack_complex_float alpha,\
	const lapack_complex_float * x, lapack_int incx, lapack_complex_float * a, lapack_int lda);

    Fptr_NL_LAPACKE_csyr csyr;

    csyr_obj = new syr_scomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);
						   

    idx = Circular_Increment_Index(idx);

    csyr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    csyr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(csyr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(csyr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*csyr library call */
    csyr = (Fptr_NL_LAPACKE_csyr)dlsym(csyr_obj->hModule, "LAPACKE_csyr");
    ASSERT_TRUE(csyr != NULL) << "failed to get the Netlib LAPACKE_csyr symbol";
    
/*Compute csyr's  o/p */
    csyr_obj->inforef = csyr( csyr_obj->matrix_layout,  csyr_obj->uplo, csyr_obj->n, csyr_obj->alpha,
								(const lapack_complex_float*)csyr_obj->xref, csyr_obj->incx, csyr_obj->Aref, csyr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    csyr_obj->info = LAPACKE_csyr( csyr_obj->matrix_layout,csyr_obj->uplo,csyr_obj->n,csyr_obj->alpha, 
										(const lapack_complex_float*)csyr_obj->x, csyr_obj->incx, csyr_obj->A,  csyr_obj->lda);

    if( csyr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_csyr is wrong\n", csyr_obj->info );
    }
    if( csyr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_csyr is wrong\n", 
        csyr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    csyr_obj->diff =  computeDiff_c( csyr_obj->bufsize_a, 
                csyr_obj->A, csyr_obj->Aref );

}

TEST_F(csyr_test, csyr1) {
    EXPECT_NEAR(0.0, csyr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyr_test, csyr2) {
    EXPECT_NEAR(0.0, csyr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyr_test, csyr3) {
    EXPECT_NEAR(0.0, csyr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(csyr_test, csyr4) {
    EXPECT_NEAR(0.0, csyr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

/* Begin dcomplex_common_parameters  class definition */
class syr_dcomplex_parameters{

   public:
	int bufsize_a;
	int bufsize_x;
	void *hModule, *dModule;
	double diff;
   /*input parameters */
	int matrix_layout;
	lapack_int n, incx;
	lapack_complex_double* A;
	lapack_complex_double* x;
	lapack_int lda;
	char uplo;	
	lapack_complex_double alpha;
	/*Output Parameter*/
	lapack_complex_double *Aref, *xref;
	/*Return Values*/	
	lapack_int info, inforef;
   public:
      syr_dcomplex_parameters (int matrix_layout, char uplo, lapack_int n, lapack_complex_double alpha, lapack_int incx);
      ~syr_dcomplex_parameters ();

};

/* Constructor definition  double_common_parameters */
syr_dcomplex_parameters:: syr_dcomplex_parameters (int matrix_layout_i, char uplo_i,  lapack_int n_i, lapack_complex_double alpha_i,lapack_int incx_i)
{
	int i;
	matrix_layout = matrix_layout_i;
	n = n_i;
	uplo = uplo_i;
	alpha = alpha_i;
	incx = incx_i;
	#if LAPACKE_TEST_VERBOSE
	printf(" \n syr dcomplex: matrix_layout = %d, uplo:%c, n: %d \n", matrix_layout, uplo, n);
	#endif
	
	lda = n;
	/*Matrix sizes*/
	bufsize_x = (1+(n-1)*abs(incx));
	bufsize_a = lda*n;
	/*Memory allocation */
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&A, &Aref, bufsize_a);
	lapacke_gtest_alloc_lapack_dcomplex_buffer_pair(&x, &xref, bufsize_x);	
	if ((A==NULL) || (Aref==NULL)||\
		(x == NULL) || (xref == NULL)){
		EXPECT_FALSE( true) << "syr_double_parameters object: malloc error.";
		syr_free();
		exit(0);
	}
	/* Initialization of input matrices */
	lapacke_gtest_init_dcomplex_buffer_pair_rand_custom_matrix( A, Aref,lda, n, uplo);
	
	for (idx =0; idx <bufsize_x; idx++)
	{
		x[idx] = x[idx]+1;
		xref[idx] = x[idx];
	}

} /* end of Constructor  */

/* Destructor definition  'double_common_parameters' class  */
syr_dcomplex_parameters :: ~syr_dcomplex_parameters ()
{
	#if LAPACKE_TEST_VERBOSE
	printf(" syr_dcomplex_parameters object: destructor invoked. \n");
	#endif
   /* De-Allocate memory for the input matrices */
   syr_free();

}
/*  Test fixture class definition */
class zsyr_test  : public  ::testing::Test {
public:
   syr_dcomplex_parameters  *zsyr_obj;
   void SetUp();
   void TearDown () { delete zsyr_obj;}
};

void zsyr_test::SetUp(){

	 /* LAPACKE zsyr prototype */
    typedef int (*Fptr_NL_LAPACKE_zsyr) (int matrix_layout, char uplo, lapack_int n, lapack_complex_double alpha,\
	const lapack_complex_double *x, lapack_int incx, lapack_complex_double *a, lapack_int lda);

    Fptr_NL_LAPACKE_zsyr zsyr;

    zsyr_obj = new syr_dcomplex_parameters ( lin_solver_paramslist[idx].matrix_layout,
                           lin_solver_paramslist[idx].Uplo,
						   lin_solver_paramslist[idx].n,
						   lin_solver_paramslist[idx].kl,
						   lin_solver_paramslist[idx].ku);
						   

    idx = Circular_Increment_Index(idx);

    zsyr_obj->dModule = dlopen(NETLIB_BLAS_LIB, RTLD_NOW | RTLD_GLOBAL);
    zsyr_obj->hModule = dlopen(NETLIB_LAPACKE_LIB, RTLD_NOW);

    ASSERT_TRUE(zsyr_obj->dModule != NULL) << "Netlib Blas handle NULL";
    ASSERT_TRUE(zsyr_obj->hModule != NULL) << "Netlib lapacke handle NULL";

	/*zsyr library call */
    zsyr = (Fptr_NL_LAPACKE_zsyr)dlsym(zsyr_obj->hModule, "LAPACKE_zsyr");
    ASSERT_TRUE(zsyr != NULL) << "failed to get the Netlib LAPACKE_zsyr symbol";
    
/*Compute zsyr's  o/p */
    zsyr_obj->inforef = zsyr( zsyr_obj->matrix_layout,  zsyr_obj->uplo, zsyr_obj->n, zsyr_obj->alpha,
								(const lapack_complex_double*)zsyr_obj->xref, zsyr_obj->incx, zsyr_obj->Aref, zsyr_obj->lda);

    /* Compute libflame's Lapacke o/p  */
    zsyr_obj->info = LAPACKE_zsyr( zsyr_obj->matrix_layout,zsyr_obj->uplo,zsyr_obj->n,zsyr_obj->alpha, 
										(const lapack_complex_double*)zsyr_obj->x, zsyr_obj->incx, zsyr_obj->A,  zsyr_obj->lda);

    if( zsyr_obj->info < 0 ) {
        printf( "\n warning: The i:%d th argument with libflame \
        LAPACKE_zsyr is wrong\n", zsyr_obj->info );
    }
    if( zsyr_obj->inforef < 0 ) {
        printf( "The i:%d th argument with Netlib LAPACKE_zsyr is wrong\n", 
        zsyr_obj->inforef );
    }

    /* Compute Difference between libflame and Netlib o/ps  */
    zsyr_obj->diff =  computeDiff_z( zsyr_obj->bufsize_a, 
                zsyr_obj->A, zsyr_obj->Aref );

}

TEST_F(zsyr_test, zsyr1) {
    EXPECT_NEAR(0.0, zsyr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyr_test, zsyr2) {
    EXPECT_NEAR(0.0, zsyr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyr_test, zsyr3) {
    EXPECT_NEAR(0.0, zsyr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}

TEST_F(zsyr_test, zsyr4) {
    EXPECT_NEAR(0.0, zsyr_obj->diff, LAPACKE_GTEST_THRESHOLD);
}
